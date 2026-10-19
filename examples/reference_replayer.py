"""Reference Replay Pack replayer — dependency-light (NumPy + safetensors only).

Implements the Gradient (C0) and Bulk-relaxation (C1) replay operations of RPK.md §6.1/§6.2
against a pack's bands (`bridge_dst`, float container): the gradient phase is the closed
form of §9.5.3, the relaxation log-weight an exact sum over the runs of the occupancy
(§9.5.7). It deliberately does NOT depend on dmipy-sim, to demonstrate that the format is
generator- and engine-agnostic.

Apache-2.0.
"""
import json
import numpy as np
from safetensors import safe_open

GAMMA = 2.6751525e8  # proton gyromagnetic ratio, rad s^-1 T^-1


def load_rpk(path):
    """Return (arrays: dict[str, np.ndarray], meta: dict). Reads the JSON metadata
    from the safetensors __metadata__ header (SPEC §12)."""
    arrays = {}
    with safe_open(path, framework="numpy") as f:
        hdr = f.metadata() or {}
        for k in f.keys():
            arrays[k] = f.get_tensor(k)
    blob = hdr.get("rpk") or hdr.get("json")   # canonical header key is "rpk" (SPEC §12)
    meta = json.loads(blob) if blob else dict(hdr)
    return arrays, meta


def _bands(arrays, meta):
    """The stored vectors [r(0), r(T) - r(0), beta_1 .. beta_K] of the three axes, (N_w, 3, K+2).
    Only the float container of `bridge_dst` is handled; a full replayer dispatches on
    meta['compression']['method'] and 'container' (registry)."""
    method = meta.get("compression", {}).get("method")
    if method != "bridge_dst" or meta["compression"].get("container"):
        raise NotImplementedError(
            f"codec {method!r} / container not implemented by this reference replayer (SPEC §9).")
    return np.stack([np.asarray(arrays[k], np.float64) for k in ("pos_x", "pos_y", "pos_z")], axis=1)


def waveform_columns(G, dt_wf, T, k):
    """The waveform's columns against the basis of §9.5.3, per gradient axis.

    G : (n, 3) T/m, sample-and-hold on ITS OWN grid dt_wf (the waveform is never resampled);
    T : the pack's duration; k : the band indices (array). Returns M0 (3,), M1 (3,), Ghat (3, len(k)):
    M0 = int G dt, M1 = int (t/T) G dt, Ghat_k = int G(t) s_k(t) dt, s_k = sqrt(2/T) sin(pi k t/T)."""
    G = np.asarray(G, np.float64)
    a = np.arange(G.shape[0]) * dt_wf
    b = a + dt_wf
    assert b[-1] <= T * (1 + 1e-9), "a waveform MUST NOT extend beyond T_max (SPEC §6.1)"
    k = np.asarray(k, np.float64)
    w = np.pi * k / T
    # int_a^b sin(w t) dt = (cos(w a) - cos(w b)) / w, in closed form for every piece
    pieces = (np.cos(np.outer(a, w)) - np.cos(np.outer(b, w))) / w
    Ghat = np.sqrt(2.0 / T) * (G.T @ pieces)
    M0 = G.sum(axis=0) * dt_wf
    M1 = (G * ((b * b - a * a) / (2.0 * T))[:, None]).sum(axis=0)
    return M0, M1, Ghat


def replay_gradient(arrays, meta, G, dt_wf, weights=None):
    """SPEC §6.1 on the bands (§9.5.3): S(G) = <w e^{i phi}>,
    phi = gamma [ r(0).M0 + (r(T) - r(0)).M1 + sum_k beta_k Ghat_k ].  Returns complex S."""
    c = _bands(arrays, meta)                          # (N_w, 3, K+2)
    K = int(meta["compression"]["K"])
    assert c.shape[2] == K + 2, "the stored width is K + 2: take K from the metadata"
    T = float(meta["walk_params"]["T_max"])
    M0, M1, Ghat = waveform_columns(G, dt_wf, T, np.arange(1, K + 1))
    phi = GAMMA * (c[:, :, 0] @ M0 + c[:, :, 1] @ M1 + np.einsum("nak,ak->n", c[:, :, 2:], Ghat))
    w = np.asarray(arrays["spin_weights"], np.float64) if "spin_weights" in arrays \
        else (np.ones(c.shape[0]) if weights is None else np.asarray(weights, np.float64))
    if weights is not None:
        w = np.asarray(weights, np.float64)
    return np.sum(w * np.exp(1j * phi)) / np.sum(w)


def occupancy_integral(arrays, meta, rate, T_e=None, column="comp"):
    """int_0^{T_e} rate[label(t)] dt per walker, exact over the runs of the occupancy (SPEC §9.5.7).
    `rate` is indexed by label; the column is a static label or an event list."""
    T_max = float(meta["walk_params"]["T_max"])
    T_e = T_max if T_e is None else float(T_e)
    rate = np.asarray(rate, np.float64)
    if f"{column}_static" in arrays:
        return rate[np.asarray(arrays[f"{column}_static"], np.int64)] * T_e
    off = np.asarray(arrays[f"{column}_events_offsets"], np.int64)
    t = np.asarray(arrays[f"{column}_events_t"], np.float64)
    lab = np.asarray(arrays[f"{column}_events_label"], np.int64)
    t_next = np.empty_like(t)
    t_next[:-1] = t[1:]
    t_next[off[1:] - 1] = T_max                       # a walker holds its last label to T_max
    overlap = np.clip(np.minimum(t_next, T_e) - t, 0.0, None)
    walker = np.repeat(np.arange(off.size - 1), np.diff(off))
    return np.bincount(walker, weights=rate[lab] * overlap, minlength=off.size - 1)


def replay_relaxation_logweight(arrays, meta, T2, T_e=None):
    """SPEC §6.2: per-walker transverse log-weight -int rate dt from the compartment channel.
    Requires the Bulk-relaxation tier. `T2` is a sequence indexed by pool id (a replay knob,
    SPEC §8.5). Returns (N_w,) log-weights (<=0)."""
    env = meta.get("replay_envelope", {})
    if not (env.get("bulk_relaxation") or env.get("relaxation") or env.get("T1T2")):
        raise ValueError("pack does not declare the Bulk-relaxation tier; "
                         "capability not present (SPEC §7/§13).")
    T2 = np.asarray(T2, float)
    rate = np.where(np.isfinite(T2) & (T2 > 0), 1.0 / np.maximum(T2, 1e-30), 0.0)
    return -occupancy_integral(arrays, meta, rate, T_e)


# ---- helper to build a demo waveform on the waveform's own grid ----
def pgse_waveform(dt_wf, delta, Delta, g_amp, direction):
    """An ideal PGSE gradient, sample-and-hold on its own grid `dt_wf`: two rectangular lobes of
    duration `delta`, separated by `Delta` (onset-to-onset), opposite sign. Returns G (n, 3) in T/m
    and the analytic b-value (s/m^2). The waveform ends at Delta + delta."""
    n = int(round((Delta + delta) / dt_wf))
    t = (np.arange(n) + 0.5) * dt_wf
    d = np.asarray(direction, float); d = d / np.linalg.norm(d)
    G = np.zeros((n, 3))
    G[t < delta] = g_amp * d
    G[(t >= Delta) & (t < Delta + delta)] = -g_amp * d
    b = (GAMMA * g_amp * delta) ** 2 * (Delta - delta / 3.0)
    return G, b
