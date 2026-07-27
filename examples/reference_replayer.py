"""Reference Replay Pack replayer — dependency-light (NumPy + safetensors only).

Implements the Gradient (T0) and Relaxation (T1) replay operations of SPEC.md §6.1/§6.2
against a raw (identity-codec) pack. This file is intentionally small: a conformant
replayer for the base tiers is ~40 lines. It deliberately does NOT depend on dmipy-sim,
to demonstrate that the format is generator- and engine-agnostic.

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
    # canonical header key is "rpk"; "json" is the accepted 1.x legacy alias (SPEC §12)
    blob = hdr.get("rpk") or hdr.get("json")
    meta = json.loads(blob) if blob else dict(hdr)
    return arrays, meta


def _positions(arrays, meta):
    """Decode the positions channel. Only the identity codec is handled here; a full
    replayer would dispatch on meta['compression']['method']."""
    method = meta.get("compression", {}).get("method", "identity")
    if method != "identity":
        raise NotImplementedError(
            f"codec {method!r} not implemented by this reference replayer; "
            f"decode to (N_w, N_t, 3) first (SPEC §9).")
    return np.asarray(arrays["positions"], np.float64)


def replay_gradient(arrays, meta, G, dt_wf=None, weights=None):
    """SPEC §6.1: S(G) = <w e^{i gamma dt sum_k G_k . r_k}>.

    G : (N_t, 3) gradient waveform on the save grid, T/m.
    Returns complex S.
    """
    r = _positions(arrays, meta)                     # (N_w, N_t, 3)
    dt = float(meta["walk_params"]["dt_traj"]) if dt_wf is None else float(dt_wf)
    G = np.asarray(G, np.float64)
    assert G.shape == r.shape[1:], f"G {G.shape} must be (N_t,3)={r.shape[1:]}"
    phi = GAMMA * dt * np.einsum("ntd,td->n", r, G)  # (N_w,)
    w = np.asarray(arrays["spin_weights"], np.float64) if "spin_weights" in arrays \
        else (np.ones(r.shape[0]) if weights is None else np.asarray(weights, np.float64))
    if weights is not None:
        w = np.asarray(weights, np.float64)
    return np.sum(w * np.exp(1j * phi)) / np.sum(w)


def replay_relaxation_logweight(arrays, meta, T2=None, T1=None, transverse=True):
    """SPEC §6.2: per-walker relaxation log-weight from the compartment channel.
    Requires the Relaxation tier. `T2`/`T1` are lists indexed by compartment id;
    default to `per_comp`. Returns (N_w,) log-weights (<=0)."""
    env = meta.get("replay_envelope", {})
    if not (env.get("bulk_relaxation") or env.get("relaxation") or env.get("T1T2")):  # +1.x aliases
        raise ValueError("pack does not declare the Relaxation tier; "
                         "capability not present (SPEC §7/§13).")
    comp = np.asarray(arrays["compartment"], np.int64)      # (N_w, N_t)
    pc = meta.get("per_comp", {})
    T2 = np.asarray(T2 if T2 is not None else pc["T2"], float)
    dt = float(meta["walk_params"]["dt_traj"])
    invT2 = np.where(np.isfinite(T2) & (T2 > 0), 1.0 / np.maximum(T2, 1e-30), 0.0)
    # transverse spin-echo: chi = 1 throughout the encoding
    rate = invT2[comp]                                       # (N_w, N_t)
    if not transverse and T1 is not None:
        invT1 = 1.0 / np.maximum(np.asarray(T1, float), 1e-30)
        rate = invT1[comp]
    return -dt * rate.sum(axis=1)


# ---- helpers to build a demo waveform on the save grid ----
def pgse_on_grid(n_t, dt, delta, Delta, g_amp, direction):
    """A trapezoid-free ideal PGSE gradient sampled on the save grid, for demos.
    Two rectangular lobes of duration `delta`, separated by `Delta` (onset-to-onset),
    opposite sign. Returns G (n_t, 3) in T/m and the analytic b-value (s/m^2)."""
    t = np.arange(n_t) * dt
    d = np.asarray(direction, float); d = d / np.linalg.norm(d)
    G = np.zeros((n_t, 3))
    lobe1 = t < delta
    lobe2 = (t >= Delta) & (t < Delta + delta)
    G[lobe1] = g_amp * d
    G[lobe2] = -g_amp * d
    b = (GAMMA * g_amp * delta) ** 2 * (Delta - delta / 3.0)
    return G, b
