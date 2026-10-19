"""A minimal, engine-agnostic Replay Pack *producer* (a "limited generator").

Demonstrates that any simulator can emit a conformant Tier-0 (Gradient) pack with nothing
but NumPy + safetensors. Here the "walk" is analytic free diffusion (independent Gaussian
increments), accumulated into the bands of SPEC §9.5 while it runs -- no positions are kept --
so replay can be checked against the exact free-diffusion signal exp(-b D).

Each batch of walkers is accumulated at block h = the step, moment order m = 0:
    beta_k = sqrt(2T) / (pi k) * sum_j cos(pi k t_j / T) dI_j ,   t_j the step's midpoint,
and the pack records the quadratic variation of the increments (QV, SPEC §9.5.5) and the
certificate terms (truncation, quadrature) for a battery of waveforms.

Run:  python write_minimal_rpk.py  ->  writes free_diffusion.rpk
Apache-2.0.
"""
import json
import math
import numpy as np
from safetensors.numpy import save_file
import reference_replayer as rr

RPK_SCHEMA_VERSION = "0.4"


def rel(m, x):
    """SPEC §9.5.5: the relative variance figure of the block quadrature at moment order m."""
    return x ** (2 * (m + 1)) / (math.factorial(m + 1) ** 2 * (2 * m + 3) * 4 ** (m + 1))


def accumulate_bands(n_walkers, T, n_steps, D, K, seed=0, batch=2000):
    """Free-diffusion walk accumulated into bands. Returns (vec (N_w, 3, K+2), QV (3,)).
    The walkers start at the origin (SPEC §4.2: continuous lab-frame positions, metres)."""
    rng = np.random.default_rng(seed)
    h = T / n_steps
    k = np.arange(1, K + 1)
    tm = (np.arange(n_steps) + 0.5) * h
    cosk = np.cos(np.pi * np.outer(tm, k) / T) * (np.sqrt(2.0 * T) / (np.pi * k))   # (n_steps, K)
    vec = np.zeros((n_walkers, 3, K + 2))
    qv = np.zeros(3)
    for i in range(0, n_walkers, batch):
        b = min(batch, n_walkers - i)
        dI = rng.normal(0.0, np.sqrt(2.0 * D * h), size=(b, 3, n_steps))
        vec[i:i + b, :, 1] = dI.sum(axis=2)           # r(T) - r(0); r(0) = 0
        vec[i:i + b, :, 2:] = dI @ cosk               # the K bands
        qv += (dI ** 2).sum(axis=(0, 2))
    return vec, qv / n_walkers


def certificate_terms(battery, T, K, h, m, D, qv):
    """SPEC §9.5.5: the largest truncation and quadrature bounds over a battery of (G, dt_wf)."""
    k_hi = np.arange(K + 1, 20001)                    # the bands beyond K, summed far enough
    k_lo = np.arange(1, K + 1)
    trunc = quad = 0.0
    for G, dt_wf in battery:
        _, _, Gh_hi = rr.waveform_columns(G, dt_wf, T, k_hi)
        _, _, Gh_lo = rr.waveform_columns(G, dt_wf, T, k_lo)
        var_free = 2.0 * D * T ** 2 / (np.pi * k_hi) ** 2
        trunc = max(trunc, 0.5 * float((rr.GAMMA ** 2 * Gh_hi ** 2 * var_free).sum()))
        w = np.pi * k_lo / T
        weight = rr.GAMMA ** 2 * Gh_lo ** 2 * rel(m, w * h) * (T / (np.pi * k_lo) ** 2)   # (3, K)
        quad = max(quad, 0.5 * float((weight * qv[:, None]).sum()))
    return trunc, quad


def write_pack(path, vec, T, h, D, K, qv, battery, seed, pack_id):
    n_w = vec.shape[0]
    trunc, quad = certificate_terms(battery, T, K, h, 0, D, qv)
    floor = 1.0 / math.sqrt(n_w)                      # MC floor of a unit-modulus ensemble mean
    meta = {
        "rpk_schema_version": RPK_SCHEMA_VERSION,
        "id": pack_id,
        "walk_params": {
            "n_walkers": int(n_w), "T_max": float(T),
            "diffusivity": float(D), "seed": int(seed),
            "block": {"h": float(h), "m": 0, "precision": "float64"},
            "qv": {"positions": [float(q) for q in qv]},
        },
        "compression": {"method": "bridge_dst", "K": int(K), "walker_preserving": True,
                        "temporal_bandwidth_hz": K / (2.0 * T)},
        "replay_envelope": {
            "gradient": True,                 # C0 only
            "bulk_relaxation": False, "surface_relaxivity": False,
            "field": False,
            "magnetization_transfer": False, "diffusivity_fixed": True,
            "acquisition": {"b_max": 5.0e9, "ogse_periods": [], "B0_list": []},
        },
        "fidelity": {"err_max": trunc + quad, "floor_max": floor,
                     "terms": {"truncation": trunc, "truncation_from": "band_variance",
                               "quadrature": {"positions": quad}},
                     "within_2x_floor": bool(trunc + quad <= 2 * floor),
                     "battery": "PGSE along x, g = 0.02 .. 0.20 T/m, delta 10 ms, Delta 40 ms",
                     "certified": "measured", "inherited_from": None},
        "provenance": {"generator": "write_minimal_rpk.py", "generator_version": "0.2",
                       "substrate": {"representation": "analytical",
                                     "description": "free diffusion (no restriction)", "ref": None}},
        "license": "CC0-1.0",
        "citation": "Replay Pack minimal example, dmrai-lab (2026).",
    }
    arrays = {f"pos_{a}": np.ascontiguousarray(vec[:, i, :], np.float32) for i, a in enumerate("xyz")}
    save_file(arrays, path, metadata={"rpk": json.dumps(meta)})
    return meta


if __name__ == "__main__":
    D, T, K, n_steps, N = 2.0e-9, 0.05, 128, 1000, 40000    # water-like D, TE = 50 ms
    h = T / n_steps
    battery = []
    for g in (0.02, 0.05, 0.10, 0.20):
        G, _ = rr.pgse_waveform(T / 500, 0.010, 0.040, g, [1, 0, 0])
        battery.append((G, T / 500))
    vec, qv = accumulate_bands(N, T, n_steps, D, K, seed=0)
    meta = write_pack("free_diffusion.rpk", vec, T, h, D, K, qv, battery, seed=0,
                      pack_id="examples/free-diffusion")
    f = meta["fidelity"]
    print(f"wrote free_diffusion.rpk  N_w={N} T={T*1e3:.0f}ms K={K} h={h*1e6:.0f}us m=0 D={D*1e9:.1f}e-9  (Tier C0)")
    print(f"QV/(2DT) = {qv / (2 * D * T)}   truncation {f['terms']['truncation']:.2e}  "
          f"quadrature {f['terms']['quadrature']['positions']:.2e}  floor {f['floor_max']:.2e}")
