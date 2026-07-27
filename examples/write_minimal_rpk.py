"""A minimal, engine-agnostic Replay Pack *producer* (a "limited generator").

Demonstrates that any simulator can emit a conformant Tier-0 (Gradient) pack with nothing
but NumPy + safetensors. Here the "walk" is analytic free diffusion (independent Gaussian
increments), so we can check replay against the exact free-diffusion signal exp(-b D).

Run:  python write_minimal_rpk.py  ->  writes free_diffusion.rpk
Apache-2.0.
"""
import json
import numpy as np
from safetensors.numpy import save_file

RPK_SCHEMA_VERSION = "1.1"


def make_free_diffusion_walk(n_walkers, n_t, dt, D, seed=0):
    """Free-diffusion trajectories: r(t_{k+1}) = r(t_k) + N(0, 2 D dt) per axis.
    Continuous, lab-frame, metres (SPEC §4/§8.2)."""
    rng = np.random.default_rng(seed)
    steps = rng.normal(0.0, np.sqrt(2.0 * D * dt), size=(n_walkers, n_t, 3)).astype(np.float32)
    steps[:, 0, :] = 0.0                       # start all walkers at the origin
    return np.cumsum(steps, axis=1)            # (N_w, N_t, 3)


def write_pack(path, positions, dt, D, seed, pack_id):
    n_w, n_t, _ = positions.shape
    meta = {
        "rpk_schema_version": RPK_SCHEMA_VERSION,
        "id": pack_id,
        "walk_params": {
            "n_walkers": int(n_w), "n_t": int(n_t),
            "dt_traj": float(dt), "T_max": float((n_t - 1) * dt),
            "diffusivity": float(D), "seed": int(seed),
        },
        "per_comp": None,
        "compression": {"method": "identity", "K": None, "walker_preserving": True},
        "replay_envelope": {
            "gradient": True, "rf": True,
            "T1T2": False, "rho": False, "mt": False,
            "B0_any": False, "orientation_any": False,
            "permeability": False, "diffusivity_fixed": True,
            "acquisition": {"b_max": 5.0e9, "ogse_periods": [], "B0_list": []},
        },
        "fidelity": None,                       # identity codec is lossless
        "provenance": {"generator": "write_minimal_rpk.py", "generator_version": "0.1",
                       "geometry": "free-diffusion (analytic)", "real_or_synthetic": "synthetic"},
        "license": "CC0-1.0",
        "citation": "Replay Pack minimal example, dmrai-lab (2026).",
    }
    arrays = {"positions": np.ascontiguousarray(positions, np.float32)}
    save_file(arrays, path, metadata={"rpk": json.dumps(meta)})
    return meta


if __name__ == "__main__":
    D, dt, n_t = 2.0e-9, 5.0e-4, 100          # water-like D, TE ~ 50 ms
    pos = make_free_diffusion_walk(n_walkers=20000, n_t=n_t, dt=dt, D=D, seed=0)
    meta = write_pack("free_diffusion.rpk", pos, dt, D, seed=0, pack_id="examples/free-diffusion")
    print(f"wrote free_diffusion.rpk  N_w={pos.shape[0]} N_t={pos.shape[1]} "
          f"dt={dt*1e3:.2f}ms D={D*1e9:.1f}e-9  (Tier 0, identity codec)")
