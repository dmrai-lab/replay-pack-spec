"""Interop: read a Replay Pack produced by the reference generator (dmipy-sim) with the
standalone spec tools in this repo — no dmipy-sim import needed on the *reader* side.

This is the whole point of the standard: dmipy-sim (or any other simulator) writes an
.rpk; the generic conformance checker and reference replayer here consume it.

Usage:
    python dmipy_sim_interop.py <pack.rpk>       # validate + replay an existing pack
    python dmipy_sim_interop.py --build          # build a tiny pack via dmipy-sim, then do both

Apache-2.0.
"""
import sys
import numpy as np
import conformance_check as cc
import reference_replayer as rr


def validate_and_replay(path):
    print(f"# validating {path} with the generic spec checker")
    rc = cc.main(path)
    print(f"\n# replaying a PGSE with the generic reference replayer")
    arrays, meta = rr.load_rpk(path)
    T = meta["walk_params"]["T_max"]
    D = meta["walk_params"].get("diffusivity")
    dt_wf = T / 500
    G, b = rr.pgse_waveform(dt_wf, min(0.010, T / 5), min(0.040, 4 * T / 5), 0.05, direction=[1, 0, 0])
    S = rr.replay_gradient(arrays, meta, G, dt_wf)
    print(f"  b={b/1e6:.0f} s/mm^2   S={abs(S):.4f}"
          + (f"   (free-limit exp(-bD)={np.exp(-b*D):.4f})" if D else ""))
    return rc


def build_small_pack(path="dmipy_interop.rpk"):
    """Build a minimal pack through the reference generator (dmipy-sim). Requires the
    dmipy-sim environment; kept tiny so it runs on CPU. The walk is a coefficient walk
    (no save grid): `simulate_trajectories` takes the duration and the band, and
    `build_replay_pack` certifies it."""
    from dmipy_sim import simulate_trajectories
    from dmipy_sim.geometries import Cylinder
    from dmipy_sim import bank
    geom = Cylinder(radius=5e-6, orientation=[0, 0, 1])
    walk = simulate_trajectories(n_walkers=4000, diffusivity=1.5e-9, geometry=geom,
                                 T_max=0.03, seed=0, require_gpu=False)
    bank.build_replay_pack(walk, id="examples/dmipy-cylinder", license="CC0-1.0",
                           citation="Replay Pack interop example, dmrai-lab (2026).",
                           out_path=path, verbose=True)
    return path


if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] == "--build":
        p = build_small_pack()
        sys.exit(validate_and_replay(p))
    elif len(sys.argv) == 2:
        sys.exit(validate_and_replay(sys.argv[1]))
    else:
        print(__doc__); sys.exit(2)
