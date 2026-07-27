"""End-to-end standalone demo (no dmipy-sim, no GPU):

  produce  ->  a Tier-0 pack (write_minimal_rpk)
  validate ->  conformance_check
  replay   ->  reference_replayer, for several b-values, vs the exact free-diffusion
               signal exp(-b D). Agreement to the MC noise floor proves the format +
               reference replay operation are self-consistent and engine-agnostic.

Run:  python demo_free_diffusion.py
Apache-2.0.
"""
import numpy as np
import write_minimal_rpk as wr
import reference_replayer as rr
import conformance_check as cc

D, dt, n_t, N = 2.0e-9, 5.0e-4, 100, 40000
pos = wr.make_free_diffusion_walk(n_walkers=N, n_t=n_t, dt=dt, D=D, seed=0)
wr.write_pack("free_diffusion.rpk", pos, dt, D, seed=0, pack_id="examples/free-diffusion")

print("== conformance ==")
assert cc.main("free_diffusion.rpk") == 0

print("\n== replay vs analytic exp(-bD) ==")
arrays, meta = rr.load_rpk("free_diffusion.rpk")
delta, Delta = 0.010, 0.040
print(f"{'g[T/m]':>8} {'b[s/mm2]':>10} {'replay':>10} {'exp(-bD)':>10} {'|Δ|':>8}")
worst = 0.0
for g in (0.02, 0.05, 0.10, 0.20):
    G, b = rr.pgse_on_grid(n_t, dt, delta, Delta, g, direction=[1, 0, 0])
    S = rr.replay_gradient(arrays, meta, G).real
    exact = np.exp(-b * D)
    worst = max(worst, abs(S - exact))
    print(f"{g:8.3f} {b/1e6:10.1f} {S:10.4f} {exact:10.4f} {abs(S-exact):8.4f}")
floor = 1.0 / np.sqrt(N)
print(f"\nmax|Δ| = {worst:.4f}   (MC noise floor ~1/sqrt(N) = {floor:.4f})")
assert worst < 3 * floor, "replay deviates from analytic free diffusion beyond the noise floor"
print("OK — standalone produce→validate→replay round-trip matches analytic free diffusion.")
