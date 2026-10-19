"""End-to-end standalone demo (no dmipy-sim, no GPU):

  produce  ->  a Tier-0 pack of bands (write_minimal_rpk)
  validate ->  conformance_check
  replay   ->  reference_replayer, for several b-values, vs the exact free-diffusion
               signal exp(-bD). Agreement to the MC noise floor proves the format +
               reference replay operation are self-consistent and engine-agnostic.

Run:  python demo_free_diffusion.py
Apache-2.0.
"""
import numpy as np
import write_minimal_rpk as wr
import reference_replayer as rr
import conformance_check as cc

D, T, K, n_steps, N = 2.0e-9, 0.05, 128, 1000, 40000
delta, Delta, dt_wf = 0.010, 0.040, T / 500
gs = (0.02, 0.05, 0.10, 0.20)
battery = [(rr.pgse_waveform(dt_wf, delta, Delta, g, [1, 0, 0])[0], dt_wf) for g in gs]
vec, qv = wr.accumulate_bands(N, T, n_steps, D, K, seed=0)
wr.write_pack("free_diffusion.rpk", vec, T, T / n_steps, D, K, qv, battery, seed=0,
              pack_id="examples/free-diffusion")

print("== conformance ==")
assert cc.main("free_diffusion.rpk") == 0

print("\n== replay vs analytic exp(-bD) ==")
arrays, meta = rr.load_rpk("free_diffusion.rpk")
print(f"{'g[T/m]':>8} {'b[s/mm2]':>10} {'replay':>10} {'exp(-bD)':>10} {'|Δ|':>8}")
worst = 0.0
for g in gs:
    G, b = rr.pgse_waveform(dt_wf, delta, Delta, g, direction=[1, 0, 0])
    S = rr.replay_gradient(arrays, meta, G, dt_wf).real
    exact = np.exp(-b * D)
    worst = max(worst, abs(S - exact))
    print(f"{g:8.3f} {b/1e6:10.1f} {S:10.4f} {exact:10.4f} {abs(S-exact):8.4f}")
floor = 1.0 / np.sqrt(N)
print(f"\nmax|Δ| = {worst:.4f}   (MC noise floor ~1/sqrt(N) = {floor:.4f})")
assert worst < 3 * floor, "replay deviates from analytic free diffusion beyond the noise floor"
print("OK — standalone produce→validate→replay round-trip matches analytic free diffusion.")
