"""Numerical check of the two band bases of SPEC §9.5 (NumPy only).

For free diffusion it walks a fine path, then checks, on the same bridge residual u of
a running integral I(t) over [0, T]:

  1. the `continuous` band  beta_k = int_0^T u(t) s_k(t) dt, s_k = sqrt(2/T) sin(pi k t/T),
     equals the cosine projection of the increments,  sqrt(2T)/(pi k) * int cos(pi k t/T) dI
     (the integration-by-parts identity of §9.5.3);
  2. the `dst1` band variance is D dt / (2 sin^2(pi k / 2N)) and the `continuous` one is
     2 D T^2 / (pi k)^2, the second being the limit of the first after the change of
     normalisation beta_cont = sqrt(dt) * beta_dst1 (§9.5.6);
  3. the series with K bands returns the endpoints exactly and the residual at an instant
     to the truncation variance (§9.5.4).

Run:  python continuous_basis_check.py          Apache-2.0.
"""
import numpy as np

D, T, N, W = 2.0e-9, 0.05, 512, 20000
dt = T / N
rng = np.random.default_rng(0)
dI = rng.normal(0.0, np.sqrt(2 * D * dt), (W, N))                 # increments of one axis
I = np.concatenate([np.zeros((W, 1)), np.cumsum(dI, axis=1)], axis=1)
t = np.arange(N + 1) * dt
u = I - (t / T) * I[:, -1:]                                       # bridge residual, §9.5.1

# dst1: orthonormal DST-I of u at the n_t - 2 = N - 1 interior saves (n_t = N + 1)
m = np.arange(1, N)
def dst1_band(k):
    return u[:, 1:-1] @ (np.sqrt(2.0 / N) * np.sin(np.pi * k * m / N))

# continuous: beta_k = int u s_k dt, by the trapezoid rule on the fine path
def cont_band(k):
    return np.trapezoid(u * np.sqrt(2.0 / T) * np.sin(np.pi * k * t / T), t, axis=1)

tm = (np.arange(N) + 0.5) * dt
for k in (1, 2, 5, 20):
    inc = dI @ np.cos(np.pi * k * tm / T)                          # int cos(pi k t/T) dI
    ratio = cont_band(k) / inc
    assert np.allclose(ratio, np.sqrt(2 * T) / (np.pi * k), rtol=1e-3), k   # §9.5.3 identity
    v_dst1 = D * dt / (2 * np.sin(np.pi * k / (2 * N)) ** 2)
    v_cont = 2 * D * T ** 2 / (np.pi * k) ** 2
    assert abs(dst1_band(k).var() / v_dst1 - 1) < 0.05, k                  # §9.5.6, dst1
    assert abs(cont_band(k).var() / v_cont - 1) < 0.05, k                  # §9.5.6, continuous
    assert abs(v_dst1 * dt / v_cont - 1) < 5e-4 * k ** 2, k                # the limit, sqrt(dt) scale
    print(f"k={k:3d}  cont/(cosine projection) = {ratio.mean() / (np.sqrt(2*T)/(np.pi*k)):.5f}  "
          f"Var ratio  dst1 {dst1_band(k).var() / v_dst1:.3f}  continuous {cont_band(k).var() / v_cont:.3f}")

# §9.5.4: the K-band series at an instant; exact at the endpoints, truncation variance elsewhere
K, t0 = 64, 0.31 * T
beta = np.stack([cont_band(k) for k in range(1, K + 1)], axis=1)
s_t = np.sqrt(2.0 / T) * np.sin(np.pi * np.arange(1, K + 1) * t0 / T)
rec = (t0 / T) * I[:, -1] + beta @ s_t                              # I(0) = 0
exact = I[:, int(round(t0 / dt))]
tail = sum(2 * D * T ** 2 / (np.pi * k) ** 2 * (2.0 / T) * np.sin(np.pi * k * t0 / T) ** 2
           for k in range(K + 1, 20000))
print(f"I(t) at t = {t0 * 1e3:.2f} ms: error variance {np.var(rec - exact):.3e}  (predicted truncation {tail:.3e})")
assert abs(np.var(rec - exact) / tail - 1) < 0.1
print("OK")
