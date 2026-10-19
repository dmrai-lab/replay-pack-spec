"""Numerical check of the basis of SPEC §9.5 (NumPy only).

For free diffusion it walks a fine path, then checks, on the bridge residual u of a running
integral I(t) over [0, T]:

  1. the band  beta_k = int_0^T u(t) s_k(t) dt,  s_k = sqrt(2/T) sin(pi k t/T),  equals the
     cosine projection of the increments,  sqrt(2T)/(pi k) * int cos(pi k t/T) dI
     (the integration-by-parts identity of §9.5.3);
  2. the band variance is 2 D T^2 / (pi k)^2 = (T / (pi k)^2) QV with QV = 2 D T (§9.5.6);
  3. the series with K bands returns the endpoints exactly and the residual at an instant
     to the truncation variance (§9.5.4);
  4. the block quadrature of §9.5.5: the error of a band accumulated at (m, h), against a
     per-step accumulation, has absolute variance rel(m, w h) * (T / (pi k)^2) * QV, for
     m = 0 .. 3, with rel(m, x) = x^(2(m+1)) / (((m+1)!)^2 (2m+3) 4^(m+1));
  5. the registry's derivation (CODEC_REGISTRY.md): the band accumulated at m = 0 and
     h = T/N equals sqrt(h) sinc(pi k / 2N) times the DST-I band of the N+1 instants.

Run:  python continuous_basis_check.py          Apache-2.0.
"""
import math
import numpy as np

D, T, N, W = 2.0e-9, 0.05, 4096, 6000
dt = T / N
rng = np.random.default_rng(0)
dI = rng.normal(0.0, np.sqrt(2 * D * dt), (W, N))                 # increments of one axis
I = np.concatenate([np.zeros((W, 1)), np.cumsum(dI, axis=1)], axis=1)
t = np.arange(N + 1) * dt
u = I - (t / T) * I[:, -1:]                                       # bridge residual, §9.5.2
tm = (np.arange(N) + 0.5) * dt


def cont_band(k):
    """beta_k = int u s_k dt, by the trapezoid rule on the fine path"""
    return np.trapezoid(u * np.sqrt(2.0 / T) * np.sin(np.pi * k * t / T), t, axis=1)


def rel(m, x):
    return x ** (2 * (m + 1)) / (math.factorial(m + 1) ** 2 * (2 * m + 3) * 4 ** (m + 1))


# 1, 2: the identity and the free-bridge variance
for k in (1, 2, 5, 20):
    inc = dI @ np.cos(np.pi * k * tm / T)                          # int cos(pi k t/T) dI
    ratio = cont_band(k) / inc
    assert np.allclose(ratio, np.sqrt(2 * T) / (np.pi * k), rtol=1e-3), k   # §9.5.3 identity
    v_cont = 2 * D * T ** 2 / (np.pi * k) ** 2
    assert abs(cont_band(k).var() / v_cont - 1) < 0.05, k                  # §9.5.6
    print(f"k={k:3d}  band/(cosine projection) = {ratio.mean() / (np.sqrt(2*T)/(np.pi*k)):.5f}  "
          f"Var ratio {cont_band(k).var() / v_cont:.3f}")

# 3: the K-band series at an instant; exact at the endpoints, truncation variance elsewhere
K, t0 = 64, 0.31 * T
beta = np.stack([cont_band(k) for k in range(1, K + 1)], axis=1)
s_t = np.sqrt(2.0 / T) * np.sin(np.pi * np.arange(1, K + 1) * t0 / T)
rec = (t0 / T) * I[:, -1] + beta @ s_t                              # I(0) = 0
exact = I[:, int(round(t0 / dt))]
tail = sum(2 * D * T ** 2 / (np.pi * k) ** 2 * (2.0 / T) * np.sin(np.pi * k * t0 / T) ** 2
           for k in range(K + 1, 20000))
print(f"I(t) at t = {t0 * 1e3:.2f} ms: error variance {np.var(rec - exact):.3e}  (predicted truncation {tail:.3e})")
assert abs(np.var(rec - exact) / tail - 1) < 0.1

# 4: the block quadrature, moments m = 0 .. 3 within blocks of n_b fine steps, against the fine accumulation
n_b = 128
h = n_b * dt
nblk = N // n_b
blk = np.arange(N) // n_b
tau = tm - (blk + 0.5) * h                                          # time from the block centre
QV = (dI ** 2).sum(axis=1).mean()                                   # = 2 D T on a free walk
for m, k in ((0, 2), (1, 12), (2, 20), (3, 28)):
    w = np.pi * k / T
    ref = dI @ np.cos(w * tm)                                       # per-step accumulation
    acc = np.zeros(W)
    for j in range(m + 1):                                          # Taylor of cos about the block centre
        Mj = np.add.reduceat(dI * tau ** j, np.arange(0, N, n_b), axis=1)      # int tau^j dI per block
        cj = (w ** j / math.factorial(j)) * np.cos(w * (np.arange(nblk) + 0.5) * h + j * np.pi / 2)
        acc += Mj @ cj
    err_var = np.var(acc - ref)
    # the band error is (sqrt(2T)/(pi k)) times this; its variance in units of (T/(pi k)^2) QV is err_var / (QV/2)
    got = err_var / (QV / 2.0)
    want = rel(m, w * h)
    print(f"m={m} k={k:2d}  w h = {w * h:.3f}  quadrature variance / ((T/(pi k)^2) QV) = {got:.3e}   rel = {want:.3e}")
    assert abs(got / want - 1) < 0.15, (m, got, want)

# 5: the registry's derivation from the DST-I of the N_s + 1 instants
Ns = 64
Is = I[:, :: N // Ns]                                               # the path at Ns + 1 instants
ts = np.arange(Ns + 1) * (T / Ns)
us = Is - (ts / T) * Is[:, -1:]
mi = np.arange(1, Ns)
dIs = np.diff(Is, axis=1)                                           # increments between the instants
hs = T / Ns
for k in (1, 3, 10):
    b_dst = us[:, 1:-1] @ (np.sqrt(2.0 / Ns) * np.sin(np.pi * k * mi / Ns))
    beta0 = np.sqrt(2 * T) / (np.pi * k) * (dIs @ np.cos(np.pi * k * (np.arange(Ns) + 0.5) / Ns))
    x = np.pi * k / (2 * Ns)
    assert np.allclose(beta0, np.sqrt(hs) * (np.sin(x) / x) * b_dst, rtol=1e-9, atol=1e-18), k
print("derivation: beta(m=0, h) = sqrt(h) sinc(pi k / 2N) * (DST-I band)   exact")
print("OK")
