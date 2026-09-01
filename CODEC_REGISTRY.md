# Replay Pack Codec Registry

**Status:** living document, versioned **independently** of the core specification.
**Registry version:** 0.4.0 (draft) — one position representation, `bridge_dst`; the four earlier position codecs (`temporal_dct`, `lowrank`, `gaussian`, `marginal`) are **retired and MUST be refused**, which is a breaking change for any pack written before it. 0.3.0 registered `susc_path_dct` for `susc_field_basis` (SPEC §6.4.1) alongside the raw-float16 baseline, the registry's first **capability-narrowing** codec (SPEC §9.4 rule 5).
**License:** CC-BY-4.0 (registry text) / Apache-2.0 (reference code).

The core [`SPEC.md`](SPEC.md) §9 defines the two storage concepts and their rules: position
**representations** (§9.1) and their exactness *domain* (§9.2), per-channel **storage codecs**
(§9.3), and the interface every method MUST satisfy (§9.4) — declared in `compression.method`,
decodes to the §5 array contract, `identity` baseline, lossy self-certifies. This registry lists
the **concrete representations and codecs** and the **exact safetensors tensor keys** each stores.
It is a **storage / inference-optimization** layer: a method changes only *how* a channel is
stored (and, for a representation, how cheaply it replays), never *what it means* — the meaning of
every channel and of the replayed signal is fixed by the core (§5–6). One method here does more
than trade accuracy: a truncated representation of a **non-linearly consumed** channel narrows
which requests the pack can serve at all, which the core requires be declared in `replay_envelope`
rather than inferred from `fidelity` (SPEC §9.4 rule 5). Such entries are marked
**capability-narrowing** and state the envelope key they oblige. The registry is versioned
independently so methods can be added, revised, or deprecated **without changing the frozen core
format**.

A conformant replayer implements the `identity` codec and MAY implement any subset of the
others; it MUST refuse (never guess) a `method` it does not recognize (SPEC §9, §13).

## Channel presence

For tiers (SPEC §7) and conformance (SPEC §13), a channel counts as **present** if the file
contains the channel's raw key **or** the codec's stored key-group listed below.

## Registered codecs

### `identity` (normative baseline — always lossless)

Every channel is stored raw under its own name.

| Decoded channel | Stored keys |
|---|---|
| `positions` | `positions` `(N_w, N_t, 3)` |
| `spin_weights` | `spin_weights` |
| `compartment` | `compartment` |
| `boundary_local_time` | `boundary_local_time` |
| `bound_fraction` | `bound_fraction` |
| `susc_field_{0,C,S}` | `susc_field_0`, `susc_field_C`, `susc_field_S` |
| `susc_field_basis` | `susc_field_basis` `(N_w, N_t, 13)` |

### Position representation

There is **one** C0 representation. It is not a menu: the codecs listed in earlier revisions of
this registry (`temporal_dct`, `lowrank`, `gaussian`, `marginal`) are **retired**, and a reader
MUST refuse them rather than decode them — see "Retired position codecs" below for why that is
a refusal and not a deprecation.

| `method` | Class | Stored keys | Notes |
|---|---|---|---|
| `bridge_dst` | walker-preserving, lossless at full rank | `pos_x`, `pos_y`, `pos_z`, each `(N_w, K+2)` | two exact endpoints, then `K` sine bands of the pinned residual |

Per axis the stored vector is `[r(0), r(T) - r(0), beta_1 .. beta_K]`, where `beta` are the
DST-I coefficients of the **Brownian bridge**

```
r(t) = r(0) + (t/T)[r(T) - r(0)] + u(t),      u(0) = u(T) = 0
```

**The leading two entries are not bands.** They are the endpoints, and `K` in `compression.K`
is the band count, so the stored width is `K + 2`. A consumer that reads the width and uses it
as `K` will compile a waveform projection of exactly the right shape to multiply and the wrong
thing to multiply by — it produces plausible numbers, not an error. Take `K` from the metadata.

Three properties motivate the split, none of them a compression claim:

- **The gradient moments are columns.** The phase is
  `gamma dt [ r(0).M0 + (r(T)-r(0)).M1 + sum_k beta_k Ghat_k ]` with `M0 = sum_n G_n` and
  `M1 = sum_n tau_n G_n`, so refocusing (`M0 = 0`) and velocity compensation (`M1 = 0`)
  annihilate the first two columns exactly. The design constraint and the storage layout are
  the same statement.
- **Both endpoints are exact at every truncation**, so a stored walk can be continued from
  where it ended rather than only replayed. A band codec's endpoint error is the size of its
  interior error.
- **The residual sits in its variance-optimal basis.** The discrete Brownian bridge has
  covariance `min(m,n) - mn/N`, whose inverse is the Dirichlet Laplacian, so its
  Karhunen-Loeve eigenvectors are exactly the DST-I vectors with eigenvalues
  `1/(4 sin^2(pi k / 2N))`.

Accuracy is indistinguishable from a cosine band truncation at equal budget, and that is a
theorem rather than a measurement: the difference operator maps one basis onto the other
exactly, `c_k(n) - c_k(n-1) = -2 sin(pi k / 2N) s_{k-1}(n)`, so a cosine expansion of a path
*is* a sine expansion of its increments and the two truncate to the same subspaces. The split
buys structure, not error.

Full rank is `K = n_t - 2`: two endpoints plus `n_t - 2` interior bands is exactly `n_t`
coefficients for `n_t` samples, so the representation is exactly rank-preserving.

### Retired position codecs

`temporal_dct`, `lowrank`, `gaussian` and `marginal` are no longer readable. This is a refusal
rather than a deprecation because the failure would otherwise be silent: `bridge_dst` puts the
two endpoints in the tensor positions a band codec fills with its two lowest bands, under the
same keys and with the same shape family. Decoding one as the other returns a trajectory that
looks like a trajectory. A pack declaring a retired method MUST be re-encoded from its master.

The per-walker physics channels each use a structure-matched **storage codec** (SPEC §9.3) —
they are **not** positions-like, so the representations above do not apply:

| Decoded channel | Codec | Stored keys | Notes |
|---|---|---|---|
| `compartment` (integer labels) | row run-length (RLE) | `comp_rle_vals`, `comp_rle_lens`, `comp_rle_counts` | lossless; long constant runs (impermeable pools) |
| `compartment` (fractional occupancy) | **quantized RLE** | same keys + meta `fractional:true, Q, scale` | **permeable** substrates: a walker crossing a membrane mid-save has a fractional time-in-compartment, so the channel is near-binary occupancy in `[0, scale]` (∈[0,1] for 2 compartments) rather than integer labels. Quantize to `Q` levels and RLE (same structure as `bound_fraction`); integer RLE would lossily int-cast the fractions. |
| `bound_fraction` | **quantized RLE** | `bfrac_rle_vals` (uint8, `Q≤256`), `bfrac_rle_lens` (uint16), `bfrac_rle_counts` | occupancy is ~binary with long dwell/free runs. Quantize `[0,1]→{0..Q-1}` (meta `Q`, default 256), RLE rows. Measured **~7×**, replay error ≪ MC floor. Lossy only to the quantization step. |
| `boundary_local_time` | **density-aware**: sparse CSR **or** dense int8 (meta `mode`) | sparse: `blt_counts`,`blt_cols`,`blt_qvals`; dense: `blt_dense_q` (int8) | not low-rank (idiosyncratic wall contacts). Density varies by substrate: isolated fibres ~15% nonzero → **sparse** (~3×); packed white matter ~55% → sparse would exceed raw float16, so **dense int8** (1 B/entry, ½ of raw f16, density-independent). Encoder picks the smaller; per-save values kept (any sequence gate, §6.6); quantization within the MC floor. |
| `susc_field_basis` | raw **float16** (lossy only to f16) | `susc_field_basis` `(N_w, N_t, 13)` | Per-walker Field basis (SPEC §6.4.1). Always valid, walker-preserving, and the only form that carries no envelope restriction. |
| `susc_field_basis` | **`susc_path_dct`** — temporal DCT-II, `K` modes, float16 | `susc_path_dct` `(N_w, n_ch, K)` | The size form; **capability-narrowing** (SPEC §9.4 rule 5) — see below. |

Codec parameters (`Q`, `scale`, `nlevels`, `n_t`) live under `compression.channels[<channel>]` in
the metadata. The reference implementation also historically stored `boundary_local_time` raw
under the legacy key `dlog_boundary_unit` and `bound_fraction` raw under `bound_frac`; readers
SHOULD accept these raw aliases across the `1.x` line.

**Why these are separate codecs.** The positions codecs (low-rank/DCT) exploit temporal
smoothness and cross-walker correlation. `bound_fraction` (a near-binary step process) and
`boundary_local_time` (a sparse sum of discrete contact events) have neither property — low-rank
needs `K>32` and still misses the noise floor on `boundary_local_time` — so RLE and sparsity,
respectively, are the structure-matched choices.

### `susc_path_dct` — the per-walker Field basis in `K` temporal modes

**Encode.** Sample the field-basis grids along each walker's **full-resolution** trajectory, then take an
orthonormal DCT-II along time of each component and keep the first `K` coefficients, stored float16 as
`susc_path_dct` `(N_w, n_ch, K)`. Sampling *before* compressing is the whole point: `φ^χ` is non-linear in
`r`, so no truncated representation of `r(t)` has a Parseval shortcut for it — but the field along a path
is far redder in time than `r(t)` is, so it compresses where the trajectory does not. This also
**decouples** the Field tier from the position representation: positions are free to be lossy, because
nothing re-samples a grid at decoded coordinates.

**`n_ch`.** `13` in general; `12` when `iso_P_zz` is dropped under the SPEC §5.3 trace identity and
reconstructed on decode; `6` when the source is purely isotropic (no `aniso_G` block). Stored order is the
§5.3 order with the omitted components removed, and `channels` lists the retained names.

**Decode to the §5.2 contract.** Inverse DCT-II back to `N_t`; re-insert
`iso_P_zz = 3·iso_local − iso_P_xx − iso_P_yy` where `iso_P_zz: "implied"`; insert an all-zero `aniso_G`
block where absent; and transpose to component-**last** `(N_w, N_t, 13)`.

**Parameters** under `compression.channels["susc_field_basis"]`: `K`, `n_t`, `n_ch`, `channels`,
`iso_P_zz` (`"implied"` | `"stored"`), `trace_residual`, `dtype`, `max_refocus_pulses`.

**Loss.** Lossless at `K = N_t` (to within `dtype`). Below that the loss is **not** a uniform accuracy
degradation but a **gate-bandwidth limit**. On the scalar-phase path this is exact: by orthonormality
`Σ_k s(t_k)·b(t_k) = Σ_m ŝ(m)·b̂(m)`, so truncation discards only the gate's out-of-band content, and a
narrowband gate is unaffected. On the **vector-Bloch** path, which decodes `b(t)` pointwise, the reference
implementation's measured requirement is `K ≥ 2·n_refocus`. A pack using this codec MUST therefore declare
`replay_envelope.acquisition.max_refocusing_pulses` (≈ `K/2`) alongside its `fidelity` (SPEC §9.4 rule 5),
and its `fidelity` battery MUST include a refocused (spin-echo / CPMG) acquisition — the small-difference
cancellation a truncated Field basis breaks first.

**Do not ship the §6.4 grid maps in the same pack as this codec over lossy positions.** The grid route
samples the field at codec-*decoded* coordinates and so is sound only while the position representation is
lossless; the path route samples the full-resolution walk at build time. Advertising both next to lossy
positions offers a replay route whose accuracy silently depends on a property the pack no longer has. Ship
the grid as a per-substrate companion artefact instead, and record which route the pack offers.

## Replay in coefficient space (the representation mechanics)

This is the mechanics of the position *representation* defined in SPEC §9.1; its exactness
*domain* — raw exact to the `Δt` grid, a full basis coincident with raw, a truncated basis exact
only within its span — is SPEC §9.2 and is not restated here.

`bridge_dst` stores positions as `r = M·c` with a basis `M` **shared across walkers** — the
constant, the ramp, and the sine functions — and per-walker coefficients `c`.
Because the gradient phase (SPEC §6.1) is **linear** in position, it factors through the basis:

```
φ_i = γΔt Σ_k G(t_k)·r_i(t_k) = c_i · (γΔt Mᵀ G_flat) + const,
```

so a replayer MAY contract the waveform against the basis **once** and evaluate every walker's
phase from its stored coefficients, never materializing the `(N_w, N_t, 3)` position array. The
per-compartment relaxation (§6.2) and surface-relaxivity (§6.3) log-weights are likewise
separable per-walker sums that need no position decode. This is a **pure performance
optimization**: the result MUST equal the decode-then-replay value (SPEC §6), and it changes
neither the stored keys nor conformance.

The **grid-map** off-resonance/susceptibility operation (§6.4) samples a field map at
`r mod cell_size` — a **nonlinear** function of position — and the **Bloch/MT** operation (§6.5)
evolves magnetization step by step; both REQUIRE the decoded positions. A replayer therefore
takes the coefficient-space fast path for the position-linear tiers (Gradient, Bulk-relaxation,
Surface) and falls back to decoded positions for the grid-map Field and Magnetization-transfer
tiers. (The **per-walker** Field basis of §6.4.1 is the exception: it is *pre-sampled* per walker,
so its phase is a per-walker separable sum — like the surface log-weight — needing **no** position
decode, and it composes with the coefficient-space fast path. Under `susc_path_dct` the sum is itself a
coefficient-space contraction, `Σ_k s(t_k)·b(t_k) = Σ_m ŝ(m)·b̂(m)`, so the scalar-phase route never
materializes `b(t)` either; the Bloch route does.) This is one reason to prefer a
linear-basis **representation**: it is not only smaller but also directly replayable — the IR
property of SPEC §9.1.

The `identity` codec stores positions raw, so it has no coefficient space. There is no coefficient space for a raw pack, so the fast path applies only to
`bridge_dst`.

## Adding a codec

A new codec entry MUST document: `method` name, class (walker-preserving vs distributional),
its stored tensor keys, its parameters, and whether it is lossless or (if lossy) how its
`fidelity` is measured (SPEC §9.4). Add it here and bump the registry version; the core
`rpk_schema_version` does **not** change for a registry-only addition.
