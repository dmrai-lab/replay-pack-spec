# Replay Pack Codec Registry

**Status:** living document, versioned **independently** of the core specification.
**Registry version:** 0.1.1 (draft).
**License:** CC-BY-4.0 (registry text) / Apache-2.0 (reference code).

The core [`SPEC.md`](SPEC.md) §9 defines the two storage concepts and their rules: position
**representations** (§9.1) and their exactness *domain* (§9.2), per-channel **storage codecs**
(§9.3), and the interface every method MUST satisfy (§9.4) — declared in `compression.method`,
decodes to the §5 array contract, `identity` baseline, lossy self-certifies. This registry lists
the **concrete representations and codecs** and the **exact safetensors tensor keys** each stores.
It is a **storage / inference-optimization** layer: a method changes only *how* a channel is
stored (and, for a representation, how cheaply it replays), never *what it means* — the meaning of
every channel and of the replayed signal is fixed by the core (§5–6). The registry is versioned
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

### Position representations (dmipy-sim)

These reference-implementation **representations** of `positions` (SPEC §9.1) MAY be revised;
each changes only how positions are stored and how cheaply they replay, never what they decode
to (§5). Params live in `compression` (e.g. `K`).

| `method` | Class | Stored keys | Notes |
|---|---|---|---|
| `temporal_dct` | walker-preserving, lossless in-band | `dct_coeffs` | temporal band-limit; keeps all per-walker channels |
| `lowrank` | walker-preserving | `modes`, `mean`, `coeffs` | KL/SVD modes + exact per-walker coefficients; lossless at full rank |
| `gaussian` | **distributional** | `modes`, `mean`, `coeff_mean`, `coeff_cov` | resamples walkers ⇒ Gradient tier only (SPEC §9.4 walker-preserving rule) |
| `marginal` | **distributional** | `modes`, `mean`, `coeff_quantiles` | resamples walkers ⇒ Gradient tier only |

The per-walker physics channels each use a structure-matched **storage codec** (SPEC §9.3) —
they are **not** positions-like, so the representations above do not apply:

| Decoded channel | Codec | Stored keys | Notes |
|---|---|---|---|
| `compartment` (integer labels) | row run-length (RLE) | `comp_rle_vals`, `comp_rle_lens`, `comp_rle_counts` | lossless; long constant runs (impermeable pools) |
| `compartment` (fractional occupancy) | **quantized RLE** | same keys + meta `fractional:true, Q, scale` | **permeable** substrates: a walker crossing a membrane mid-save has a fractional time-in-compartment, so the channel is near-binary occupancy in `[0, scale]` (∈[0,1] for 2 compartments) rather than integer labels. Quantize to `Q` levels and RLE (same structure as `bound_fraction`); integer RLE would lossily int-cast the fractions. |
| `bound_fraction` | **quantized RLE** | `bfrac_rle_vals` (uint8, `Q≤256`), `bfrac_rle_lens` (uint16), `bfrac_rle_counts` | occupancy is ~binary with long dwell/free runs. Quantize `[0,1]→{0..Q-1}` (meta `Q`, default 256), RLE rows. Measured **~7×**, replay error ≪ MC floor. Lossy only to the quantization step. |
| `boundary_local_time` | **density-aware**: sparse CSR **or** dense int8 (meta `mode`) | sparse: `blt_counts`,`blt_cols`,`blt_qvals`; dense: `blt_dense_q` (int8) | not low-rank (idiosyncratic wall contacts). Density varies by substrate: isolated fibres ~15% nonzero → **sparse** (~3×); packed white matter ~55% → sparse would exceed raw float16, so **dense int8** (1 B/entry, ½ of raw f16, density-independent). Encoder picks the smaller; per-save values kept (any sequence gate, §6.6); quantization within the MC floor. |

Codec parameters (`Q`, `scale`, `nlevels`, `n_t`) live under `compression.channels[<channel>]` in
the metadata. The reference implementation also historically stored `boundary_local_time` raw
under the legacy key `dlog_boundary_unit` and `bound_fraction` raw under `bound_frac`; readers
SHOULD accept these raw aliases across the `1.x` line.

**Why these are separate codecs.** The positions codecs (low-rank/DCT) exploit temporal
smoothness and cross-walker correlation. `bound_fraction` (a near-binary step process) and
`boundary_local_time` (a sparse sum of discrete contact events) have neither property — low-rank
needs `K>32` and still misses the noise floor on `boundary_local_time` — so RLE and sparsity,
respectively, are the structure-matched choices.

## Distributional codecs are Gradient-only

A codec marked **distributional** (recognizable by the presence of `coeff_mean`/`coeff_cov`/
`coeff_quantiles`) resamples walkers and breaks per-walker channel alignment. Per SPEC §9.4, a
pack using one MUST NOT carry any per-walker channel and MUST declare only the Gradient tier.

## Replay in coefficient space (the representation mechanics)

This is the mechanics of the position *representation* defined in SPEC §9.1; its exactness
*domain* — raw exact to the `Δt` grid, a full basis coincident with raw, a truncated basis exact
only within its span — is SPEC §9.2 and is not restated here.

The `lowrank` and `temporal_dct` representations store positions as `r = M·c + μ` with a basis `M`
(KL/SVD modes or the DCT matrix) **shared across walkers** and per-walker coefficients `c`.
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

The **off-resonance/susceptibility** operation (§6.4) samples a field map at `r mod cell_size`
— a **nonlinear** function of position — and the **Bloch/MT** operation (§6.5) evolves
magnetization step by step; both REQUIRE the decoded positions. A replayer therefore takes the
coefficient-space fast path for the position-linear tiers (Gradient, Bulk-relaxation, Surface)
and falls back to decoded positions for the Field and Magnetization-transfer tiers. This is one reason to
prefer a linear-basis **representation**: it is not only smaller but also directly replayable — the IR property of SPEC §9.1.

The `identity` codec stores positions raw, so it has no coefficient space; the fast path applies
only to codecs that expose a shared linear basis (`modes`/`dct_coeffs`). Distributional codecs
already carry only the Gradient tier, and their coefficient statistics are consumed by §6.1
directly.

## Adding a codec

A new codec entry MUST document: `method` name, class (walker-preserving vs distributional),
its stored tensor keys, its parameters, and whether it is lossless or (if lossy) how its
`fidelity` is measured (SPEC §9.4). Add it here and bump the registry version; the core
`rpk_schema_version` does **not** change for a registry-only addition.
