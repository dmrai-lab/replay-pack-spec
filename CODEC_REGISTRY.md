# Replay Pack Codec Registry

**Status:** living document, versioned **independently** of the core specification.
**Registry version:** 0.1.0 (draft).
**License:** CC-BY-4.0 (registry text) / Apache-2.0 (reference code).

The core [`SPEC.md`](SPEC.md) §9 defines only the codec *interface* — a codec is declared in
`compression.method`, decodes to the §5 array contract, `identity` is the lossless baseline,
and lossy codecs self-certify against the MC noise floor. This registry lists the **concrete
codecs** and the **exact safetensors tensor keys** each stores. It is intentionally separate
so new compression methods (several still being developed and described in the methodology
paper) can be added, revised, or deprecated **without changing the frozen core format**.

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

### Reference-implementation codecs (dmipy-sim)

These are **experimental** and MAY change while the methodology paper is in preparation. Each
applies to the `positions` channel unless noted; params live in `compression` (e.g. `K`).

| `method` | Class | Stored keys | Notes |
|---|---|---|---|
| `temporal_dct` | walker-preserving, lossless in-band | `dct_coeffs` | temporal band-limit; keeps all per-walker channels |
| `lowrank` | walker-preserving | `modes`, `mean`, `coeffs` | KL/SVD modes + exact per-walker coefficients; lossless at full rank |
| `gaussian` | **distributional** | `modes`, `mean`, `coeff_mean`, `coeff_cov` | resamples walkers ⇒ Gradient tier only (SPEC §9 walker-preserving rule) |
| `marginal` | **distributional** | `modes`, `mean`, `coeff_quantiles` | resamples walkers ⇒ Gradient tier only |

The dense per-walker physics channels each have their own structure-matched codec (they are
**not** positions-like, so the position codecs above do not apply):

| Decoded channel | Codec | Stored keys | Notes |
|---|---|---|---|
| `compartment` | row run-length (RLE) | `comp_rle_vals`, `comp_rle_lens`, `comp_rle_counts` | lossless; long constant runs (impermeable pools) |
| `bound_fraction` | **quantized RLE** | `bfrac_rle_vals` (uint8, `Q≤256`), `bfrac_rle_lens` (uint16), `bfrac_rle_counts` | occupancy is ~binary with long dwell/free runs. Quantize `[0,1]→{0..Q-1}` (meta `Q`, default 256), RLE rows. Measured **~7×**, replay error ≪ MC floor. Lossy only to the quantization step. |
| `boundary_local_time` | **sparse CSR** | `blt_counts` (int32, nonzeros/row), `blt_cols` (int16 col indices), `blt_qvals` (int16 quantized values) | channel is ~85% zeros and **not** low-rank (idiosyncratic wall contacts). Store nonzeros only, values quantized to a signed global `scale`/`nlevels` (meta). Per-save values kept (any sequence gate, §6.6). Measured **~3×**, exact at `nlevels=4096`. |

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
`coeff_quantiles`) resamples walkers and breaks per-walker channel alignment. Per SPEC §9, a
pack using one MUST NOT carry any per-walker channel and MUST declare only the Gradient tier.

## Adding a codec

A new codec entry MUST document: `method` name, class (walker-preserving vs distributional),
its stored tensor keys, its parameters, and whether it is lossless or (if lossy) how its
`fidelity` is measured (SPEC §9). Add it here and bump the registry version; the core
`rpk_schema_version` does **not** change for a registry-only addition.
