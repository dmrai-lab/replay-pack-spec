# Replay Pack Codec Registry

**Status:** living document, versioned **independently** of the core specification.
**Registry version:** 0.7.0 (draft). A pack is bands on **one** basis, the continuous sines of `[0, T]` (SPEC §9.5.3), plus exact occupancy in time (SPEC §9.5.7). The registry defines the stored tensor keys of each channel, the integer band container (per-band scales with a leading **block** axis and an optional per-walker `band_block`, so shards of a distributed fill merge without re-quantisation), the `static` and event-list forms of C1, and the forms a reader refuses.
**License:** CC-BY-4.0 (registry text) / Apache-2.0 (reference code).

The core [`RPK.md`](RPK.md) §9 defines the two storage concepts and their rules: the position
**representation** (§9.1) and its exactness *domain* (§9.2), per-channel **storage codecs**
(§9.3), and the interface every method MUST satisfy (§9.4) — declared in `compression.method`,
decodes to the §5 channel contract, every pack self-certifies. This registry lists
the **concrete representation and codecs** and the **exact safetensors tensor keys** each stores.
It is a **storage / inference-optimization** layer: a method changes only *how* a channel is
stored (and, for a representation, how cheaply it replays), never *what it means* — the meaning of
every channel and of the replayed signal is fixed by the core (§5–6). One method here does more
than trade accuracy: a truncated representation of a **non-linearly consumed** channel narrows
which requests the pack can serve at all, which the core requires be declared in `replay_envelope`
rather than inferred from `fidelity` (SPEC §9.4 rule 5). Such entries are marked
**capability-narrowing** and state the envelope key they oblige. The registry is versioned
independently so methods can be added or revised **without changing the frozen core
format**.

A conformant replayer implements `bridge_dst` and MAY implement any subset of the
others; it MUST refuse (never guess) a `method` it does not recognize (SPEC §9, §13).

## Channel presence

For tiers (SPEC §7) and conformance (SPEC §13), a channel counts as **present** if the file
contains the codec's stored key-group listed below. A channel has no raw key: there is no
tensor of samples on a grid.

## Registered codecs

### Stored keys by channel

| Decoded channel | Stored keys |
|---|---|
| `positions` | `pos_x`, `pos_y`, `pos_z` `(N_w, K+2)` (`bridge_dst`, float container) or the integer band container below |
| `spin_weights` | `spin_weights` `(N_w,)` |
| `compartment` | per column: `comp_static` `(N_w,)`, or `comp_events_offsets`, `comp_events_t`, `comp_events_label`; the `bound` column the same under `bound_` |
| `boundary_local_time` | `blt_bridge_dst`, `blt_start`, `blt_endpoint`, or the integer band container |
| `susc_field_basis` | `susc_path_dct` `(N_w, n_ch, K)` and `susc_path_scale` |

### Position representation

There is **one** C0 representation.

| `method` | Class | Stored keys | Notes |
|---|---|---|---|
| `bridge_dst` | walker-preserving | **float container**: `pos_x`, `pos_y`, `pos_z`, each `(N_w, K+2)` float32 | two exact endpoints, then `K` sine bands of the pinned residual |
| `bridge_dst` (integer band container) | walker-preserving; the bands quantised per band | per axis `pos_x_ends` `(N_w, 2)` float32 and, per band range `i` of `bits` bits, `pos_x_b<i>` `(N_w, k1-k0)` int16 or int8 when `bits` is 16 or 8, else the **bit-packed** row `pos_x_b<i>u<bits>` `(N_w, ceil((k1-k0) * bits / 8))` uint8; `pos_band_scale` `(n_blocks, n_axes, K)` float32; optional `band_block` `(N_w,)` uint16 / int32 | `compression.container` lists the ranges `[{bands: [k0, k1], bits}, ...]` covering `0 .. K-1`, `bits` an integer in `3 .. 16`. Band `k` of walker `w` on axis `a` is `q * pos_band_scale[block(w), a, k]`, with `block(w) = band_block[w]` (0 when absent). The scale of a band is `max_w |beta_k| / (2^(bits-1) - 1)` over the block's walkers. A bit-packed range holds its bands side by side in one row per walker, each as the offset-binary `q + 2^(bits-1)` in `bits` bits, most significant bit first, the last byte padded with zero bits; the depth is in the tensor's name so that a reader without the metadata unpacks it. A reader MUST refuse a pack whose scale table has more than one block and no `band_block`. |

The method name `bridge_dst` is the registered identifier of this representation. Per axis the stored
vector is `[r(0), r(T) - r(0), beta_1 .. beta_K]`, where `beta` are the sine-basis coefficients of the
**Brownian bridge** residual (SPEC §9.5.2, §9.5.3)

```
r(t) = r(0) + (t/T)[r(T) - r(0)] + u(t),      u(0) = u(T) = 0,
s_k(t) = sqrt(2/T) sin(pi k t / T),           beta_k = int_0^T u(t) s_k(t) dt .
```

**The leading two entries are not bands.** They are the endpoints, and `K` in `compression.K`
is the band count, so the stored width is `K + 2`. A consumer that reads the width and uses it
as `K` will compile a waveform projection of exactly the right shape to multiply and the wrong
thing to multiply by — it produces plausible numbers, not an error. Take `K` from the metadata.

**`K` is a count; the band is a frequency.** `K` sine bands over `T_max` resolve the residual up to
`K / (2 · T_max)` hertz (SPEC §9.4 rule 6), recorded as `compression.temporal_bandwidth_hz`. The same
`K` on a walk twice as long is half the band; a bank compares walks and requests packs by that frequency,
and a prefix of the walk (SPEC §3) keeps it by re-encoding at `K' = ⌈K · TE' / T_max⌉` or more --
the certificate, not the ratio, fixes `K'`.

Three properties motivate the split, none of them a compression claim:

- **The gradient moments are columns.** The phase is
  `gamma [ r(0).M0 + (r(T)-r(0)).M1 + sum_k beta_k Ghat_k ]` with `M0 = int G dt` and
  `M1 = int (t/T) G dt`, so refocusing (`M0 = 0`) and velocity compensation (`M1 = 0`)
  annihilate the first two columns exactly. The design constraint and the storage layout are
  the same statement.
- **Both endpoints are exact at every truncation**, so a stored walk can be continued from
  where it ended rather than only replayed. A band codec's endpoint error is the size of its
  interior error.
- **The residual sits in its variance-optimal basis.** The Brownian bridge has covariance
  `2 D (min(s,t) - s t / T)`, whose Karhunen-Loeve eigenfunctions are exactly the `s_k`, with
  variances `2 D T^2 / (pi k)^2` (SPEC §9.5.6).

A cosine expansion of a path *is* a sine expansion of its increments: `beta_k = (sqrt(2T) / (pi k)) int cos(pi k t / T) dI`
(SPEC §9.5.3), which is what the walk accumulates (SPEC §9.5.5). The split buys structure, not error.

#### Derivation from the DST-I of saves

*Informative; the derivation of the accumulation identity, not a readable kind.* A walk saved at `n_t` instants
`t_m = m h`, `h = T/N`, `N = n_t - 1`, has the orthonormal DST-I bands of its bridge residual,
`b_k = sum_{m=1}^{N-1} sqrt(2/N) sin(pi k m / N) u_m`. By summation by parts the DST-I band of the saves is
`b_k = sqrt(2/N) / (2 sin(pi k / 2N)) * sum_j cos(pi k (j + 1/2) / N) dI_j`, and the sine-basis band of the
same increments accumulated at moment order `m = 0` and block `h` (cosines at the block midpoint) is

```
beta_k = sqrt(h) * sinc(pi k / (2N)) * b_k ,          sinc(x) = sin(x) / x .
```

The free-bridge variance of the DST-I band, `D h / (2 sin^2(pi k / 2N))`, then maps to `2 D T^2 / (pi k)^2`
at `m = 0` as `N -> infinity`. `examples/continuous_basis_check.py` verifies the identity numerically.

### Refused forms

A reader MUST refuse, by name and with the instruction to re-walk (SPEC §9.5.1), a pack carrying any of:

- the metadata keys `n_t` or `dt_traj` (in `walk_params`, its window table, or anywhere under `compression`) or `compression.basis`;
- the codec `identity`, and the tensors that go with it: `positions`, `compartment`, `boundary_local_time`, `susc_field_basis`, `dlog_boundary_unit`;
- the occupancy and contact forms on a grid: `comp_rle_vals`/`comp_rle_lens`/`comp_rle_counts`, `bound_rle_*`, `blt_counts`/`blt_cols`/`blt_qvals`, `blt_dense_q`, and any `<name>_rle_*`;
- the field grid as a field route: the tensors `susc_field_0`, `susc_field_C`, `susc_field_S` and the metadata key `walk_params.cell_size`;
- the retired position codecs `temporal_dct`, `lowrank`, `gaussian` and `marginal` (`dct_coeffs`, `modes`), the retired `blt_dct_coeffs`, and the retired `bfrac_rle_*` keys of MT as a channel of its own (it is the `bound` column of `compartment`, SPEC §5.2, §6.5).

This is a refusal rather than a conversion because the failure would otherwise be silent: the band tensors
of a pack on a save grid have the same keys and the same shape family as those of this revision, and a retired
position codec puts its two lowest bands where `bridge_dst` puts its endpoints. Decoding one as the other
returns a trajectory that looks like a trajectory. Such a pack is re-walked.

## Occupancy and the physics channels

The per-walker physics channels each use a structure-matched **storage form** (SPEC §9.3) —
they are **not** positions-like, so the representation above does not apply:

| Decoded channel | Codec | Stored keys | Notes |
|---|---|---|---|
| `compartment` → a column, one label per walker | **static** | `comp_static` `(N_w,)` int8 (`bound_static` for the `bound` column) | The label of the walker over the whole walk: nothing crossed (an impermeable substrate, no binding). Descriptor `kind: "static"`. A producer MUST write this form when it applies, and MUST NOT when any walker changes label. |
| `compartment` → a column, labels changing in time | **events** | `<c>_events_offsets` `(N_w + 1,)` int64, `<c>_events_t` `(N_ev,)` float64, `<c>_events_label` `(N_ev,)` int8 (int16 where a pool id exceeds 127), with `<c>` = `comp` or `bound` | Descriptor `kind: "events"`. Walker `i`'s events are the entries `offsets[i] .. offsets[i+1]-1`. `t` is in **seconds from the walk's start**, `label` is the label **entered** at `t`: a pool id for `comp`, `0` free / `1` bound for `bound`. Every walker's first event is at `t = 0` and gives its initial label; one walker's times strictly increase within `[0, T_max]`; no event repeats the label before it. The walker holds `label_j` on `[t_j, t_{j+1})` and its last label to `T_max`. A crossing, a binding and an unbinding are one event each. In a window table each window's list holds its own events, times still from the walk's start, opening with an event at the window's start. A producer that caps the events per walker MUST refuse on overflow rather than drop an event. The `bound` run intervals are the dwell times. |
| `boundary_local_time` | **`bridge_dst`** — two exact endpoints + sine bands of the pinned cumulative `B(t)` (meta `mode: "bridge_dst"`); the bands in the **float container** `blt_bridge_dst` `(N_w, K)` (float16 or float32, meta `dtype`) **or the integer band container** `blt_b<i>` per range with `blt_band_scale` `(n_blocks, K)` and the shared `band_block` (meta `container`, `dtype: "bands"`), the same rule as positions | `blt_bridge_dst` or `blt_b<i>` + `blt_band_scale`; `blt_start` (f32), `blt_endpoint` (f32) | The cumulative `B(t)` is smooth and nondecreasing, so it takes the same split as C0: hold both endpoints and expand only the doubly-pinned residual. A truncated *cosine* residual does not vanish at `t = T`, so replay would miss the stored total by percent at **every** `K` and segment chaining would drift; the sine form is identically zero there. Endpoints stay float32 (they are read by subtraction, §6.3); bands may be float16. The walk records `walk_params.qv.local_time` for the certificate. |
| `susc_field_basis` | **`susc_path_dct`** — `K` temporal modes of the field `b(t)`, integer container (meta `bits`, default **8**, one depth over the modes; or meta `container` listing ranges of modes `[{modes: [k0, k1], bits}, ...]` at their own depths, the same rule as the position bands) | `susc_path_dct` `(N_w, n_ch, K)` int8/int16, or per range `i` `susc_path_dct_b<i>` (int8/int16) / `susc_path_dct_b<i>u<bits>` (bit-packed uint8, `(N_w, n_ch, ceil((k1-k0) * bits / 8))`); `susc_path_scale` `(n_ch, K)` float32 | The only form of this channel; **capability-narrowing** (SPEC §9.4 rule 5) — see below. The modes are `c_k = int_0^T e_k(t) b(t) dt`, `e_0 = T^(-1/2)`, `e_k = sqrt(2/T) cos(pi k t / T)`, the sine bands of `F = int b dt` (SPEC §9.5.3). Above its mean the sampled field is nearly white (band variance falls ~1 decade over 63 bands, slope ≈ −0.6, against ≈ −2 for a trajectory), so the transform buys no energy compaction — it buys a basis in which the **gate** is sparse, and the distortion weight is `ε̂_k²`. Bits, not `K`, is therefore the graded knob: int8 at full `K` is half of float16 and beats `K/2` at float16 at equal bytes **while keeping** the refocusing capability truncation forfeits. `bits: null` stores raw floats. `susc_path_scale` is not walker-indexed and so is read whole. The walk records `walk_params.qv.field` and `walk_params.field_read_interval` for the certificate. |

Codec parameters (`Q`, `scale`, `nlevels`) live under `compression.channels[<channel>]` in
the metadata, and no descriptor carries `n_t`. Per-column descriptors (`name`, `kind`) live under
`compression.channels.compartment.columns`. There is **no** separate `bound_fraction` channel: the MT bound pool is
a `bound` column of `compartment` (SPEC §5.2, §6.5).

**Why these are separate codecs.** The positions codec exploits temporal smoothness. Occupancy columns
(step processes in time) and `boundary_local_time` (a sparse sum of discrete contact events) have a different
structure, so an event list and a bridge of the cumulative are the structure-matched choices.

**Not registered, and why.** Reducing the `n_ch` components to a few combinations of them was
measured and **rejected**: they are strongly correlated (rank 4 holds 97% of the variance) but a
rank-6 KLT at matched bytes is two to three orders worse than band truncation. The distortion
weight `ε̂_k²` lives on `k` and is blind to `c`, so a discarded *combination* removes signal the
gate passes at full weight, whereas a discarded *mode* is one the gate barely reaches. Bits belong
in frequency, never across channels.

### `susc_path_dct` — the per-walker Field basis in `K` temporal modes

**Encode.** The walk accumulates `c_k = int e_k b dt` for each component of the field basis at the walker's own
position, at its reads (SPEC §9.5.5), and the pack stores the first `K` as integers with a per-(channel, band) scale
as `susc_path_dct` `(N_w, n_ch, K)`. Reading the field at the walk's position, never along the stored bands, is
the whole point: `φ^χ` is non-linear in `r`, so no truncated representation of `r(t)` has a Parseval shortcut for
it — but the field along a path is far redder in time than `r(t)` is, so it compresses where the trajectory does
not. This also **decouples** the Field tier from the position representation: positions are band-limited, and
nothing re-samples a grid at decoded coordinates.

**`n_ch`.** `13` in general; `12` when `iso_P_zz` is dropped under the SPEC §5.3 trace identity and
reconstructed on decode; `6` when the source is purely isotropic (no `aniso_G` block). Stored order is the
§5.3 order with the omitted components removed, and `channels` lists the retained names.

**Decode to the §5.2 contract.** Evaluate `b(t) = sum_k c_k e_k(t)` at an instant; re-insert
`iso_P_zz = 3·iso_local − iso_P_xx − iso_P_yy` where `iso_P_zz: "implied"`; insert an all-zero `aniso_G`
block where absent; and order the components as §5.3.

**Parameters** under `compression.channels["susc_field_basis"]`: `K`, `n_ch`, `channels`,
`iso_P_zz` (`"implied"` | `"stored"`), `trace_residual`, `dtype`, `bits` or `container` (the mode ranges and their depths), `max_refocus_pulses`.

**Loss.** There is no full rank: the series is a band-limited `b(t)`. The loss is **not** a uniform accuracy
degradation but a **gate-bandwidth limit**. On the scalar-phase path this is exact: by orthonormality
`∫ s(t)·b(t) dt = Σ_k ŝ_k c_k`, so truncation discards only the gate's out-of-band content, and a
narrowband gate is unaffected. On the **vector-Bloch** path, which decodes `b(t)` pointwise, the reference
implementation's measured requirement is `K ≥ 2·n_refocus`. A pack using this codec MUST therefore declare
`replay_envelope.acquisition.max_refocusing_pulses` (≈ `K/2`) alongside its `fidelity` (SPEC §9.4 rule 5),
and its `fidelity` battery MUST include a refocused (spin-echo / CPMG) acquisition — the small-difference
cancellation a truncated Field basis breaks first.

**No field grid is a field route.** The field is accumulated at the walk's own position and shipped as `susc_path_dct` only. A field grid MAY accompany a pack as a companion artefact or an `x_` channel, as provenance or for visualisation; no replay reads it.

## Replay in coefficient space (the representation mechanics)

This is the mechanics of the position *representation* defined in SPEC §9.1; its exactness
*domain* — exact endpoints, exact within the span of the `K` modes, bounded and certified outside it — is SPEC §9.2
and is not restated here.

`bridge_dst` stores positions as `r = M·c` with a basis `M` **shared across walkers** — the
constant, the ramp, and the sine functions — and per-walker coefficients `c`.
Because the gradient phase (SPEC §6.1) is **linear** in position, it factors through the basis:

```
φ_i = γ ∫ G(t)·r_i(t) dt = c_i · (γ Mᵀ G) + const,
```

so a replayer contracts the waveform against the basis **once**, in closed form over the waveform's pieces
(SPEC §9.5.3), and evaluates every walker's phase from its stored coefficients, never materializing the
trajectory. The per-compartment relaxation (§6.2) and surface-relaxivity (§6.3) log-weights are likewise
separable per-walker quantities (a sum over the runs of the occupancy; a difference of two reads of `B(t)`)
that need no position decode.

The **Bloch/MT** operation (§6.5) evolves magnetization over time and REQUIRES the position evaluated at
instants (SPEC §9.5.4). A replayer therefore takes the coefficient-space path for the position-linear tiers
(Gradient, Bulk-relaxation, Surface) and evaluates the function at instants for the Magnetization-transfer tier.
The **per-walker** Field basis of §6.4.1 is *accumulated* per walker, so its phase is a per-walker separable sum --
like the surface log-weight -- needing **no** position decode, and it composes with the coefficient-space path.
Under `susc_path_dct` the sum is itself a coefficient-space contraction, `∫ s(t)·b(t) dt = Σ_m ŝ_m c_m`, so the
scalar-phase route never materializes `b(t)` either; the Bloch route does. This is one reason for the
representation: it is not only compact but directly replayable (SPEC §9.1).

## Adding a codec

A new codec entry MUST document: `method` name, class (walker-preserving vs distributional),
its stored tensor keys, its parameters, and how its `fidelity` is measured (SPEC §9.4). Add it here and bump the
registry version; the core `rpk_schema_version` does **not** change for a registry-only addition.
