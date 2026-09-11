# The Acquisition Specification (`ScannerSequence`)

**Status:** draft, versioned independently of the core `.rpk` specification.
**Version:** 0.2.0 (draft for comment). Nothing is numbered `1.0` before publication.
**License:** CC-BY-4.0 (text) / Apache-2.0 (reference code).

*Companion to the Replay Pack Specification ([`RPK.md`](RPK.md)), the Replay Phantom Specification
([`RPH.md`](RPH.md)) and the Substrate Specification ([`SUBSTRATE.md`](SUBSTRATE.md)). `RPK.md` §6.6 says
that the acquisition is a replay knob and that a pack stores no sequence; this document says what that knob
**is**.*

An **acquisition** is what the scanner does from the excitation to the readout: the gradient it plays, the
RF pulses it fires and when, the dead times it needs around them, the emergent crushing it applies, and the
sample at which it reads. It is one object. Every consumer -- the fused Monte-Carlo walk, the replay of a
pack (scalar-phase tiers), the vector-Bloch replay, an analytical signal model, a Pulseq exporter -- reads
that one object and nothing beside it, and everything a consumer might want to know that is not stored
(the effective gradient, the coherence state, where the echoes form, the b-value) is **derived** from it by
the rules of §4.

---

## 1. Scope and purpose

This specification defines:

- the **data model** of an acquisition (§3): the stored fields, their units and shapes, and the serialised
  record;
- the **derivation rules** (§4): how the effective gradient, the coherence state, the echo times and the
  b-value follow from the stored gradient and the RF schedule;
- the **invariants** a conformant acquisition satisfies (§5) and that a consumer MAY assume;
- the **families** (§6): how the standard encodings -- PGSE, PGSTE, OGSE, CPMG, gradient echo, spherical and
  planar tensor encoding -- are assembled from a small set of gradient shapes around a small set of RF
  schedules, so that a family is a thin name for one mechanics rather than its own code path;
- **what each consumer reads** (§7) and the **scanner-limit catalogue** an acquisition is built against (§8);
- the **Pulseq interchange** (§9) and **conformance** (§10).

It does not standardise how a walk is integrated, how a pack is stored (`RPK.md`), or how a design
optimiser searches over acquisitions: an optimiser is a producer of this object under the limits of §8.

### 1.1 Why a standard

Before this document the same acquisition existed in several shapes at once: a gradient array with an
implied ideal spin echo, a gradient array with an explicit list of finite pulses, an analytical parameter
set (b, δ, Δ, TE), and a set of constructor flags (`chi_perp`, `TM`, `stimulated_echo`, `echo_indices`)
restating what the pulses already said. Consumers took the gradient and its RF schedule as separate
arguments and could be handed two that disagreed; some stored the *effective* gradient (the 180 already
folded in), others the *physical* one, and a consumer that folded twice, or never, produced a plausible
number. A timing budget (pulse durations, readout dead time) lived in two other packages and in neither
container. Scanner limits were catalogued five times with three values for the same slew rate.

The fix is one object with one rule per fact: the gradient stored is the one the scanner plays; the RF
schedule is the one list of pulses; every coherence fact is derived from those two; a timing budget is a
field; and no consumer takes the RF, the echo or the refocus time as a separate argument.

---

## 2. Terminology

- **Acquisition** (`ScannerSequence`) -- the object of §3: one RF schedule, one grid, one echo time.
- **Physical gradient** `G` -- what the gradient amplifiers play, per measurement, in T/m. Stored.
- **Effective gradient** `G_eff` -- the gradient the transverse phase integrates, `G_eff(t) = G(t)·s(t)`.
  Derived (§4.2).
- **Gate** `s(t)` -- the transverse-phase sign of `RPK.md` §6.6: `+1`, flipping at every refocusing pulse
  and at a stimulated echo's recall. Derived from the schedule (§4.2).
- **RF event** -- one pulse: its instant, flip angle, role, B1 axis, duration, carrier offset, and
  optionally its played envelope (§3.2).
- **RF schedule** -- the time-ordered tuple of RF events of an acquisition.
- **Role** -- what a pulse does to the coherence: `excite` (z → transverse), `store` (transverse → z),
  `recall` (z → transverse again), `refocus` (invert the transverse phase). A pulse's **label** states its
  role; an unlabelled pulse is read by its flip (§4.3).
- **Coherence mask** `χ(t)` -- the transverse fraction of the encoded magnetisation over the grid: `1`
  while transverse, `0` while stored along z, fractional across a finite pulse (§4.3).
- **Readout** -- the grid samples at which the signal is read; the last of them is **the echo sample**.
- **Echo time** `TE` -- the time of the last readout sample, `TE = (n_t − 1)·dt`. The grid runs from the
  excitation (`t = 0`) to `TE`.
- **Timing budget** -- the dead times of a real scanner: the excitation lead-in, the refocusing window, the
  readout tail, a preparation offset (§3.3).
- **Encoding** -- the per-measurement parameters an analytical model reads (`b`, direction, δ, Δ, TE, …):
  what the gradient *means*, not what it is (§3.4).
- **Shape** -- a unit-amplitude gradient block as a function of its amplitude (§6.1).
- **Assembler** -- an RF schedule and the rule placing blocks around it (§6.2).
- **Family** -- a named shape × assembler with its literature parameters (§6.3).
- **Protocol** -- a tuple of acquisitions, one echo time each: a multi-TE scheme.

RFC-2119 keywords (MUST, SHOULD, MAY) have their usual meaning.

---

## 3. Data model

An acquisition has the following fields. Stored fields are the whole state; everything in §4 is a function
of them and MUST NOT be stored alongside (a stored copy can disagree with its source).

| field | type / shape | unit | meaning |
|---|---|---|---|
| `G` | float32 `(n_meas, n_t, 3)` | T/m | the PHYSICAL gradient, one row per measurement, on the grid |
| `dt` | float | s | the grid step; `T = TE = (n_t − 1)·dt` |
| `rf` | RF schedule (§3.2) | -- | the pulses, in time order; MAY be empty (a bare gradient echo) |
| `readout` | tuple of int | samples | the samples read; defaults to the last sample, or to every echo of a train |
| `timing` | budget (§3.3) or null | -- | the dead times the acquisition was built to |
| `encoding` | encoding (§3.4) or null | -- | the analytical parameters, `n_meas` rows |
| `crusher` | crusher (§3.5) or null | -- | the emergent voxel-scale crusher the vector-Bloch route models |
| `family` | string | -- | the family name (§6.3), `"waveform"` for an arbitrary played gradient |
| `build_spec` | `(name, kwargs)` or null | -- | provenance: the builder and its arguments; rebuilding reproduces `G` |
| `prescription` | prescription (§3.6) or null | -- | where in the bore and on what voxels the acquisition images; nothing is derived from it |

`n_meas ≥ 1`; `n_t ≥ 2`. A measurement is one row of `G`: one direction and amplitude of the same
sequence. **One acquisition has one RF schedule and one echo time**; a scheme with several echo times is a
**Protocol** (a tuple of acquisitions), never a per-measurement `TE` array with a single schedule.

### 3.1 The grid

The grid is `t_k = k·dt`, `k = 0 .. n_t − 1`. Sample `k` acts over the step `[k·dt, (k+1)·dt)`: a walker's
phase over that step is `γ·dt·G_eff[k]·r_k` (`RPK.md` §6.1). Consequently **the readout sample `n_t − 1`
acts over nothing** -- it is the instant the signal is read -- and a conformant builder leaves it at zero.
A shape (§6.1) is sampled at the **middle** of each step, `(k + ½)·dt`, so that a lobe rasterises
symmetrically however short its ramps.

### 3.2 The RF event and the schedule

An RF event is the record

| key | type | unit | meaning |
|---|---|---|---|
| `t_s` | float | s | the pulse instant (its centre) |
| `flip_deg` | float | degrees | the nominal flip angle |
| `label` | string | -- | the role: `"Mz→Mxy"` / `"excitation"` (excite), `"store"`, `"recall"`, `"refocus"` / `"refocusing"`; `""` for unlabelled |
| `axis_deg` | float, default 0 | degrees | the B1 phase (0 = x, 90 = y) |
| `duration_s` | float, default 0 | s | the pulse duration; `0` is an instantaneous hard pulse |
| `offset_hz` | float, default 0 | Hz | the carrier offset over the pulse |
| `envelope` | optional `{b1_re, b1_im, dt}` | T, s | the played complex `B1(t)`; when present its length is the duration and `γ ∫|B1| dt` the flip |

A finite pulse **occupies the window** `[t_s − duration_s/2, t_s + duration_s/2]`. A hard pulse occupies its
instant and constrains nothing (§5.2).

The schedule is the tuple of events sorted by `t_s`. Its serialised record is the list of event records
(JSON; the Pulseq `dmipy_rf_events` definition of §9). A reader MUST reject a record with unknown keys or
a missing `t_s` / `flip_deg`; it MUST NOT accept a pulse spelled as anything but this record.

### 3.3 The timing budget

| key | unit | meaning |
|---|---|---|
| `t_prep` | s | the excitation centre (a preparation before it shifts the whole schedule); `0` by default |
| `t_excite` | s | the excitation duration; encoding MAY begin at the **lead-in** `t_lead = t_prep + t_excite` |
| `t_refocus` | s | the refocusing pulse duration: a window of that width centred on the pulse carries no gradient |
| `t_readout_pre_echo` | s | the readout tail before every readout sample that carries no gradient |
| `TE` | s or null | the echo time the budget was resolved for, when fixed |

`min_TE = max(2·(t_lead + t_refocus/2), 2·(t_readout_pre_echo + t_refocus/2))` is the smallest spin-echo
`TE` for which both encoding windows exist. A budget that names finite durations makes the pulses of the
schedule finite: the events of an acquisition built to a budget carry `duration_s = t_excite` (excitation,
store, recall) and `t_refocus` (refocusing).

### 3.4 The encoding

The per-measurement analytical view, all arrays of length `n_meas` unless noted: `bvalues` (s/m²),
`gradient_directions` (`(n_meas, 3)`, unit vectors), `TE` (s), `qvalues` (1/m), `gradient_strengths`
(T/m), `delta`, `Delta` (s), `minimum_te` (s), `te_auto` (bool), `tau_perp_SE` (s: time transverse),
`ste_flip_angles` (deg, 3-tuple), `ramp_time` (s), `oscillation_frequency` (Hz), `gradient_rise_time` (s),
`n_oscillation_cycles`, `gradient_duration` (s), `cpmg_n_echoes`, `cpmg_TE` (s), `cpmg_beta_deg`,
`n_t_per_echo`, `refocused` (bool). A field a family does not define is null.

`bvalues` MUST equal the b-value derived from the acquisition (§4.5) to relative `1e-6`; the encoding
restates, it never overrides.

### 3.5 The crusher

The emergent voxel-scale crusher of the vector-Bloch route: `{windows_s: [(t0, t1), ...], n_cycles}` -- a
dephasing across the voxel of `n_cycles` full cycles applied over each window, so the transverse
magnetisation that was not stored along z is destroyed there. Only the vector-Bloch route reads it.

### 3.6 The prescription

The acquisition in **space**, as the rest of the object is the acquisition in **time**:

```jsonc
{"isocenter_m": [0, 0, 0],            // the point the scanner is focused on
 "axes": "RAS",                       // the scanner direction each voxel index runs along (i -> +x, j -> +y, k -> +z)
 "voxel_size_m": [1.5e-3, 1.5e-3, 1.5e-3],
 "matrix": [40, 40, 1],
 "origin_m": [-0.02925, -0.02925, 0]} // scanner coordinate of the CENTRE of voxel (0, 0, 0); default: the FOV centred on the isocenter
```

It is OPTIONAL and **derives nothing**: `G`, the schedule, the b-value and the echo are the same with or
without it (§5.7 holds). What it fixes is the frame of §4.1 -- the gradient and B0 directions are given along
`axes` -- and the voxels: a consumer that bins a sample into voxels (a partitioned replay phantom whose grid is
attached to the bore, `RPH.md`) MUST use this prescription when the phantom declares no grid of its own, and a
consumer holding its own grid MUST refuse an acquisition prescribed on other `axes` rather than rotate either.
A Pulseq export carries it as the file's `FOV` definition (`matrix · voxel_size_m`) plus the full object under
`dmipy_prescription`; an import without either has no prescription.

---

## 4. Derivation rules (normative)

### 4.1 Units and frame

SI everywhere (`RPK.md` §4.1): T/m, s, Hz, m. `γ = 267.513·10⁶ rad/s/T` (proton). The acquisition is in the
**scanner frame**; a substrate's pose rotates the gradient and the field direction into the substrate frame
at consumption (`RPH.md` §4), never the stored `G`.

### 4.2 The gate and the effective gradient

The gate over the grid is

    s(t_k) = ∏_{p ∈ refocus ∪ recall} (−1)^{[ t_k ≥ t_p − τ ]},      τ = 10⁻¹² s,

the product over every pulse whose role is `refocus` or `recall` (an unlabelled pulse with
`|flip − 180°| < 20°` counts as refocusing). A grid time equal to a pulse instant up to floating rounding
IS that instant: the sample at the pulse is after it. The effective gradient is

    G_eff[k] = G[k] · s(t_k).

`s` is `±1` and its own inverse: the same rule un-folds an effective gradient into the physical one. A
consumer of the scalar-phase tiers (`RPK.md` §6.1, §6.4) reads `G_eff`; the vector-Bloch route reads `G`
and applies the pulses itself (§7).

### 4.3 The coherence state, the echoes, the mixing time

Magnetisation starts along z. Walking the schedule in time order with the state `transverse ∈ {no, yes}`:

- an `excite` (or a `recall`) while not transverse makes it transverse and sets the reference time
  `t_ref` to the pulse instant; a `recall` adds `t − t_store` to the **mixing time** `TM`;
- a `store` while transverse makes it longitudinal and records `t_store`;
- a `refocus` while transverse forms an **echo** at `2·t − t_ref` and sets `t_ref` to that echo.

A labelled pulse plays its label's role whatever its flip (a stimulated echo's store may be 60°); an
unlabelled pulse is read by its flip (90° excites, stores or recalls by state; 180° refocuses); other flips
are not tracked (they are the vector-Bloch route's). Each transition happens at the pulse's instant.

The **coherence mask** `χ(t_k)` is `1` while transverse and `0` while stored, binary for hard pulses. Across
a **finite** pulse's window it is the transverse fraction of the pathway averaged over the ensemble's
azimuth, with `θ` running across the pulse: an excitation or a recall tips z into the plane as `sin²θ`
(`θ: 0 → π/2`), a store tips the plane onto z as `cos²θ`, a 180 keeps the component along B1 transverse and
swings the perpendicular one through z, `½ + ½cos²θ` (`θ: 0 → π`) -- a quarter of the pulse spent
longitudinal. `χ ≡ 1` (no storage anywhere) is reported as *no mask*.

An acquisition is a **stimulated echo** when its schedule stores and recalls; its `TM` is the total
longitudinal storage time; the readout then carries the idealised **0.5** amplitude factor of the stored
half (§7.1). An empty schedule is transverse throughout and forms no echo: a bare gradient echo.

### 4.4 The readout and the echo

`readout` defaults to `(n_t − 1,)` -- the grid ends at the readout -- or, for a schedule forming two or more
echoes (a train), to the sample of every echo. A stated readout MUST lie in `[0, n_t)`, MUST agree with the
schedule's echoes within **2 samples** when the schedule forms any, and MUST list every echo of a train. The
**echo sample** is the last readout sample.

### 4.5 The b-value and the B-tensor

With `q_k = γ · Σ_{j ≤ k} G_eff[j]·dt` (rectangular accumulation, as the walk's phase),

    b   = ∫ |q|² dt   (trapezoidal rule over the grid),      B_ij = ∫ q_i q_j dt,

so `trace(B) = b`. These are THE integrals: every builder's declared `bvalues` is scaled to them exactly
and an analytical layer reads them.

### 4.6 The refocusing residual

The **net moment at a readout** is `q` accumulated over the steps *before* it (§3.1),
`q(t_i) = Σ_{k < i} G_eff[k]·dt`. The residual `max_i |q(t_i)| / max_k |q_k|` over the readout samples and
measurements is the refocusing residual.

---

## 5. Invariants (normative)

A conformant acquisition satisfies all of the following; a builder MUST verify them on every build and
refuse an acquisition that fails one, naming the failure.

### 5.1 Refocused at every readout

For every family that declares a refocusing pulse or a stimulated echo, and for a self-refocusing gradient
echo, the refocusing residual (§4.6) MUST be `≤ 10⁻³`. A reference builder realises it to rounding
(`≤ 10⁻⁹`) by construction: a spin echo plays the **same block** on both sides of its 180, so `G_eff` is
`B, −B` sample for sample.

### 5.2 No gradient through a finite pulse; a hard pulse constrains nothing

Every sample whose time lies in a finite pulse's window (§3.2) MUST carry zero gradient. A hard pulse
(`duration_s = 0`) imposes nothing: a constant gradient through an ideal 180 train is Carr-Purcell and is
legal.

### 5.3 The budget's dead times

With a timing budget, no **whole step** `[t_k, t_k + dt)` MAY lie inside the lead-in `[0, t_lead]` or a
readout tail `[t_ro − t_readout_pre_echo, t_ro]` before any readout sample while carrying gradient (a step
straddling a window edge is grid rounding, not a violation), and `TE ≥ min_TE`.

### 5.4 The readout is where the echo forms

§4.4. A readout off the schedule's echo is refused, not shifted.

### 5.5 The declared b is the played b

`encoding.bvalues` equals §4.5 to relative `10⁻⁶`; `encoding.gradient_directions` has `n_meas` rows.

### 5.6 One schedule, one TE

`rf` is one schedule and `TE = (n_t − 1)·dt` is one time. A multi-TE scheme is a Protocol.

### 5.7 Nothing is stored that §4 derives

`G_eff`, `s`, `χ`, `TM`, the echo times, the stimulated-echo state, the b-value and the refocusing
residual are derived on read and MUST NOT appear as stored fields or constructor flags.

---

## 6. Families: shapes × assemblers

### 6.1 Shapes

A shape is a unit-amplitude block sampled mid-step (§3.1). Its ramps follow the amplitude `g` and the slew
limit: `ε = g / slew_rate` on every edge, vertical (`ε = 0`) at `slew_rate = ∞`. **Slew is a limit, never a
fork**: the structure of a block is the same at every slew rate. A block that cannot reach its amplitude
inside its own span (a ramp longer than the lobe, a train lobe shorter than two ramps, a cosine whose own
slope `2πf·g` exceeds the limit) MUST be refused, never clipped.

| shape | parameters | block |
|---|---|---|
| `trapezoid` | `δ` (half-amplitude width), `ε` | ramp 0→1 over `ε`, flat to `δ`, ramp 1→0 over `ε`; span `δ + ε`; `ε = 0` is the square lobe |
| `trapezoid_train` | `N` lobes, lobe span `L`, `ε` | `N` adjacent lobes alternating `+, −, +, …`, each ramping over `ε` at both ends (Drobnjak et al. 2016); one lobe is a PGSE lobe |
| `cosine` | `n` whole periods, `f`, `ε` | `cos(2πft)` over `n/f`, DC-free in `q`, its edges ramped over `ε` |
| `bipolar` | `δ`, `Δ`, `ε` | a lobe and its negative `Δ` later (centre to centre): self-refocusing |
| `axis_pairs` | axes `{a_i}`, duration `σ`, `ε` | for each axis in turn, `+a_i` then `−a_i` over `σ / (2·n_axes)` each: each axis's `q` returns to zero inside its own pair, so `B` has no off-diagonal term |

### 6.2 Assemblers

An assembler is an RF schedule and where the blocks sit relative to its pulses. It computes the smallest
`TE` that fits (the blocks, their gaps, and the budget's windows) and accepts a longer one, adding the extra
time as dead time as stated below. Blocks are placed in whole steps: a block **ends** against a pulse (its
last step ends by the pulse instant) or **starts** after one (its first step starts at or after it).

| assembler | schedule | placement | TE |
|---|---|---|---|
| **spin echo** | 90 at `t_prep`, 180 at `TE/2`, echo at `TE` | the same block before and after the 180, each row's pair **centred on the 180** with the row's own gap between the blocks; the extra TE is dead time split equally outside the pair | `min TE = max_rows(2·span + gap) + 2·max(t_lead, t_readout_pre_echo)`; the gap MUST hold `t_refocus` |
| **stimulated echo** | 90, store, recall; no 180 | the block before the store (rows with shorter ramps end at the same store) and the same block after the recall; the time transverse before the store equals the time after the recall, `TE = 2·t_store + TM` | `min TE = 2·max span + t_excite + TM + 2·max(t_lead, t_readout_pre_echo)` |
| **gradient echo** | 90 only | one self-refocusing block starting at the lead-in; the rest of TE is free precession | `min TE = t_lead + span + t_readout_pre_echo` |
| **echo train** | 90, then `n` pulses of `β` at `(k + ½)·TE_echo`, an echo read at every `k·TE_echo` | the gradient fills every stretch the pulses and readouts leave, one lobe per stretch, a stretch never crossing an echo, at **constant** polarity or **alternating** per interval; the same dead time (the larger of lead-in and readout tail) on both sides of every echo so each interval's halves balance | `TE = n·TE_echo`, fixed |

### 6.3 The families

Every family takes the measurement axis -- `gradient_directions` and exactly one of `bvalues` (the b to
realise: the amplitude is iterated, the ramps following it, then scaled so §4.5 holds exactly) or
`gradient_strengths` (the amplitude to play: the b follows) -- an optional `TE`, the grid `n_t`, the
`slew_rate` limit and an optional timing budget.

| family | assembler | shape | parameters | notes |
|---|---|---|---|---|
| `pgse` | spin echo | trapezoid | `δ`, `Δ` (centre to centre) | gap `= Δ − δ − ε`; both lobes the same sign (the 180 folds the second) |
| `pgste` | stimulated echo | trapezoid | `δ`, `TM`, `ste_flip_angles` | `Δ = δ + TM`; the recall's sign flip folds the second lobe; `tau_perp_SE = TE − TM` |
| `ogse` | spin echo | trapezoid_train or cosine | `f`, `σ` (block duration), `shape`, optional `Δ` | `2fσ` lobes (trapezoid) or `fσ` periods (cosine) MUST be whole -- refused, not snapped; gap `= Δ − σ` when `Δ` is given, else the refocusing window; the same block on both sides so `G_eff` continues the oscillation |
| `cpmg` | echo train | (fills the stretches) | `n_echoes`, `TE_echo`, `polarity`, `β` | `bvalues` is the b of the whole train; with neither `bvalues` nor `gradient_strengths` the train carries no gradient (a pure-T2 train) |
| `gre` | gradient echo | bipolar | `TE`, optional `δ`, `Δ` | no 180: `G_eff = G`, a static field is not refocused, the readout is complex |
| `ste` | gradient echo | axis_pairs (x, y, z) | `σ` | spherical tensor encoding, `b_δ = 0`: `B = (b/3)·I` |
| `pte` | gradient echo | axis_pairs (two in-plane axes) | plane normal, `σ` | planar tensor encoding, `b_δ = −½` |
| `waveform` | -- | as played | `G`, `dt`, directions | an arbitrary played gradient with no declared pulses: MUST refocus on its own (§5.1) |

Readers of a played gradient that declare its schedule exist for the spin echo (`from_btensor_waveform`: the
180 at `TE/2`, the only instant at which the static field refocuses at the echo) and the stimulated echo
(`from_pgste_waveform`: the two transverse periods MUST match, `τ₁ = τ₃`).

---

## 7. What a consumer reads

### 7.1 Scalar-phase tiers (the fused walk, the pack's `replay`, `pose_response`)

Read `G_eff` (§4.2) for the gradient phase and `χ` (§4.3) for relaxation: `T2` acts while transverse, `T1`
while stored (`RPK.md` §6.2); the static field term of `RPK.md` §6.4 is gated by `s`. A stimulated echo's
signal carries the factor `0.5`. The readout is the echo sample; a train returns one value per echo.

### 7.2 The vector-Bloch route (`simulate_bloch`, the pack's `replay_bloch`, a phantom's `replay_bloch`)

Reads the PHYSICAL `G` and applies the schedule itself: every pulse as a rotation about its B1 axis (a
finite pulse over its window with its carrier; an envelope when given), the crusher (§3.5) over its
windows, relaxation throughout. A transmit scale multiplies every flip. It returns the transverse
magnetisation at the readout sample, or at every echo of a train. This route MUST be used for anything the
gate cannot express -- a non-nominal flip, a finite pulse's pathways, an MT saturation -- and a scalar-tier
consumer MUST refuse such an acquisition rather than approximate it (`RPK.md` §6.6).

### 7.3 Analytical layers

Read `encoding` (§3.4). An analytical model MUST NOT re-derive an encoding from `G`; it reads the declared
one, whose `bvalues` are guaranteed by §5.5.

### 7.4 No side channels

A consumer takes the acquisition and nothing beside it: not its RF schedule, not its echo samples, not
its refocus time as a separate argument. (A kernel that operates on bare arrays -- a trajectory, a gradient
array, a list of events -- is not a consumer of an acquisition and is exempt.)

---

## 8. Scanner limits (data, not solver)

An acquisition is built against a scanner's limits, which are **data** with a citation. The catalogue is
one JSON document with the sections

| section | content |
|---|---|
| `citations` | `key → {title, url / doi, accessed}` |
| `scanners` | per model: `vendor`, `model`, `field_T`, and the leaves `gradient.max_amplitude`, `gradient.max_slew_rate`, `gradient.gradient_raster_time`, `rf.*`, each a **cited leaf** |
| `safety` | the SAR / `B1_rms` / dB/dt (PNS) limits of the standard cited, and the coefficients of a PNS model (SAFE) where published |
| `envelopes` | declared limit points that are not a machine (a certificate's clinical / insert envelope, a Pulseq example system) |
| `classes` | the short names of scanner classes and which model they resolve to |
| `aliases` | short names for models |

A **cited leaf** is `{value, unit, field_T, context, source_key, location, confidence}` with `confidence ∈
{cited, inferred, assumed}`. A typed SI view (`ScannerLimits.of(name, regime)`) resolves a name through the
aliases and classes and converts the leaves; `regime = "diffusion"` returns the PNS-derated slew where one
is catalogued. No other module carries a scanner number. A PNS solver that reads these coefficients is a
designer's concern and lives outside this specification.

---

## 9. Pulseq interchange

**Export.** The PHYSICAL gradient is written as arbitrary-gradient blocks on the acquisition's own raster,
the pulses as block pulses of one raster at their instants, the ADC at the readout. A pulse whose instant
falls in a live raster is **inserted** -- the gradient pauses and resumes unchanged, the sequence lengthens
by one raster per inserted pulse -- and the exporter MUST say so; a pulse on the boundary between a free
raster and a live one is played in the free raster before it. An acquisition built to a budget exports
without insertion. The `[DEFINITIONS]` carry `dmipy_dt`, `dmipy_n_t`, `dmipy_echo_idx`, `dmipy_gradient =
"physical"`, `dmipy_rf_events` (the schedule's exact record, §3.2) and `dmipy_timing` (the budget, §3.3).

**Import.** The gradient is rasterised exactly onto a uniform grid anchored at the excitation; between
gradient events the gradient is zero. The schedule is read from the RF blocks; a file this specification
exported also carries the exact schedule in its metadata and the reader MUST prefer it (the blocks are its
raster). A foreign spin echo (one 90, one 180, one ADC) MAY have a timing budget derived from its blocks.
An older file whose gradient was written folded (no `dmipy_gradient` definition and more pulses in the
metadata than in the blocks) is un-folded through §4.2 on read.

---

## 10. Conformance

A **conformant producer** (a builder, a designer, an importer) MUST: store the physical gradient; declare
every pulse as an RF event with its role; leave the readout sample at zero; verify §5 on every build and
refuse rather than repair; take exactly one of `bvalues` / `gradient_strengths`; treat the slew rate as a
limit (§6.1); refuse a shape that cannot be played rather than clip it.

A **conformant consumer** MUST: read the acquisition and nothing beside it (§7.4); derive `G_eff`, `χ`, the
echoes and `TM` by §4 and never from a stored flag; read the encoding rather than re-derive it; refuse an
acquisition whose RF the scalar gate cannot express instead of approximating it; refuse a bare gradient
array (it says nothing about its pulses).

### 10.1 Reference implementation

`dmipy-sim` (Apache-2.0): `dmipy_sim/acquisition/scanner_sequence.py` (the object, §3--§5),
`dmipy_sim/acquisition/rf.py` (the RF event and schedule, §3.2, §4.2--§4.3),
`dmipy_sim/acquisition/timing.py` (the budget, §3.3), `dmipy_sim/sequences/assemble.py` (shapes and
assemblers, §6.1--§6.2), `dmipy_sim/sequences/builders.py` (the families, §6.3),
`dmipy_sim/acquisition/scanner_constants.json` and `scanners.py` (§8), `dmipy_sim/sequences/pulseq.py` (§9).
Its `tests/test_api_surface.py` locks §5.7 and §7.4 by inspection of the package.

---

## Appendix A. Change log

- **0.1.0** -- first draft: one object, the derivation rules, the invariants, shapes × assemblers, the
  scanner catalogue, Pulseq interchange.
