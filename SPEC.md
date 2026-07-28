# The Replay Pack Specification (`.rpk`)

**Version:** 0.1.0 (draft for comment)
**Status:** Draft — stable enough to implement against; field names and metadata keys are frozen for the `1.x` container schema, semantics may be clarified.
**Container schema described:** `rpk_schema_version = "1.1"`
**License of this document:** CC-BY-4.0. **License of the reference code:** Apache-2.0.

The key words **MUST**, **MUST NOT**, **REQUIRED**, **SHALL**, **SHOULD**, **SHOULD NOT**, **MAY**, and **OPTIONAL** in this document are to be interpreted as described in RFC 2119.

---

## 1. Scope and purpose

A **Replay Pack** (`.rpk`) is a portable, self-describing, code-free file that stores the *state of a Monte-Carlo random walk through a diffusion-MRI substrate* in a form that lets a **replayer** reconstruct the measured signal for **many different acquisitions and physical settings** without re-running the walk.

This specification defines:

1. the **replay model** — the physical invariant that makes replay valid (§3);
2. the **data model** — the named channels a pack MAY contain, their units, shapes, and physical definitions (§5);
3. the **replay operations** — the exact, reference math each channel feeds, so that *any* conformant replayer reproduces *any* conformant pack (§6);
4. **capability tiers** — how a generator that models only some physics declares a partial but fully interoperable pack (§7);
5. the **generator invariants** a producer MUST satisfy (§8);
6. **codecs / compression** as an *orthogonal* per-channel encoding layer — a pack is well-defined with raw, uncompressed channels (§9);
7. the **metadata schema**, **declared envelope**, and **self-certification** (§10–11);
8. the **container format** and **conformance rules** (§12–13).

This specification does **not** standardize: how a generator *produces* trajectories (any Monte-Carlo, lattice, or SDE integrator is admissible), the substrate geometry format, or the fitting/analysis performed on replayed signals.

### 1.1 Why a standard

Producing a converged random walk through realistic microstructure is the expensive step in Monte-Carlo diffusion MRI; evaluating a pulse sequence against a *stored* walk is cheap. The economically and scientifically useful unit is therefore **"one walk, replayed many times."** A shared on-disk contract lets:

- different simulators emit interoperable packs (even partial ones — see tiers, §7);
- a single reference replayer serve packs from any producer;
- packs be deposited, cited, and pooled in a public **substrate bank**;
- lossy compression be introduced (or not) without changing the format or the replay semantics.

---

## 2. Terminology

- **Walker / spin.** One random-walk realization. Index `i ∈ [0, N_w)`.
- **Trajectory.** The stored position of a walker at each save time: `r_i(t_k)`.
- **Save grid.** The `N_t` times `t_k = k · Δt`, `k ∈ [0, N_t)`, at which positions are stored. `Δt` is the **save interval** (`dt_traj`), not necessarily the integrator's internal sub-step.
- **Channel.** A named array in the pack (e.g. `positions`, `compartment`). Some are REQUIRED, most OPTIONAL.
- **Replay knob.** A physical parameter applied *at replay time* — gradient waveform `G(t)`, static field `B0`, substrate orientation, `T2`/`T1`, surface relaxivity `ρ`, susceptibility, magnetization-transfer pools. Knobs are **not** stored per value; they parameterize the replay operation.
- **Replay operation.** A deterministic function mapping (channels, knobs) → signal (§6).
- **Capability tier.** A named group of channels that unlocks one class of replay (§7).
- **Declared envelope.** The region of knob space over which the producer certifies the pack's replay to a stated fidelity (§11).
- **Codec.** A per-channel, reversible-or-lossy encoding of a channel's array (§9). The identity codec (raw array) is always valid.
- **MC noise floor.** The irreducible variance of a finite-`N_w` ensemble estimate; the reference against which any lossy codec's error is judged (§9.3).

---

## 3. The replay model (the invariant)

Replay rests on one physical fact.

> **Replay invariant.** For a fixed substrate geometry, a fixed (possibly per-compartment) free diffusivity `D`, and a fixed pseudo-random seed, the ensemble of walker trajectories `{r_i(t_k)}` and the boundary-contact record are **independent of** the gradient waveform `G(t)`, the static field `B0`, the substrate orientation, the relaxation times `T2`/`T1`, the surface relaxivity `ρ`, the susceptibility distribution, and the magnetization-transfer parameters.

All of those quantities enter the *signal* only through a phase or a (log-)weight accumulated **along** an already-determined trajectory. Therefore they are **replay knobs**: the trajectory is computed once and stored; the signal for any knob setting is a cheap functional of the stored channels (§6).

Two consequences are used throughout:

- **(TE-prefix property.)** A walk stored to `T_max` contains every shorter echo time as a leading prefix: the first `N_t' = ⌈TE'/Δt⌉` saves are a valid, converged walk to `TE' ≤ T_max`. A replayer MAY slice.
- **(Separation of state and knobs.)** A pack stores only *walk state*. It never stores a signal, a b-value, or a waveform. Signals are always derived at replay time.

The invariant defines the boundary of the format. Anything a producer cannot recover by replay from stored channels — most importantly **a change of geometry, `D`, or seed** — is **out of scope**: it requires a new walk and hence a new pack. Diffusivity is therefore recorded as a fixed pack property (§10), not a knob.

---

## 4. Conventions

### 4.1 Units — SI everywhere, no implicit scale factors

Two kinds of quantity carry units in this standard, both fixed to SI so that a pack from one generator reproduces **identically** under another implementation's replayer. A producer MUST convert from any internal units on write; a replayer MUST assume SI.

**(a) Quantities stored in the pack** — the format proper.

| Stored quantity | Symbol | Unit | Channel / field |
|---|---|---|---|
| position | `r` | **metre (m)** | `positions` (continuous lab frame, §4.2) |
| save interval, echo time | `Δt`, `T_max` | second (s) | `dt_traj`, `T_max` |
| diffusivity | `D` | m²·s⁻¹ | `walk_params.diffusivity` (fixed pack property, not a knob) |
| boundary local time | `ℓ` | metre (m) | `boundary_local_time`, at `ρ/D = 1` (§6.3) |
| relaxation times | `T1`, `T2` | second (s) | `per_comp` |
| bound fraction, spin weight | `b_i`, `w_i` | dimensionless | `bound_fraction`, `spin_weights` |
| field-map grid spacing | — | metre (m) | `walk_params.cell_size` |

A custom `x_` channel (§5.2) MUST declare its own unit in metadata.

**(b) Replay-operation inputs (knobs)** — **not stored**; the caller supplies these to the normative replay operations (§6). Their units are fixed here *because the operations are normative*: without a fixed convention the same pack and the "same" acquisition would give different signals across implementations.

| Knob | Symbol | Unit | Used in |
|---|---|---|---|
| gradient waveform | `G(t)` | tesla·metre⁻¹ (T·m⁻¹) | §6.1 |
| gyromagnetic ratio | `γ` | rad·s⁻¹·T⁻¹ (¹H = 2.675 2219 × 10⁸) | §6.1, §6.4 |
| static field | `B₀` | tesla (T) | §6.4 |
| **volume susceptibility** | `χ`, `Δχ` | **dimensionless, SI — NOT CGS (differ by 4π)** | §6.4 |
| orientation, azimuth | `θ`, `α` | radian (rad) | §6.4 |
| surface relaxivity | `ρ` | metre·s⁻¹ | §6.3 |
| *(derived)* off-resonance, phase | `Δω`, `φ` | rad·s⁻¹, rad | §6.1, §6.4 |

**The b-value is neither stored nor a replay input.** It is a quantity the *caller* derives from its own `G(t)`; the replay operation (§6.1) consumes `G` directly and never sees a b-value, so no b-value convention is imposed on producers or replayers. A b-value appears in this standard only as the optional human-readable descriptor `replay_envelope.acquisition.b_max`, whose unit is **s·m⁻² (SI, not s·mm⁻²)** for that field alone.

### 4.2 Frame and positions
- **Frame.** Positions are in a fixed **laboratory frame**. Substrate *orientation* is a replay knob applied by rotating the **waveform** (and field direction), never by rotating stored positions (which would invalidate periodic and field-map channels).
- **Continuous.** Stored positions MUST be the *unwrapped* lab-frame trajectory (§8.2). Periodicity, if any, is a property of auxiliary field-map channels, not of `positions`.

### 4.3 Time
The save grid is uniform with interval `Δt = dt_traj`. `N_t` and `T_max = (N_t − 1)·Δt` are recorded. The integrator's internal sub-step is a producer concern and is neither stored nor standardized.

### 4.4 Signal and sign convention
The canonical signal is complex, `S = ⟨w · e^{+iφ}⟩`, with phase `φ = γ ∫ G·r dt` (§6.1), a **positive** gyromagnetic sign (`γ > 0`), and the `e^{+iφ}` rotation sense. Magnitude DWI is insensitive to this choice, but it is fixed here so the **complex** signal and the susceptibility phase (§6.4) are reproducible across implementations. A real spin-echo value is `|S|` or `⟨w·cos φ⟩`. Averages `⟨·⟩` are weighted means over walkers (weights `w_i`, default `1`).

### 4.5 Grids
Field-map channels (§5.2) are sampled at the walker's position **wrapped** into the field cell, `r mod cell_size`, on a uniform grid whose axes are in metres and whose spacing is `walk_params.cell_size`.

---

## 5. Data model — channels

A pack is a set of named arrays (**channels**) plus a metadata object (§10). Array element types are given as NumPy dtypes; a codec (§9) MAY change the *stored* representation but the **decoded** array MUST match the type/shape/units below.

### 5.1 Required channel

| Channel | Shape | Dtype (decoded) | Units | Definition |
|---|---|---|---|---|
| `positions` | `(N_w, N_t, 3)` | float32 (float64 permitted) | metres | Continuous lab-frame position `r_i(t_k)` of walker `i` at save `k`. |

`positions` is the **irreducible** channel. A pack with only `positions` (plus metadata) is valid and supports the Gradient tier (§7).

### 5.2 Optional channels

| Channel | Shape | Dtype | Units | Definition | Unlocks |
|---|---|---|---|---|---|
| `spin_weights` | `(N_w,)` | float32 | — | Per-walker statistical weight `w_i` (e.g. compartment volume weighting). Default `1`. | all tiers |
| `compartment` | `(N_w, N_t)` | int16 | label | Compartment id of walker `i` at save `k`. Id **0 is the extra-cellular / free compartment by convention**; positive ids are producer-defined and described in metadata `per_comp`. | Relaxation |
| `boundary_local_time` | `(N_w, N_t)` | float32 | m (see below) | Per-step accumulated wall-contact measure with **`ρ/D = 1`** (dimension of length ÷ diffusivity is folded so that replay multiplies by the desired `ρ/D`). See §6.3. | Surface |
| `bound_fraction` | `(N_w, N_t)` | float32 | — ∈ [0,1] | Fraction of walker `i` bound to the restricted (macromolecular) pool at save `k`. | Exchange/MT |
| `susc_field_0` | producer grid | float32 | per unit susceptibility | m=0 (isotropic / mean) component of the substrate's normalized off-resonance field map. | Field |
| `susc_field_C` | producer grid | float32 | " | ℓ=2 cosine (anisotropic) component. | Field |
| `susc_field_S` | producer grid | float32 | " | ℓ=2 sine (anisotropic) component. | Field |

Field-map channels are **normalized off-resonance basis maps** (§6.4): `susc_field_0` carries the isotropic part and `susc_field_{C,S}` the axially-anisotropic part, so the tier supports isotropic and/or anisotropic susceptibility (a fully general susceptibility tensor is a §14 extension). They are accompanied by `walk_params.cell_size` (grid spacing, m), the susceptibility scale(s) the maps are normalized to (the reference implementation stores the anisotropic scale as `walk_params.delta_chi_a`), and, where the substrate is oriented, per-compartment rotations in `per_comp.R` (each a 3×3 mapping the compartment's local frame to the lab frame).

A producer **MAY** define **additional** channels prefixed `x_` (e.g. `x_temperature`). Replayers **MUST** ignore unknown channels they do not implement (§14).

### 5.3 Channel invariants

- Every per-walker/per-save channel MUST share `N_w` and `N_t` with `positions`.
- `compartment` id `0` MUST denote the extra-cellular/free pool; other ids MUST be described in `per_comp`.
- `bound_fraction` MUST lie in `[0, 1]`, and at `k = 0` it MUST be the equilibrium occupancy (the walk is pre-burned-in; §8.8).
- `boundary_local_time` MUST be non-negative and expressed in the `ρ/D = 1` normalization of §6.3.

---

## 6. Replay operations (normative reference semantics)

These define the meaning of the channels. A **conformant replayer** MUST reproduce these operations (to the decoded arrays' numerical precision). `γ` is the proton gyromagnetic ratio. Sums over `k` run on the save grid with weight `Δt` (trapezoidal or left-rule is producer-agnostic *provided the producer used the same rule to converge the walk*; left-rule is the reference).

The operations split by their dependence on position: the gradient phase (§6.1) and the relaxation/surface log-weights (§6.2–6.3) are **linear** in the stored positions, whereas the off-resonance field-map lookup (§6.4) and the Bloch/MT evolution (§6.5) are **nonlinear**. This is informative only — the signal is identical either way — but it lets a replayer evaluate the position-linear operations directly in a linear codec's coefficient space without decoding the full trajectory (see `CODEC_REGISTRY.md`, "Replay in coefficient space").

### 6.1 Gradient phase (Gradient tier — always available)

For a gradient waveform `G(t_k)` (T·m⁻¹, resampled to the save grid) the accumulated phase of walker `i` is

```
φ_i = γ · Δt · Σ_k G(t_k) · r_i(t_k)
```

and the signal is the weighted ensemble mean

```
S(G) = Σ_i w_i · e^{i φ_i} / Σ_i w_i .
```

This yields any b-value, any b-tensor, OGSE/PGSE/arbitrary `G(t)`, and (by sweeping `q`) the ensemble average propagator. **No stored quantity depends on `G`** — this is the whole point.

### 6.2 Bulk relaxation `T2`/`T1` (Bulk-relaxation tier)

This tier replays the **intrinsic per-compartment** relaxation of each pool's water — distinct from surface relaxivity at walls, which is the separate Surface tier (§6.3).


With per-compartment transverse/longitudinal rates and the `compartment` channel `c_i(t_k)`, accumulate a log-weight

```
log W_i = − Δt · Σ_k [ χ_k / T2[c_i(t_k)]  +  (1 − χ_k) / T1[c_i(t_k)] ]
```

where `χ_k ∈ {0,1}` encodes transverse vs longitudinal periods of the sequence (a spin-echo is transverse throughout its encoding). The relaxation-weighted signal replaces `w_i → w_i · W_i` in §6.1. `T2[·]`, `T1[·]` are replay knobs supplied from `per_comp` defaults or overridden by the caller.

### 6.3 Surface relaxivity (Surface tier)

With the `boundary_local_time` channel `ℓ_i(t_k)` (stored at `ρ/D = 1`) and a chosen surface relaxivity `ρ` and diffusivity `D`,

```
log W_i += − (ρ / D) · Σ_k ℓ_i(t_k) .
```

Because the channel is stored normalized, one walk serves every `ρ`. This is a replay knob; `ρ = 0` recovers §6.2.

### 6.4 Static off-resonance / susceptibility (Field tier)

The Field tier stores a set of **normalized off-resonance field-basis maps** the producer precomputes for the substrate. The reference set is `{Φ_0, Φ_C, Φ_S}`: `Φ_0` is the isotropic (m=0) component and `Φ_C, Φ_S` are the ℓ=2 (cylindrically anisotropic) components. Together they represent an **isotropic and/or axially-anisotropic** susceptibility source — not only anisotropic. A fully general rank-2 susceptibility-tensor field is a documented extension requiring a richer basis set (§14).

At replay, for a field strength `B0`, a substrate orientation `θ` relative to `B0`, an azimuth `α`, and the susceptibility scale(s) `Δχ` the maps are normalized to, the off-resonance seen by walker `i` is sampled from the maps at the *wrapped* position `r_i(t_k) mod cell_size`:

```
Δω_i(t_k) = γ · B0 · [ χ_iso · Φ_0(r)  +  χ_aniso · sin²θ · (cos2α·Φ_C(r) + sin2α·Φ_S(r)) ]
```

`χ_iso` and `χ_aniso` are the isotropic and anisotropic susceptibility scales (the reference implementation exposes the anisotropic scale as `walk_params.delta_chi_a`; a purely isotropic source uses `Φ_0` alone).

The off-resonance contributes to the signal only through the **acquisition's phase history**, which is a replay knob (§6.6). The replayer supplies a per-step transverse-phase gate `s(t_k)` derived from the pulse sequence, and

```
φ^χ_i = Σ_k s(t_k) · Δω_i(t_k) · Δt ,      φ_i → φ_i + φ^χ_i
```

The susceptibility phase adds to the gradient phase inside the *same* exponential, so the diffusion×susceptibility covariance (the cross-term) is preserved. Field strength, orientation, susceptibility scale, and sequence are all knobs; one walk serves any combination.

### 6.5 Magnetization transfer / exchange (Exchange tier)

With `bound_fraction` `b_i(t_k)`, replay blends per-step relaxation and off-resonance toward a bound-pool set `(T2_b, T1_b, Δω_b)` by occupancy, within a **vector-Bloch replay** — RF pulses act as rotations of the magnetization vector, which the scalar log-weight model cannot represent (it has no `M_z` reservoir). Saturation transfer is *emergent*. The bound-pool parameters are replay knobs. This same vector-Bloch path is the general route for *arbitrary RF* (§6.6). See the reference implementation for the Bloch–McConnell blend.

**Equilibrium start.** An Exchange pack's `t=0` is the **bound-pool equilibrium**: the walk is generated pre-burned-in (§8.8), so replay begins from a fully-relaxed steady-state occupancy and never sees the fill-up transient.

### 6.6 The acquisition is a replay knob (no sequence is baked into the pack)

A pack stores **no** assumption about the pulse sequence. The full acquisition — the gradient waveform `G(t)` **and** the RF (flip angles, refocusing, storage) — is supplied at replay. This is what makes one walk serve GRE, spin-echo, CPMG, stimulated-echo/PGSTE, and arbitrary sequences alike.

For the **scalar-phase tiers** (§6.1, §6.4) the sequence enters through a per-step **transverse-phase gate** `s(t_k)` that the replayer derives from the sequence:

- `s ≡ +1` — gradient echo (GRE), no refocusing;
- `s` flips sign at each 180° refocusing pulse — spin echo, CPMG;
- `s = 0` while magnetization is stored along `z` — stimulated echo (STE/PGSTE);

so `s(t_k) ∈ {−1, 0, +1}` in the common cases. The **gradient phase (§6.1)** already encodes any bipolar/refocused *gradient* structure in the waveform's polarity and needs no gate; the gate applies to the **static off-resonance** term (§6.4), which a refocusing pulse inverts and a storage period freezes. (The transverse/longitudinal indicator `χ_k` of §6.2 is the same idea for relaxation.)

Arbitrary RF — non-180° flips, adiabatic pulses, MT saturation — cannot be reduced to a scalar sign and requires the **vector-Bloch replay** (§6.5). A conformant replayer MAY implement only the scalar-gate model; it MUST then **refuse** acquisitions that need full Bloch evolution rather than approximate them (§13).

---

## 7. Capability tiers

A producer that models only part of the physics still emits a **fully interoperable** pack; it simply declares fewer tiers. Tiers are **independent** — a pack MAY declare Field without Relaxation.

Metadata flags use **explicit, self-describing names** (§10).

| Tier | Required channels (beyond `positions`) | Replay unlocked | Metadata flag(s) |
|---|---|---|---|
| **T0 Gradient** | *(none)* | §6.1 — any `G(t)`, b-tensor, EAP | `gradient: true` (always) |
| **T1 Bulk relaxation** | `compartment` (+ `per_comp.T2`,`per_comp.T1`) | §6.2 — any `T2`/`T1` | `bulk_relaxation: true` |
| **T2 Surface** | `boundary_local_time` | §6.3 — any surface relaxivity | `surface_relaxivity: true` |
| **T3 Field** | `susc_field_{C,S,0}` (+ `cell_size`, susceptibility scale) | §6.4 — any `B0`, orientation, susceptibility | `field_offresonance: true`, `field_orientation: true` |
| **T4 Exchange** | `bound_fraction` | §6.5 — MT/exchange (vector-Bloch) | `magnetization_transfer: true` |

The declared tier set lives in `replay_envelope` (§10). A replayer asked for a knob outside the declared tiers MUST refuse with a clear "capability not present" error, **not** silently return an approximate result (§11, §13).

*On membrane permeability (out of scope in `1.x`).* Permeability changes the **trajectory** (a crossing is a per-event random draw), so it is not a replay knob at all — a pack is walked at one permeability and varying it requires new walks. The `1.x` envelope therefore does **not** carry a permeability flag; the fixed value the walk used SHOULD be recorded under `provenance`. If a future version finds a way to make crossing a replayable quantity, it will be added as a new tier (§14).

---

## 8. Generator invariants (what a producer MUST guarantee)

A pack is only replayable if the producer respected the model. A conformant producer:

### 8.1 Determinism
MUST record the `seed` and produce byte-reproducible trajectories from `(geometry, D, seed)` on a fixed implementation. The seed is metadata, not a knob.

### 8.2 Continuous positions
MUST store **unwrapped** lab-frame positions. If the walk used a periodic cell, the producer MUST unwrap before storing (undo per-save jumps `> L/2` on periodic axes; leave non-periodic axes untouched), so that the gradient phase in §6.1 is correct. Storing wrapped positions is a conformance failure: each wrap injects a spurious `q·L` phase.

### 8.3 TE-prefix consistency
MUST ensure the save grid is uniform and that any leading prefix is itself a converged walk to that shorter time (§3). Producers MUST NOT reorder or subsample walkers across the save axis.

### 8.4 No sequence assumption
A producer MUST NOT bake any pulse-sequence choice into the pack. The acquisition (gradients and RF) is applied at replay (§6.6); the producer's obligation is only to store the trajectory (and any tier channels) over the full walk to `T_max`. In particular there is no "refocusing time" a producer must record — refocusing is a property of the replayed sequence, not of the walk.

### 8.5 Compartment convention
MUST use id `0` for the extra-cellular/free pool and MUST list every other id with its `(T2, T1[, R])` in `per_comp`.

### 8.6 Units and frame
MUST honor §4 exactly. A producer using non-SI internal units MUST convert on write.

### 8.7 Honest envelope
MUST declare a `replay_envelope` it can defend, and, for any lossy codec, MUST attach a `fidelity` self-report (§9.3, §11). A producer MUST NOT declare a tier whose channel it did not actually populate.

### 8.8 Equilibrium start for Exchange packs
An Exchange-tier (T4) pack MUST begin from the **bound-pool equilibrium**. A Monte-Carlo binding walk started from an arbitrary state (e.g. all spins free) shows a transient while the bound fraction relaxes to its steady state `f_b = k_f/(k_f + k_r)`; that transient does not represent a fully-relaxed sample and must not appear in a replay. The producer MUST therefore **equilibrate** the binding dynamics and **discard** that preamble, saving the walk from `t=0 =` the equilibrated state (the `bound_fraction` channel accordingly starts at equilibrium). Replayers assume `t=0` is equilibrium (§6.5). The producer SHOULD verify the occupancy has reached steady state before saving. This applies only to the Exchange tier; the other tiers have no such initial-condition transient.

*Reference implementation.* dmipy-sim exposes this as `equilibrate_binding = 'auto' | 'burnin' | 'fast' | 'off'`: `'burnin'` (the `'auto'` default when MT is on) runs an **adaptive RF-off burn-in in ~dwell-sized chunks until the ensemble occupancy plateaus** (geometry-agnostic, with a convergence flag); `'fast'` seeds the equilibrium occupancy analytically from a known `S/V` when that is position-invariant; `'off'` keeps the legacy all-free start. A producer MAY use any method that satisfies the equilibrium-start requirement above.

---

## 9. Codecs and compression (orthogonal layer)

**The format does not depend on compression.** Every channel is defined by its *decoded* array (§5). How that array is *stored* is a per-channel **codec** choice. The identity codec — the raw array — is the normative baseline, and a pack stored entirely with identity codecs is fully conformant. This section defines only the codec **interface**; the concrete codecs and their exact stored-key layouts live in a separate, independently versioned document, [`CODEC_REGISTRY.md`](CODEC_REGISTRY.md), so that new compression methods (some still being developed) can be added without touching this frozen core.

A conformant codec MUST satisfy four rules:

1. **Declared.** The codec is named (with any parameters) in metadata `compression`, e.g. `{"method": "identity"}` or `{"method": "lowrank", "K": 64, "walker_preserving": true}`. A replayer that does not recognize the `method` MUST refuse (never guess). Recognized method names and their stored tensor keys are defined in the registry.
2. **Decodes to the contract.** Decoding MUST yield the exact channel array of §5 (dtype/shape/units). Channel *presence* for tiers (§7) and conformance (§13) is satisfied by the channel's raw key **or** the registry-defined stored keys of its codec.
3. **Identity is the baseline.** `identity` (the raw array under its channel name) is always valid and lossless, and every replayer MUST support it.
4. **Lossy ⇒ self-certified.** A pack using any lossy codec MUST include a `fidelity` object reporting the maximum decoded-vs-raw replay error over a declared **acquisition battery**, measured against the **split-half Monte-Carlo noise floor** (the ensemble split into two *random* halves — random because walkers are typically seeded in compartment order). The convention: "lossless to the noise floor" means `err_max ≤ tol · floor_max` (reference `tol = 2`). Lossless codecs SHOULD report `err_max = 0`.

One structural rule belongs in the core because it constrains tiers, not any particular algorithm:

> **Walker-preserving requirement.** A codec that **resamples walkers** (a *distributional* codec — it stores a coefficient *distribution*, not per-walker coefficients) breaks per-walker channel alignment. Such a pack **MUST NOT** carry any per-walker channel (`compartment`, `boundary_local_time`, `bound_fraction`, `spin_weights`) and **MUST** declare only the Gradient tier. Per-walker channels REQUIRE a walker-preserving codec.

Compression is thus a *quality claim about a channel*, never a change to the replay contract. See `CODEC_REGISTRY.md` for the current registered codecs (`identity`, and reference-implementation methods) and their stored-key tables.

---

## 10. Metadata schema

Metadata is a JSON object embedded in the container (§12) and validated by `schema/rpk_metadata.schema.json`. Required and optional keys:

```jsonc
{
  "rpk_schema_version": "1.1",          // REQUIRED, semver of the container schema
  "id": "canonical-wm/g070-f055-3T",    // REQUIRED, stable pack identifier
  "walk_params": {                       // REQUIRED
    "n_walkers": 120000, "n_t": 200,
    "dt_traj": 2.65e-4, "T_max": 0.053,  // s
    "diffusivity": 0.6e-9,               // m^2/s, FIXED (not a knob)
    "seed": 0,
    "cell_size": 5e-6,                   // m, field-map grid spacing; present iff Field tier
    "delta_chi_a": -0.1e-6               // anisotropic susceptibility scale the Field maps
                                         //   are normalized to; present iff Field tier
  },
  "per_comp": {                          // REQUIRED iff Bulk-relaxation tier; else null
    "T2": [0.08, 0.05],                  // s, index = compartment id (0 = extra/free)
    "T1": [1.0, 0.8],                    // s
    "R": null                            // OPTIONAL 3x3 rotation(s) mapping each compartment's
                                         //   local frame -> lab frame, used to orient the
                                         //   anisotropic Field maps (§6.4); null if unoriented
  },
  "compression": {                       // REQUIRED — codec descriptor (see CODEC_REGISTRY.md)
    "method": "lowrank", "K": 64, "walker_preserving": true
  },
  "replay_envelope": {                   // REQUIRED — declared capabilities + domain
    "gradient": true,                    // T0 — always true
    "bulk_relaxation": true,             // T1 — intrinsic per-compartment T1/T2
    "surface_relaxivity": false,         // T2 — any surface relaxivity
    "field_offresonance": true,          // T3 — any B0 / susceptibility off-resonance
    "field_orientation": true,           // T3 — any substrate orientation w.r.t. B0
    "magnetization_transfer": false,     // T4 — MT/exchange (vector-Bloch)
    "diffusivity_fixed": true,           // D is a fixed pack property, not a knob
    "acquisition": { "b_max": 3.0e9, "ogse_periods": [], "B0_list": [3.0, 7.0] }
  },
  "fidelity": {                          // REQUIRED iff any lossy codec (§9)
    "err_max": 0.0015, "floor_max": 0.0074, "within_2x_floor": true,
    "battery": "..." },
  "provenance": {                        // RECOMMENDED — how the walk was made
    "generator": "dmipy-sim", "generator_version": "...",
    "geometry": "PackedMyelinatedCylinders", "real_or_synthetic": "synthetic",
    "permeability": 0.0 },               // the FIXED permeability the walk used (§7), if any
  "license": "CC-BY-4.0",                // REQUIRED — the SOURCE substrate's license
  "citation": "Fick RHJ, dmrai-lab (2026) ..."   // REQUIRED
}
```

Rules: `diffusivity` and `seed` are fixed pack properties, not knobs. `license` records the **source substrate's** license and MUST NOT relicense upstream geometry. Any tier flag set `true` in `replay_envelope` MUST have its channels present (§7) and its `per_comp`/`walk_params` fields populated. The envelope flags use explicit, self-describing names; *compatibility:* across the `1.x` line a reader SHOULD also accept the pre-rename aliases `T1T2`/`relaxation→bulk_relaxation`, `rho→surface_relaxivity`, `B0_any→field_offresonance`, `orientation_any→field_orientation`, `mt→magnetization_transfer` (and ignore the retired `rf`/`permeability` flags).

---

## 11. Declared envelope and self-certification

A pack is a **certificate**, not a black box. Two objects make its guarantees explicit:

- **`replay_envelope`** — the *domain of validity*: which tiers, and the acquisition/field range (`b_max`, OGSE periods, `B0_list`, …) over which the producer stands behind the pack. A request inside the envelope is guaranteed; a request outside is a *known limit*, and the replayer MUST surface it as such (§13), never silently extrapolate.
- **`fidelity`** — the *measured* error of any lossy encoding against the noise floor (§9.3).

The substrate bank (a separate service) accepts a pack by validating it against this specification (§13) and MAY re-measure `fidelity` independently.

---

## 12. Container format

- **Arrays:** a single **safetensors** file. Safetensors is REQUIRED because it (a) carries **no executable code** — a pack is safe to download and memory-map from untrusted sources — (b) is strongly typed, and (c) is zero-copy. Under the identity codec each channel is one tensor under its channel name (§5); under another codec a channel's tensor keys are those the registry defines for that method (`CODEC_REGISTRY.md`).
- **Metadata:** the JSON object of §10, serialized as a single string and stored in the safetensors `__metadata__` header under the key **`"rpk"`**. Producers MUST write it there; a sibling `.json` copy is OPTIONAL for human inspection. *Compatibility:* across the `1.x` line a reader MUST also accept the legacy header key `"json"` (used by the reference implementation before this spec fixed the canonical key); producers SHOULD migrate to `"rpk"`.
- **Dataset sidecar (OPTIONAL):** a Croissant (`schema.org/Dataset` + MLCommons) `.croissant.jsonld` carrying `license`, `citation`, `provenance`, and `replay_envelope` for dataset-catalog interoperability.
- **Extension:** `.rpk`. **Naming:** the file basename SHOULD equal the last path segment of `id`.

A pack MUST contain no code and MUST be openable without executing anything.

---

## 13. Conformance

A file is a **conformant Replay Pack** iff:

1. it is a valid safetensors file with a `__metadata__` block that validates against `schema/rpk_metadata.schema.json`;
2. it contains the REQUIRED `positions` channel (raw or under a declared, decodable codec), with `(N_w, N_t, 3)` shape after decode;
3. every tier flagged `true` in `replay_envelope` has all its §7 channels present and its §10 fields populated;
4. all per-walker/per-save channels share `(N_w, N_t)` with `positions`;
5. if any codec is lossy, a `fidelity` object is present (§9.3);
6. `license` and `citation` are present.

A **conformant replayer**:

- implements §6.1 and at least refuses cleanly for tiers it does not implement;
- decodes every codec it claims to support to the §5 array contract before applying §6;
- for any requested knob outside the pack's `replay_envelope`, returns a **capability/domain error**, never a silent approximation (§11).

A **conformant producer** satisfies §8 and §9.3.

Conformance is layered: a producer MAY implement only T0; a replayer MAY implement only T0; they still interoperate on T0.

---

## 14. Versioning and extensibility

- `rpk_schema_version` is semantic. **Minor** bumps add OPTIONAL channels/keys; a replayer MUST ignore channels and metadata keys it does not recognize (forward-compatible). **Major** bumps may change required semantics.
- New channels intended for later standardization SHOULD be introduced under the `x_` prefix first.
- New tiers extend the table in §7 without altering existing tiers.

---

## 15. Security considerations

A pack carries data only. Safetensors executes no code on load, and the metadata and Croissant sidecar are inert JSON. Replayers MUST NOT `eval`/`exec` any pack content. Field-map and position arrays are bounded numeric data; a replayer SHOULD sanity-check shapes against `walk_params` before allocating.

---

## 16. References

- RFC 2119 — Key words for requirement levels.
- safetensors — https://github.com/huggingface/safetensors
- Croissant (MLCommons) — http://mlcommons.org/croissant/
- The Replay Pack methodology paper — *(forthcoming; this specification claims priority on the infrastructure/format.)*

---

## Appendix A. Minimal conformant pack (Tier 0)

Arrays: `positions` `(N_w, N_t, 3)` float32.
Metadata: `rpk_schema_version`, `id`, `walk_params{n_walkers,n_t,dt_traj,T_max,diffusivity,seed}`, `compression{method:"identity"}`, `replay_envelope{gradient:true, ...all others false...}`, `license`, `citation`.
Replay: §6.1 only. This is the floor of interoperability — a walk and nothing else — and every richer pack is this plus declared channels. See `examples/`.
