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

- **Units are SI.** Length in metres (m), time in seconds (s), magnetic field in tesla (T), gyromagnetic ratio `γ` in rad·s⁻¹·T⁻¹, gradient in T·m⁻¹.
- **Frame.** Positions are in a fixed **laboratory frame**. Substrate *orientation* is a replay knob applied by rotating the **waveform** (and field direction), never by rotating stored positions (which would invalidate periodic and field-map channels).
- **Positions are continuous.** Stored positions MUST be the *unwrapped* lab-frame trajectory (see §8.2). Periodicity, if any, is a property of auxiliary field-map channels, not of `positions`.
- **Time.** The save grid is uniform with interval `Δt = dt_traj`. `N_t` and `T_max = (N_t − 1)·Δt` are recorded. The integrator's internal sub-step is a producer concern and is neither stored nor standardized.
- **Complex signal.** The canonical signal is complex, `S = ⟨w · e^{iφ}⟩`; a real spin-echo magnitude is `|S|` or `⟨w·cos φ⟩` as appropriate. Averages `⟨·⟩` are weighted means over walkers with weights `w_i` (§5, `spin_weights`; default `w_i = 1`).

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
| `susc_field_C` | producer grid | float32 | dimensionless (per unit `Δχ`) | ℓ=2 cosine component of the substrate's normalized off-resonance field map. | Field |
| `susc_field_S` | producer grid | float32 | " | ℓ=2 sine component. | Field |
| `susc_field_0` | producer grid | float32 | " | m=0 (monopole/mean) component. | Field |

Field-map channels are accompanied by the metadata keys `walk_params.cell_size` (grid spacing, m) and `walk_params.delta_chi_a` (the anisotropic susceptibility scale `Δχ`) and, where the substrate is oriented, a per-compartment rotation in `per_comp.R`.

A producer **MAY** define **additional** channels prefixed `x_` (e.g. `x_temperature`). Replayers **MUST** ignore unknown channels they do not implement (§14).

### 5.3 Channel invariants

- Every per-walker/per-save channel MUST share `N_w` and `N_t` with `positions`.
- `compartment` id `0` MUST denote the extra-cellular/free pool; other ids MUST be described in `per_comp`.
- `bound_fraction` MUST lie in `[0, 1]`.
- `boundary_local_time` MUST be non-negative and expressed in the `ρ/D = 1` normalization of §6.3.

---

## 6. Replay operations (normative reference semantics)

These define the meaning of the channels. A **conformant replayer** MUST reproduce these operations (to the decoded arrays' numerical precision). `γ` is the proton gyromagnetic ratio. Sums over `k` run on the save grid with weight `Δt` (trapezoidal or left-rule is producer-agnostic *provided the producer used the same rule to converge the walk*; left-rule is the reference).

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

### 6.2 Relaxation `T2`/`T1` (Relaxation tier)

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

### 6.4 Susceptibility / static off-resonance (Field tier)

Given the normalized field-map components `{Φ_C, Φ_S, Φ_0}`, a field strength `B0`, a substrate orientation `θ` relative to `B0`, and an anisotropic scale `Δχ = delta_chi_a`, the off-resonance experienced by walker `i` at `t_k` is sampled from the maps at the walker's *wrapped* position `r_i(t_k) mod cell_size`:

```
Δω_i(t_k) = γ · B0 · Δχ · [ Φ_0(r) + sin²θ·(cos2α·Φ_C(r) + sin2α·Φ_S(r)) ]      (reference form)
φ^χ_i     = Σ_k ε_k · Δω_i(t_k) · Δt
```

where `ε_k ∈ {+1, −1}` is the spin-echo refocusing sign, **which MUST flip at `TE/2`** (a static-field refocusing requirement; see §8.4). The susceptibility phase adds to the gradient phase inside the same exponential (`φ_i → φ_i + φ^χ_i`), so their covariance — the diffusion×susceptibility cross-term — is preserved. Field strength and orientation are replay knobs; one walk serves any `(B0, θ)`.

### 6.5 Magnetization transfer / exchange (Exchange tier)

With `bound_fraction` `b_i(t_k)`, replay blends per-step relaxation and off-resonance toward a bound-pool set `(T2_b, T1_b, Δω_b)` by occupancy, within a vector-Bloch replay (the scalar log-weight engine cannot represent an `M_z` reservoir). Saturation transfer is *emergent* (RF rotates bound spins). The bound-pool parameters are replay knobs. See the reference implementation for the Bloch–McConnell blend.

---

## 7. Capability tiers

A producer that models only part of the physics still emits a **fully interoperable** pack; it simply declares fewer tiers. Tiers are **independent** — a pack MAY declare Field without Relaxation.

| Tier | Required channels (beyond `positions`) | Replay unlocked | Metadata flag |
|---|---|---|---|
| **T0 Gradient** | *(none)* | §6.1 — any `G(t)`, b-tensor, EAP | `gradient: true` (always) |
| **T1 Relaxation** | `compartment` (+ `per_comp.T2`,`per_comp.T1`) | §6.2 — any `T2`/`T1` | `T1T2: true` |
| **T2 Surface** | `boundary_local_time` | §6.3 — any `ρ` | `rho: true` |
| **T3 Field** | `susc_field_{C,S,0}` (+ `cell_size`,`delta_chi_a`) | §6.4 — any `B0`, orientation | `B0_any: true`, `orientation_any: true` |
| **T4 Exchange** | `bound_fraction` | §6.5 — MT/exchange | `mt: true` |

The declared tier set lives in `replay_envelope` (§10). A replayer asked for a knob outside the declared tiers MUST refuse with a clear "capability not present" error, **not** silently return an approximate result (§11, §13).

*Note on `permeability`.* Membrane permeability changes the trajectory (a crossing depends on a per-event random draw), so it is **not** a pure replay knob: a pack is walked at one permeability, recorded as fixed in metadata (`replay_envelope.permeability: false` meaning "not a free knob"). Varying permeability requires new walks.

---

## 8. Generator invariants (what a producer MUST guarantee)

A pack is only replayable if the producer respected the model. A conformant producer:

### 8.1 Determinism
MUST record the `seed` and produce byte-reproducible trajectories from `(geometry, D, seed)` on a fixed implementation. The seed is metadata, not a knob.

### 8.2 Continuous positions
MUST store **unwrapped** lab-frame positions. If the walk used a periodic cell, the producer MUST unwrap before storing (undo per-save jumps `> L/2` on periodic axes; leave non-periodic axes untouched), so that the gradient phase in §6.1 is correct. Storing wrapped positions is a conformance failure: each wrap injects a spurious `q·L` phase.

### 8.3 TE-prefix consistency
MUST ensure the save grid is uniform and that any leading prefix is itself a converged walk to that shorter time (§3). Producers MUST NOT reorder or subsample walkers across the save axis.

### 8.4 Refocusing time
For Field-tier packs, the producer MUST document the encoding so that a replayer can place the `TE/2` sign flip (§6.4); the reference convention is a single 180° at `TE/2`.

### 8.5 Compartment convention
MUST use id `0` for the extra-cellular/free pool and MUST list every other id with its `(T2, T1[, R])` in `per_comp`.

### 8.6 Units and frame
MUST honor §4 exactly. A producer using non-SI internal units MUST convert on write.

### 8.7 Honest envelope
MUST declare a `replay_envelope` it can defend, and, for any lossy codec, MUST attach a `fidelity` self-report (§9.3, §11). A producer MUST NOT declare a tier whose channel it did not actually populate.

---

## 9. Codecs and compression (orthogonal layer)

**The format does not depend on compression.** Every channel is defined by its *decoded* array (§5). How that array is *stored* is a per-channel **codec** choice recorded in metadata `compression`. The identity codec — the raw array — is always valid, and a pack stored entirely with identity codecs is fully conformant.

### 9.1 Codec descriptor
Each compressed channel carries a codec name and parameters in `compression` (e.g. `{"method": "lowrank", "K": 64, "walker_preserving": true}`). A replayer MUST decode to the array contract of §5 before applying §6. Unknown codecs → the replayer MUST refuse (not guess).

### 9.2 Codec classes (informative)
Reference codecs in the dmipy-sim implementation:
- **identity** — raw positions (default; lossless).
- **temporal_dct** — temporal band-limit; walker-preserving; lossless in-band.
- **lowrank** — KL/SVD modes + exact coefficients; walker-preserving; the lossless-to-floor tier.
- **gaussian / marginal** — distributional (resample walkers); *not* walker-preserving, therefore **Gradient/EAP tier only** (they break per-walker channel alignment and MUST set the tier flags accordingly).

Only **walker-preserving** codecs may accompany the per-walker channels (`compartment`, `boundary_local_time`, `bound_fraction`, `spin_weights`); a distributional codec forces a Gradient-only pack.

### 9.3 Fidelity self-report (REQUIRED for lossy codecs)
A pack using any lossy codec MUST include a `fidelity` object reporting the maximum replay error of decoded-vs-raw channels over a declared **acquisition battery**, measured against the **split-half MC noise floor** (the ensemble split into two random halves — a *random* permutation, because walkers are typically seeded in compartment order). The convention: a codec is "lossless to the noise floor" when `err_max ≤ tol · floor_max` (reference `tol = 2`). Lossless (identity/dct/lowrank-at-full-rank) codecs SHOULD report `err_max = 0`.

Compression is thus a *quality claim about a channel*, never a change to the replay contract.

### 9.4 Reserved storage keys

A codec MAY replace a channel's single decoded array with several **reserved sub-arrays** under fixed tensor keys. **Channel presence** (for tiers, §7, and conformance, §13) is satisfied by the raw key **or** the codec's sub-arrays. The reference codecs use:

| Decoded channel | Codec | Stored tensor keys |
|---|---|---|
| `positions` | identity | `positions` |
| `positions` | `temporal_dct` | `dct_coeffs` |
| `positions` | `lowrank` (walker-preserving) | `modes`, `mean`, `coeffs` |
| `positions` | `gaussian` (distributional) | `modes`, `mean`, `coeff_mean`, `coeff_cov` |
| `positions` | `marginal` (distributional) | `modes`, `mean`, `coeff_quantiles` |
| `compartment` | RLE (row run-length) | `comp_rle_vals`, `comp_rle_lens`, `comp_rle_counts` |
| `spin_weights` | identity | `spin_weights` |
| `boundary_local_time` | identity | `dlog_boundary_unit` |
| `bound_fraction` | identity | `bound_frac` |
| `susc_field_{C,S,0}` | identity | `susc_field_C`, `susc_field_S`, `susc_field_0` |

A **distributional** positions codec (`gaussian`, `marginal` — recognizable by `coeff_mean`/`coeff_cov`/`coeff_quantiles`) resamples walkers and therefore breaks per-walker alignment: such a pack **MUST NOT** carry the per-walker channels (`compartment`, `boundary_local_time`, `bound_fraction`, `spin_weights`) and **MUST** declare only the Gradient tier. New codecs MUST document their reserved keys before use.

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
    "cell_size": 5e-6,                   // m, present iff Field tier
    "delta_chi_a": -0.1e-6               // present iff Field tier
  },
  "per_comp": {                          // REQUIRED iff Relaxation tier; else nulls
    "T2": [0.08, 0.05], "T1": [1.0, 0.8], "R": null   // index = compartment id
  },
  "compression": {                       // REQUIRED
    "method": "lowrank", "K": 64, "walker_preserving": true
  },
  "replay_envelope": {                   // REQUIRED — the declared capability + domain
    "gradient": true, "rf": true,
    "T1T2": true, "rho": false, "mt": false,
    "B0_any": true, "orientation_any": true,
    "permeability": false, "diffusivity_fixed": true,
    "acquisition": { "b_max": 3.0e9, "ogse_periods": [], "B0_list": [3.0, 7.0] }
  },
  "fidelity": {                          // REQUIRED iff any lossy codec
    "err_max": 0.0015, "floor_max": 0.0074, "within_2x_floor": true,
    "battery": "..." },
  "provenance": {                        // RECOMMENDED — how the walk was made
    "generator": "dmipy-sim", "generator_version": "...",
    "geometry": "PackedMyelinatedCylinders", "real_or_synthetic": "synthetic" },
  "license": "CC-BY-4.0",                // REQUIRED — the SOURCE geometry's license
  "citation": "Fick RHJ, dmrai-lab (2026) ..."   // REQUIRED
}
```

Rules: `diffusivity` and `seed` are fixed pack properties. `license` records the **source substrate's** license and MUST NOT relicense upstream geometry. Any tier flag set `true` in `replay_envelope` MUST have its channels present and its `per_comp`/`walk_params` fields populated.

---

## 11. Declared envelope and self-certification

A pack is a **certificate**, not a black box. Two objects make its guarantees explicit:

- **`replay_envelope`** — the *domain of validity*: which tiers, and the acquisition/field range (`b_max`, OGSE periods, `B0_list`, …) over which the producer stands behind the pack. A request inside the envelope is guaranteed; a request outside is a *known limit*, and the replayer MUST surface it as such (§13), never silently extrapolate.
- **`fidelity`** — the *measured* error of any lossy encoding against the noise floor (§9.3).

The substrate bank (a separate service) accepts a pack by validating it against this specification (§13) and MAY re-measure `fidelity` independently.

---

## 12. Container format

- **Arrays:** a single **safetensors** file. Safetensors is REQUIRED because it (a) carries **no executable code** — a pack is safe to download and memory-map from untrusted sources — (b) is strongly typed, and (c) is zero-copy. Channel names are the safetensors tensor keys; codec-produced sub-arrays use the reserved suffixes documented by the codec (e.g. `positions` may be replaced by `pos_modes`, `pos_coeff` for `lowrank`).
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
