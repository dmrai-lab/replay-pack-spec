# The Replay Pack Specification (`.rpk`)

**Version:** 0.2.3 (draft for comment)
**Status:** Draft, pre-publication — **everything stays below `1.0` until the format is published**, and the container line is `0.x` for the same reason. Under semantic versioning a `0.x` minor bump MAY break, which is the honest description of where this is: the previous `1.x` numbering promised that "field names and metadata keys are frozen", a promise `0.3` breaks by removing a channel. Stable enough to implement against, but pin a version and expect to re-encode. `0.3.0` is the first **breaking** revision: the `bound_fraction` channel is removed and magnetization transfer's bound pool becomes a `bound` **occupancy column** of `compartment` (§5.2, §5.3, §6.5) — binding is an occupancy, so it needed no channel of its own, and C4 remains a tier because a tier is defined by replay capability, not storage. The equilibrium-start requirement moves with it, anchored to the column rather than to the tier flag (§8.8), since bound occupancy biases the scalar replay on its own. A pack carrying the retired channel MUST be refused, not reinterpreted. `0.2.3` separates **what a channel contains** from **how it is stored** where the two had blurred: it pins the previously-unstated component order and degrees of freedom of `susc_field_basis` (§5.3), adds interface rule 5 — a codec that narrows *reach* rather than accuracy declares it in `replay_envelope`, not `fidelity` (§9.4, §6.4.1) — and with it two OPTIONAL envelope keys, `acquisition.max_refocusing_pulses` and `acquisition.field_modes`. No new required fields. `0.2.2` is clarification-only (no new required fields): it states **time-additivity / resume** — a walk MAY be extended in time by resuming from its final state, the forward complement of the TE-prefix property (§3) — and RECOMMENDS **archiving the geometry of non-reproducible generated substrates** so that resume and walker-shards remain possible (§8.9). `0.2.1` added the OPTIONAL orientation frame `walk_params.substrate_frame` (§4.2, §10). `0.2.0` added two OPTIONAL *representations* (a per-walker Field basis, §6.4.1; a parametric two-pool MT model, §6.5.1) and stated the additive-shard property (§3). Packs written against `0.1.x`/`0.2.x` are **not** conformant to `0.3` and must be re-encoded — that is what a `0.x` minor bump is allowed to mean, and the retired keys are refused by name so a stale pack fails loudly rather than decoding to plausible wrong numbers.
**Container schema described:** `rpk_schema_version = "0.3"`
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
6. **representations and codecs** — how `positions` are stored (raw samples or a linear *representation*) and how channels are encoded, both orthogonal to a channel's meaning (§9);
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
- **Replay knob.** A physical parameter applied *at replay time* — gradient waveform `G(t)`, static field `B0`, substrate orientation, `T2`/`T1`, surface relaxivity `ρ`, susceptibility, and the bound-pool relaxation/off-resonance of magnetization transfer. Knobs are **not** stored per value; they parameterize the replay operation.
- **Replay operation.** A deterministic function mapping (channels, knobs) → signal (§6).
- **Capability tier.** A named group of channels that unlocks one class of replay (§7).
- **Declared envelope.** The region of knob space over which the producer certifies the pack's replay to a stated fidelity (§11).
- **Representation.** The *form* `positions` are stored in — raw samples (the *source*) or a linear basis (an intermediate form, cf. a compiler IR) chosen so replay is cheap on it (§9.1). It changes replay cost and the exactness *domain* (§9.2), never a pack's meaning.
- **Codec.** A per-channel, reversible-or-lossy **storage** encoding of a channel's array (§9.4). The identity codec (raw array) is always valid.
- **MC noise floor.** The irreducible variance of a finite-`N_w` ensemble estimate; the reference against which any lossy codec's error is judged (§9.4).
- **Shard.** One pack among several that walk the *same* substrate/physics with *different* seeds; their walker ensembles combine into a larger effective ensemble (§3, additive property).

---

## 3. The replay model (the invariant)

Replay rests on one physical fact.

> **Replay invariant.** For a fixed substrate geometry, a fixed (possibly per-compartment) free diffusivity `D`, and a fixed pseudo-random seed, the ensemble of walker trajectories `{r_i(t_k)}` and the boundary-contact record are **independent of** the gradient waveform `G(t)`, the static field `B0`, the substrate orientation, the relaxation times `T2`/`T1`, the surface relaxivity `ρ`, the susceptibility distribution, and the bound-pool relaxation/off-resonance of magnetization transfer.

All of those quantities enter the *signal* only through a phase or a (log-)weight accumulated **along** an already-determined trajectory. Therefore they are **replay knobs**: the trajectory is computed once and stored; the signal for any knob setting is a cheap functional of the stored channels (§6).

*Magnetization transfer is split.* The MT **binding kinetics** (binding rate and dwell time) are **not** knobs: a bound spin is frozen for its dwell, so they change the trajectory itself and are fixed at walk time — recorded through the C1 `bound` occupancy column (§5.2), exactly as membrane permeability (§7) is a walk-time property. Only the **bound-pool relaxation/off-resonance** (`T2_b`, `T1_b`, `Δω_b`) are replay knobs (§6.5). A pack MAY instead carry MT in a coarser but cheaper **parametric two-pool** representation (§6.5.1) whose pool size and exchange rate are derived from the substrate's surface-to-volume ratio — no per-walker binding channel, and no binding walk — with the RF saturation as the replay knob; the two representations are the emergent and mean-field forms of the same C4 tier.

Three consequences are used throughout:

- **(TE-prefix property.)** A walk stored to `T_max` contains every shorter echo time as a leading prefix: the first `N_t' = ⌈TE'/Δt⌉` saves are a valid, converged walk to `TE' ≤ T_max`. A replayer MAY slice.
- **(Separation of state and knobs.)** A pack stores only *walk state*. It never stores a signal, a b-value, or a waveform. Signals are always derived at replay time.
- **(Resolution limit.)** The save step `Δt` sets the temporal bandwidth of every replay: a pack reproduces any `G(t)` resolvable on the `Δt` grid (Nyquist, up to `1/(2Δt)`); structure finer than `Δt` was never stored and cannot be replayed. This is a property of the stored walk, independent of how positions are represented (§9).
- **(Additivity — walkers.)** Because walkers are independent realizations, two packs that share the *same* substrate geometry, diffusivity, save grid `Δt`, and channel set but use **different seeds** are **shards** of one larger ensemble: a replayer MAY pool them by concatenating their per-walker channels (or averaging their signals weighted by `N_w`), and the pooled MC noise floor falls as `1/√(Σ N_w)`. A bank therefore need not commit to a final `N_w` up front — it MAY publish a first shard and append more later to drive the floor down, with no re-walk and no change to any published shard. Shards MUST agree on the fixed-at-walk-time properties below; only the seed (and hence the walkers) differs.
- **(Additivity — time; resume.)** Because the walk is **Markov** — a walker's position is its full state — a walk to `T_max` MAY be **extended in time** to any `T_max' > T_max` by resuming each walker from its `T_max` position, continuing the process for the additional `T_max' − T_max`, and concatenating the per-walker channels along the time axis. The `[0, T_max]` segment is preserved **exactly**; the extension is a valid (fresh-seed) continuation of the same process. This is the forward complement of the TE-prefix property: a bank MAY publish a **short** walk and lengthen it later — with no re-walk of the existing segment — to reach longer echo/diffusion times. Two conditions: (i) the extension MUST use the identical fixed-at-walk-time properties below (geometry, diffusivity, `Δt`, per-compartment `D`, permeability, MT kinetics); and (ii) any per-walker walk-time state that persists across the seam — notably the MT bound/free occupancy (§6.5) — MUST be carried from the `T_max` state (a parametric-MT or no-MT pack has none, so it resumes freely; an emergent-MT pack MUST supply it or MUST NOT be resumed). A pack is therefore extendable along **both** axes — more walkers (shards) and more time (resume) — provided the substrate geometry is available to re-walk (§8.9).

The invariant defines the boundary of the format. The following are **fixed at walk time — not replay knobs**; changing any requires a new walk and a new pack:

- the **substrate geometry**;
- the **diffusivity** `D`;†
- the **random seed**;
- **membrane permeability / crossing** — a crossing is a per-event random draw that alters the trajectory (§7);
- the **MT binding kinetics** (rate, dwell) — recorded in the C1 `bound` column (§6.5);
- the **save step** `Δt` — the temporal resolution of the stored walk (resolution limit, above).

† `D` sets the trajectory, so it is recorded as a fixed pack property (§10). The sole exception is the trivial Brownian rescaling — a walk at `D` over `t` equals one at `D/α` over `αt` — which holds **only** in free diffusion, or when the geometry length-scales *and* the save grid `Δt` are co-scaled by the same `α`. In restricted geometry the fixed wall scales break it, so `D` is **not** a general replay knob.

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
| bound occupancy, spin weight | `b_i`, `w_i` | dimensionless | `compartment` (`bound` column), `spin_weights` |
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
- **Substrate frame (intrinsic orientation).** Because orientation is applied by rotating the acquisition, a replayer needs to know how the substrate is oriented *within* the stored lab frame — otherwise "parallel" vs "perpendicular" B₀ (very different susceptibility dephasing) and each diffusion-gradient direction are undefined relative to the tissue. A pack SHOULD therefore declare `walk_params.substrate_frame` (§10): a **3×3 orthonormal, right-handed matrix `R`** whose columns are a canonical substrate basis expressed in stored coordinates, with **column 3 (`z`) the primary structural axis** (e.g. the fibre-bundle direction) and columns 1–2 a **fixed perpendicular basis**. A single axis is insufficient: it leaves a free rotation about itself, so the perpendicular plane — and hence any non-axisymmetric structure, multi-direction scheme, or b-tensor — would be gauge-dependent from run to run. The producer MUST pin `x`/`y` deterministically, aligning them to genuine in-plane structure where one exists. For a **single, axisymmetric** population, `x`/`y` are a fixed but arbitrary basis (e.g. Gram–Schmidt from a fixed reference). For a **multi-population (e.g. crossing)** substrate the producer MUST anchor the frame to a **primary** population, never to a principal-axis/PCA average over all positions — that average is the bundles' *bisector*, which is no real population's axis and rotates as the crossing angle opens, leaving a crossing-angle series with no shared reference. The RECOMMENDED convention: `z` = the primary population's mean axis (chosen by a deterministic rule declared in provenance — e.g. largest by volume, tie-broken by index, independent of packing RNG); `y` = the secondary population's mean axis projected perpendicular to `z` (so the crossing opens into `+y`); `x = y × z`. The primary population then sits at `z` in every pack of a crossing-angle series and the secondary rotates `z → +y` (one bundle held straight, the other rotating), so replays are directly comparable. Two populations fully pin the frame; further populations are merely described by their axes within it. A direction given in the **canonical frame** (`z` along the primary axis) maps to stored coordinates as `d_stored = R · d_canonical`; a replayer orients an acquisition onto the pack by this rotation. Every pack SHOULD carry `substrate_frame` for uniform behaviour, **including isotropic substrates**, for which `R` MAY be the identity but MUST still be fixed. A pack that omits it MUST be treated as `R = I` (positions already in the canonical frame); producers of oriented substrates are strongly encouraged not to rely on this default.

### 4.3 Time
The save grid is uniform with interval `Δt = dt_traj`. `N_t` and `T_max = (N_t − 1)·Δt` are recorded. The integrator's internal sub-step is a producer concern and is neither stored nor standardized.

### 4.4 Signal and sign convention
The canonical signal is complex, `S = ⟨w · e^{+iφ}⟩`, with phase `φ = γ ∫ G·r dt` (§6.1), a **positive** gyromagnetic sign (`γ > 0`), and the `e^{+iφ}` rotation sense. Magnitude DWI is insensitive to this choice, but it is fixed here so the **complex** signal and the susceptibility phase (§6.4) are reproducible across implementations. A real spin-echo value is `|S|` or `⟨w·cos φ⟩`. Averages `⟨·⟩` are weighted means over walkers (weights `w_i`, default `1`).

### 4.5 Grids
Shared Field **grid maps** (§5.2, §6.4) are sampled at the walker's position on a uniform grid whose axes are in metres and whose spacing is `walk_params.cell_size`. For a **periodic** substrate the position is **wrapped** into the field cell, `r mod cell_size`; for a **finite, non-periodic** substrate it is **clamped** to the grid bounds. The producer declares which via `walk_params.boundary` (§10); wrapped is the default for back-compatibility. (The per-walker Field **basis** representation of §6.4.1 is sampled at walk time and carries no replay-time grid, so this convention does not apply to it.)

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
| `compartment` | `(N_w, N_t)` per column | int16 / float32 | label or ∈ [0,1] | **Occupancy over the declared pools**, as one or more named columns. The `comp` column is the exclusive geometric axis: a compartment id of walker `i` at save `k`, id **0 the extra-cellular / free compartment by convention**, positive ids producer-defined and described in `per_comp` (or a fraction where a permeable crossing splits a save). Further columns are occupancies on **independent** axes — `bound` for the MT macromolecular pool (§6.5). | Relaxation (`comp`); Magnetization transfer (`bound`) |
| `boundary_local_time` | `(N_w, N_t)` | float32 | m (see below) | Per-step accumulated wall-contact measure with **`ρ/D = 1`** (dimension of length ÷ diffusivity is folded so that replay multiplies by the desired `ρ/D`). See §6.3. | Surface |
| `susc_field_0` | producer grid | float32 | per unit susceptibility | m=0 (isotropic / mean) component of the substrate's normalized off-resonance field map. | Field |
| `susc_field_C` | producer grid | float32 | " | ℓ=2 cosine (anisotropic) component. | Field |
| `susc_field_S` | producer grid | float32 | " | ℓ=2 sine (anisotropic) component. | Field |
| `susc_field_basis` | `(N_w, N_t, 13)` | float32 | per unit susceptibility | Per-walker off-resonance field basis **sampled along the walk** — the *general* Field representation (§6.4.1). The 13 components are `iso_local` (1), `iso_P` (6), and `aniso_G` (6): the isotropic-source and rank-2 anisotropic-source responses, projected onto the 6-vector `Q(H)` of the main-field direction at replay. The component axis is **last**, matching `positions`, and its order is fixed by §5.3; the 13 components carry 12 degrees of freedom (§5.3). Per-walker ⇒ needs a walker-preserving codec (§9.4). | Field |

The Field tier has **two representations of the same off-resonance physics** (§6.4), a producer choosing one:
- **Shared grid maps** `susc_field_{0,C,S}` — **normalized off-resonance basis maps** on a grid, sampled at replay (§6.4). `susc_field_0` carries the isotropic part and `susc_field_{C,S}` the axially-anisotropic (ℓ=2) part. Compact and reusable, but assumes a field **translationally invariant along the fibre axis** (the straight-cylinder idealization); best for periodic/analytic substrates. Accompanied by `walk_params.cell_size` (grid spacing, m), the normalization scale(s) (reference: `walk_params.delta_chi_a`), and, where oriented, `per_comp.R`.
- **Per-walker basis** `susc_field_basis` — the field basis **sampled along each trajectory** (§6.4.1), the analog of `boundary_local_time` for the field tier. It carries no invariance assumption and so serves **arbitrary 3-D morphology** (undulating/dispersed meshes) and **any** main-field direction via the full rank-2 `Q(H)` contraction; the cost is per-walker storage. Accompanied by `walk_params.delta_chi_a` and the field-grid metadata used to build it (`walk_params.field_grid`, §10).

A fully general per-voxel susceptibility *tensor* field remains a §14 extension; the `susc_field_basis` set covers an isotropic source plus a single rank-2 anisotropic director (the myelin case).

A producer **MAY** define **additional** channels prefixed `x_` (e.g. `x_temperature`). Replayers **MUST** ignore unknown channels they do not implement (§14).

### 5.3 Channel invariants

- Every per-walker/per-save channel MUST share `N_w` and `N_t` with `positions`.
- `compartment` MUST carry the exclusive `comp` column; its id `0` MUST denote the extra-cellular/free pool and other ids MUST be described in `per_comp`. Every column MUST share `N_t`.
- At each save the occupancy over all declared pools sums to 1. Columns are the **compact encoding of that joint simplex**, not independent quantities: with a bound occupancy `b` and geometric label `g`, the pool occupancies are `f_bound = b` and `f_g = (1 − b)·[comp = g]`. Replay weights the per-pool rates by them, `R(t) = Σ_c f_c(t) R_c` — which is what makes an integer label (one-hot), a permeable crossing (a fraction on the geometric axis), and MT binding (a fraction on an independent axis) the same object.
- The `bound` column MUST lie in `[0, 1]`, and at `k = 0` it MUST be the equilibrium occupancy (the walk is pre-burned-in; §8.8).
- `boundary_local_time` MUST be non-negative and expressed in the `ρ/D = 1` normalization of §6.3.
- `susc_field_basis` MUST order its component axis (the **last** axis) as `iso_local`, then `iso_P` as `(xx, yy, zz, xy, xz, yz)`, then `aniso_G` in the same six-component order — so that slices `[..., 1:7]` and `[..., 7:13]` contract directly against `Q(H)` as written in §6.4.1. A substrate with no anisotropic source stores an all-zero `aniso_G` block.
- `susc_field_basis` carries **13 components but 12 degrees of freedom**: the isotropic block obeys the exact trace identity `iso_P_xx + iso_P_yy + iso_P_zz = 3 · iso_local`, so `iso_P_zz` is recoverable from the other three. This holds to machine precision when the isotropic source grid the basis was built from is un-apodised, and is **broken** by a k-space window applied to `iso_P` but not to `iso_local`. A producer MUST verify the identity rather than assume it; a codec MAY exploit it (dropping and reconstructing `iso_P_zz`) only where it verifies.

---

## 6. Replay operations (normative reference semantics)

These define the meaning of the channels. A **conformant replayer** MUST reproduce these operations (to the decoded arrays' numerical precision). `γ` is the proton gyromagnetic ratio. Sums over `k` run on the save grid with weight `Δt` (trapezoidal or left-rule is producer-agnostic *provided the producer used the same rule to converge the walk*; left-rule is the reference).

The operations split by their dependence on position: the gradient phase (§6.1) and the relaxation/surface log-weights (§6.2–6.3) are **linear** in the stored positions, whereas the off-resonance field-map lookup (§6.4) and the Bloch/MT evolution (§6.5) are **nonlinear**. This is informative only — the signal is identical either way — but it lets a replayer evaluate the position-linear operations directly in a linear *representation*'s coefficient space without decoding the full trajectory (§9.1; exactness domain §9.2).

### 6.1 Gradient phase (Gradient tier — always available)

For a gradient waveform `G(t_k)` (T·m⁻¹, resampled to the save grid) the accumulated phase of walker `i` is

```
φ_i = γ · Δt · Σ_k G(t_k) · r_i(t_k)
```

and the signal is the weighted ensemble mean

```
S(G) = Σ_i w_i · e^{i φ_i} / Σ_i w_i .
```

This yields any b-value, any b-tensor, OGSE/PGSE/arbitrary `G(t)` **resolvable on the `Δt` grid** (Nyquist; the resolution limit of §3), and (by sweeping `q`) the ensemble average propagator. **No stored quantity depends on `G`** — this is the whole point.

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

#### 6.4.1 Per-walker field basis (general representation)

When the field is **not** translationally invariant — an undulating, dispersed, or otherwise irregular 3-D substrate — the ℓ=2 grid maps do not apply, and the producer instead stores the field basis **sampled along each trajectory** in `susc_field_basis` (§5.2). Its 13 components per `(i, k)` are the geometry-only responses `iso_local`, `iso_P[6]`, `aniso_G[6]`. At replay, for a main-field **unit direction** `H = (H_x,H_y,H_z)`, form the symmetric-tensor contraction 6-vector

```
Q(H) = (H_x², H_y², H_z², 2 H_x H_y, 2 H_x H_z, 2 H_y H_z)
```

and the off-resonance seen by walker `i` is

```
Δω_i(t_k) = γ · B0 · [ χ_iso · ( iso_local − Q·iso_P )  +  χ_aniso · ( Q·aniso_G ) ]_{i,k}
```

which feeds the **same** transverse-phase gate and phase accumulation as §6.4 (`φ^χ_i = Σ_k s(t_k)·Δω_i(t_k)·Δt`), preserving the cross-term identically. This is the *general* form: the grid-map operation (§6.4) is its special case for a field invariant along a fixed fibre axis, and `Q(H)` supersedes the `sin²θ·(cos2α,sin2α)` factor with the full main-field direction — so **any** B0 orientation is a replay knob, not only the polar sweep. Because the basis is pre-sampled, this representation decodes and contracts without a field grid at replay; the trade is per-walker storage versus the grid's shared footprint.

**A truncated temporal storage of this channel costs reach, not accuracy.** The phase reaches the field only through the gate, `φ^χ_i = Σ_k s(t_k)·Δω_i(t_k)·Δt`. A codec that stores the channel as `K` temporal modes (registry) is therefore band-limited **in the gate**, not in the field: the retained modes reproduce `φ^χ` for any `s(t)` inside their band while reproducing `Δω_i(t_k)` itself poorly, because the field along a trajectory is broadband — the walker crosses μm-scale field structure every save step. `K` is thus a **capability**: it bounds the refocusing-pulse density the pack can serve (reference implementation: `K ≥ 2·n_refocus`), and MUST be declared in `replay_envelope` (§9.4 rule 5, §11), never read as "the field is accurate to `K` modes".

### 6.5 Magnetization transfer (Magnetization-transfer tier)

With the C1 `bound` occupancy column `b_i(t_k)` (§5.2), replay blends per-step relaxation and off-resonance toward a bound-pool set `(T2_b, T1_b, Δω_b)` by occupancy, within a **vector-Bloch replay** — RF pulses act as rotations of the magnetization vector, which the scalar log-weight model cannot represent (it has no `M_z` reservoir). Saturation transfer is *emergent*. The bound-pool parameters are replay knobs. This same vector-Bloch path is the general route for *arbitrary RF* (§6.6). See the reference implementation for the Bloch–McConnell blend.

**MT stores nothing of its own.** The bound pool is a column of the C1 occupancy channel, because binding *is* an occupancy: a spin bound for a fraction of a save has that fraction of the step at the bound pool's relaxation, which is exactly the occupancy weighting C1 already defines. What makes magnetization transfer a distinct **capability** tier is the replay side — RF as vector-Bloch rotations (the scalar log-weight path cannot express it), the bound-pool knobs, the equilibrium start below, and the certification duty of §9.2 — not the storage. Two tiers sharing a storage form is normal here: C0 and C2 share a byte-identical bridge layout (§9.3) and remain separate tiers because their knobs differ. A useful consequence of storing binding as C1 rather than as a summary parameter: the run lengths of the `bound` column **are** the dwell times, so an emergent pack carries the realised exchange statistics, not just a mean rate.

**Equilibrium start.** A magnetization-transfer pack's `t=0` is the **bound-pool equilibrium**: the walk is generated pre-burned-in (§8.8), so replay begins from a fully-relaxed steady-state occupancy and never sees the fill-up transient. *(This applies to the emergent per-walker representation; the parametric representation of §6.5.1 carries no `t=0` occupancy transient.)*

#### 6.5.1 Parametric two-pool MT (mean-field representation)

The emergent representation above requires the `bound` occupancy column, which a producer obtains from a binding walk that resolves near-wall motion at a fine sub-step. On sub-micron restricting features that sub-step makes the binding walk intractable. A producer MAY therefore declare C4 in a **mean-field two-pool** form instead: rather than a per-walker channel, it stores a small **pool descriptor** in metadata `mt` (§10) — the bound-pool fraction `f_b`, the forward exchange rate `k_f`, and the bound-pool relaxation `(T2_b, T1_b)` — derived from the substrate's myelin **surface-to-volume ratio** `S/V` (with the binding reactivity `κ` and dwell `τ`): `k_f = κ·S/V`, `f_b = k_f·τ / (1 + k_f·τ)`. These are **fixed-at-walk-time geometry/tissue properties** (like `D`), not knobs.

At replay the observable free-pool signal follows the two-pool **Bloch–McConnell** equations under an RF saturation of nutation `w1` and offset `Δ` for duration `t_sat`; a replayer offering this representation MUST expose at least the **Z-spectrum** (free-pool `M_z` vs `Δ`) and the derived **MT ratio**. RF offset, power, and duration are the replay knobs. This form is validated to reproduce the emergent equilibrium `f_b` (the mean matches `κ·S/V·τ`), but it deliberately does **not** model the binding↔diffusion coupling the emergent walk captures — a bound spin's frozen path — so it is the correct choice for MT saturation / qMT observables and **not** for diffusion–MT coupling. The two are the mean-field and emergent forms of the **same C4 tier**; a pack declares exactly one.

### 6.6 The acquisition is a replay knob (no sequence is baked into the pack)

A pack stores **no** assumption about the pulse sequence. The full acquisition — the gradient waveform `G(t)` **and** the RF / `B1` pulses (their flip angle, phase, shape, and **timing**) — is supplied at replay. Moving the refocusing pulse when the diffusion time or `TE` changes, or switching GRE↔spin-echo↔CPMG↔stimulated-echo/PGSTE, is a replay-time choice against the one stored walk, never a re-walk. (Spatially *uniform* `B1` pulses need nothing stored; `B1` **inhomogeneity** is a scalar flip-angle knob per pack — becoming a per-voxel scalar only when packs are tiled into a mesoscopic grid — not a within-substrate field; §14.)

For the **scalar-phase tiers** (§6.1, §6.4) the sequence enters through a per-step **transverse-phase gate** `s(t_k)` that the replayer derives from the sequence:

- `s ≡ +1` — gradient echo (GRE), no refocusing;
- `s` flips sign at each 180° refocusing pulse — spin echo, CPMG;
- `s = 0` while magnetization is stored along `z` — stimulated echo (STE/PGSTE);

so `s(t_k) ∈ {−1, 0, +1}` in the common cases. The **gradient phase (§6.1)** already encodes any bipolar/refocused *gradient* structure in the waveform's polarity and needs no gate; the gate applies to the **static off-resonance** term (§6.4), which a refocusing pulse inverts and a storage period freezes. (The transverse/longitudinal indicator `χ_k` of §6.2 is the same idea for relaxation.)

Arbitrary RF — non-180° flips, adiabatic pulses, MT saturation — cannot be reduced to a scalar sign and requires the **vector-Bloch replay** (§6.5). A conformant replayer MAY implement only the scalar-gate model; it MUST then **refuse** acquisitions that need full Bloch evolution rather than approximate them (§13).

### 6.7 Orientation dispersion (informative, non-normative — analytical overlay)

*This subsection is **informative**: it defines no new channel, tier, or metadata flag, and a pack is fully conformant without it.* It records a common, useful **derived** operation that a replayer MAY offer on top of the substrate-orientation knob (§4.2).

The orientation knob already evaluates one coherent bundle at any orientation by rotating the waveform. An orientation *distribution* is therefore the fibre-response ⊗ ODF superposition of that same bundle over orientations: for a single-fibre response `R` (the stored bundle) and an orientation distribution `p(Ω)`,

```
S_voxel(G) = ∫_{S²} p(Ω) · S(R_Ω⁻¹ · G) dΩ .
```

For an axially-symmetric fibre response this is the standard rotational-harmonic (Funk–Hecke) convolution: expand the response in zonal harmonics `R_l` and the ODF in real SH `f_lm`, then `S_lm = 2π·√(4π/(2l+1))·R_l·f_lm`. It is evaluated at replay from the **single stored bundle** — sample the acquisition at a set of polar angles, fit `R_l`, convolve — so it costs nothing extra to store.

**It is an analytical model, not simulated dispersion. The caveats are the point:**

- **No exchange.** It superposes **independent, non-exchanging** fibre populations: water never crosses between differently-oriented fascicles. Genuinely dispersed physics (crossing geometry, inter-fascicle exchange) changes the *trajectory* and requires a new walk in a dispersed substrate — it is not recoverable by replay.
- **Homogeneous copies.** Every fascicle is the *same* stored bundle rotated — one radius distribution, `D`, packing, and relaxation. Per-population heterogeneity is not represented.
- **Scope.** The zonal convolution covers the **diffusion signal** (§6.1) and separable per-compartment relaxation (§6.2–6.3). It does **not** cover off-resonance/susceptibility dispersion (§6.4 would additionally require rotating the field basis per orientation) or MT dispersion (§6.5, which compounds the no-exchange caveat); a replayer MUST refuse those rather than convolve them.
- **Convergence.** The SH order `lmax` must resolve both the response and the ODF; an under-resolved order biases the convolution.

Because it stores nothing and changes no replay semantics, it needs no envelope flag; a producer MAY note an intended dispersion in `provenance`. A replayer offering it **MUST** present the result as an analytical dispersion model, not Monte-Carlo ground truth.

---

## 7. Capability tiers

A producer that models only part of the physics still emits a **fully interoperable** pack; it simply declares fewer tiers. Tiers are **independent** — a pack MAY declare Field without Relaxation.

Tiers are named **`C0`–`C4`** (for *capability*) — deliberately not `T#`, which would collide with the relaxation times `T1`/`T2`. Metadata flags use **explicit, self-describing names** (§10).

| Tier | Required channels (beyond `positions`) | Replay unlocked | Metadata flag(s) |
|---|---|---|---|
| **C0 Gradient** | *(none)* | §6.1 — any `G(t)`, b-tensor, EAP | `gradient: true` (always) |
| **C1 Bulk relaxation** | `compartment` (+ `per_comp.T2`,`per_comp.T1`) | §6.2 — any `T2`/`T1` | `bulk_relaxation: true` |
| **C2 Surface** | `boundary_local_time` | §6.3 — any surface relaxivity | `surface_relaxivity: true` |
| **C3 Field** | `susc_field_{C,S,0}` grid maps (+ `cell_size`) **or** `susc_field_basis` (+ scale) | §6.4 / §6.4.1 — any `B0`, susceptibility, orientation | `field: true` |
| **C4 Magnetization transfer** | C1 `bound` occupancy column (emergent) **or** parametric `mt{}` pool (§6.5.1) | §6.5 / §6.5.1 — magnetization transfer | `magnetization_transfer: true` |

The declared tier set lives in `replay_envelope` (§10). A replayer asked for a knob outside the declared tiers MUST refuse with a clear "capability not present" error, **not** silently return an approximate result (§11, §13).

*On membrane permeability (out of scope in `0.x`).* Permeability changes the **trajectory** (a crossing is a per-event random draw), so it is not a replay knob at all — a pack is walked at one permeability and varying it requires new walks. The `0.x` envelope therefore does **not** carry a permeability flag; the fixed value the walk used SHOULD be recorded under `provenance`. If a future version finds a way to make crossing a replayable quantity, it will be added as a new tier (§14).

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
A producer MUST NOT bake any pulse-sequence choice into the pack. The acquisition (gradients and RF) is applied at replay (§6.6); the producer's obligation is only to store the trajectory (and any tier channels) over the full walk to `T_max`. In particular there is no "refocusing time" a producer must record — refocusing is a property of the replayed sequence, not of the walk. (A producer that chooses a capability-narrowing codec does declare a refocusing-pulse *density* it can serve — `replay_envelope.acquisition.max_refocusing_pulses`, §9.4 rule 5. That is a limit of the chosen **storage**, not a sequence assumption baked into the walk: the raw channel carries none.)

### 8.5 Compartment convention
MUST use id `0` for the extra-cellular/free pool and MUST list every other id with its `(T2, T1[, R])` in `per_comp`.

### 8.6 Units and frame
MUST honor §4 exactly. A producer using non-SI internal units MUST convert on write.

### 8.7 Honest envelope
MUST declare a `replay_envelope` it can defend, and, for any lossy codec, MUST attach a `fidelity` self-report whose battery exercises **every declared tier** (§9.4, §11). A producer MUST NOT declare a tier whose channel it did not actually populate.

### 8.8 Equilibrium start for magnetization-transfer packs
A pack carrying the C1 `bound` occupancy column MUST begin from the **bound-pool equilibrium**. A Monte-Carlo binding walk started from an arbitrary state (e.g. all spins free) shows a transient while the bound fraction relaxes to its steady state `f_b = k_f/(k_f + k_r)`; that transient does not represent a fully-relaxed sample and must not appear in a replay. The producer MUST therefore **equilibrate** the binding dynamics and **discard** that preamble, saving the walk from `t=0 =` the equilibrated state (the `bound` column accordingly starts at equilibrium). Replayers assume `t=0` is equilibrium (§6.5). The producer SHOULD verify the occupancy has reached steady state before saving.

The requirement attaches to the **column, not to the tier flag** — the column's presence is itself the assertion, and a consumer needs no C4 declaration to depend on it. The reason is that bound occupancy biases the *scalar* replay on its own: a spin bound for part of a save carries the bound pool's short `T2_b` for that fraction under the ordinary C1 occupancy weighting `R(t) = Σ_c f_c(t) R_c` (§5.3), with no RF and no `M_z` reservoir involved. A replayer using only bulk relaxation therefore reads the column too, so anchoring the guarantee to the magnetization-transfer tier would leave that consumer relying on a promise nobody had made. No other channel has an initial-condition requirement: the transient is a property of the binding kinetics, not of the walk.

### 8.9 Reproducibility of generated substrates
A substrate from a **stochastic generator that is not deterministically reproducible from its configuration** (e.g. a mesher seeded from wall-clock time or process id) cannot be recovered from the config alone. For such substrates the producer SHOULD **archive the exact geometry used for the walk** — as a hashed sidecar — and reference it in `provenance.substrate.ref` (§10). This is needed not only so the substrate is reproducible and inspectable, but because the additive operations of §3 — **resuming/extending the walk in time** and **adding walker shards** — both *re-walk the same geometry*; without an archived geometry they become impossible once the generator's working directory is gone. The geometry *format* is out of scope for the container (§1), but for interoperability the archive SHOULD use an **open, widely-readable geometry format** (e.g. PLY/OBJ/STL for meshes; NIfTI for label volumes) rather than a producer-specific binary, so any consumer — in any language — can read it. Compression is orthogonal: a plain `xz`/`gzip`/`zstd` wrapper keeps it a one-line decompress; a container SHOULD state which. A substrate that *is* bit-reproducible from its config (config + version + seed) MAY instead reference that recipe.

*Reference implementation.* dmipy-sim exposes this as `equilibrate_binding = 'auto' | 'burnin' | 'fast' | 'off'`: `'burnin'` (the `'auto'` default when MT is on) runs an **adaptive RF-off burn-in in ~dwell-sized chunks until the ensemble occupancy plateaus** (geometry-agnostic, with a convergence flag); `'fast'` seeds the equilibrium occupancy analytically from a known `S/V` when that is position-invariant; `'off'` keeps the legacy all-free start. A producer MAY use any method that satisfies the equilibrium-start requirement above.

---

## 9. Representations and codecs (orthogonal to meaning)

**A channel's *meaning* is its decoded array (§5); how it is stored is orthogonal to that meaning.** A pack stored entirely as raw arrays (the `identity` baseline) is fully conformant. Two independent storage choices exist, and neither changes what a pack *means* — only its size and its replay *cost*:

- a **representation** of `positions` — the *form* the trajectory is stored in (§9.1);
- a per-channel **storage codec** — a reversible-or-lossy encoding of any channel's array (§9.4).

Both are declared in metadata `compression`, both MUST decode to the §5 contract, and the concrete methods with their exact stored-key layouts live in the separately versioned [`CODEC_REGISTRY.md`](CODEC_REGISTRY.md), so they evolve without touching this frozen core.

### 9.1 Representation of positions — source vs. intermediate form

Positions may be stored as **raw samples** — the *source* — or projected onto a **linear basis** `r = M·c + μ` (a shared basis `M`, per-walker coefficients `c`). This is **not mere compression**: it is a change of representation, closer to a compiler's **intermediate representation (IR)** than to a zip file. The raw samples are the source; the basis is the form chosen because the operation that runs on it — replay — is *cheaper* on it. Because the gradient phase (§6.1) and the separable relaxation/surface log-weights (§6.2–6.3) are **linear** in position, replay in the basis is a projection onto `K` coefficients instead of an `O(N_t)` integral per walker, evaluated without ever reconstructing the trajectory (§6 preamble; mechanics in the registry). The nonlinear operations — susceptibility (§6.4) and Bloch/MT (§6.5) — decode to raw first.

### 9.2 Exactness and its domain

A representation's exactness has a **domain** — the single sharpest point to state plainly:

- **Raw samples** are exact for any `G(t)` resolvable on the `Δt` grid (the resolution limit of §3).
- **A linear basis** is exact for any `G(t)` whose induced phase lies in the **span** of its `K` modes. A full, in-band basis spans everything `Δt` can resolve and therefore **coincides with raw** — exact wherever the source is. A **truncated** basis (`K` below full rank) is exact only **within its span**; outside it the error is bounded and MUST be certified by `fidelity` (§9.4).

So raw is exact to the grid, a full basis matches raw, and a truncated basis trades exactness outside its span for size and replay speed — and declares that trade in `fidelity`.

### 9.3 Storage codecs

The non-position channels (`compartment`, `boundary_local_time`, …) use per-channel **storage codecs** — run-length, quantization, sparse/dense — matched to each channel's structure. These are pure storage: they decode to the §5 array and do not change replay cost.

### 9.4 Interface rules

A conformant representation or codec MUST satisfy five rules:

1. **Declared.** The codec is named (with any parameters) in metadata `compression`, e.g. `{"method": "identity"}` or `{"method": "lowrank", "K": 64, "walker_preserving": true}`. A replayer that does not recognize the `method` MUST refuse (never guess). Recognized method names and their stored tensor keys are defined in the registry.
2. **Decodes to the contract.** Decoding MUST yield the exact channel array of §5 (dtype/shape/units). Channel *presence* for tiers (§7) and conformance (§13) is satisfied by the channel's raw key **or** the registry-defined stored keys of its codec.
3. **Identity is the baseline.** `identity` (the raw array under its channel name) is always valid and lossless, and every replayer MUST support it.
4. **Lossy ⇒ self-certified.** A pack using any lossy codec MUST include a `fidelity` object reporting the maximum decoded-vs-raw replay error over a declared **acquisition battery**, measured against the **split-half Monte-Carlo noise floor** (the ensemble split into two *random* halves — random because walkers are typically seeded in compartment order). The convention: "lossless to the noise floor" means `err_max ≤ tol · floor_max` (reference `tol = 2`). Lossless codecs SHOULD report `err_max = 0`. The battery MUST exercise **every tier the pack declares** — in particular the *nonlinear* Field (C3) and magnetization-transfer (C4) operations: a truncated linear representation's error on those is generally **larger than, and uncorrelated with,** its error on the gradient phase (§9.2), so a gradient-only battery does not certify them.
5. **Capability-narrowing ⇒ declared in the envelope, not in `fidelity`.** Most codecs cost only accuracy, and rule 4 is the whole story. A truncated representation of a channel that replay consumes **non-linearly** can instead cost *reach* — leaving some requests unservable at any accuracy — and `fidelity` cannot express that: an `err_max` measured over a battery that never leaves the retained span reports a small number for a request the pack cannot serve at all. Such a codec MUST declare the resulting limit in `replay_envelope` (§11) as well as its `fidelity`. The registered case is the per-walker Field basis stored as `K` temporal modes, which bounds the servable gate bandwidth rather than the field accuracy (§6.4.1) and declares `replay_envelope.acquisition.max_refocusing_pulses`.

The `fidelity` object is **codec-agnostic** — it reports the decoded-vs-raw replay error whatever method produced it — so it lives in the **core** metadata (§10) and is read identically by any replayer or bank; only the codec *algorithms* live in the registry. How to read it: `err_max` is the worst-case replay error over the battery, `floor_max` the irreducible MC noise floor for the same ensemble, and `within_2x_floor: true` means the loss sits below `2×` that floor (scientifically negligible).

One structural rule belongs in the core because it constrains tiers, not any particular algorithm:

> **Walker-preserving requirement.** A codec that **resamples walkers** (a *distributional* codec — it stores a coefficient *distribution*, not per-walker coefficients) breaks per-walker channel alignment. Such a pack **MUST NOT** carry any per-walker channel (`compartment`, `boundary_local_time`, `spin_weights`, `susc_field_basis`) and **MUST** declare only the Gradient tier. Per-walker channels REQUIRE a walker-preserving codec.

A channel's storage never changes what a pack *means* — the decoded array of §5 and the operations of §6 are fixed. For channels replay consumes **linearly** it is a pure *size/speed/quality choice*, and §9.2's span argument bounds the trade. For a channel consumed **non-linearly** — the Field basis (§6.4.1), Bloch/MT (§6.5) — a truncated representation may additionally narrow the *domain of validity*, which is a matter for `replay_envelope` (rule 5), never a silent property of the stored bytes. See [`CODEC_REGISTRY.md`](CODEC_REGISTRY.md) for the registered representations and storage codecs and their stored-key tables.

---

## 10. Metadata schema

Metadata is a JSON object embedded in the container (§12) and validated by `schema/rpk_metadata.schema.json`. Required and optional keys:

```jsonc
{
  "rpk_schema_version": "0.3",          // REQUIRED, semver of the container schema
  "id": "canonical-wm/g070-f055-3T",    // REQUIRED, stable pack identifier
  "walk_params": {                       // REQUIRED
    "n_walkers": 120000, "n_t": 200,
    "dt_traj": 2.65e-4, "T_max": 0.053,  // s
    "diffusivity": 0.6e-9,               // m^2/s, FIXED (not a knob)
    "seed": 0,
    "boundary": "periodic",              // OPTIONAL "periodic" | "reflecting" | "finite";
                                         //   sets grid-map sampling (§4.5) + bounds the envelope
                                         //   (a finite box restricts long-diffusion-time replays);
                                         //   default "periodic" for back-compatibility
    "substrate_frame": [[1,0,0],[0,1,0],[0,0,1]],  // OPTIONAL 3x3 orthonormal, right-handed frame
                                         //   R (§4.2): columns = canonical substrate basis in stored
                                         //   coords; column 3 (z) = primary axis (fibre), columns 1-2
                                         //   = a FIXED perpendicular basis (no free rotation about the
                                         //   fibre). A canonical direction maps to stored as R·d.
                                         //   Omit => identity (positions already canonical).
    "cell_size": 5e-6,                   // m, field-map grid spacing; present iff Field grid maps
    "delta_chi_a": -0.1e-6,              // anisotropic susceptibility scale the Field basis/maps
                                         //   are normalized to; present iff Field tier
    "field_grid": {                      // present iff `susc_field_basis` (§6.4.1); how it was built
      "origin": [0,0,0], "voxel_size": [2e-7,2e-7,2e-7],   // m
      "components": ["iso_local","iso_P","aniso_G"],        // order of the 13 basis channels
      "refocus": "gate" }                //   phase gate s(t_k) supplied at replay (§6.6)
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
    "gradient": true,                    // C0 — always true
    "bulk_relaxation": true,             // C1 — intrinsic per-compartment T1/T2
    "surface_relaxivity": false,         // C2 — any surface relaxivity
    "field": true,                       // C3 — any B0 / susceptibility off-resonance; orientation too when the ℓ=2 maps are present
    "magnetization_transfer": false,     // C4 — magnetization transfer (vector-Bloch)
    "diffusivity_fixed": true,           // D is a fixed pack property, not a knob
    "acquisition": { "b_max": 3.0e9, "ogse_periods": [], "B0_list": [3.0, 7.0],
                     "field_modes": 32,              // K temporal modes kept for susc_field_basis
                     "max_refocusing_pulses": 16 }   // ~K/2 — the gate bandwidth K buys (§9.4 r5)
  },
  "fidelity": {                          // REQUIRED iff any lossy codec (§9); codec-agnostic
    "err_max": 0.0015,                   // worst decoded-vs-raw replay error over `battery`
    "floor_max": 0.0074,                 // irreducible MC noise floor (same ensemble)
    "within_2x_floor": true,             // err_max <= 2*floor_max -> loss is negligible
    "battery": "..." },
  "mt": {                                // REQUIRED iff C4 parametric two-pool (§6.5.1); else absent
    "model": "two_pool",                 //   fixed-at-walk-time pool descriptor (not knobs)
    "f_bound": 0.081, "k_forward": 88.3, //   dimensionless, s^-1  (from S/V, kappa, dwell)
    "T2_bound": 1.0e-5, "T1_bound": 0.44,//   s
    "S_over_V": 4.4e6 },                 //   m^-1, the geometry the pool was derived from
  "provenance": {                        // RECOMMENDED — how the walk was made
    "generator": "dmipy-sim", "generator_version": "...",
    "substrate": {                       // WHAT was walked (descriptive; geometry format is out of scope, §1)
      "representation": "analytical",    //   "analytical" | "mesh" | "voxel" | ...
      "description": "PackedMyelinatedCylinders, gamma radii, f=0.55",
      "ref": null },                     //   optional URI / DOI / hash of the source substrate
    "permeability": 0.0 },               // the FIXED permeability the walk used (§7), if any
  "license": "CC-BY-4.0",                // REQUIRED — the SOURCE substrate's license
  "citation": "Fick RHJ, dmrai-lab (2026) ..."   // REQUIRED
}
```

Rules: `diffusivity` and `seed` are fixed pack properties, not knobs. `license` records the **source substrate's** license and MUST NOT relicense upstream geometry. Any tier flag set `true` in `replay_envelope` MUST have its channels present (§7) and its `per_comp`/`walk_params` fields populated. A tier with two representations is satisfied by **either** one: C3 (`field`) by the `susc_field_{0,C,S}` grid maps **or** the `susc_field_basis` channel (+ `walk_params.field_grid`); C4 (`magnetization_transfer`) by the C1 `bound` occupancy column **or** the `mt` pool descriptor. The envelope flags use explicit, self-describing names; *compatibility:* across the `0.x` line a reader SHOULD also accept the pre-rename aliases `T1T2`/`relaxation→bulk_relaxation`, `rho→surface_relaxivity`, `B0_any`/`orientation_any`/`field_offresonance`/`field_orientation→field`, `mt→magnetization_transfer` (and ignore the retired `rf`/`permeability` flags). The substrate bank derives its catalog card from `provenance` and the optional Croissant sidecar (§12); the card's layout is the bank's concern, outside this format.

---

## 11. Declared envelope and self-certification

A pack is a **certificate**, not a black box. Two objects make its guarantees explicit:

- **`replay_envelope`** — the *domain of validity*: which tiers, and the acquisition/field range (`b_max`, OGSE periods, `B0_list`, `max_refocusing_pulses`, …) over which the producer stands behind the pack. Its limits have two sources — what the *walk* covers (diffusion time, box size) and what a capability-narrowing codec leaves servable (§9.4 rule 5) — and a replayer treats them identically. A request inside the envelope is guaranteed; a request outside is a *known limit*, and the replayer MUST surface it as such (§13), never silently extrapolate.
- **`fidelity`** — the *measured* error of any lossy encoding against the noise floor (§9.4).

The substrate bank (a separate service) accepts a pack by validating it against this specification (§13) and MAY re-measure `fidelity` independently.

---

## 12. Container format

- **Arrays:** a single **safetensors** file. Safetensors is REQUIRED because it (a) carries **no executable code** — a pack is safe to download and memory-map from untrusted sources — (b) is strongly typed, and (c) is zero-copy. Under the identity codec each channel is one tensor under its channel name (§5); under another codec a channel's tensor keys are those the registry defines for that method (`CODEC_REGISTRY.md`).
- **Metadata:** the JSON object of §10, serialized as a single string and stored in the safetensors `__metadata__` header under the key **`"rpk"`**. Producers MUST write it there; a sibling `.json` copy is OPTIONAL for human inspection. *Compatibility:* across the `0.x` line a reader MUST also accept the legacy header key `"json"` (used by the reference implementation before this spec fixed the canonical key); producers SHOULD migrate to `"rpk"`.
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
5. if any codec is lossy, a `fidelity` object is present (§9.4);
6. `license` and `citation` are present.

A **conformant replayer**:

- implements §6.1 and at least refuses cleanly for tiers it does not implement;
- decodes every codec it claims to support to the §5 array contract before applying §6;
- for any requested knob outside the pack's `replay_envelope`, returns a **capability/domain error**, never a silent approximation (§11).

A **conformant producer** satisfies §8 and §9.4.

Conformance is layered: a producer MAY implement only C0; a replayer MAY implement only C0; they still interoperate on C0.

---

## 14. Versioning and extensibility

- `rpk_schema_version` is semantic. **Minor** bumps add OPTIONAL channels/keys; a replayer MUST ignore channels and metadata keys it does not recognize (forward-compatible). **Major** bumps may change required semantics.
- New channels intended for later standardization SHOULD be introduced under the `x_` prefix first.
- New tiers extend the table in §7 without altering existing tiers.
- **Recognized extensions** (permitted by the invariant, not yet standardized):
  - **Transmit-field (`B1+`) inhomogeneity.** Within a single substrate `B1` is uniform — it varies on the RF-wavelength (cm) scale, far larger than a μm–mm pack — so it is **not** a stored map but a **scalar flip-angle knob** `α → κ_B1·α` on the vector-Bloch path (§6.5), the transmit-side analog of the `B0` *field-strength* scalar (not of the susceptibility map). Uniform pulses are already replay knobs (§6.6); the knob adds only the inhomogeneity scale. When packs are tiled into a **mesoscopic voxel grid** (a brain-sized volume), `B1` inhomogeneity becomes a **per-voxel scalar map** — one `κ_B1` per voxel position, applied uniformly inside each voxel's substrate — which is a property of the *assembly*, not of a single pack.
  - a fully general per-voxel rank-2 **susceptibility-tensor** field (§5.2, §6.4). The `0.2.0` `susc_field_basis` (§6.4.1) already covers an isotropic source plus a single rank-2 anisotropic director (the myelin case) for arbitrary morphology; the open extension is a spatially-varying full tensor;
  - membrane **permeability/crossing**, should it ever become a replayable quantity (§7).

---

## 15. Security considerations

A pack carries data only. Safetensors executes no code on load, and the metadata and Croissant sidecar are inert JSON. Replayers MUST NOT `eval`/`exec` any pack content. Field-map and position arrays are bounded numeric data; a replayer SHOULD sanity-check shapes against `walk_params` before allocating.

---

## 16. References

- RFC 2119 — Key words for requirement levels.
- safetensors — https://github.com/huggingface/safetensors
- Croissant (MLCommons) — http://mlcommons.org/croissant/
- This document is the citable, versioned reference for the format; its Zenodo DOI ([10.5281/zenodo.21692505](https://doi.org/10.5281/zenodo.21692505)) establishes priority on the format and infrastructure.

---

## Appendix A. Minimal conformant pack (Tier C0)

Arrays: `positions` `(N_w, N_t, 3)` float32.
Metadata: `rpk_schema_version`, `id`, `walk_params{n_walkers,n_t,dt_traj,T_max,diffusivity,seed}`, `compression{method:"identity"}`, `replay_envelope{gradient:true, ...all others false...}`, `license`, `citation`.
Replay: §6.1 only. This is the floor of interoperability — a walk and nothing else — and every richer pack is this plus declared channels. See `examples/`.
