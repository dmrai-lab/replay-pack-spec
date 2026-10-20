# The Replay Pack Specification (`.rpk`)

**Version:** 0.4.4 (draft for comment)
**Status:** Draft, pre-publication — **everything stays below `1.0` until the format is published**, and the container line is `0.x` for the same reason. Under semantic versioning a `0.x` minor bump MAY break, which is the honest description of where this is: the previous `1.x` numbering promised that "field names and metadata keys are frozen", a promise `0.3` breaks by removing a channel. Stable enough to implement against, but pin a version and expect to re-encode. `0.4.0` is the second **breaking** revision and states the rule the format was built on: **a pack stores channels and no tissue value**. The `per_comp` block (per-pool `T2` / `T1`) and the bound pool's `T2_bound` / `T1_bound` are removed from the metadata -- they were replay knobs stored as pack properties -- and a pack RECOMMENDS embedding the **substrate specification** it was walked from (`substrate`, SUBSTRATE.md) so that its compartment ids have a definition and a replayer can read the nominal tissue values from the specification, explicitly, rather than from the pack. A pack carrying `per_comp` MUST be refused. `0.3.0` is the first **breaking** revision: the `bound_fraction` channel is removed and magnetization transfer's bound pool becomes a `bound` **occupancy column** of `compartment` (§5.2, §5.3, §6.5) — binding is an occupancy, so it needed no channel of its own, and C4 remains a tier because a tier is defined by replay capability, not storage. The equilibrium-start requirement moves with it, anchored to the column rather than to the tier flag (§8.8), since bound occupancy biases the scalar replay on its own. A pack carrying the retired channel MUST be refused, not reinterpreted. `0.2.3` separates **what a channel contains** from **how it is stored** where the two had blurred: it pins the previously-unstated component order and degrees of freedom of `susc_field_basis` (§5.3), adds interface rule 5 — a codec that narrows *reach* rather than accuracy declares it in `replay_envelope`, not `fidelity` (§9.4, §6.4.1) — and with it two OPTIONAL envelope keys, `acquisition.max_refocusing_pulses` and `acquisition.field_modes`. No new required fields. `0.4.4` states **one basis**: a pack is bands on the continuous sine basis of §9.5.3 and nothing else, with `walk_params.block`, `walk_params.qv` and `fidelity.terms` REQUIRED (§9.5.5); occupancy is exact in time, as a static label or an event list (§9.5.7); and **storage windows** -- a walk stored as consecutive windows of one duration, each encoded and certified on its own and sharing its boundary instant with the next -- with the window table `walk_params.windows = {n, T, walks}`, the per-window certificates `fidelity.windows` and the `s{i}/` tensor prefix (§4.3, §10). A pack on a save grid is refused by name and re-walked, never converted (§9.5.1). `0.4.2` is clarification-only (no new required fields): it states the two **similarity views** under which one walk is read as another situation -- in time, `(a D, a kappa)` on a grid divided by `a`, `a >= 1`; in space, every length times `s` at `s^2 D` and `s kappa`, served up to `b_max / s^2` -- with what each keeps, moves and refuses (§6.8), corrects §7 and §10 where they called the diffusivity and the permeability fixed without qualification, and adds the OPTIONAL provenance records `diffusivity_scaled` and `length_scaled` a view carries (§10). `0.4.1` is clarification-only (no new required fields): it states that **an acquisition ends at its echo** and that **a band-limited prefix is a re-encoding** with its own certificate and a `provenance.prefix` naming its parent (§3), and adds interface rule 6 -- **a temporal band is a frequency**, `K / (2 · T_max)`, RECOMMENDED as `compression.temporal_bandwidth_hz` (§9.4, §10). `0.2.2` is clarification-only (no new required fields): it states **time-additivity / resume** — a walk MAY be extended in time by resuming from its final state, the forward complement of the TE-prefix property (§3) — and RECOMMENDS **archiving the geometry of non-reproducible generated substrates** so that resume and walker-shards remain possible (§8.9). `0.2.1` added the OPTIONAL orientation frame `walk_params.substrate_frame` (§4.2, §10). `0.2.0` added two OPTIONAL *representations* (a per-walker Field basis, §6.4.1; a parametric two-pool MT model, §6.5.1) and stated the additive-shard property (§3). Packs written against `0.1.x`/`0.2.x` are **not** conformant to `0.3` and must be re-encoded — that is what a `0.x` minor bump is allowed to mean, and the retired keys are refused by name so a stale pack fails loudly rather than decoding to plausible wrong numbers.
**Container schema described:** `rpk_schema_version = "0.4"`
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
6. **representation and storage** — how `positions` and the other channels are stored (bands on one basis, and occupancy as events), orthogonal to a channel's meaning (§9);
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
- **Trajectory.** The position of a walker as a function of time on `[0, T_max]`: `r_i(t)`.
- **Running integral.** The function `I(t)` of a walker that a channel stores: the position `r(t)`, the cumulative boundary local time, the integral of the sampled field (§9.5.2).
- **Band.** One coefficient `β_k` of a running integral on the sine basis of §9.5.3. A pack stores `K` of them per axis or component, beside the two exact endpoints.
- **Block.** The length `h` over which the producer summed a channel's increments while the walk ran, with the moment order `m` it kept within a block (§9.5.5).
- **Channel.** A named array in the pack (e.g. `positions`, `compartment`). Some are REQUIRED, most OPTIONAL.
- **Replay knob.** A physical parameter applied *at replay time* — gradient waveform `G(t)`, static field `B0`, substrate orientation, `T2`/`T1`, surface relaxivity `ρ`, susceptibility, and the bound-pool relaxation/off-resonance of magnetization transfer. Knobs are **not** stored per value; they parameterize the replay operation.
- **Replay operation.** A deterministic function mapping (channels, knobs) → signal (§6).
- **Capability tier.** A named group of channels that unlocks one class of replay (§7).
- **Declared envelope.** The region of knob space over which the producer certifies the pack's replay to a stated fidelity (§11).
- **Representation.** The *form* every channel of a running integral is stored in: bands on the one sine basis of §9.5. It sets replay cost and the exactness *domain* (§9.2), never a pack's meaning.
- **Codec.** A per-channel **storage** encoding of a channel (§9.4): the band tensors, the static label or event list of the occupancy, the field modes.
- **MC noise floor.** The irreducible variance of a finite-`N_w` ensemble estimate; the reference against which a pack's error is judged (§9.4).
- **Shard.** One pack among several that walk the *same* substrate/physics with *different* seeds; their walker ensembles combine into a larger effective ensemble (§3, additive property).

---

## 3. The replay model (the invariant)

Replay rests on one physical fact.

> **Replay invariant.** For a fixed substrate geometry, a fixed (possibly per-compartment) free diffusivity `D`, and a fixed pseudo-random seed, the ensemble of walker trajectories `{r_i(t)}` and the boundary-contact record are **independent of** the gradient waveform `G(t)`, the static field `B0`, the substrate orientation, the relaxation times `T2`/`T1`, the surface relaxivity `ρ`, the susceptibility distribution, and the bound-pool relaxation/off-resonance of magnetization transfer.

All of those quantities enter the *signal* only through a phase or a (log-)weight accumulated **along** an already-determined trajectory. Therefore they are **replay knobs**: the trajectory is computed once and stored; the signal for any knob setting is a cheap functional of the stored channels (§6).

*Magnetization transfer is split.* The MT **binding kinetics** (binding rate and dwell time) are **not** knobs: a bound spin is frozen for its dwell, so they change the trajectory itself and are fixed at walk time — recorded through the C1 `bound` occupancy column (§5.2), exactly as membrane permeability (§7) is a walk-time property. Only the **bound-pool relaxation/off-resonance** (`T2_b`, `T1_b`, `Δω_b`) are replay knobs (§6.5). A pack MAY instead carry MT in a coarser but cheaper **parametric two-pool** representation (§6.5.1) whose pool size and exchange rate are derived from the substrate's surface-to-volume ratio — no per-walker binding channel, and no binding walk — with the RF saturation as the replay knob; the two representations are the emergent and mean-field forms of the same C4 tier.

Three consequences are used throughout:

- **(TE-prefix property.)** A walk stored to `T_max` contains every shorter echo time as a leading prefix: the walk on `[0, TE']`, `TE' ≤ T_max`, is a valid, converged walk to `TE'`. A replayer MAY read it directly from the bands (§9.5.4). Two consequences for a replayer and a producer. **An acquisition ends at its echo:** every replay operation of §6 -- the gradient phase, the relaxation and surface sums, the field gate -- is evaluated over the acquisition's own duration, and a walk stored beyond a shorter acquisition's echo contributes nothing to that acquisition (a 10 ms spin echo on a 20 ms pack carries `exp(-10 ms / T2)`, never `exp(-20 ms / T2)`). **A band-limited prefix is a re-encoding**, not a slice: a truncated representation (§9) stores `K` temporal modes over `T_max`, so its prefix to `TE'` is the decoded channels cut at `TE'` and re-projected onto `K'` modes, starting from `K' = ⌈K · TE' / T_max⌉` -- the same **bands per second** (§9.4 rule 6) -- and certified as any lossy pack is (rule 4), against the parent's decoded prefix. The certificate decides, not the ratio: the truncation error of a pinned residual has a floor in the band count that does not shrink with the duration, so a producer raises `K'` (up to the parent's `K`) until the certificate passes and records every band tried (reference measurement: a 100 ms walk at `K = 48` prefixes to 25 ms at `K' = 12` but needs `K' = 12`, not 6, at 12.5 ms). The child MUST name its parent in `provenance.prefix` (the parent's `id` and content digest, `T_max`, `K`, and the cut `TE'`, `K'`), and a child whose certificate fails MUST be refused rather than written.
- **(Separation of state and knobs.)** A pack stores only *walk state*. It never stores a signal, a b-value, or a waveform. Signals are always derived at replay time.
- **(Resolution limit.)** The band count `K` sets the temporal bandwidth of every replay: a pack reproduces any `G(t)` whose spectrum lies below `K / (2 · T_max)` hertz (§9.4 rule 6); structure above it was never stored and cannot be replayed. The certificate (§9.5.5) bounds what the bands beyond `K` and the block quadrature cost.
- **(Additivity — walkers.)** Because walkers are independent realizations, two packs that share the *same* substrate geometry, diffusivity, duration `T_max`, band count `K`, block `(h, m)` and channel set but use **different seeds** are **shards** of one larger ensemble: a replayer MAY pool them by concatenating their per-walker channels (or averaging their signals weighted by `N_w`), and the pooled MC noise floor falls as `1/√(Σ N_w)`. A bank therefore need not commit to a final `N_w` up front — it MAY publish a first shard and append more later to drive the floor down, with no re-walk and no change to any published shard. Shards MUST agree on the fixed-at-walk-time properties below; only the seed (and hence the walkers) differs.
- **(Additivity — time; resume.)** Because the walk is **Markov** — a walker's position is its full state — a walk to `T_max` MAY be **extended in time** to any `T_max' > T_max` by resuming each walker from its `T_max` position, continuing the process for the additional `T_max' − T_max`, and storing the continuation as further storage windows (§4.3). The `[0, T_max]` segment is preserved **exactly**; the extension is a valid (fresh-seed) continuation of the same process. This is the forward complement of the TE-prefix property: a bank MAY publish a **short** walk and lengthen it later — with no re-walk of the existing segment — to reach longer echo/diffusion times. Two conditions: (i) the extension MUST use the identical fixed-at-walk-time properties below (geometry, diffusivity, block `(h, m)`, per-compartment `D`, permeability, MT kinetics); and (ii) any per-walker walk-time state that persists across the seam — notably the MT bound/free occupancy (§6.5) — MUST be carried from the `T_max` state (a parametric-MT or no-MT pack has none, so it resumes freely; an emergent-MT pack MUST supply it or MUST NOT be resumed). A pack is therefore extendable along **both** axes — more walkers (shards) and more time (resume) — provided the substrate geometry is available to re-walk (§8.9).

The invariant defines the boundary of the format. The following are **fixed at walk time — not replay knobs**; changing any requires a new walk and a new pack:

- the **substrate geometry**;
- the **diffusivity** `D`;†
- the **random seed**;
- **membrane permeability / crossing** — a crossing is a per-event random draw that alters the trajectory (§7);
- the **MT binding kinetics** (rate, dwell) — recorded in the C1 `bound` column (§6.5);
- the **block** `(h, m)` — how the walk was accumulated, which sets the quadrature term of the certificate (§9.5.5).

† `D` sets the trajectory, so it is recorded as a fixed pack property (§10). The sole exception is the trivial Brownian rescaling — a walk at `D` over `t` equals one at `D/α` over `αt` — which holds **only** in free diffusion, or when the geometry's lengths and the block length `h` are co-scaled consistently (§6.8). In restricted geometry the fixed wall scales break it, so `D` is **not** a general replay knob.

---

## 4. Conventions

### 4.1 Units — SI everywhere, no implicit scale factors

Two kinds of quantity carry units in this standard, both fixed to SI so that a pack from one generator reproduces **identically** under another implementation's replayer. A producer MUST convert from any internal units on write; a replayer MUST assume SI.

**(a) Quantities stored in the pack** — the format proper.

| Stored quantity | Symbol | Unit | Channel / field |
|---|---|---|---|
| position | `r` | **metre (m)** | `positions` (continuous lab frame, §4.2) |
| duration of the stored walk, block length | `T_max`, `h` | second (s) | `T_max`, `walk_params.block.h` |
| diffusivity | `D` | m²·s⁻¹ | `walk_params.diffusivity` (fixed pack property, not a knob) |
| boundary local time | `ℓ` | metre (m) | `boundary_local_time`, at `ρ/D = 1` (§6.3) |
| relaxation times | `T1`, `T2` | second (s) | *not stored* -- replay knobs applied over the `compartment` channel (§6.2) |
| bound occupancy, spin weight | `b_i`, `w_i` | dimensionless | `compartment` (`bound` column), `spin_weights` |

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
- **Continuous.** Stored positions MUST be the *unwrapped* lab-frame trajectory (§8.2). Periodicity, if any, is a property of the walk's box (§4.5), not of `positions`.
- **Substrate frame (intrinsic orientation).** Because orientation is applied by rotating the acquisition, a replayer needs to know how the substrate is oriented *within* the stored lab frame — otherwise "parallel" vs "perpendicular" B₀ (very different susceptibility dephasing) and each diffusion-gradient direction are undefined relative to the tissue. A pack SHOULD therefore declare `walk_params.substrate_frame` (§10): a **3×3 orthonormal, right-handed matrix `R`** whose columns are a canonical substrate basis expressed in stored coordinates, with **column 3 (`z`) the primary structural axis** (e.g. the fibre-bundle direction) and columns 1–2 a **fixed perpendicular basis**. A single axis is insufficient: it leaves a free rotation about itself, so the perpendicular plane — and hence any non-axisymmetric structure, multi-direction scheme, or b-tensor — would be gauge-dependent from run to run. The producer MUST pin `x`/`y` deterministically, aligning them to genuine in-plane structure where one exists. For a **single, axisymmetric** population, `x`/`y` are a fixed but arbitrary basis (e.g. Gram–Schmidt from a fixed reference). For a **multi-population (e.g. crossing)** substrate the producer MUST anchor the frame to a **primary** population, never to a principal-axis/PCA average over all positions — that average is the bundles' *bisector*, which is no real population's axis and rotates as the crossing angle opens, leaving a crossing-angle series with no shared reference. The RECOMMENDED convention: `z` = the primary population's mean axis (chosen by a deterministic rule declared in provenance — e.g. largest by volume, tie-broken by index, independent of packing RNG); `y` = the secondary population's mean axis projected perpendicular to `z` (so the crossing opens into `+y`); `x = y × z`. The primary population then sits at `z` in every pack of a crossing-angle series and the secondary rotates `z → +y` (one bundle held straight, the other rotating), so replays are directly comparable. Two populations fully pin the frame; further populations are merely described by their axes within it. A direction given in the **canonical frame** (`z` along the primary axis) maps to stored coordinates as `d_stored = R · d_canonical`; a replayer orients an acquisition onto the pack by this rotation. Every pack SHOULD carry `substrate_frame` for uniform behaviour, **including isotropic substrates**, for which `R` MAY be the identity but MUST still be fixed. A pack that omits it MUST be treated as `R = I` (positions already in the canonical frame); producers of oriented substrates are strongly encouraged not to rely on this default.

### 4.3 Time
The stored walk is a function of time on `[0, T_max]`, `T_max = walk_params.T_max` (s). There is no save grid: a channel is the function of §9.5.4, and the temporal resolution of a pack is its band `K / (2 · T_max)` (§9.4 rule 6). The integrator's step is a producer concern and is neither stored nor standardized; what is stored of how the walk was accumulated is the block `(h, m)` (§9.5.5).

**Storage windows.** A producer MAY store a walk as `n` consecutive **storage windows** of one duration `T`, consecutive windows sharing their boundary instant, so that `T_max = n T`. A storage window is a unit of storage, not of acquisition: each is encoded on its own (its own endpoints, its own bands) and certified on its own, and the walk is extended (§3, additivity in time), cut to a prefix of whole windows, or joined at a window boundary without re-encoding the others. A pack stored this way declares the window table `walk_params.windows = {n, T, walks}`, `walks` listing the walk runs that produced the windows as `{first, last, seed}` in window indices (a resumed walk adds a run with its own seed); window 0's tensors carry the channel names and window `i ≥ 1`'s the same names under the prefix `s{i}/`, beside the tensors every window shares (the weights, a static label), which appear once; `fidelity.windows` lists each window's certificate (§9.4) and the whole's `fidelity` is the bound over them. A replay of an acquisition spanning several windows is the sum of each window's contraction over the part of the acquisition it holds, the shared boundary instant carrying no weight of its own. A pack without the table is one window, `n = 1`, `T = T_max`; the reference implementation writes the table on every pack.

### 4.4 Signal and sign convention
The canonical signal is complex, `S = ⟨w · e^{+iφ}⟩`, with phase `φ = γ ∫ G·r dt` (§6.1), a **positive** gyromagnetic sign (`γ > 0`), and the `e^{+iφ}` rotation sense. Magnitude DWI is insensitive to this choice, but it is fixed here so the **complex** signal and the susceptibility phase (§6.4) are reproducible across implementations. A real spin-echo value is `|S|` or `⟨w·cos φ⟩`. Averages `⟨·⟩` are weighted means over walkers (weights `w_i`, default `1`).

### 4.5 Boundary
`walk_params.boundary` (§10) declares the walk's box: `periodic`, `reflecting` or `finite`; the default is `periodic`. It bounds the envelope (a finite box restricts long-diffusion-time replays). Stored positions are unwrapped either way (§8.2).

---

## 5. Data model — channels

A pack is a set of named arrays (**channels**) plus a metadata object (§10). A channel is a function of time per walker, stored as the tensors the registry defines (§9); its **decoded** form, read at an instant `t` (§9.5.4), MUST match the type and units below. Element types are given as NumPy dtypes.

### 5.1 Required channel

| Channel | Shape | Dtype (decoded) | Units | Definition |
|---|---|---|---|---|
| `positions` | `(N_w, 3)` functions of `t` | float32 (float64 permitted) | metres | Continuous lab-frame position `r_i(t)` of walker `i`, `t ∈ [0, T_max]`; stored as bands (§9.5). |

`positions` is the **irreducible** channel. A pack with only `positions` (plus metadata) is valid and supports the Gradient tier (§7).

### 5.2 Optional channels

| Channel | Shape | Dtype | Units | Definition | Unlocks |
|---|---|---|---|---|---|
| `spin_weights` | `(N_w,)` | float32 | — | Per-walker statistical weight `w_i` (e.g. compartment volume weighting). Default `1`. | all tiers |
| `compartment` | `(N_w,)` piecewise-constant functions of `t`, per column | int8 / int16 | label | **Occupancy over the declared pools**, as one or more named columns. The `comp` column is the exclusive geometric axis: the compartment id of walker `i` at time `t`, id **0 the extra-cellular / free compartment by convention**, positive ids the pool ids of the embedded `substrate` specification. Further columns are occupancies on **independent** axes — `bound` for the MT macromolecular pool (§6.5), `0` free and `1` bound. A column is a static label or an event list (§9.5.7). | Relaxation (`comp`); Magnetization transfer (`bound`) |
| `boundary_local_time` | `(N_w,)` functions of `t` | float32 | m (see below) | Cumulative wall-contact measure `B_i(t)` of walker `i`, nondecreasing, with **`ρ/D = 1`** (dimension of length ÷ diffusivity is folded so that replay multiplies by the desired `ρ/D`); stored as bands (§9.5.2). See §6.3. | Surface |
| `susc_field_basis` | `(N_w, 13)` functions of `t` | float32 | per unit susceptibility (dimensionless) | Per-walker off-resonance field basis `b(t)` **accumulated along the walk at the walker's own position** — the *general* Field representation (§6.4.1). The 13 components are `iso_local` (1), `iso_P` (6), and `aniso_G` (6): the isotropic-source and rank-2 anisotropic-source responses, projected onto the 6-vector `Q(H)` of the main-field direction at replay. The component order is fixed by §5.3; the 13 components carry 12 degrees of freedom (§5.3). Stored as `K` temporal modes of each component (§9.5.2). Per-walker ⇒ needs a walker-preserving codec (§9.4). | Field |

The Field tier has **one representation**, `susc_field_basis` — the field basis **accumulated along each trajectory** (§6.4.1), the analog of `boundary_local_time` for the field tier. It carries no invariance assumption and so serves **arbitrary 3-D morphology** (undulating/dispersed meshes) and **any** main-field direction via the full rank-2 `Q(H)` contraction; the cost is per-walker storage. Accompanied by `walk_params.delta_chi_a`, the field-grid metadata used to build it (`walk_params.field_grid`, §10) and the read interval `walk_params.field_read_interval` (§9.5.5). No replay reads a field grid (§6.4).

A fully general per-voxel susceptibility *tensor* field remains a §14 extension; the `susc_field_basis` set covers an isotropic source plus a single rank-2 anisotropic director (the myelin case).

A producer **MAY** define **additional** channels prefixed `x_` (e.g. `x_temperature`). Replayers **MUST** ignore unknown channels they do not implement (§14).

### 5.3 Channel invariants

- Every per-walker channel MUST share `N_w` with `positions` and the duration `T_max`.
- `compartment` MUST carry the exclusive `comp` column; its id `0` MUST denote the extra-cellular/free pool and the other ids are the pool ids of the embedded `substrate` specification (§10) when the pack carries one, and MUST otherwise be described in `provenance`. A `comp` column whose every walker keeps one label over the walk MUST be stored as one label per walker (the registry's `static` kind); any other column is an event list (§9.5.7). A reader takes the label at an instant `t` as that of the last event at or before `t`.
- At every instant the occupancy over all declared pools sums to 1. Columns are the **compact encoding of that joint simplex**, not independent quantities: with a bound occupancy `b` and geometric label `g`, the pool occupancies are `f_bound = b` and `f_g = (1 − b)·[comp = g]`. Replay weights the per-pool rates by them, `R(t) = Σ_c f_c(t) R_c` — which is what makes a geometric label (one-hot) and MT binding (an occupancy on an independent axis) the same object.
- The `bound` column MUST be `0` or `1`, and at `t = 0` it MUST be the equilibrium occupancy (the walk is pre-burned-in; §8.8).
- `boundary_local_time` MUST be non-negative and nondecreasing in `t`, and expressed in the `ρ/D = 1` normalization of §6.3.
- `susc_field_basis` MUST order its component axis (the **last** axis) as `iso_local`, then `iso_P` as `(xx, yy, zz, xy, xz, yz)`, then `aniso_G` in the same six-component order — so that slices `[..., 1:7]` and `[..., 7:13]` contract directly against `Q(H)` as written in §6.4.1. A substrate with no anisotropic source stores an all-zero `aniso_G` block.
- `susc_field_basis` carries **13 components but 12 degrees of freedom**: the isotropic block obeys the exact trace identity `iso_P_xx + iso_P_yy + iso_P_zz = 3 · iso_local`, so `iso_P_zz` is recoverable from the other three. This holds to machine precision when the isotropic source grid the basis was built from is un-apodised, and is **broken** by a k-space window applied to `iso_P` but not to `iso_local`. A producer MUST verify the identity rather than assume it; a codec MAY exploit it (dropping and reconstructing `iso_P_zz`) only where it verifies.

---

## 6. Replay operations (normative reference semantics)

These define the meaning of the channels. A **conformant replayer** MUST reproduce these operations (to the decoded arrays' numerical precision). `γ` is the proton gyromagnetic ratio. Every contraction is an integral over time on `[0, T_e]`, `T_e ≤ T_max` the duration of the acquisition (§3), taken over the pieces of the waveform or of the occupancy in closed form; no grid of the pack enters.

The operations split by their dependence on position: the gradient phase (§6.1) and the relaxation/surface log-weights (§6.2–6.3) are **linear** in the stored positions, whereas the Bloch/MT evolution (§6.5) is **nonlinear** and the off-resonance phase (§6.4) is a per-walker contraction of the field the walk accumulated. This is informative only — the signal is identical either way — but it lets a replayer evaluate the position-linear operations directly in the bands' coefficient space without evaluating the trajectory (§9.1; exactness domain §9.2).

### 6.1 Gradient phase (Gradient tier — always available)

The stored walk is the function `r_i(t)` on `[0, T_max]` (§9.5.4). For a gradient waveform `G(t)` (T·m⁻¹),
sample-and-hold on **its own** grid, the accumulated phase of walker `i` is the exact integral

```
φ_i = γ ∫₀^{T_e} G(t) · r_i(t) dt ,
```

which on the bands is the closed form of §9.5.3, each piece of the waveform integrated against the
basis. A replayer MUST evaluate the integral over the waveform's own pieces and MUST NOT
resample `G` onto a grid of its own. A waveform MUST NOT extend beyond `T_max` (the path ends there).

The signal is the weighted ensemble mean

```
S(G) = Σ_i w_i · e^{i φ_i} / Σ_i w_i .
```

This yields any b-value, any b-tensor, OGSE/PGSE/arbitrary `G(t)` whose spectrum the retained bands carry (the band limit of §9.4 rule 6), and (by sweeping `q`) the ensemble average propagator. **No stored quantity depends on `G`** — this is the whole point.

### 6.2 Bulk relaxation `T2`/`T1` (Bulk-relaxation tier)

This tier replays the **intrinsic per-compartment** relaxation of each pool's water — distinct from surface relaxivity at walls, which is the separate Surface tier (§6.3).


With per-compartment transverse/longitudinal rates and the `compartment` channel `c_i(t)`, accumulate a log-weight

```
log W_i = − ∫₀^{T_e} [ χ(t) / T2[c_i(t)]  +  (1 − χ(t)) / T1[c_i(t)] ] dt
```

where `χ(t) ∈ {0,1}` encodes transverse vs longitudinal periods of the sequence (a spin-echo is transverse throughout its encoding). The integral is exact over the runs of the occupancy (§9.5.7) and the pieces of `χ`. The relaxation-weighted signal replaces `w_i → w_i · W_i` in §6.1. `T2[·]`, `T1[·]` are replay knobs supplied by the caller, whose default MAY be the nominal values of the embedded substrate specification (§8.5).

### 6.3 Surface relaxivity (Surface tier)

With the `boundary_local_time` channel `B_i(t)` (stored at `ρ/D = 1`) and a chosen surface relaxivity `ρ` and diffusivity `D`,

```
log W_i += − (ρ / D) · [ B_i(T_e) − B_i(0) ] .
```

Because the channel is stored normalized, one walk serves every `ρ`. This is a replay knob; `ρ = 0` recovers §6.2.

### 6.4 Static off-resonance / susceptibility (Field tier)

The Field tier stores the off-resonance field basis `b(t)` accumulated along each walker's trajectory by the walk, at the walker's own position (§6.4.1, §9.5.5). **No replay reads a field grid**, and the field is never sampled at the stored positions, which are band-limited (§9.5.4). A field grid MAY accompany a pack as provenance or for visualisation, as a companion artefact or an `x_` channel; it is not a field route, and a pack that carries `susc_field_0`, `susc_field_C` or `susc_field_S` as its field is refused (§9.5.1). A fully general rank-2 susceptibility-tensor field is a documented extension requiring a richer basis set (§14).

The off-resonance contributes to the signal only through the **acquisition's phase history**, which is a replay knob (§6.6). The replayer supplies a transverse-phase gate `s(t)` derived from the pulse sequence, and

```
φ^χ_i = ∫₀^{T_e} s(t) · Δω_i(t) dt ,      φ_i → φ_i + φ^χ_i
```

with `Δω_i(t)` of §6.4.1. The susceptibility phase adds to the gradient phase inside the *same* exponential, so the diffusion×susceptibility covariance (the cross-term) is preserved. Field strength, orientation, susceptibility scale, and sequence are all knobs; one walk serves any combination.

#### 6.4.1 Per-walker field basis (general representation)

The producer stores the field basis accumulated along each trajectory in `susc_field_basis` (§5.2). It carries no invariance assumption: an undulating, dispersed, or otherwise irregular 3-D substrate is served as a straight one. Its 13 components per walker and instant are the geometry-only responses `iso_local`, `iso_P[6]`, `aniso_G[6]`. At replay, for a main-field **unit direction** `H = (H_x,H_y,H_z)`, form the symmetric-tensor contraction 6-vector

```
Q(H) = (H_x², H_y², H_z², 2 H_x H_y, 2 H_x H_z, 2 H_y H_z)
```

and the off-resonance seen by walker `i` is

```
Δω_i(t) = γ · B0 · [ χ_iso · ( iso_local − Q·iso_P )  +  χ_aniso · ( Q·aniso_G ) ]_i(t)
```

which feeds the **same** transverse-phase gate and phase accumulation as §6.4 (`φ^χ_i = ∫ s(t)·Δω_i(t) dt`), preserving the cross-term identically. `Q(H)` carries the full main-field direction, so **any** B0 orientation is a replay knob. Because the basis is accumulated by the walk, it contracts without a field grid at replay; the cost is per-walker storage.

**The temporal storage of this channel costs reach, not accuracy.** The phase reaches the field only through the gate, `φ^χ_i = ∫ s(t)·Δω_i(t) dt`. The channel is stored as `K` temporal modes (§9.5.3, registry), so it is band-limited **in the gate**, not in the field: by Parseval `∫ s b dt = Σ_k ŝ_k c_k`, the retained modes reproduce `φ^χ` for any `s(t)` inside their band while reproducing `Δω_i(t)` itself poorly, because the field along a trajectory is broadband — the walker crosses μm-scale field structure every step. `K` is thus a **capability**: it bounds the refocusing-pulse density the pack can serve (reference implementation: `K ≥ 2·n_refocus`), and MUST be declared in `replay_envelope` (§9.4 rule 5, §11), never read as "the field is accurate to `K` modes".

### 6.5 Magnetization transfer (Magnetization-transfer tier)

With the C1 `bound` occupancy column `b_i(t)` (§5.2), replay blends per-step relaxation and off-resonance toward a bound-pool set `(T2_b, T1_b, Δω_b)` by occupancy, within a **vector-Bloch replay** — RF pulses act as rotations of the magnetization vector, which the scalar log-weight model cannot represent (it has no `M_z` reservoir). Saturation transfer is *emergent*. The bound-pool parameters are replay knobs. This same vector-Bloch path is the general route for *arbitrary RF* (§6.6). See the reference implementation for the Bloch–McConnell blend.

**MT stores nothing of its own.** The bound pool is a column of the C1 occupancy channel, because binding *is* an occupancy: a spin bound for a fraction of an interval has that fraction of it at the bound pool's relaxation, which is exactly the occupancy weighting C1 already defines. What makes magnetization transfer a distinct **capability** tier is the replay side — RF as vector-Bloch rotations (the scalar log-weight path cannot express it), the bound-pool knobs, the equilibrium start below, and the certification duty of §9.2 — not the storage. Two tiers sharing a storage form is normal here: C0 and C2 share a byte-identical bridge layout (§9.3) and remain separate tiers because their knobs differ. A useful consequence of storing binding as C1 rather than as a summary parameter: the intervals between the events of the `bound` column **are** the dwell times, so an emergent pack carries the realised exchange statistics, not just a mean rate.

**Equilibrium start.** A magnetization-transfer pack's `t=0` is the **bound-pool equilibrium**: the walk is generated pre-burned-in (§8.8), so replay begins from a fully-relaxed steady-state occupancy and never sees the fill-up transient. *(This applies to the emergent per-walker representation; the parametric representation of §6.5.1 carries no `t=0` occupancy transient.)*

#### 6.5.1 Parametric two-pool MT (mean-field representation)

The emergent representation above requires the `bound` occupancy column, which a producer obtains from a binding walk that resolves near-wall motion at a fine sub-step. On sub-micron restricting features that sub-step makes the binding walk intractable. A producer MAY therefore declare C4 in a **mean-field two-pool** form instead: rather than a per-walker channel, it stores a small **pool descriptor** in metadata `mt` (§10) — the bound-pool fraction `f_b`, the forward exchange rate `k_f`, and the bound-pool relaxation `(T2_b, T1_b)` — derived from the substrate's myelin **surface-to-volume ratio** `S/V` (with the binding reactivity `κ` and dwell `τ`): `k_f = κ·S/V`, `f_b = k_f·τ / (1 + k_f·τ)`. These are **fixed-at-walk-time geometry/tissue properties** (like `D`), not knobs.

At replay the observable free-pool signal follows the two-pool **Bloch–McConnell** equations under an RF saturation of nutation `w1` and offset `Δ` for duration `t_sat`; a replayer offering this representation MUST expose at least the **Z-spectrum** (free-pool `M_z` vs `Δ`) and the derived **MT ratio**. RF offset, power, and duration are the replay knobs, and so are the bound pool's `T2_bound` / `T1_bound`: the descriptor stores what the walk fixed (`f_bound`, `k_forward`, `S_over_V`) and no relaxation time. This form is validated to reproduce the emergent equilibrium `f_b` (the mean matches `κ·S/V·τ`), but it deliberately does **not** model the binding↔diffusion coupling the emergent walk captures — a bound spin's frozen path — so it is the correct choice for MT saturation / qMT observables and **not** for diffusion–MT coupling. The two are the mean-field and emergent forms of the **same C4 tier**; a pack declares exactly one.

### 6.6 The acquisition is a replay knob (no sequence is baked into the pack)

A pack stores **no** assumption about the pulse sequence. The full acquisition — the gradient waveform `G(t)` **and** the RF / `B1` pulses (their flip angle, phase, shape, and **timing**) — is supplied at replay. Moving the refocusing pulse when the diffusion time or `TE` changes, or switching GRE↔spin-echo↔CPMG↔stimulated-echo/PGSTE, is a replay-time choice against the one stored walk, never a re-walk. (Spatially *uniform* `B1` pulses need nothing stored; `B1` **inhomogeneity** is a scalar flip-angle knob per pack — becoming a per-voxel scalar only when packs are tiled into a mesoscopic grid — not a within-substrate field; §14.)

For the **scalar-phase tiers** (§6.1, §6.4) the sequence enters through a **transverse-phase gate** `s(t)` that the replayer derives from the sequence:

- `s ≡ +1` — gradient echo (GRE), no refocusing;
- `s` flips sign at each 180° refocusing pulse — spin echo, CPMG;
- `s = 0` while magnetization is stored along `z` — stimulated echo (STE/PGSTE);

so `s(t) ∈ {−1, 0, +1}` in the common cases. The **gradient phase (§6.1)** already encodes any bipolar/refocused *gradient* structure in the waveform's polarity and needs no gate; the gate applies to the **static off-resonance** term (§6.4), which a refocusing pulse inverts and a storage period freezes. (The transverse/longitudinal indicator `χ_k` of §6.2 is the same idea for relaxation.)

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

### 6.8 Similarity views (normative)

A pack is walked at one diffusivity `D`, one set of wall rates (permeability, surface relaxivity, MT reactivity; m/s) and one size. Two similarities of diffusion let a replayer read the stored walk as a walk of another situation **without decoding and without a new certificate of its own**, because the dimensionless walk -- `sqrt(D Δt) / L`, `kappa L / D`, `rho L / D`, `D T / L^2`, with `L` the smallest feature -- is unchanged. A conformant replayer MAY offer either or both; one that does MUST move every quantity the similarity couples, MUST state the view's own values (its `walk_params.diffusivity`, its embedded `substrate`, its `replay_envelope`), and MUST record the view in `provenance` (§10). A view that moves one quantity alone is a different walk and MUST NOT be presented as the stored one.

- **In time (`a = D' / D`, `a >= 1`).** The bands read on a time axis divided by `a` (`T_max / a`, `h / a`, the band in hertz times `a`) are the walk at `a D`; every permeability and MT reactivity is `a` times the walked value, since the per-encounter probability `kappa sqrt(dt / D)` of the walk's step `dt` is what the path realised; the surface term's rate divides by the new `D` (§6.3). Every length stays. Only `a >= 1` is served: the walk's in-step error goes as `D dt^2` for the step `dt` its step law fixed, which in the view is `D dt^2 / a`, so a stretched time axis (`a < 1`) grows it past what the step law kept below the floor, and a replayer MUST refuse it rather than return it uncertified. The view is shorter, `T_max / a`.
- **In space (`s = L' / L`, both directions).** Every stored length times `s` -- the positions, the boundary local time (§6.3), the field basis's grid metadata `field_grid` (origin and voxel size, §6.4.1), the embedded substrate's lengths -- is the walk at `s^2 D` with every wall rate at `s` times the walked value; the block, the duration, the gates, the occupancy and the field values stay (a dipole field depends on shape, not size). The gradient phase of the view at `G` is the stored walk's at `s G` (§6.1), and the surface term at `rho` is the stored walk's at `rho / s`. The certificate was measured over the battery up to `replay_envelope.acquisition.b_max`, so the view's envelope is `b_max / s^2` (and `tissue.rho2_over_D_max / s`), which it MUST state and MUST hold its replays to (§11).
- **Composed.** A target `(L', D', T')` is reached from the walk `(L, D, T_max)` by the space view then the time view when `D' / L'^2 >= D / L^2` and `D' T' / L'^2 <= D T_max / L^2`; the prefix property (§3) serves a shorter acquisition.

---

## 7. Capability tiers

A producer that models only part of the physics still emits a **fully interoperable** pack; it simply declares fewer tiers. Tiers are **independent** — a pack MAY declare Field without Relaxation.

Tiers are named **`C0`–`C4`** (for *capability*) — deliberately not `T#`, which would collide with the relaxation times `T1`/`T2`. Metadata flags use **explicit, self-describing names** (§10).

| Tier | Required channels (beyond `positions`) | Replay unlocked | Metadata flag(s) |
|---|---|---|---|
| **C0 Gradient** | *(none)* | §6.1 — any `G(t)`, b-tensor, EAP | `gradient: true` (always) |
| **C1 Bulk relaxation** | `compartment` (the `comp` column) | §6.2 — any `T2`/`T1`, given per pool at replay | `bulk_relaxation: true` |
| **C2 Surface** | `boundary_local_time` | §6.3 — any surface relaxivity | `surface_relaxivity: true` |
| **C3 Field** | `susc_field_basis` (stored as `susc_path_dct` modes; + `delta_chi_a`, `field_grid`, `field_read_interval`) | §6.4 / §6.4.1 — any `B0`, susceptibility, orientation | `field: true` |
| **C4 Magnetization transfer** | C1 `bound` occupancy column (emergent) **or** parametric `mt{}` pool (§6.5.1) | §6.5 / §6.5.1 — magnetization transfer | `magnetization_transfer: true` |

The declared tier set lives in `replay_envelope` (§10). A replayer asked for a knob outside the declared tiers MUST refuse with a clear "capability not present" error, **not** silently return an approximate result (§11, §13).

*On membrane permeability.* Permeability changes the **trajectory** (a crossing is a per-event random draw), so it is not a replay knob: a pack is walked at one permeability, and the only other permeabilities the same walk serves are the pairs `(a D, a kappa)` of the time view and `(s^2 D, s kappa)` of the space view (§6.8), which keep the realised crossing probability per encounter. Any other value requires a new walk. The `0.x` envelope therefore carries no permeability flag; the fixed value the walk used SHOULD be recorded under `provenance` and in the embedded substrate specification.

---

## 8. Generator invariants (what a producer MUST guarantee)

A pack is only replayable if the producer respected the model. A conformant producer:

### 8.1 Determinism
MUST record the `seed` and produce byte-reproducible trajectories from `(geometry, D, seed)` on a fixed implementation. The seed is metadata, not a knob.

### 8.2 Continuous positions
MUST store **unwrapped** lab-frame positions. If the walk used a periodic cell, the producer MUST unwrap before storing (undo per-step jumps `> L/2` on periodic axes; leave non-periodic axes untouched), so that the gradient phase in §6.1 is correct. Storing wrapped positions is a conformance failure: each wrap injects a spurious `q·L` phase.

### 8.3 TE-prefix consistency
MUST ensure that any leading prefix of the walk is itself a converged walk to that shorter time (§3). Producers MUST NOT reorder walkers between a channel and `positions`.

### 8.4 No sequence assumption
A producer MUST NOT bake any pulse-sequence choice into the pack. The acquisition (gradients and RF) is applied at replay (§6.6); the producer's obligation is only to store the trajectory (and any tier channels) over the full walk to `T_max`. In particular there is no "refocusing time" a producer must record — refocusing is a property of the replayed sequence, not of the walk. (A producer that chooses a capability-narrowing codec does declare a refocusing-pulse *density* it can serve — `replay_envelope.acquisition.max_refocusing_pulses`, §9.4 rule 5. That is a limit of the chosen **storage**, not a sequence assumption baked into the walk: the walk carries none.)

### 8.5 Compartment convention
MUST use id `0` for the extra-cellular/free pool. Every other id MUST be defined: by the pool of that id in the embedded `substrate` specification (§10, SUBSTRATE.md), or, for a pack without one, in `provenance`. A producer MUST NOT store a tissue value -- `T2`, `T1`, surface relaxivity, susceptibility, bound-pool relaxation -- as a pack property: those are replay knobs (§2), and the place for their nominal values is the substrate specification. A replayer's default MAY be the **nominal replay** -- the specification's declared values, its `nominal_field_T` as `B0` -- so that a published pack reproduces its paper by firing a pulse at it; every value stays a knob the caller overrides, and a pack without a specification replays the gradient alone.

### 8.6 Units and frame
MUST honor §4 exactly. A producer using non-SI internal units MUST convert on write.

### 8.7 Honest envelope
MUST declare a `replay_envelope` it can defend, and MUST attach a `fidelity` self-report whose battery exercises **every declared tier** (§9.4, §11). A producer MUST NOT declare a tier whose channel it did not actually populate.

### 8.8 Equilibrium start for magnetization-transfer packs
A pack carrying the C1 `bound` occupancy column MUST begin from the **bound-pool equilibrium**. A Monte-Carlo binding walk started from an arbitrary state (e.g. all spins free) shows a transient while the bound fraction relaxes to its steady state `f_b = k_f/(k_f + k_r)`; that transient does not represent a fully-relaxed sample and must not appear in a replay. The producer MUST therefore **equilibrate** the binding dynamics and **discard** that preamble, storing the walk from `t=0 =` the equilibrated state (the `bound` column accordingly starts at equilibrium). Replayers assume `t=0` is equilibrium (§6.5). The producer SHOULD verify the occupancy has reached steady state before saving.

The requirement attaches to the **column, not to the tier flag** — the column's presence is itself the assertion, and a consumer needs no C4 declaration to depend on it. The reason is that bound occupancy biases the *scalar* replay on its own: a spin bound for part of a save carries the bound pool's short `T2_b` for that fraction under the ordinary C1 occupancy weighting `R(t) = Σ_c f_c(t) R_c` (§5.3), with no RF and no `M_z` reservoir involved. A replayer using only bulk relaxation therefore reads the column too, so anchoring the guarantee to the magnetization-transfer tier would leave that consumer relying on a promise nobody had made. No other channel has an initial-condition requirement: the transient is a property of the binding kinetics, not of the walk.

### 8.9 Reproducibility of generated substrates
A substrate from a **stochastic generator that is not deterministically reproducible from its configuration** (e.g. a mesher seeded from wall-clock time or process id) cannot be recovered from the config alone. For such substrates the producer SHOULD **archive the exact geometry used for the walk** — as a hashed sidecar — and reference it in `provenance.substrate.ref` (§10). This is needed not only so the substrate is reproducible and inspectable, but because the additive operations of §3 — **resuming/extending the walk in time** and **adding walker shards** — both *re-walk the same geometry*; without an archived geometry they become impossible once the generator's working directory is gone. The geometry *format* is out of scope for the container (§1), but for interoperability the archive SHOULD use an **open, widely-readable geometry format** (e.g. PLY/OBJ/STL for meshes; NIfTI for label volumes) rather than a producer-specific binary, so any consumer — in any language — can read it. Compression is orthogonal: a plain `xz`/`gzip`/`zstd` wrapper keeps it a one-line decompress; a container SHOULD state which. A substrate that *is* bit-reproducible from its config (config + version + seed) MAY instead reference that recipe.

*Reference implementation.* dmipy-sim exposes this as `equilibrate_binding = 'auto' | 'burnin' | 'fast' | 'off'`: `'burnin'` (the `'auto'` default when MT is on) runs an **adaptive RF-off burn-in in ~dwell-sized chunks until the ensemble occupancy plateaus** (geometry-agnostic, with a convergence flag); `'fast'` seeds the equilibrium occupancy analytically from a known `S/V` when that is position-invariant; `'off'` keeps the all-free start. A producer MAY use any method that satisfies the equilibrium-start requirement above.

---

## 9. Representation and storage (orthogonal to meaning)

**A channel's *meaning* is its decoded function (§5); how it is stored is orthogonal to that meaning.** A pack stores every channel of a running integral as bands on the one basis of §9.5, the occupancy as exact events in time, and nothing as samples on a grid. Two storage choices remain, and neither changes what a pack *means* -- only its size and its replay *cost*:

- the **representation** of `positions` -- the bands (§9.1);
- a per-channel **storage codec** -- the encoding of a channel's tensors (§9.3, §9.4).

Both are declared in metadata `compression`, both MUST decode to the §5 contract, and the concrete methods with their exact stored-key layouts live in the separately versioned [`CODEC_REGISTRY.md`](CODEC_REGISTRY.md), so they evolve without touching this frozen core.

### 9.1 Representation of positions — bands

Positions are stored as the two exact endpoints and `K` bands of the bridge residual on a linear basis, `r = M·c + μ` (a shared basis `M`, per-walker coefficients `c`, §9.5). Because the gradient phase (§6.1) and the separable relaxation/surface log-weights (§6.2–6.3) are **linear** in position, replay is a projection onto `K` coefficients instead of an integral along the path, evaluated without ever reconstructing the trajectory (§6 preamble; mechanics in the registry). The nonlinear operation -- Bloch/MT (§6.5) -- evaluates the function at instants first (§9.5.4).

### 9.2 Exactness and its domain

A pack's exactness has a **domain**:

- the endpoints are exact at every `K` (§9.5.3);
- a pack is exact for any `G(t)` whose induced phase lies in the **span** of its `K` modes;
- there is no full-rank pack, since a diffusing path has bands at every `k`; outside its span the error is bounded and MUST be certified by `fidelity` (§9.4, §9.5.5).

So a pack trades exactness outside its span for size and replay speed, and declares that trade in `fidelity`.

### 9.3 Storage codecs

The channels other than `positions` each use a **storage form** matched to their structure, with the stored keys defined in the registry: the occupancy columns as a static label or an event list (§9.5.7); the contact `B(t)` as bands, as positions; the field as `K` temporal modes. These are pure storage: they decode to the §5 contract and do not change replay cost.

A storage codec MAY offer more than one **container** for the same coefficients -- for the bridge form, the float32 tensor per axis or an integer band container with a per-band scale (registry 0.6). A container is a choice of bytes, not of meaning: the decoded channel and every replay contraction are the same to the container's quantisation, the pack's fidelity certificate (§10) is computed through the container it is written in, and a reader that meets a container it does not implement MUST refuse the pack rather than read its tensors under another layout.

### 9.4 Interface rules

A conformant representation or codec MUST satisfy six rules:

1. **Declared.** The codec is named (with any parameters) in metadata `compression`, e.g. `{"method": "bridge_dst", "K": 256, "walker_preserving": true}`. A replayer that does not recognize the `method` MUST refuse (never guess). Recognized method names and their stored tensor keys are defined in the registry.
2. **Decodes to the contract.** Decoding MUST yield the channel of §5 (type/units, read at an instant). Channel *presence* for tiers (§7) and conformance (§13) is satisfied by the registry-defined stored keys of its codec.
3. **Bands are the baseline.** `bridge_dst` is the codec of `positions` and every replayer MUST support it. A channel has no raw form: a pack is bands (§9.5).
4. **Every pack is self-certified.** A pack is band-limited and accumulated in blocks, so it MUST include a `fidelity` object reporting the maximum replay error over a declared **acquisition battery**, measured against a reference replay of the same acquisition (a per-step accumulation of the same walk, or the forward engine's signal) and against the **split-half Monte-Carlo noise floor** (the ensemble split into two *random* halves -- random because walkers are typically seeded in compartment order). The convention: "lossless to the noise floor" means `err_max ≤ tol · floor_max` (reference `tol = 2`). The battery MUST exercise **every tier the pack declares** -- in particular the *nonlinear* Field (C3) and magnetization-transfer (C4) operations: a truncated linear representation's error on those is generally **larger than, and uncorrelated with,** its error on the gradient phase (§9.2), so a gradient-only battery does not certify them. The certificate's bounded terms are `fidelity.terms` (§9.5.5).
5. **Capability-narrowing ⇒ declared in the envelope, not in `fidelity`.** Most of a pack's storage costs only accuracy, and rule 4 is the whole story. A truncated representation of a channel that replay consumes **non-linearly** can instead cost *reach* -- leaving some requests unservable at any accuracy -- and `fidelity` cannot express that: an `err_max` measured over a battery that never leaves the retained span reports a small number for a request the pack cannot serve at all. Such a codec MUST declare the resulting limit in `replay_envelope` (§11) as well as its `fidelity`. The registered case is the per-walker Field basis stored as `K` temporal modes, which bounds the servable gate bandwidth rather than the field accuracy (§6.4.1) and declares `replay_envelope.acquisition.max_refocusing_pulses`.
6. **A temporal band is a frequency.** A codec that stores a channel as `K` temporal modes over the walk's `T_max` resolves that channel up to `K / (2 · T_max)` hertz, and it is that frequency -- not the count `K`, whose meaning changes with every walk length -- that a scanner envelope is compared against and that a prefix (§3) preserves. A pack SHOULD record it as `compression.temporal_bandwidth_hz = K / (2 · T_max)` (§10) -- a reader derives the same number from `K` and `T_max` when it is absent -- and a bank SHOULD request a pack by bands per second. The frequency is the codec's **reach**: a gradient waveform whose spectrum lies below it is inside the retained span. It is not the accuracy for a given waveform -- that is the product of the waveform's out-of-band spectrum with the channel's residual (the walk's own spectrum falls only as `1/m` in mode index), and remains the business of the `fidelity` certificate (rule 4) and of a consumer's per-waveform check.

   **Inherited certification.** Packs that are blocks of one *fill* -- the same substrate, walk parameters, block `(h, m)`, codec and containers, differing only in the walkers they hold -- MAY inherit the battery measurement instead of repeating it: one or more **certifying** packs of the fill carry the measured `fidelity` (`certified: "measured"`); every other pack carries `certified: "inherited"` and `inherited_from {id, err_max, floor_max}` naming a certifying pack, reports that pack's `err_max` (whole and per family and per tier), and reads its **own** split-half `floor_max` -- and, when partitioned, its own per-voxel floors -- over the same battery in the coded domain (the phases from the stored coefficients, no path decoded). A pack whose codec parameters differ from the pack it cites is non-conformant; a pack that inherited a certificate cannot certify another; a bank MAY re-measure any pack.

The `fidelity` object is **codec-agnostic** — it reports the replay error whatever method produced it — so it lives in the **core** metadata (§10) and is read identically by any replayer or bank; only the codec *algorithms* live in the registry. How to read it: `err_max` is the worst-case replay error over the battery, `floor_max` the irreducible MC noise floor for the same ensemble, and `within_2x_floor: true` means the loss sits below `2×` that floor (scientifically negligible).

One structural rule belongs in the core because it constrains tiers, not any particular algorithm:

> **Walker-preserving requirement.** A codec that **resamples walkers** (a *distributional* codec — it stores a coefficient *distribution*, not per-walker coefficients) breaks per-walker channel alignment. Such a pack **MUST NOT** carry any per-walker channel (`compartment`, `boundary_local_time`, `spin_weights`, `susc_field_basis`) and **MUST** declare only the Gradient tier. Per-walker channels REQUIRE a walker-preserving codec.

A channel's storage never changes what a pack *means* — the decoded channel of §5 and the operations of §6 are fixed. For channels replay consumes **linearly** it is a pure *size/speed/quality choice*, and §9.2's span argument bounds the trade. For a channel consumed **non-linearly** — the Field basis (§6.4.1), Bloch/MT (§6.5) — a truncated representation may additionally narrow the *domain of validity*, which is a matter for `replay_envelope` (rule 5), never a silent property of the stored bytes. See [`CODEC_REGISTRY.md`](CODEC_REGISTRY.md) for the registered representation and storage codecs and their stored-key tables.

### 9.5 The one basis

A pack stores every running integral on one basis, the continuous sines of §9.5.3, accumulated by the walk while it runs (§9.5.5). A coefficient means something only against its basis, and there is exactly one.

#### 9.5.1 One basis, and the refusal

Every band-stored channel -- C0 positions, C2 contact and C3 field (§9.5.2) -- is on the sine basis of `[0, T]` (`T = T_max`, or the window duration of a window table, §4.3). No metadata field names the basis, because there is no other.

**Refusal.** A reader MUST refuse a pack that carries `n_t` or `dt_traj` (in `walk_params`, in its window table, or anywhere under `compression`), `compression.basis`, `walk_params.cell_size` (a field grid offered as a field route), the codec `identity`, or any tensor key of the registry's refused list, and MUST name the key it found and give this instruction: *"This pack carries `<key>`, which this revision does not read; re-walk the pack with a conformant producer. A pack is never converted."* A reader MUST NOT decode any tensor of such a pack: a pack on a save grid (`n_t` saves of spacing `dt_traj`, bands on the saves, or positions as raw samples) stores band tensors under the same keys and in the same widths as a pack of this revision, so a coefficient read against the wrong basis is a plausible number and not an error. The metadata schema rejects such a pack, and `examples/conformance_check.py` names the key and gives the instruction.

#### 9.5.2 Running integrals and the stored vector

A band-stored channel stores the coefficients of one **running integral** `I(t)` per axis or component. `I` is the position `r(t)` for C0, the cumulative boundary local time `B(t)` of §6.3 for C2, and `F(t) = ∫₀ᵗ b dt'` for C3, the running integral of the field `b` the walk accumulated at the walker's own position (§9.5.5). The **bridge residual** is

```
u(t) = I(t) − I(0) − (t/T)·[I(T) − I(0)],        u(0) = u(T) = 0,
```

and the stored vector per axis is `[I(0), I(T) − I(0), β_1, …, β_K]`: the two exact endpoints, then `K` bands of `u`.

#### 9.5.3 The sine basis on `[0, T]`

The basis is

```
s_k(t) = sqrt( 2 / T ) · sin( π k t / T ),     k = 1, 2, …     (orthonormal in L²(0, T)),
β_k = ∫₀ᵀ u(t) · s_k(t) dt ,     I(t) = I(0) + (t/T)·[I(T) − I(0)] + Σ_k β_k · s_k(t).
```

A band has the units of `I` times `sqrt(s)`. `K` is a band count only: there is no full rank, since a diffusing path has bands at every `k`, and what a truncation at `K` costs is the business of the certificate (§9.5.5).

*The bands are the projections of the increments.* Let `e_k(t) = sqrt(2/T)·cos(π k t / T)`. For `k ≥ 1`,

```
β_k = (T / (π k)) · ∫₀ᵀ e_k(t) dI(t) = ( sqrt(2 T) / (π k) ) · ∫₀ᵀ cos(π k t / T) dI(t) .
```

By parts, with the antiderivative `−(T/(π k))·e_k` of `s_k`: the boundary term `[u · (−T e_k/(π k))]₀ᵀ` vanishes because `u(0) = u(T) = 0`; what remains is `(T/(π k)) ∫ e_k du`, and `du = dI − ([I(T) − I(0)]/T) dt` with `∫₀ᵀ e_k dt = 0` for `k ≥ 1`, which drops the ramp. For a path that is not of bounded variation the integral against `dI` is defined by this identity (for a diffusing path it is the Wiener integral). The cosine projection of the increments at `k = 0` is the displacement `I(T) − I(0)`, which the second stored entry holds exactly. This identity is what a walk accumulates (§9.5.5).

*The endpoints are exact.* Every `s_k` vanishes at `0` and `T`, so the series returns `I(0)` and `I(T)` at every `K`, and a stored walk is continued from where it ended (§3, additivity in time) to the stored precision of the endpoints.

*One basis for every channel.* For C3 the stored modes are `c_k = ∫₀ᵀ e_k(t)·b(t) dt` (`e_0 = T^{−1/2}`, `k = 0 … K − 1`), that is `∫ e_k dF`, and by the identity above the sine bands of `F`: C0, C2 and C3 are one basis.

*Gradient phase.* §6.1 on a pack reads

```
φ = γ · [ r(0)·M0 + (r(T) − r(0))·M1 + Σ_{k=1}^{K} β_k · Ĝ_k ],     M0 = ∫G dt,   M1 = ∫ (t/T) G dt,   Ĝ_k = ∫ G(t) s_k(t) dt,
```

each integral over the waveform's own pieces in closed form and over the acquisition's own duration (§3, an acquisition ends at its echo), per gradient axis against the matching axis of `β`. Refocusing (`M0 = 0`) and velocity compensation (`M1 = 0`) remove the first two columns exactly.

#### 9.5.4 Evaluating a pack at an instant

A reader needing the channel at a time `t` -- to draw a path, to continue a walk, to cut a prefix inside a window (§3) -- evaluates the basis at `t`:

```
I(t) = I(0) + (t/T)·[I(T) − I(0)] + Σ_{k=1}^{K} β_k · sqrt(2/T) · sin(π k t / T),        0 ≤ t ≤ T.
```

`t` is a time in seconds. There is no index to round it to and no grid it must lie on; the value at `t = 0` and `t = T` is exact, and in between it is the reconstruction at `K` bands. An array on a grid of the reader's choosing is a sampling of this function made by the reader, not the stored channel. The C2 contact of an acquisition ending at `t_e` is `B(t_e) − B(0)` read the same way, and the field gate (§6.6) contracts against the modes `c_k` of §9.5.3 by Parseval, `∫ s(t)·b(t) dt = Σ_k ŝ_k c_k` with `ŝ_k = ∫ s e_k dt`. The band-limited prefix of §3 is the function cut at `TE'` and re-projected onto the basis of `[0, TE']`; its error is certified as any pack's is.

#### 9.5.5 `walk_params.block`, the quadratic variation, and the certificate

A pack is accumulated while the walk runs, not read off stored samples, so how it was accumulated is part of its error. `walk_params.block = {h, m, precision}` is REQUIRED:

| key | type | meaning |
|---|---|---|
| `h` | number, s | the block length: the producer summed each channel's increments within blocks of this length |
| `m` | integer `≥ 0` | the moment order: within a block it kept `∫ τ^j dI` for `j = 0 … m`, `τ` from the block centre |
| `precision` | `"float32"` or `"float64"` | the floating-point type in which block contributions were accumulated across blocks |

`(h, m)` is a constant of the pack, fixed at walk time like `D`, never a replay knob. `m = 0` is a legal value: the block's displacement alone, at the block's midpoint.

*The quadratic variation.* For a channel with running integral `I`, accumulated over the increments `dI` of the walk's own steps (for the field, of its reads), the **quadratic variation** is

```
QV(I) = Σ (dI)² ,
```

the sum of the squared increments over the walk, averaged over the walkers. The walk records it per channel in `walk_params.qv`, REQUIRED:

| key | channel, `I` | unit | meaning |
|---|---|---|---|
| `positions` | C0, `r(t)` | m² | `[QV_x, QV_y, QV_z]`, one per axis. Free diffusion has `2 D T`; a wall only lowers it |
| `local_time` | C2, `B(t)` | m² | `Σ dℓ²` over the contact increments, with `ρ/D = 1` (§6.3). A point process, not Brownian; present iff `boundary_local_time` is |
| `field` | C3, `F(t) = ∫ b dt` | s² | one `Σ (b dt_f)²` per component of `susc_field_basis`, in the order of `compression.channels.susc_field_basis.channels`, over the field reads (below); smooth and small. Present iff `susc_field_basis` is |
| `field_increments` | C3, the read series `b_i` | — (dimensionless²) | one `Σ (b_{i+1} − b_i)²` per component of `susc_field_basis`, same order, over the field reads; it sets `fidelity.terms.field_read` (below). Present iff `susc_field_basis` is |

*The relative figure.* The block quadrature perturbs band `k` of a channel by an error independent of the others, whose variance is at most

```
rel(m, ω_k h) · ( T / (π k)² ) · QV(I),      ω_k = π k / T,
rel(m, x) = x^{2(m+1)} / ( ((m+1)!)² · (2m+3) · 4^{m+1} ),
```

with `T / (π k)² · QV` the band's variance under the free bridge (§9.5.6, where `QV = 2 D T`), scaled by `rel`:

| `m` | `rel(m, ω h)` |
|---|---|
| 0 | `(ωh)² / 12` |
| 1 | `(ωh)⁴ / 320` |
| 2 | `(ωh)⁶ / 16128` |
| 3 | `(ωh)⁸ / 1327104` |

It is an **absolute** error, not a fraction of the band's own variance: a restricted walk's high bands are far smaller than the free value while the quadrature error is not, so the same error is a large fraction of a restricted band and below the table in absolute terms.

*The quadrature term, per channel.* A measurement's phase (or log-weight) is `Σ_c Σ_k g_{c,k} β_{c,k}` plus endpoint terms, which are exact, with `g_{c,k}` the weight of the replay (§6) on band `k` of channel `c`: for C0, `γ Ĝ_k` per axis (§9.5.3); for C2, `ρ/D` times the band's projection of the acquisition's contact window; for C3, `γ B0 χ` times the band's projection of the gate. Block quadrature adds to that measurement an independent term of variance at most

```
v_c = Σ_k g_{c,k}² · rel(m, ω_k h) · ( T / (π k)² ) · QV_c ,
```

and the signal is biased by at most `½ v_c` times `|S|`. On a free walk (C0) this is at most `½ · rel(m, π K h / T) · 2 b D`, with `b` the b-value (s·m⁻²) of the measurement and `D` the largest free diffusivity of the walk: it does not depend on the walker count or on the substrate, only on the band and on the `b` the pack can reach (`replay_envelope.acquisition.b_max`). The certificate carries `v_c / 2` per channel (below). The quadrature term does not include the rounding of the accumulation, which a producer with `precision = "float32"` MUST measure and carry in `err_max` as the term `rounding` (below).

*The field read.* The field channel is accumulated at the walk's own position, at the step, and is never read from the stored bands: the bands are a band-limited path, and a read of a micron-scale field along them is a quadrature at the band's resolution. `walk_params.field_read_interval = dt_f` (s) records the interval between the walk's reads of the field; it is REQUIRED iff `susc_field_basis` is present. Between reads the walker moves `sqrt(2 D dt_f)`, and the trapezoid of two consecutive reads, `(b_i + b_{i+1}) dt_f / 2`, misses the path's deviation from the chord. Over a walk that deviation adds to the field phase an independent term whose variance is set by the field the walker ACTUALLY crossed, the quadratic variation of the read series, which the walk records per component in `walk_params.qv.field_increments`:

```
QV_Δb  = Σ_i (b_{i+1} − b_i)²            ≈ 2 D ∫ |∇b(x(t))|² dt ,
v_read = (γ B0 χ)² · dt_f² / 12 · QV_Δb ,
```

with `QV_Δb` combined over the channel's components as the square of the summed roots, `(Σ_c √QV_Δb,c)²` (the components are functionals of one path and the replay contracts them with weights of modulus at most one in any orientation, so this bounds every orientation), `B0` (T) and `χ` (dimensionless, SI) the field strength and susceptibility scale of the replay and the `1/12` the trapezoid's constant for a field linear over one read (the chord's curvature across a sheath raises it; a producer that measures a larger constant on a pilot states it). The gate bounds `|s| ≤ 1`. `QV_Δb` is a measured property of the walk, like the positions' `QV`: it costs one register per component (the previous read, which the trapezoid keeps anyway) and is independent of any assumed maximum of the field. `dt_f` is chosen from a pilot's `QV_Δb` per pool, the coarsest interval holding `v_read` below the square of the floor at the envelope's `B0` and `χ`; it cannot usefully be finer than the walk's step (a read along the chord between steps carries the same deviation), and when the floor would ask for less the read is the step and the envelope states the `B0 · χ` product under which the term stays below the floor (§9.4 rule 5). A merge of shards with their own intervals states the coarsest and charges the term at it.

*The certificate.* A pack has no raw walk to be decoded against, so its `fidelity` (§9.4 rule 4) reports its bounded errors as separate terms in `fidelity.terms`, each a bound on the replay error over the battery in the units of `err_max`, the largest over the battery:

| key | meaning |
|---|---|
| `truncation` | the bands beyond `K`. Walk-free, from the free-bridge band variances of §9.5.6 weighted by the battery waveforms' projections; or measured on a pilot subset as `K` against `2K` bands accumulated in the same walk, which `truncation_from` (`"band_variance"` or `"pilot"`, OPTIONAL) records. REQUIRED |
| `quadrature` | an object with one number per channel, `½ v_c`: `positions` always, `local_time` iff C2, `field` iff the field channel is present. Each is `½ Σ_k g_{c,k}² rel(m, ω_k h) (T/(π k)²) QV_c` with the `QV_c` of `walk_params.qv`. REQUIRED |
| `field_read` | `½ v_read` of the field read above, from `walk_params.qv.field_increments` at `walk_params.field_read_interval`. REQUIRED iff `walk_params.field_read_interval` is present |
| `rounding` | the rounding of the accumulation across blocks: with `r` the measured relative residual of the chain (the producer's `k = 0` mode against the exact displacement, the largest over the walkers), a phase error correlated with the phase, bounded by `r · sqrt(2 · b_max · D)` with `b_max` the battery's largest b-value and `D` the largest free diffusivity. REQUIRED iff `walk_params.block.precision` is `"float32"` |

`err_max` is at least the sum of the terms. Quantisation by a container is part of `err_max` as for any pack, and the split-half floors (`floor_max`, per voxel where partitioned) are read from the coefficients as ever. The terms are measured once against a per-step accumulation (h = the step) -- the slab for C0, a reflecting pore at the envelope's `ρ_max` for C2, a sheath for C3 -- and the signal against the forward engine within the floor; agreement of the bands alone is not the certificate.

#### 9.5.6 The free bridge's band variances

Free diffusion in one axis has displacement variance `2 D t`. Its bridge residual has covariance `2 D (min(s,t) − s t / T)`, whose eigenfunctions are the `s_k` and whose bands are uncorrelated, with variances

```
Var β_k = 2 D T² / (π k)² = ( T / (π k)² ) · QV ,        QV = 2 D T .
```

The phase variance a truncation at `K` leaves for a waveform is `γ² Σ_{k>K} Ĝ_k² Var β_k` per axis, with `Ĝ_k` as in §9.5.3. These are the free-walk values: the reference scale against which a restricted walk's measured band variance is read, and the scale of the quadrature term (§9.5.5). The temporal resolution of a pack is the band `K / (2 · T_max)` (§9.4 rule 6).

#### 9.5.7 Occupancy as an event list

The occupancy columns of C1 (§5.2) are piecewise constant in time and are stored exactly in time, never sampled. A column is stored in one of two forms, keyed in the registry:

- **static**: one label per walker (`comp_static`, `(N_w,)`), when no walker changes label over the walk. A producer MUST write this form when it applies and MUST NOT when any walker changes label.
- **events**: three tensors per column `<c>` (`comp` for the geometric axis, `bound` for the MT pool), `<c>_events_offsets` `(N_w + 1,)` int64, `<c>_events_t` `(N_ev,)` float64 and `<c>_events_label` `(N_ev,)` int8 (int16 where a pool id exceeds 127). Walker `i`'s events are the entries `offsets[i] … offsets[i+1] − 1`.

An event `(t_j, label_j)` says that the walker **enters** `label_j` at the time `t_j`, in **seconds from the walk's start**; the label is a pool id for `comp` and `0` (free) or `1` (bound) for `bound`. The first event of every walker is at `t = 0` and gives its initial label; the times of one walker strictly increase and lie in `[0, T_max]`; no event repeats the label before it. The walker holds `label_j` on `[t_j, t_{j+1})`, its last label to `T_max`. A crossing of a permeable wall, a binding and an unbinding are each one event, at the walk's own time, not at an instant of any grid of the pack. In a window table (§4.3) each window's list holds the events inside that window, with times still from the walk's start, and opens with an event at the window's start carrying the label held there.

A reader takes the label at an instant as that of the last event at or before it, and the integral of a rate over an interval as the sum of the rate of each label times its overlap with the interval -- exact. A producer that caps the events per walker MUST refuse the walk on overflow rather than drop an event.

---

## 10. Metadata schema

Metadata is a JSON object embedded in the container (§12) and validated by `schema/rpk_metadata.schema.json`. Required and optional keys:

```jsonc
{
  "rpk_schema_version": "0.4",          // REQUIRED, semver of the container schema
  "id": "canonical-wm/g070-f055-3T",    // REQUIRED, stable pack identifier
  "walk_params": {                       // REQUIRED
    "n_walkers": 120000,
    "T_max": 0.1,                        // s, the duration of the stored walk (§4.3)
    "diffusivity": 0.6e-9,               // m^2/s, FIXED (not a knob)
    "seed": 0,
    "block": { "h": 2.0e-5, "m": 1,      // REQUIRED (§9.5.5): block length (s), moment order,
               "precision": "float64" }, //   accumulator type
    "qv": { "positions": [1.2e-10, 1.2e-10, 1.2e-10],   // REQUIRED (§9.5.5): the quadratic variation per
            "local_time": 3.0e-14,                      //   channel, m^2 (positions, local_time) and s^2 (field);
            "field": [4.0e-14, 1.0e-14],                //   local_time iff C2, field iff susc_field_basis
            "field_increments": [2.0e-9, 5.0e-10] },     //   Σ (Δb)² per component, iff susc_field_basis
    "field_read_interval": 2.0e-6,       // s; REQUIRED iff susc_field_basis (§9.5.5)
    "boundary": "periodic",              // OPTIONAL "periodic" | "reflecting" | "finite";
                                         //   the walk's box (§4.5) + bounds the envelope
                                         //   (a finite box restricts long-diffusion-time replays);
                                         //   default "periodic"
    "substrate_frame": [[1,0,0],[0,1,0],[0,0,1]],  // OPTIONAL 3x3 orthonormal, right-handed frame
                                         //   R (§4.2): columns = canonical substrate basis in stored
                                         //   coords; column 3 (z) = primary axis (fibre), columns 1-2
                                         //   = a FIXED perpendicular basis (no free rotation about the
                                         //   fibre). A canonical direction maps to stored as R·d.
                                         //   Omit => identity (positions already canonical).
    "delta_chi_a": -0.1e-6,              // anisotropic susceptibility scale the Field basis
                                         //   are normalized to; present iff Field tier
    "field_grid": {                      // present iff `susc_field_basis` (§6.4.1); how it was built
      "origin": [0,0,0], "voxel_size": [2e-7,2e-7,2e-7],   // m
      "components": ["iso_local","iso_P","aniso_G"],        // order of the 13 basis channels
      "refocus": "gate" },               //   phase gate s(t_k) supplied at replay (§6.6)
    "windows": { "n": 1, "T": 0.1,       // OPTIONAL — the storage windows (§4.3): n windows of
      "walks": [ { "first": 0, "last": 0, "seed": 0 } ] }  //   T s each; the runs that walked them
  },
  "substrate": { ... },                  // RECOMMENDED — the substrate specification the walk was
                                         //   driven by (SUBSTRATE.md, one JSON document): defines the
                                         //   compartment ids (its pool ids), the walls, the seeding and
                                         //   the NOMINAL tissue values a replayer may choose to apply
  "compression": {                       // REQUIRED — codec descriptor (see CODEC_REGISTRY.md)
    "method": "bridge_dst", "K": 256, "walker_preserving": true,   // the bands of §9.5; no basis key
    "temporal_bandwidth_hz": 1280.0       // RECOMMENDED — K / (2 T_max): the band as a frequency (§9.4 rule 6)
  },
  "replay_envelope": {                   // REQUIRED — declared capabilities + domain
    "gradient": true,                    // C0 — always true
    "bulk_relaxation": true,             // C1 — intrinsic per-compartment T1/T2
    "surface_relaxivity": false,         // C2 — any surface relaxivity
    "field": true,                       // C3 — any B0 / susceptibility off-resonance and orientation
    "magnetization_transfer": false,     // C4 — magnetization transfer (vector-Bloch)
    "diffusivity_fixed": true,           // D is a fixed pack property, not a knob
    "acquisition": { "b_max": 3.0e9, "ogse_periods": [], "B0_list": [3.0, 7.0],
                     "field_modes": 32,              // K temporal modes kept for susc_field_basis
                     "max_refocusing_pulses": 16 }   // ~K/2 — the gate bandwidth K buys (§9.4 r5)
  },
  "fidelity": {                          // REQUIRED (§9.4 rule 4); codec-agnostic
    "err_max": 0.0015,                   // worst replay error over `battery`, against the reference replay
    "floor_max": 0.0074,                 // irreducible MC noise floor (same ensemble)
    "terms": {                           // REQUIRED (§9.5.5): bounds on the replay error, in err_max's units
      "truncation": 4.0e-4, "truncation_from": "band_variance",
      "quadrature": { "positions": 2.0e-6, "local_time": 1.0e-6, "field": 3.0e-7 },   // one per channel
      "field_read": 5.0e-5 },            //   iff field_read_interval
    "within_2x_floor": true,             // err_max <= 2*floor_max -> loss is negligible
    "battery": "...",
    "certified": "measured",             // "measured" here, or "inherited" from a certifying pack of the same fill
    "inherited_from": null,              //   then {id, err_max, floor_max} of that pack (§9.4 rule 4)
    "windows": null },                   // OPTIONAL — one certificate per storage window (§4.3); the
                                         //   whole's numbers are then the bound over them
  "mt": {                                // REQUIRED iff C4 parametric two-pool (§6.5.1); else absent
    "model": "two_pool",                 //   fixed-at-walk-time pool descriptor (not knobs)
    "f_bound": 0.081, "k_forward": 88.3, //   dimensionless, s^-1  (from S/V, kappa, dwell)
    "S_over_V": 4.4e6 },                 //   m^-1, the geometry the pool was derived from
  "provenance": {                        // RECOMMENDED — how the walk was made
    "generator": "dmipy-sim", "generator_version": "...",
    "substrate": {                       // WHAT was walked (descriptive; geometry format is out of scope, §1)
      "representation": "analytical",    //   "analytical" | "mesh" | "voxel" | ...
      "description": "PackedMyelinatedCylinders, gamma radii, f=0.55",
      "ref": null },                     //   optional URI / DOI / hash of the source substrate
    "permeability": 0.0,                 // the FIXED permeability the walk used (§7), if any
    "diffusivity_scaled": null,          // OPTIONAL: a time view (§6.8): {walked, replayed, a, T_max, permeability}
    "length_scaled": null,               // OPTIONAL: a space view (§6.8): {s, diffusivity, b_max: {walked, served}}
    "prefix": null },                    // REQUIRED on a re-encoded prefix (§3): the parent's id and content
                                         //   digest, its T_max and K, and this child's TE' and K'; null otherwise
  "license": "CC-BY-4.0",                // REQUIRED — the SOURCE substrate's license
  "pack_license": "CC-BY-NC-SA-4.0 OR LicenseRef-Commercial",   // OPTIONAL — the terms the PACK is distributed under
  "pack_license_url": "https://huggingface.co/datasets/<owner>/<name>/blob/main/LICENSING.md",   // OPTIONAL
  "citation": "Fick RHJ, dmrai-lab (2026) ..."   // REQUIRED
}
```

Rules: `pack_license` is the licence of the pack file itself, as an SPDX expression, and `pack_license_url` where its full text lives; absent, the pack is distributed under `license`, the source substrate's. `diffusivity` and `seed` are fixed pack properties, not knobs; a similarity view (§6.8) states the diffusivity it realises and records the walked one under `provenance`. `license` records the **source substrate's** license and MUST NOT relicense upstream geometry. Any tier flag set `true` in `replay_envelope` MUST have its channels present (§7) and its `walk_params` fields populated. C3 (`field`) is carried by the `susc_field_basis` channel (+ `walk_params.field_grid`, `field_read_interval`); C4 (`magnetization_transfer`) by the C1 `bound` occupancy column **or** the `mt` pool descriptor, either one satisfying the flag. The envelope flags are `gradient`, `bulk_relaxation`, `surface_relaxivity`, `field`, `magnetization_transfer` and `diffusivity_fixed`, each under that name only. The substrate bank derives its catalog card from `provenance` and the optional Croissant sidecar (§12); the card's layout is the bank's concern, outside this format.

---

## 11. Declared envelope and self-certification

A pack is a **certificate**, not a black box. Two objects make its guarantees explicit:

- **`replay_envelope`** — the *domain of validity*: which tiers, and the acquisition/field range (`b_max`, OGSE periods, `B0_list`, `max_refocusing_pulses`, …) over which the producer stands behind the pack. Its limits have two sources — what the *walk* covers (diffusion time, box size) and what a capability-narrowing codec leaves servable (§9.4 rule 5) — and a replayer treats them identically. A request inside the envelope is guaranteed; a request outside is a *known limit*, and the replayer MUST surface it as such (§13), never silently extrapolate.
- **`fidelity`** — the *certified* error of the pack's encoding against the noise floor (§9.4, §9.5.5).

The substrate bank (a separate service) accepts a pack by validating it against this specification (§13) and MAY re-measure `fidelity` independently.

---

## 12. Container format

- **Arrays:** a single **safetensors** file. Safetensors is REQUIRED because it (a) carries **no executable code** — a pack is safe to download and memory-map from untrusted sources — (b) is strongly typed, and (c) is zero-copy. A channel's tensor keys are those the registry defines for its method (`CODEC_REGISTRY.md`).
- **Metadata:** the JSON object of §10, serialized as a single string and stored in the safetensors `__metadata__` header under the key **`"rpk"`**. Producers MUST write it there and readers read it from there; a sibling `.json` copy is OPTIONAL for human inspection.
- **Dataset sidecar (OPTIONAL):** a Croissant (`schema.org/Dataset` + MLCommons) `.croissant.jsonld` carrying `license`, `citation`, `provenance`, and `replay_envelope` for dataset-catalog interoperability.
- **Extension:** `.rpk`. **Naming:** the file basename SHOULD equal the last path segment of `id`.

A pack MUST contain no code and MUST be openable without executing anything.

---

## 13. Conformance

A file is a **conformant Replay Pack** iff:

1. it is a valid safetensors file with a `__metadata__` block that validates against `schema/rpk_metadata.schema.json`;
2. it contains the REQUIRED `positions` channel as the band tensors the registry defines (`pos_x`, `pos_y`, `pos_z`, each of width `K + 2`);
3. every tier flagged `true` in `replay_envelope` has all its §7 channels present and its §10 fields populated;
4. all per-walker channels share `N_w` with `positions`, and every channel the duration `T_max`;
5. `walk_params.block`, `walk_params.qv` and `fidelity` with `fidelity.terms` (§9.5.5) are present;
6. `license` and `citation` are present;
7. it carries none of the forms §9.5.1 refuses.

A **conformant replayer**:

- implements §6.1 and at least refuses cleanly for tiers it does not implement;
- decodes every codec it claims to support to the §5 channel contract before applying §6;
- refuses a pack carrying a form of §9.5.1 by name, with the instruction to re-walk it;
- for any requested knob outside the pack's `replay_envelope`, returns a **capability/domain error**, never a silent approximation (§11).

A **conformant producer** satisfies §8 and §9.4.

Conformance is layered: a producer MAY implement only C0; a replayer MAY implement only C0; they still interoperate on C0.

---

## 14. Versioning and extensibility

- `rpk_schema_version` is semantic. **Minor** bumps add OPTIONAL channels/keys; a replayer MUST ignore channels and metadata keys it does not recognize. **Major** bumps may change required semantics.
- New channels intended for later standardization SHOULD be introduced under the `x_` prefix first.
- New tiers extend the table in §7 without altering existing tiers.
- **Recognized extensions** (permitted by the invariant, not yet standardized):
  - **Transmit-field (`B1+`) inhomogeneity.** Within a single substrate `B1` is uniform — it varies on the RF-wavelength (cm) scale, far larger than a μm–mm pack — so it is **not** a stored map but a **scalar flip-angle knob** `α → κ_B1·α` on the vector-Bloch path (§6.5), the transmit-side analog of the `B0` *field-strength* scalar (not of the susceptibility map). Uniform pulses are already replay knobs (§6.6); the knob adds only the inhomogeneity scale. When packs are tiled into a **mesoscopic voxel grid** (a brain-sized volume), `B1` inhomogeneity becomes a **per-voxel scalar map** — one `κ_B1` per voxel position, applied uniformly inside each voxel's substrate — which is a property of the *assembly*, not of a single pack.
  - a fully general per-voxel rank-2 **susceptibility-tensor** field (§5.2, §6.4). The `0.2.0` `susc_field_basis` (§6.4.1) already covers an isotropic source plus a single rank-2 anisotropic director (the myelin case) for arbitrary morphology; the open extension is a spatially-varying full tensor;
  - membrane **permeability/crossing**, should it ever become a replayable quantity (§7).

---

## 15. Security considerations

A pack carries data only. Safetensors executes no code on load, and the metadata and Croissant sidecar are inert JSON. Replayers MUST NOT `eval`/`exec` any pack content. Field-basis and band arrays are bounded numeric data; a replayer SHOULD sanity-check shapes against `walk_params` before allocating.

---

## 16. References

- RFC 2119 — Key words for requirement levels.
- safetensors — https://github.com/huggingface/safetensors
- Croissant (MLCommons) — http://mlcommons.org/croissant/
- This document is the citable, versioned reference for the format; its Zenodo DOI ([10.5281/zenodo.21692505](https://doi.org/10.5281/zenodo.21692505)) establishes priority on the format and infrastructure.

---

## Appendix A. Minimal conformant pack (Tier C0)

Arrays: `pos_x`, `pos_y`, `pos_z`, each `(N_w, K + 2)` float32.
Metadata: `rpk_schema_version`, `id`, `walk_params{n_walkers,T_max,diffusivity,seed,block,qv}`, `compression{method:"bridge_dst",K,walker_preserving:true}`, `replay_envelope{gradient:true, ...all others false...}`, `fidelity{err_max,floor_max,terms{truncation,quadrature{positions}}}`, `license`, `citation`.
Replay: §6.1 only. This is the floor of interoperability — a walk and nothing else — and every richer pack is this plus declared channels. See `examples/`.
