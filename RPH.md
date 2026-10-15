# The Replay Phantom Specification (`.rph`)

**Status:** draft, versioned independently of the core `.rpk` specification.
**Version:** 0.5.0 (draft). `0.5.0` adds the second route to a phantom, **partition** (§9): one walk of a substrate larger than a voxel, cut into voxels by where each walker started, with membership derived rather than stored, fractions emergent from the walkers' weights, a grid attached to the tissue or to the bore, and a rigid pose. `0.4.0` says what a pose is -- a **rotation**, not an axis -- and adds the two modes that follow from it: `frames`, a rotation per slot, and `bingham`, a frame with two concentrations, so an anisotropically dispersed population (a fan) is expressible. `0.3.0` places the grid in the scanner (`grid.origin_m`, `grid.isocenter_m`, the meaning of `frame`), turns `scalars` into a registry of **macroscopic layers** with stated replay semantics (§5.1), lets a phantom declare the tissue values a pack substrate replays at (§3.2), and fixes the conformance rule to what §3 always said: fractions sum to **one**. `0.2.0` added the ODF mode and the analytic substrate.
**License:** CC-BY-4.0 (text) / Apache-2.0 (reference code).

A **replay phantom** is a spatial arrangement of solved substrates -- an assembly of replay
packs ([`RPK.md`](RPK.md)) -- together with the per-voxel information needed to compose them
into a signal. It stores **no walkers and no trajectories of its own**. Everything physical
lives in the packs; the phantom adds only *which substrates are where, how they are oriented,
in what proportion, and with what proton density*.

This is a **secondary** format. It extends nothing in the replay invariant (SPEC §3) and adds
no capability tier (SPEC §7). A phantom is exactly as replayable as the packs it carries, and
its tier is the **intersection** of theirs.

## 1. Why a phantom is its own format

There are **two ways to arrive at a phantom**, and this document specifies both. **Composition** (§3-§7) places
solved packs into voxels the producer declares: fractions, poses and proton density per voxel, on a grid that
is the producer's. **Partition** (§9) bins the walkers of one walk of a substrate larger than a voxel into the
voxels they started in: the grid is free, the poses are the tissue's, the fractions are measured. At the pack
layer the two are identical -- walkers over tiers, compressed -- and walker provenance is invisible to the
replay mechanics, so partition adds no capability tier and changes nothing in the replay invariant. Composition
is the voxel-averaged projection of a partition; the projection is not invertible, which is why both exist.

A pack answers *what does this microstructure do to the magnetization*: one substrate at one
pose. A voxel is a distribution of poses, a volume is a field of such distributions, and both
are properties of the arrangement rather than of any substrate in it. SPEC §14 already draws
that line -- transmit-field inhomogeneity becomes a per-voxel scalar map once packs are tiled,
belonging to the arrangement and not to any pack.

The split earns itself three times over. The expensive object stays shared: one solved pack
serves every voxel and every orientation that cites it, so a whole-brain phantom over a handful
of tissue types is small next to the packs it draws on. Packs and phantoms then evolve
independently -- re-solving a substrate touches no phantom, re-arranging tissue re-solves
nothing. And provenance stays honest, because a phantom names what it is made of rather than
absorbing and anonymising it.

## 2. Referenced and embedded substrates

A phantom cites each substrate either **by reference** or **embedded**, and MAY mix the two.

*Referenced* keeps the phantom small and is right for a bank, where the packs are already
published and addressable. *Embedded* makes the phantom a **standalone artifact**: a brain
phantom with its white-matter, grey-matter and CSF substrates inside it can be shared, cited
and replayed as one file, with no resolution step and nothing to go missing. That is the mode
to prefer for anything archival, since a reference is only as durable as what it points at.

Embedded packs are stored **as tensors, not as opaque blobs**: substrate `i`'s arrays appear
under the prefix `substrate{i}/`, so `substrate0/pos_x` is the `pos_x` channel of the first
substrate. This keeps the properties SPEC §12 requires safetensors for -- strong typing and
zero-copy access -- and lets a reader memory-map one channel of one substrate without
materialising the rest. The prefix also removes any collision between substrates that share
channel names, which all of them do.

Each embedded pack's own metadata object is carried verbatim under `substrates[i].pack_meta`.
A reader MUST treat it exactly as it would that pack's `"rpk"` header: the embedding changes
where the bytes live, never what they mean. `sha256` is REQUIRED in both modes and pins the
identity of the solved physics either way.

## 3. Data model

One safetensors file; SPEC §12 conventions apply unchanged. The grid is stored **sparsely** --
only occupied voxels appear -- so an anatomy that fills a fraction of its bounding box costs
only what it occupies.

| Array | Shape | dtype | Meaning |
|---|---|---|---|
| `voxel_index` | `(N_v, 3)` | int32 | voxel coordinates on the grid |
| `substrate_id` | `(N_v, P)` | int16 | index into `substrates` (§5); `-1` marks an unused slot |
| `geometric_fraction` | `(N_v, P)` | float32 | fraction of the voxel volume occupied by that substrate |
| `peak_dir` | `(N_v, P, 3)` | float32 | *peaks mode*: unit direction of each slot, its azimuth unstated |
| `odf_sh` | `(N_v, P, n_c)` | float32 | *ODF mode*: distribution over directions, even-order real SH |
| `pose_quat` | `(N_v, P, 4)` | float32 | *frames / bingham mode*: the slot's rotation, `(x, y, z, w)`, unit |
| `bingham_kappa` | `(N_v, P, 2)` | float32 | *bingham mode*: concentrations about the frame's first two axes |
| `roll_kappa` | `(N_v, P)` | float32 | *bingham mode*, OPTIONAL: concentration of the azimuth about the frame |
| `scalars` | `(N_v, S)` | float32 | OPTIONAL per-voxel scalars, named in metadata |

`P` is the number of slots per voxel. A voxel using fewer pads with `substrate_id = -1` and
zero fraction.

**`geometric_fraction` is geometry, not signal.** It is the share of the voxel volume the
substrate occupies, before any relaxation or proton density is applied. Rows MUST sum to
**one**, and a replayer MUST reject a row that does not, rather than normalising or accepting
it.

A voxel is always full. There is no vacuum in a sample, so a row summing to less than one is
not a voxel with a void in it -- it is a voxel whose remainder was not modelled, and composing
it returns a signal that is quietly too low while looking entirely legitimate. Anything that is
not tissue is therefore declared as a substrate rather than left as slack. §3.1 gives the
three kinds a substrate may be.

Partial volume needs no slack in the sum, because it is carried by the *ratio* of the
fractions. A voxel that is 60% white matter and 30% grey matter with the remainder outside the
modelled object is `0.6 / 0.3 / 0.1-inert`; one made only of those two tissues is `2/3 / 1/3`.
Both sum to one, and they are different voxels -- which is exactly what a boundary needs to
express.

### 3.1 Three kinds of substrate

| `kind` | Response from | Needs a pack | Needs an orientation |
|---|---|---|---|
| `pack` | stored walkers (`RPK.md`) | yes | yes |
| `analytic` | a declared closed form | no | only if the model is anisotropic |
| `inert` | nothing | no | no |

**`analytic` exists because for some substrates a pack is not merely wasteful but unusable.**
Free water is the case. Its signal decays exponentially in `b`, while the Monte-Carlo error of
a walker ensemble decays only as `1/sqrt(N_w)` and is roughly `b`-independent, so the *relative*
error grows like `exp(+bD)/sqrt(N_w)`. At `D = 3.0e-9 m^2/s` with 4000 walkers the noise floor
is 0.1x the signal at `b = 1000 s/mm^2`, **113x** at 3000 and **8.5e4** at 5000; reaching 1%
relative accuracy would need of order `5e11` walkers at `b = 3000` and `3e17` at 5000. No pack
can be built. The closed form `exp(-bD)` is exact, has no walkers to store, and -- being
isotropic -- carries no orientation, so it also needs no ODF.

This version defines one analytic model:

```jsonc
{"id": "csf/free-water", "kind": "analytic", "m0": 1.00,
 "model": "free_water", "params": {"diffusivity": 3.0e-9}, "T2_s": 2.0}
```

with response `E = exp(-b D) · exp(-TE / T2_s)`, independent of gradient direction and of `B0`
(`T2_s`, and likewise `T1_s` over a mixing time, are OPTIONAL: absent, no relaxation applies). A replayer MUST
refuse an analytic `model` it does not recognise rather than guessing, exactly as SPEC §9
requires for codecs. Further models are additions to this table, not changes to the format.

**Namespaced models.** A `model` of the form `"<package>:<Name>"` is a closed form another package
defines (a compartment model of dmipy-fit: `"dmipy_fit:C1Stick"`, parameters by the model's own names
and units under `params`). A reader resolves it by importing `<package>.phantom` and calling its
`analytic_substrate(entry)`; without that package it MUST refuse, naming the package, never guess.

**A closed form with an axis takes a pose.** An analytic entry that carries `"oriented": true` has an
axis (a stick, a cylinder) and is composed exactly as a pack is (§6): its response at a pose -- the
form evaluated with its axis along the third column of the pose -- is expanded over SO(3) and contracted
with the voxel's orientation distribution, so it MUST have an orientation field like a pack and MUST NOT
disperse itself (a model that carries its own dispersion is not citable: the phantom's field would
disperse it twice). An entry without the flag (free water) is isotropic and MUST NOT be given one.
A closed form is **full-tier**: every knob a pack takes has an exact value for it. A form with no
susceptibility source has a field of zero at any `B0`; with no wall it has no surface relaxivity; with no
bound pool no magnetisation transfer; under an RF train it is a static spin's response. What it carries it
evaluates exactly: free water's bulk relaxation, `exp(-TE / T2)` (and `exp(-TM / T1)`), when the entry
declares `T2_s` (`T1_s`). A replayer MUST NOT warn about, skip, or fill in physics for a closed form: the
zeros are the physics.

**`inert` is not air.** It is defined by contributing nothing, which is a modelling statement,
not a material. Air is the opposite of inert magnetically: the air--tissue susceptibility step
is of order 9 ppm, roughly two orders of magnitude larger than the sub-ppm anisotropy the packs
carry, and it perturbs the field of *neighbouring* voxels over centimetres. A phantom composes
packs voxel-by-voxel and has no mechanism for that inter-voxel field, so a genuine air cavity
MUST NOT be represented by an inert substrate -- doing so would put an air interface in the
geometry and silently omit the dominant effect it has. Representing air needs a macroscopic
`B0` field map over the grid, which is an assembly property in the sense of SPEC §14 and is not
in this version.

### 3.2 What a pack substrate replays at

A pack carries the substrate specification it was walked from (`RPK.md` §10) and no tissue value; the values a
replay applies are knobs. A phantom fixes them per substrate: `substrates[i].tissue`, when present, is an
object of the knobs of a pack's replay -- per-pool `T2` / `T1` (s), the wall relaxivity `rho` (m/s), the field
source's `chi_iso` / `chi_aniso` -- and MUST be applied by the replayer to that substrate; a knob not listed
takes the pack's **nominal** value (the embedded specification's declared value, its `nominal_field_T` as
`B0` when the replay gives none). The override is part of the phantom's declaration and travels with the
file: the same phantom replays the same way everywhere.

## 4. Orientation: a pose is a rotation

Exactly one orientation mode is declared per phantom. What they have in common is the thing to state first: the
pose of a substrate is a **rotation**, not a direction. Naming the direction a substrate points along fixes two
of its three degrees of freedom and leaves the spin about that axis unstated, and a substrate's response depends
on that spin unless the substrate happens to be axially symmetric -- which a finite bundle of tortuous strands
is not, and neither is a fanned population. The four modes differ in how much of the rotation they pin down.

**Peaks** (`peak_dir`) give each slot a direction and leave the azimuth about it unstated. A replayer MUST then
integrate over that azimuth rather than choose a value for it: leaving it to a convention would make the signal
depend on a private choice of the producer's frame. This is the representation for discrete crossings -- a voxel
with two fibre populations is two slots with two directions and two fractions -- and `P <= 3` covers the
configurations resolvable in practice.

**ODF** (`odf_sh`) gives each slot a distribution over directions in the even-order real spherical harmonic
basis, for dispersion and fanning about a mean direction, again with the azimuth unstated and integrated away.
Peaks are the zero-dispersion limit of this mode.

**Frames** (`pose_quat`) give each slot a whole rotation, as a unit quaternion `(x, y, z, w)`. This is what a
substrate whose response is not axially symmetric needs in order to be placed unambiguously, and what any
operation acting on the magnetisation vector needs -- an RF pulse at a scaled flip angle is applied to a pose,
not to an axis.

**Bingham** (`pose_quat` + `bingham_kappa`, optionally `roll_kappa`) gives each slot a frame and two
concentrations, one about each of the frame's first two axes: a population dispersed **anisotropically**, wide
in one plane and narrow in the other. Equal concentrations are a Watson cone of the same width, so this mode
contains the isotropic case rather than replacing it. `roll_kappa`, when present, concentrates the substrate's
own azimuth about the frame; absent, that azimuth is free and is integrated away as in the modes above.

A conformant replayer MUST produce the same signal from a peak set as from ODFs concentrated on those
directions with the same weights, to within the SH truncation, and the same signal from a `bingham` slot with
equal concentrations as from the matching Watson ODF -- acceptance tests, not remarks. The modes are stored
separately because each states exactly what it knows: evaluating a response at a pose is cheaper and more exact
than contracting a near-singular distribution, and a mode that pins the azimuth cannot be recovered from one
that does not.

In peaks and frames mode a slot's `geometric_fraction` is the weight of that population, so a voxel with two
populations of one substrate is two slots citing the same substrate with two poses and two fractions.

### 4.1 The spherical-harmonic basis is normative

Phantoms will be built from ODFs produced elsewhere -- MRtrix, DIPY, dmipy -- and those tools
do not share a real spherical-harmonic convention. The differences are not cosmetic here: the
composition of §6 goes through the SH addition theorem, which holds only for an **orthonormal**
basis, and one widely used convention is not orthonormal. `odf_sh` MUST therefore be stored in
the basis below, and `orientation.convention` MUST name it.

**The required basis.** Orthonormal real spherical harmonics, even orders only, one contiguous
block per order with `m` ascending from `-l` to `+l`, `Y_{l,0} = sqrt((2l+1)/4pi) P_l(cos t)`
and `Y_{l,±m}` the `sqrt(2)`-scaled cosine (`+m`) and sine (`-m`) terms. This is DIPY's
`real_sh_tournier(..., legacy=False)`.

**Converting from the common alternatives.** Both relations below are exact and per-coefficient;
neither is a resampling.

| Source | Orthonormal | Coefficient conversion to the required basis |
|---|---|---|
| `tournier` (DIPY, `legacy=False`) | yes | identity |
| `mrtrix` / `tournier_legacy` (DIPY default; MRtrix `.mif` FODs) | **no** | `c ← c/sqrt(2)` for `m ≠ 0`; `m = 0` unchanged |
| `descoteaux` (DIPY) | yes | `c_{l,m} ← s_m · c_{l,-m}`, with `s_m = (-1)^m` for `m > 0` and `+1` otherwise |

The two failure modes differ in how loudly they fail, which is the reason to be strict. The
MRtrix basis differs by a **scale** on `m ≠ 0` and breaks the addition theorem, so a phantom
imported without conversion is wrong by an amount that **vanishes exactly when the gradient is
parallel to `B0`** -- the one geometry a cursory check would test. The Descoteaux basis is
orthonormal but is a different basis, related by a signed `m → -m` permutation within each
band; unconverted it is wrong everywhere, which at least announces itself.

A producer that cannot establish which convention its source used MUST NOT declare one. There
is no safe default: guessing `tournier` for MRtrix output silently rescales every `m ≠ 0`
coefficient.

## 5. Proton density

Every substrate carries its own **`m0`** -- equilibrium proton density, in whatever units the
phantom declares -- in `substrates[i].m0`. This is required, not optional: white matter, grey
matter and CSF differ in proton density by tens of per cent, and a phantom that omits it
silently asserts they do not. It is a property of the substrate as used here, so it lives with
the substrate rather than in a per-voxel array; genuine spatial variation of proton density
within one tissue belongs in `scalars`.

### 5.1 Macroscopic layers: the `scalars` registry

The packs carry the microstructure. What a scanner adds on top of it varies over centimetres, not microns,
and is a property of the voxel, not of the substrate: a transmit field, a macroscopic off-resonance, a
proton-density gradient. These are **per-voxel scalars**, `scalars[v, s]`, named in `metadata.scalars`,
each name in the registry below with its meaning in the replay. A replayer MUST refuse a name it does not
know rather than ignore it: a layer that is silently dropped is a phantom that replays wrong while looking
right. Layers are added to a phantom one at a time; a phantom that declares none composes exactly per §6.

| name | unit | meaning | how it enters the replay |
|---|---|---|---|
| `kappa_B1` | – | transmit (B1+) scale at the voxel, 1 = nominal | multiplies every RF flip angle of the acquisition; needs the RF-aware (vector-Bloch) replay of each cited pack at the slot's pose, so it applies in peaks mode and to a magnitude-only gradient replay it MUST be refused, not dropped |
| `delta_B0_T` | T | macroscopic off-resonance at the voxel (a field map value), relative to the nominal `B0` | a uniform precession `gamma * delta_B0_T` over the voxel: a phase `gamma delta_B0_T int s(t) dt` with `s` the acquisition's coherence gate (zero for a spin echo whose 180 sits at TE/2); it does not dephase within the voxel -- the intra-voxel gradient of a field map is an acquisition-side term and is not this layer |
| `m0_scale` | – | proton-density variation within one tissue, 1 = the substrate's `m0` | multiplies `m0` of every slot of the voxel |

A layer that needs the field of *neighbouring* voxels -- an air or bone interface, a susceptibility step at a
tissue boundary -- is not a per-voxel scalar; it is an assembly property (`RPK.md` §14) and is out of this
version, exactly as §3.1 says of air.

## 6. Replay operation

For voxel `v` and acquisition `q`, with `E_i` the response of substrate `i` and `F_{v,p}` the
orientation of slot `p`,

```
S_v(q) = sum_p  geometric_fraction[v,p] * m0[substrate_id[v,p]]
                * INT_{S^2} F_{v,p}(n) E_{substrate_id[v,p]}(n; q) dn
```

with the integral replaced by `E(peak_dir[v,p]; q)` in peaks mode, and by the closed form
itself for an `analytic` substrate whose model is isotropic -- there the orientation carries no
information and MUST be ignored rather than applied. Where the response depends
on the field direction as well as the gradient -- susceptibility -- the two axes MUST be
composed jointly; the reduction to independent one-dimensional convolutions does not hold
there.

Composition is linear in the orientation, in `geometric_fraction` and in `m0`, so a phantom
introduces no new physics and cannot supply a tier its packs lack.

## 7. Metadata

JSON under the safetensors header key **`"rph"`**:

```jsonc
{
  "rph_schema_version": "0.2.0",
  "id": "phantoms/brain/hcp-like-1mm",
  "grid": {"shape": [180, 216, 180], "voxel_size_m": [1e-3, 1e-3, 1e-3],
           "origin_m": [-0.0895, -0.1075, -0.0895],       // scanner coordinates of the CENTRE of voxel (0,0,0)
           "isocenter_m": [0.0, 0.0, 0.0],                 // where the scanner is focused (default: the grid centre)
           "frame": "RAS",                                 // grid axes in the scanner: i -> +x (R), j -> +y (A), k -> +z (S)
           "attach": "substrate"},                         // what the grid is welded to under a pose (§9.2); inert for a composition
  "orientation": {"mode": "peaks", "max_peaks": 3},   // or {"mode": "odf_sh", "lmax": 8,
                                                      //      "basis": "real",
                                                      //      "convention": "orthonormal"}
                                                      // or {"mode": "frames"} / {"mode": "bingham"}
  "substrates": [
    {"id": "canonical/wm/g070-f055", "kind": "pack", "m0": 0.70, "embedded": true,
     "sha256": "…", "pack_meta": {…},                 // arrays under substrate0/
     "tissue": {"T2": [0.055, 0.05, 0.01], "chi_iso": -0.1e-6}},   // OPTIONAL knobs this substrate replays at (§3.2)
    {"id": "canonical/gm/…", "kind": "pack", "m0": 0.85, "embedded": true,
     "sha256": "…", "pack_meta": {…}},
    {"id": "csf/free-water",   "kind": "analytic", "m0": 1.00,
     "model": "free_water", "params": {"diffusivity": 3.0e-9}},
    {"id": "background/inert", "kind": "inert", "m0": 0.00}
  ],
  "scalars": ["kappa_B1"],                              // names from the registry of §5.1, columns of `scalars`
  "license": "…", "citation": "…", "provenance": {…}
}
```

`grid.frame` names the scanner axes the grid indices run along; the acquisition's gradient directions and
`B0` direction are given in the scanner frame and a replayer rotates them into each slot's orientation through
it. `origin_m` and `isocenter_m` are what make a macroscopic layer a function of position in the bore; a
replayer MUST default a missing `isocenter_m` to the grid centre and MUST NOT default a missing `origin_m`
when any layer is declared.

In ODF mode `orientation.convention` MUST be `orthonormal`. The composition goes through the
spherical-harmonic addition theorem, which is false for the non-orthonormal real conventions in
common use, and the error it introduces vanishes exactly when the gradient is parallel to `B0`
-- the one geometry a cursory check would test.

## 8. Conformance

A file is a conformant `.rph` when it is a safetensors container carrying the arrays of §3 and
the metadata of §7; every `substrate_id` resolves; `geometric_fraction` rows sum to **one**
(§3); every declared scalar is in the registry of §5.1; every substrate declares `m0` and `sha256`; embedded substrates carry their arrays under
`substrate{i}/` and their `pack_meta`; and, in ODF mode, `odf_sh` is in the declared orthonormal
basis.

A conformant replayer resolves or reads every cited substrate, refuses (never guesses) one it
cannot, composes per §6 over both axes where susceptibility is present, integrates the azimuth away in the
modes that leave it unstated rather than choosing a value for it (§4), reproduces a peak set from concentrated
ODFs and a Watson from an equal-concentration `bingham`, applies every declared layer as §5.1 states or refuses
the acquisition that cannot carry it, applies `substrates[i].tissue` where given, and reports the phantom's tier
as the intersection of its substrates' tiers.

A partition file (§9) is conformant when it carries the metadata of §9.4 and, for every declared slot, its
`declared{i}/` tensors; it carries **no** `voxel_index`, `substrate_id`, `geometric_fraction` or orientation array
for its pack slots, since a reader derives them. A conformant replayer of a partition derives membership from
the stored `r(0)` (§9.1), verifies any membership cache against that derivation rather than trusting it, takes
pack fractions from the weights (§9.3), applies the pose with the declared attachment (§9.2), and never applies a
per-voxel orientation to a partitioned pack (§9.2).

## 9. Partition: one walk, cut into voxels

A **partitioned** phantom cites one walk of a substrate larger than a voxel -- an anatomy, a strand phantom of
a cubic millimetre -- and assigns each walker to the voxel its **starting position** fell in. Everything
below is an *addressing mode* on top of the pack layer: nothing new is stored per walker, no tier is added.

### 9.1 Membership is derived

`voxel(i) = floor((pose · r_i(0) − corner) / voxel_size)`, with `r_i(0)` the pack's stored start position
(`RPK.md`: the position codec keeps the two endpoints exactly, so `r(0)` is a stored value and identical at
every `K`; measured, membership changes for 0.0000 % of 20,000 walkers between `K = 16` and `K = 64`, the
residual being float32 rounding of the coordinate), `corner` the low corner of voxel `(0, 0, 0)`
(`origin_m − voxel_size_m / 2`), and `pose` the rotation of §9.2 when the grid is the bore's, the identity when
it is the tissue's. A walker outside the grid belongs to no voxel. Membership MUST NOT be stored as a
normative array; a producer MAY store a cache (`walker_voxel (N_w,) int32` per pack), and a reader MUST verify
it against the derivation rather than trust it.

Because membership is derived, **the grid is free**: the same walk yields any voxel size, any shift, any
parcellation, without a re-walk or a re-encode -- refinement included, which a composition cannot do, having
stored only voxel summaries.

### 9.2 Rigid scope, and what the grid is attached to

A partitioned anatomy is one physical object. Its voxels are artificial subdivisions of it and cannot be posed
one by one, so orientation has **rigid scope**: one rotation, `pose.rotation` (substrate frame → scanner
frame), for the whole phantom. Per-voxel `peak_dir` / `odf_sh` MAY be carried as *descriptive* ground truth
(`orientation.role: "descriptive"`) and MUST NOT be applied in the replay: the trajectories already contain
the local orientation, dispersion and crossings, and applying an orientation operator would count them twice.

"Rotating the phantom" is two independent things, and `grid.attach` says which:

| `grid.attach` | the grid is welded to | a pose changes |
|---|---|---|
| `"substrate"` | the tissue | the physics only: the acquisition is rotated into the tissue frame (`g · R r = (Rᵀ g) · r`); membership is invariant |
| `"lab"` | the bore | the physics **and** the binning: the posed starts are rebinned, fractions re-emerge |

The two coincide at the identity pose. A partition with **no grid of its own** (`grid: null`) is bore-attached
by definition and takes its voxels from the acquisition's prescription (`ACQUISITION.md` §3.6) at replay time;
a replayer MUST refuse an acquisition without one. `attach` is inert for a composed phantom, which never poses.

Translation collapses to one degree of freedom -- only the relative offset of specimen and grid matters -- and
is the grid's (`origin_m`); `pose` carries a rotation only.

### 9.3 Fractions are emergent

A pack slot of a partition declares **no** `geometric_fraction`: under a bore-attached grid and a non-identity
pose a stored fraction is only valid at the pose it was measured at, and recomputing it needs the geometry,
which no pack carries. Instead each walker carries its statistical weight `spin_weights` (`RPK.md` §5.2), an
importance weight `(true volume density) / (sampling density)` of its compartment, and for voxel `v`

```
W_p(v) = sum_{i in v, pack p} w_i          f_p(v) = (1 − d(v)) · W_p(v) / sum_q W_q(v)
```

with `d(v)` the sum of the **declared** fractions of the voxel's walker-less slots (a myelin the walk excluded)
-- the only fractions a partition states, because nothing weighs them. A voxel with no walkers is the
**outside** substrate (`analytic` or `inert`, §3.1) at `1 − d(v)`, or is absent from the phantom when none is
declared. Rows sum to one by construction; §3's rule becomes a check on the derived quantity.

The trap this avoids is worth stating: with stratified seeding, **unweighted counts are not fractions**, and
a phantom inferring them from counts is wrong in proportion to the stratification. With the weights, weighted
counts *are* the fractions, at any grid and any pose.

### 9.4 Metadata and tensors

```jsonc
{
  "rph_schema_version": "0.5.0",
  "addressing": "partition",
  "grid": {"shape": [40, 40, 40], "voxel_size_m": [25e-6, 25e-6, 25e-6], "origin_m": [...], "isocenter_m": [0, 0, 0],
           "frame": "RAS", "attach": "substrate"},          // or null: the voxels are the acquisition's prescription
  "pose": {"rotation": [[1, 0, 0], [0, 1, 0], [0, 0, 1]]},  // substrate frame -> scanner frame
  "orientation": {"mode": "rigid", "scope": "rigid"},
  "substrates": [
    {"id": "disco/intra", "kind": "pack", "m0": 1.0, "addressing": "partition", "embedded": true, "sha256": "…", "pack_meta": {…}},
    {"id": "disco/extra", "kind": "pack", "m0": 1.0, "addressing": "partition", "uri": "disco_extra.rpk", "sha256": "…"},
    {"id": "myelin", "kind": "inert", "m0": 0.0, "declared": true},           // tensors declared2/voxel_index (N_d, 3), declared2/fraction (N_d,)
    {"id": "csf/free-water", "kind": "analytic", "m0": 1.0, "model": "free_water", "params": {"diffusivity": 0.6e-9}, "outside": true}
  ],
  "license": "…", "citation": "…"
}
```

Several pack slots citing **the same walk** (one pack per compartment, each with its weights) share every
voxel through §9.3; a producer MUST make them one walk, since binning two independent walks into the same
voxels is a composition wearing a partition's name.

### 9.5 Replay

For voxel `v`, with `E_i(q)` the complex signal of walker `i` under acquisition `q` in the tissue frame
(the acquisition rotated by `poseᵀ`), `ew_i` its weight with the relaxation and surface terms of the tiers
applied, and the outside and declared slots as in §6,

```
S_v(q) = sum_p m0_p f_p(v) · [ sum_{i in v, p} ew_i E_i(q) ] / W_p(v)  +  sum_{declared, outside} m0_s f_s(v) E_s(q)
```

Every pack is replayed **once** per acquisition -- the sum over walkers is regrouped by voxel, not recomputed
-- so the cost is the walk, not the voxel count. A macroscopic layer given per voxel (§5.1) is looked up
through membership and applied **per walker**: a transmit scale becomes a per-walker `B1+` scale through one
propagation of the whole walk. A walker's signal is attributed to the voxel it **started** in, not the one it
is in at the readout; walkers leave their voxel, and that is what makes the slicing exact.

### 9.6 The two routes agree, and where they part

A partition and its composition projection -- the fractions measured per voxel, an orientation distribution
fitted per voxel, one pack per tissue cited at it -- MUST agree at the grid the projection was taken on, to
within the packs' Monte-Carlo floor. They part wherever a sub-voxel or grid-changing question is asked:

| operation | composition | partition |
|---|---|---|
| coarsen the grid | merge slots; refit once `P` binds | exact rebin |
| refine the grid | impossible | exact rebin |
| shift the grid | resampling | exact rebin |
| rigid rotation, bore-attached grid | resampling | exact: rotate `r(0)`, rebin |
| reuse one pack across many voxels and poses | yes -- the whole-brain economy | no: the anatomy must be walked |

The last row is why composition is not merely a lossy partition, and why both routes are in this document.

