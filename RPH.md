# The Replay Phantom Specification (`.rph`)

**Status:** draft, versioned independently of the core `.rpk` specification.
**Version:** 0.2.0 (draft).
**License:** CC-BY-4.0 (text) / Apache-2.0 (reference code).

A **replay phantom** is a spatial arrangement of solved substrates -- an assembly of replay
packs ([`SPEC.md`](SPEC.md)) -- together with the per-voxel information needed to compose them
into a signal. It stores **no walkers and no trajectories of its own**. Everything physical
lives in the packs; the phantom adds only *which substrates are where, how they are oriented,
in what proportion, and with what proton density*.

This is a **secondary** format. It extends nothing in the replay invariant (SPEC §3) and adds
no capability tier (SPEC §7). A phantom is exactly as replayable as the packs it carries, and
its tier is the **intersection** of theirs.

## 1. Why a phantom is its own format

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
| `peak_dir` | `(N_v, P, 3)` | float32 | *peaks mode*: unit orientation of each slot |
| `odf_sh` | `(N_v, P, n_c)` | float32 | *ODF mode*: orientation distribution, even-order real SH |
| `scalars` | `(N_v, S)` | float32 | OPTIONAL per-voxel scalars, named in metadata |

`P` is the number of slots per voxel. A voxel using fewer pads with `substrate_id = -1` and
zero fraction.

**`geometric_fraction` is geometry, not signal.** It is the share of the voxel volume the
substrate occupies, before any relaxation or proton density is applied. Rows MUST sum to at
most one, and a replayer MUST NOT normalise a row that sums to less; a row summing to more is
an error in the phantom and MUST be rejected rather than rescaled.

That rule is what makes **partial volume** expressible. A voxel that is 60% white matter and
30% grey matter returns 90% of the signal the same voxel returns when those two fill it, because
the remaining 10% is not there. Normalising by the row sum would assert that the cited
substrates fill every voxel, turning every tissue boundary into pure tissue -- and partial
volume at boundaries is one of the main things a phantom is built to exercise. A short row is
therefore a statement, not a defect: it says what is modelled, and leaves what the remainder is
to the consumer.

## 4. Orientation: peaks or an ODF

Exactly one orientation mode is declared per phantom.

**Peaks** (`peak_dir`) give each slot a single direction. This is the representation for
discrete crossings: a voxel with two fibre populations is two slots with two directions, and
`P ≤ 3` covers the crossing configurations that are resolvable in practice. The two populations
may cite the *same* substrate at different orientations -- a pure crossing of one solved
microstructure -- or different substrates.

**ODF** (`odf_sh`) gives each slot an orientation distribution in the even-order real spherical
harmonic basis, for dispersion, fanning, and anything a discrete peak set cannot express.

The two are not separate physics. A peak is the zero-dispersion limit of an ODF, and a
conformant replayer MUST produce the same signal from a peak set as from ODFs concentrated on
those directions with the same weights, to within the SH truncation. Peaks are stored
separately because evaluating the response at a direction is exact and cheaper than contracting
a near-singular ODF, not because they mean something different.

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

## 6. Replay operation

For voxel `v` and acquisition `q`, with `E_i` the response of substrate `i` and `F_{v,p}` the
orientation of slot `p`,

```
S_v(q) = sum_p  geometric_fraction[v,p] * m0[substrate_id[v,p]]
                * INT_{S^2} F_{v,p}(n) E_{substrate_id[v,p]}(n; q) dn
```

with the integral replaced by `E(peak_dir[v,p]; q)` in peaks mode. Where the response depends
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
  "grid": {"shape": [180, 216, 180], "voxel_size_m": [1e-3, 1e-3, 1e-3], "frame": "RAS"},
  "orientation": {"mode": "peaks", "max_peaks": 3},   // or {"mode": "odf_sh", "lmax": 8,
                                                      //      "basis": "real",
                                                      //      "convention": "orthonormal"}
  "substrates": [
    {"id": "canonical/wm/g070-f055", "m0": 0.70, "embedded": true,
     "sha256": "…", "pack_meta": {…}},                // arrays under substrate0/
    {"id": "canonical/gm/…",  "m0": 0.85, "embedded": true,  "sha256": "…", "pack_meta": {…}},
    {"id": "canonical/csf/…", "m0": 1.00, "embedded": false, "sha256": "…", "uri": "hf://…"}
  ],
  "scalars": ["kappa_B1"],
  "license": "…", "citation": "…", "provenance": {…}
}
```

In ODF mode `orientation.convention` MUST be `orthonormal`. The composition goes through the
spherical-harmonic addition theorem, which is false for the non-orthonormal real conventions in
common use, and the error it introduces vanishes exactly when the gradient is parallel to `B0`
-- the one geometry a cursory check would test.

## 8. Conformance

A file is a conformant `.rph` when it is a safetensors container carrying the arrays of §3 and
the metadata of §7; every `substrate_id` resolves; `geometric_fraction` rows sum to at most
one; every substrate declares `m0` and `sha256`; embedded substrates carry their arrays under
`substrate{i}/` and their `pack_meta`; and, in ODF mode, `odf_sh` is in the declared orthonormal
basis.

A conformant replayer resolves or reads every cited substrate, refuses (never guesses) one it
cannot, composes per §6 over both axes where susceptibility is present, and reports the
phantom's tier as the intersection of its substrates' tiers.
