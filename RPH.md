# The Replay Phantom Specification (`.rph`)

**Status:** draft, versioned independently of the core `.rpk` specification.
**Version:** 0.1.0 (draft).
**License:** CC-BY-4.0 (text) / Apache-2.0 (reference code).

A **replay phantom** is a spatial arrangement of solved substrates -- an assembly of replay
packs ([`SPEC.md`](SPEC.md)) -- together with the per-voxel information needed to compose them
into a signal. It stores **no walkers and no trajectories**. Everything physical lives in the
packs it references; the phantom adds only *where they are, how they are oriented, and in what
proportion*.

This is a **secondary** format. It extends nothing in the replay invariant (SPEC §3) and adds
no capability tier (SPEC §7). A phantom is exactly as replayable as the packs it cites, and a
reader that can replay those packs plus the composition of §3 below can replay the phantom.

## 1. Why a phantom is its own format

A pack answers *what does this microstructure do to the magnetization*. It is one substrate at
one pose. A voxel is a distribution of poses, and a volume is a field of such distributions.
Those are properties of an arrangement, not of any substrate in it, and SPEC §14 already
separates them: transmit-field inhomogeneity becomes a per-voxel scalar map once packs are
tiled into a voxel grid, belonging to the arrangement rather than to any pack in it.

Keeping the phantom in its own file has three consequences worth the split. The expensive
object stays shared — one solved pack serves every voxel and every orientation that cites it,
so a phantom is small and a bank of phantoms costs little beyond the packs. The two evolve
independently: re-solving a substrate replaces a pack without touching any phantom, and
re-arranging tissue replaces a phantom without re-solving anything. And provenance stays
honest, because a phantom names the packs it is made of rather than absorbing and anonymising
them.

## 2. Data model

Arrays, one safetensors file (SPEC §12 conventions apply unchanged):

| Array | Shape | dtype | Meaning |
|---|---|---|---|
| `voxel_index` | `(N_v, D)` | int32 | voxel coordinates on the grid; `D` is 2 or 3 |
| `pack_id` | `(N_v, P)` | int16 | index into `packs` (§4); `-1` marks an empty slot |
| `pack_fraction` | `(N_v, P)` | float32 | volume fraction per cited pack; row sums ≤ 1 |
| `odf_sh` | `(N_v, P, n_c)` | float32 | orientation distribution per (voxel, pack), even-order real SH |
| `scalars` | `(N_v, S)` | float32 | OPTIONAL per-voxel phantom scalars, named in metadata |

`P` is the maximum number of packs cited by any voxel. A voxel citing fewer pads with
`pack_id = -1` and zero fraction. A row summing to less than one leaves the remainder
unmodelled — a phantom MUST NOT silently normalise it, since "the rest is free water" is a
modelling choice belonging to the consumer.

The single-orientation case is the `n_c = 1` degenerate ODF; storing a direction instead of a
distribution is not a separate mode. That keeps one composition path rather than two.

## 3. Replay operation

For voxel `v` and acquisition `q`, with `E_p` the pack response and `F_{v,p}` the orientation
distribution,

```
S_v(q) = sum_p  pack_fraction[v,p] * INT_{S^2}  F_{v,p}(n) E_p(n; q)  dn
```

The inner integral is the composition of the coupled angular spectrum: linear in the SH
coefficients of `F`, so it is a contraction rather than a quadrature per voxel. Where the
response depends on the field direction as well as the gradient — susceptibility — the two axes
must be composed jointly, and the reduction to independent one-dimensional convolutions does
not hold. A conformant replayer MUST compose both axes.

Composition is **linear in `F` and in `pack_fraction`**, so a phantom introduces no new physics
and cannot repair a pack that does not carry the tier an acquisition needs: a phantom's tier is
the *intersection* of the tiers of the packs it cites.

## 4. Metadata

JSON under the safetensors header key **`"rph"`** (SPEC §12 conventions):

```jsonc
{
  "rph_schema_version": "0.1.0",
  "id": "phantoms/circular-wm/annulus-30",
  "packs": [                          // resolved by id, NOT embedded
    {"id": "winther/g6/axon06", "sha256": "…", "uri": "hf://…"}
  ],
  "grid": {"shape": [30, 30], "voxel_size_m": [1e-3, 1e-3], "frame": "…"},
  "sh": {"basis": "real", "convention": "orthonormal", "lmax": 8},
  "scalars": ["kappa_B1"],            // names for the columns of `scalars`
  "license": "…", "citation": "…", "provenance": {…}
}
```

Two fields carry weight. `packs[].sha256` pins *which* solved physics a result used, so a
phantom is reproducible even as a bank re-issues packs. And `sh.convention` MUST be stated and
MUST be `orthonormal`: the composition goes through the spherical-harmonic addition theorem,
which is false for the non-orthonormal real conventions in common use, and the error it causes
vanishes exactly when the gradient is parallel to `B0` — the one geometry a cursory check
would test.

## 5. Conformance

A file is a conformant `.rph` when it is a safetensors container carrying the arrays of §2 and
the metadata of §4; every `pack_id` resolves to a cited pack; `pack_fraction` rows sum to at
most one; and `odf_sh` is in the declared orthonormal basis.

A conformant replayer resolves the cited packs, refuses (never guesses) an unresolvable one,
composes per §3 over both axes where susceptibility is present, and reports the phantom's tier
as the intersection of its packs' tiers.
