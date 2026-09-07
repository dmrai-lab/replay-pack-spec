# The Substrate Specification (`.sub.json`)

*Companion to the Replay Pack Specification (`SPEC.md`) and the Replay Phantom Specification (`RPH.md`).
Draft `0.1.1` for comment; nothing is numbered `1.0` before publication. `0.1.1`: `susceptibility.chi_iso` / `chi_aniso` MAY be `null` -- the producer declares the field source and not its values.*

## 1. Scope and purpose

A **substrate spec** is a portable, declarative, code-free description of *the situation a random
walker is in*: the domain and what its faces do, the water pools and their bulk properties, the walls
between pools and what each wall does to a spin that hits it from either side, which pools are
populated at `t = 0`, which pools carry a magnetic susceptibility, and what the substrate was asked to
be versus what it turned out to be. It is the **only** input a conformant Monte-Carlo generator
accepts, and every Replay Pack (`.rpk`, §10 of `SPEC.md`) embeds the spec of the substrate it walked.

It exists because a substrate that is a program (a class, a loader, a notebook) cannot be audited,
reproduced or exchanged, and because the same ambiguity re-appears at every stage: whether the space
outside a cylinder is a pool, whether a box face wraps or reflects, which pool the walkers start in,
which pool has a susceptibility, whether a packing achieved the fraction it was asked for. This
document makes each of those a field.

This specification defines the data model (§3), the meaning of each field for the walk (§4), the
request/realisation contract for generated substrates (§5), the relation to the pack (§6), producers
(§7) and conformance (§8). It does not standardise how a generator integrates the walk, nor any mesh
file format: surfaces are referenced as files in open formats (§3.4).

## 2. Terminology

- **Pool** — a region of space with one set of bulk properties, identified by a small integer id.
  Id `0` is always the free / extra-cellular pool. Ids are dense: `0 .. n_pools - 1`.
- **Wall** — a surface separating two pools (or a pool from *void*, §3.4). Walls carry the wall
  physics: permeability per crossing direction, surface relaxivity and MT reactivity per side.
- **Domain** — the box the substrate is defined in and the boundary condition on each of its six faces
  (three axes): `periodic` (the face wraps), `reflect` (the face is a specular wall), `open` (there is
  no face; walkers may leave the box and the substrate is undefined beyond it).
- **Void** — the complement of every pool. A walker never enters void: a wall whose outside is void
  is impermeable by definition.
- **Request / realisation** — for a generated substrate, what the generator was asked for and what it
  produced, recorded separately so a consumer can tell the two apart.

Units are SI throughout (metres, seconds, m²/s, m/s, tesla), as in `SPEC.md` §4.1. Positions and
directions are in the **substrate frame**; an acquisition rotation is applied by the engine and never
changes the spec.

## 3. Data model

A spec is one JSON object (`schema/substrate.schema.json` is normative for types and presence).

### 3.1 Header

| key | type | meaning |
|---|---|---|
| `substrate_spec_version` | string | `"0.1"` |
| `id` | string | stable identifier; a pack's `provenance.substrate.id` equals it |
| `description` | string | free text |

### 3.2 `domain`

| key | type | meaning |
|---|---|---|
| `box_min`, `box_max` | number[3] | the box (m); `box_min < box_max` on every axis |
| `boundary` | string[3] | per axis, one of `periodic`, `reflect`, `open` |

A packed cell is `["periodic", "periodic", "open"]`; an isolated analytic object is
`["open", "open", "open"]` with a box large enough to seed and to hold any field grid; a voxelised
mesh substrate declares whatever its box does. `reflect` faces are walls of the domain and take the
domain-wide `domain.surface_relaxivity` (default 0) on their inside.

### 3.3 `frame`

| key | type | meaning |
|---|---|---|
| `axis` | number[3] | unit vector: the substrate's principal axis in the substrate frame (fibre axis) |
| `in_plane` | number[3] or null | optional second axis fixing the in-plane orientation |

### 3.4 `pools`

An array; `pools[i].id == i` and `pools[0].name == "extra"` (or `"free"`). Each pool:

| key | type | meaning |
|---|---|---|
| `id` | integer | dense id |
| `name` | string | `extra`, `intra`, `myelin`, or another unique name |
| `D` | number or null | free diffusivity in the pool (m²/s); `0` for a stuck pool; `null` means the diffusivity the walk is driven with (an analytic geometry that carries none) |
| `T2`, `T1` | number or null | **nominal** relaxation times (s); replay knobs, never copied into a pack |
| `water_fraction` | number in [0, 1] | proton density relative to free water; the seeding weight |
| `susceptibility` | object or null | `{chi_iso, chi_aniso, director}`; `director` is `"none"` (isotropic), `"radial"` (a sheath: the local outward normal of the nearest wall) or `"file"` (a per-voxel director grid, `file`). The values are **nominal**: the field basis is derived from the pool's occupancy and director alone, and `chi_iso` / `chi_aniso` are applied at replay, and are `null` when the producer declares the source but not its values (an analytic sheath) |

A pool with a `susceptibility` object is a **field source**; the set of field-source pools answers
"which pools generate susceptibility fields". A substrate with no field-source pool has no Field tier.

### 3.5 `walls`

An array. Each wall:

| key | type | meaning |
|---|---|---|
| `name` | string | unique |
| `surface` | object | see below |
| `inside_pool` | integer | the pool on the surface's inside |
| `outside_pool` | integer or null | the pool outside, or `null` for void (then `permeability` MUST be zero) |
| `permeability` | `{in_to_out, out_to_in}` | κ (m/s) per crossing direction; `0` is impermeable. Unequal values are a pump, not passive exchange, and MUST be declared as such in `description` |
| `surface_relaxivity` | `{inside, outside}` | ρ₂ (m/s) seen from each side |
| `mt_reactivity` | `{inside, outside}` | κ_MT (m/s) per side, `0` off |

`surface.kind` is one of `sphere`, `cylinder`, `ellipsoid`, `plane`, `swept_polyline`, `mesh`, with:

- `sphere`: `center[3]`, `radius`; `cylinder`: `center[3]`, `axis[3]`, `radius` (infinite along `axis`
  unless `length` is given); `ellipsoid`: `center[3]`, `semiaxes[3]`, `rotation[3][3]` optional;
  `plane`: `point[3]`, `normal[3]` (inside = the half-space the normal points away from);
  `swept_polyline`: `centerline[n][3]`, `radius` (a sphere-swept polyline); `mesh`: `file`,
  `format` (`ply`, `obj`, `stl`), `scale` (file units → metres), `sha256`.
- `instances` (optional): arrays of per-instance parameters (`centers`, `radii`, …) so a packing of
  `N` identical-role objects is **one** wall entry with `N` instances sharing pools and properties;
  instance `k` is the `k`-th object where a consumer needs object ids (`SPEC.md` §8.5).

Two walls MAY share a pool on one side (a myelin sheath is `inner: intra|myelin` and
`outer: myelin|extra`). A pool MUST be bounded by walls or by `reflect`/`periodic` faces unless it is
pool `0` in an `open` domain.

### 3.6 `seeding`

| key | type | meaning |
|---|---|---|
| `pools` | integer[] | the pools populated at `t = 0` |
| `rule` | string | `uniform_by_volume` (density ∝ volume × `water_fraction`) or `explicit` (positions supplied by the driver) |
| `weights` | string | `water_fraction` (spins carry the pool's water fraction) or `thin` (spins are dropped to the water fraction so weights are uniform) |

### 3.7 `request` and `realisation` (generated substrates)

Both optional and free in form, but a producer that sets one SHOULD set both and SHOULD use these keys
when they apply: `generator {name, version}`, `seed`, `n_objects`, `packing_fraction`,
`diameter_law {family, shape, scale, d_min}`, `g_ratio`, `min_gap`, `cell_side`, `smallest_feature`.
`realisation` records what was **achieved**; a consumer MUST read the fraction, the cell and the
count from `realisation`, never from `request`.

### 3.8 `validity`

| key | type | meaning |
|---|---|---|
| `smallest_feature` | number | the smallest length the walk must resolve (smallest radius, thinnest sheath); the sub-step rules divide it |
| `min_gap` | number or null | narrowest passage between walls |
| `mesh_edge_feature_ratio` | number or null | median edge over `smallest_feature`, meshes only |
| `tiers` | string[] | subset of `gradient`, `relaxation`, `surface`, `field`, `exchange`: what a walk on this substrate CAN record |

### 3.9 `provenance`

`source` (dataset / generator), `files [{path, sha256}]`, `scale`, `transformations` (ordered list of
what was done to the source: dropped surfaces, box choice, unit scale, g-ratio measurement),
`created` (ISO date), `software {name, version}`.

## 4. Meaning for the walk (normative)

1. A walker's pool at `t = 0` is the pool it was seeded in (§3.6). Its pool changes only by a granted
   crossing at a wall with non-zero `permeability` in that direction. The engine MUST reject any step
   that would change the pool otherwise (`SPEC.md` §8.5).
2. A wall hit applies, from the side hit, `surface_relaxivity` as a boundary local time increment and
   `mt_reactivity` as a binding probability; a crossing decision uses `permeability` for the direction
   of travel. A wall with `outside_pool: null` reflects always.
3. A `periodic` face wraps the position (the stored trajectory stays continuous); a `reflect` face is
   a specular wall of the domain; an `open` face does nothing.
4. The pool's `D` sets the step length of a walker in it; unequal `D` across a permeable wall is a
   diffusivity-discontinuity interface and MUST be refused unless the engine implements it.
5. The Field tier's basis is computed from the occupancy of every field-source pool and its director,
   on the domain box (periodic where the domain is): geometry only. `chi_iso`, `chi_aniso`, `B0` and
   its direction are replay knobs.
6. `T2`, `T1`, `water_fraction`-as-weight, `surface_relaxivity`, `permeability`, `mt_reactivity` and
   `susceptibility` values are **nominal** in the spec. Of these only `permeability`, `mt_reactivity`
   (when walked emergently), `water_fraction` and `D` shape the walk; the rest are replay knobs and
   MUST NOT appear in a pack as values (§6).

## 5. Request and realisation

A generator MUST check feasibility **before** it walks or packs: a packing fraction above what its
placement method can reach, a `min_gap` below the step floor, a `diameter_law` whose `d_min` yields a
`smallest_feature` that the engine's sub-step cap would refuse — each is a refusal that names the
limit, not a failure after the attempts run out. `realisation` MUST be filled from the produced
geometry, not copied from `request`.

## 6. Relation to the Replay Pack

A pack embeds the spec verbatim under `substrate` in its metadata, and `provenance.substrate.id`
equals `spec.id`. The pack's compartment channel uses the spec's pool ids; its object ids, where
present, are wall instance indices. A pack MUST NOT carry `T2`, `T1`, `ρ`, `χ`, `B0` or MT rate
*values*: it carries channels (positions, occupancy, boundary local time, bound fraction, field basis)
and the spec that says what they mean. A replayer takes the physical values from the caller, who MAY
take them from the spec's nominal fields — explicitly, never as a silent default.

## 7. Producers

Anything that constructs a substrate is a **producer** and emits a spec (plus surface files): an
analytic constructor, a random packer, a mesh pipeline (CACTUS, CATERPillar, Winther), a strand
phantom (DiSCo). A producer MUST NOT construct engine objects; the engine reads the spec. The
producer records `provenance.transformations` for every decision it took on the user's behalf
(dropping an open surface, choosing a box, measuring a g-ratio).

## 8. Conformance

A spec is conformant when it validates against `schema/substrate.schema.json` and satisfies the
invariants the schema cannot express: dense pool ids with `0` the free pool; every wall's pools exist;
`outside_pool: null` implies zero permeability; every seeded pool exists; `box_min < box_max`;
`smallest_feature > 0`; `tiers` consistent with the content (`surface` needs a wall, `field` a
field-source pool, `relaxation` more than one pool or a `T2`, `exchange` a non-zero `mt_reactivity`
or `permeability`). The reference validator is `dmipy_sim.spec.validate`.
