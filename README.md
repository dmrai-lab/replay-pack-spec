# Replay Pack (`.rpk`) — an open format for replayable Monte-Carlo diffusion-MRI walks

A **Replay Pack** is a portable, code-free file that stores the *state of a Monte-Carlo
random walk* through a diffusion-MRI substrate so that the measured signal can be
**replayed** for many different acquisitions and physical settings — any gradient
waveform, b-tensor, field strength, orientation, `T2`/`T1`, surface relaxivity, or
magnetization-transfer setting — **without re-running the expensive walk**.

This repository is the **normative specification** of that format, plus a small
standalone reference implementation and worked examples. Its goal is to let *any*
Monte-Carlo simulator emit interoperable packs (even partial ones), served by a single
reference replayer, and pooled in a public substrate bank.

- **The spec:** [`SPEC.md`](SPEC.md) — the definitive, RFC-2119 document.
- **Metadata schema:** [`schema/rpk_metadata.schema.json`](schema/rpk_metadata.schema.json).
- **Examples:** [`examples/`](examples/) — a standalone writer + replayer (no dependencies
  beyond NumPy + safetensors), a conformance checker, and dmipy-sim interop.

## The one idea

> For a fixed geometry, diffusivity, and seed, the walker trajectories are **independent of**
> the gradient waveform, field, relaxation, and susceptibility. Those are **replay knobs**:
> compute the walk once, store it, and evaluate any acquisition as a cheap functional of the
> stored positions.

Everything in the format follows from this invariant. In particular **the format does not
depend on compression**: a pack storing raw, uncompressed positions is fully conformant.
Compression is an *orthogonal* per-channel codec layer (§9 of the spec).

## Capability tiers (a "limited" generator is still interoperable)

| Tier | Adds channel | Replays |
|---|---|---|
| **T0 Gradient** | `positions` (required) | any `G(t)`, b-tensor, EAP |
| **T1 Relaxation** | `compartment` | any `T2`/`T1` |
| **T2 Surface** | `boundary_local_time` | any surface relaxivity `ρ` |
| **T3 Field** | `susc_field_{C,S,0}` | any `B0`, orientation |
| **T4 Exchange** | `bound_fraction` | magnetization transfer / exchange |

A producer declares the tiers it populates; a replayer refuses (never fakes) tiers a pack
does not carry.

## Status

**v0.1.0 (draft for comment).** The `1.x` container schema (channel names, metadata keys) is
frozen; semantics may be clarified. The reference implementation
([dmipy-sim](https://github.com/) — private for now) emits `rpk_schema_version = "1.1"`.

Intended trajectory: this spec is deposited (Zenodo) to establish priority on the
infrastructure/format; the methodology paper follows (arXiv); the substrate bank and example
substrates open thereafter.

## Licensing

Dual-licensed by artifact type:

- **Specification & documentation** (`SPEC.md`, this README, schema docs): **CC-BY-4.0** — see [`LICENSE-SPEC`](LICENSE-SPEC).
- **Code** (`examples/`, `schema/` machinery): **Apache-2.0** — see [`LICENSE`](LICENSE).

Adopt, implement, and extend freely under those terms. A pack's *own* `license` field records
the license of the **source substrate**, which this format never relicenses.

## Citing

See [`CITATION.cff`](CITATION.cff). Please cite the specification when you implement or emit
Replay Packs.

## Implementing a replayer

You need ~40 lines: read safetensors, validate metadata, and implement the gradient phase
`φ = γ Δt Σ G(t)·r(t)` (spec §6.1). See [`examples/reference_replayer.py`](examples/reference_replayer.py)
for a complete, dependency-light T0/T1 replayer, and [`examples/conformance_check.py`](examples/conformance_check.py)
to validate a pack against the spec.
