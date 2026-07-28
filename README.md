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
- **Codec registry:** [`CODEC_REGISTRY.md`](CODEC_REGISTRY.md) — the concrete position *representations* and storage codecs and their stored keys, versioned separately from the frozen core so methods can be added without touching the format.
- **Metadata schema:** [`schema/rpk_metadata.schema.json`](schema/rpk_metadata.schema.json).
- **Examples:** [`examples/`](examples/) — a standalone writer + replayer (no dependencies
  beyond NumPy + safetensors), a conformance checker, and dmipy-sim interop.

## The one idea

> For a fixed geometry, diffusivity, and seed, the walker trajectories are **independent of**
> the gradient waveform, field, relaxation, and susceptibility. Those are **replay knobs**:
> compute the walk once, store it, and evaluate any acquisition as a cheap functional of the
> stored positions.

Everything follows from this invariant. Two things the spec states up front:

- **Resolution limit.** One walk replays any `G(t)` resolvable on its save step `Δt` (Nyquist) — temporal structure finer than `Δt` isn't stored.
- **Representation ≠ meaning.** A pack's meaning is the raw trajectory. Storing positions in a linear basis is a *change of representation* — closer to a compiler's intermediate representation (IR) than to a zip file — that turns replay into a projection onto `K` coefficients instead of an integral over every timestep. A full basis is exact wherever raw is; a truncated one is exact within its span and certifies the rest. So the format is representation-agnostic: a raw pack is fully conformant, and compression never changes what a pack *means* (§9 of the spec).

## Capability tiers (a "limited" generator is still interoperable)

| Tier | Adds channel | Replays |
|---|---|---|
| **C0 Gradient** | `positions` (required) | any `G(t)` resolvable at `Δt`, b-tensor, EAP |
| **C1 Bulk relaxation** | `compartment` | any intrinsic `T2`/`T1` |
| **C2 Surface** | `boundary_local_time` | any surface relaxivity `ρ` |
| **C3 Field** | `susc_field_{C,S,0}` | any `B0`, orientation |
| **C4 Magnetization transfer** | `bound_fraction` | magnetization transfer |

A producer declares the tiers it populates; a replayer refuses (never fakes) tiers a pack
does not carry. (Tiers are named **C0–C4** — not `T#`, which would collide with the
relaxation times `T1`/`T2`.)

## Status

**v0.1.1 (draft for comment).** The `1.x` container schema (channel names, metadata keys) is
frozen; semantics may be clarified. The reference implementation (dmipy-sim, private for now)
emits `rpk_schema_version = "1.1"`.

Intended trajectory: this spec is deposited (Zenodo) to **fix the format's definition and
date**; reference tooling and example substrates open thereafter.

## Feedback

Found an ambiguity, or want to propose a channel, tier, or codec? Open a **GitHub Issue** —
that is where specification discussion happens.

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
for a complete, dependency-light C0/C1 replayer, and [`examples/conformance_check.py`](examples/conformance_check.py)
to validate a pack against the spec.
