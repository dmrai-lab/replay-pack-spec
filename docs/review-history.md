# Review history — Replay Pack Specification

Preserved record of the pre-1.0 review of this specification, originally conducted on
`dmrai-lab/replay-pack-spec` PR #1 ("Replay Pack Specification v0.1.0-draft"). Comments
are reproduced here so the review provenance survives; the DOI'd history lives in git.

## Summary of the review
First version of the **Replay Pack (`.rpk`)** open specification, opened against an empty `main` so we can review and iterate **inline** before merging the first version.

## What to review
The whole diff is the first version. Suggested reading order:
1. **`SPEC.md`** — the normative document. Key sections to scrutinize:
   - §3 the **replay invariant** (the founding axiom) and the state-vs-knobs boundary.
   - §5 **channels** + §6 **replay operations** (the exact math each channel feeds).
   - §7 **capability tiers** (T0 Gradient → T4 Exchange), independent + partial-producer friendly.
   - §9 **codecs** — the "format is independent of compression" claim + §9.4 reserved storage keys.
   - §10–§11 metadata schema + declared envelope / self-certification.
   - §12 container (safetensors; canonical `"rpk"` header key, `"json"` legacy alias).
2. **`schema/rpk_metadata.schema.json`** — machine-readable metadata contract.
3. **`examples/`** — standalone producer + reference replayer + codec-aware conformance checker.

## Validated
- Standalone produce → validate → replay matches analytic `exp(-bD)` to the MC noise floor.
- The conformance checker validates real dmipy-sim packs (DiSCo, canonical-WM), incl. the RLE
  compartment + susceptibility-field codec keys and the legacy `json` header alias.

## Open questions for inline discussion
- Tier granularity: is T2 Surface / T4 Exchange the right split, or should surface-relaxivity fold into Relaxation?
- Should `permeability` ever be expressible (currently a fixed pack property, not a knob)?
- Field-tier map parameterization (§6.4): is the ℓ=2 + m=0 decomposition the right normative minimum?
- Metadata: which `provenance` fields should be REQUIRED vs RECOMMENDED for bank submission?
- Naming: `.rpk` / "Replay Pack" / tier names — lock before Zenodo deposit.

Licenses: spec **CC-BY-4.0**, code **Apache-2.0**. Intended trajectory: merge → Zenodo deposit (format/infrastructure priority) → replay paper (arXiv) → open the substrate bank.

## Inline review comments (chronological)

### `schema/rpk_metadata.schema.json:40`  — @rutgerfick  ·  2026-07-27

rotation of what?

### `schema/rpk_metadata.schema.json:64`  — @rutgerfick  ·  2026-07-27

i think these properties must have more explicit names, or at least a explicit description lik above. "rho" can mean a billion things but here you mean surface relaxivity right.

### `schema/rpk_metadata.schema.json:68`  — @rutgerfick  ·  2026-07-27

right now we don't have permeability as a replay knob no? if we figure out a way to do it later then we can add it at a next version of the spec

### `SPEC.md:161`  — @rutgerfick  ·  2026-07-27

does it have to be an *anisotropic* scale? can't it be isotropic? or can a generalized suscpetibility tensor be provided?

### `SPEC.md:168`  — @rutgerfick  ·  2026-07-27

why focus on spin echo refocussing sign only? we can do a GRE (no 180) or PGSTE (2x90) as well. epx_k doesn't make sense for pgste for example only being -1 and 1 no?

### `SPEC.md:208`  — @rutgerfick  ·  2026-07-27

again this doesn't make sense. any gradien tor RF pulse can be fired so this sign flip is PGSE only special casing that is irrelevant for the replay

### `SPEC.md:320`  — @rutgerfick  ·  2026-07-27

this codec and compression section is very nicely fleshed out but it also imagines a lot of normal "keys" for compression parameters that completely don't exist yet right? like we are coming up with these compression algorithms in the paper, and not even all these variables mean anything without the contex tof the compression algorithm iself. 

### `schema/rpk_metadata.schema.json:40`  — @rutgerfick  ·  2026-07-27

Clarified: `R` is an optional per-compartment 3×3 rotation mapping the compartment's local frame → lab frame, used to orient the anisotropic Field maps (SPEC §6.4); `null` if unoriented. Schema description + SPEC §5.2/§10 updated.

### `schema/rpk_metadata.schema.json:64`  — @rutgerfick  ·  2026-07-27

Agreed — renamed all envelope flags to explicit, self-describing names with descriptions: `rho`→`surface_relaxivity`, `T1T2`→`relaxation`, `B0_any`→`field_offresonance`, `orientation_any`→`field_orientation`, `mt`→`magnetization_transfer`. Coordinated rename in `bank.py` too; the old names are accepted as 1.x read-aliases so existing packs still validate.

### `schema/rpk_metadata.schema.json:68`  — @rutgerfick  ·  2026-07-27

Agreed — removed `permeability` from `replay_envelope` entirely (it's not a replay knob: crossing changes the trajectory). The fixed value the walk used now goes under `provenance`. If we ever make crossing replayable it's a new tier in a later version (SPEC §7 note + §14). Also dropped the now-meaningless `rf` flag.

### `SPEC.md:161`  — @rutgerfick  ·  2026-07-27

Good catch — no longer hardcoded to anisotropic. §6.4 now frames the Field channels as normalized off-resonance *basis maps*: `Φ0` = isotropic (m=0), `ΦC/ΦS` = ℓ=2 anisotropic, so the tier covers isotropic and/or axially-anisotropic susceptibility. A fully general rank-2 susceptibility *tensor* field is flagged as a §14 extension (needs a richer basis). `delta_chi_a` is now described as just the anisotropic scale the maps are normalized to.

### `SPEC.md:168`  — @rutgerfick  ·  2026-07-27

You're right, this was PGSE-shaped and `ε_k∈{+1,−1}` can't express PGSTE. Reworked: new §6.6 makes the *acquisition* (gradients **and** RF) a replay knob. Scalar tiers take a general per-step transverse-phase gate `s(t)`: `+1` for GRE, sign-flip at each 180° for SE/CPMG, and **`0` during z-storage for STE/PGSTE** (so `{−1,0,+1}`). Arbitrary RF (non-180° flips, adiabatic, MT sat) → the vector-Bloch replay (§6.5). No sequence is stored in the pack.

### `SPEC.md:208`  — @rutgerfick  ·  2026-07-27

Agreed — removed. §8.4 no longer mentions refocusing time; it now just says a producer MUST NOT bake any sequence choice into the pack. Refocusing is a property of the replayed sequence (§6.6), not of the walk.

### `SPEC.md:320`  — @rutgerfick  ·  2026-07-27

Agreed — this was coupling the frozen format to WIP compression research. §9 now defines only the codec *interface* (declared name+params, decodes-to-§5-contract, `identity` baseline, lossy self-certifies) plus the one structural rule (distributional codecs are Gradient-only). The concrete codecs and their reserved storage keys moved to a new, independently versioned **CODEC_REGISTRY.md** — algorithms land there as they solidify in the paper, without touching the core spec or bumping `rpk_schema_version`.

### `CODEC_REGISTRY.md:11`  — @rutgerfick  ·  2026-07-28

no need to say new method are being develoed. just stick to the compression methods being independent from the spec.

### `CODEC_REGISTRY.md:39`  — @rutgerfick  ·  2026-07-28

again no need to say any methodology paper in preparation anywhere in this PR. people can still come up with better ways to compress after it's prepared too. it just means that this is a derivative spec document that is more geared towards either memory or inference optimization, but not specification of the *meaning* of the input or output.

### `README.md:39`  — @rutgerfick  ·  2026-07-28

I guess I'm not such a huge fan of the "Exchange" part as an alternative naming for magnetization transfer. Exchange in my ears can also refer to permeability, which is something we don't replay. can we just call it magnetization transfer and not exchange?

### `CODEC_REGISTRY.md:4`  — @rutgerfick  ·  2026-07-28

why is this version 0.2 when the other spec is 0.1?

### `schema/rpk_metadata.schema.json:65`  — @rutgerfick  ·  2026-07-28

how can we have two T3 envelopes, do we treat them sepearately?

### `schema/rpk_metadata.schema.json:79`  — @rutgerfick  ·  2026-07-28

how is a user supposed to inperpret this fidelity field for lossy codecs? didn't we split off codecs?

### `schema/rpk_metadata.schema.json:89`  — @rutgerfick  ·  2026-07-28

here in the provenance we also need to actually describe the substrate (mesh or analytical substrate) itself, not just whether it was real or not. in fact whether it's real is irrelevant we just need to know what was the actual thing.

actuallhy ohw does this transfer to the substrate card on the substrate bank?

### `CODEC_REGISTRY.md:11`  — @rutgerfick  ·  2026-07-28

Done (86f9259). Dropped the "being developed / methodology paper" phrasing. Your comment-2 framing is now the registry's opening line: it's a **storage/inference-optimization** layer that governs only *how* a channel is stored, never its *meaning* — which is fixed by the core spec (§5–6) and independent of any codec.


### `CODEC_REGISTRY.md:39`  — @rutgerfick  ·  2026-07-28

Done (86f9259). Removed the "experimental / paper in preparation" line here, the "(some still being developed)" in SPEC §9, and the forthcoming-paper bullet in §16. The registry now reads as a derivative optimization layer for memory/inference — not a specification of the meaning of the input or output.


### `README.md:39`  — @rutgerfick  ·  2026-07-28

Agreed (86f9259). Renamed the T4 tier "Exchange" → "Magnetization transfer" everywhere (README, SPEC §7 / §6.5 / §8.8 / §2). "Exchange" wrongly implied permeability / diffusive exchange, which we do **not** replay. The flag stays `magnetization_transfer`.


### `CODEC_REGISTRY.md:4`  — @rutgerfick  ·  2026-07-28

Fair — reset to **0.1.0** (86f9259) so registry and spec ship at parity at launch. Independent versioning is by design (that's why the registry is a separate doc), but the divergence should start the first time a codec actually changes post-release, not on day one.


### `schema/rpk_metadata.schema.json:65`  — @rutgerfick  ·  2026-07-28

Collapsed to a single `field` flag (86f9259) — one tier, one flag, like the others. Whether orientation is a knob is now self-described by which basis maps are present: `susc_field_0` alone = isotropic (no θ dependence); the ℓ=2 `susc_field_C/S` add orientation. Old `field_offresonance`/`field_orientation` are kept as 1.x read-aliases, and the conformance checker now requires only `susc_field_0` for the tier so isotropic-only packs validate.


### `schema/rpk_metadata.schema.json:79`  — @rutgerfick  ·  2026-07-28

Good tension to flag. We split off the codec *algorithms* (→ registry), but `fidelity` is the codec-**agnostic** trust certificate — the decoded-vs-raw error whatever method produced it — so it's part of the meaning/trust contract and stays in the **core** (86f9259). Clarified the schema + SPEC §9 on *why* it stays and *how* to read it: `err_max` = worst replay error over `battery`, `floor_max` = the irreducible MC noise floor, `within_2x_floor` = loss below 2× that floor (negligible).


### `schema/rpk_metadata.schema.json:89`  — @rutgerfick  ·  2026-07-28

Agreed (86f9259). Replaced `real_or_synthetic`/`geometry` with a descriptive `substrate` block — `{representation: "analytical"|"mesh"|"voxel", description, ref}` — i.e. *what was actually walked*. Kept descriptive, not a normative geometry schema, since the geometry **format** is explicitly out of scope (§1).

On the bank: the substrate-bank catalog card is a downstream service, fed by `provenance` + the optional Croissant sidecar (§12); the card's layout is the bank's concern, not this format. Noted that in §10.



## Discussion (issue-level comments)

**@rutgerfick** · 2026-07-27

Pushed a revision round (`b41f1bf`) addressing all seven comments — replied inline on each thread. Summary of what changed:

- **RF/sequence is now a replay knob, not baked physics.** New §6.6: a general per-step transverse-phase gate `s(t)` for the scalar tiers (GRE/SE/CPMG/STE, incl. `s=0` during PGSTE storage); arbitrary RF → vector-Bloch (§6.5). Removed the spin-echo/TE-2 special-casing from §6.4 and the "refocusing time" generator invariant (§8.4).
- **Compression decoupled.** §9 is now just the codec *interface*; concrete codecs + reserved storage keys live in the new, independently versioned **`CODEC_REGISTRY.md`**.
- **Permeability out** of the envelope (→ provenance; future tier if ever replayable). `rf` flag dropped.
- **Explicit envelope names** (`relaxation`, `surface_relaxivity`, `field_offresonance`, `field_orientation`, `magnetization_transfer`) + schema descriptions; old names kept as 1.x read-aliases. Coordinated rename landed in `bank.py`.
- **Susceptibility generalized** to isotropic (`Φ0`) and/or anisotropic (ℓ=2), with the general tensor as a §14 extension; clarified `per_comp.R`.

Re-validated: standalone demo matches `exp(-bD)` to the MC floor; the conformance checker still passes the real DiSCo and canonical-WM packs via the alias + codec-key paths. Ready for another pass.

---

