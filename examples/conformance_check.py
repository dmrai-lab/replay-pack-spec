"""Conformance checker for a Replay Pack against SPEC.md §13.

Validates: safetensors readability, metadata against schema/rpk_metadata.schema.json,
required `positions` channel + shape, tier-flag ⇒ channel-present consistency, shared
(N_w, N_t) across per-walker channels, and the lossy-codec ⇒ fidelity rule. Prints the
declared tiers and envelope.

Usage:  python conformance_check.py <pack.rpk>
Exit code 0 = conformant, 1 = not. Apache-2.0.
"""
import json
import os
import sys
import numpy as np
from safetensors import safe_open

HERE = os.path.dirname(os.path.abspath(__file__))
SCHEMA = os.path.join(HERE, "..", "schema", "rpk_metadata.schema.json")

# Channel presence is codec-aware (SPEC §9): a channel is "present" if its raw key OR any of
# its codec sub-array key-groups are in the file. These key-groups mirror CODEC_REGISTRY.md
# (the living registry) — update them there and here together when a codec is added.
CHANNEL_KEYS = {
    # `bridge_dst` writes pos_{x,y,z}; the four retired position codecs (dct_coeffs, modes+*)
    # are NOT listed -- a pack carrying them is refused, not read (registry 0.4.0).
    "positions": [["positions"], ["pos_x", "pos_y", "pos_z"]],
    # C1 is one or more occupancy columns: `comp` is required, others (e.g. `bound`) optional.
    "compartment": [["compartment"], ["comp_rle_vals", "comp_rle_lens", "comp_rle_counts"]],
    "boundary_local_time": [["boundary_local_time"], ["dlog_boundary_unit"],
                            ["blt_bridge_dst", "blt_start", "blt_endpoint"],
                            ["blt_counts", "blt_cols", "blt_qvals"], ["blt_dense_q"]],
    # MT's bound pool is a C1 column, not a channel: bfrac_rle_* is retired (registry 0.5.0).
    "bound_fraction": [["bound_fraction"], ["bound_frac"],
                       ["bound_rle_vals", "bound_rle_lens", "bound_rle_counts"]],
    "spin_weights": [["spin_weights"]],
    "susc_field_C": [["susc_field_C"]], "susc_field_S": [["susc_field_S"]],
    "susc_field_0": [["susc_field_0"]],
}
TIER_CHANNELS = {          # explicit replay_envelope flag -> channel(s) that MUST be present
    "bulk_relaxation": ["compartment"],
    "surface_relaxivity": ["boundary_local_time"],
    "magnetization_transfer": ["bound_fraction"],
    "field": ["susc_field_0"],   # isotropic Phi_0 is the minimum; l=2 C/S add orientation
}
# pre-rename 0.x aliases accepted on read (SPEC §10)
ENV_ALIASES = {"T1T2": "bulk_relaxation", "relaxation": "bulk_relaxation",
               "rho": "surface_relaxivity", "B0_any": "field", "orientation_any": "field",
               "field_offresonance": "field", "field_orientation": "field",
               "mt": "magnetization_transfer"}
DISTRIBUTIONAL_KEYS = {"coeff_mean", "coeff_cov", "coeff_quantiles"}
# `bridge_dst` is lossless only at K = N_t - 2 (the interior dimension), so it is not
# unconditionally lossless and does not belong here; `identity` is the only codec that is.
LOSSLESS_CODECS = {"identity"}
# Refused on sight rather than read: each stores different quantities under names a current
# reader would reinterpret, yielding plausible wrong numbers (SPEC §9.4, registry 0.4.0/0.5.0).
RETIRED_CODECS = {"temporal_dct", "lowrank", "gaussian", "marginal"}
RETIRED_KEYS = {"dct_coeffs": "positions as cosine bands (retired: use bridge_dst / pos_*)",
                "blt_dct_coeffs": "C2 as detrended cosine bands (retired: use blt_bridge_dst)",
                "bfrac_rle_vals": "MT as its own channel (retired: use the C1 `bound` column)"}


def _present(channel, keys):
    """True iff `channel` is realized by any codec key-group present in `keys`."""
    return any(all(k in keys for k in group) for group in CHANNEL_KEYS.get(channel, [[channel]]))


def _norm_env(meta):
    """replay_envelope with 0.x aliases folded onto the explicit names (SPEC §10)."""
    env = dict(meta.get("replay_envelope", {}))
    for old, new in ENV_ALIASES.items():
        if old in env and new not in env:
            env[new] = env[old]
    return env


def check(path):
    errs, warns = [], []
    try:
        with safe_open(path, framework="numpy") as f:
            keys = set(f.keys())
            hdr = f.metadata() or {}
            shapes = {k: f.get_slice(k).get_shape() for k in keys}
    except Exception as e:
        return [f"not a readable safetensors file: {e}"], []

    blob = hdr.get("rpk") or hdr.get("json")   # "rpk" canonical, "json" = 0.x legacy alias
    if not blob:
        return ["__metadata__ has no 'rpk' (or legacy 'json') JSON blob (SPEC §12)"], []
    meta = json.loads(blob)
    if "rpk" not in hdr and "json" in hdr:
        warns.append("metadata under legacy header key 'json'; SHOULD migrate to 'rpk' (SPEC §12)")

    # 1. metadata schema (best-effort; jsonschema optional)
    try:
        import jsonschema
        jsonschema.validate(meta, json.load(open(SCHEMA)))
    except ImportError:
        warns.append("jsonschema not installed; skipped full metadata schema validation")
    except Exception as e:
        errs.append(f"metadata fails schema: {getattr(e, 'message', e)}")

    # 2. retired representations are REFUSED, not read. Each stores different quantities under
    # names a current reader would reinterpret, so a stale pack decodes to plausible wrong numbers
    # rather than failing. The pack must be re-encoded from its master.
    method = meta.get("compression", {}).get("method", "identity")
    if method in RETIRED_CODECS:
        errs.append(f"position codec {method!r} is retired and MUST be refused; re-encode as "
                    f"'bridge_dst' (SPEC §9.4, registry 0.4.0)")
    for k, why in RETIRED_KEYS.items():
        if k in keys:
            errs.append(f"retired key {k!r} present — {why}; re-encode from the master")
    c1 = ((meta.get("compression", {}).get("channels") or {}).get("compartment") or {})
    if c1 and "columns" not in c1:
        errs.append("C1 metadata predates the occupancy-column layout (no 'columns'); the arrays "
                    "are readable but declared differently — re-encode (registry 0.5.0)")

    # 2b. required positions channel (codec-aware; shape only checkable for identity)
    if not _present("positions", keys):
        errs.append("missing required 'positions' channel / codec sub-arrays (SPEC §5.1/§9)")
    elif method == "identity":
        if len(shapes["positions"]) != 3 or shapes["positions"][2] != 3:
            errs.append(f"positions shape {shapes['positions']} is not (N_w, N_t, 3)")
    else:
        warns.append(f"codec {method!r}: positions shape checked only after decode")

    # 3. tier flag ⇒ channels present (codec-aware; explicit names + 0.x aliases)
    env = _norm_env(meta)
    for flag, chans in TIER_CHANNELS.items():
        if env.get(flag):
            missing = [c for c in chans if not _present(c, keys)]
            if missing:
                errs.append(f"replay_envelope.{flag}=true but missing channel(s) {missing} (SPEC §7/§13)")

    # 3b. distributional positions codec MUST NOT carry per-walker channels (SPEC §9)
    if keys & DISTRIBUTIONAL_KEYS:
        bad = [c for c in ("compartment", "boundary_local_time", "bound_fraction", "spin_weights")
               if _present(c, keys)]
        if bad:
            errs.append(f"distributional codec carries per-walker channel(s) {bad}; "
                        f"forbidden — such packs are Gradient-only (SPEC §9)")
        if any(env.get(f) for f in ("bulk_relaxation", "surface_relaxivity", "magnetization_transfer")):
            errs.append("distributional codec must declare Gradient tier only (SPEC §9)")
    if env.get("bulk_relaxation") and not (meta.get("per_comp") or {}).get("T2"):
        errs.append("bulk_relaxation tier declared but per_comp.T2 absent (SPEC §10)")
    if env.get("field"):
        wp = meta.get("walk_params", {})
        if wp.get("cell_size") is None:
            errs.append("field tier declared but walk_params.cell_size absent (SPEC §7)")

    # 4. shared (N_w, N_t) across per-walker/per-save channels
    if "positions" in shapes and len(shapes["positions"]) == 3:
        nw, nt = shapes["positions"][0], shapes["positions"][1]
        for c in ("compartment", "boundary_local_time", "bound_fraction"):
            if c in shapes and tuple(shapes[c][:2]) != (nw, nt):
                errs.append(f"channel {c} shape {shapes[c]} != (N_w={nw}, N_t={nt}, ...)")
        if "spin_weights" in shapes and shapes["spin_weights"][0] != nw:
            errs.append(f"spin_weights length {shapes['spin_weights'][0]} != N_w={nw}")

    # 5. lossy codec ⇒ fidelity present
    if method not in LOSSLESS_CODECS and not meta.get("fidelity"):
        errs.append(f"lossy codec {method!r} requires a 'fidelity' self-report (SPEC §9.3)")

    # 6. license + citation
    for k in ("license", "citation"):
        if not meta.get(k):
            errs.append(f"missing required metadata '{k}' (SPEC §13)")

    return errs, warns


def main(path):
    errs, warns = check(path)
    with safe_open(path, framework="numpy") as f:
        _h = f.metadata() or {}
        meta = json.loads(_h.get("rpk") or _h.get("json") or "{}")
    env = _norm_env(meta)
    tiers = [name for name, on in [
        ("C0-Gradient", env.get("gradient")), ("C1-BulkRelax", env.get("bulk_relaxation")),
        ("C2-Surface", env.get("surface_relaxivity")), ("C3-Field", env.get("field")),
        ("C4-MT", env.get("magnetization_transfer"))] if on]
    print(f"pack id:    {meta.get('id')}")
    print(f"schema:     {meta.get('rpk_schema_version')}  codec: {meta.get('compression', {}).get('method')}")
    print(f"tiers:      {', '.join(tiers) or '(none)'}")
    print(f"license:    {meta.get('license')}")
    for w in warns:
        print(f"  warn: {w}")
    if errs:
        print(f"NOT CONFORMANT ({len(errs)}):")
        for e in errs:
            print(f"  - {e}")
        return 1
    print("CONFORMANT ✓")
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__); sys.exit(2)
    sys.exit(main(sys.argv[1]))
