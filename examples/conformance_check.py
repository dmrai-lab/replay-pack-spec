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
    "positions": [["positions"], ["dct_coeffs"], ["modes", "coeffs"],
                  ["modes", "coeff_mean", "coeff_cov"], ["modes", "coeff_quantiles"]],
    "compartment": [["compartment"], ["comp_rle_vals", "comp_rle_lens", "comp_rle_counts"]],
    "boundary_local_time": [["boundary_local_time"], ["dlog_boundary_unit"],
                            ["blt_counts", "blt_cols", "blt_qvals"]],
    "bound_fraction": [["bound_fraction"], ["bound_frac"],
                       ["bfrac_rle_vals", "bfrac_rle_lens", "bfrac_rle_counts"]],
    "spin_weights": [["spin_weights"]],
    "susc_field_C": [["susc_field_C"]], "susc_field_S": [["susc_field_S"]],
    "susc_field_0": [["susc_field_0"]],
}
TIER_CHANNELS = {          # explicit replay_envelope flag -> channel(s) that MUST be present
    "bulk_relaxation": ["compartment"],
    "surface_relaxivity": ["boundary_local_time"],
    "magnetization_transfer": ["bound_fraction"],
    "field_offresonance": ["susc_field_C", "susc_field_S", "susc_field_0"],
}
# pre-rename 1.x aliases accepted on read (SPEC §10)
ENV_ALIASES = {"T1T2": "bulk_relaxation", "relaxation": "bulk_relaxation",
               "rho": "surface_relaxivity", "B0_any": "field_offresonance",
               "orientation_any": "field_orientation", "mt": "magnetization_transfer"}
DISTRIBUTIONAL_KEYS = {"coeff_mean", "coeff_cov", "coeff_quantiles"}
LOSSLESS_CODECS = {"identity", "temporal_dct", "lowrank"}


def _present(channel, keys):
    """True iff `channel` is realized by any codec key-group present in `keys`."""
    return any(all(k in keys for k in group) for group in CHANNEL_KEYS.get(channel, [[channel]]))


def _norm_env(meta):
    """replay_envelope with 1.x aliases folded onto the explicit names (SPEC §10)."""
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

    blob = hdr.get("rpk") or hdr.get("json")   # "rpk" canonical, "json" = 1.x legacy alias
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

    # 2. required positions channel (codec-aware; shape only checkable for identity)
    method = meta.get("compression", {}).get("method", "identity")
    if not _present("positions", keys):
        errs.append("missing required 'positions' channel / codec sub-arrays (SPEC §5.1/§9)")
    elif method == "identity":
        if len(shapes["positions"]) != 3 or shapes["positions"][2] != 3:
            errs.append(f"positions shape {shapes['positions']} is not (N_w, N_t, 3)")
    else:
        warns.append(f"codec {method!r}: positions shape checked only after decode")

    # 3. tier flag ⇒ channels present (codec-aware; explicit names + 1.x aliases)
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
    if env.get("field_offresonance"):
        wp = meta.get("walk_params", {})
        if wp.get("cell_size") is None:
            errs.append("field_offresonance tier declared but walk_params.cell_size absent (SPEC §7)")

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
        ("T0-Gradient", env.get("gradient")), ("T1-BulkRelax", env.get("bulk_relaxation")),
        ("T2-Surface", env.get("surface_relaxivity")), ("T3-Field", env.get("field_offresonance")),
        ("T4-Exchange", env.get("magnetization_transfer"))] if on]
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
