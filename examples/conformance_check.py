"""Conformance checker for a Replay Pack against RPK.md §13.

Validates: safetensors readability, metadata against schema/rpk_metadata.schema.json, the
refusals of §9.5.1 (a pack carrying n_t, dt_traj, compression.basis, the identity codec or a
save-grid tensor is refused BY NAME with the instruction to re-walk it), the required band
tensors of `positions`, tier-flag => channel-present consistency, shared N_w across per-walker
channels, the event lists of the occupancy (§9.5.7), and the certificate (§9.5.5): block, qv,
fidelity.terms. Prints the declared tiers and envelope.

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

# Channel presence is the registry's stored key-groups (CODEC_REGISTRY.md): a channel is
# "present" if all keys of any one group are in the file. There is no raw key.
CHANNEL_KEYS = {
    "positions": [["pos_x", "pos_y", "pos_z"],
                  ["pos_x_ends", "pos_y_ends", "pos_z_ends", "pos_band_scale"]],
    # C1 is one or more occupancy columns: `comp` is required, others (e.g. `bound`) optional.
    "compartment": [["comp_static"], ["comp_events_offsets", "comp_events_t", "comp_events_label"]],
    "bound": [["bound_static"], ["bound_events_offsets", "bound_events_t", "bound_events_label"]],
    "boundary_local_time": [["blt_bridge_dst", "blt_start", "blt_endpoint"],
                            ["blt_band_scale", "blt_start", "blt_endpoint"]],
    "spin_weights": [["spin_weights"]],
    "field": [["susc_path_dct", "susc_path_scale"]],
}
TIER_CHANNELS = {          # explicit replay_envelope flag -> channel(s) that MUST be present
    "bulk_relaxation": ["compartment"],
    "surface_relaxivity": ["boundary_local_time"],
    "magnetization_transfer": ["bound"],      # or the parametric `mt` descriptor
    "field": ["field"],                       # the path modes (SPEC §6.4.1)
}
DISTRIBUTIONAL_KEYS = {"coeff_mean", "coeff_cov", "coeff_quantiles"}
RETIRED_CODECS = {"temporal_dct", "lowrank", "gaussian", "marginal"}
# Refused tensors (registry, "Refused forms"): each belongs to a form this revision does not read.
REFUSED_KEYS = {
    "susc_field_0": "a field grid as a field route", "susc_field_C": "a field grid as a field route",
    "susc_field_S": "a field grid as a field route",
    "positions": "raw positions on a save grid",
    "compartment": "the occupancy sampled on a save grid",
    "boundary_local_time": "the contact sampled on a save grid",
    "susc_field_basis": "the field sampled on a save grid",
    "dlog_boundary_unit": "the contact sampled on a save grid",
    "blt_counts": "the contact sampled on a save grid", "blt_cols": "the contact sampled on a save grid",
    "blt_qvals": "the contact sampled on a save grid", "blt_dense_q": "the contact sampled on a save grid",
    "comp_rle_vals": "the occupancy as runs on a save grid", "comp_rle_lens": "the occupancy as runs on a save grid",
    "comp_rle_counts": "the occupancy as runs on a save grid", "bound_rle_vals": "the occupancy as runs on a save grid",
    "bound_rle_lens": "the occupancy as runs on a save grid", "bound_rle_counts": "the occupancy as runs on a save grid",
    "dct_coeffs": "positions as cosine bands", "blt_dct_coeffs": "the contact as detrended cosine bands",
    "bfrac_rle_vals": "MT as a channel of its own (it is the `bound` column)",
}
INSTRUCTION = ("this revision does not read it: re-walk the pack with a conformant producer "
               "(SPEC §9.5.1). A pack is never converted.")


def _present(channel, keys):
    """True iff `channel` is realized by any registry key-group present in `keys`."""
    return any(all(k in keys for k in group) for group in CHANNEL_KEYS.get(channel, [[channel]]))


def _walk(obj, path=""):
    """Yield (path, key) for every dict key under obj."""
    if isinstance(obj, dict):
        for k, v in obj.items():
            yield f"{path}.{k}" if path else k, k
            yield from _walk(v, f"{path}.{k}" if path else k)
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            yield from _walk(v, f"{path}[{i}]")


def refusals(meta, keys):
    """The forms SPEC §9.5.1 refuses, each named, with the instruction to re-walk."""
    out = []
    wp, comp = meta.get("walk_params", {}), meta.get("compression", {})
    for k in ("n_t", "dt_traj"):
        if k in wp:
            out.append(f"walk_params.{k} is present: {INSTRUCTION}")
    for k in ("n_t",):
        if k in (wp.get("windows") or {}):
            out.append(f"walk_params.windows.{k} is present: {INSTRUCTION}")
    for path, k in _walk(comp, "compression"):
        if k == "n_t" or path == "compression.basis":
            out.append(f"{path} is present: {INSTRUCTION}")
    if "cell_size" in wp:
        out.append(f"walk_params.cell_size is present (a field grid as a field route): {INSTRUCTION}")
    if "per_comp" in meta:
        out.append("per_comp is present (a pack stores no tissue value, SPEC §8.5): re-walk the pack")
    method = comp.get("method")
    if method == "identity":
        out.append(f"compression.method 'identity' is present: {INSTRUCTION}")
    if method in RETIRED_CODECS:
        out.append(f"compression.method {method!r} is present, a retired codec: {INSTRUCTION}")
    for k, why in REFUSED_KEYS.items():
        if k in keys:
            out.append(f"tensor {k!r} is present ({why}): {INSTRUCTION}")
    for k in sorted(keys):
        if k.endswith(("_rle_vals", "_rle_lens", "_rle_counts")) and k not in REFUSED_KEYS:
            out.append(f"tensor {k!r} is present (the occupancy as runs on a save grid): {INSTRUCTION}")
    return out


def check(path):
    errs, warns = [], []
    try:
        with safe_open(path, framework="numpy") as f:
            keys = set(f.keys())
            hdr = f.metadata() or {}
            shapes = {k: f.get_slice(k).get_shape() for k in keys}
            events = {k: f.get_tensor(k) for k in keys if "_events_" in k}
    except Exception as e:
        return [f"not a readable safetensors file: {e}"], []

    blob = hdr.get("rpk")
    if not blob:
        return ["__metadata__ has no 'rpk' JSON blob (SPEC §12)"], []
    meta = json.loads(blob)

    # 1. the refusals of §9.5.1 come first: nothing of a refused pack is read.
    ref = refusals(meta, keys)
    if ref:
        return ref, warns

    # 2. metadata schema (best-effort; jsonschema optional)
    try:
        import jsonschema
        jsonschema.validate(meta, json.load(open(SCHEMA)))
    except ImportError:
        warns.append("jsonschema not installed; skipped full metadata schema validation")
    except Exception as e:
        errs.append(f"metadata fails schema: {getattr(e, 'message', e)}")

    comp = meta.get("compression", {})
    method = comp.get("method")
    wp = meta.get("walk_params", {})
    K = comp.get("K")

    # 3. the pack is bands, accumulated in blocks and certified (SPEC §9.5.5)
    if not wp.get("block"):
        errs.append("walk_params.block {h, m, precision} is missing (SPEC §9.5.5)")
    if not (wp.get("qv") or {}).get("positions"):
        errs.append("walk_params.qv.positions is missing (SPEC §9.5.5)")
    fid = meta.get("fidelity") or {}
    terms = fid.get("terms") or {}
    if "truncation" not in terms or "positions" not in (terms.get("quadrature") or {}):
        errs.append("fidelity.terms {truncation, quadrature.positions} is missing (SPEC §9.5.5)")
    if "local_time" in (wp.get("qv") or {}) and "local_time" not in (terms.get("quadrature") or {}):
        errs.append("walk_params.qv.local_time without fidelity.terms.quadrature.local_time (SPEC §9.5.5)")
    if "field" in (wp.get("qv") or {}) and "field" not in (terms.get("quadrature") or {}):
        errs.append("walk_params.qv.field without fidelity.terms.quadrature.field (SPEC §9.5.5)")
    if "field_read_interval" in wp and "field_read" not in terms:
        errs.append("walk_params.field_read_interval without fidelity.terms.field_read (SPEC §9.5.5)")
    if ("field" in (wp.get("qv") or {})) and "field_read_interval" not in wp:
        errs.append("walk_params.qv.field without walk_params.field_read_interval (SPEC §9.5.5)")
    if fid and "err_max" in fid and terms.get("truncation") is not None:
        s = terms["truncation"] + sum((terms.get("quadrature") or {}).values()) + terms.get("field_read", 0.0)
        if fid["err_max"] < s * (1 - 1e-9):
            errs.append(f"fidelity.err_max {fid['err_max']:.3e} is below the sum of its terms {s:.3e} (SPEC §9.5.5)")

    # 4. required positions channel: the band tensors of width K + 2
    if not _present("positions", keys):
        errs.append("missing required 'positions' channel: pos_x, pos_y, pos_z (SPEC §5.1/§9, registry)")
        n_w = None
    elif "pos_x" in keys:
        n_w = shapes["pos_x"][0]
        for a in "xyz":
            sh = shapes.get(f"pos_{a}")
            if sh is None or len(sh) != 2 or (K is not None and sh[1] != K + 2) or sh[0] != n_w:
                errs.append(f"pos_{a} shape {sh} is not (N_w, K+2) with K = {K}")
    else:
        n_w = shapes["pos_x_ends"][0]
    if n_w is not None and wp.get("n_walkers") not in (None, n_w):
        errs.append(f"walk_params.n_walkers {wp.get('n_walkers')} != N_w {n_w} of positions")

    # 5. tier flag => channels present
    env = meta.get("replay_envelope", {})
    for flag, chans in TIER_CHANNELS.items():
        if env.get(flag):
            if flag == "magnetization_transfer" and meta.get("mt"):
                continue
            missing = [c for c in chans if not _present(c, keys)]
            if missing:
                errs.append(f"replay_envelope.{flag}=true but missing channel(s) {missing} (SPEC §7/§13)")
    if env.get("field") and "susc_path_dct" in keys:
        if "field" not in (wp.get("qv") or {}) or "field_read_interval" not in wp:
            errs.append("susc_path_dct without walk_params.qv.field and field_read_interval (SPEC §9.5.5)")
        if (env.get("acquisition") or {}).get("max_refocusing_pulses") is None:
            errs.append("susc_path_dct without replay_envelope.acquisition.max_refocusing_pulses (SPEC §9.4 rule 5)")

    # 6. distributional codec MUST NOT carry per-walker channels (SPEC §9.4)
    if keys & DISTRIBUTIONAL_KEYS:
        bad = [c for c in ("compartment", "boundary_local_time", "bound", "spin_weights") if _present(c, keys)]
        if bad:
            errs.append(f"distributional codec carries per-walker channel(s) {bad}; "
                        f"forbidden — such packs are Gradient-only (SPEC §9.4)")
        if any(env.get(f) for f in ("bulk_relaxation", "surface_relaxivity", "magnetization_transfer")):
            errs.append("distributional codec must declare Gradient tier only (SPEC §9.4)")

    # 7. the occupancy columns (SPEC §9.5.7): static labels or event lists, N_w shared
    T_max = wp.get("T_max")
    for col in ("comp", "bound"):
        if f"{col}_static" in shapes and n_w is not None and shapes[f"{col}_static"][0] != n_w:
            errs.append(f"{col}_static length {shapes[f'{col}_static'][0]} != N_w={n_w}")
        if f"{col}_events_offsets" in events:
            off = np.asarray(events[f"{col}_events_offsets"], np.int64)
            t = np.asarray(events[f"{col}_events_t"], np.float64)
            lab = np.asarray(events[f"{col}_events_label"])
            if n_w is not None and off.size != n_w + 1:
                errs.append(f"{col}_events_offsets has {off.size} entries, not N_w+1={n_w + 1}")
                continue
            if off[0] != 0 or np.any(np.diff(off) < 1) or off[-1] != t.size or t.size != lab.size:
                errs.append(f"{col}_events: offsets must start at 0, give every walker >= 1 event, and end at "
                            f"N_ev={t.size}={lab.size}")
                continue
            first = t[off[:-1]]
            if np.any(first != 0.0):
                errs.append(f"{col}_events: every walker's first event must be at t = 0")
            walker = np.repeat(np.arange(off.size - 1), np.diff(off))
            same = walker[1:] == walker[:-1]
            if np.any(np.diff(t)[same] <= 0):
                errs.append(f"{col}_events: times of a walker must strictly increase")
            if np.any(lab[1:][same] == lab[:-1][same]):
                errs.append(f"{col}_events: an event repeats the label before it")
            if T_max is not None and np.any(t > T_max * (1 + 1e-9)):
                errs.append(f"{col}_events: a time exceeds T_max")
            if len(np.unique(walker)) == off.size - 1 and np.all(np.diff(off) == 1):
                errs.append(f"{col}: every walker has one event; MUST be stored as {col}_static")
    for c in ("boundary_local_time",):
        if n_w is not None and "blt_start" in shapes and shapes["blt_start"][0] != n_w:
            errs.append(f"blt_start length {shapes['blt_start'][0]} != N_w={n_w}")
    if "spin_weights" in shapes and n_w is not None and shapes["spin_weights"][0] != n_w:
        errs.append(f"spin_weights length {shapes['spin_weights'][0]} != N_w={n_w}")

    # 8. fidelity present (every pack is certified, SPEC §9.4 rule 4)
    if not fid:
        errs.append("a 'fidelity' self-report is REQUIRED on every pack (SPEC §9.4 rule 4)")

    # 9. license + citation
    for k in ("license", "citation"):
        if not meta.get(k):
            errs.append(f"missing required metadata '{k}' (SPEC §13)")

    return errs, warns


def main(path):
    errs, warns = check(path)
    with safe_open(path, framework="numpy") as f:
        _h = f.metadata() or {}
        meta = json.loads(_h.get("rpk") or "{}")
    env = meta.get("replay_envelope", {})
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
