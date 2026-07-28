# Examples

Minimal, dependency-light demonstrations of the Replay Pack format. Requirements:
`numpy`, `safetensors` (and `jsonschema` for full metadata validation). The `dmipy_sim_interop`
`--build` path additionally needs the dmipy-sim environment; everything else is standalone.

| File | What it shows |
|---|---|
| `write_minimal_rpk.py` | A **producer** in ~50 lines: emits a conformant Tier-0 pack (raw positions, identity codec) from an analytic free-diffusion walk. Proves any simulator can write `.rpk` with only NumPy + safetensors. |
| `reference_replayer.py` | A **replayer** in ~40 lines: `load_rpk` + the Gradient (§6.1) and Relaxation (§6.2) replay operations. No dmipy-sim. |
| `conformance_check.py` | A **validator** against SPEC §13 + the JSON schema. `python conformance_check.py pack.rpk`. |
| `demo_free_diffusion.py` | **End-to-end**: produce → validate → replay at several b-values, checked against exact `exp(-bD)` to the MC noise floor. Run this first. |
| `dmipy_sim_interop.py` | **Interop**: consume a pack written by dmipy-sim with the generic tools here (`python dmipy_sim_interop.py pack.rpk`), or `--build` a tiny one via dmipy-sim then validate+replay it. |

## Quick start

```bash
pip install numpy safetensors jsonschema
python demo_free_diffusion.py
```

Expected: `CONFORMANT ✓` and replayed signals matching `exp(-bD)` within ~`1/sqrt(N)`.
