# Examples

Minimal, dependency-light demonstrations of the Replay Pack format. Requirements:
`numpy`, `safetensors` (and `jsonschema` for full metadata validation). The `dmipy_sim_interop`
`--build` path additionally needs the dmipy-sim environment; everything else is standalone.

| File | What it shows |
|---|---|
| `write_minimal_rpk.py` | A **producer** in ~100 lines: accumulates an analytic free-diffusion walk into the bands of SPEC §9.5 while it runs, records the quadratic variation and the certificate terms, and emits a conformant Tier-0 pack. Proves any simulator can write `.rpk` with only NumPy + safetensors. |
| `reference_replayer.py` | A **replayer**: `load_rpk` + the Gradient (§6.1, closed form on the bands) and Relaxation (§6.2, exact over the occupancy events) replay operations. No dmipy-sim. |
| `conformance_check.py` | A **validator** against SPEC §13 + the JSON schema; it refuses a pack carrying `n_t`, `dt_traj`, `compression.basis`, the identity codec or a save-grid tensor by name, with the instruction to re-walk it (§9.5.1). `python conformance_check.py pack.rpk`. |
| `demo_free_diffusion.py` | **End-to-end**: produce → validate → replay at several b-values, checked against exact `exp(-bD)` to the MC noise floor. Run this first. |
| `continuous_basis_check.py` | **Numerical check of SPEC §9.5**: the integration-by-parts identity between the bands and the cosine projection of the increments, the free-bridge band variance, the K-band series at an instant, the block-quadrature table `rel(m, ωh)` for `m = 0 … 3`, and the registry's derivation from the DST-I of instants. NumPy only. |
| `dmipy_sim_interop.py` | **Interop**: consume a pack written by dmipy-sim with the generic tools here (`python dmipy_sim_interop.py pack.rpk`), or `--build` a tiny one via dmipy-sim then validate+replay it. |

## Quick start

```bash
pip install numpy safetensors jsonschema
python demo_free_diffusion.py
```

Expected: `CONFORMANT ✓` and replayed signals matching `exp(-bD)` within ~`1/sqrt(N)`.

## Standalone checks

With the dependencies above installed, run from the repository root:

```bash
python -m unittest discover -s examples -v
python examples/continuous_basis_check.py
python examples/demo_free_diffusion.py
```

These checks also run in `.github/workflows/checks.yml`. The battery regression tests
use repository-local temporary packs, require no GPU, git commands or network access,
and can use the same command on a GPU runner. They require `jsonschema` so schema
validation cannot silently be skipped. `fidelity.battery` is optional for standalone packs and required for inherited
certification and merge inputs. When present it is a parameter object (§9.4),
including SI b-values in s/m²; strings and malformed parameters fail with a named path.
