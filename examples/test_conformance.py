"""Battery contract regression tests. Run: python -m unittest discover -s examples."""
import copy
import json
from pathlib import Path
import tempfile
import unittest

import jsonschema
import numpy as np
from safetensors.numpy import save_file

import conformance_check as cc


BATTERY = {
    "bvals": [0.0, 0.5e9, 1e9, 2e9, 3e9],
    "dirs": [[0, 0, 1], [1, 0, 0], [1, 0, 1]],
    "ogse_periods": [1, 2, 3, 5], "shortd_b": 1e9,
    "shortd_deltas_frac": [0.2, 0.1, 0.05, 0.025],
    "delta_frac": 0.2, "Delta_frac": 0.5,
}


class BatteryConformance(unittest.TestCase):
    def setUp(self):
        self.meta = {
            "rpk_schema_version": "0.4", "id": "tests/battery",
            "walk_params": {
                "n_walkers": 2, "T_max": 0.1, "diffusivity": 2e-9, "seed": 0,
                "block": {"h": 1e-5, "m": 1, "precision": "float64"},
                "qv": {"positions": [4e-10, 4e-10, 4e-10]},
            },
            "compression": {"method": "bridge_dst", "K": 2, "walker_preserving": True},
            "replay_envelope": {"gradient": True, "diffusivity_fixed": True},
            "fidelity": {"err_max": 0, "floor_max": 0.1,
                         "terms": {"truncation": 0, "quadrature": {"positions": 0}},
                         "battery": copy.deepcopy(BATTERY)},
            "license": "CC0-1.0", "citation": "Battery regression fixture",
        }
        # Keep generated files within the repository, including on CI/GPU runners.
        self.temp = tempfile.TemporaryDirectory(dir=Path(__file__).parent)
        self.addCleanup(self.temp.cleanup)
        self.path = str(Path(self.temp.name) / "battery.rpk")

    def check_metadata(self):
        save_file({f"pos_{a}": np.zeros((2, 4), dtype=np.float32) for a in "xyz"},
                  self.path, metadata={"rpk": json.dumps(self.meta)})
        errors, warnings = cc.check(self.path)
        self.assertEqual(warnings, [])
        return errors

    def assert_battery_error(self):
        errors = self.check_metadata()
        self.assertTrue(errors)
        self.assertTrue(any("fidelity.battery" in e for e in errors), errors)

    def test_reference_object_passes(self):
        self.assertEqual(self.check_metadata(), [])

    def test_extension_passes(self):
        self.meta["fidelity"]["battery"]["theta_deg"] = [0, 45, 90]
        self.assertEqual(self.check_metadata(), [])

    def test_string_fails_by_name(self):
        self.meta["fidelity"]["battery"] = "default"
        self.assert_battery_error()

    def test_nonobject_batteries_fail(self):
        for value in (None, [], 1, True):
            with self.subTest(value=value):
                self.meta["fidelity"]["battery"] = value
                self.assert_battery_error()

    def test_missing_battery_passes(self):
        del self.meta["fidelity"]["battery"]
        self.assertEqual(self.check_metadata(), [])

    def test_missing_parameters(self):
        for key in BATTERY:
            with self.subTest(key=key):
                self.meta["fidelity"]["battery"] = copy.deepcopy(BATTERY)
                del self.meta["fidelity"]["battery"][key]
                self.assert_battery_error()

    def test_malformed_parameters(self):
        cases = [("bvals", "1e9"), ("bvals", [-1]), ("bvals", []),
                 ("dirs", [[1, 0]]), ("dirs", [[1, "x", 0]]), ("dirs", []),
                 ("ogse_periods", [1.5]), ("ogse_periods", [0]),
                 ("shortd_b", -1), ("shortd_b", "1e9"),
                 ("shortd_deltas_frac", [0]), ("shortd_deltas_frac", [1.1]),
                 ("delta_frac", 0), ("delta_frac", "0.2"), ("Delta_frac", 1.1)]
        for key, value in cases:
            with self.subTest(key=key, value=value):
                self.meta["fidelity"]["battery"] = copy.deepcopy(BATTERY)
                self.meta["fidelity"]["battery"][key] = value
                self.assert_battery_error()

    def test_unused_families_pass(self):
        self.meta["fidelity"]["battery"]["ogse_periods"] = []
        self.meta["fidelity"]["battery"]["shortd_deltas_frac"] = []
        self.assertEqual(self.check_metadata(), [])

    def test_schema_is_valid(self):
        with open(cc.SCHEMA) as source:
            jsonschema.Draft202012Validator.check_schema(json.load(source))


if __name__ == "__main__":
    unittest.main()
