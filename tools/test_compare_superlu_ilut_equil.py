"""Smoke tests for the optional SciPy/SuperLU ablation utility."""

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

try:
    import numpy as np
    from scipy import sparse
    from scipy.io import mmwrite
except ImportError:
    np = None

SCRIPT = Path(__file__).with_name("compare_superlu_ilut_equil.py")


@unittest.skipIf(np is None, "NumPy and SciPy required for the optional benchmark")
class TestSuperLUAblation(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.directory = Path(self.tmp.name)
        n = 40
        a = sparse.diags(
            [-0.2 * np.ones(n - 1), 4.0 * np.ones(n),
             -0.4 * np.ones(n - 1)], [-1, 0, 1], format="lil"
        )
        # A genuine nontrivial column correction is required so a no-op
        # SuperLU Equil option cannot accidentally pass this regression.
        a[4, 5] = 12.0
        a = sparse.diags(np.logspace(-4, 4, n)) @ a.tocsc()
        self.a = self.directory / "snapshot_A.mtx"
        self.b = self.directory / "snapshot_b.mtx"
        mmwrite(self.a, a)
        mmwrite(self.b, (a @ np.linspace(1, 2, n)).reshape(-1, 1))
        self.csv = self.directory / "ablation.csv"

    def invoke(self, *flags):
        return subprocess.run(
            [sys.executable, str(SCRIPT), str(self.a), "--csv", str(self.csv), *flags],
            capture_output=True, text=True,
        )

    def test_four_modes_with_original_rhs(self):
        run = self.invoke("--rcm", "--row-max-scale")
        self.assertEqual(run.returncode, 0, run.stderr)
        with self.csv.open(newline="") as fp:
            rows = list(csv.DictReader(fp))
        self.assertEqual(len(rows), 4)
        self.assertEqual(
            {(r["permc_spec"], r["equil"]) for r in rows},
            {("NATURAL", "False"), ("NATURAL", "True"),
             ("COLAMD", "False"), ("COLAMD", "True")},
        )
        for row in rows:
            self.assertEqual(row["status"], "converged", row)
            self.assertLess(float(row["rel_residual"]), 1e-8)
            self.assertEqual(row["rhs_source"], str(self.b))
            self.assertEqual(row["row_max_scale"], "True")
            if row["equil"] == "True":
                self.assertGreater(int(row["col_scaled_count"]), 0, row)
                self.assertGreater(float(row["col_scale_max"]), 1.0, row)
            else:
                self.assertEqual(int(row["col_scaled_count"]), 0, row)

    def test_missing_rhs_requires_explicit_manufacturing(self):
        self.b.unlink()
        run = self.invoke("--rcm")
        self.assertNotEqual(run.returncode, 0)
        run = self.invoke("--manufactured-rhs")
        self.assertEqual(run.returncode, 0, run.stderr)
        with self.csv.open(newline="") as fp:
            rows = list(csv.DictReader(fp))
        self.assertTrue(all(r["rhs_source"].startswith("manufactured") for r in rows))


if __name__ == "__main__":
    unittest.main()
