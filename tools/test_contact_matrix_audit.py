"""No-dependency synthetic regressions for contact_matrix_audit.py."""
import json
import math
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import contact_matrix_audit as auditor


class ContactSnapshotAuditTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def write(self, name, contents):
        file = self.root / name
        file.write_text(contents, encoding="utf-8")
        return file

    def symmetric_example(self):
        return self.write("stiffness.mtx", """%%MatrixMarket matrix coordinate real symmetric
% Three mechanical active DOFs
3 3 5
1 1 4
2 1 -1
2 2 4
3 2 -1
3 3 3
""")

    def test_symmetric_operator_true_residual_modes_and_partial_nodes(self):
        a = self.symmetric_example()
        b = self.write("rhs.mtx", """%%MatrixMarket matrix array real general
3 1
2
4
7
""")
        x = self.write("solution.mtx", """%%MatrixMarket matrix array real general
3 1
1
2
3
""")
        modes = self.write("rbm.mtx", """%%MatrixMarket matrix array real general
3 2
1
1
1
1
0
-1
""")
        dof = self.write("dofs.csv", """eq,node,dof,body
0,10,1,part-A
1,10,2,part-A
2,11,1,part-B
""")
        report = auditor.audit(a, b, None, x, modes, dof)
        op = report["operator"]
        self.assertEqual(op["n"], 3)
        self.assertEqual(op["nnz"], 7)
        self.assertEqual(op["zero_diagonal_count"], 0)
        self.assertEqual(op["min_nonzero_abs_diagonal"], 3.0)
        self.assertEqual(op["frobenius_relative_asymmetry"], 0.0)
        self.assertAlmostEqual(report["solution_relative_residual"], 0.0)
        self.assertAlmostEqual(
            report["initial_relative_residual"], 1.0
        )
        self.assertEqual(len(report["near_nullspace_modes"]), 2)
        self.assertAlmostEqual(
            report["near_nullspace_modes"][0]["rayleigh_quotient"], 7.0 / 3.0
        )
        self.assertEqual(report["dof_map"]["physical_nodes"], 2)
        self.assertEqual(report["dof_map"]["structural_bodies"], 2)
        self.assertEqual(report["dof_map"]["nodes_with_incomplete_dofs"], 2)

    def test_nonsymmetric_operator_and_initial_guess(self):
        a = self.write("unsymmetric.mtx", """%%MatrixMarket matrix coordinate real general
2 2 3
1 1 2
1 2 3
2 2 1
""")
        b = self.write("b.mtx", """%%MatrixMarket matrix array real general
2 1
2
1
""")
        x0 = self.write("x0.mtx", """%%MatrixMarket matrix array real general
2 1
1
0
""")
        report = auditor.audit(a, rhs_path=b, initial_path=x0)
        self.assertAlmostEqual(
            report["operator"]["frobenius_relative_asymmetry"],
            math.sqrt(18.0 / 14.0)
        )
        self.assertAlmostEqual(report["initial_relative_residual"], 1 / math.sqrt(5))

    def test_duplicate_coordinate_entries_are_summed(self):
        a = self.write("duplicates.mtx", """%%MatrixMarket matrix coordinate real general
1 1 2
1 1 2
1 1 3
""")
        report = auditor.audit(a)
        self.assertEqual(report["operator"]["nnz"], 1)
        self.assertEqual(report["operator"]["max_abs_diagonal"], 5.0)

    def test_invalid_nonfinite_incomplete_and_out_of_range(self):
        inputs = [
            """%%MatrixMarket matrix coordinate real general
1 1 1
1 1 nan
""",
            """%%MatrixMarket matrix coordinate real general
1 1 2
1 1 1
""",
            """%%MatrixMarket matrix coordinate real general
1 1 1
2 1 1
""",
        ]
        for n, body in enumerate(inputs):
            with self.subTest(case=n):
                f = self.write("bad%d.mtx" % n, body)
                with self.assertRaises(ValueError):
                    auditor.audit(f)

    def test_dof_map_requires_unique_complete_equations(self):
        a = self.symmetric_example()
        bad = self.write("bad.csv", """eq,node,dof,body
0,10,1,A
0,10,2,A
2,11,1,B
""")
        with self.assertRaises(ValueError):
            auditor.audit(a, dof_map_path=bad)
        bad.write_text("""eq,node,dof,body
0,10,1,A
1,10,2,A
""", encoding="utf-8")
        with self.assertRaises(ValueError):
            auditor.audit(a, dof_map_path=bad)

    def test_cli_success_and_nonzero_status_on_bad_input(self):
        a = self.symmetric_example()
        script = Path(auditor.__file__).resolve()
        ok = subprocess.run(
            [sys.executable, str(script), "--matrix", str(a)],
            check=False, capture_output=True, text=True
        )
        self.assertEqual(ok.returncode, 0, ok.stderr)
        self.assertEqual(json.loads(ok.stdout)["operator"]["n"], 3)
        bad = self.write("bad.mtx", "invalid\n")
        fail = subprocess.run(
            [sys.executable, str(script), "--matrix", str(bad)],
            check=False, capture_output=True, text=True
        )
        self.assertEqual(fail.returncode, 2)
        self.assertIn("contact-matrix-audit:", fail.stderr)


if __name__ == "__main__":
    unittest.main()
