"""Pure-Python regression for LIS native contact benchmark audit."""
from pathlib import Path
import argparse
import tempfile
import unittest

from run_native_ilutp_contact import (
    integer_nonnegative, parse_full_row, parse_time_file, sha256_file,
    positive_number, nonnegative_number,
)


class NativeContactAuditTest(unittest.TestCase):
    def test_full_row_selection_is_strict(self):
        csv = 'variant,status,physical_relres,api_error,solver_status\n'
        full = 'ilutp_full,converged,3e-14,0,0\n'
        row, error = parse_full_row(csv + full)
        self.assertIsNone(error)
        self.assertEqual(row['physical_relres'], '3e-14')
        _, error = parse_full_row(csv + full + full)
        self.assertIn('expected one', error)
        _, error = parse_full_row(csv + 'legacy_ilut,converged,1e-4,0,0\n')
        self.assertIn('expected one', error)

    def test_integer_limit_rejects_fraction_and_overflow(self):
        self.assertEqual(integer_nonnegative('20000000'), 20000000)
        self.assertEqual(integer_nonnegative('2147483648'), 2147483648)
        self.assertEqual(integer_nonnegative('9223372036854775807'), 2**63 - 1)
        for item in ('-1', '+1', '1.2', '1e3', 'nan', '',
                     '9223372036854775808'):
            with self.assertRaises(argparse.ArgumentTypeError):
                integer_nonnegative(item)

    def test_finite_numeric_validation(self):
        self.assertEqual(positive_number('1e-5'), 1e-5)
        for item in ('nan', 'inf', '0', '-1'):
            with self.assertRaises(argparse.ArgumentTypeError):
                positive_number(item)
        self.assertEqual(nonnegative_number('0'), 0)

    def test_hash_and_time_report(self):
        with tempfile.TemporaryDirectory() as d:
            file = Path(d) / 'data.mtx'
            file.write_bytes(b'%%MatrixMarket\n1 1 1\n1 1 7\n')
            first = sha256_file(file)
            self.assertEqual(first, sha256_file(file))
            file.write_bytes(b'%%MatrixMarket\n1 1 1\n1 1 9\n')
            self.assertNotEqual(first, sha256_file(file))
            time_log = Path(d) / 'time'
            time_log.write_text('MAX_RSS_KIB=57376\nWALL_SECONDS=0.11\n')
            self.assertEqual(parse_time_file(time_log)['peak_rss_kib'], 57376)
            self.assertAlmostEqual(parse_time_file(time_log)['gnu_wall_seconds'], 0.11)


if __name__ == '__main__':
    unittest.main()