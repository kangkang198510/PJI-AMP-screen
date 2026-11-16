"""Synthetic unit fixtures only; these are not study observations."""
import csv
import math
import random
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from shared.benchmark_io import (calculate_metrics, empirical_roc, load_predictions, load_counts,
                                 read_csv, write_csv, write_audit, ensure_distinct_input_output)


class BenchmarkTests(unittest.TestCase):
    def test_known_confusion_counts(self):
        m = calculate_metrics(8, 7, 2, 3)
        self.assertEqual(m["n"], 20)
        self.assertAlmostEqual(m["sensitivity"], 8/11)
        self.assertAlmostEqual(m["specificity"], 7/9)
        self.assertAlmostEqual(m["accuracy"], .75)
        self.assertAlmostEqual(m["precision"], .8)
        self.assertAlmostEqual(m["F1"], 16/21)
        self.assertAlmostEqual(m["MCC"], 50/math.sqrt(10*11*9*10))

    def test_undefined_metrics_not_fabricated(self):
        m = calculate_metrics(0, 10, 0, 0)
        self.assertIsNone(m["sensitivity"])
        self.assertIsNone(m["precision"])
        self.assertIsNone(m["F1"])
        self.assertIsNone(m["MCC"])
        self.assertEqual(m["accuracy"], 1)

    def test_known_auc_and_ties(self):
        rows = [{"y_true": y, "score": s} for y, s in zip([0, 0, 1, 1], [.1, .4, .35, .8])]
        area, points = empirical_roc(rows)
        self.assertEqual(area, .75)
        self.assertEqual(points[-1]["TPR"], 1)
        self.assertEqual(points[-1]["FPR"], 1)
        self.assertEqual(empirical_roc([{"y_true": 0, "score": .5}, {"y_true": 1, "score": .5}])[0], .5)

    def test_binary_missing_nonfinite_and_single_class_scores_refused(self):
        fixtures = [([0, 1], [0, 1]), ([0, 1], [.2, ""]), ([0, 1], [.2, "nan"]), ([1, 1], [.2, .8])]
        for labels, scores in fixtures:
            with self.subTest(labels=labels, scores=scores), self.assertRaises(ValueError):
                empirical_roc([{"y_true": y, "score": s} for y, s in zip(labels, scores)])

    def test_roc_integral_matches_independent_pairwise_auc(self):
        rng = random.Random(17)
        for _ in range(50):
            pos = [rng.choice([.2, .4, .6, .8]) for _ in range(9)]
            neg = [rng.choice([.2, .4, .6, .8]) for _ in range(11)]
            expected = sum((p > n) + .5 * (p == n) for p in pos for n in neg) / (len(pos) * len(neg))
            rows = [{"y_true": 1, "score": x} for x in pos] + [{"y_true": 0, "score": x} for x in neg]
            self.assertAlmostEqual(empirical_roc(rows)[0], expected)

    def make_csv(self, directory, rows, fields):
        path = Path(directory) / "synthetic.csv"
        with path.open("w", newline="", encoding="utf-8") as handle:
            w = csv.DictWriter(handle, fieldnames=fields)
            w.writeheader()
            w.writerows(rows)
        return path

    def test_duplicates_and_mismatched_cohorts_refused(self):
        base = {"peptide_id": "synthetic_1", "model": "A", "y_true": "1", "y_pred": "0"}
        with tempfile.TemporaryDirectory() as d:
            for other in [base, {**base, "model": "B", "peptide_id": "synthetic_2"}, {**base, "model": "B", "y_true": "0"}]:
                path = self.make_csv(d, [base, other], list(base))
                with self.assertRaises(ValueError):
                    load_predictions(path)

    def test_model_order_does_not_change_cohort(self):
        rows = [{"peptide_id": str(i), "model": model, "y_true": str(i % 2), "y_pred": "0"}
                for model, ids in [("A", [0, 1]), ("B", [1, 0])] for i in ids]
        with tempfile.TemporaryDirectory() as d:
            data = load_predictions(self.make_csv(d, rows, list(rows[0])))
            self.assertEqual(len(data), 2)

    def test_count_values_must_be_nonnegative_integers(self):
        with self.assertRaises(ValueError):
            calculate_metrics(-1, 2, 3, 4)
        with tempfile.TemporaryDirectory() as d:
            row = dict(model="A", TP="1.5", TN="1", FP="1", FN="1")
            with self.assertRaises(ValueError):
                load_counts(self.make_csv(d, [row], list(row)))

    def test_malformed_csv_not_silently_parsed(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "bad.csv"
            path.write_text('a,b\n"unterminated,b\n', encoding="utf-8")
            with self.assertRaises(ValueError):
                read_csv(path)

    def test_preflight_protects_input_sidecars_and_existing_files(self):
        with tempfile.TemporaryDirectory() as d:
            original = Path(d) / "source.csv"
            original.write_text("original", encoding="utf-8")
            output = Path(d) / "result.csv"
            sidecar = Path(d) / "result.csv.audit.json"
            sidecar.write_text("keep audit", encoding="utf-8")
            for paths in [(output, original), (output, sidecar), (output, output)]:
                with self.assertRaises(ValueError):
                    ensure_distinct_input_output(original, *paths)
            self.assertEqual(original.read_text(), "original")
            self.assertEqual(sidecar.read_text(), "keep audit")
            self.assertFalse(output.exists())

    def test_writers_exclusively_create_files(self):
        with tempfile.TemporaryDirectory() as d:
            original = Path(d) / "source.csv"
            original.write_text("original", encoding="utf-8")
            with self.assertRaises(FileExistsError):
                write_csv(original, [{"value": 1}])
            with self.assertRaises(FileExistsError):
                write_audit(original, original, {})
            self.assertEqual(original.read_text(), "original")

    def test_cli_sidecar_collision_fails_before_primary_write(self):
        root = Path(__file__).resolve().parents[1]
        with tempfile.TemporaryDirectory() as d:
            output = Path(d) / "result.csv"
            sidecar = Path(str(output) + ".audit.json")
            sidecar.write_text("keep audit", encoding="utf-8")
            run = subprocess.run([sys.executable, str(root / "03_model_benchmarking/calculate_metrics.py"),
                                  "--input", str(root / "example_data/demo_benchmark.csv"), "--output", str(output)],
                                 capture_output=True, text=True)
            self.assertNotEqual(run.returncode, 0)
            self.assertIn("Output already exists", run.stderr)
            self.assertFalse(output.exists())
            self.assertEqual(sidecar.read_text(), "keep audit")


if __name__ == "__main__":
    unittest.main()
