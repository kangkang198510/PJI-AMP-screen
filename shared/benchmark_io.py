"""CSV contracts, metrics and deterministic empirical ROC for AMP benchmarks.

Scores must be oriented so larger values mean stronger AMP predictions. No
threshold is inferred, missing predictions are not discarded, and reported
metrics are fractions rather than percentages. No real study data are embedded.
"""
from __future__ import annotations

import csv
import hashlib
import importlib.metadata
import json
import math
import platform
from datetime import datetime, timezone
from pathlib import Path


def read_csv(path):
    path = Path(path)
    if not path.is_file():
        raise ValueError(f"Input file does not exist: {path.name}")
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, strict=True)
        fields = reader.fieldnames
        if not fields or len(fields) != len(set(fields)) or any(not x.strip() for x in fields):
            raise ValueError("CSV header must contain distinct nonempty names.")
        rows = []
        try:
            for row in reader:
                if None in row or any(v is None for v in row.values()):
                    raise ValueError("CSV row width differs from header.")
                rows.append({k: v.strip() for k, v in row.items()})
        except csv.Error as exc:
            raise ValueError(f"Malformed CSV: {exc}") from exc
    if not rows:
        raise ValueError("Input CSV contains no data rows.")
    return rows


def binary(value, name):
    if value not in ("0", "1"):
        raise ValueError(f"{name} must be encoded exactly as 0 or 1, with no missing values.")
    return int(value)


def load_predictions(path, require_predictions=True):
    rows = read_csv(path)
    required = {"peptide_id", "model", "y_true"}
    if require_predictions:
        required.add("y_pred")
    if not required.issubset(rows[0]):
        raise ValueError("Required columns: " + ", ".join(sorted(required)))
    by_model = {}
    seen = set()
    for row in rows:
        model, peptide = row["model"], row["peptide_id"]
        if not model or not peptide:
            raise ValueError("model and peptide_id cannot be empty.")
        key = (model, peptide)
        if key in seen:
            raise ValueError(f"Duplicate model/peptide_id pair: {model}, {peptide}")
        seen.add(key)
        row["y_true"] = binary(row["y_true"], "y_true")
        if require_predictions:
            row["y_pred"] = binary(row["y_pred"], "y_pred")
        by_model.setdefault(model, []).append(row)
    reference = None
    for model, model_rows in by_model.items():
        model_rows.sort(key=lambda x: x["peptide_id"])
        cohort = {r["peptide_id"]: r["y_true"] for r in model_rows}
        if reference is None:
            reference = cohort
        elif cohort != reference:
            raise ValueError(f"Cohort/ground-truth mismatch for {model}; identical IDs and labels required.")
    return dict(sorted(by_model.items()))


def load_counts(path):
    rows = read_csv(path)
    if not {"model", "TP", "TN", "FP", "FN"}.issubset(rows[0]):
        raise ValueError("Count CSV requires model,TP,TN,FP,FN.")
    seen = set()
    for row in rows:
        if not row["model"] or row["model"] in seen:
            raise ValueError("Count CSV model names must be nonempty and unique.")
        seen.add(row["model"])
        for field in ("TP", "TN", "FP", "FN"):
            try:
                value = int(row[field])
            except ValueError as exc:
                raise ValueError(f"{field} must be a nonnegative integer.") from exc
            if value < 0:
                raise ValueError(f"{field} cannot be negative.")
            row[field] = value
        if sum(row[k] for k in ("TP", "TN", "FP", "FN")) == 0:
            raise ValueError("Count table must contain at least one observation per model.")
    return rows


def counts_from_predictions(rows):
    counts = dict(TP=0, TN=0, FP=0, FN=0)
    names = {(1, 1): "TP", (0, 0): "TN", (0, 1): "FP", (1, 0): "FN"}
    for row in rows:
        counts[names[row["y_true"], row["y_pred"]]] += 1
    return counts


def calculate_metrics(tp, tn, fp, fn):
    values = (tp, tn, fp, fn)
    if any(isinstance(x, bool) or not isinstance(x, int) or x < 0 for x in values):
        raise ValueError("Confusion counts must be nonnegative integers.")
    ratio = lambda numerator, denominator: numerator / denominator if denominator else None
    mcc_denom = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    metrics = {
        "n": sum(values), "TP": tp, "TN": tn, "FP": fp, "FN": fn,
        "sensitivity": ratio(tp, tp + fn), "specificity": ratio(tn, tn + fp),
        "accuracy": ratio(tp + tn, sum(values)), "precision": ratio(tp, tp + fp),
        "F1": ratio(2 * tp, 2 * tp + fp + fn),
        "MCC": ratio(tp * tn - fp * fn, mcc_denom),
    }
    metrics["undefined_metrics"] = ";".join(k for k, v in metrics.items() if v is None)
    return metrics


def metrics_table(by_model):
    output = []
    for model, rows in by_model.items():
        counts = counts_from_predictions(rows)
        output.append({"model": model, **calculate_metrics(*(counts[x] for x in ("TP", "TN", "FP", "FN")))})
    return output


def empirical_roc(rows):
    """Return AUC and full tied-threshold ROC. Binary-only scores are refused.

    The trapezoidal empirical ROC integral equals the probability a randomly
    selected positive outranks a negative, with half credit for score ties.
    No bootstrap CI, significance test, or interpolation is invented.
    """
    if {r["y_true"] for r in rows} != {0, 1}:
        raise ValueError("ROC requires both positive and negative ground-truth labels.")
    pairs = []
    for row in rows:
        try:
            score = float(row.get("score", ""))
        except (TypeError, ValueError) as exc:
            raise ValueError("ROC requires a numeric, nonmissing score for every peptide.") from exc
        if not math.isfinite(score):
            raise ValueError("ROC score must be finite for every peptide.")
        pairs.append((score, row["y_true"]))
    if set(s for s, _ in pairs).issubset({0.0, 1.0}):
        raise ValueError("Binary-only score values are not accepted as continuous ROC predictions.")
    positive = sum(label for _, label in pairs)
    negative = len(pairs) - positive
    points = [{"threshold": "inf", "FPR": 0.0, "TPR": 0.0}]
    grouped = {}
    for score, label in pairs:
        grouped.setdefault(score, [0, 0])[label] += 1
    tp = fp = 0
    for threshold in sorted(grouped, reverse=True):
        n_negative, n_positive = grouped[threshold]
        tp += n_positive
        fp += n_negative
        points.append({"threshold": threshold, "FPR": fp / negative, "TPR": tp / positive})
    area = sum((b["FPR"] - a["FPR"]) * (a["TPR"] + b["TPR"]) / 2
               for a, b in zip(points, points[1:]))
    return area, points


def write_csv(path, rows):
    path = Path(path)
    if not rows:
        raise ValueError("No output rows available.")
    if any(set(row) != set(rows[0]) for row in rows):
        raise ValueError("Output rows must have identical columns.")
    path.parent.mkdir(parents=True, exist_ok=True)
    # Exclusive creation protects existing files even if created after preflight.
    with path.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_audit(path, input_path, details, packages=()):
    versions = {"python": platform.python_version()}
    for package in packages:
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = "not installed"
    audit = {"timestamp_utc": datetime.now(timezone.utc).isoformat(),
             "input_sha256": hashlib.sha256(Path(input_path).read_bytes()).hexdigest(),
             "seed": 20261116, "random_sampling": False, "versions": versions, **details}
    payload = json.dumps(audit, ensure_ascii=False, indent=2, allow_nan=False) + "\n"
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with Path(path).open("x", encoding="utf-8") as handle:
        handle.write(payload)


def ensure_distinct_input_output(input_path, *outputs):
    """Preflight every artifact before any writes; never replace existing data."""
    original = Path(input_path).resolve()
    resolved = [Path(p).resolve() for p in outputs]
    if any(p == original for p in resolved):
        raise ValueError("Output paths must not overwrite the input file.")
    if len(resolved) != len(set(resolved)):
        raise ValueError("Output artifact paths must be distinct from each other.")
    for path in outputs:
        candidate = Path(path)
        if candidate.exists() or candidate.is_symlink():
            raise ValueError(f"Output already exists: {candidate.name}; choose a new output path/directory.")
