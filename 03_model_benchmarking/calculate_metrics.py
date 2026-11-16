"""Calculate AMP benchmark metrics from exported predictions or confusion counts."""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from shared.benchmark_io import (calculate_metrics, load_counts, load_predictions,
                                 metrics_table, write_csv, write_audit, ensure_distinct_input_output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Long prediction CSV; or count CSV with --counts.")
    parser.add_argument("--output", required=True, help="Output metrics CSV (undefined values are blank).")
    parser.add_argument("--counts", action="store_true", help="Input columns model,TP,TN,FP,FN; cohort identity cannot be verified.")
    args = parser.parse_args()
    try:
        ensure_distinct_input_output(args.input, args.output, args.output + ".audit.json")
        if args.counts:
            rows = [{"model": r["model"], **calculate_metrics(*(r[k] for k in ("TP", "TN", "FP", "FN")))}
                    for r in load_counts(args.input)]
        else:
            rows = metrics_table(load_predictions(args.input))
        write_csv(args.output, rows)
        write_audit(args.output + ".audit.json", args.input,
                    {"analysis": "confusion metrics", "input_mode": "counts" if args.counts else "predictions",
                     "cohort_alignment_verified": not args.counts, "models": len(rows),
                     "undefined_metrics_policy": "blank, with undefined_metrics column", "AUROC_computed": False})
    except (ValueError, OSError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
