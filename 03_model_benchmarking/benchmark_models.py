"""Benchmark already-exported AMP predictions; does not train or run AMP models."""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from shared.benchmark_io import load_predictions, empirical_roc, metrics_table, write_csv, write_audit, ensure_distinct_input_output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Long CSV peptide_id,model,y_true,y_pred[,score].")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--score-direction", choices=["higher-is-AMP"],
                        help="Required to evaluate any score column. Otherwise AUROC is left blank.")
    args = parser.parse_args()
    try:
        out = Path(args.output_dir)
        ensure_distinct_input_output(args.input, out / "model_metrics.csv", out / "roc_points.csv", out / "benchmark_audit.json")
        by_model = load_predictions(args.input)
        result = metrics_table(by_model)
        all_points = []
        for record in result:
            model = record["model"]
            record.update(AUROC=None, eligible_AUROC_gt_0_90=None, AUROC_status="not evaluated: score orientation not confirmed")
            if args.score_direction:
                try:
                    area, points = empirical_roc(by_model[model])
                    record.update(AUROC=area, eligible_AUROC_gt_0_90=area > 0.90, AUROC_status="computed")
                    all_points.extend({"model": model, **p} for p in points)
                except ValueError as exc:
                    record["AUROC_status"] = "not evaluated: " + str(exc)
        out.mkdir(parents=True, exist_ok=True)
        write_csv(out / "model_metrics.csv", result)
        if all_points:
            write_csv(out / "roc_points.csv", all_points)
        write_audit(out / "benchmark_audit.json", args.input,
                    {"analysis": "benchmark of exported model outputs", "models": len(by_model),
                     "score_direction": args.score_direction, "cohort_alignment_verified": True,
                     "selection_rule": "AUROC > 0.90; unavailable AUROC yields no selection decision",
                     "missing_rows_dropped": False, "CI_computed": False, "P_values_computed": False})
    except (ValueError, OSError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
