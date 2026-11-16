"""Compute empirical AUROC only from complete continuous exported model scores."""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from shared.benchmark_io import load_predictions, empirical_roc, write_csv, write_audit, ensure_distinct_input_output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Long CSV: peptide_id,model,y_true,score.")
    parser.add_argument("--output", required=True, help="AUROC CSV; ROC points saved beside it.")
    parser.add_argument("--score-direction", choices=["higher-is-AMP"], required=True,
                        help="Explicitly confirm orientation from each model's documentation.")
    args = parser.parse_args()
    try:
        points_path = str(Path(args.output).with_suffix(".points.csv"))
        ensure_distinct_input_output(args.input, args.output, points_path, args.output + ".audit.json")
        by_model = load_predictions(args.input, require_predictions=False)
        metrics, all_points = [], []
        for model, rows in by_model.items():
            auc, points = empirical_roc(rows)
            metrics.append({"model": model, "n": len(rows), "AUROC": auc,
                            "eligible_AUROC_gt_0_90": auc > 0.90})
            all_points.extend({"model": model, **p} for p in points)
        write_csv(args.output, metrics)
        write_csv(points_path, all_points)
        write_audit(args.output + ".audit.json", args.input,
                    {"analysis": "empirical ROC trapezoidal AUC, half credit for ties",
                     "score_direction": args.score_direction, "models": len(by_model),
                     "cohort_alignment_verified": True, "missing_scores_dropped": False,
                     "CI_computed": False, "P_values_computed": False})
    except (ValueError, OSError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
