"""Plot empirical ROC curves from complete scores using Python only.

Figure contract: compare discrimination of exported AMP model scores on one
shared labeled cohort, using ROC curves plus empirical AUROC. This is a single
quantitative panel, 180 x 140 mm, with an external legend, complete axes, sample
counts and redundant line styles. It makes no claim of external validity or
statistical superiority. Outputs PDF/SVG, 600 dpi TIFF, PNG preview and source points.
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from shared.benchmark_io import load_predictions, empirical_roc, write_csv, write_audit, ensure_distinct_input_output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Long CSV peptide_id,model,y_true,score.")
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--score-direction", choices=["higher-is-AMP"], required=True)
    args = parser.parse_args()
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        parser.error("ROC plotting requires matplotlib and Pillow; install repository requirements.txt.")
    try:
        ensure_distinct_input_output(args.input, *[args.output_prefix + x for x in
            (".pdf", ".svg", ".tiff", ".png", ".points.csv", ".audit.json")])
        by_model = load_predictions(args.input, require_predictions=False)
        # Validate every model before opening a graphics device or writing output.
        rocs = {model: empirical_roc(rows) for model, rows in by_model.items()}
        plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
                             "font.size": 7, "svg.fonttype": "none", "pdf.fonttype": 42,
                             "axes.spines.top": False, "axes.spines.right": False,
                             "axes.linewidth": 0.8, "legend.frameon": False})
        palette = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#333333"]
        styles = ["-", "--", "-.", ":"]
        markers = [None, "o", "s", "^", "D"]
        fig, ax = plt.subplots(figsize=(180 / 25.4, 140 / 25.4))
        all_points = []
        for index, (model, (area, points)) in enumerate(rocs.items()):
            ax.plot([p["FPR"] for p in points], [p["TPR"] for p in points],
                    color=palette[index % len(palette)], linestyle=styles[index % len(styles)],
                    marker=markers[(index // len(styles)) % len(markers)], markersize=2.5,
                    markevery=max(1, len(points) // 12),
                    linewidth=1.2, label=f"{model} (AUROC {area:.3f})")
            all_points.extend({"model": model, **p} for p in points)
        ax.plot([0, 1], [0, 1], color="0.65", linestyle="--", linewidth=0.7)
        first = next(iter(by_model.values()))
        n_positive = sum(r["y_true"] for r in first)
        ax.set(xlim=(0, 1), ylim=(0, 1), xlabel="False-positive rate (1 − specificity)",
               ylabel="True-positive rate (sensitivity)",
               title=f"AMP benchmark: n = {len(first)} peptides ({n_positive} AMP / {len(first)-n_positive} non-AMP)")
        ax.set_aspect("equal", adjustable="box")
        ax.legend(loc="center left", bbox_to_anchor=(1.03, 0.5), fontsize=7)
        fig.tight_layout()
        Path(args.output_prefix).parent.mkdir(parents=True, exist_ok=True)
        for extension in ("pdf", "svg", "tiff"):
            with Path(args.output_prefix + "." + extension).open("xb") as handle:
                fig.savefig(handle, format=extension, dpi=600, bbox_inches="tight")
        with Path(args.output_prefix + ".png").open("xb") as handle:
            fig.savefig(handle, format="png", dpi=150, bbox_inches="tight")
        plt.close(fig)
        write_csv(args.output_prefix + ".points.csv", all_points)
        write_audit(args.output_prefix + ".audit.json", args.input,
                    {"analysis": "empirical ROC figure", "models": len(by_model), "n_peptides": len(first),
                     "score_direction": args.score_direction, "CI_computed": False, "P_values_computed": False,
                     "statistical_method": "empirical ROC; trapezoidal AUROC; tied scores grouped"},
                    packages=("matplotlib", "Pillow"))
    except (ValueError, OSError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
