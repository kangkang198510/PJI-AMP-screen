# Figure contract: display EBS exploratory or SDPS FDR-controlled results with matching gates.
# A single quantitative panel; no claim of mechanism. PDF/SVG + 600 dpi TIFF.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
source(file.path(dirname(script), "..", "shared", "r_utils.R"))
a <- parse_cli(c("input", "out", "workflow"), c("input", "out", "workflow"),
 "Rscript volcano_plot.R --input results.csv --out results/volcano --workflow SDPS|EBS\nRequires peptide_id, log2FC, and q_value (SDPS) or p_value and FC (EBS). Includes all evaluable peptides, not just selected rows.")
if (!a$workflow %in% c("EBS", "SDPS")) stop("Workflow must be EBS or SDPS")
for (pkg in c("ggplot2", "svglite")) if (!requireNamespace(pkg, quietly = TRUE)) stop("Install R package: ", pkg)
set.seed(20261116)
x <- read_csv_strict(a$input); pcol <- if (a$workflow == "SDPS") "q_value" else "p_value"
need_cols(x, c("peptide_id", "log2FC", pcol))
if (anyNA(x$peptide_id) || anyDuplicated(x$peptide_id)) stop("peptide_id must be unique and present")
if ("workflow" %in% names(x) && any(is.na(x$workflow) | x$workflow != a$workflow)) stop("Input workflow differs from requested plot")
x$log2FC <- numeric_strict(x$log2FC, "log2FC"); x$probability <- numeric_strict(x[[pcol]], pcol)
if (any(x$probability < 0 | x$probability > 1, na.rm = TRUE)) stop("P/q values must be in [0,1]")
selected <- !is.na(x$probability) & x$probability < 0.05
if (a$workflow == "EBS") {
 need_cols(x, "FC"); x$FC <- numeric_strict(x$FC, "FC", nonnegative = TRUE)
 selected <- selected & !is.na(x$FC) & (x$FC >= 5 | x$FC <= 0.2)
}
x$category <- ifelse(selected, ifelse(x$log2FC > 0, "Higher in PJI", ifelse(x$log2FC < 0, "Lower in PJI", "Equal mean")), "Not selected")
x$plotted <- is.finite(x$log2FC) & is.finite(x$probability)
# A numerical floor only for rendering values exported as exact zero; not a changed P/q.
x$neg_log10 <- -log10(pmax(x$probability, .Machine$double.xmin))
if (!any(x$plotted)) stop("No finite points to plot")
shown <- x[x$plotted, , drop = FALSE]
caption <- paste0(a$workflow, ": ", if (a$workflow == "SDPS") "two-sided Mann-Whitney U; BH q < 0.05" else "vendor two-sided Student t test; P < 0.05; FC >= 5 or <= 0.2",
 "; peptides plotted = ", nrow(shown), "; excluded nonfinite = ", sum(!x$plotted), ".")
if (all(c("n_PJI", "n_OA") %in% names(shown))) {
 for (nm in c("n_PJI", "n_OA")) shown[[nm]] <- numeric_strict(shown[[nm]], nm, nonnegative = TRUE, allow_na = FALSE)
 caption <- paste0(caption, "\nAvailable biological observations per peptide: PJI n = ",
  paste(range(shown$n_PJI), collapse = "-"), "; OA n = ", paste(range(shown$n_OA), collapse = "-"), ".")
} else {
 caption <- paste0(caption, "\nBiological sample counts: not encoded in the supplied result table; confirm in source metadata.")
}
p <- ggplot2::ggplot(shown, ggplot2::aes(log2FC, neg_log10, colour = category, shape = category)) +
 ggplot2::geom_point(size = 1.1, alpha = 0.65) +
 ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.4) +
 ggplot2::scale_colour_manual(values = c("Higher in PJI" = "#D55E00", "Lower in PJI" = "#0072B2", "Equal mean" = "#333333", "Not selected" = "#999999")) +
 ggplot2::scale_shape_manual(values = c("Higher in PJI" = 17, "Lower in PJI" = 15, "Equal mean" = 4, "Not selected" = 16)) +
 ggplot2::labs(x = "log2 fold change (PJI / OA)", y = paste0("-log10(", pcol, ")"), colour = NULL, shape = NULL, caption = caption) +
 ggplot2::theme_classic(base_size = 8) + ggplot2::theme(legend.position = "top", plot.caption = ggplot2::element_text(size = 6, hjust = 0))
if (a$workflow == "EBS") p <- p + ggplot2::geom_vline(xintercept = c(-log2(5), log2(5)), linetype = "dotted", linewidth = 0.4)
d <- out_dir(a$out)
audit(d, a, a$input, c(caption, "Color plus shape encodes category. P/q=0 is floored to machine minimum only for plotting; source values remain unchanged.",
 "Input should be the full result table. Peptide count is not biological sample n. Per-peptide biological n is reported by SDPS analysis."))
write_table(x, file.path(d, "volcano_source_data.csv"))
save_plot(p, file.path(d, paste0("volcano_", a$workflow)))
