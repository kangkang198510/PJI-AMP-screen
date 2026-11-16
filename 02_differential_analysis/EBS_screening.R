# Exploratory screen of existing vendor Student t-test statistics; independent of SDPS.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
source(file.path(dirname(script), "..", "shared", "r_utils.R"))
a <- parse_cli(c("input", "out", "p-column", "fc-column", "confirm-test"),
 c("input", "out", "confirm-test"),
 "Rscript EBS_screening.R --input vendor_stats.csv --out results/ebs --confirm-test two-sided-student [--p-column p_value] [--fc-column FC]\nInput: peptide_id, sequence, p_value, FC (PJI/OA ratio). p < 0.05 AND (FC >= 5 OR FC <= 0.2). Do not supply SDPS-prefiltered statistics.")
if (a[["confirm-test"]] != "two-sided-student") stop("Vendor P values must be confirmed as two-sided Student's t-test results")
set.seed(20261116)
x <- read_csv_strict(a$input); pc <- if (is.null(a[["p-column"]])) "p_value" else a[["p-column"]]
fc <- if (is.null(a[["fc-column"]])) "FC" else a[["fc-column"]]
need_cols(x, c("peptide_id", "sequence", pc, fc))
if (anyNA(x$peptide_id) || anyDuplicated(x$peptide_id) || anyNA(x$sequence)) stop("peptide_id must be unique/nonmissing; sequence required")
x$p_value <- numeric_strict(x[[pc]], pc); x$FC <- numeric_strict(x[[fc]], fc, nonnegative = TRUE)
if (any(x$p_value < 0 | x$p_value > 1, na.rm = TRUE)) stop("P values must be between 0 and 1")
x$log2FC <- ifelse(x$FC > 0, log2(x$FC), NA_real_)
x$selected <- !is.na(x$p_value) & !is.na(x$FC) & x$p_value < 0.05 & (x$FC >= 5 | x$FC <= 0.2)
x$status <- ifelse(is.na(x$p_value) | is.na(x$FC), "not_evaluable", "evaluated")
x$workflow <- "EBS"; d <- out_dir(a$out)
audit(d, a, a$input, c("EBS uses unadjusted P < 0.05 and FC >= 5 or <= 0.2; exploratory, not FDR-controlled.",
 "Screen only: original vendor matrix preprocessing and t-test input scale remain unverified.",
 "Zero FC is retained for screening but log2FC is undefined for plotting; missing P/FC marked not_evaluable."))
write_table(x, file.path(d, "ebs_all_results.csv"))
write_table(x[x$selected, , drop = FALSE], file.path(d, "ebs_candidates.csv"))
