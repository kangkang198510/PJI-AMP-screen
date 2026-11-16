# Vendor normalization algorithm has not been supplied. This is an import/validation gate.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
source(file.path(dirname(script), "..", "shared", "r_utils.R"))
a <- parse_cli(c("input", "samples", "out", "mode", "demo"), c("input", "samples", "out", "mode"),
  "Rscript normalization.R --input vendor_normalized.csv --samples samples.csv --mode vendor-normalized --out results/normalized [--demo]\nRequires an ALREADY NORMALIZED, nonnegative linear abundance matrix. Performs validation/export only; normalization parameters are unavailable.")
if (a$mode != "vendor-normalized") stop("Only vendor-normalized import is supported. Supply the confirmed vendor algorithm before implementing raw-data normalization.")
set.seed(20261116)
s <- read_samples(a$samples, isTRUE(a$demo)); x <- read_matrix(a$input, s)
mat <- as.matrix(x[, s$sample_id, drop = FALSE]); d <- out_dir(a$out)
audit(d, a, c(a$input, a$samples), c("NO normalization, imputation, or log transform performed.",
 "The supplied table must already contain vendor-normalized LINEAR abundances. Vendor method remains to be confirmed.",
 "Metabolomics internal-standard normalization and half-minimum imputation are not transferred to peptidomics."))
write_table(x, file.path(d, "normalized_matrix.csv"))
qc <- data.frame(sample_id = s$sample_id, group = s$group, n_observed = colSums(!is.na(mat)),
 median = apply(mat, 2, function(v) if (all(is.na(v))) NA_real_ else median(v, na.rm = TRUE)),
 zero_count = colSums(mat == 0, na.rm = TRUE))
write_table(qc, file.path(d, "normalization_qc.csv"))
