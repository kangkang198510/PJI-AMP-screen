# Apply missingness only to the ORIGINAL quantitative matrix, before imputation.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
source(file.path(dirname(script), "..", "shared", "r_utils.R"))
a <- parse_cli(c("input", "samples", "out", "demo"), c("input", "samples", "out"),
  "Rscript missingness_filter.R --input raw.csv --samples samples.csv --out results/missingness [--demo]\nRetain missing fraction <= 0.50 across all manifest samples. NA/blank means missing; zero remains measured zero.")
set.seed(20261116)
s <- read_samples(a$samples, isTRUE(a$demo)); x <- read_matrix(a$input, s)
mat <- as.matrix(x[, s$sample_id, drop = FALSE]); fraction <- rowMeans(is.na(mat))
qc <- data.frame(peptide_id = x$peptide_id, n_samples = ncol(mat), n_observed = rowSums(!is.na(mat)),
                 missing_fraction = fraction, retained = fraction <= 0.5)
d <- out_dir(a$out)
audit(d, a, c(a$input, a$samples), "Fixed source-backed threshold: exclude missingness > 50%; no imputation; zeros not reclassified.")
write_table(qc, file.path(d, "missingness_qc.csv"))
write_table(x[qc$retained, , drop = FALSE], file.path(d, "filtered_matrix.csv"))
write_table(data.frame(sample_id = s$sample_id, group = s$group, missing_fraction = colMeans(is.na(mat))),
            file.path(d, "sample_qc.csv"))
message(sum(qc$retained), " / ", nrow(x), " peptides retained; no expected result count enforced.")
