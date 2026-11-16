# Synthetic regression checks only; no patient data or manuscript-result expectations.
# Run: Rscript tests/test_r_pipeline.R
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
root <- file.path(dirname(script), "..")
rscript <- file.path(R.home("bin"), "Rscript.exe")
if (!file.exists(rscript)) rscript <- file.path(R.home("bin"), "Rscript")
work <- tempfile("pji_r_validation_"); dir.create(work)
set.seed(20261116)
write_input <- function(x, name) {
 p <- file.path(work, name); write.csv(x, p, row.names = FALSE, na = ""); p
}
run <- function(entry, opts, success = TRUE) {
 log <- tempfile(tmpdir = work, fileext = ".log")
 status <- system2(rscript, c(shQuote(file.path(root, entry)), vapply(opts, shQuote, character(1))),
                   stdout = log, stderr = log)
 if ((status == 0) != success) stop(paste(readLines(log, warn = FALSE), collapse = "\n"))
 invisible(status)
}
samples <- data.frame(sample_id = c(paste0("P", 1:4), paste0("O", 1:4)), group = rep(c("PJI", "OA"), each = 4))
sp <- write_input(samples, "samples.csv")
mat <- rbind(c(10,11,12,13,1,2,3,4), rep(2,8), c(1,2,NA,NA,3,4,NA,NA), c(1,NA,NA,NA,2,3,NA,NA), rep(0,8))
raw <- data.frame(peptide_id = paste0("p", 1:5), sequence = c("SYNTHA", "SYNTHB", "SYNTHC", "SYNTHD", "SYNTHA"), mat, check.names = FALSE)
names(raw)[3:10] <- samples$sample_id
rp <- write_input(raw, "raw.csv"); np <- write_input(raw, "norm.csv")
run("01_peptidomics_preprocessing/missingness_filter.R", c("--input", rp, "--samples", sp, "--out", file.path(work,"missing"), "--demo"))
qc <- read.csv(file.path(work, "missing", "missingness_qc.csv"))
stopifnot(identical(qc$retained, c(TRUE,TRUE,TRUE,FALSE,TRUE)), qc$missing_fraction[3] == 0.5)
run("01_peptidomics_preprocessing/missingness_filter.R", c("--input", rp, "--samples", sp, "--out", file.path(work,"wrong_n")), success = FALSE)
run("01_peptidomics_preprocessing/normalization.R", c("--input", np, "--samples", sp, "--mode", "vendor-normalized", "--out", file.path(work,"normalization"), "--demo"))
norm <- read.csv(file.path(work,"normalization","normalized_matrix.csv"), check.names = FALSE)
stopifnot(identical(norm, read.csv(np, check.names = FALSE)))
common <- c("--raw", rp, "--input", np, "--samples", sp, "--method", "asymptotic", "--continuity", "yes", "--bootstrap", "100", "--seed", "42", "--demo")
run("02_differential_analysis/SDPS_analysis.R", c(common, "--out", file.path(work,"sdps")))
run("02_differential_analysis/SDPS_analysis.R", c(common, "--out", file.path(work,"sdps_repeat")))
res <- read.csv(file.path(work,"sdps","sdps_all_results.csv"))
repeat_res <- read.csv(file.path(work,"sdps_repeat","sdps_all_results.csv"))
stopifnot(identical(res, repeat_res), nrow(res)==4, res$rank_biserial[1] == 1,
 res$effect_CI_low[1] == 1, res$effect_CI_high[1] == 1,
 abs(res$p_value[1] - 0.0303828219765775) < 1e-12,
 res$p_value[2] == 1, res$rank_biserial[2] == 0,
 all.equal(res$q_value, p.adjust(res$p_value, "BH")) == TRUE,
 is.na(res$FC[4]), !"p4" %in% res$peptide_id)
bad <- raw; bad$P1[1] <- -1; bp <- write_input(bad, "bad.csv")
run("01_peptidomics_preprocessing/normalization.R", c("--input", bp, "--samples", sp, "--mode", "vendor-normalized", "--out", file.path(work,"negative"), "--demo"), success = FALSE)
bad$P1[1] <- "oops"; bp <- write_input(bad, "bad_text.csv")
run("01_peptidomics_preprocessing/normalization.R", c("--input", bp, "--samples", sp, "--mode", "vendor-normalized", "--out", file.path(work,"text"), "--demo"), success = FALSE)
ebs <- data.frame(peptide_id = paste0("e",1:5), sequence = paste0("SYNTH",1:5), p_value = c(0.049,0.049,0.05,0.001,NA), FC = c(5,0.2,6,4,0.1))
ep <- write_input(ebs, "ebs.csv")
run("02_differential_analysis/EBS_screening.R", c("--input", ep, "--out", file.path(work,"ebs"), "--confirm-test", "two-sided-student"))
er <- read.csv(file.path(work,"ebs","ebs_all_results.csv"))
stopifnot(identical(er$selected, c(TRUE,TRUE,FALSE,FALSE,FALSE)))
cat("PASS: raw missingness boundary; study group sizes; unchanged vendor-normalized values;\n")
cat("Mann-Whitney reference P; rank-biserial direction; ties/zero means; BH; bootstrap reproducibility;\n")
cat("negative/nonnumeric rejection; EBS P/FC boundaries. Synthetic fixture only.\n")
cat(R.version.string, "\n")
