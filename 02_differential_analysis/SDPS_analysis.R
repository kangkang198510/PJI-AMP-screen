# Reconstruction of manuscript SDPS. R implementation choices are explicit, not claimed historical.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
source(file.path(dirname(script), "..", "shared", "r_utils.R"))
a <- parse_cli(c("raw", "input", "samples", "out", "method", "continuity", "bootstrap", "seed", "demo"),
 c("raw", "input", "samples", "out", "method", "continuity"),
 "Rscript SDPS_analysis.R --raw raw.csv --input vendor_normalized.csv --samples samples.csv --out results/sdps --method asymptotic --continuity yes [--bootstrap 2000] [--seed 20261116] [--demo]\nRaw matrix sets <=50% missingness gate; normalized nonnegative LINEAR values are tested. Two-sided Mann-Whitney, BH q<0.05, no FC gate. Method choices must be checked against historical settings.")
if (!a$method %in% c("asymptotic", "auto")) stop("--method must be asymptotic or auto (exact only if no ties and min(n)<=8)")
if (!a$continuity %in% c("yes", "no")) stop("--continuity must be yes or no")
B <- if (is.null(a$bootstrap)) 2000L else numeric_strict(a$bootstrap, "bootstrap", allow_na = FALSE)
seed <- if (is.null(a$seed)) 20261116L else numeric_strict(a$seed, "seed", allow_na = FALSE)
if (length(B) != 1 || B < 100 || B != as.integer(B)) stop("--bootstrap must be integer >=100")
if (seed < 0 || seed > .Machine$integer.max || seed != as.integer(seed)) stop("Invalid seed")
set.seed(seed)
s <- read_samples(a$samples, isTRUE(a$demo)); raw <- read_matrix(a$raw, s); x <- read_matrix(a$input, s)
if (!setequal(raw$peptide_id, x$peptide_id)) stop("Raw and normalized inputs must have identical peptide_id sets (unfiltered matrices)")
x <- x[match(raw$peptide_id, x$peptide_id), , drop = FALSE]
if (!identical(raw$sequence, x$sequence)) stop("Sequence mismatch between raw and normalized tables")
fraction <- rowMeans(is.na(as.matrix(raw[, s$sample_id, drop = FALSE])))
keep <- fraction <= 0.5
if (!any(keep)) stop("No peptides pass raw missingness filter")
qc <- data.frame(peptide_id = raw$peptide_id, missing_fraction_raw = fraction, retained = keep)
x <- x[keep, , drop = FALSE]
pj <- s$sample_id[s$group == "PJI"]; oa <- s$sample_id[s$group == "OA"]
# Rank-biserial correlation: 2 * P(PJI > OA, ties half-weight) - 1.
rank_effect <- function(p, o) 2 * mean(outer(p, o, ">") + 0.5 * outer(p, o, "==")) - 1
result <- lapply(seq_len(nrow(x)), function(i) {
 p <- as.numeric(x[i, pj]); p <- p[!is.na(p)]
 o <- as.numeric(x[i, oa]); o <- o[!is.na(o)]
 mp <- if (length(p)) mean(p) else NA_real_; mo <- if (length(o)) mean(o) else NA_real_
 fc <- if (is.finite(mo) && mo > 0) mp / mo else NA_real_
 z <- data.frame(peptide_id = x$peptide_id[i], sequence = x$sequence[i],
  n_PJI = length(p), n_OA = length(o), mean_PJI = mp, mean_OA = mo, FC = fc,
  log2FC = if (is.finite(fc) && fc > 0) log2(fc) else NA_real_, p_value = NA_real_,
  rank_biserial = NA_real_, effect_CI_low = NA_real_, effect_CI_high = NA_real_,
  test_method = NA_character_, status = "insufficient_observations")
 if (min(length(p), length(o)) < 2) return(z)
 exact <- a$method == "auto" && min(length(p), length(o)) <= 8 && !anyDuplicated(c(p, o))
 z$test_method <- if (exact) "exact" else paste0("asymptotic_continuity_", a$continuity)
 if (length(unique(c(p, o))) == 1L) { z$p_value <- 1; z$status <- "all_values_equal" } else {
  test <- wilcox.test(p, o, alternative = "two.sided", exact = exact, correct = a$continuity == "yes")
  z$p_value <- test$p.value; z$status <- "tested"
 }
 z$rank_biserial <- rank_effect(p, o)
 boot <- replicate(B, rank_effect(sample(p, length(p), replace = TRUE), sample(o, length(o), replace = TRUE)))
 ci <- quantile(boot, c(0.025, 0.975), names = FALSE, type = 7)
 z$effect_CI_low <- ci[1]; z$effect_CI_high <- ci[2]; z
})
res <- do.call(rbind, result)
res$q_value <- NA_real_; ok <- is.finite(res$p_value)
res$q_value[ok] <- p.adjust(res$p_value[ok], method = "BH")
res$selected <- !is.na(res$q_value) & res$q_value < 0.05; res$workflow <- "SDPS"
d <- out_dir(a$out)
audit(d, a, c(a$raw, a$input, a$samples), c(
 "BH family: all finite P values among peptides passing <=50% raw missingness; no FC significance gate.",
 "Normalized input: existing vendor linear abundances; no imputation or transformation by this script.",
 paste("Bootstrap:", B, "independent within-group resamples; percentile 95% CI; seed", seed),
 "Rank-biserial effect and percentile bootstrap algorithm are explicit reconstruction choices; original effect-size definition/bootstrap settings unconfirmed.",
 "Missing samples omitted per peptide with n reported; all-equal peptide assigned P=1. Zero denominator mean -> FC missing.",
 "435 is a manuscript-reported count, not an enforced output target."))
write_table(qc, file.path(d, "missingness_qc.csv"))
write_table(res, file.path(d, "sdps_all_results.csv"))
write_table(res[res$selected, , drop = FALSE], file.path(d, "sdps_deps.csv"))
message(sum(res$selected), " DEPs; ", sum(ok), " tested peptides; no expected result count enforced.")
