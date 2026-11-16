# Figure contract: count peptide-set overlap, not prediction accuracy or biological activity.
# Adaptation of supplied UpSetR scripts: explicit sets, per-set deduplication, exact intersection audit.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
source(file.path(dirname(script), "..", "shared", "r_utils.R"))
a <- parse_cli(c("input", "out", "sets", "max-intersections"), c("input", "out", "sets"),
 "Rscript upset_plot.R --input peptide_sets.csv --sets ModelA,ModelB,ModelC --out results/upset [--max-intersections 30]\nCSV: one set per named column; each nonblank cell is a peptide identifier/sequence. NOT a binary prediction matrix. Repeated peptides within a set count once.")
for (pkg in c("UpSetR", "svglite")) if (!requireNamespace(pkg, quietly = TRUE)) stop("Install R package: ", pkg)
set.seed(20261116)
x <- read_csv_strict(a$input); sets <- strsplit(a$sets, ",", fixed = TRUE)[[1]]
if (length(sets) < 2 || anyDuplicated(sets) || any(!nzchar(sets))) stop("Provide at least two distinct set columns")
need_cols(x, sets)
nmax <- if (is.null(a[["max-intersections"]])) 30L else numeric_strict(a[["max-intersections"]], "max-intersections", allow_na = FALSE)
if (nmax < 1 || nmax != as.integer(nmax)) stop("max-intersections must be a positive integer")
lst <- lapply(x[, sets, drop = FALSE], function(v) unique(trimws(v[!is.na(v) & nzchar(trimws(v))])))
if (!length(unique(unlist(lst)))) stop("All sets are empty")
universe <- sort(unique(unlist(lst)))
membership <- do.call(cbind, lapply(lst, function(v) as.integer(universe %in% v)))
colnames(membership) <- sets
combinations <- apply(membership, 1, function(v) paste(sets[v == 1], collapse = ";"))
counts <- as.data.frame(table(combinations), stringsAsFactors = FALSE)
names(counts) <- c("exclusive_combination", "intersection_size")
counts <- counts[order(-counts$intersection_size, counts$exclusive_combination), , drop = FALSE]
membership_export <- data.frame(peptide = universe, membership, check.names = FALSE)
size <- data.frame(set = sets, n_unique = lengths(lst),
 n_input_nonmissing = vapply(x[, sets, drop = FALSE], function(v) sum(!is.na(v) & nzchar(trimws(v))), integer(1)))
size$n_removed_within_set_duplicates <- size$n_input_nonmissing - size$n_unique
d <- out_dir(a$out)
audit(d, a, a$input, c("Peptide strings are trimmed, NA/blanks removed, and duplicates removed WITHIN each set. Input case is preserved.",
 "UpSet bars show exclusive intersections. Full intersection table is exported even if figure shows top 30.",
 "Counts refer to peptide membership, not independent biological samples or validated antimicrobial activity."))
write_table(membership_export, file.path(d, "upset_membership.csv"))
write_table(counts, file.path(d, "exclusive_intersections.csv"))
write_table(size, file.path(d, "set_qc.csv"))
draw <- function() {
 print(UpSetR::upset(UpSetR::fromList(lst), sets = sets, nsets = length(sets),
  nintersects = nmax, keep.order = TRUE, order.by = "freq", main.bar.color = "#333333",
  sets.bar.color = "#0072B2", matrix.color = "#333333", shade.color = "#EAEAEA",
  mainbar.y.label = "Exclusive intersection size (peptides)", sets.x.label = "Set size (peptides)",
  text.scale = 1.1, point.size = 2.5, line.size = 0.6))
}
prefix <- file.path(d, "upset")
grDevices::cairo_pdf(paste0(prefix, ".pdf"), width = 8, height = 5); draw(); grDevices::dev.off()
svglite::svglite(paste0(prefix, ".svg"), width = 8, height = 5); draw(); grDevices::dev.off()
grDevices::tiff(paste0(prefix, ".tiff"), width = 8, height = 5, units = "in", res = 600, compression = "lzw"); draw(); grDevices::dev.off()
