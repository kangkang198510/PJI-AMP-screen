# Shared validation and audit helpers. Source from an entry-point script.
parse_cli <- function(allowed, required = character(), usage) {
  argv <- commandArgs(trailingOnly = TRUE)
  if (!length(argv) || "--help" %in% argv) { cat(usage, "\n"); quit(status = 0) }
  out <- list(); i <- 1L
  while (i <= length(argv)) {
    key <- sub("^--", "", argv[i])
    if (!startsWith(argv[i], "--") || !key %in% allowed) stop("Unknown argument: ", argv[i])
    if (!is.null(out[[key]])) stop("Duplicate argument: ", key)
    if (key == "demo") { out[[key]] <- TRUE; i <- i + 1L; next }
    if (i == length(argv) || startsWith(argv[i + 1L], "--")) stop("Missing value: ", key)
    out[[key]] <- argv[i + 1L]; i <- i + 2L
  }
  absent <- setdiff(required, names(out))
  if (length(absent)) stop("Missing options: ", paste(absent, collapse = ", "))
  out
}
read_csv_strict <- function(path) {
  if (!file.exists(path)) stop("Input file not found: ", path)
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
                colClasses = "character", na.strings = c("", "NA"), fileEncoding = "UTF-8-BOM")
  if (anyDuplicated(names(x))) stop("Duplicate column names: ", path)
  if (!nrow(x)) stop("Empty input table: ", path)
  x
}
need_cols <- function(x, cols) {
  absent <- setdiff(cols, names(x))
  if (length(absent)) stop("Missing columns: ", paste(absent, collapse = ", "))
}
numeric_strict <- function(x, label, nonnegative = FALSE, allow_na = TRUE) {
  y <- suppressWarnings(as.numeric(x))
  if (any(is.na(y) & !is.na(x))) stop("Nonnumeric value in ", label)
  if (any(!is.finite(y) & !is.na(y))) stop("Nonfinite value in ", label)
  if (!allow_na && anyNA(y)) stop("Missing value in ", label)
  if (nonnegative && any(y < 0, na.rm = TRUE)) stop("Negative value in ", label)
  y
}
read_samples <- function(path, demo = FALSE) {
  s <- read_csv_strict(path); need_cols(s, c("sample_id", "group"))
  if (anyNA(s$sample_id) || anyDuplicated(s$sample_id)) stop("Sample IDs must be present and unique")
  if (anyNA(s$group) || !all(s$group %in% c("PJI", "OA"))) stop("Groups must be PJI or OA")
  counts <- table(factor(s$group, levels = c("PJI", "OA")))
  if (any(counts < 2)) stop("At least two independent samples per group required")
  if (!demo && !identical(as.integer(counts), c(24L, 17L))) stop("Study requires 24 PJI and 17 OA samples; --demo is for synthetic examples only")
  s
}
read_matrix <- function(path, samples) {
  x <- read_csv_strict(path); need_cols(x, c("peptide_id", "sequence", samples$sample_id))
  if (anyNA(x$peptide_id) || anyDuplicated(x$peptide_id)) stop("peptide_id must be present and unique; do not deduplicate sequences silently")
  if (anyNA(x$sequence)) stop("Missing peptide sequence")
  for (nm in samples$sample_id) x[[nm]] <- numeric_strict(x[[nm]], nm, nonnegative = TRUE)
  x
}
out_dir <- function(path) {
  if (file.exists(path) && !dir.exists(path)) stop("Output path is a file")
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(path)) stop("Cannot create output directory")
  path
}
write_table <- function(x, path) {
  if (file.exists(path)) stop("Refusing to overwrite: ", path)
  write.csv(x, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
}
audit <- function(folder, args, files, notes = character()) {
  path <- file.path(folder, "run_log.txt")
  if (file.exists(path)) stop("Output folder already has run_log.txt; use a fresh output folder")
  hashes <- tools::md5sum(files)
  lines <- c("New implementation reconstructed from supplied materials; not an original historical executable.",
             paste("Time UTC:", format(Sys.time(), tz = "UTC")),
             "Arguments:", capture.output(str(args)), "Notes:", notes,
             "Input MD5 (integrity only):", paste(basename(names(hashes)), unname(hashes)),
             "Runtime/package versions:", capture.output(sessionInfo()))
  writeLines(lines, path, useBytes = TRUE)
}
save_plot <- function(p, prefix, width = 7.2, height = 4.8) {
  for (ext in c("pdf", "svg", "tiff")) if (file.exists(paste0(prefix, ".", ext))) stop("Output already exists: ", prefix)
  ggplot2::ggsave(paste0(prefix, ".pdf"), p, width = width, height = height, device = grDevices::cairo_pdf)
  ggplot2::ggsave(paste0(prefix, ".svg"), p, width = width, height = height, device = svglite::svglite)
  ggplot2::ggsave(paste0(prefix, ".tiff"), p, width = width, height = height, dpi = 600,
                 device = "tiff", compression = "lzw")
}
