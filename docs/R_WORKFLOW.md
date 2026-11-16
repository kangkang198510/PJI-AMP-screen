# R workflow and validation

Run these commands from the repository root. Install R and the packages listed in
`requirements-R.txt`; `Rscript` must be on PATH (or use its full executable path).
The source manuscript reports R 4.4.2. The reconstruction was tested in R 4.6.1;
these versions are not represented as interchangeable historical environments.

```r
install.packages(c("ggplot2", "svglite", "UpSetR"))
```

Every pipeline command writes a `run_log.txt` containing arguments, input-file MD5
checksums, statistical choices and `sessionInfo()`. Use a fresh output directory
for each run: existing output files are never deliberately overwritten.

## Input contract

- Matrix CSV: `peptide_id,sequence,<sample columns>`. `peptide_id` is unique and
  stable; repeated sequences are allowed if they represent distinct source rows.
- Manifest CSV: `sample_id,group`, with group exactly `PJI` or `OA`.
- Actual study inputs require 24 PJI and 17 OA samples. `--demo` explicitly
  permits synthetic smaller cohorts, with at least two samples per group.
- Empty cells and `NA` represent missingness. Numeric zero is an observed value.
  Negative, infinite and malformed abundance values are rejected.
- Abundance inputs must be on the nonnegative **linear** scale. The scripts do
  not silently exponentiate log-scale data, impute missing values, infer groups,
  collapse sequence duplicates or select sample columns by position.
- Raw and normalized SDPS matrices must have the same peptide IDs and manifest
  samples. Raw data determine the missingness gate; normalized data supply the
  observations for testing. A script cannot infer whether a vendor already
  imputed the normalized table; retain the vendor processing record.

## Synthetic end-to-end example

The included demo matrix has three synthetic samples per group. Reusing this
single table as raw and normalized is only a software demonstration, not a
normalization procedure or a representation of the study cohort.

```sh
Rscript 01_peptidomics_preprocessing/missingness_filter.R --input example_data/demo_matrix.csv --samples example_data/demo_samples.csv --out results/demo_missing --demo
Rscript 01_peptidomics_preprocessing/normalization.R --input example_data/demo_matrix.csv --samples example_data/demo_samples.csv --mode vendor-normalized --out results/demo_norm --demo
Rscript 02_differential_analysis/EBS_screening.R --input example_data/demo_vendor_statistics.csv --out results/demo_ebs --confirm-test two-sided-student
Rscript 02_differential_analysis/SDPS_analysis.R --raw example_data/demo_matrix.csv --input example_data/demo_matrix.csv --samples example_data/demo_samples.csv --out results/demo_sdps --method auto --continuity yes --bootstrap 2000 --seed 20261116 --demo
Rscript 05_figures/volcano_plot.R --input results/demo_sdps/sdps_all_results.csv --out results/demo_volcano --workflow SDPS
Rscript 05_figures/upset_plot.R --input example_data/demo_sets.csv --sets Demo_A,Demo_B,Demo_C --out results/demo_upset
Rscript tests/test_r_pipeline.R
```

`--help` is available on each entry point. A real-data invocation omits `--demo`
and supplies separately exported raw and vendor-normalized matrices.

## What is implemented

`missingness_filter.R` retains peptides missing in no more than 50% of the raw
sample columns. It exports every peptide's missingness and per-sample QC.

`normalization.R` validates and re-exports an **already vendor-normalized** matrix.
It performs no normalization, imputation or log transform. The actual vendor
normalization algorithm is not recoverable from the supplied code; the filename
does not imply that the missing algorithm has been reconstructed. Metabolomics
internal-standard normalization and half-minimum imputation must not be imported
into this peptidomics pipeline without source evidence.

`EBS_screening.R` screens existing, author-confirmed vendor two-sided Student
t-test statistics: unadjusted `p_value < 0.05` and `FC >= 5` or `FC <= 0.2`, where
FC means PJI/OA. The explicit `--confirm-test two-sided-student` flag records the
input contract, not an independent validation of the vendor's method. This
script does not recompute t tests on an unknown processing scale. EBS is
exploratory and is independent of SDPS; do not feed an SDPS-filtered list to EBS.

`SDPS_analysis.R` first applies raw missingness, then tests each retained peptide
using two-sided Mann–Whitney U statistics and nonmissing normalized observations.
At least two observations per group are required. BH correction uses all finite
P values in the retained test family, and `q_value < 0.05` defines selection.
There is **no additional FC threshold**. Missing/untestable results remain
explicit, and the number 435 is never used as a target or gate.

`--method auto` uses an exact test when at least one group has at most eight
observations and there are no pooled ties; otherwise it uses the asymptotic
test. `--method asymptotic` always uses the asymptotic test. The continuity
correction is explicit (`--continuity yes|no`). These options must be matched to
the original analysis; the old supplied Python code relied on the SciPy default
method without recording the package version. All-equal observations receive
P=1. Arithmetic means of observed normalized abundances determine FC; zero OA
mean gives undefined FC, which remains missing instead of silently adding a
pseudocount.

The effect is **rank-biserial correlation**,
`2 * [P(PJI > OA) + 0.5 P(PJI = OA)] - 1`, oriented positive for larger PJI values.
The 95% interval uses independent within-group nonparametric bootstrap samples
and percentile quantiles. The seed, number of resamples (default 2000), effect
definition and interval algorithm are explicit reconstruction choices. The
original manuscript's precise effect/CI implementation remains to be confirmed;
these intervals must not be described as the original historical calculations.

## Figures

Volcano plots use the full result table: EBS shows raw P with FC boundaries,
whereas SDPS shows BH q with no FC significance boundaries. Shape and color
encode categories; undefined log2FC/P/q values are excluded from plotting and
counted in source-data QC. Exact zero P/q is floored only for log rendering.

UpSet expects **columns of peptide identifiers**, not a 0/1 membership matrix.
It trims blank padding, removes missing cells, and deduplicates within each set.
Bars represent exclusive intersections. Full membership, intersection counts
and deduplication QC are exported even when only the largest 30 intersections
are shown. It does not equate overlap with prediction accuracy or activity.

Both scripts export PDF, SVG and 600 dpi TIFF. Actual sample size belongs in the
associated legend/source table; the number of plotted peptide features is not
the number of biological replicates. UpSetR emits ggplot2 deprecation warnings
with the tested packages; these warnings did not prevent rendering.

## Completed checks and limits

- R parsing and synthetic regression checks passed: the 50% boundary, study
  group counts, unchanged normalized import, a reference Mann–Whitney P value,
  tie handling, rank-biserial direction, zero denominator handling, BH, seeded
  bootstrap repeatability, malformed-data rejection and EBS threshold edges.
- All six R entry points ran on the included synthetic examples. R-generated
  PNG previews of the demo volcano/UpSet figures were visually inspected; PDF,
  SVG and 600 dpi TIFF exports were generated.
- A private validation on the 1,050 retained study rows yielded 435 SDPS DEPs.
  After the audited sequence-ID suffix mapping, the selected sequence set
  matched the supplied reported 435 sequences. This verifies the retained
  matrix and the tested settings, not the unknown normalization algorithm.
- The study validation used 100 bootstrap resamples to test execution; those
  CI estimates are validation outputs, not final manuscript intervals. The
  complete raw-matrix filtering audit and sequence-ID mapping are documented
  separately. Private matrices and validation result tables are not bundled
  in the public-ready example directory.

For Windows installations with an inherited invalid `C.UTF-8` locale, use a
UTF-8 locale supported by Windows. If temporary validation libraries are used,
include both that directory and the existing user library in `R_LIBS_USER`;
replacing the user-library path can hide installed ggplot2 dependencies.
