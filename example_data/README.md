# Input contracts and synthetic examples

Every `demo_*` file is deliberately synthetic, for software tests only. It contains
no patient data, no study results, and no verified AMP activity. The 3+3 demo
sample manifest is NOT the study cohort; the study has 24 PJI and 17 OA samples.

## Peptidomics

CSV matrix: `peptide_id,sequence,<sample columns>`. Use a stable, unique source-row
identifier for peptide_id when peptide sequences repeat across protein mappings.
Do not silently collapse duplicate sequences. Values are nonnegative LINEAR
abundances; empty cell or NA means missing, zero means observed zero.
CSV manifest: `sample_id,group` with group exactly PJI or OA. Explicitly map the
vendor's Normal columns to OA using this manifest; never infer clinical groups.
The raw and vendor-normalized matrices must share all peptide IDs and sample IDs.
SDPS uses RAW missingness, then VENDOR-NORMALIZED observations for testing. Do not
replace unavailable raw data with normalized data in the actual study.
The demo reuses one matrix for both solely to test software behavior.

EBS input: `peptide_id,sequence,p_value,FC` from the vendor's reported two-sided
Student t tests, with FC=PJI/OA. This is an explicit adapter schema, not an assertion
that source workbook columns already have these exact names.

## Benchmarking

Long CSV: `peptide_id,model,y_true,y_pred,score`. True AMP=1, non-AMP=0.
`y_pred` must be original class calls under documented model-specific thresholds.
For score-only ROC analysis, y_pred may be omitted. Scores must be continuous,
finite, complete, and oriented so larger means AMP. They need not be calibrated
probabilities. Never fill scores using y_pred or confusion-matrix totals.
Models must use the same peptide IDs and truth labels for direct comparison.
Counts-only schema: `model,TP,TN,FP,FN` (cannot recover AUROC).

## Consensus and candidate review

Long CSV: `peptide_id,sequence,model,prediction` with explicit binary prediction for
every selected model and peptide. Missing model outputs are NOT negative predictions.
The code rejects duplicate peptide/model rows instead of inflating vote counts.
Model names are exact and case-sensitive; no unverified aliases are silently merged.
The --models list, --family-map CSV (model,family) and --min-votes threshold must
be stated for each run. Votes count TOOL FAMILIES, not individual variants.
CAMP_Bert, CAMP_LSTM and CAMP_Attention share one CAMP-family vote. If multiple
variants are supplied, --family-rules (family,rule) must explicitly specify
all/any/majority; the historical aggregation rule is not assumed. Alternatively,
supply the verified existing family-level output alone. See docs/MODEL_FAMILIES.md.
This script aggregates existing predictions; it does not download/train/run AMP models.

Manual review CSV: `peptide_id,decision,reason`; decision is include/exclude/pending.
Unreviewed peptides remain pending. Consensus alone does not establish activity,
safety, synthesis feasibility, or selection for animal experiments.

Original workbook adapters require an explicit sheet and column mapping. Clinical
workbooks and model source tables stay outside this repository until the author
has reviewed what may be publicly released.
