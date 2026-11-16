# Evaluated outputs and consensus tool families

Author clarification, 2026-11-17:

> Although seven evaluated outputs exceeded AUROC 0.90, five tool families were
> used for final consensus prediction; CAMP_Bert, CAMP_LSTM, and CAMP_Attention
> are three model variants belonging to the same CAMP framework/tool family.

| Evaluated outputs exceeding AUROC 0.90 in the available score table | Consensus family |
|---|---|
| AI4AMP | AI4AMP |
| AMPScanner | AMPScanner |
| macrel | macrel |
| amplify | AMPlify |
| CAMP_Bert; CAMP_LSTM; CAMP_Attention | CAMP (one family vote) |

Thus seven eligible outputs map to five families. This distinction resolves the
7-versus-5 interpretation; it does not make the models statistically independent,
nor resolve score-table/summary-table numerical discrepancies.
The available table contains ten evaluated outputs in total. If a manuscript
reports eight “models”, specify whether that count refers to outputs or families
and document the complete taxonomy; this file only confirms the three CAMP
variants' shared family as clarified by the author. CAMP_R4 is kept separate in
the supplied map, consistent with the clarified three-variant grouping; its
relationship to the broader CAMP naming still needs explicit author documentation.

`config/model_families.csv` is a visible, editable mapping for observed spellings.
It maps family identity only: it does not assert that variant scores or outputs
are interchangeable, and does not merge or average their AUROC values.

Consensus requires `--family-map` and counts at most one positive vote per family.
Two supported input routes:

1. One verified output per family (including an existing family-level CAMP output).
   Supply that output as its own model column/rows; no invented aggregation is used.
2. Multiple variant outputs within a family. Supply an explicit `--family-rules`
   CSV with columns `family,rule`, using `all`, `any`, or `majority` only if it
   matches the author's actual method. Strict majority requires more than half;
   ties are negative. No rule is the default.

The current source labels a prediction set `LSTM_BERT+Attention`; how this output
was generated remains pending confirmation. The code refuses to collapse three
CAMP variants without a supplied rule. The manuscript must ultimately state the
actual CAMP vote construction and the minimum number of positive families.

Suggested wording (does not add an unverified aggregation rule):

Seven evaluated outputs achieved AUROC > 0.90 and represented five tool families:
AI4AMP, AMPScanner, macrel, AMPlify, and CAMP. CAMP_Bert, CAMP_LSTM, and
CAMP_Attention were treated as variants within the CAMP family rather than as
three independent prediction frameworks. Final consensus prediction was
performed across the five tool families. [Specify how the CAMP-family prediction
was obtained and the required number of positive families.]
