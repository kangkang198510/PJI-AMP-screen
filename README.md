# PJI-AMP-analysis

PJI 滑液多肽组、AMP 模型评价与共识筛选代码整理版。

**状态：有来源依据的重构代码，尚非完整历史分析复现或正式公开发行版。**
研究队列为 24 例 PJI、17 例 OA；稿件报告 435 条 DEPs。代码不会把这些
结果数量硬编码为筛选目标。示例文件全部为软件测试用合成数据。

## 代码与方法来源

本版根据作者提供的绘图/交集脚本、差异分析幻灯片，以及当前 ChatGPT
项目本地可用的稿件文本和补充数据整理。没有访问作者全部 ChatGPT 历史对话。
已存在的脚本经过路径参数化、输入校验和错误处理；缺失脚本依据明确方法重构。
这里的“重构”不表示找回了原始执行脚本或模型训练代码。

| 模块 | 当前实现与边界 |
|---|---|
| `01_peptidomics_preprocessing/missingness_filter.R` | 原始矩阵缺失率 >50% 排除；零不自动视为缺失 |
| `01_peptidomics_preprocessing/normalization.R` | 检查并导出已有 vendor-normalized 数据；归一化公式未核实，不自行再归一化 |
| `02_differential_analysis/EBS_screening.R` | 筛选已有 Student t 检验统计表：P<0.05 且 FC≥5 或≤0.2；不重新猜测检验尺度 |
| `02_differential_analysis/SDPS_analysis.R` | 原始缺失率筛选后，对供应商归一化值做双侧 Mann–Whitney U、BH q<0.05；无 FC 门槛 |
| `03_model_benchmarking/calculate_metrics.py` | TP/TN/FP/FN、sensitivity/specificity/accuracy/precision/F1/MCC |
| `03_model_benchmarking/roc_auc_analysis.py` | 从连续评分计算经验 ROC/AUROC；不从硬分类或混淆矩阵生成 AUC |
| `03_model_benchmarking/benchmark_models.py` | 评价已导出的预测；不是模型训练或批量推理程序 |
| `04_consensus_screening/consensus_prediction.py` | 按工具家族计票；CAMP 三个变体不会贡献三票；家族内规则须显式提供 |
| `04_consensus_screening/candidate_selection.py` | 合并作者审阅决定；无审阅则 pending，不自动断言安全性或可合成性 |
| `05_figures/volcano_plot.R` | EBS 与 SDPS 各自使用正确阈值，不混用 P 与 q |
| `05_figures/ROC_plot.py` | 连续评分 ROC，颜色及线型区分；导出 PDF/SVG/600 dpi TIFF |
| `05_figures/upset_plot.R` | 集合内去重，输出精确互斥交集和输入 QC |

辅助目录 `shared/` 为校验函数，`tests/` 为回归测试，`docs/` 为方法边界与检查说明。

## 环境

Python 3.10+；核心评价和共识代码只需标准库。绘图安装：

```bash
python -m pip install -r requirements.txt
```

R 依赖见 `requirements-R.txt`。本次实际测试环境另见 `docs/VALIDATION.md`；
稿件所述 R 4.4.2 是历史报告信息，不代表本次测试运行版本。

## 输入与快速验证

从仓库根目录执行。输入 CSV 约定、缺失值语义和模型覆盖要求见
[example_data/README.md](example_data/README.md)。所有路径均为参数；不依赖桌面文件路径。
请使用新的输出目录，保留每次运行日志。

```bash
python -m unittest discover -s tests -p "test_*.py" -v
Rscript tests/test_r_pipeline.R

python 03_model_benchmarking/benchmark_models.py --input example_data/demo_benchmark.csv --output-dir results/demo_benchmark --score-direction higher-is-AMP
python 05_figures/ROC_plot.py --input example_data/demo_benchmark.csv --output-prefix results/demo_roc --score-direction higher-is-AMP
python 04_consensus_screening/consensus_prediction.py --input example_data/demo_predictions.csv --models Demo_A Demo_B --family-map example_data/demo_family_map.csv --min-votes 2 --output results/demo_consensus.csv
python 04_consensus_screening/candidate_selection.py --input results/demo_consensus.csv --review example_data/demo_review.csv --output results/demo_candidates.csv

Rscript 01_peptidomics_preprocessing/missingness_filter.R --input example_data/demo_matrix.csv --samples example_data/demo_samples.csv --out results/demo_missing --demo
Rscript 01_peptidomics_preprocessing/normalization.R --input example_data/demo_matrix.csv --samples example_data/demo_samples.csv --out results/demo_normalized --mode vendor-normalized --demo
Rscript 02_differential_analysis/SDPS_analysis.R --raw example_data/demo_matrix.csv --input example_data/demo_matrix.csv --samples example_data/demo_samples.csv --out results/demo_sdps --method auto --continuity yes --bootstrap 100 --seed 20261116 --demo
Rscript 02_differential_analysis/EBS_screening.R --input example_data/demo_vendor_statistics.csv --out results/demo_ebs --confirm-test two-sided-student
Rscript 05_figures/volcano_plot.R --input results/demo_sdps/sdps_all_results.csv --workflow SDPS --out results/demo_volcano
Rscript 05_figures/upset_plot.R --input example_data/demo_sets.csv --sets Demo_A,Demo_B,Demo_C --out results/demo_upset
```

`--demo` 仅用于明确的合成测试，不用于绕过真实研究的样本量检查。SDPS 示例
的 bootstrap=100 是快速软件验证参数，不是原研究参数或正式 CI 分析推荐。
`--method auto` 使用无 ties 且较小组≤8时精确检验，否则使用渐近检验；历史
软件版本及连续性校正仍需核实。rank-biserial 和组内百分位 bootstrap 为本版
明确实现的效应量算法，不能声称已复现历史效应量/CI。

## 必须解决的复现差异

现有数据支持原始缺失率筛选保留 1,050 条；在供应商归一化值上重跑 SDPS，
435 条序列集合与结果表一致，P/q/FC 差异仅在浮点精度范围。
作者已确认：7 个 AUROC>0.90 的输出对应 5 个工具家族，CAMP_Bert、
CAMP_LSTM、CAMP_Attention 属于同一 CAMP 家族。详见
[docs/MODEL_FAMILIES.md](docs/MODEL_FAMILIES.md)。
模型评分表与汇总表仍存在 AUROC 数值和样本总数差异。详见
[docs/REPRODUCIBILITY_GAPS.md](docs/REPRODUCIBILITY_GAPS.md)。
这些差异不会由代码自动“修正”或用目标结果倒推参数。

## 公开发布前

1. 核对模型身份、真实预测阈值、训练/验证数据来源及候选审阅规则。
2. 确认供应商归一化方法、缺失处理和历史统计实现。
3. 选择软件许可证并确认改编代码的权利；当前 LICENSE 不授予开源使用许可。
4. 确认作者/版本/引用信息后再创建 CITATION.cff、GitHub Release 和 Zenodo DOI。

原始材料已另行备份，不纳入本仓库；临床数据、账户凭据、Word 和整本补充
工作簿不在发布包内。本任务未上传 GitHub、未发布 Zenodo、未修改论文或引文域。
