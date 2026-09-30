# OTOMO: Outlier-Conscious Transfer-Oriented Minority Oversampling

**OTOMO (Outlier-Conscious Transfer-Oriented Minority Oversampling)** is an outlier-aware and target-oriented oversampling method developed for **Cross-Project Defect Prediction (CPDP)**.

OTOMO generates informative synthetic defective instances by jointly considering the **source minority structure**, **target-domain relevance**, and **outlier information**.

Unlike conventional oversampling methods that mainly rely on the distribution of the source minority class, OTOMO incorporates information from the **unlabeled target project** to determine which source minority instances should contribute more strongly to synthetic sample generation.

---

## Key Features

OTOMO is designed around three main principles:

### 1. Source Minority Structure

OTOMO identifies minority-class instances in the labeled source project and analyzes their relationships with neighboring minority instances. This allows synthetic instances to be generated while considering the local structure of the defective class.

### 2. Target-Oriented Relevance

Information from the unlabeled target project is incorporated into the oversampling process. Source minority instances that are more relevant to the target-domain distribution receive greater importance during synthetic sample generation.

### 3. Outlier-Aware Sampling

OTOMO identifies potential outliers in the source minority class and incorporates this information into the sampling process to reduce the influence of atypical instances during synthetic sample generation.


## Contents

- [Requirements](#requirements)
- [Repository structure](#repository-structure)
- [Quick start](#quick-start)
- [Full reproduction](#full-reproduction)
- [Understanding the outputs](#understanding-the-outputs)
- [Reproducibility notes](#reproducibility-notes)
- [Using OTOMO on your own data](#using-otomo-on-your-own-data)

---

## Requirements

- **MATLAB.** All results in this repository were produced and re-checked with **MATLAB R2025a**. The earliest release that can run every script is **R2021a**: `friedman_nemenyi_figure.m` uses `tiledlayout(...,'TileSpacing','tight')`, which was added in R2021a (other recent features used are `readtable(...,'VariableNamingRule',...)` from R2020b and `exportgraphics` from R2020a). Releases before R2025a have not been tested.
- **Statistics and Machine Learning Toolbox.** No other toolbox is needed; to check that it is installed, run `ver('stats')` in MATLAB (it should list `Statistics and Machine Learning Toolbox`). Microsoft Excel is not needed: MATLAB reads the `.xlsx` data files itself.
- Disk space: the clone is under 20 MB, and the scripts write only a few MB.

Toolbox functions used by the code:

| Function | Used for | Used in |
|---|---|---|
| `TreeBagger`, `predict` | Random Forest classifier | `main.m`, `sensitivity_figure.m` |
| `fitclinear`, `predict` | Logistic Regression classifier | `main.m`, `sensitivity_figure.m` |
| `perfcurve` | ROC-AUC | `main.m`, `sensitivity_figure.m` |
| `kmeans` | Clustering of the target project | `otomo.m`, `tomo.m`, `target_clustering_summary.m`, `otomo_clustering.m` |
| `pdist2` | Distances between instances | `otomo.m`, `tomo.m`, `smote_oversample.m`, `adasyn_oversample.m`, `dpp_oversample.m` |
| `pdist` | TCA+ normalization choice | `normal.m` |
| `chi2inv` | Mahalanobis outlier threshold | `otomo.m`, `outlier_plot.m` |
| `mahal` | Mahalanobis distances | `outlier_plot.m` |
| `tiedrank`, `chi2cdf`, `finv`, `fcdf`, `normpdf`, `normcdf` | Friedman, Iman–Davenport and Nemenyi tests | `friedman_nemenyi_figure.m` |

The optional Wilcoxon check in [`docs/SIGNIFICANCE_TESTS.md`](docs/SIGNIFICANCE_TESTS.md) also uses `signrank`. Everything else (`readmatrix`, `readtable`, `writecell`, `writetable`, `exportgraphics`, `tiledlayout`, `savefig`, `normalize`, `integral`, `fzero`, ...) is base MATLAB.

---

## Repository structure

```text
OTOMO/
├── main.m                          # Experiment driver: 7 methods x 66 pairs x 30 seeds, RF and LR -> Results/Metric_Summaries/
├── otomo.m                         # OTOMO oversampling (the proposed method)
├── coral.m                         # CORAL alignment: the CORAL baseline, and applied after every oversampler
├── tca_plus.m                      # TCA+ baseline (Nam et al., 2013)
├── tca.m                           # Transfer Component Analysis, called by tca_plus.m
├── normal.m                        # TCA+ normalization rules, called by tca_plus.m
├── smote_oversample.m              # SMOTE baseline
├── adasyn_oversample.m             # ADASYN baseline
├── dpp_oversample.m                # DPP baseline
├── tomo.m                          # TOMO baseline
├── friedman_nemenyi_figure.m       # Friedman / Iman-Davenport / Nemenyi tests and CD diagrams
├── sensitivity_figure.m            # Gamma sensitivity experiment (OTOMO, gamma = 0, 0.1, ..., 1)
├── sensitivity_figure_from_csv.m   # Redraws the sensitivity figures from the saved summary CSV
├── plot_sensitivity_figures.m      # Plotting function used by the two sensitivity scripts
├── target_clustering_summary.m     # Agreement of OTOMO's k-means step with the true labels (purity, ARI, NMI)
├── otomo_clustering.m              # Stand-alone copy of that clustering step (no script calls it)
├── outlier_plot.m                  # Mahalanobis outlier plots and counts for every project
├── data/                           # 18 projects in 4 families, one .xlsx file per project
│   ├── AEEEM/                      #   EQ, JDT, LC, ML, PDE
│   ├── NASA/                       #   CM1, MW1, PC1, PC3, PC4
│   ├── RELINK/                     #   A, S, Z
│   └── SOFTLAB/                    #   ar1, ar3, ar4, ar5, ar6
├── Results/
│   ├── Metric_Summaries/           # 8 CSVs (2 classifiers x 4 metrics), written by main.m
│   └── Significance_Tests/
│       └── gamma_001_0.1/          # Summary CSV and CD diagrams, written by friedman_nemenyi_figure.m
├── OTOMO_Gamma_Sensitivity/        # Summary CSV, .mat file and Sensitivity_Figures/, written by sensitivity_figure.m
├── target clustering summary/      # OTOMO_Project_Clustering_Summary.csv, written by target_clustering_summary.m
├── outlier_plots/                  # Per-project plots (.png/.eps/.fig) and 2 CSVs, written by outlier_plot.m
├── Result.xlsx                     # Workbook assembled by hand from the CSVs above (not written by any script)
├── docs/
│   └── SIGNIFICANCE_TESTS.md       # Details of the statistical tests and their current results
└── README.md
```

### Data

Each file `data/<Family>/<Project>.xlsx` has one header row, one column per software metric, and the defect label in the **last column** (`1` = defective, `0` = clean). Rows with missing or non-numeric values are dropped when the file is read (none of the included files has any).

| Family | Projects (instances, defective) | Features | Source → target pairs |
|---|---|---:|---:|
| AEEEM | EQ (324, 129), JDT (997, 206), LC (691, 64), ML (1862, 245), PDE (1497, 209) | 61 | 20 |
| NASA | CM1 (327, 42), MW1 (251, 25), PC1 (696, 55), PC3 (1073, 132), PC4 (1276, 176) | 37 | 20 |
| RELINK | A (194, 98), S (56, 22), Z (399, 118) | 26 | 6 |
| SOFTLAB | ar1 (121, 9), ar3 (63, 8), ar4 (107, 20), ar5 (36, 8), ar6 (101, 15) | 29 | 20 |

Every ordered pair of distinct projects in the same family is one cross-project task (the source is labelled, the target is used as the test set): 66 pairs in total.

---

## Quick start

Clone the repository, start MATLAB, and make the repository root the current folder. **All scripts must be run from the repository root**, because they build every path from `pwd`.

```bash
git clone https://github.com/suravi0827/OTOMO.git
cd OTOMO
```

These instructions describe the `main` branch. The older tag `v1.0.0` has a different layout (`main_otomo.m`, no baselines or analysis scripts) and does not match this README.

```matlab
>> cd path/to/OTOMO     % the folder that contains main.m and data/
```

### 1. Checks that only need the published CSVs (about a minute each)

These scripts do not train any model. They rebuild figures and tables from the CSVs already in the repository, and overwrite the committed copies in `Results/Significance_Tests/gamma_001_0.1/` and `OTOMO_Gamma_Sensitivity/Sensitivity_Figures/` (`git checkout -- Results OTOMO_Gamma_Sensitivity` restores them):

```matlab
>> friedman_nemenyi_figure        % Friedman / Nemenyi tests and CD diagrams, from Results/Metric_Summaries/
>> sensitivity_figure_from_csv    % gamma sensitivity figures, from OTOMO_Gamma_Sensitivity/OTOMO_Gamma_Sensitivity_Summary.csv
```

With R2025a the regenerated `Friedman_Nemenyi_Summary.csv` is byte-identical to the committed one, so `git status` lists only the 20 figure files (10 `.eps`, 10 `.png`) as modified. This is expected: the `.eps` files differ only in their `%%CreationDate` line and the `.png` files only in embedded metadata. Their pixels are identical, except for a one-pixel shift of the legend in `LR_Gamma_MCC.png`.

### 2. Reproduce part of the main experiment (under 10 minutes)

The full experiment takes several hours. For a quick check, run `main.m` on the smallest family, RELINK (3 projects, 6 pairs, all 7 methods, seeds 1–30), and compare with the published numbers. The snippet below writes a copy of `main.m` whose only changes are the `datasets` line and a separate output folder, so the published results are not touched:

```matlab
% Make main_check.m: one family only, output to Results_check/ instead of Results/
fam = 'RELINK';   % family to check: 'RELINK' (about 6 min) or 'SOFTLAB' (about 9 min)
src = fileread('main.m');
src = strrep(src,'datasets = ["NASA","RELINK","SOFTLAB","AEEEM"];',['datasets = "' fam '";']);
src = strrep(src,'outputRoot = fullfile(pwd,''Results'');','outputRoot = fullfile(pwd,''Results_check'');');
assert(contains(src,['datasets = "' fam '";']) && contains(src,'Results_check'),'Config lines not found: edit main_check.m by hand.');
fid = fopen('main_check.m','w'); fwrite(fid,src); fclose(fid);
main_check
```

(You can also copy `main.m` to `main_check.m` and change the two lines by hand.) It prints about 800 lines of progress (OTOMO and TOMO report every run). It is finished when `ALL EXPERIMENTS COMPLETE` appears, and it should print no warnings. Then compare the family block of each new CSV (header, family row, pair rows, `W/T/L` and `Average (%)`) with the same block of the published CSV. The `Overall W/T/L` and `Overall Average (%)` rows of the new files cover this one family only, so they are not compared:

```matlab
files = dir(fullfile('Results_check','Metric_Summaries','*.csv'));
assert(numel(files) == 8,'Expected 8 CSVs in Results_check/Metric_Summaries; did main_check finish?');
for f = files'
    new = splitlines(strtrim(string(fileread(fullfile(f.folder,f.name)))));
    old = splitlines(strtrim(string(fileread(fullfile('Results','Metric_Summaries',f.name)))));
    i = find(old == new(2),1);                                 % first row of the family block, e.g. "RELINK,,,,,,,,"
    j = i + find(startsWith(old(i+1:end),"Average (%)"),1);    % the block's "Average (%)" row
    a = new(1:j-i+2);  b = [old(1); old(i:j)];
    same = isequal(a,b);
    fprintf('%-62s identical: %d\n',f.name,same);
    if ~same, fprintf('    yours:     %s\n    published: %s\n',[a(a ~= b), b(a ~= b)].'); end
end
```

On MATLAB R2025a all eight files print `identical: 1`: every `mean +/- sd` cell matches to the last printed digit. On other MATLAB versions or platforms, small differences in the last digits are possible (see [Reproducibility notes](#reproducibility-notes)). Clean up afterwards with `delete main_check.m` and `rmdir('Results_check','s')`.

To check SOFTLAB instead, set `fam = 'SOFTLAB';` in the first block (20 pairs instead of 6; about 9–10 minutes, also identical on R2025a). The comparison code needs no change. NASA and AEEEM have larger projects and take much longer.

---

## Full reproduction

Run the steps below in this order, from the repository root. Each script has its settings at the top of the file; the defaults are the settings used for the published results.

> **Every script overwrites its published output folder.** To keep the committed results for comparison, change the output setting named in each step (`outputRoot`, `outputFolder` or `figureFolder`) before running, or work in a separate clone. If you have already overwritten them, `git checkout -- Results OTOMO_Gamma_Sensitivity outlier_plots "target clustering summary"` restores every published output (`git status` shows what changed). Several scripts start with `clear; close all; clc`, so save anything you need from the MATLAB workspace first.

### Step 1: `main.m`, the main experiment

**What it computes.** For every source → target pair in every family, for every method and every seed:

1. Source and target are z-scored separately, each with its own mean and standard deviation.
2. The method is applied to the source:
   - `OTOMO`, `SMOTE`, `ADASYN`, `DPP`, `TOMO`: oversample the defective class of the source until it has as many instances as the clean class (`samplingRatio = 1.0`), then align the balanced source to the target with CORAL (`coral.m`).
   - `CORAL`: CORAL alignment without oversampling.
   - `TCA+`: `tca_plus.m` (TCA+ normalization, then TCA) maps source and target to a shared space.
3. A Random Forest (`TreeBagger`, `nTrees` trees, `floor(sqrt(p))` predictors per split, minimum leaf size 1) and a ridge-regularized Logistic Regression (`fitclinear`, `lbfgs`, penalty `lambdaLR`) are trained on the source and evaluated on the labelled target.
4. ROC-AUC, F1-score, G-mean and MCC are recorded (defective = positive class).

Per pair, the 30 seed values are summarized as `mean +/- sd`, and one CSV per classifier and metric is written.

**Command.**

```matlab
>> main
```

**Settings** (top of `main.m`):

| Variable | Default | Meaning |
|---|---|---|
| `datasets` | `["NASA","RELINK","SOFTLAB","AEEEM"]` | Families to run (sub-folders of `data/`) |
| `gammas` | `0.1` | OTOMO's `gamma` (weight of source-minority distances against target distances). A vector runs every value and writes one file set per value (`gamma_001_<value>`, `gamma_002_<value>`, ...) |
| `seeds` | `1:30` | Random seeds; each run calls `rng(seed,'twister')` |
| `nTrees` | `100` | Number of trees in the Random Forest |
| `lambdaLR` | `1e-4` | Ridge penalty of the Logistic Regression |
| `methods` | `["OTOMO","CORAL","TCA+","ADASYN","DPP","SMOTE","TOMO"]` | Methods to run. In the CSVs, the selected methods always appear in the order OTOMO, CORAL, TCA+, DPP, ADASYN, SMOTE, TOMO; without OTOMO, no `W/T/L` rows are written |
| `samplingRatio` | `1.0` | Target defective/clean ratio after oversampling (1.0 = balanced) |
| `kSMOTE`, `kADASYN` | `7`, `7` | Number of nearest neighbours for SMOTE and ADASYN |
| `tomoLambda` | `0.1` | Trade-off parameter of TOMO |
| `otomoConf` | `0.995` | Chi-square confidence level of OTOMO's Mahalanobis outlier threshold |
| `otomoBeta` | `0.1` | Not used by the code |
| `outputRoot` | `fullfile(pwd,'Results')` | Output folder. **Change it (for example to `fullfile(pwd,'Results_rerun')`) to keep the published results** |

Other values are fixed in the code, for example: DPP uses a kernel width of `0.1` (in `main.m`); OTOMO and TOMO pick each neighbour among the 7 nearest source-minority instances (`otomo.m`, `tomo.m`); and their k-means step uses 20 replicates.

**Outputs.** Eight CSVs in `<outputRoot>/Metric_Summaries/`:

```text
Random_Forest_All_Datasets_ROC_AUC_gamma_001_0.1.csv
Random_Forest_All_Datasets_F1_Score_gamma_001_0.1.csv
Random_Forest_All_Datasets_G_Measure_gamma_001_0.1.csv
Random_Forest_All_Datasets_MCC_Score_gamma_001_0.1.csv
Logistic_Regression_All_Datasets_ROC_AUC_gamma_001_0.1.csv
Logistic_Regression_All_Datasets_F1_Score_gamma_001_0.1.csv
Logistic_Regression_All_Datasets_G_Measure_gamma_001_0.1.csv
Logistic_Regression_All_Datasets_MCC_Score_gamma_001_0.1.csv
```

The files are written only **after all families have finished**; an interrupted run saves nothing. Per-seed values are kept in memory only and are not saved.

**Runtime.** Several hours for all four families (under 10 minutes for RELINK alone). The console output is long, because OTOMO and TOMO print a few lines for every run.

**Checking a full rerun.** With `outputRoot = fullfile(pwd,'Results_rerun')`, compare all eight files line by line:

```matlab
for f = dir(fullfile('Results_rerun','Metric_Summaries','*.csv'))'
    new = splitlines(strtrim(string(fileread(fullfile(f.folder,f.name)))));
    old = splitlines(strtrim(string(fileread(fullfile('Results','Metric_Summaries',f.name)))));
    fprintf('%-62s identical: %d\n',f.name,isequal(new,old));
end
```

### Step 2: `friedman_nemenyi_figure.m`, significance tests

**What it computes.** For each classifier and metric, the per-pair means of the seven methods are ranked on every pair. A Friedman test with the Iman–Davenport F-approximation checks whether the methods differ, and, if they do, the Nemenyi post-hoc test with its critical difference (CD) shows which ones differ. This is done over all 66 pairs and within each family. The method and the current results are described in [`docs/SIGNIFICANCE_TESTS.md`](docs/SIGNIFICANCE_TESTS.md).

**Command.**

```matlab
>> friedman_nemenyi_figure
```

**Settings** (top of the file):

| Variable | Default | Meaning |
|---|---|---|
| `summaryFolder` | `fullfile(pwd,'Results','Metric_Summaries')` | Input CSVs from Step 1. Point it to `Results_rerun/Metric_Summaries` to test your own rerun |
| `outputFolder` | `fullfile(pwd,'Results','Significance_Tests')` | Output folder; a `<gammaTag>` sub-folder is added. **Change it to keep the published results** |
| `gammaTag` | `""` | Which `gamma_*` file set to use; empty means the only one present (an error if there are several) |
| `datasets` | `["AEEEM","RELINK","NASA","SOFTLAB"]` | Families tested separately and pooled into `All` |
| `alpha` | `0.05` | Significance level |
| `tieCorrection` | `false` | `false`: Friedman statistic without tie correction, as in the paper; `true`: tie-corrected statistic |
| `highlightMethod` | `"OTOMO"` | Method drawn in blue and used as the reference in the verdict columns |

**Outputs** in `<outputFolder>/gamma_001_0.1/`:

- `Friedman_Nemenyi_Summary.csv`: 40 rows (4 metrics x 2 classifiers x 5 groups: the four families and `All`).
- `RF_Overall_Friedman_Nemenyi_CD.eps` / `.png` and `LR_Overall_Friedman_Nemenyi_CD.eps` / `.png`: CD diagrams over all 66 pairs, one panel per metric (PNG at 600 dpi).

**Runtime.** Under a minute.

### Step 3: `sensitivity_figure.m`, gamma sensitivity

**What it computes.** Runs the OTOMO pipeline of Step 1 (OTOMO, then CORAL, then RF and LR) for `gamma = 0, 0.1, ..., 1.0` on all 66 pairs with seeds 1–30. For each family and gamma, it averages the per-pair means over the family's pairs and then draws one figure per classifier and metric.

**Command.**

```matlab
>> sensitivity_figure
```

**Settings** (top of the file):

| Variable | Default | Meaning |
|---|---|---|
| `datasets` | `["NASA","RELINK","SOFTLAB","AEEEM"]` | Families to run |
| `gammas` | `0:0.1:1.0` | Gamma values to test |
| `seeds` | `1:30` | Random seeds |
| `nTrees`, `lambdaLR` | `100`, `1e-4` | Classifier settings, as in `main.m` |
| `samplingRatio`, `otomoConf` | `1.0`, `0.995` | OTOMO settings, as in `main.m` |
| `outputRoot` | `fullfile(pwd,'OTOMO_Gamma_Sensitivity')` | Output folder. **Change it to keep the published results** |

**Outputs** in `<outputRoot>/`:

- `OTOMO_Gamma_Sensitivity_Summary.csv`: one row per family and gamma (44 rows).
- `OTOMO_Gamma_Sensitivity.mat`: the variables `datasets`, `gammas`, `RF_ROC`, `RF_G`, `RF_F1`, `RF_MCC`, `LR_ROC`, `LR_G`, `LR_F1`, `LR_MCC` (each families x gammas) and `ValidPairs`.
- `Sensitivity_Figures/<RF|LR>_Gamma_<ROC_AUC|G_Mean|F1|MCC>.eps` and `.png` (8 figures, PNG at 600 dpi).

**Runtime.** Several hours (11 gamma values x 66 pairs x 30 seeds).

**Consistency check.** The `Gamma = 0.1` rows use exactly the setting of Step 1. In the published files, every `Gamma = 0.1` value, multiplied by 100 and rounded to two decimals, equals the OTOMO entry of the matching family's `Average (%)` row in `Results/Metric_Summaries/` (all 32 values: 4 families x 2 classifiers x 4 metrics).

**Figures only: `sensitivity_figure_from_csv.m`.** Redraws the eight figures from the summary CSV in about a minute, without running any experiment:

```matlab
>> sensitivity_figure_from_csv
```

| Variable | Default | Meaning |
|---|---|---|
| `csvFile` | `OTOMO_Gamma_Sensitivity/OTOMO_Gamma_Sensitivity_Summary.csv` | Summary CSV to plot |
| `figureFolder` | `OTOMO_Gamma_Sensitivity/Sensitivity_Figures` | Output folder (overwritten) |
| `datasetOrder` | `["NASA","RELINK","SOFTLAB","AEEEM"]` | Legend and line order |

### Step 4: `target_clustering_summary.m`, quality of the target clustering

**What it computes.** OTOMO splits the unlabelled target project into two clusters with k-means (k = 2, squared Euclidean distance, 20 replicates, `rng(seed,'twister')`) and treats the smaller cluster as the likely defective one. This script applies the same step to each of the 18 projects (z-scored) for seeds 1–30 and measures how well the two clusters agree with the true labels: purity, adjusted Rand index (ARI) and normalized mutual information (NMI).

**Command.**

```matlab
>> target_clustering_summary
```

**Settings:** `datasets` (default `["NASA","RELINK","SOFTLAB","AEEEM"]`), `seeds` (`1:30`) and `outputFolder` (`fullfile(pwd,'target clustering summary')`; change it to keep the published file).

**Output:** `target clustering summary/OTOMO_Project_Clustering_Summary.csv` (the folder name contains spaces).

**Runtime.** Under a minute. With R2025a the regenerated CSV is byte-identical to the committed one.

### Step 5: `outlier_plot.m`, Mahalanobis outliers in each project

**What it computes.** A descriptive analysis of the data, separate from the experiment. For each project (all instances, both classes), it drops constant features, z-scores the rest, computes each instance's squared Mahalanobis distance with `mahal`, and marks as outliers the instances above the chi-square threshold `chi2inv(1 - alpha, p)`, where `p` is the number of features. It saves one plot per project and two summary tables.

**Command.**

```matlab
>> outlier_plot
```

**Settings:** `datasets` (default `["AEEEM","RELINK","SOFTLAB","NASA"]`), `alpha` (`0.005`, that is, the 99.5% quantile, the same level as `otomoConf`) and `outputRoot` (`fullfile(pwd,'outlier_plots')`; change it to keep the published files).

**Outputs** in `outlier_plots/`:

- `<Family>/<Project>_mahalanobis_outliers.png` (300 dpi), `.eps` and `.fig` for each of the 18 projects.
- `mahalanobis_outlier_summary_all_datasets.csv` and `mahalanobis_outlier_summary_paper.csv`.

**Runtime.** About 1.5–2 minutes. MATLAB prints `Matrix is close to singular or badly scaled` warnings from `mahal` for 11 of the 18 projects; they are expected and do not stop the script. With R2025a both regenerated CSVs are byte-identical to the committed ones. `git status` still lists all 54 plot files as modified. This is expected: the `.png` files are pixel-identical, the `.eps` files differ only in their `%%CreationDate` line, and the `.fig` files embed their creation time but contain the same plotted data.

### Step 6: `Result.xlsx`, the combined workbook

`Result.xlsx` is **assembled by hand** from the CSVs above; no script writes it. To rebuild it, copy the tables into a workbook as follows:

| Sheet | Source file | Content |
|---|---|---|
| `RF_AUC`, `RF_GM`, `RF_F1`, `RF_MCC` | `Results/Metric_Summaries/Random_Forest_All_Datasets_<ROC_AUC, G_Measure, F1_Score, MCC_Score>_gamma_001_0.1.csv` | The CSV with source and target merged into one `Source→Target` column (columns A–H), and the same table in percent (columns I–P) |
| `LR_AUC`, `LR_GM`, `LR_F1`, `LR_MCC` | `Results/Metric_Summaries/Logistic_Regression_All_Datasets_<...>_gamma_001_0.1.csv` | As above, for Logistic Regression |
| `F1` | `W/T/L` and `Average (%)` rows of the two `F1_Score` CSVs | Per-family and overall summary, Random Forest then Logistic Regression |
| `MCC` | `W/T/L` and `Average (%)` rows of the two `MCC_Score` CSVs | As above, for MCC |
| `target_clustering` | `target clustering summary/OTOMO_Project_Clustering_Summary.csv` | Purity, ARI and NMI per project |

All 66 pair rows of the eight metric sheets match the CSVs exactly (the SOFTLAB project names are written in upper case there). In the hand-made `F1` and `MCC` summary sheets, 8 averages differ from the CSV `Average (%)` values by 0.01; the CSVs are the reference. The significance tests, sensitivity results and outlier counts are not in the workbook; use their CSVs directly.

---

## Understanding the outputs

### Metric summaries (`Results/Metric_Summaries/*.csv`)

File names follow `<Classifier>_All_Datasets_<Metric>_gamma_<index>_<gamma>.csv`, where `<Classifier>` is `Random_Forest` or `Logistic_Regression`, `<Metric>` is `ROC_AUC`, `F1_Score`, `G_Measure` or `MCC_Score`, and `gamma_001_0.1` means the first (and only) gamma value, 0.1.

```text
Source,Target,OTOMO,CORAL,TCA+,DPP,ADASYN,SMOTE,TOMO
NASA,,,,,,,,
CM1,MW1,0.7554 +/- 0.0103,0.7272 +/- 0.0249,0.7075 +/- 0.0174,...
...
W/T/L,,,16/0/4,14/0/6,9/0/11,8/0/12,5/0/15,19/0/1
Average (%),,72.65,68.94,68.87,71.82,71.26,71.72,61.23
RELINK,,,,,,,,
...
Overall W/T/L,,,36/2/28,48/0/18,...
Overall Average (%),,73.83,72.13,...
```

- **Family blocks** appear in the order NASA, RELINK, SOFTLAB, AEEEM. A row with only the family name starts a block.
- **Pair rows**: `mean +/- sd` of the metric over the 30 seeds (sample standard deviation), with 4 decimals.
- **`W/T/L`**: for each baseline column, the number of pairs in the family where OTOMO's mean is higher / equal / lower than that baseline's mean, after rounding both to 4 decimals. The OTOMO column is empty. All four metrics are higher-is-better.
- **`Average (%)`**: 100 x the mean of the per-pair means over the family's pairs, with 2 decimals.
- **`Overall W/T/L`** and **`Overall Average (%)`**: the same over all 66 pairs pooled.
- If all runs of a pair failed for some method, extra rows `Valid tasks (average)` and (when OTOMO is selected) `Valid comparisons vs OTOMO` report how many pairs have values. They do not appear in the published files.

The metrics treat defective (label 1) as the positive class:

| File tag | Metric |
|---|---|
| `ROC_AUC` | Area under the ROC curve (`perfcurve`) of the classifier's score for class 1 |
| `F1_Score` | F1-score of the defective class |
| `G_Measure` | **G-mean**, the geometric mean `sqrt(TPR * TNR)` of recall (true-positive rate) and specificity (true-negative rate) |
| `MCC_Score` | Matthews correlation coefficient |

F1 and MCC are set to 0 when they are undefined (no predicted or no actual positives), and ROC-AUC is 0.5 when all scores are equal.

### Significance tests (`Results/Significance_Tests/gamma_001_0.1/`)

`Friedman_Nemenyi_Summary.csv` has one row per metric, classifier and group, with the Friedman and Iman–Davenport statistics, the Nemenyi CD, the average rank of each method and a verdict for OTOMO against each baseline. The CD diagrams show the average ranks (rank 1 = best, on the right); methods joined by a bar are not significantly different. See [`docs/SIGNIFICANCE_TESTS.md`](docs/SIGNIFICANCE_TESTS.md) for the column definitions, how to read the diagrams and the current results.

### Gamma sensitivity (`OTOMO_Gamma_Sensitivity/`)

`OTOMO_Gamma_Sensitivity_Summary.csv` has the columns `Dataset, Gamma, Valid_Pairs, RF_ROC_AUC, RF_G_Mean, RF_F1, RF_MCC, LR_ROC_AUC, LR_G_Mean, LR_F1, LR_MCC`. Each metric value is the mean, over the family's pairs, of the per-pair mean over the 30 seeds. Values are fractions (0–1), not percentages. `Valid_Pairs` is the number of pairs with a result for all four RF metrics (20, or 6 for RELINK). The figures plot these values against gamma, one line per family.

### Target clustering (`target clustering summary/OTOMO_Project_Clustering_Summary.csv`)

Columns: `Dataset`, `Project`, `Purity (%) mean +/- SD`, `ARI (%) mean +/- SD` and `NMI (%) mean +/- SD`, each over the 30 seeds, in percent. The SD is often 0.00 because k-means with 20 replicates usually finds the same partition for every seed. Low ARI and NMI mean that the two clusters agree only weakly with the true defective/clean labels.

### Outliers (`outlier_plots/`)

- `mahalanobis_outlier_summary_all_datasets.csv`: `Dataset, Project, Samples, Features` (after dropping constant features), `Defect_Count, Defect_Percentage, MD_Threshold, Outlier_Count, Outlier_Percentage`.
- `mahalanobis_outlier_summary_paper.csv`: the same counts formatted as `n (x%)`.
- Each plot shows the squared Mahalanobis distance of every instance (by index), the outliers as crosses and the threshold as a dashed line.

---

## Reproducibility notes

- **Fixed seeds.** Every run uses one of the seeds 1–30 through `rng(seed,'twister')`, reset before the method and again before each classifier. Reruns on the same MATLAB version give the same numbers: on R2025a, the RELINK check above matches the published CSVs exactly. Other MATLAB versions or platforms may give small floating-point differences, which can occasionally change a last digit or a W/T/L count.
- **Run from the repository root.** All paths are built from `pwd` (`data/...`, `Results/...`). Running from another folder does not always stop with an error: `main.m` prints `No files found` warnings, creates an empty `Results/` folder in the current folder and still ends with `ALL EXPERIMENTS COMPLETE`, and `outlier_plot.m` writes empty summary CSVs into a new `outlier_plots/` folder there. `friedman_nemenyi_figure`, `sensitivity_figure_from_csv` and `target_clustering_summary` stop with an error.
- **Close Excel before running.** While a data file is open in Excel, Excel creates a lock file such as `data/NASA/~$CM1.xlsx`. `main.m`, `sensitivity_figure.m` and `outlier_plot.m` read every `*.xlsx` in the family folder, so the lock file is taken as a project and loading fails (`Unable to open file ... as a workbook`). Close Excel and delete any leftover `~$*.xlsx` files.
- **Do not open and re-save the result CSVs in Excel or Numbers.** Spreadsheet programs can turn cells such as `2/2/2` into dates and can change number formats, delimiters or encoding, which breaks `friedman_nemenyi_figure.m` and the comparisons above. View them in a text editor, or restore them with `git checkout -- Results`.
- **RELINK source A is not oversampled.** Project A has more defective (98) than clean (96) instances, so with `samplingRatio = 1.0` no synthetic instances are created. OTOMO, DPP, ADASYN, SMOTE and TOMO then train on exactly the same data as CORAL, and for A → S and A → Z all six give identical results (only TCA+ differs). This is expected, and it is why most RELINK `W/T/L` entries contain 2 ties.
- **TCA+ needs `tca.m` and `normal.m`**, both included. `normal.m` implements the TCA+ normalization rules of Nam et al. (2013).
- **Failed runs do not stop `main.m`.** An error in one run is reported as a warning (`... seed=... failed: ...`) and that run is left out of the mean; a `writeMethodSummary:MissingRuns` warning then says that some pairs have fewer than 30 runs, and pairs with no run at all are counted in extra `Valid tasks` rows. The RELINK check above runs without such warnings, and the published CSVs have no `Valid tasks` rows.

---

## Using OTOMO on your own data

OTOMO is a single function, [`otomo.m`](otomo.m):

```matlab
[Xs_bal,Ys_bal,purity,ARI,NMI] = otomo(Xs,Ys,Xt,Yt,gamma,ratio,conf,seed)
```

| Argument | Meaning |
|---|---|
| `Xs` | Source features, n_s x p numeric matrix |
| `Ys` | Source labels, n_s x 1, with `1` = defective (the minority class) and `0` = clean |
| `Xt` | Target features, n_t x p, same p features as the source. Target labels are not needed |
| `Yt` | Target labels, used **only** to print and return the clustering diagnostics (`purity`, `ARI`, `NMI`); they do not affect the oversampling. Pass `[]` if the target is unlabelled |
| `gamma` | Weight in [0, 1] between source-minority distances (`gamma`) and distances to the target (`1 - gamma`) when choosing neighbours. The published results use `0.1`; `[]` gives the default `0.5` |
| `ratio` | Target defective/clean ratio after oversampling; `1.0` (default) balances the classes |
| `conf` | Confidence level of the chi-square threshold for Mahalanobis outliers; default `0.995` |
| `seed` | Random seed, **required** (calls `rng(seed,'twister')`) |

| Output | Meaning |
|---|---|
| `Xs_bal`, `Ys_bal` | The original source followed by the synthetic defective instances (label 1) |
| `purity`, `ARI`, `NMI` | Agreement between the target clusters and `Yt` (fractions, not percent); `NaN` when `Yt` is empty or no oversampling was done |

When oversampling is needed, the source must contain at least 7 defective instances: OTOMO picks each neighbour among the 7 nearest, and with 2 to 6 defective instances it usually stops with an index error (it can succeed when only a few synthetic instances are needed); with 7 or fewer, some synthetic instances can be copies of existing ones. `otomo` returns the source unchanged when no synthetic instances are needed (the source already has at least `floor(ratio * clean)` defective instances), when the source has fewer than 2 defective or no clean instances, or when the target has fewer than 2 instances. It prints a few diagnostic lines to the console whenever it creates synthetic instances.

In `main.m`, OTOMO is used as follows: both projects are z-scored separately, the source is oversampled, and the balanced source is aligned to the target with CORAL (`coral(Xs_bal,Xt)`) before a classifier is trained on it. A minimal example, run from the repository root:

```matlab
% Source and target projects: features in all columns but the last, 0/1 label in the last
S = readmatrix(fullfile('data','NASA','CM1.xlsx'));
T = readmatrix(fullfile('data','NASA','PC1.xlsx'));
Xs = S(:,1:end-1);  Ys = S(:,end);
Xt = T(:,1:end-1);  Yt = T(:,end);          % Yt is used only to evaluate at the end

% z-score each project on its own statistics
zs = @(X) (X - mean(X,1)) ./ (std(X,0,1) + (std(X,0,1) < eps));   % SD 1 for constant features, as in main.m
Xs = zs(Xs);  Xt = zs(Xt);

% OTOMO oversampling followed by CORAL alignment (the OTOMO pipeline of main.m)
gamma = 0.1;  ratio = 1.0;  conf = 0.995;  seed = 1;
rng(seed,'twister');
[XsBal,YsBal] = otomo(Xs,Ys,Xt,[],gamma,ratio,conf,seed);
Xtrain = real(coral(XsBal,Xt));

% Train on the balanced, aligned source; predict the target
rng(seed,'twister');
rf = TreeBagger(100,Xtrain,YsBal,'Method','classification', ...
    'NumPredictorsToSample',floor(sqrt(size(Xtrain,2))),'MinLeafSize',1);
[label,score] = predict(rf,Xt);
pred = str2double(label);
[~,~,~,auc] = perfcurve(Yt,score(:,strcmp(rf.ClassNames,'1')),1);
tpr = mean(pred(Yt == 1) == 1);  tnr = mean(pred(Yt == 0) == 0);
fprintf('AUC = %.4f | G-mean = %.4f\n',auc,sqrt(tpr*tnr));
```

With R2025a this prints `AUC = 0.7414 | G-mean = 0.6411` for seed 1, after OTOMO's four diagnostic lines. This is exactly the seed-1 run of `main.m`: repeating the lines from `gamma = ...` to `perfcurve` for `seed = 1:30` and taking `mean` and `std` of the 30 AUC values gives `0.7377 +/- 0.0074`, the published value for CM1 → PC1 (OTOMO, Random Forest, ROC-AUC).

To use your own projects, replace the two `readmatrix` calls, or build `Xs`, `Ys` and `Xt` any other way. First drop rows with missing values, as `main.m` does (for example `S = S(all(isfinite(S),2),:);`). If the target file is unlabelled and has no label column, use all its columns as features (`Xt = T;`) and drop the evaluation lines.

---

## Citation

If you use OTOMO or this implementation in your research, please cite:

> **OTOMO: Outlier-Conscious Transfer-Oriented Minority Oversampling for Cross-Project Defect Prediction**

The complete bibliographic information will be added after publication.

---



## Contact

For questions regarding OTOMO, please contact suravi.akhter@ulab.edu.bd or golam.kibria@ulab.edu.bd.
