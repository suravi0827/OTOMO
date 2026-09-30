# Statistical Significance Tests

This note explains the Friedman–Nemenyi analysis in [`friedman_nemenyi_figure.m`](../friedman_nemenyi_figure.m), the results it currently produces, and a supplementary head-to-head check. The notation follows the paper: $\eta$ is the number of cross-project pairs, $L$ the number of methods, and $R_j$ the average rank of method $j$.

## 1. Overview

OTOMO is compared with six methods (CORAL, TCA+, DPP, ADASYN, SMOTE, TOMO) on 66 cross-project pairs. The comparison uses two classifiers, Random Forest (RF) and Logistic Regression (LR), and four metrics: ROC-AUC, G-mean, F1-score and MCC. No method wins on every pair, so we test whether the differences are larger than chance would explain. For each classifier–metric combination, a Friedman test checks the null hypothesis that all seven methods perform equivalently. If that hypothesis is rejected, the Nemenyi post-hoc test identifies which pairs of methods differ. The script draws the critical-difference (CD) diagrams used in the paper and writes a summary table.

## 2. Data used

- **Pairs.** There are 66 source → target pairs, one for every ordered pair of distinct projects within a dataset family: AEEEM (5 projects, 20 pairs), NASA (5, 20), RELINK (3, 6) and SOFTLAB (5, 20).
- **Methods.** $L = 7$: OTOMO, CORAL, TCA+, DPP, ADASYN, SMOTE and TOMO.
- **Values.** For each pair and method, the test uses the mean of the metric over 30 random seeds. The metric summary CSVs store this as `mean +/- sd` with four decimals. The standard deviation is not used.
- **Tests.** One test is run per classifier × metric, giving 2 × 4 = 8 tests. Each test uses a $66 \times 7$ matrix: the rows (the blocks) are the pairs and the columns are the methods.
- **Nothing aggregated enters a test.** The script skips the `W/T/L`, `Average (%)`, `Overall W/T/L` and `Overall Average (%)` rows of the CSVs. No averages over pairs and no win/tie/loss counts are used.

The script also repeats each test within each family (the `AEEEM`, `NASA`, `RELINK` and `SOFTLAB` rows of the summary table). The figures and Section 6 use all 66 pairs.

## 3. Procedure (as in the paper)

The settings are $\eta = 66$ pairs, $L = 7$ methods and $\alpha = 0.05$.

**Ranking.** On every pair, the methods are ranked by metric value, with rank 1 for the best (all four metrics are higher-is-better). Tied methods receive the average of the ranks they occupy; for example, two methods tied for second place both get 2.5. $R_j$ is the average rank of method $j$ over the $\eta$ pairs.

**Friedman test.** The null hypothesis $H_0$ is that the methods perform equivalently. Under $H_0$, every $R_j$ should be close to $(L+1)/2 = 4$. The statistic is

$$\chi^2_F = \frac{12\eta}{L(L+1)}\left[\sum_{j=1}^{L} R_j^2 - \frac{L(L+1)^2}{4}\right]$$

and under $H_0$ it approximately follows a $\chi^2$ distribution with $L-1 = 6$ degrees of freedom.

**Iman–Davenport decision.** The decision uses the F-approximation

$$F_F = \frac{(\eta-1)\,\chi^2_F}{\eta(L-1) - \chi^2_F},$$

which is compared with an F-distribution with $L-1$ and $(L-1)(\eta-1)$ degrees of freedom, here $F(6, 390)$. $H_0$ is rejected when $F_F > F_{\text{crit}}$, the upper-tail critical value at $\alpha$: $F_{\text{crit}}(0.05;\,6,\,390) = 2.1218$. The summary table also reports the $\chi^2$ p-value (`Friedman_p`), but the decision is based on $F_F$, as in the paper.

**Nemenyi post-hoc test.** This test is applied only when $H_0$ is rejected. Methods $i$ and $j$ differ significantly when $\lvert R_i - R_j \rvert > \mathrm{CD}$, where

$$\mathrm{CD} = q_{\alpha,L}\sqrt{\frac{L(L+1)}{6\eta}}.$$

$q_{\alpha,L}$ is the upper-$\alpha$ point of the studentized range for $L$ groups and infinite degrees of freedom, divided by $\sqrt{2}$. The script computes it numerically as $q_{0.05,7} = 2.9483$. Demšar (2006, Table 5) lists 2.949, which comes from a rounded tabulated value. With $L = 7$ and $\eta = 66$, $\mathrm{CD} = 2.9483 \times \sqrt{56/396} = 1.1087$. When $H_0$ is not rejected, the script gives no post-hoc verdicts (`not tested (Friedman n.s.)`).

## 4. Why this design is reasonable

- **No normality assumption.** The per-pair metric values are bounded, spread differently across families, and often not normally distributed. A rank test avoids the normality and equal-variance assumptions of a repeated-measures ANOVA.
- **Pairs of different difficulty become comparable.** Under RF, for example, the best ROC-AUC achieved on a pair ranges from 0.58 (ar3 → ar6) to 0.95 (ar4 → ar5). Ranking within each pair removes this level and scale, so every pair has equal weight and a few extreme pairs cannot dominate the result, as they can in an average.
- **Pairs as blocks.** All methods are evaluated on the same 66 pairs. This is the complete block design the Friedman test requires; independence of the blocks holds only approximately (Section 8).
- **Error control over all comparisons.** Testing the 21 method pairs one by one at $\alpha$ would inflate false positives. An omnibus test followed by Nemenyi keeps the family-wise error rate at $\alpha$ across all 21 comparisons.
- **Iman–Davenport.** The $\chi^2$ approximation of the Friedman statistic is known to be conservative; the F-approximation is more accurate (Iman & Davenport, 1980).
- **Standard practice.** Demšar (2006) recommends Friedman + Nemenyi with CD diagrams for comparing several methods over multiple data sets. The procedure is widely used in defect prediction, so the analysis follows the convention of prior work.

## 5. Reading the critical-difference diagrams

The script produces one figure per classifier, `RF_Overall_Friedman_Nemenyi_CD` and `LR_Overall_Friedman_Nemenyi_CD`. Each has four panels, (a) ROC-AUC, (b) G-mean, (c) F1-score and (d) MCC, and every panel uses all 66 pairs.

- **Axis.** The axis shows average rank from 7 (left) to 1 (right), so better methods are further right. Each method's line meets the axis at its $R_j$. The better half of the methods is labelled on the right and the rest on the left.
- **CD bar.** The bar labelled `CD=1.1087` above the axis shows the critical difference on the rank scale.
- **Group bars.** The thick bars below the axis join groups of methods whose average ranks differ by at most CD. Each bar is a maximal group. Methods joined by a bar are not significantly different; two methods that share no bar are significantly different. The bar colours only separate overlapping groups and carry no meaning.
- **OTOMO** is labelled in blue (`highlightMethod`).
- **Friedman n.s.** If a Friedman test did not reject $H_0$, the panel caption would read `(Friedman n.s., p = …)`, with the Iman–Davenport p-value, and its bars should not be interpreted. None of the current 66-pair panels is in this situation.

Example: in RF (a) ROC-AUC, OTOMO ($R_j = 3.05$) shares a bar with DPP, CORAL and SMOTE, so it is not significantly different from them. It shares no bar with ADASYN (4.27), TCA+ (4.86) or TOMO (5.61), so it is significantly better than those three.

## 6. Current results over all 66 pairs

These results were produced by `friedman_nemenyi_figure.m` with its default settings ($\alpha = 0.05$, `tieCorrection = false`) from the `gamma_001_0.1` metric summaries on 2026-09-30. They are the `Dataset = All` rows of `Friedman_Nemenyi_Summary.csv`. In all eight tests, $\eta = 66$, $L = 7$, $F_{\text{crit}} = 2.1218$ and $\mathrm{CD} = 1.1087$.

| Classifier | Metric | $\chi^2_F$ | $F_F$ | $p$ of $F_F$ | OTOMO $R_j$ (position) | OTOMO significantly better than | Not significant |
|---|---|---:|---:|---:|---|---|---|
| RF | ROC-AUC | 77.06 | 15.71 | 3.6e-16 | 3.05 (1st) | TCA+, ADASYN, TOMO | CORAL, DPP, SMOTE |
| RF | G-mean | 183.70 | 56.24 | 6.7e-50 | 2.47 (1st) | CORAL, TCA+, DPP | ADASYN, SMOTE, TOMO |
| RF | F1-score | 129.13 | 31.45 | 7.9e-31 | 2.71 (1st) | CORAL, TCA+, DPP | ADASYN, SMOTE, TOMO |
| RF | MCC | 65.95 | 12.99 | 2.1e-13 | 2.92 (1st) | TCA+, TOMO | CORAL, DPP, ADASYN, SMOTE |
| LR | ROC-AUC | 101.58 | 22.43 | 1.0e-22 | 2.15 (1st) | all six | none |
| LR | G-mean | 115.55 | 26.78 | 1.0e-26 | 2.89 (2nd; DPP 2.85) | CORAL, TCA+, TOMO | DPP, ADASYN, SMOTE |
| LR | F1-score | 69.96 | 13.95 | 2.2e-14 | 2.53 (1st) | CORAL, TCA+, ADASYN, TOMO | DPP, SMOTE |
| LR | MCC | 64.00 | 12.53 | 6.3e-13 | 2.41 (1st) | all six | none |

- $H_0$ is rejected in all eight tests.
- OTOMO has the best average rank in 7 of 8 tests. The exception is LR G-mean, where DPP is ahead by 0.04, which is not significant.
- Across the 48 OTOMO-vs-baseline comparisons, OTOMO is significantly better in 30 and never significantly worse. By baseline: TCA+ 8/8, CORAL 6/8, TOMO 6/8, DPP 4/8, ADASYN 4/8 and SMOTE 2/8.
- **Per-family rows.** The summary table also has rows for each family: AEEEM, NASA and SOFTLAB ($\eta = 20$, $F_{\text{crit}} = 2.1791$, $\mathrm{CD} = 2.0141$) and RELINK ($\eta = 6$, $F_{\text{crit}} = 2.4205$, $\mathrm{CD} = 3.6772$). With this few pairs the CD is large (2.01 and 3.68 ranks), so the family-level post-hoc tests have little power; the omnibus test rejects $H_0$ in 28 of the 32 family-level tests. For LR on RELINK, the Friedman test does not reject $H_0$ for any of the four metrics.

## 7. Supplementary check (not in the paper): Wilcoxon signed-rank test

This analysis is supplementary. It is not produced by `friedman_nemenyi_figure.m` and is not part of the procedure described in the paper. It asks a narrower question: taking one baseline at a time, are OTOMO's 66 per-pair values systematically higher than that baseline's?

- **Test.** A two-sided Wilcoxon signed-rank test on the 66 paired differences (OTOMO − baseline). Unlike the Friedman ranks, it uses the size of each difference, and it ignores the other five methods. Zero differences are dropped. RELINK A → S and A → Z are exact ties: every method except TCA+ gets the same value there, in every file. Each comparison therefore has 63–66 non-zero differences, 64 in most cases.
- **Multiple comparisons.** Holm's step-down correction is applied over the six baselines within each classifier × metric.
- **Effect size.** Cliff's $\delta = P(X > Y) - P(X < Y)$, computed over all $66 \times 66$ combinations of an OTOMO value $X$ and a baseline value $Y$. The thresholds of Romano et al. (2006) are used: $\lvert\delta\rvert < 0.147$ is negligible, $< 0.33$ small, $< 0.474$ medium, and anything larger is large. $\delta$ ignores the pairing, so it is dominated by the differences between pairs (Section 4): it indicates practical size, not how consistently OTOMO wins on the same pair.

**Agreement with Nemenyi.** The two procedures agree on 37 of 48 comparisons: in 30 both find OTOMO significantly better, and in 7 neither finds a significant difference. In the other 11, Wilcoxon + Holm finds OTOMO significantly better and Nemenyi does not; the reverse never happens. All 48 median differences are positive. Of the 11 extra results, 8 have a negligible Cliff's $\delta$ and 3 a small one, although OTOMO wins on 39–49 of the 63–64 non-tied pairs in each (paired, matched-pairs rank-biserial effect 0.30–0.55):

| Classifier | Metric | Baseline | Holm $p$ | Median difference | Cliff's $\delta$ | Magnitude |
|---|---|---|---:|---:|---:|---|
| RF | G-mean | SMOTE | 1.63e-02 | +0.0085 | +0.048 | negligible |
| RF | G-mean | TOMO | 1.59e-02 | +0.0147 | +0.232 | small |
| RF | F1-score | ADASYN | 3.57e-02 | +0.0059 | +0.069 | negligible |
| RF | F1-score | SMOTE | 3.04e-02 | +0.0081 | +0.058 | negligible |
| RF | F1-score | TOMO | 1.25e-02 | +0.0188 | +0.201 | small |
| RF | MCC | CORAL | 6.72e-03 | +0.0315 | +0.225 | small |
| RF | MCC | DPP | 1.78e-03 | +0.0228 | +0.109 | negligible |
| RF | MCC | ADASYN | 6.72e-03 | +0.0132 | +0.109 | negligible |
| RF | MCC | SMOTE | 6.72e-03 | +0.0109 | +0.062 | negligible |
| LR | F1-score | DPP | 9.89e-04 | +0.0145 | +0.071 | negligible |
| LR | F1-score | SMOTE | 2.95e-04 | +0.0134 | +0.083 | negligible |

The Wilcoxon test is more sensitive for three reasons. It adjusts for 6 comparisons (OTOMO against each baseline) instead of all 21, its result does not depend on the other five methods (Benavoli et al., 2016), and it ranks the sizes of the 66 differences instead of using only each pair's ordering of seven methods. This lets it detect small but consistent advantages that the Nemenyi test does not flag. The following MATLAB code reproduces this analysis when run from a folder that contains `Results/Metric_Summaries/`:

```matlab
folder = fullfile('Results','Metric_Summaries'); tag = "gamma_001_0.1";
for clf = ["Random_Forest","Logistic_Regression"]
    for met = ["ROC_AUC","G_Measure","F1_Score","MCC_Score"]
        lines = splitlines(string(fileread(fullfile(folder,clf + "_All_Datasets_" + met + "_" + tag + ".csv"))));
        methods = strtrim(split(lines(1),",")).'; methods = methods(3:end);
        X = zeros(0,numel(methods));
        for i = 2:numel(lines)
            c = strtrim(split(lines(i),",")).';
            if numel(c) < numel(methods) + 2 || c(1) == "" || c(2) == "", continue; end   % family, W/T/L, average rows
            X(end+1,:) = str2double(extractBefore(c(3:end) + " "," "));                    %#ok<SAGROW> mean of 'mean +/- sd'
        end
        o = X(:,methods == "OTOMO"); base = methods(methods ~= "OTOMO");
        p = zeros(1,numel(base)); d = p; md = p;
        for b = 1:numel(base)
            y = X(:,methods == base(b));
            p(b) = signrank(o,y);                    % two-sided, zero differences dropped
            md(b) = median(o - y);                   % median paired difference
            d(b) = mean(sign(o - y.'),'all');        % Cliff's delta over all 66 x 66 combinations
        end
        [ps,ix] = sort(p); pHolm(ix) = min(1,cummax(ps .* (numel(p):-1:1)));   % Holm over the 6 baselines
        for b = 1:numel(base)
            fprintf('%-19s %-9s OTOMO vs %-6s p_Holm=%.2e  median diff=%+.4f  Cliff d=%+.3f  significant=%d\n', ...
                clf,met,base(b),pHolm(b),md(b),d(b),pHolm(b) < 0.05);
        end
    end
end
```

The Holm $p$ values in the table above are the output of this code. MATLAB's `signrank` uses the normal approximation here, because each comparison has more than 15 non-zero differences. An independent SciPy implementation gives p-values within 2% of these and the same 48 decisions.

## 8. Limitations and assumptions

- **The pairs are not independent.** The 66 pairs come from 18 projects. Each project appears in 8 pairs (4 as source, 4 as target), and each RELINK project in 4. The tests treat the pairs as independent blocks, so the effective sample size is likely below 66, and the p-values and CD are likely too optimistic by an unknown amount. Treating pairs as blocks is standard practice in cross-project defect prediction (CPDP). As a sensitivity check, if the 66 pairs carried only the information of 18 independent pairs (one per project: $\chi^2_F$ scaled by 18/66, CD = 2.12), all eight omnibus tests would still reject $H_0$ (smallest $F_F$ = 3.28 against $F_{\text{crit}}$ = 2.19), but only 17 of OTOMO's 30 significant post-hoc wins would remain.
- **Nemenyi is conservative.** It adjusts for all 21 pairwise comparisons, although the main interest is OTOMO against each baseline. For comparisons against a control method, Demšar (2006) describes the more powerful Bonferroni–Dunn test. It uses the same CD formula with $q_{0.05} = 2.6383$ for $L = 7$ (Demšar lists 2.638), which gives $\mathrm{CD}_{\mathrm{BD}} = 0.9921$. With this test, OTOMO would also be significantly better than TOMO for RF F1-score (rank difference 1.09) and than SMOTE for LR F1-score (1.00); no other verdict would change. Both comparisons are among the 11 extra Wilcoxon results. We report Nemenyi to match the paper.
- **Seed variability is not modelled.** Each pair enters the test as a mean over 30 seeds, so the tests address variation across pairs, not across seeds within a pair. Because the means are stored with four decimals, methods whose means agree to four decimals count as tied.
- **Ties and the `tieCorrection` switch.** Ties receive average ranks, and the paper's $\chi^2_F$ has no tie correction. Ties reduce the variance of the ranks, so the uncorrected statistic is smaller (slightly for the 66 pairs, markedly within RELINK), which makes the test conservative. In each test, 3–6 of the 66 pairs contain ties. Setting `tieCorrection = true` switches to the tie-corrected statistic used by MATLAB's `friedman`, R's `friedman.test` and SciPy's `friedmanchisquare`, where $r_{ij}$ is the rank of method $j$ on pair $i$:

  $$\chi^2_{F,\text{ties}} = \frac{(L-1)\,\eta^2 \sum_{j}\left(R_j - \frac{L+1}{2}\right)^2}{\sum_{i,j} r_{ij}^2 - \eta L (L+1)^2/4}$$

  On the 66 pairs, the correction raises $\chi^2_F$ slightly (for example, from 77.06 to 78.60 for RF ROC-AUC) and changes no decision; CD is unaffected. Within RELINK ($\eta = 6$; A → S and A → Z are six-way ties), it would reject $H_0$ for LR ROC-AUC, G-mean and MCC, which the paper's formula does not.
- **Several separate tests.** There are eight tests (2 classifiers × 4 metrics), each at $\alpha = 0.05$ with no correction across them. The metrics are computed from the same models and are correlated. The eight results should therefore be read as consistency across evaluation settings, not as eight independent confirmations.

## 9. How to run

In MATLAB R2021a or later with the Statistics and Machine Learning Toolbox, run the script from the repository root. The script uses `pwd` and starts with `clear; close all; clc`.

```matlab
friedman_nemenyi_figure
```

**Inputs.** The script reads eight CSVs in `Results/Metric_Summaries/` named `<Classifier>_All_Datasets_<Metric>_<gammaTag>.csv`, where `<Classifier>` is `Random_Forest` or `Logistic_Regression` and `<Metric>` is `ROC_AUC`, `G_Measure`, `F1_Score` or `MCC_Score`; for example, `Random_Forest_All_Datasets_ROC_AUC_gamma_001_0.1.csv`. These files are written by the experiment driver. The CSVs used in Sections 6 and 7 are not included in this repository; any CSVs in this format can be used:

```text
Source,Target,OTOMO,CORAL,TCA+,DPP,ADASYN,SMOTE,TOMO
NASA,,,,,,,,
CM1,MW1,0.7554 +/- 0.0103,0.7272 +/- 0.0249,0.7075 +/- 0.0174,...
...
W/T/L,,,9/0/11,11/0/9,...
```

The format rules are:

- The header starts with `Source,Target`, followed by the method names.
- A row with only its first cell filled starts a family.
- Each pair row holds a source, a target and one cell per method. Only the leading number of each cell (the mean) is read.
- Rows with an empty second cell (W/T/L, averages) are ignored.
- Pairs with a missing value are dropped with a warning, and duplicate pairs raise an error.
- All eight files must have the same method columns.

**Outputs** go to `Results/Significance_Tests/<gammaTag>/`:

- `RF_Overall_Friedman_Nemenyi_CD.eps` / `.png` (600 dpi) and `LR_Overall_Friedman_Nemenyi_CD.eps` / `.png`: the CD figures described in Section 5.
- `Friedman_Nemenyi_Summary.csv`: 40 rows, one per metric × classifier × group (the four families plus `All`). After `Metric`, `Classifier` and `Dataset`, the columns are `N_Pairs` ($\eta$), `K_Methods` ($L$), `Friedman_Chi2`, `Friedman_p`, `Friedman_Significant` (χ² decision, for reference only), `ImanDavenport_F`, `ImanDavenport_Fcrit`, `ImanDavenport_p`, `Reject_H0` (the Iman–Davenport decision that gates Nemenyi), `Nemenyi_q`, `CD`, `Best_Method`, `AvgRank_<method>`, and one `OTOMO_vs_<baseline>` verdict per baseline: `significantly better`, `significantly worse`, `no significant difference`, or `not tested (Friedman n.s.)` (or `n/a` if a method is missing). The console shows $\eta$, $L$, $\chi^2_F$ and its p-value, $F_F$, $F_{\text{crit}}$, the decision and CD.

**Settings** are at the top of the script:

| Variable | Default | Effect |
|---|---|---|
| `alpha` | `0.05` | Significance level for $F_{\text{crit}}$ and $q_{\alpha,L}$ |
| `tieCorrection` | `false` | `false`: the paper's formula; `true`: tie-corrected $\chi^2_F$ (Section 8) |
| `datasets` | `["AEEEM","RELINK","NASA","SOFTLAB"]` | Families tested separately and pooled into `All`; a listed family missing from a CSV is an error, and other families are excluded with a warning |
| `highlightMethod` | `"OTOMO"` | Method drawn in blue and used as the reference in the verdict columns |
| `gammaTag` | `""` | Which `gamma_*` file set to use; empty means the only one present (error if there are several) |
| `summaryFolder`, `outputFolder` | `Results/Metric_Summaries`, `Results/Significance_Tests` | Input and output folders |

## 10. References

- Benavoli, A., Corani, G., & Mangili, F. (2016). Should we really use post-hoc tests based on mean-ranks? *Journal of Machine Learning Research*, 17(5), 1–10.
- Cliff, N. (1993). Dominance statistics: Ordinal analyses to answer ordinal questions. *Psychological Bulletin*, 114(3), 494–509.
- Demšar, J. (2006). Statistical comparisons of classifiers over multiple data sets. *Journal of Machine Learning Research*, 7, 1–30.
- Friedman, M. (1937). The use of ranks to avoid the assumption of normality implicit in the analysis of variance. *Journal of the American Statistical Association*, 32(200), 675–701.
- Friedman, M. (1940). A comparison of alternative tests of significance for the problem of m rankings. *The Annals of Mathematical Statistics*, 11(1), 86–92.
- Holm, S. (1979). A simple sequentially rejective multiple test procedure. *Scandinavian Journal of Statistics*, 6(2), 65–70.
- Iman, R. L., & Davenport, J. M. (1980). Approximations of the critical region of the Friedman statistic. *Communications in Statistics – Theory and Methods*, 9(6), 571–595.
- Nemenyi, P. B. (1963). *Distribution-free multiple comparisons*. PhD thesis, Princeton University.
- Romano, J., Kromrey, J. D., Coraggio, J., & Skowronek, J. (2006). Appropriate statistics for ordinal level data: Should we really be using t-test and Cohen's d for evaluating group differences on the NSSE and other surveys? Paper presented at the annual meeting of the Florida Association of Institutional Research.
- Wilcoxon, F. (1945). Individual comparisons by ranking methods. *Biometrics Bulletin*, 1(6), 80–83.
