# SAS Healthcare Data Pipeline — PBC Analysis

**Team:** Team Techies
**Members:** Rucha Keni, Satvik Nayak, Likitha Sri Kode, Tye Jamian

---

A structured SAS pipeline for clinical data analysis using the **Primary Biliary Cirrhosis (PBC)** dataset from the Mayo Clinic (418 patients, 312 randomized to trial). The pipeline follows a modular use-case design, progressing from raw data ingestion through statistical comparison.

---

## Dataset

**PBC (Primary Biliary Cirrhosis)** — `pbc.csv`

| Attribute | Detail |
|-----------|--------|
| Total rows | 418 patients |
| Trial patients | 312 (randomized to D-Penicillamine vs. Placebo) |
| Registry-only | 106 (non-trial, missing 8+ clinical columns) |
| Key outcome | Survival status (0=Censored, 1=Transplant, 2=Dead) |

---

## Current Progress

### Use Case 1 — Data Upload and Validation (`01_upload_validate.sas`)

Ingests the raw CSV and performs comprehensive data quality checks before any transformation.

**Steps:**
1. Import `pbc.csv` via `PROC IMPORT` with full row guessing (`guessingrows=418`)
2. `PROC CONTENTS` — verify variable names (dot-notation auto-converted to underscores)
3. `PROC FREQ` — frequency distributions for all categorical variables (`status`, `trt`, `sex`, `ascites`, `hepato`, `spiders`, `edema`, `stage`)
4. `PROC MEANS` — descriptive statistics (n, nmiss, min, mean, median, max, std) for all numeric lab values
5. Missing pattern analysis — cross-tab of `trt × ascites` to surface the 106 non-trial rows
6. Range validation — flag physiologically impossible values for `bili`, `albumin`, `protime`, `age`, `chol`
7. Duplicate ID check via `PROC SORT NODUPKEY`
8. `PROC UNIVARIATE` — histograms with normal/kernel overlay for key lab distributions

---

### Use Case 2 — Cleaning and Preprocessing (`02_clean_preprocess.sas`)

Transforms the raw data into two analysis-ready datasets with full imputation and derived features.

**Steps:**
1. Drop redundant `rownames` column (duplicate of `id`)
2. Flag 106 registry-only patients (`non_trial = 1`) for separation
3. Impute missing values using **trial-patient statistics only**:
   - Continuous variables (`chol`, `trig`, `copper`, `alk_phos`, `ast`, `platelet`, `protime`) → **median**
   - Ordinal variable (`stage`) → **mode**
   - Imputation flags created per variable to track which rows were filled
4. Apply variable labels across all 20 columns
5. Create derived variables:
   - `death` — binary death indicator (status = 2)
   - `agegroup` — decade buckets (1=<40 through 5=70+)
   - `bili_cat` — bilirubin severity (1=Normal → 4=Severe, clinical thresholds)
   - `low_albumin` — hypoalbuminemia flag (< 3.5 g/dL)
   - `female` — numeric sex indicator for regression
6. Save two output datasets:
   - `pbc.pbc_all` — all 418 rows (imputed, flagged)
   - `pbc.pbc_trial` — 312 trial patients only

**Post-clean validation:** zero-missing confirmation, treatment arm balance, imputation volume summary.

---

### Use Case 3 — ANOVA Comparison (`03_anova.sas`)

Tests biomarker differences across treatment groups and disease stages using parametric and non-parametric methods.

**Input:** `pbc.pbc_trial` (312 patients)

**Sections:**

| Section | Analysis |
|---------|----------|
| 1 | Baseline balance — means/SD for all biomarkers by treatment arm |
| 2 | Normality testing — Shapiro-Wilk per arm per biomarker |
| 3 | One-way ANOVA by treatment group — Levene variance test + Tukey HSD post-hoc (macro loop over 9 biomarkers) |
| 4 | Non-parametric Wilcoxon rank-sum — for skewed biomarkers (`bili`, `chol`, `copper`, `alk_phos`, `ast`, `trig`) |
| 5 | One-way ANOVA by histologic stage (4 levels) — Tukey HSD with confidence limit differences |
| 6 | Two-way ANOVA — Treatment × Stage interaction for bilirubin and albumin (interaction plots) |
| 7 | Summary SQL table — biomarker means side-by-side by treatment group |

---

## Repository Structure

```
SAS_Health_Report/
├── pbc.csv                    # Raw dataset (418 patients)
├── 01_upload_validate.sas     # Use Case 1: Ingestion & QC
├── 02_clean_preprocess.sas    # Use Case 2: Cleaning & Feature Engineering
├── 03_anova.sas               # Use Case 3: ANOVA Statistical Testing
├── build_docs.py              # Documentation build script
└── build_progress_report.py   # Progress report generator
```

---

## SAS Environment

- **Platform:** SAS OnDemand for Academics (SAS 9.4)
- **Library path:** `/home/u64462473/SPL_Project`
- **ODS Graphics:** enabled for diagnostic and interaction plots

---

## Planned Use Cases

- Use Case 4 — Survival Analysis (Kaplan-Meier, Cox regression)
- Use Case 5 — Logistic Regression (mortality prediction)
- Use Case 6 — Reporting & Output (ODS PDF/RTF export)
