# SAS Visual Analytics - PBC Clinical Trial Dashboard

## Dashboard Link

[View the Interactive Dashboard](https://vfl-072.engage.sas.com/SASVisualAnalytics/)

## Project Overview

This interactive dashboard presents the statistical analysis of a randomized, double-blind, placebo-controlled clinical trial evaluating D-Penicillamine for the treatment of Primary Biliary Cirrhosis (PBC). The dataset contains 312 trial patients from the Mayo Clinic.

## Dashboard Pages

### Page 1 - Patient Overview
- Treatment group distribution (158 Drug, 154 Placebo)
- Disease stage distribution across all patients
- Age distribution by treatment group

### Page 2 - Clinical Analysis
- Mean bilirubin and albumin trends across disease stages (Dual Axis Chart)
- Ascites rates by treatment group
- Average bilirubin heat map by treatment group and disease stage

### Page 3 - Outcomes
- Death rates by treatment group and disease stage
- Key findings summary

## Interactive Features

- **Global Filter:** Toggle between D-Penicillamine, Placebo, or All patients across every page
- **Chart-to-Chart Linking:** Click on any chart element to filter other charts on the same page
- **Page Navigation:** Click on chart elements to navigate to related pages with filters applied
- **Hover Tooltips:** Hover over any data point for detailed information

## Key Findings

1. Ascites rates were higher in D-Penicillamine patients compared to Placebo, raising potential safety concerns about the drug worsening fluid retention
2. Disease stage is the strongest predictor of liver function decline (p < 0.0001 for bilirubin, albumin, copper, and prothrombin time)
3. Stage 4 patients had significantly worse outcomes than all other stages across every biomarker
4. Adverse event rates were comparable between Drug and Placebo groups with no additional safety concerns
5. Treatment groups were well balanced at baseline confirming proper randomization
6. Conclusion: D-Penicillamine did not demonstrate efficacy in treating Primary Biliary Cirrhosis

## Dataset

- **Source:** PBC Trial Dataset (Mayo Clinic)
- **Patients:** 312 randomized trial participants
- **Treatment Arms:** D-Penicillamine (n=158) vs Placebo (n=154)
- **Variables:** 35 including demographics, lab values, clinical indicators, and survival outcomes
- **Disease:** Primary Biliary Cirrhosis (chronic liver disease)

## Technology Stack

- **Data Processing:** SAS Studio (SAS OnDemand for Academics)
- **Statistical Analysis:** SAS 9.4 (PROC GLM, PROC LIFETEST, PROC PHREG, PROC FREQ)
- **Report Generation:** SAS ODS (PDF output)
- **Visualization:** SAS Viya Visual Analytics
- **Dashboard:** SAS Visual Analytics (interactive, web-based)

## Team Members

- Rucha Keni
- Satvik Nayak
- Likhita Sri Kode
- Tye Jamian
