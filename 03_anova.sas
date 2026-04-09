/*============================================================
  USE CASE 3 — ANOVA Comparison
  Input  : pbc.pbc_trial (312 trial patients)
  Goal   : Compare lab biomarkers across treatment groups
           and histologic stages; test group differences
============================================================*/

%let proj = /home/u64462473/SPL_Project;
libname pbc "&proj.";

ods graphics on;

/* ============================================================
   SECTION 1 — Baseline Balance: Treatment Groups
   Simple check that D-pen and Placebo arms are comparable
   at enrollment before testing treatment effects.
   ============================================================ */

proc means data=pbc.pbc_trial mean std median min max maxdec=2;
    class trt;
    var age bili albumin chol copper alk_phos ast trig platelet protime;
    title "Section 1 — Baseline Biomarkers by Treatment Group (1=D-Pen, 2=Placebo)";
run;

/* ============================================================
   SECTION 2 — Normality Testing
   Shapiro-Wilk per treatment arm for each biomarker.
   Result informs whether to use parametric ANOVA or
   non-parametric Wilcoxon rank-sum below.
   ============================================================ */

proc univariate data=pbc.pbc_trial normal;
    class trt;
    var bili albumin chol copper alk_phos ast trig platelet protime;
    ods select TestsForNormality;
    title "Section 2 — Shapiro-Wilk Normality Tests by Treatment Group";
run;

/* ============================================================
   SECTION 3 — One-Way ANOVA: Treatment Group (trt)
   Levene test for equal variances included.
   Tukey HSD post-hoc (though only 2 groups here, confirms
   the pairwise p-value matches the overall F-test).
   ============================================================ */

%macro anova_trt(var=, label=);
    proc glm data=pbc.pbc_trial plots=(diagnostics residuals);
        class trt;
        model &var. = trt;
        means trt / hovtest=levene(type=abs) tukey;
        title "ANOVA — &label. by Treatment Group";
    quit;
%mend;

%anova_trt(var=bili,     label=Serum Bilirubin (mg/dL));
%anova_trt(var=albumin,  label=Serum Albumin (g/dL));
%anova_trt(var=chol,     label=Serum Cholesterol (mg/dL));
%anova_trt(var=copper,   label=Urine Copper (ug/day));
%anova_trt(var=alk_phos, label=Alkaline Phosphatase (U/L));
%anova_trt(var=ast,      label=SGOT/AST (U/mL));
%anova_trt(var=trig,     label=Triglycerides (mg/dL));
%anova_trt(var=platelet, label=Platelet Count (x1000/mL));
%anova_trt(var=protime,  label=Prothrombin Time (sec));

/* ============================================================
   SECTION 4 — Non-Parametric Alternative (Wilcoxon)
   Used for skewed biomarkers (bili, copper, alk_phos, ast).
   Wilcoxon rank-sum = Mann-Whitney U for 2-group comparison.
   ============================================================ */

proc npar1way data=pbc.pbc_trial wilcoxon;
    class trt;
    var bili chol copper alk_phos ast trig;
    title "Section 4 — Wilcoxon Rank-Sum: Skewed Biomarkers by Treatment";
run;

/* ============================================================
   SECTION 5 — One-Way ANOVA: Histologic Stage (4 levels)
   Tukey HSD post-hoc identifies which stage pairs differ.
   ============================================================ */

%macro anova_stage(var=, label=);
    proc glm data=pbc.pbc_trial
             plots=(diagnostics meanplot(connect));
        class stage;
        model &var. = stage;
        means stage / hovtest=levene(type=abs) tukey cldiff;
        title "ANOVA — &label. by Histologic Stage (1–4)";
    quit;
%mend;

%anova_stage(var=bili,     label=Serum Bilirubin (mg/dL));
%anova_stage(var=albumin,  label=Serum Albumin (g/dL));
%anova_stage(var=protime,  label=Prothrombin Time (sec));
%anova_stage(var=copper,   label=Urine Copper (ug/day));
%anova_stage(var=alk_phos, label=Alkaline Phosphatase (U/L));

/* ============================================================
   SECTION 6 — Two-Way ANOVA: Treatment × Stage Interaction
   Tests whether treatment effect on bilirubin and albumin
   differs across disease stages.
   ============================================================ */

proc glm data=pbc.pbc_trial plots=(intplot);
    class trt stage;
    model bili = trt | stage;           /* trt stage trt*stage */
    lsmeans trt*stage / pdiff adjust=tukey slice=trt;
    title "Section 6 — Two-Way ANOVA: Bilirubin ~ Treatment × Stage";
quit;

proc glm data=pbc.pbc_trial plots=(intplot);
    class trt stage;
    model albumin = trt | stage;
    lsmeans trt*stage / pdiff adjust=tukey slice=trt;
    title "Section 6 — Two-Way ANOVA: Albumin ~ Treatment × Stage";
quit;

/* ============================================================
   SECTION 7 — Summary Table
   Collect means + p-values in one SQL table for reporting.
   ============================================================ */

proc sql;
    title "Section 7 — Biomarker Means by Treatment Group";
    select
        "bili"     as variable length=12,
        mean(case when trt=1 then bili     end) as mean_trt1 format=8.2,
        mean(case when trt=2 then bili     end) as mean_trt2 format=8.2
    from pbc.pbc_trial
    outer union corr
    select "albumin",
        mean(case when trt=1 then albumin  end),
        mean(case when trt=2 then albumin  end)
    from pbc.pbc_trial
    outer union corr
    select "chol",
        mean(case when trt=1 then chol     end),
        mean(case when trt=2 then chol     end)
    from pbc.pbc_trial
    outer union corr
    select "copper",
        mean(case when trt=1 then copper   end),
        mean(case when trt=2 then copper   end)
    from pbc.pbc_trial
    outer union corr
    select "alk_phos",
        mean(case when trt=1 then alk_phos end),
        mean(case when trt=2 then alk_phos end)
    from pbc.pbc_trial
    outer union corr
    select "ast",
        mean(case when trt=1 then ast      end),
        mean(case when trt=2 then ast      end)
    from pbc.pbc_trial
    outer union corr
    select "trig",
        mean(case when trt=1 then trig     end),
        mean(case when trt=2 then trig     end)
    from pbc.pbc_trial
    outer union corr
    select "platelet",
        mean(case when trt=1 then platelet end),
        mean(case when trt=2 then platelet end)
    from pbc.pbc_trial
    outer union corr
    select "protime",
        mean(case when trt=1 then protime  end),
        mean(case when trt=2 then protime  end)
    from pbc.pbc_trial;
quit;

ods graphics off;
title;
