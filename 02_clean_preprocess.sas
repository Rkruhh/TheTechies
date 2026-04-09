/*
  USE CASE 2 — Cleaning and Preprocessing
  Cleaning steps applied:
    1. Drop redundant rownames column
    2. Exclude / flag the 106 non-trial patients
    3. Impute missing values (median for continuous, mode for stage)
    4. Label all variables + create derived variables
*/

%let proj = /home/u64462473/SPL_Project;
libname pbc "&proj.";

/* 
   STEP 1 — Drop rownames (identical to id, adds no value)
*/
data work.pbc_s1;
    set pbc.pbc_raw (drop=rownames);
run;

%put NOTE: Step 1 complete — rownames column dropped.;

/* 
   STEP 2 — Flag non-trial patients (trt = missing)
*/
data work.pbc_s2;
    set work.pbc_s1;
    non_trial = (trt = .);
    label non_trial = "Non-Trial Patient (1=Registry-only, 0=RCT)";
run;

proc freq data=work.pbc_s2;
    tables non_trial / nocum nopercent;
    title "Step 2 — Trial vs Non-Trial Patient Counts";
run;

/* 
   STEP 3 — Impute missing values
*/

/* 3a. Compute medians from trial patients only */
proc means data=work.pbc_s2 median noprint;
    where non_trial = 0;
    var chol trig copper alk_phos ast platelet protime;
    output out=work.medians(drop=_type_ _freq_)
           median = med_chol med_trig med_copper
                    med_alkphos med_ast med_platelet med_protime;
run;

/* 3b. Push medians into macro variables */
data _null_;
    set work.medians;
    call symputx('med_chol',     put(med_chol,     best12.));
    call symputx('med_trig',     put(med_trig,     best12.));
    call symputx('med_copper',   put(med_copper,   best12.));
    call symputx('med_alkphos',  put(med_alkphos,  best12.));
    call symputx('med_ast',      put(med_ast,      best12.));
    call symputx('med_platelet', put(med_platelet, best12.));
    call symputx('med_protime',  put(med_protime,  best12.));
run;

%put NOTE: Imputation medians (from trial patients):;
%put NOTE:   chol=&med_chol.  trig=&med_trig.  copper=&med_copper.;
%put NOTE:   alk_phos=&med_alkphos.  ast=&med_ast.;
%put NOTE:   platelet=&med_platelet.  protime=&med_protime.;

/* 3c. Compute mode of stage from trial patients */
proc freq data=work.pbc_s2 noprint;
    where non_trial = 0;
    tables stage / out=work.stage_freq(drop=percent);
run;

proc sort data=work.stage_freq; by descending count; run;

data _null_;
    set work.stage_freq (obs=1);
    call symputx('mode_stage', stage);
run;

%put NOTE: Stage mode (from trial patients) = &mode_stage.;

/* 3d. Apply imputation across all 418 rows */
data work.pbc_s3;
    set work.pbc_s2;

    if chol     = . then do; chol     = &med_chol.;     chol_imp     = 1; end;
    if trig     = . then do; trig     = &med_trig.;     trig_imp     = 1; end;
    if copper   = . then do; copper   = &med_copper.;   copper_imp   = 1; end;
    if alk_phos = . then do; alk_phos = &med_alkphos.;  alkphos_imp  = 1; end;
    if ast      = . then do; ast      = &med_ast.;      ast_imp      = 1; end;
    if platelet = . then do; platelet = &med_platelet.; platelet_imp = 1; end;
    if protime  = . then do; protime  = &med_protime.;  protime_imp  = 1; end;
    if stage    = . then do; stage    = &mode_stage.;   stage_imp    = 1; end;

    /* Default imputation flags to 0 */
    array impflags {8} chol_imp trig_imp copper_imp alkphos_imp
                       ast_imp platelet_imp protime_imp stage_imp;
    do i = 1 to 8;
        if impflags{i} = . then impflags{i} = 0;
    end;
    drop i;

    label
        chol_imp     = "Cholesterol imputed flag"
        trig_imp     = "Triglycerides imputed flag"
        copper_imp   = "Copper imputed flag"
        alkphos_imp  = "Alk. Phosphatase imputed flag"
        ast_imp      = "AST imputed flag"
        platelet_imp = "Platelet imputed flag"
        protime_imp  = "Protime imputed flag"
        stage_imp    = "Stage imputed flag";
run;

/* 
   STEP 4 — Label variables and create derived columns
*/
data work.pbc_s4;
    set work.pbc_s3;

    label
        id       = "Patient ID"
        time     = "Follow-up Time (days)"
        status   = "Status (0=Censored, 1=Transplant, 2=Dead)"
        trt      = "Treatment (1=D-Penicillamine, 2=Placebo)"
        age      = "Age (years)"
        sex      = "Sex (f=Female, m=Male)"
        ascites  = "Ascites Present (0=No, 1=Yes)"
        hepato   = "Hepatomegaly (0=No, 1=Yes)"
        spiders  = "Spider Angiomata (0=No, 1=Yes)"
        edema    = "Edema (0=None, 0.5=Treated/Resolved, 1=Persistent)"
        bili     = "Serum Bilirubin (mg/dL)"
        chol     = "Serum Cholesterol (mg/dL)"
        albumin  = "Serum Albumin (g/dL)"
        copper   = "Urine Copper (ug/day)"
        alk_phos = "Alkaline Phosphatase (U/liter)"
        ast      = "SGOT / AST (U/mL)"
        trig     = "Triglycerides (mg/dL)"
        platelet = "Platelet Count (per cubic mL x1000)"
        protime  = "Prothrombin Time (seconds)"
        stage    = "Histologic Stage (1=earliest, 4=advanced)";

    death = (status = 2);
    label death = "Death Indicator (0=Alive/Transplant, 1=Died)";

    if      age <  40 then agegroup = 1;
    else if age <  50 then agegroup = 2;
    else if age <  60 then agegroup = 3;
    else if age <  70 then agegroup = 4;
    else                   agegroup = 5;
    label agegroup = "Age Group (1=<40, 2=40s, 3=50s, 4=60s, 5=70+)";

    if      bili <  1.2  then bili_cat = 1;   
    else if bili <  3.0  then bili_cat = 2;  
    else if bili < 10.0  then bili_cat = 3;   
    else                      bili_cat = 4;   
    label bili_cat = "Bilirubin Category (1=Normal ─> 4=Severe)";

    low_albumin = (albumin < 3.5);
    label low_albumin = "Low Albumin Flag (<3.5 g/dL)";

    female = (sex = 'f');
    label female = "Female (1=Yes, 0=No)";

run;

/* 
   STEP 5 — Save final datasets
   pbc_all   : all 418 rows (imputed, flagged)
   pbc_trial : 312 trial patients only
*/
data pbc.pbc_all;
    set work.pbc_s4;
run;

data pbc.pbc_trial;
    set work.pbc_s4;
    where non_trial = 0;
run;

%put NOTE: pbc_all   saved — %sysfunc(attrn(%sysfunc(open(pbc.pbc_all)),   nobs)) rows.;
%put NOTE: pbc_trial saved — %sysfunc(attrn(%sysfunc(open(pbc.pbc_trial)), nobs)) rows.;

/* 
   POST-CLEAN VALIDATION
*/

proc means data=pbc.pbc_trial nmiss maxdec=0;
    var time age bili chol albumin copper alk_phos ast
        trig platelet protime stage;
    title "Post-Clean: Missing Count Check — pbc_trial (expect all zeros)";
run;

proc freq data=pbc.pbc_trial;
    tables trt sex stage agegroup bili_cat / nocum nopercent;
    title "Post-Clean: Distributions in Trial Dataset (n=312)";
run;

proc means data=pbc.pbc_trial sum maxdec=0;
    var chol_imp trig_imp copper_imp alkphos_imp
        ast_imp platelet_imp protime_imp stage_imp;
    title "Imputation Counts — pbc_trial (rows that were filled in)";
run;

title;
