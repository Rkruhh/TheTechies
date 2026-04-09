/*
  USE CASE 1 — Data Upload and Validation
*/

/* 1. Global path macro  */
%let proj = /home/u64462473/SPL_Project;
libname pbc "&proj.";

/* 2. Import raw CSV  */
proc import
    datafile = "&proj./pbc.csv"
    out      = pbc.pbc_raw
    dbms     = csv
    replace;
    getnames = yes;
    guessingrows = 418;
run;

data pbc.pbc_raw;
    set pbc.pbc_raw;
    alk_phos = 'alk.phos'n;
    drop 'alk.phos'n;
run;

/* 3. Rename dot-notation variable  */
/* SAS converts dots to underscores on import automatically  */
proc contents data=pbc.pbc_raw varnum;
    title "PBC Raw — Variable Listing";
run;

/* 4. Frequency check — all categorical variables */
proc freq data=pbc.pbc_raw;
    tables status trt sex ascites hepato spiders edema stage
           / missing nocum nopercent;
    title "PBC — Categorical Variable Distributions";
run;

/* 5. Descriptive stats — all numeric variables */
proc means data=pbc.pbc_raw
    n nmiss min mean median max std
    maxdec=2;
    var time age bili chol albumin copper alk_phos ast trig platelet protime;
    title "PBC — Numeric Variable Summary (n, nmiss, range)";
run;

/* 6. Missing value heat-map by variable */
proc freq data=pbc.pbc_raw;
    tables trt*ascites / missing list nocum nopercent;
    title "PBC — Missing Pattern: trt vs ascites (106 non-trial rows)";
run;

/* 7. Range validation — flag physiologically impossible values */
data pbc.pbc_flags;
    set pbc.pbc_raw;

    flag_bili    = (bili    <  0 | bili    > 50);
    flag_albumin = (albumin <  1 | albumin >  6);
    flag_protime = (protime <  5 | protime > 30);
    flag_age     = (age     < 18 | age     > 100);
    flag_chol    = (chol    <  0 | chol    > 2000);

    any_flag = max(flag_bili, flag_albumin, flag_protime,
                   flag_age, flag_chol);
run;

proc freq data=pbc.pbc_flags;
    tables any_flag flag_bili flag_albumin flag_protime
           flag_age flag_chol / nocum nopercent;
    title "PBC — Out-of-Range Flag Counts";
run;

/* 8. Duplicate ID check  */
proc sort data=pbc.pbc_raw out=_sort nodupkey dupout=_dups;
    by id;
run;

data _null_;
    if 0 then set _dups nobs=n;
    call symputx('ndups', n);
    stop;
run;

%put NOTE: Duplicate IDs found = &ndups.;

/* 9. Quick visual — distribution of key lab values */
proc univariate data=pbc.pbc_raw noprint;
    var bili chol albumin protime;
    histogram / normal kernel;
    inset n nmiss mean std / position=ne;
    title "PBC — Lab Value Distributions";
run;

title;
