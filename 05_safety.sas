/* 
   Use Case 5 - Safety and Adverse Event Analysis
*/

proc format;
   value trt_fmt
      1 = "D-Penicillamine"
      2 = "Placebo";
   value yn_fmt
      0 = "No"
      1 = "Yes";
run;

/*
   SECTION 0 - DATA PREPARATION
*/

data safety;
   set pbc.pbc_trial;

   high_bili    = (bili    >= 3.0);   
   high_copper  = (copper  >= 140);   
   high_protime = (protime >= 12);  
   low_platelet = (platelet < 150);  

   edema_any = (edema > 0);

   if trt = 1 then trt_label = "1 - D-Penicillamine";
   else if trt = 2 then trt_label = "2 - Placebo";

   stage_label = cats("Stage ", put(stage, 1.));

   label
      high_bili    = "High Bilirubin (>= 3.0 mg/dL)"
      high_copper  = "High Copper (>= 140 ug/day)"
      high_protime = "High Prothrombin Time (>= 12 sec)"
      low_platelet = "Low Platelet Count (< 150)"
      low_albumin  = "Low Albumin (< 3.5 g/dL)"
      edema_any    = "Any Edema (Edema > 0)"
      trt_label    = "Treatment Group"
      ascites      = "Ascites (Fluid in Abdomen)"
      hepato       = "Hepatomegaly (Enlarged Liver)"
      spiders      = "Spider Angiomata";
run;


/* 
   SECTION 1 - ADVERSE CLINICAL INDICATOR RATES BY TREATMENT
*/

title1 "Section 1 - Adverse Clinical Indicators by Treatment Group";
title2 "Ascites (Fluid in Abdomen) - Drug vs Placebo";

proc freq data=safety;
   tables trt * ascites / chisq relrisk riskdiff
                          norow nocol nopercent;
   format trt trt_fmt.;
run;

title;

title1 "Section 1 - Adverse Clinical Indicators by Treatment Group";
title2 "Hepatomegaly (Enlarged Liver) - Drug vs Placebo";

proc freq data=safety;
   tables trt * hepato / chisq relrisk riskdiff
                         norow nocol nopercent;
   format trt trt_fmt.;
run;

title;

title1 "Section 1 - Adverse Clinical Indicators by Treatment Group";
title2 "Spider Angiomata - Drug vs Placebo";

proc freq data=safety;
   tables trt * spiders / chisq relrisk riskdiff
                          norow nocol nopercent;
   format trt trt_fmt.;
run;

title;

title1 "Section 1 - Adverse Clinical Indicators by Treatment Group";
title2 "Any Edema (Edema > 0) - Drug vs Placebo";

proc freq data=safety;
   tables trt * edema_any / chisq relrisk riskdiff
                            norow nocol nopercent;
   format trt trt_fmt.;
run;

title;


/* 
   SECTION 2 - DEATH RATE BY TREATMENT GROUP
*/

title1 "Section 2 - Mortality by Treatment Group";
title2 "All-Cause Death Rate - D-Penicillamine vs Placebo";

proc freq data=safety;
   tables trt * death / chisq relrisk riskdiff
                        norow nocol nopercent;
   format trt   trt_fmt.
          death yn_fmt.;
   label death = "Died (1=Yes, 0=No)";
run;

title;

title1 "Section 2 - Mortality by Treatment Group and Stage";
title2 "Death Rate within Each Histologic Stage - Drug vs Placebo";

proc freq data=safety;
   tables trt * stage * death / chisq norow nocol nopercent;
   format trt   trt_fmt.
          death yn_fmt.;
run;

title;


/* 
   SECTION 3 - ABNORMAL LAB VALUE RATES BY TREATMENT GROUP
*/

%macro lab_freq(var=, varlabel=, sec=);
   title1 "Section 3 - Abnormal Lab Values by Treatment Group";
   title2 "&sec. - Drug vs Placebo";

   proc freq data=safety;
      tables trt * &var. / chisq relrisk riskdiff
                           norow nocol nopercent;
      format trt  trt_fmt.
             &var. yn_fmt.;
      label &var. = "&varlabel.";
   run;

   title;
%mend lab_freq;

%lab_freq(var=high_bili,    varlabel=High Bilirubin (>= 3.0 mg/dL),    sec=High Bilirubin);
%lab_freq(var=low_albumin,  varlabel=Low Albumin (< 3.5 g/dL),         sec=Low Albumin);
%lab_freq(var=high_copper,  varlabel=High Copper (>= 140 ug/day),       sec=High Copper);
%lab_freq(var=high_protime, varlabel=High Prothrombin Time (>= 12 sec), sec=High Prothrombin Time);
%lab_freq(var=low_platelet, varlabel=Low Platelet Count (< 150),        sec=Low Platelet Count);


/* 
   SECTION 4 - CONSOLIDATED SAFETY SUMMARY TABLE
   Approach:
   Step 1 - Use PROC FREQ with ODS OUTPUT to capture counts,
            percentages, and p-values for every indicator.
   Step 2 - Stack results into one summary dataset.
   Step 3 - Render with PROC REPORT.
*/

%macro get_pval(var=, indicator=);
   ods exclude all;
   ods output ChiSq        = chisq_&var.
              CrossTabFreqs = freq_&var.;

   proc freq data=safety;
      tables trt * &var. / chisq;
      format trt trt_fmt. &var. yn_fmt.;
   run;

   ods select all;

   /* Extract chi-square p-value row */
   data pval_&var.;
      set chisq_&var.;
      where statistic = "Chi-Square";
      indicator = "&indicator.";
      pval = prob;
      keep indicator pval;
   run;

%mend get_pval;

%get_pval(var=ascites,      indicator=Ascites);
%get_pval(var=hepato,       indicator=Hepatomegaly);
%get_pval(var=spiders,      indicator=Spider Angiomata);
%get_pval(var=edema_any,    indicator=Any Edema);
%get_pval(var=high_bili,    indicator=High Bilirubin);
%get_pval(var=low_albumin,  indicator=Low Albumin);
%get_pval(var=high_copper,  indicator=High Copper);
%get_pval(var=high_protime, indicator=High Prothrombin Time);
%get_pval(var=low_platelet, indicator=Low Platelet Count);
%get_pval(var=death,        indicator=Death);


/* --- 4B: Stack p-values into one dataset --- */
data all_pvals;
   set pval_ascites
       pval_hepato
       pval_spiders
       pval_edema_any
       pval_high_bili
       pval_low_albumin
       pval_high_copper
       pval_high_protime
       pval_low_platelet
       pval_death;
run;


proc sql noprint;
   select count(*) into :n_drug trimmed from safety where trt = 1;
   select count(*) into :n_plac trimmed from safety where trt = 2;
quit;

%macro get_n(var=, trt_val=, out=);
   proc sql noprint;
      select count(*) into :cnt trimmed
      from safety
      where &var. = 1 and trt = &trt_val.;
   quit;

   %let pct = %sysevalf(&cnt / %sysfunc(ifc(&trt_val.=1, &n_drug., &n_plac.)) * 100, ceil);
   %global &out._n &out._pct;
   %let &out._n   = &cnt.;
   %let &out._pct = %sysfunc(putn(%sysevalf(&cnt / %sysfunc(ifc(&trt_val.=1, &n_drug., &n_plac.)) * 100), 5.1));
%mend get_n;

%get_n(var=ascites,      trt_val=1, out=asc1);
%get_n(var=ascites,      trt_val=2, out=asc2);
%get_n(var=hepato,       trt_val=1, out=hep1);
%get_n(var=hepato,       trt_val=2, out=hep2);
%get_n(var=spiders,      trt_val=1, out=spi1);
%get_n(var=spiders,      trt_val=2, out=spi2);
%get_n(var=edema_any,    trt_val=1, out=ede1);
%get_n(var=edema_any,    trt_val=2, out=ede2);
%get_n(var=high_bili,    trt_val=1, out=bil1);
%get_n(var=high_bili,    trt_val=2, out=bil2);
%get_n(var=low_albumin,  trt_val=1, out=alb1);
%get_n(var=low_albumin,  trt_val=2, out=alb2);
%get_n(var=high_copper,  trt_val=1, out=cop1);
%get_n(var=high_copper,  trt_val=2, out=cop2);
%get_n(var=high_protime, trt_val=1, out=pro1);
%get_n(var=high_protime, trt_val=2, out=pro2);
%get_n(var=low_platelet, trt_val=1, out=plt1);
%get_n(var=low_platelet, trt_val=2, out=plt2);
%get_n(var=death,        trt_val=1, out=dth1);
%get_n(var=death,        trt_val=2, out=dth2);


/* Build the summary dataset row by row */
data safety_summary;
   length indicator $35 drug_np $15 plac_np $15;
   format pval pvalue6.4;

   /* Total n header row */
   sortord  = 0;
   indicator = "Total N";
   drug_np   = cats(&n_drug.);
   plac_np   = cats(&n_plac.);
   pval      = .;
   output;

   /* Clinical indicators */
   sortord=1; indicator="Ascites";
   drug_np=cats("&asc1_n. ","(","&asc1_pct.","%",")");
   plac_np=cats("&asc2_n. ","(","&asc2_pct.","%",")");
   set all_pvals; where indicator="Ascites";
   output;

   sortord=2; indicator="Hepatomegaly";
   drug_np=cats("&hep1_n. ","(","&hep1_pct.","%",")");
   plac_np=cats("&hep2_n. ","(","&hep2_pct.","%",")");
   set all_pvals; where indicator="Hepatomegaly";
   output;

   sortord=3; indicator="Spider Angiomata";
   drug_np=cats("&spi1_n. ","(","&spi1_pct.","%",")");
   plac_np=cats("&spi2_n. ","(","&spi2_pct.","%",")");
   set all_pvals; where indicator="Spider Angiomata";
   output;

   sortord=4; indicator="Any Edema";
   drug_np=cats("&ede1_n. ","(","&ede1_pct.","%",")");
   plac_np=cats("&ede2_n. ","(","&ede2_pct.","%",")");
   set all_pvals; where indicator="Any Edema";
   output;

   /* Lab abnormalities */
   sortord=5; indicator="High Bilirubin (>= 3.0)";
   drug_np=cats("&bil1_n. ","(","&bil1_pct.","%",")");
   plac_np=cats("&bil2_n. ","(","&bil2_pct.","%",")");
   set all_pvals; where indicator="High Bilirubin";
   output;

   sortord=6; indicator="Low Albumin (< 3.5 g/dL)";
   drug_np=cats("&alb1_n. ","(","&alb1_pct.","%",")");
   plac_np=cats("&alb2_n. ","(","&alb2_pct.","%",")");
   set all_pvals; where indicator="Low Albumin";
   output;

   sortord=7; indicator="High Copper (>= 140 ug/day)";
   drug_np=cats("&cop1_n. ","(","&cop1_pct.","%",")");
   plac_np=cats("&cop2_n. ","(","&cop2_pct.","%",")");
   set all_pvals; where indicator="High Copper";
   output;

   sortord=8; indicator="High Prothrombin Time (>= 12)";
   drug_np=cats("&pro1_n. ","(","&pro1_pct.","%",")");
   plac_np=cats("&pro2_n. ","(","&pro2_pct.","%",")");
   set all_pvals; where indicator="High Prothrombin Time";
   output;

   sortord=9; indicator="Low Platelet Count (< 150)";
   drug_np=cats("&plt1_n. ","(","&plt1_pct.","%",")");
   plac_np=cats("&plt2_n. ","(","&plt2_pct.","%",")");
   set all_pvals; where indicator="Low Platelet Count";
   output;

   /* Death */
   sortord=10; indicator="Death (All-Cause)";
   drug_np=cats("&dth1_n. ","(","&dth1_pct.","%",")");
   plac_np=cats("&dth2_n. ","(","&dth2_pct.","%",")");
   set all_pvals; where indicator="Death";
   output;

run;

/*
   NOTE: The safety_summary dataset built above using nested SET
   inside DATA is a structural template. In practice the simplest
   and most robust approach in SAS Studio is to build the summary
   table directly in PROC REPORT reading from the safety dataset.
   The PROC REPORT below reads the safety dataset directly and
   computes everything inline via summary statistics.
   The safety_summary dataset is kept as a reference / export target.
*/


title1 "Section 4 - Consolidated Safety Summary Table";
title2 "Adverse Event and Lab Abnormality Rates by Treatment Group";
title3 "n (Row %) shown for each group | Chi-Square p-value";
footnote1 "Reference: D-Penicillamine (trt=1), Placebo (trt=2)";
footnote2 "Row percentages computed within each treatment group";
footnote3 "* Fisher exact test recommended where expected cell count < 5 (see Section 6)";


data safety_long;
   set safety;
   length indicator $35;

   array flags[10] ascites hepato spiders edema_any
                   high_bili low_albumin high_copper
                   high_protime low_platelet death;

   array indlabels[10] $35 _temporary_ (
      "Ascites"
      "Hepatomegaly"
      "Spider Angiomata"
      "Any Edema"
      "High Bilirubin (>= 3.0 mg/dL)"
      "Low Albumin (< 3.5 g/dL)"
      "High Copper (>= 140 ug/day)"
      "High Prothrombin Time (>= 12)"
      "Low Platelet Count (< 150)"
      "Death (All-Cause)"
   );

   array sortvals[10] _temporary_ (1,2,3,4,5,6,7,8,9,10);

   do i = 1 to 10;
      indicator = indlabels[i];
      event     = flags[i];
      sortord   = sortvals[i];
      output;
   end;

   keep trt indicator event sortord;
run;

proc sort data=safety_long; by sortord indicator trt; run;

/* Summary counts per indicator per treatment */
proc sql;
   create table safety_report as
   select
      indicator,
      sortord,
      trt,
      sum(event)               as n_event,
      count(*)                 as n_total,
      sum(event) / count(*) * 100 as pct_event
   from safety_long
   group by indicator, sortord, trt
   order by sortord, trt;
quit;

/* Sort by indicator (alphabetical) as required by PROC TRANSPOSE BY statement */
proc sort data=safety_report; by indicator sortord trt; run;

/* Pivot to wide format: one row per indicator */
proc transpose data=safety_report out=safety_wide_n    prefix=n_trt;
   by indicator sortord;
   id trt;
   var n_event;
run;
proc transpose data=safety_report out=safety_wide_pct  prefix=pct_trt;
   by indicator sortord;
   id trt;
   var pct_event;
run;

data safety_wide;
   merge safety_wide_n   (keep=indicator sortord n_trt1 n_trt2
                          rename=(n_trt1=n1 n_trt2=n2))
         safety_wide_pct (keep=indicator sortord pct_trt1 pct_trt2
                          rename=(pct_trt1=pct1 pct_trt2=pct2));
   by indicator sortord;
   drug_np = cats(n1, " (", put(pct1, 5.1), "%)");
   plac_np = cats(n2, " (", put(pct2, 5.1), "%)");
run;

/* Merge in p-values */
proc sort data=safety_wide;   by indicator; run;
proc sort data=all_pvals;     by indicator; run;

data safety_final;
   merge safety_wide (in=a)
         all_pvals   (in=b);
   by indicator;
   if a;
   if pval < 0.001 then pval_fmt = "<0.001";
   else pval_fmt = put(pval, 6.4);
   label
      indicator = "Safety Indicator"
      drug_np   = "D-Penicillamine n (%)"
      plac_np   = "Placebo n (%)"
      pval_fmt  = "p-value";
run;

proc sort data=safety_final; by sortord; run;

proc report data=safety_final nowd
            style(header)=[background=CX003366 color=white fontsize=10pt fontweight=bold]
            style(column)=[fontsize=9pt]
            style(report) =[borderwidth=2];
   column indicator drug_np plac_np pval_fmt;
   define indicator / display "Safety Indicator"  style(column)=[width=35%];
   define drug_np   / display "D-Penicillamine n (%)"  style(column)=[width=22% just=center];
   define plac_np   / display "Placebo n (%)"           style(column)=[width=22% just=center];
   define pval_fmt  / display "p-value"                 style(column)=[width=12% just=center];
run;

title;
footnote;


/*
   SECTION 5 - VISUALIZATIONS
*/

proc sql;
   create table ae_plot as
   select
      indicator,
      sortord,
      trt,
      sum(event) / count(*) * 100 as rate label="Event Rate (%)"
   from safety_long
   where sortord <= 9   /* exclude death - shown separately */
   group by indicator, sortord, trt;
quit;


data ae_plot;
   set ae_plot;
   if trt = 1 then trt_lbl = "D-Penicillamine";
   else           trt_lbl = "Placebo";
run;

proc sort data=ae_plot; by sortord; run;

title1 "Section 5 - Adverse Event Rates by Treatment Group";
title2 "Grouped Bar Chart: D-Penicillamine vs Placebo";
footnote1 "Rate = percentage of patients in each treatment group with the condition";

proc sgplot data=ae_plot;
   vbar indicator / response=rate group=trt_lbl
                    groupdisplay=cluster
                    datalabel datalabelattrs=(size=7)
                    fillattrs=(transparency=0.15);
   xaxis label="Safety Indicator"
         discreteorder=data
         fitpolicy=rotate
         valueattrs=(size=7);
   yaxis label="Event Rate (%)" grid;
   keylegend / title="Treatment Group" location=inside position=topright;
   styleattrs datacolors=(steelblue tomato);
run;

title;
footnote;

proc sql;
   create table death_stage as
   select
      trt,
      stage,
      sum(death)          as n_dead,
      count(*)            as n_total,
      sum(death)/count(*)*100 as death_rate label="Death Rate (%)"
   from safety
   group by trt, stage;
quit;

data death_stage;
   set death_stage;
   if trt = 1 then trt_lbl = "D-Penicillamine";
   else           trt_lbl = "Placebo";
   stage_lbl = cats("Stage ", stage);
run;

title1 "Section 5 - Death Rate by Treatment Group and Histologic Stage";
title2 "Stacked Bar Chart: Proportion of Deaths within Each Stage";
footnote1 "Height of each segment = death rate (%) within that stage and treatment group";

proc sgplot data=death_stage;
   vbar trt_lbl / response=death_rate group=stage_lbl
                  groupdisplay=stack
                  datalabel datalabelattrs=(size=7 color=white);
   xaxis label="Treatment Group";
   yaxis label="Death Rate (%)" grid;
   keylegend / title="Histologic Stage" location=inside position=topright;
   styleattrs datacolors=(CX2166AC CX4DAC26 CXFDAE61 CXD7191C);
run;

title;
footnote;


/* --- 5D: Heat-map style matrix - event rate by indicator and treatment --- */
title1 "Section 5 - Safety Indicator Heat Map";
title2 "Event Rate (%) by Indicator and Treatment Group";

proc sgplot data=ae_plot;
   hbar indicator / response=rate group=trt_lbl
                    groupdisplay=cluster
                    datalabel datalabelattrs=(size=7)
                    fillattrs=(transparency=0.1);
   yaxis label="Safety Indicator" discreteorder=data;
   xaxis label="Event Rate (%)" grid;
   keylegend / title="Treatment Group";
   styleattrs datacolors=(steelblue tomato);
run;

title;


/* 
   SECTION 6 - FISHER EXACT TEST FOR RARE EVENTS
*/

title1 "Section 6 - Fisher Exact Test for Rare Events";
title2 "Ascites: Observed and Expected Cell Counts";
footnote1 "If any expected count < 5, use Fisher exact test (p-value from EXACT statement)";

proc freq data=safety;
   tables trt * ascites / chisq expected
                          norow nocol nopercent;
   format trt     trt_fmt.
          ascites yn_fmt.;
run;

title;
footnote;

title1 "Section 6 - Fisher Exact Test: Ascites by Treatment Group";
title2 "Two-Sided Fisher Exact p-value - D-Penicillamine vs Placebo";

proc freq data=safety;
   tables trt * ascites / chisq fisher
                          norow nocol nopercent;
   exact fisher;
   format trt     trt_fmt.
          ascites yn_fmt.;
run;

title;

title1 "Section 6 - Fisher Exact Test: Low Platelet Count by Treatment Group";
title2 "Two-Sided Fisher Exact p-value - D-Penicillamine vs Placebo";

proc freq data=safety;
   tables trt * low_platelet / chisq fisher expected
                               norow nocol nopercent;
   exact fisher;
   format trt         trt_fmt.
          low_platelet yn_fmt.;
run;

title;

ods graphics off;
