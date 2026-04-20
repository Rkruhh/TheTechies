/* 
   Use Case 6 - Automated Report Generation (ODS)
 */


/* 
   PART A - GLOBAL SETUP
*/

%let proj = /home/u64462473/SPL_Project;
libname pbc "&proj.";

ods graphics on / width=6in height=4in attrpriority=color;


%let rpt_date = %sysfunc(today(), worddate.);

proc format;
   value trt_fmt
      1 = "D-Penicillamine"
      2 = "Placebo";

   value stage_fmt
      1 = "Stage 1"
      2 = "Stage 2"
      3 = "Stage 3"
      4 = "Stage 4";

   value sex_fmt
      0 = "Male"
      1 = "Female";

   value yn_fmt
      0 = "No"
      1 = "Yes";

   value edema_fmt
      0   = "None"
      0.5 = "Treated"
      1   = "Untreated";
run;

/* 
   PART B - MACRO: %demographics
*/

%macro demographics(data=);

   data _demog;
      set &data.;
      format trt   trt_fmt.
             stage stage_fmt.
             female sex_fmt.;
      label trt     = "Treatment Group"
            age     = "Age (Years)"
            bili    = "Serum Bilirubin (mg/dL)"
            albumin = "Serum Albumin (g/dL)"
            copper  = "Urine Copper (ug/day)"
            protime = "Prothrombin Time (Seconds)"
            female  = "Sex"
            stage   = "Histologic Stage";
   run;

   title2 "Table 1.1 - Baseline Continuous Variables by Treatment Group";
   title3 "Mean (SD) and Median shown for each variable";

   proc tabulate data=_demog noseps;
      class  trt;
      var    age bili albumin copper protime;
      table  age bili albumin copper protime,
             trt * (mean*f=8.2 std*f=8.2 median*f=8.2)
             / box="Variable" rts=30;
      keylabel mean   = "Mean"
               std    = "SD"
               median = "Median";
      format trt trt_fmt.;
   run;

   title2;

   title2 "Table 1.2 - Baseline Categorical Variables by Treatment Group";
   title3 "n (column %) within each treatment group";

   proc freq data=_demog;
      tables female * trt / nocum norow nopercent
                            plots=none;
      format trt    trt_fmt.
             female sex_fmt.;
      label female = "Sex";
   run;

   proc freq data=_demog;
      tables stage * trt / nocum norow nopercent
                           plots=none;
      format trt   trt_fmt.
             stage stage_fmt.;
      label stage = "Histologic Stage";
   run;

   title2;
   title3;

   proc datasets lib=work nolist;
      delete _demog;
   quit;

%mend demographics;


/* 
   PART C - MACRO: %run_anova
*/

%macro run_anova(data=);

   %let avars = bili albumin chol copper alk_phos ast protime;
   %let nvars = 7;
   %let varlabels = Bilirubin|Albumin|Cholesterol|Copper|
                    Alk Phosphatase|AST (SGOT)|Prothrombin Time;

   %macro one_anova(var=, lab=, grp=, out=);
      ods select none;
      ods output OverallANOVA = _anova_&out.;
      proc glm data=&data.;
         class &grp.;
         model &var. = &grp.;
      run; quit;
      ods select all;

      data _anova_&out.;
         set _anova_&out.;
         where source = "Model";
         variable = "&lab.";
         groupvar = "&grp.";
         keep variable groupvar fvalue probf;
      run;
   %mend one_anova;

   %one_anova(var=bili,     lab=Bilirubin (mg/dL),        grp=trt, out=t1);
   %one_anova(var=albumin,  lab=Albumin (g/dL),           grp=trt, out=t2);
   %one_anova(var=chol,     lab=Cholesterol (mg/dL),      grp=trt, out=t3);
   %one_anova(var=copper,   lab=Copper (ug/day),          grp=trt, out=t4);
   %one_anova(var=alk_phos, lab=Alkaline Phosphatase,     grp=trt, out=t5);
   %one_anova(var=ast,      lab=AST SGOT (U/mL),          grp=trt, out=t6);
   %one_anova(var=protime,  lab=Prothrombin Time (sec),   grp=trt, out=t7);

   %one_anova(var=bili,     lab=Bilirubin (mg/dL),        grp=stage, out=s1);
   %one_anova(var=albumin,  lab=Albumin (g/dL),           grp=stage, out=s2);
   %one_anova(var=chol,     lab=Cholesterol (mg/dL),      grp=stage, out=s3);
   %one_anova(var=copper,   lab=Copper (ug/day),          grp=stage, out=s4);
   %one_anova(var=alk_phos, lab=Alkaline Phosphatase,     grp=stage, out=s5);
   %one_anova(var=ast,      lab=AST SGOT (U/mL),          grp=stage, out=s6);
   %one_anova(var=protime,  lab=Prothrombin Time (sec),   grp=stage, out=s7);

   proc sql noprint;
      create table _means_trt as
      select
         mean(bili)     as bili_1     label="Bili",
         mean(albumin)  as albumin_1  label="Albumin",
         mean(chol)     as chol_1     label="Cholesterol",
         mean(copper)   as copper_1   label="Copper",
         mean(alk_phos) as alkp_1     label="Alk Phos",
         mean(ast)      as ast_1      label="AST",
         mean(protime)  as pt_1       label="Prothrombin"
      from &data. where trt = 1;

      create table _means_plac as
      select
         mean(bili)     as bili_2,
         mean(albumin)  as albumin_2,
         mean(chol)     as chol_2,
         mean(copper)   as copper_2,
         mean(alk_phos) as alkp_2,
         mean(ast)      as ast_2,
         mean(protime)  as pt_2
      from &data. where trt = 2;
   quit;

   data _trt_anova;
      set _anova_t1 _anova_t2 _anova_t3 _anova_t4
          _anova_t5 _anova_t6 _anova_t7;
      sortord = _n_;
   run;

   data _means_trt_long;
      set _means_trt;
      array mvals[7] bili_1 albumin_1 chol_1 copper_1 alkp_1 ast_1 pt_1;
      do i = 1 to 7;
         mean_drug = mvals[i];
         sortord   = i;
         output;
      end;
      keep mean_drug sortord;
   run;

   data _means_plac_long;
      set _means_plac;
      array mvals[7] bili_2 albumin_2 chol_2 copper_2 alkp_2 ast_2 pt_2;
      do i = 1 to 7;
         mean_plac = mvals[i];
         sortord   = i;
         output;
      end;
      keep mean_plac sortord;
   run;

   data _table21;
      merge _trt_anova _means_trt_long _means_plac_long;
      by sortord;
      if probf < 0.001 then pval_fmt = "<0.001";
      else pval_fmt = put(probf, 6.4);
      label variable  = "Lab Variable"
            mean_drug = "Mean - D-Penicillamine"
            mean_plac = "Mean - Placebo"
            fvalue    = "F-Statistic"
            pval_fmt  = "p-value";
      format mean_drug mean_plac 8.2
             fvalue 8.2;
   run;

   title2 "Table 2.1 - Lab Values: D-Penicillamine vs Placebo (One-Way ANOVA)";

   proc report data=_table21 nowd
               style(header)=[background=CX003366 color=white
                               fontweight=bold fontsize=9pt]
               style(column)=[fontsize=9pt]
               style(report) =[borderwidth=2];
      column variable mean_drug mean_plac fvalue pval_fmt;
      define variable  / display "Lab Variable"           style=[width=28%];
      define mean_drug / display "Mean (Drug)"            style=[width=17% just=center];
      define mean_plac / display "Mean (Placebo)"         style=[width=17% just=center];
      define fvalue    / display "F-Statistic"            style=[width=17% just=center];
      define pval_fmt  / display "p-value"                style=[width=14% just=center];
   run;

   title2;

   proc sql noprint;
      create table _stg_means as
      select
         "Bilirubin (mg/dL)"      as variable length=30, 1 as sortord,
         mean(case when stage=1 then bili else . end) as s1_mean,
         mean(case when stage=2 then bili else . end) as s2_mean,
         mean(case when stage=3 then bili else . end) as s3_mean,
         mean(case when stage=4 then bili else . end) as s4_mean
      from &data.
      union all
      select "Albumin (g/dL)",   2,
         mean(case when stage=1 then albumin else . end),
         mean(case when stage=2 then albumin else . end),
         mean(case when stage=3 then albumin else . end),
         mean(case when stage=4 then albumin else . end)
      from &data.
      union all
      select "Cholesterol (mg/dL)", 3,
         mean(case when stage=1 then chol else . end),
         mean(case when stage=2 then chol else . end),
         mean(case when stage=3 then chol else . end),
         mean(case when stage=4 then chol else . end)
      from &data.
      union all
      select "Copper (ug/day)", 4,
         mean(case when stage=1 then copper else . end),
         mean(case when stage=2 then copper else . end),
         mean(case when stage=3 then copper else . end),
         mean(case when stage=4 then copper else . end)
      from &data.
      union all
      select "Alkaline Phosphatase", 5,
         mean(case when stage=1 then alk_phos else . end),
         mean(case when stage=2 then alk_phos else . end),
         mean(case when stage=3 then alk_phos else . end),
         mean(case when stage=4 then alk_phos else . end)
      from &data.
      union all
      select "AST SGOT (U/mL)", 6,
         mean(case when stage=1 then ast else . end),
         mean(case when stage=2 then ast else . end),
         mean(case when stage=3 then ast else . end),
         mean(case when stage=4 then ast else . end)
      from &data.
      union all
      select "Prothrombin Time (sec)", 7,
         mean(case when stage=1 then protime else . end),
         mean(case when stage=2 then protime else . end),
         mean(case when stage=3 then protime else . end),
         mean(case when stage=4 then protime else . end)
      from &data.;
   quit;

   data _stage_anova;
      set _anova_s1 _anova_s2 _anova_s3 _anova_s4
          _anova_s5 _anova_s6 _anova_s7;
      sortord = _n_;
   run;

   data _table22;
      merge _stg_means (in=a) _stage_anova (keep=fvalue probf sortord);
      by sortord;
      if a;
      if probf < 0.001 then pval_fmt = "<0.001";
      else pval_fmt = put(probf, 6.4);
      label variable = "Lab Variable"
            s1_mean  = "Stage 1 Mean"
            s2_mean  = "Stage 2 Mean"
            s3_mean  = "Stage 3 Mean"
            s4_mean  = "Stage 4 Mean"
            fvalue   = "F-Statistic"
            pval_fmt = "p-value";
      format s1_mean s2_mean s3_mean s4_mean fvalue 8.2;
   run;

   title2 "Table 2.2 - Lab Values by Histologic Stage (One-Way ANOVA)";

   proc report data=_table22 nowd
               style(header)=[background=CX003366 color=white
                               fontweight=bold fontsize=9pt]
               style(column)=[fontsize=9pt]
               style(report) =[borderwidth=2];
      column variable s1_mean s2_mean s3_mean s4_mean fvalue pval_fmt;
      define variable  / display "Variable"     style=[width=24%];
      define s1_mean   / display "Stage 1 Mean" style=[width=11% just=center];
      define s2_mean   / display "Stage 2 Mean" style=[width=11% just=center];
      define s3_mean   / display "Stage 3 Mean" style=[width=11% just=center];
      define s4_mean   / display "Stage 4 Mean" style=[width=11% just=center];
      define fvalue    / display "F-Stat"       style=[width=11% just=center];
      define pval_fmt  / display "p-value"      style=[width=11% just=center];
   run;

   title2;

   title2 "Table 2.3 - Tukey Pairwise Comparisons (Significant Variables Only, p < 0.05)";
   title3 "Stage groups compared where overall ANOVA p-value < 0.05";

   %macro tukey_if_sig(var=, lab=, fout=);
      data _null_;
         set _anova_&fout.;
         call symputx("pcheck", probf);
      run;
      %if %sysevalf(&pcheck. < 0.05) %then %do;
         title3 "Tukey Groupings: &lab. by Stage (overall p = %sysfunc(putn(&pcheck., 6.4)))";
         title4 "Groups sharing a letter are not significantly different (alpha=0.05)";
         ods select MCLines;
         proc glm data=&data.;
            class stage;
            model &var. = stage;
            means stage / tukey;
            format stage stage_fmt.;
         run; quit;
         ods select all;
      %end;
      title3;
      title4;
   %mend tukey_if_sig;

   %tukey_if_sig(var=bili,     lab=Bilirubin,         fout=s1);
   %tukey_if_sig(var=albumin,  lab=Albumin,           fout=s2);
   %tukey_if_sig(var=chol,     lab=Cholesterol,       fout=s3);
   %tukey_if_sig(var=copper,   lab=Copper,            fout=s4);
   %tukey_if_sig(var=alk_phos, lab=Alkaline Phosph.,  fout=s5);
   %tukey_if_sig(var=ast,      lab=AST SGOT,          fout=s6);
   %tukey_if_sig(var=protime,  lab=Prothrombin Time,  fout=s7);

   title2;

   /* Clean up */
   proc datasets lib=work nolist;
      delete _anova_t1-_anova_t7 _anova_s1-_anova_s7
             _means_trt _means_plac _means_trt_long _means_plac_long
             _trt_anova _table21 _stg_means _stage_anova _table22;
   quit;

%mend run_anova;


/* 
   PART D - MACRO: %run_survival
*/

%macro run_survival(data=);

   title2 "Figure 3.1 - Kaplan-Meier Survival Curve by Treatment Group";
   title3 "Outcome: All-Cause Mortality | Log-Rank Test";

   ods select SurvivalPlot HomTests;
   proc lifetest data=&data.
                 plots=survival(atrisk test cl nocensor)
                 method=km;
      time   time * death(0);
      strata trt / test=logrank;
      format trt trt_fmt.;
      label  time = "Follow-up Time (Days)";
   run;
   ods select all;

   title2;
   title3;

   title2 "Figure 3.2 - Kaplan-Meier Survival Curve by Histologic Stage";
   title3 "Outcome: All-Cause Mortality | Log-Rank Test";

   ods select SurvivalPlot HomTests;
   proc lifetest data=&data.
                 plots=survival(atrisk test cl nocensor)
                 method=km;
      time   time * death(0);
      strata stage / test=logrank;
      format stage stage_fmt.;
      label  time  = "Follow-up Time (Days)"
             stage = "Histologic Stage";
   run;
   ods select all;

   title2;
   title3;

   title2 "Table 3.1 - Cox Regression: Hazard Ratios for All-Cause Mortality";
   title3 "Covariates: Treatment, Age, Sex, Bilirubin, Albumin, Copper, Prothrombin Time, Stage, Edema";

   ods select ParameterEstimates;
   proc phreg data=&data.;
      class trt   (ref='1')
            stage (ref='1');
      model time * death(0) = trt age female bili albumin copper
                              protime stage edema
                              / ties=efron rl;
      label trt     = "Treatment (Ref: D-Penicillamine)"
            age     = "Age (Years)"
            female  = "Female Sex (Ref: Male)"
            bili    = "Serum Bilirubin (mg/dL)"
            albumin = "Serum Albumin (g/dL)"
            copper  = "Urine Copper (ug/day)"
            protime = "Prothrombin Time (Seconds)"
            stage   = "Histologic Stage (Ref: Stage 1)"
            edema   = "Edema (0/0.5/1)";
   run;
   ods select all;

   title2;
   title3;

%mend run_survival;


/* 
   PART E - MACRO: %run_safety
*/

%macro run_safety(data=);

   data _safety;
      set &data.;
      high_bili    = (bili    >= 3.0);
      high_copper  = (copper  >= 140);
      high_protime = (protime >= 12);
      low_platelet = (platelet < 150);
      edema_any    = (edema    > 0);
   run;

   proc sql noprint;
      select count(*) into :n_drug trimmed from _safety where trt = 1;
      select count(*) into :n_plac trimmed from _safety where trt = 2;
   quit;

   %macro safety_row(var=, lab=, sortord=);
      proc sql noprint;
         select count(*) into :cnt1 trimmed from _safety
            where &var. = 1 and trt = 1;
         select count(*) into :cnt2 trimmed from _safety
            where &var. = 1 and trt = 2;
      quit;

      %let pct1 = %sysfunc(putn(%sysevalf(&cnt1./&n_drug.*100), 5.1));
      %let pct2 = %sysfunc(putn(%sysevalf(&cnt2./&n_plac.*100), 5.1));

      ods select none;
      ods output ChiSq = _cs_&var.;
      proc freq data=_safety;
         tables trt * &var. / chisq;
      run;
      ods select all;

      data _row_&var.;
         set _cs_&var.;
         where statistic = "Chi-Square";
         indicator = "&lab.";
         sortord   = &sortord.;
         drug_np   = cats("&cnt1."," ","(","&pct1.","%",")");
         plac_np   = cats("&cnt2."," ","(","&pct2.","%",")");
         if prob < 0.001 then pval_fmt = "<0.001";
         else pval_fmt = put(prob, 6.4);
         keep indicator sortord drug_np plac_np pval_fmt;
      run;
   %mend safety_row;

   %safety_row(var=ascites,      lab=Ascites,                          sortord=1);
   %safety_row(var=hepato,       lab=Hepatomegaly,                     sortord=2);
   %safety_row(var=spiders,      lab=Spider Angiomata,                 sortord=3);
   %safety_row(var=edema_any,    lab=Any Edema,                        sortord=4);
   %safety_row(var=high_bili,    lab=High Bilirubin (>= 3.0 mg/dL),   sortord=5);
   %safety_row(var=low_albumin,  lab=Low Albumin (< 3.5 g/dL),        sortord=6);
   %safety_row(var=death,        lab=Death (All-Cause),                sortord=7);

   data _safety_table;
      set _row_ascites _row_hepato _row_spiders _row_edema_any
          _row_high_bili _row_low_albumin _row_death;
      label indicator = "Safety Indicator"
            drug_np   = "D-Penicillamine n (%)"
            plac_np   = "Placebo n (%)"
            pval_fmt  = "p-value";
   run;

   title2 "Table 4.1 - Consolidated Safety Summary";
   title3 "n (%) = proportion of patients in each treatment group with the condition";
   footnote1 "Chi-Square p-value shown; Fisher exact recommended where expected cell count < 5";
   footnote2 "Total: D-Penicillamine n = &n_drug. | Placebo n = &n_plac.";

   proc report data=_safety_table nowd
               style(header)=[background=CX003366 color=white
                               fontweight=bold fontsize=9pt]
               style(column)=[fontsize=9pt]
               style(report) =[borderwidth=2];
      column indicator drug_np plac_np pval_fmt;
      define indicator / display "Safety Indicator"       style=[width=36%];
      define drug_np   / display "D-Penicillamine n (%)"  style=[width=22% just=center];
      define plac_np   / display "Placebo n (%)"          style=[width=22% just=center];
      define pval_fmt  / display "p-value"                style=[width=13% just=center];
   run;

   title2;
   title3;
   footnote;

   proc sql;
      create table _ae_plot as
      select trt,
         "Ascites"           as indicator length=35, 1 as sortord,
         sum(ascites)/count(*)*100 as rate
      from _safety group by trt
      union all
      select trt, "Hepatomegaly",   2, sum(hepato)/count(*)*100    from _safety group by trt
      union all
      select trt, "Spider Angiomata", 3, sum(spiders)/count(*)*100 from _safety group by trt
      union all
      select trt, "Any Edema",      4, sum(edema_any)/count(*)*100 from _safety group by trt
      union all
      select trt, "High Bilirubin", 5, sum(high_bili)/count(*)*100 from _safety group by trt
      union all
      select trt, "Low Albumin",    6, sum(low_albumin)/count(*)*100 from _safety group by trt
      union all
      select trt, "Death",          7, sum(death)/count(*)*100     from _safety group by trt;
   quit;

   data _ae_plot;
      set _ae_plot;
      if trt = 1 then trt_lbl = "D-Penicillamine";
      else           trt_lbl = "Placebo";
   run;

   proc sort data=_ae_plot; by sortord; run;

   title2 "Figure 4.1 - Adverse Event Rates by Treatment Group";
   title3 "Percentage of patients with each safety indicator";
   footnote1 "D-Penicillamine (n=&n_drug.) and Placebo (n=&n_plac.)";

   proc sgplot data=_ae_plot;

      format rate 5.1;

      vbar indicator / response=rate
                       group=trt_lbl
                       groupdisplay=cluster
                       datalabel
                       datalabelattrs=(size=7)
                       fillattrs=(transparency=0.1);
      xaxis label="Safety Indicator"
            discreteorder=data
            fitpolicy=rotate
            valueattrs=(size=8);
      yaxis label="Event Rate (%)" grid max=80;
      keylegend / title="Treatment Group" location=inside position=topright;
      styleattrs datacolors=(steelblue tomato);
   run;

   title2;
   title3;
   footnote;

   proc datasets lib=work nolist;
      delete _safety _ae_plot _safety_table
             _cs_ascites _cs_hepato _cs_spiders _cs_edema_any
             _cs_high_bili _cs_low_albumin _cs_death
             _row_ascites _row_hepato _row_spiders _row_edema_any
             _row_high_bili _row_low_albumin _row_death;
   quit;

%mend run_safety;


/* 
   PART F - MACRO: %generate_report
*/

%macro generate_report(data=, outfile=);

   /* 
      Open ODS PDF destination
   */
   ods pdf file="&outfile."
           style=htmlblue
           pdftoc=2
           bookmarkgen=yes
           notoc;

   ods pdf text="^S={just=c font_size=18pt font_weight=bold
                      topmargin=3in}
                 PBC Clinical Trial Statistical Report";

   ods pdf text="^S={just=c font_size=13pt topmargin=0.4in}
                 D-Penicillamine vs Placebo for Primary Biliary Cirrhosis";

   ods pdf text="^S={just=c font_size=10pt topmargin=0.3in color=gray}
                 Report Date: &rpt_date.";

   ods pdf text="^S={just=c font_size=10pt topmargin=0.2in}
                 Rucha Keni | Satvik Nayak | Likhita Sri Kode | Tye Jamian";

   ods pdf text="^S={just=c font_size=8pt topmargin=0.3in color=gray}
                 Dataset: PBC Clinical Trial | N = 312 Trial Patients";

   ods pdf startpage=now;


   /* 
      SECTION 1 - DEMOGRAPHICS
   */
   ods pdf startpage=now;

   title1 height=14pt bold "Section 1 - Patient Demographics";
   title2 height=11pt "Baseline Characteristics by Treatment Group";

   %demographics(data=&data.);

   title;
   ods pdf startpage=now;


   /* 
      SECTION 2 - ANOVA
   */
   ods pdf startpage=now;

   title1 height=14pt bold "Section 2 - Analysis of Variance (ANOVA)";
   title2 height=11pt "Lab Value Comparisons by Treatment Group and Histologic Stage";

   %run_anova(data=&data.);

   title;
   ods pdf startpage=now;

   /* 
      SECTION 3 - SURVIVAL ANALYSIS
   */
   ods pdf startpage=now;

   title1 height=14pt bold "Section 3 - Survival Analysis";
   title2 height=11pt "Kaplan-Meier Curves and Cox Proportional Hazards Model";

   %run_survival(data=&data.);

   title;
   ods pdf startpage=now;


   /* 
      SECTION 4 - SAFETY SUMMARY
   */
   ods pdf startpage=now;

   title1 height=14pt bold "Section 4 - Safety and Adverse Event Summary";
   title2 height=11pt "Clinical Indicators and Lab Abnormalities by Treatment Group";

   %run_safety(data=&data.);

   title;

   /* 
      Close ODS PDF
   */
   ods pdf close;

   %put NOTE: Report written to -- &outfile.;

%mend generate_report;

/* 
   PART G - MASTER EXECUTION
*/

%generate_report(
   data    = pbc.pbc_trial,
   outfile = /home/u64462473/SPL_Project/PBC_Clinical_Trial_Report.pdf
);

ods graphics off;
