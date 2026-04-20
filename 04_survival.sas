/*
   Purpose  : Use Case 4 - Survival Analysis
*/

%let proj = /home/u64462473/SPL_Project;
libname pbc "&proj.";

ods graphics on;


/* 
   SECTION 1 - KAPLAN-MEIER SURVIVAL ANALYSIS
*/

/* KM Curves by Treatment Group (Drug vs Placebo) */
title1 "Kaplan-Meier Survival Analysis";
title2 "Survival by Treatment Group: D-Penicillamine vs Placebo";

proc lifetest data=pbc.pbc_trial
              plots=survival(atrisk test cl)
              method=km;
   time   time * death(0);
   strata trt / test=logrank;
   label  time  = "Follow-up Time (Days)"
          death = "Death Indicator (1=Died)";
run;

title;


/* KM Curves by Histologic Stage (4 groups) */
title1 "Kaplan-Meier Survival Analysis";
title2 "Survival by Histologic Stage (1 to 4)";

proc lifetest data=pbc.pbc_trial
              plots=survival(atrisk test cl)
              method=km;
   time   time * death(0);
   strata stage / test=logrank;
   label  time  = "Follow-up Time (Days)"
          stage = "Histologic Stage";
run;

title;


/* KM Curves by Bilirubin Severity Category (bili_cat) */
title1 "Kaplan-Meier Survival Analysis";
title2 "Survival by Bilirubin Severity Category";
title3 "(1=Normal to 4=Very High)";

proc lifetest data=pbc.pbc_trial
              plots=survival(atrisk test cl)
              method=km;
   time   time * death(0);
   strata bili_cat / test=logrank;
   label  time     = "Follow-up Time (Days)"
          bili_cat = "Bilirubin Severity Category";
run;

title;


/* 
   SECTION 2 - COX PROPORTIONAL HAZARDS REGRESSION
*/

/* 2A: Multivariate Cox Model */
title1 "Cox Proportional Hazards Regression";
title2 "Multivariate Model - All Trial Patients";

proc phreg data=pbc.pbc_trial;
   class trt   (ref='1')   
         stage (ref='1');  
   model time * death(0) = trt age female bili albumin copper
                           protime stage edema
                           / ties=efron rl;
   hazardratio trt   / cl=wald;
   hazardratio stage / cl=wald;
   label trt     = "Treatment (1=Drug, 2=Placebo)"
         age     = "Age (Years)"
         female  = "Sex (1=Female, 0=Male)"
         bili    = "Serum Bilirubin (mg/dL)"
         albumin = "Serum Albumin (g/dL)"
         copper  = "Urine Copper (ug/day)"
         protime = "Prothrombin Time (Seconds)"
         stage   = "Histologic Stage"
         edema   = "Edema (0/0.5/1)";
run;

title;


/* Assess Proportional Hazards Assumption */
title1 "Cox Model - Proportional Hazards Assumption";
title2 "Assessment via Cumulative Martingale Residuals";

proc phreg data=pbc.pbc_trial;
   class trt   (ref='1')
         stage (ref='1');
   model time * death(0) = trt age female bili albumin copper
                           protime stage edema
                           / ties=efron;
   assess ph / resample seed=12345;
run;

title;


/* 
   SECTION 3 - SUBGROUP COX MODELS
*/

/* Cox Model by Histologic Stage */

%macro cox_by_stage(stg=);
   title1 "Subgroup Cox Regression - Stage and stg.";
   proc phreg data=pbc.pbc_trial;
      where stage = &stg.;
      class trt (ref='1');
      model time * death(0) = trt age female bili albumin copper
                              protime edema
                              / ties=efron rl;
      hazardratio trt / cl=wald;
   run;
   title;
%mend cox_by_stage;

%cox_by_stage(stg=1);
%cox_by_stage(stg=2);
%cox_by_stage(stg=3);
%cox_by_stage(stg=4);


/* --- 3B: Cox Model by Sex --- */

title1 "Subgroup Cox Regression - Females Only";
proc phreg data=pbc.pbc_trial;
   where female = 1;
   class trt   (ref='1')
         stage (ref='1');
   model time * death(0) = trt age bili albumin copper
                           protime stage edema
                           / ties=efron rl;
   hazardratio trt / cl=wald;
run;
title;

title1 "Subgroup Cox Regression - Males Only";
proc phreg data=pbc.pbc_trial;
   where female = 0;
   class trt   (ref='1')
         stage (ref='1');
   model time * death(0) = trt age bili albumin copper
                           protime stage edema
                           / ties=efron rl;
   hazardratio trt / cl=wald;
run;
title;


/* 
   SECTION 4 - FOREST PLOT DATA PREPARATION and VISUALIZATION
*/

/* Build the forest plot summary dataset  */
data forest_data;
   length subgroup $30 label $40;

   subgroup = "Stage 1";
   label    = "Stage 1 (Mildest)";
   hr  = 1.00; lcl = 0.50; ucl = 2.00;  
   output;

   subgroup = "Stage 2";
   label    = "Stage 2";
   hr  = 1.00; lcl = 0.50; ucl = 2.00;  
   output;

   subgroup = "Stage 3";
   label    = "Stage 3";
   hr  = 1.00; lcl = 0.50; ucl = 2.00;   
   output;

   subgroup = "Stage 4";
   label    = "Stage 4 (Most Severe)";
   hr  = 1.00; lcl = 0.50; ucl = 2.00;  
   output;

   /* Subgroups by sex */
   subgroup = "Female";
   label    = "Female";
   hr  = 1.00; lcl = 0.50; ucl = 2.00;  
   output;

   subgroup = "Male";
   label    = "Male";
   hr  = 1.00; lcl = 0.50; ucl = 2.00;   
   output;

   /* Overall model (full cohort) */
   subgroup = "Overall";
   label    = "Overall";
   hr  = 1.00; lcl = 0.50; ucl = 2.00;  
   output;

   /* Log-transform CIs for display on natural scale */
   log_hr  = log(hr);
   log_lcl = log(lcl);
   log_ucl = log(ucl);

run;

data forest_data;
   set forest_data;
   select (subgroup);
      when ("Overall") sortord = 1;
      when ("Female")  sortord = 2;
      when ("Male")    sortord = 3;
      when ("Stage 1") sortord = 4;
      when ("Stage 2") sortord = 5;
      when ("Stage 3") sortord = 6;
      when ("Stage 4") sortord = 7;
      otherwise        sortord = 99;
   end;
run;

proc sort data=forest_data; by sortord; run;


/* Render the forest plot */
title1 "Forest Plot - Treatment Hazard Ratios by Subgroup";
title2 "Outcome: All-Cause Mortality | Reference: D-Penicillamine (trt=1)";
title3 "(Update HR/LCL/UCL values from HAZARDRATIO output before presenting)";

proc sgplot data=forest_data noautolegend;

   refline 1 / axis=x lineattrs=(color=gray pattern=shortdash);
   highlow y=label low=lcl high=ucl / type=line
           lineattrs=(color=steelblue thickness=2)
           clipcap;
   scatter y=label x=hr / markerattrs=(symbol=squarefilled
           color=steelblue size=10);

   xaxis label="Hazard Ratio (Placebo vs D-Penicillamine)"
         logbase=e logstyle=linear
         values=(0.25 0.5 1 2 4)
         valuesdisplay=("0.25" "0.50" "1.00" "2.00" "4.00");
   yaxis label="Subgroup" discreteorder=data reverse;

   inset "HR greater than 1 favors D-Penicillamine" / position=topright
         textattrs=(size=8);
run;

title;

ods graphics off;
