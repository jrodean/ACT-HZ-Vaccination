* ==============================================================================
* HTE_all.do
*
* Every heterogeneity-of-treatment-effect (HTE) result in the paper: Figure 4
* panels A-E and the 8 supplementary subgroup-split figures/tables (Figures S30-S37, with Tables S17-S21 for panels A-E). The ONLY
* input is data/analysis_panel.dta, written by R/00c_lasso_risk_score.R at
* the end of the R data-construction pipeline -- so this script, together
* with the R scripts, still replicates everything starting from the raw ACT
* CSV extracts; this file is just the Stata half of that pipeline.
*
* Each panel is produced by a single call to the `hte_panel` program (below),
* which:
*   1. fits two separate-sample event-study regressions (group==1, group==0)
*   2. draws the coefplot in the exact style of the project's own HTE do-files
*      (Xi/v0326/HTE_0926.do, Xi/v0326/HTE_0326.do) -- plotplain scheme,
*      filled circle vs. open triangle, dashed grey CIs, red vaccination-year
*      line, bottom legend box, default (gridded) y-axis -- confirmed against
*      the actual published Figure 4A image embedded in the manuscript
*   3. writes the Table S "difference in event-study treatment effects"
*      companion table (Estimate/CI/P/BH), using the exact two-model method
*      confirmed against cleaned_lilac/cleaned_data/Stata_code/HTE_difference.do:
*      displayed columns from the two separate-sample models, but the
*      pretrend joint Wald p-value from a SEPARATE pooled regression with
*      group x event-time interactions (not from combining the separate
*      models' covariances -- those are not the same test).
*
* Run from the "replication package" directory:
*   stata -b do stata/HTE_all.do
* ==============================================================================

cd "`c(pwd)'"
set scheme plotplain

cap mkdir output
cap mkdir output/figures
cap mkdir output/tables

use "data/analysis_panel.dta", clear
keep if !missing(zost1_date) & !missing(casi_total)

* -------- Event time: 2-year bins, cap [-14,12] (same formula everywhere) --------
cap drop event_time
gen int event_time = floor((visit_date - zost1_date)/(2*365))*2
replace event_time =  12 if event_time >=  12
replace event_time = -14 if event_time < -14

cap drop trt_m* trt_p*
forvalues i = 1/6 {
    gen byte trt_m`i' = (event_time == -`i'*2)
}
forvalues i = 0/6 {
    gen byte trt_p`i' = (event_time ==  `i'*2)
}
replace trt_m1 = 0

label var trt_m6 "12 yrs before"
label var trt_m5 "10 yrs before"
label var trt_m4 "8 yrs before"
label var trt_m3 "6 yrs before"
label var trt_m2 "4 yrs before"
label var trt_m1 "2 yrs before, ref"
label var trt_p0 "ZVL vaccination year"
label var trt_p1 "2 yrs after"
label var trt_p2 "4 yrs after"
label var trt_p3 "6 yrs after"
label var trt_p4 "8 yrs after"
label var trt_p5 "10 yrs after"
label var trt_p6 "12 yrs after"

* ------------------------------------------------------------------------
* hte_panel: one HTE panel (plot + difference table) for a 0/1 splitvar
*   args: splitvar  label1  label0  legendtitle  outname
* ------------------------------------------------------------------------
cap program drop hte_panel
program define hte_panel
    args splitvar label1 label0 legendtitle outname

    quietly count if !missing(`splitvar')
    if r(N) == 0 {
        display as error "[skip] `splitvar' is entirely missing"
        exit
    }

    preserve
    keep if !missing(`splitvar')

    * ---- separate-sample event studies ----
    quietly reghdfe casi_total trt_m* trt_p* if `splitvar'==1, absorb(subject visit_year visit_age) vce(robust)
    estimates store p1
    quietly reghdfe casi_total trt_m* trt_p* if `splitvar'==0, absorb(subject visit_year visit_age) vce(robust)
    estimates store p0

    * ---- plot ----
    coefplot ///
        (p1, label("`label1'") msymbol(O) mcolor(black) ciopts(recast(rspike) lcolor(gs10) lpattern(dash))) ///
        (p0, label("`label0'") msymbol(Th) mcolor(gs6)  ciopts(recast(rspike) lcolor(gs10) lpattern(dash))), omitted ///
        keep(trt_m6 trt_m5 trt_m4 trt_m3 trt_m2 trt_m1 trt_p0 trt_p1 trt_p2 trt_p3 trt_p4 trt_p5 trt_p6) ///
        order(trt_m6 trt_m5 trt_m4 trt_m3 trt_m2 trt_m1 trt_p0 trt_p1 trt_p2 trt_p3 trt_p4 trt_p5 trt_p6) ///
        vertical xline(7, lp(dash) lcolor(red)) yline(0, lp(dash) lcolor(gs10)) ///
        xlabel(1 "12 yrs before" 2 "10 yrs before" 3 "8 yrs before" 4 "6 yrs before" ///
               5 "4 yrs before"  6 "2 yrs before"  7 "ZVL vaccination year" ///
               8 "2 yrs after"   9 "4 yrs after"  10 "6 yrs after" ///
               11 "8 yrs after" 12 "10 yrs after" 13 "12 yrs after", angle(45) labsize(small) nogrid) ///
        xtitle("Time relative to the year of ZVL vaccination", size(small) margin(t=5)) ///
        ytitle("Effect on CASI Total Score") ///
        legend(title("`legendtitle'", size(vsmall) color(black)) pos(6) rows(1) size(vsmall)) ///
        graphregion(color(white)) bgcolor(white)

    graph export "output/figures/HTE_`outname'.png", replace width(1800) height(1400)

    * ---- displayed difference table: from the two SEPARATE models ----
    estimates restore p1
    matrix b1 = e(b)
    matrix V1 = e(V)
    estimates restore p0
    matrix b0 = e(b)
    matrix V0 = e(V)

    local evlist trt_m6 trt_m5 trt_m4 trt_m3 trt_m2 trt_m1 trt_p0 trt_p1 trt_p2 trt_p3 trt_p4 trt_p5 trt_p6

    tempname P
    tempfile out
    postfile `P' str25 EventTime str18 Estimate str18 LB str18 UB str18 Pvalue str18 BHp using `out', replace

    foreach d of local evlist {
        local lbl ""
        if "`d'"=="trt_m6" local lbl "12 years before"
        if "`d'"=="trt_m5" local lbl "10 years before"
        if "`d'"=="trt_m4" local lbl "8 years before"
        if "`d'"=="trt_m3" local lbl "6 years before"
        if "`d'"=="trt_m2" local lbl "4 years before"
        if "`d'"=="trt_m1" local lbl "2 years before"
        if "`d'"=="trt_p0" local lbl "Vaccination year"
        if "`d'"=="trt_p1" local lbl "2 years after"
        if "`d'"=="trt_p2" local lbl "4 years after"
        if "`d'"=="trt_p3" local lbl "6 years after"
        if "`d'"=="trt_p4" local lbl "8 years after"
        if "`d'"=="trt_p5" local lbl "10 years after"
        if "`d'"=="trt_p6" local lbl "12 years after"

        if "`d'"=="trt_m1" {
            post `P' ("`lbl'") ("Ref") ("") ("") ("") ("")
            continue
        }

        local c1 = colnumb(b1, "`d'")
        local c0 = colnumb(b0, "`d'")
        if missing(`c1') | missing(`c0') {
            post `P' ("`lbl'") ("") ("") ("") ("") ("")
            continue
        }

        local b_1  = b1[1,`c1']
        local se_1 = sqrt(V1[`c1',`c1'])
        local b_0  = b0[1,`c0']
        local se_0 = sqrt(V0[`c0',`c0'])

        local diff = `b_1' - `b_0'
        local sed  = sqrt(`se_1'^2 + `se_0'^2)
        local z    = `diff'/`sed'
        local p    = 2*normal(-abs(`z'))
        local lbnd = `diff' - 1.96*`sed'
        local ubnd = `diff' + 1.96*`sed'

        post `P' ("`lbl'") (string(`diff',"%9.3f")) (string(`lbnd',"%9.3f")) ///
                 (string(`ubnd',"%9.3f")) (string(`p',"%9.3f")) ("")
    }
    postclose `P'

    * ---- pretrend joint Wald p-value: SEPARATE pooled interaction model ----
    * (run BEFORE switching the loaded dataset to the difference table --
    * Stata does not support true nested preserve/restore, so this program
    * never calls preserve more than once; the single outer preserve/restore
    * pair brackets the whole program)
    quietly reghdfe casi_total trt_m6 trt_m5 trt_m4 trt_m3 trt_m2 trt_m1 ///
        trt_p0 trt_p1 trt_p2 trt_p3 trt_p4 trt_p5 trt_p6 ///
        c.`splitvar'#c.trt_m6 c.`splitvar'#c.trt_m5 c.`splitvar'#c.trt_m4 ///
        c.`splitvar'#c.trt_m3 c.`splitvar'#c.trt_m2 ///
        c.`splitvar'#c.trt_p0 c.`splitvar'#c.trt_p1 c.`splitvar'#c.trt_p2 ///
        c.`splitvar'#c.trt_p3 c.`splitvar'#c.trt_p4 c.`splitvar'#c.trt_p5 c.`splitvar'#c.trt_p6, ///
        absorb(subject visit_year visit_age) vce(robust)
    quietly test c.`splitvar'#c.trt_m6 c.`splitvar'#c.trt_m5 c.`splitvar'#c.trt_m4 ///
                 c.`splitvar'#c.trt_m3 c.`splitvar'#c.trt_m2
    local pretrend_p = r(p)

    * ---- BH-adjust the displayed p-values (now safe to switch datasets) ----
    use `out', clear
    tempvar pnum rank m
    gen double `pnum' = real(Pvalue)
    gen byte ok = !missing(`pnum')
    egen long `m' = total(ok)
    gsort -ok `pnum'
    gen long `rank' = sum(ok) if ok
    gen double praw = `pnum' * `m' / `rank' if ok
    gsort -ok -praw
    gen double q = praw
    quietly {
        local N = _N
        forvalues i = 2/`N' {
            if ok[`i'-1]==1 & q[`i'-1] < q[`i'] & q[`i']<. {
                replace q = q[`i'-1] in `i'
            }
        }
    }
    replace q = 1 if q > 1 & q < .
    replace BHp = string(q, "%9.3f") if ok
    drop `pnum' ok `m' `rank' praw q

    export delimited using "output/tables/HTE_`outname'.csv", replace

    file open fh using "output/tables/HTE_`outname'.csv", write append text
    file write fh _n "# `legendtitle': `label1' vs `label0'" _n
    file write fh "# Pretrend joint Wald test p-value (pooled interaction model): " ///
        (string(`pretrend_p',"%6.4f")) _n
    file write fh "# Notes: Same event-study as primary analysis, estimated separately by `legendtitle'; reported coefficients are differences between subgroup estimates. Robust SEs from subgroup regressions. BH-adjusted p-values (Benjamini-Hochberg corrected) are reported to correct for multiple testing." _n
    file close fh

    display as result "[done] `outname': pretrend p = " %6.4f `pretrend_p'

    restore
end

* ==============================================================================
* Figure 4 panels A-E
* ==============================================================================
hte_panel firstdeg_memprob "Family history: memory problems" "No family history: memory problems" "Family history" panelA_TableS17_firstdeg_memprob
hte_panel firstdeg_dem     "Family history: any dementia" "No family history: any dementia" "Family history" panelB_TableS18_firstdeg_dem

cap drop frail_pre
bys subject: egen frail_pre = max(cond(event_time < 0, adl_flag, .))
label define frail 0 "No ADL limitation" 1 "Any ADL limitation", replace
label values frail_pre frail
hte_panel frail_pre "Any ADL limitation" "No ADL limitation" "Frailty status (ADL-based)" panelC_TableS19_frailty_ADL

hte_panel fastdec "Faster baseline decline" "Slower baseline decline" "Baseline rate of cognition decline" panelD_TableS20_baseline_decline

hte_panel high_baseline_vulnerability "Above-median dementia risk" "Below-median dementia risk" "LASSO-predicted baseline dementia risk" panelE_TableS21_LASSO_risk

* ==============================================================================
* Supplementary subgroup-split figures/tables
* ==============================================================================
hte_panel apoe_any4 "APOE e4 carrier" "Non-APOE e4 carrier" "APOE e4 carrier status" S34_apoe4
hte_panel white     "White" "Non-White" "Ethnicity" S35_ethnicity
hte_panel college   "College graduated" "Not college graduated" "College graduation" S37_college

* Figure S36: baseline household income, split at the sample median.
* income is coded on two different scales across survey waves: a pre-2016
* numeric 1-6 scale (9 = don't know/refused, treated as missing) and a
* post-2016 letter A-F scale on a different range. Per Table S22's
* documented restriction, this panel is limited to participants with
* pre-2016 (numeric) income coding.
cap drop income_num income_pre2016 income_med high_income
destring income, generate(income_num) force
replace income_num = . if income_num == 9 | income_num > 6
bys subject (visit_date): gen byte _first_visit = (_n == 1)
bys subject: egen income_pre2016 = max(cond(_first_visit == 1, income_num, .))
drop _first_visit
sum income_pre2016, detail
local income_med = r(p50)
* The income scale is heavily ceiling-effected (median = max = 6), so a
* strict ">" split around the median would be degenerate (empty "above"
* group). Use ">=" so the at-median group falls into "above".
gen byte high_income = (income_pre2016 >= `income_med') if !missing(income_pre2016)
label define incomelbl 0 "Below-median income" 1 "Above-median income", replace
label values high_income incomelbl
hte_panel high_income "Above-median income" "Below-median income" "Baseline household income" S36_income

cap drop female
gen byte female = (gender==2) if !missing(gender)
hte_panel female "Female" "Male" "Gender" S31_gender

* prior other vaccination history (flu/pneumococcal/Td/Tdap), year-level comparison
cap drop zost_year vax_flu_year vax_pneu_year vax_td_year vax_dtap_year
gen int zost_year     = year(zost1_date)
gen int vax_flu_year  = year(flu1_date)  if !missing(flu1_date)
gen int vax_pneu_year = year(pneu1_date) if !missing(pneu1_date)
gen int vax_td_year   = year(td1_date)   if !missing(td1_date)
gen int vax_dtap_year = year(dtap1_date) if !missing(dtap1_date)

cap drop flu_pre pneu_pre td_pre dtap_pre any_vax_pre
gen byte flu_pre  = (vax_flu_year  < zost_year) if !missing(vax_flu_year,  zost_year)
gen byte pneu_pre = (vax_pneu_year < zost_year) if !missing(vax_pneu_year, zost_year)
gen byte td_pre   = (vax_td_year   < zost_year) if !missing(vax_td_year,   zost_year)
gen byte dtap_pre = (vax_dtap_year < zost_year) if !missing(vax_dtap_year, zost_year)
bys subject: egen flu_any_pre  = max(flu_pre)
bys subject: egen pneu_any_pre = max(pneu_pre)
bys subject: egen td_any_pre   = max(td_pre)
bys subject: egen dtap_any_pre = max(dtap_pre)
gen byte any_vax_pre = (flu_any_pre==1 | pneu_any_pre==1 | td_any_pre==1 | dtap_any_pre==1)
hte_panel any_vax_pre "Prior vaccination (flu/pneumococcal/Td/Tdap)" "No prior vaccination" "Prior other vaccination history" S32_prior_vax

* prior herpes simplex history (exact-date comparison, matches HTE.do's hs_prior)
cap drop prior_hsv
gen byte prior_hsv = (!missing(hs1_date) & hs1_date < zost1_date)
hte_panel prior_hsv "Prior herpes simplex history" "No prior herpes simplex history" "Prior herpes simplex history" S33_hsv

* shingles after vaccination (year-level, both-doses comparison, matches HTE.do)
cap drop vax_noshingle_pre shingles_followup hz_year1 hz_year2
gen int hz_year1 = year(hz1_date) if !missing(hz1_date)
gen int hz_year2 = year(hz2_date) if !missing(hz2_date)
gen byte vax_noshingle_pre = !missing(zost_year) & (missing(hz_year1) | hz_year1 >= zost_year) & (missing(hz_year2) | hz_year2 >= zost_year)
gen byte shingles_followup = .
replace shingles_followup = 0 if vax_noshingle_pre==1 & missing(hz_year1) & missing(hz_year2)
replace shingles_followup = 1 if vax_noshingle_pre==1 & !missing(hz_year1) & hz_year1 > zost_year
preserve
keep if vax_noshingle_pre==1
hte_panel shingles_followup "Shingles during follow-up" "No shingles during follow-up" "Shingles after vaccination" S30_shingles_followup
restore

* ==============================================================================
* CASI>60-restricted companion panels for Figures S30 (shingles-during-
* follow-up) and S31 (gender), matching cleaned_lilac/cleaned_data/R code/
* hetergeneity-difference-plot_XI.R (commented-out line: `filter(all(casi_
* total > 60))`): restrict to subjects whose CASI Total Score never fell to
* or below 60 at any visit, as a robustness check excluding participants
* with a period of severe cognitive impairment.
* ==============================================================================
cap drop casi_always_above60
bys subject: egen byte casi_always_above60 = min(casi_total > 60)

preserve
keep if casi_always_above60==1
hte_panel female "Female" "Male" "Gender (CASI>60 restricted)" S31_gender_casi60
restore

preserve
keep if casi_always_above60==1 & vax_noshingle_pre==1
hte_panel shingles_followup "Shingles during follow-up" "No shingles during follow-up" "Shingles after vaccination (CASI>60 restricted)" S30_shingles_followup_casi60
restore

display as result "ALL HTE PANELS DONE"
