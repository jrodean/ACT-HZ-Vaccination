* ==============================================================================
* HTE_all.do
*
* Every heterogeneity-of-treatment-effect (HTE) result: Figure 4 panels A-E
* and the 8 supplementary subgroup-split figures/tables (S30-S37, with
* Tables S17-S21 for panels A-E). Ported section-by-section from the
* original Xi/v0926/HTE_1026.do, using the original's own frozen
* data/raw/event_all_120125.dta (+ event_all_120125_with_risk.csv,
* decline01.csv, EHR.dta, ACT_cleaned.csv) directly -- the original rebuilds
* its working sample independently per section (several standalone `use`
* + section-specific merges), so this file does the same instead of sharing
* one filtered panel across all panels.
*
* Run from the "replication package" directory: stata -b do stata/HTE_all.do
* ==============================================================================

cd "`c(pwd)'"
set scheme plotplain

cap mkdir output
cap mkdir output/figures
cap mkdir output/tables

* ------------------------------------------------------------------------
* hz_beforeafter (HTE_1026.do lines 46-62): shingles before/after ZVL,
* from the EHR event file -- used by the S30 shingles-follow-up panel.
* ------------------------------------------------------------------------
use "data/raw/EHR.dta", clear
keep if event == "ep_hz"
rename days_to_event days_to_hz
tempfile hz_long
save `hz_long'

use "data/raw/event_all_120125.dta", clear
keep if zost1_date != .
keep zost1_date subject days_to_zost1
bysort subject: keep if _n == 1
merge 1:m subject using `hz_long'
drop if _merge == 2
drop _merge
gen hz_after_zost1  = days_to_hz > days_to_zost1 & days_to_hz != .
gen hz_before_zost1 = days_to_hz < days_to_zost1 & days_to_hz != .
collapse (max) hz_after_zost1 hz_before_zost1, by(subject)
tempfile hz_beforeafter
save `hz_beforeafter'

* ------------------------------------------------------------------------
* Base panel (HTE_1026.do lines 65-130): merge ML risk score + baseline
* decline rate, keep vaccinated, build event_time/trt_* [-14,12] and fastdec.
* ------------------------------------------------------------------------
import delimited "data/raw/event_all_120125_with_risk.csv", clear
tempfile vulnfile
save `vulnfile'

import delimited "data/raw/decline01.csv", clear
destring decline01, replace force
tempfile declfile
save `declfile'

use "data/raw/event_all_120125.dta", clear
merge m:1 subject using `vulnfile', keep(match master) nogen
merge m:1 subject using `declfile', keep(match master) nogen

keep if !missing(zost1_date)

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

cap drop fastdec
gen byte fastdec = (decline01 >= 0) if !missing(decline01)

save `"`c(tmpdir)'base_panel.dta"', replace

* ------------------------------------------------------------------------
* hte_panel: one HTE panel (plot + difference table) for a 0/1 splitvar.
* Matches HTE_1026.do's baselineplot/coefplot style exactly.
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

    keep if !missing(`splitvar')

    quietly reghdfe casi_total trt_m* trt_p* if `splitvar'==1, absorb(subject visit_year visit_age) vce(robust)
    estimates store p1
    quietly reghdfe casi_total trt_m* trt_p* if `splitvar'==0, absorb(subject visit_year visit_age) vce(robust)
    estimates store p0

    estimates restore p1
    matrix b1 = e(b)
    matrix V1 = e(V)
    estimates restore p0
    matrix b0 = e(b)
    matrix V0 = e(V)

    local evlist trt_m6 trt_m5 trt_m4 trt_m3 trt_m2 trt_m1 trt_p0 trt_p1 trt_p2 trt_p3 trt_p4 trt_p5 trt_p6

    * A term absorbed/collinear in either subgroup regression (e.g. too few
    * observations that far from the vaccination date -- reghdfe keeps such
    * a column but pins it at coefficient 0 with 0 variance, rather than
    * dropping it) is excluded from the panel entirely here, rather than
    * shown as a spurious zero-effect point with no CI. This is also why
    * the published figure for some panels (e.g. S36, baseline income)
    * stops short of +/-12 years.
    local keeplist ""
    local xlab ""
    local pos = 0
    local vaxpos = .
    foreach d of local evlist {
        local lbl ""
        if "`d'"=="trt_m6" local lbl "12 yrs before"
        if "`d'"=="trt_m5" local lbl "10 yrs before"
        if "`d'"=="trt_m4" local lbl "8 yrs before"
        if "`d'"=="trt_m3" local lbl "6 yrs before"
        if "`d'"=="trt_m2" local lbl "4 yrs before"
        if "`d'"=="trt_m1" local lbl "2 yrs before"
        if "`d'"=="trt_p0" local lbl "ZVL vaccination year"
        if "`d'"=="trt_p1" local lbl "2 yrs after"
        if "`d'"=="trt_p2" local lbl "4 yrs after"
        if "`d'"=="trt_p3" local lbl "6 yrs after"
        if "`d'"=="trt_p4" local lbl "8 yrs after"
        if "`d'"=="trt_p5" local lbl "10 yrs after"
        if "`d'"=="trt_p6" local lbl "12 yrs after"

        local keepit = 0
        if "`d'"=="trt_m1" local keepit = 1
        if `keepit' == 0 {
            local c1 = colnumb(b1, "`d'")
            local c0 = colnumb(b0, "`d'")
            if !missing(`c1') & !missing(`c0') {
                if V1[`c1',`c1'] > 1e-12 & V0[`c0',`c0'] > 1e-12 local keepit = 1
            }
        }
        if `keepit' {
            local pos = `pos' + 1
            local keeplist `"`keeplist' `d'"'
            local xlab `" `xlab' `pos' "`lbl'" "'
            if "`d'"=="trt_p0" local vaxpos = `pos'
        }
    }

    coefplot ///
        (p1, label("`label1'") msymbol(O) mcolor(black) ciopts(recast(rspike) lcolor(gs10) lpattern(dash))) ///
        (p0, label("`label0'") msymbol(Th) mcolor(gs6)  ciopts(recast(rspike) lcolor(gs10) lpattern(dash))), ///
        keep(`keeplist') order(`keeplist') ///
        vertical xline(`vaxpos', lp(dash) lcolor(red)) yline(0, lp(dash) lcolor(gs10)) ///
        xlabel(`xlab', angle(45) labsize(small) nogrid) ///
        ylabel(, nogrid angle(0)) ///
        xtitle("Time relative to the year of ZVL vaccination", size(small) margin(t=5)) ///
        ytitle("Effect on CASI Total Score") ///
        legend(title("`legendtitle'", size(vsmall) color(black)) pos(6) rows(1) size(vsmall)) ///
        graphregion(color(white) lcolor(none)) plotregion(lcolor(none) lwidth(none)) bgcolor(white)

    graph export "output/figures/HTE_`outname'.png", replace width(1800) height(1400)

    tempname P
    tempfile out
    postfile `P' str25 EventTime str18 Estimate str18 LB str18 UB str18 Pvalue str18 BHp double Praw using `out', replace

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
            post `P' ("`lbl'") ("Ref") ("") ("") ("") ("") (.)
            continue
        }

        local c1 = colnumb(b1, "`d'")
        local c0 = colnumb(b0, "`d'")
        if missing(`c1') | missing(`c0') {
            post `P' ("`lbl'") ("") ("") ("") ("") ("") (.)
            continue
        }
        if V1[`c1',`c1'] <= 1e-12 | V0[`c0',`c0'] <= 1e-12 {
            post `P' ("`lbl'") ("") ("") ("") ("") ("") (.)
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
                 (string(`ubnd',"%9.3f")) (string(`p',"%9.3f")) ("") (`p')
    }
    postclose `P'

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

    use `out', clear
    tempvar pnum rank m
    gen double `pnum' = Praw
    gen byte ok = !missing(`pnum')
    egen long `m' = total(ok)
    gsort -ok `pnum'
    gen long `rank' = sum(ok) if ok
    gen double praw = `pnum' * `m' / `rank' if ok
    gsort -ok -`rank'
    gen double q = praw
    quietly {
        local N = _N
        forvalues i = 2/`N' {
            if ok[`i'-1]==1 & ok[`i']==1 & q[`i'-1] < q[`i'] {
                replace q = q[`i'-1] in `i'
            }
        }
    }
    replace q = 1 if q > 1 & q < .
    replace BHp = string(q, "%9.3f") if ok
    drop `pnum' ok `m' `rank' praw q Praw

    export delimited using "output/tables/HTE_`outname'.csv", replace

    file open fh using "output/tables/HTE_`outname'.csv", write append text
    file write fh _n "# `legendtitle': `label1' vs `label0'" _n
    file write fh "# Pretrend joint Wald test p-value (pooled interaction model): " ///
        (string(`pretrend_p',"%6.4f")) _n
    file write fh "# Notes: Same event-study as primary analysis, estimated separately by `legendtitle'; reported coefficients are differences between subgroup estimates. Robust SEs from subgroup regressions. BH-adjusted p-values (Benjamini-Hochberg corrected) are reported to correct for multiple testing." _n
    file close fh

    display as result "[done] `outname': pretrend p = " %6.4f `pretrend_p'
end

* ==============================================================================
* Figure 4 panels A-E (HTE_1026.do "Main Text" section)
* ==============================================================================
use `"`c(tmpdir)'base_panel.dta"', clear
hte_panel firstdeg_memprob "Memory problems" "No memory problems" "Family history" panelA_TableS17_firstdeg_memprob

use `"`c(tmpdir)'base_panel.dta"', clear
hte_panel firstdeg_dem "Any dementia" "No dementia" "Family history" panelB_TableS18_firstdeg_dem

use `"`c(tmpdir)'base_panel.dta"', clear
cap drop frail_pre
bys subject: egen frail_pre = max(cond(event_time < 0, adl_flag, .))
hte_panel frail_pre "Any ADL limitation" "No ADL limitation" "Frailty status (ADL-based)" panelC_TableS19_frailty_ADL

use `"`c(tmpdir)'base_panel.dta"', clear
hte_panel fastdec "Faster baseline decline" "Slower baseline decline" "Baseline rate of cognition decline" panelD_TableS20_baseline_decline

use `"`c(tmpdir)'base_panel.dta"', clear
hte_panel high_baseline_vulnerability "Above-median dementia risk" "Below-median dementia risk" "LASSO-predicted baseline dementia risk" panelE_TableS21_LASSO_risk

* ==============================================================================
* Supplementary subgroup-split figures/tables (HTE_1026.do "APPENDIX" section)
* ==============================================================================
use `"`c(tmpdir)'base_panel.dta"', clear
label var apoe_any4 "APOE e4 carrier"
hte_panel apoe_any4 "APOE e4 carrier" "Non-APOE e4 carrier" "APOE e4 carrier status" S34_apoe4

* white, college: recomputed inline exactly as HTE_1026.do lines 854-901,
* not from the dta's own (differently-sourced) white/college columns.
use `"`c(tmpdir)'base_panel.dta"', clear
cap drop white
gen byte white = (raceeth_c==7) if !missing(raceeth_c)
hte_panel white "White" "Non-White" "Ethnicity" S35_ethnicity

use `"`c(tmpdir)'base_panel.dta"', clear
cap drop college
gen byte college = (educat4 == 3) if !missing(educat4)
hte_panel college "College graduated" "Not college graduated" "College graduation" S37_college

* female: HTE_1026.do line 725 uses sex_c (ACT coding: 2=female), not gender.
use `"`c(tmpdir)'base_panel.dta"', clear
cap drop female
gen byte female = (sex_c==2) if !missing(sex_c)
hte_panel female "Female" "Male" "Gender" S31_gender

* any_vax_pre: HTE_1026.do lines 600-629, using the dta's own vax_flu/
* vax_pneu/vax_td/vax_dtap (already present -- no merge needed).
use `"`c(tmpdir)'base_panel.dta"', clear
cap drop zost_year vax_flu_year vax_pneu_year vax_td_year vax_dtap_year
gen int zost_year     = year(zost1_date)
gen int vax_flu_year  = year(vax_flu)  if !missing(vax_flu)
gen int vax_pneu_year = year(vax_pneu) if !missing(vax_pneu)
gen int vax_td_year   = year(vax_td)   if !missing(vax_td)
gen int vax_dtap_year = year(vax_dtap) if !missing(vax_dtap)
cap drop flu_pre pneu_pre td_pre dtap_pre
gen byte flu_pre  = (vax_flu_year  < zost_year) if !missing(vax_flu_year,  zost_year)
gen byte pneu_pre = (vax_pneu_year < zost_year) if !missing(vax_pneu_year, zost_year)
gen byte td_pre   = (vax_td_year   < zost_year) if !missing(vax_td_year,   zost_year)
gen byte dtap_pre = (vax_dtap_year < zost_year) if !missing(vax_dtap_year, zost_year)
cap drop flu_any_pre pneu_pre_any td_pre_any dtap_pre_any any_vax_pre
bys subject: egen flu_any_pre   = max(flu_pre)
bys subject: egen pneu_pre_any  = max(pneu_pre)
bys subject: egen td_pre_any    = max(td_pre)
bys subject: egen dtap_pre_any  = max(dtap_pre)
gen byte any_vax_pre = (flu_any_pre==1 | pneu_pre_any==1 | td_pre_any==1 | dtap_pre_any==1)
hte_panel any_vax_pre "Prior vaccination (flu/pneumococcal/Td/Tdap)" "No prior vaccination" "Prior other vaccination history" S32_prior_vax

* ------------------------------------------------------------------------
* prior_hsv (S33): HTE_1026.do lines 941-1101. HS date comes from
* ACT_cleaned.csv's date_ep_hs1 (the dta has no hs1_date), merged by subject.
* ------------------------------------------------------------------------
import delimited "data/raw/ACT_cleaned.csv", clear varnames(1)
keep subject date_ep_hs1 date_av_acycl1 date_av_acycl2 date_av_acycl3 ///
     date_av_valac1 date_av_valac2 date_av_valac3
bysort subject: keep if _n==1
tempfile act_dates
save `act_dates', replace

use "data/raw/event_all_120125.dta", clear
merge m:1 subject using `act_dates'
drop if _merge==2
drop _merge

format zost1_date %td
capture confirm numeric variable visit_date
if _rc {
    gen double visit_date_d = date(visit_date, "YMD") if !missing(visit_date) & visit_date!=""
    format visit_date_d %td
    drop visit_date
    rename visit_date_d visit_date
}
format visit_date %td

keep if !missing(zost1_date)

cap drop hs_date
gen double hs_date = .
capture confirm string variable date_ep_hs1
if _rc==0 {
    replace hs_date = date(date_ep_hs1, "YMD") if !missing(date_ep_hs1) & date_ep_hs1!=""
    replace hs_date = date(date_ep_hs1, "DMY") if missing(hs_date) & !missing(date_ep_hs1) & date_ep_hs1!=""
}
else {
    replace hs_date = date_ep_hs1
    replace hs_date = . if hs_date==0
}
format hs_date %td

cap drop prior_hsv
gen byte prior_hsv = (!missing(hs_date) & hs_date < zost1_date)

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

hte_panel prior_hsv "Yes" "No" "Prior herpes simplex history" S33_hsv

* ------------------------------------------------------------------------
* shingles_followup (S30): HTE_1026.do lines 800-847, via the
* EHR-based hz_beforeafter merge (hz_after_zost1/hz_before_zost1), not the
* dta's own hz_year1/hz_year2 columns.
* ------------------------------------------------------------------------
use `"`c(tmpdir)'base_panel.dta"', clear
merge m:1 subject using `hz_beforeafter'
drop _merge
cap drop shingles_followup vax_noshingle_pre
gen byte shingles_followup = hz_after_zost1==1
gen byte vax_noshingle_pre = hz_before_zost1==0
keep if vax_noshingle_pre==1
hte_panel shingles_followup "Shingles during follow-up" "No shingles during follow-up" "Shingles after vaccination" S30_shingles_followup

* ------------------------------------------------------------------------
* income (S36): HTE_1026.do lines 1104-1237. act_cleaned socio-demo columns
* merged by subject+visit_year (keep _merge==3 -- changes the row sample),
* pre-2016 numeric 1-6 income scale only, median split.
* ------------------------------------------------------------------------
import delimited "data/raw/ACT_cleaned.csv", clear varnames(1)
destring visit_year, replace ignore(" ,")
tostring subject, replace
keep sex_c race_c hispanic_c raceeth_c education degree educat4 college ///
     marital marital_change ever_worked current_work current_work_hours ///
     retired_bl retired_age_bl retired_reason_bl ///
     retired_vis retired_age_vis retired_reason_vis ///
     income basicneeds_met basicneeds_met_food ///
     basicneeds_met_housing basicneeds_met_clothing ///
     basicneeds_met_medical afford_luxury public_assist ///
     income_sufficient subject visit_year
bysort subject visit_year: keep if _n==1
tempfile act_tmp
save `act_tmp'

use "data/raw/event_all_120125.dta", clear
destring visit_year, replace ignore(" ,")
tostring subject, replace
merge m:1 subject visit_year using `act_tmp'
keep if _merge == 3
drop _merge
destring subject, replace

keep if !missing(zost1_date)
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

drop if inlist(income,"A","B","C","D","E","F","9")
cap drop income_ord
gen income_ord = .
replace income_ord = 1 if inlist(income,"1","2")
replace income_ord = 2 if inlist(income,"3","4")
replace income_ord = 3 if income=="5"
replace income_ord = 4 if income=="6"
cap drop inc_ord_base
bys subject: egen inc_ord_base = mean(cond(event_time < 0 & visit_year < 2016, income_ord, .))
summ inc_ord_base, detail
scalar inc_med = r(p50)
cap drop income_above_base
gen byte income_above_base = inc_ord_base >= inc_med if !missing(inc_ord_base)

hte_panel income_above_base "At/above-median baseline income" "Below-median baseline income" "Baseline household income" S36_income

* ==============================================================================
* CASI>60-restricted companion panels for S30 (shingles-during-follow-up)
* and S31 (gender): subjects whose CASI Total Score never fell below 60.
* ==============================================================================
use `"`c(tmpdir)'base_panel.dta"', clear
cap drop casi_always_above60
bys subject: egen byte casi_always_above60 = min(casi_total >= 60)
cap drop female
gen byte female = (sex_c==2) if !missing(sex_c)
keep if casi_always_above60==1
hte_panel female "Female" "Male" "Gender (CASI>60 restricted)" S31_gender_casi60

use `"`c(tmpdir)'base_panel.dta"', clear
merge m:1 subject using `hz_beforeafter'
drop _merge
bys subject: egen byte casi_always_above60 = min(casi_total >= 60)
cap drop shingles_followup vax_noshingle_pre
gen byte shingles_followup = hz_after_zost1==1
gen byte vax_noshingle_pre = hz_before_zost1==0
keep if casi_always_above60==1 & vax_noshingle_pre==1
hte_panel shingles_followup "Shingles during follow-up" "No shingles during follow-up" "Shingles after vaccination (CASI>60 restricted)" S30_shingles_followup_casi60

* ==============================================================================
* Sample sizes for the four subgroup splits built only on the Stata side
* (frailty, prior other vaccination, prior herpes simplex, income), for
* Table S22. Read back by R/27_finalize_table_s22.R.
* ==============================================================================
tempname ssz
file open `ssz' using "output/tables/TableS22_stata_sample_sizes.csv", write replace
file write `ssz' "Panel,Level,N_individuals,N_person_visits" _n

use `"`c(tmpdir)'base_panel.dta"', clear
cap drop frail_pre
bys subject: egen frail_pre = max(cond(event_time < 0, adl_flag, .))
foreach v in 1 0 {
    quietly count if frail_pre == `v'
    local nvis = r(N)
    quietly bys subject: gen byte _f = (_n == 1) if frail_pre == `v'
    quietly count if _f == 1
    local nind = r(N)
    drop _f
    file write `ssz' `""Any ADL limitation pre-vaccination (Table S19)",`v',`nind',`nvis'"' _n
}

use `"`c(tmpdir)'base_panel.dta"', clear
cap drop zost_year vax_flu_year vax_pneu_year vax_td_year vax_dtap_year
gen int zost_year     = year(zost1_date)
gen int vax_flu_year  = year(vax_flu)  if !missing(vax_flu)
gen int vax_pneu_year = year(vax_pneu) if !missing(vax_pneu)
gen int vax_td_year   = year(vax_td)   if !missing(vax_td)
gen int vax_dtap_year = year(vax_dtap) if !missing(vax_dtap)
gen byte flu_pre  = (vax_flu_year  < zost_year) if !missing(vax_flu_year,  zost_year)
gen byte pneu_pre = (vax_pneu_year < zost_year) if !missing(vax_pneu_year, zost_year)
gen byte td_pre   = (vax_td_year   < zost_year) if !missing(vax_td_year,   zost_year)
gen byte dtap_pre = (vax_dtap_year < zost_year) if !missing(vax_dtap_year, zost_year)
bys subject: egen flu_any_pre   = max(flu_pre)
bys subject: egen pneu_pre_any  = max(pneu_pre)
bys subject: egen td_pre_any    = max(td_pre)
bys subject: egen dtap_pre_any  = max(dtap_pre)
gen byte any_vax_pre = (flu_any_pre==1 | pneu_pre_any==1 | td_pre_any==1 | dtap_pre_any==1)
foreach v in 1 0 {
    quietly count if any_vax_pre == `v'
    local nvis = r(N)
    quietly bys subject: gen byte _f = (_n == 1) if any_vax_pre == `v'
    quietly count if _f == 1
    local nind = r(N)
    drop _f
    file write `ssz' `""Prior other vaccination (flu/pneumococcal/Td/Tdap)",`v',`nind',`nvis'"' _n
}

import delimited "data/raw/ACT_cleaned.csv", clear varnames(1)
keep subject date_ep_hs1 date_av_acycl1 date_av_acycl2 date_av_acycl3 ///
     date_av_valac1 date_av_valac2 date_av_valac3
bysort subject: keep if _n==1
tempfile act_dates2
save `act_dates2', replace
use "data/raw/event_all_120125.dta", clear
merge m:1 subject using `act_dates2'
drop if _merge==2
drop _merge
keep if !missing(zost1_date)
gen double hs_date = .
capture confirm string variable date_ep_hs1
if _rc==0 {
    replace hs_date = date(date_ep_hs1, "YMD") if !missing(date_ep_hs1) & date_ep_hs1!=""
    replace hs_date = date(date_ep_hs1, "DMY") if missing(hs_date) & !missing(date_ep_hs1) & date_ep_hs1!=""
}
else {
    replace hs_date = date_ep_hs1
    replace hs_date = . if hs_date==0
}
gen byte prior_hsv = (!missing(hs_date) & hs_date < zost1_date)
foreach v in 1 0 {
    quietly count if prior_hsv == `v'
    local nvis = r(N)
    quietly bys subject: gen byte _f = (_n == 1) if prior_hsv == `v'
    quietly count if _f == 1
    local nind = r(N)
    drop _f
    file write `ssz' `""Prior herpes simplex episode",`v',`nind',`nvis'"' _n
}

import delimited "data/raw/ACT_cleaned.csv", clear varnames(1)
destring visit_year, replace ignore(" ,")
tostring subject, replace
keep income subject visit_year
bysort subject visit_year: keep if _n==1
tempfile act_tmp2
save `act_tmp2'
use "data/raw/event_all_120125.dta", clear
destring visit_year, replace ignore(" ,")
tostring subject, replace
merge m:1 subject visit_year using `act_tmp2'
keep if _merge == 3
drop _merge
destring subject, replace
keep if !missing(zost1_date)
gen int event_time = floor((visit_date - zost1_date)/(2*365))*2
drop if inlist(income,"A","B","C","D","E","F","9")
gen income_ord = .
replace income_ord = 1 if inlist(income,"1","2")
replace income_ord = 2 if inlist(income,"3","4")
replace income_ord = 3 if income=="5"
replace income_ord = 4 if income=="6"
bys subject: egen inc_ord_base = mean(cond(event_time < 0 & visit_year < 2016, income_ord, .))
summ inc_ord_base, detail
scalar inc_med2 = r(p50)
gen byte income_above_base = inc_ord_base >= inc_med2 if !missing(inc_ord_base)
foreach v in 1 0 {
    quietly count if income_above_base == `v'
    local nvis = r(N)
    quietly bys subject: gen byte _f = (_n == 1) if income_above_base == `v'
    quietly count if _f == 1
    local nind = r(N)
    drop _f
    file write `ssz' `""At/above-median baseline household income",`v',`nind',`nvis'"' _n
}

file close `ssz'

display as result "ALL HTE PANELS DONE"
