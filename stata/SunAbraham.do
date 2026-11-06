* ==============================================================================
* SunAbraham.do
*
* Figure S20. Sun & Abraham (2021) interaction-weighted event-study
* estimator, via Stata's own `eventstudyinteract` -- ported line-for-line
* from the "SUN" section of the original analysis_Xi_1225.do (found in
* Xi/v0326/ and Xi/v1225/; the first of four eventstudyinteract blocks in
* that file, the one on zost1_date/casi_total -- not the later Shingrix/RSV
* blocks further down the same file): cohort(vax_age) (age at vaccination),
* control_cohort(last_vax), last_vax = vaccinated in 2017, sample dropped
* to visits before 2017 so the 2017 cohort's surviving rows are all
* pre-vaccination (its own "last treated" control role), event time capped
* to [-14, +8] (the original code's own cap -- not [-12,+8], which is what
* the Supplement's prose paraphrase of this figure states; the code is
* followed exactly here per instruction, not reconciled with the prose).
*
* Requires: ssc install eventstudyinteract avar
* Input: data/analysis_panel.dta (written by R/00c_lasso_risk_score.R).
* Run from the "replication package" directory: stata -b do stata/SunAbraham.do
* ==============================================================================
clear all
set more off
cap ssc install eventstudyinteract
local plotplain "plotregion(lcolor(none)) graphregion(color(white) lcolor(none) ilcolor(none)) bgcolor(white)"

use "data/analysis_panel.dta", clear
keep if !missing(zost1_date)

bys subject: egen birthyear = max(cond(days_since_baseline == 0, visit_year - b_visit_age, .))
gen vax_year = year(zost1_date)
gen vax_age = vax_year - birthyear
gen last_vax = (vax_year == 2017)
drop if visit_year >= 2017

gen event_time = floor((visit_date - zost1_date)/(2*365)) * 2
replace event_time = 8   if event_time > 8  & event_time != .
replace event_time = -14 if event_time < -14 & event_time != .

forvalues i = 1/6 {
    gen trt_m`i' = (event_time == -`i'*2)
}
forvalues i = 0/4 {
    gen trt_p`i' = (event_time == `i'*2)
}
replace trt_m1 = 0

eventstudyinteract casi_total trt_m1-trt_m6 trt_p0 trt_p1 trt_p2 trt_p3 trt_p4, ///
    absorb(subject visit_year visit_age) cohort(vax_age) control_cohort(last_vax) vce(robust)

matrix b_iw = e(b_iw)
matrix V_iw = e(V_iw)

tempname P
tempfile out
postfile `P' str20 term double estimate double se using `out', replace
local cnames : colnames b_iw
local i = 1
foreach c of local cnames {
    post `P' ("`c'") (b_iw[1, `i']) (sqrt(V_iw[`i', `i']))
    local ++i
}
postclose `P'
use `out', clear
gen ci_low = estimate - 1.96*se
gen ci_high = estimate + 1.96*se
export delimited using "output/tables/TableS20_sunab_estimator.csv", replace
list

coefplot, b(b_iw) v(V_iw) baselevels omitted ///
    mcolor(black) m(circle) msize(medium) lcolor(black) ///
    ciopts(recast(rspike) lcolor(gs10) lpattern(dash)) levels(95) vertical ///
    keep(trt_*) ///
    order(trt_m6 trt_m5 trt_m4 trt_m3 trt_m2 trt_m1 trt_p0 trt_p1 trt_p2 trt_p3 trt_p4) ///
    xlabel(1 "12 yrs before" 2 "10 yrs before" 3 "8 yrs before" 4 "6 yrs before" ///
           5 "4 yrs before" 6 "2 yrs before" 7 "ZVL vaccination year" ///
           8 "2 yrs after" 9 "4 yrs after" 10 "6 yrs after" 11 "8 yrs after", ///
           labsize(small) angle(45) nogrid) ///
    ylabel(, angle(0) labsize(medium) nogrid) ///
    xtitle("Time relative to the year of ZVL vaccination", size(medium) margin(small)) ///
    ytitle("Effect on CASI Total Score", size(medium) margin(small)) ///
    xline(7, lp(dash) lcolor(red) lwidth(medium)) ///
    yline(0, lp(dash) lcolor(gs10) lwidth(medium)) ///
    xscale(range(-0.1 11.5)) ///
    text(0.5 0.1 "Benefit", placement(e) size(small) color(black)) ///
    text(-0.5 0.1 "Harm", placement(e) size(small) color(black)) ///
    addplot( ///
    pcarrowi 0.00 0.3 0.3 0.3, lcolor(black) lwidth(thin) mcolor(black) ///
    || pcarrowi 0.00 0.3 -0.3 0.3, lcolor(black) lwidth(thin) mcolor(black) ///
    ) `plotplain'
graph export "output/figures/FigureS20_sunab_estimator.png", replace width(1800) height(1200)

display as result "SunAbraham.do DONE"
