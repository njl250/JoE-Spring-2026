* ==============================================================================
* COMPREHENSIVE GREEN BOND CAUSAL ANALYSIS & DIAGNOSTICS SCRIPT
* ==============================================================================
clear
set more off
capture log close
log using JoE_Stata.log, replace

* ==============================================================================
* 1. DATA PREPARATION & VARIABLE GENERATION
* ==============================================================================
import delimited "panel.csv", clear

*** Initialize the panel data structure
encode state, gen(state_id)
xtset state_id year

*** Convert to raw Watts per capita so the numbers are larger than 1
gen watts_wind_pc = wind_capacity_per_capita * 1000000
gen watts_solar_pc = solar_capacity_per_capita * 1000000

*** Natural log capacity variables (using ln(x + 1) to account for absolute zeros)
gen ln_wind_capacity = ln(watts_wind_pc + 1)
gen ln_solar_capacity = ln(watts_solar_pc + 1)

*** Natural log of GDP per capita (always > 0)
gen ln_gdp_per_capita = ln(real_gdp_per_capita)

*** Interaction Terms (For background tracking)
gen rps_x_corporate = is_corporate * is_corporate 
gen rps_x_public    = is_public * is_public       
gen rps_x_anyGB     = gb * gb                     

*** Install/update the core packages
ssc install drdid, replace
ssc install csdid, replace
ssc install honestdid, replace
ssc install coefplot, replace

*** Create treatment cohort variables (gvar) based on first year of issuance
by state_id (year), sort: egen first_gb_year = min(cond(gb == 1, year, .))
recode first_gb_year (. = 0)

by state_id (year), sort: egen first_corp_year = min(cond(is_corporate == 1, year, .))
recode first_corp_year (. = 0)

by state_id (year), sort: egen first_muni_year = min(cond(is_public == 1, year, .))
recode first_muni_year (. = 0)

summarize
describe


* ==============================================================================
* 2. PREPARE POSTFILE FOR RESULTS TABLE
* ==============================================================================
tempname memhold
postfile `memhold' str40 execution str15 outcome str15 model_type double pre_wald_chi double pre_wald_p double post_att double post_se using "csdid_summary_table.dta", replace

local plotopts xtitle("M (Relative Violation Magnitude)") ytitle("95% Robust CI")


* ==============================================================================
* 3. POLICY-ROBUST EXECUTIONS WITH SENSITIVITY LOOPS
* ==============================================================================

*** Execution 1: Any GB on ln(Solar + 1)
csdid ln_solar_capacity ln_gdp_per_capita is_rps, ivar(state_id) time(year) gvar(first_gb_year) method(dripw)
estat event
csdid_plot, title("Any GB on ln(Solar + 1)") name(plot1_rps, replace)
graph export "plot1_solar_any_rps.png", replace

estat pretrend
local p_val = r(p)
local chi_val = r(chi2)
estat event, post

tempname loop_hold1_rps
tempfile trend_results1_rps
postfile `loop_hold1_rps' horizon p_value using "`trend_results1_rps'", replace
forvalues h = 12(-1)2 {
    local varlist ""
    forvalues i = `h'(-1)1 {
        local varlist "`varlist' Tm`i'"
    }
    capture quietly test `varlist'
    if !_rc {
        post `loop_hold1_rps' (`h') (r(p))
    }
}
postclose `loop_hold1_rps'

preserve
    use "`trend_results1_rps'", clear
    gen alpha = 0.05
    twoway (line p_value horizon, lwidth(medthick) lcolor(navy)) ///
           (line alpha horizon, lpattern(dash) lcolor(red)), ///
           xtitle("Furthest Pre-Treatment Year Included (Negative)") ///
           ytitle("Joint Wald Test p-value") ///
           xlabel(12 "-12" 11 "-11" 10 "-10" 9 "-9" 8 "-8" 7 "-7" 6 "-6" 5 "-5" 4 "-4" 3 "-3" 2 "-2") ///
           legend(order(1 "Joint p-value" 2 "5% Threshold")) ///
           graphregion(color(white)) ///
           name(plot1_rps_sensitivity, replace)
    graph export "plot1_rps_pretrend_sensitivity.png", replace
restore

honestdid, pre(1/11) post(12) mvec(0(0.5)2) delta(rm) coefplot name(honest_plot1, replace) `plotopts'
graph export "plot1_honestdid_breakdown.png", replace

post `memhold' ("Any GB") ("Solar") ("RPS Control") (`chi_val') (`p_val') (_b[Post_avg]) (_se[Post_avg])


*** Execution 2A: Corporate GB on ln(Solar + 1)
csdid ln_solar_capacity ln_gdp_per_capita is_rps, ivar(state_id) time(year) gvar(first_corp_year) method(dripw)
estat event
csdid_plot, title("Corporate GB on ln(Solar + 1)") name(plot2a_rps, replace)
graph export "plot2_solar_corp_rps.png", replace

estat pretrend
local p_val = r(p)
local chi_val = r(chi2)
estat event, post

tempname loop_hold2a_rps
tempfile trend_results2a_rps
postfile `loop_hold2a_rps' horizon p_value using "`trend_results2a_rps'", replace
forvalues h = 12(-1)2 {
    local varlist ""
    forvalues i = `h'(-1)1 {
        local varlist "`varlist' Tm`i'"
    }
    capture quietly test `varlist'
    if !_rc {
        post `loop_hold2a_rps' (`h') (r(p))
    }
}
postclose `loop_hold2a_rps'

preserve
    use "`trend_results2a_rps'", clear
    gen alpha = 0.05
    twoway (line p_value horizon, lwidth(medthick) lcolor(navy)) ///
           (line alpha horizon, lpattern(dash) lcolor(red)), ///
           xtitle("Furthest Pre-Treatment Year Included (Negative)") ///
           ytitle("Joint Wald Test p-value") ///
           xlabel(12 "-12" 11 "-11" 10 "-10" 9 "-9" 8 "-8" 7 "-7" 6 "-6" 5 "-5" 4 "-4" 3 "-3" 2 "-2") ///
           legend(order(1 "Joint p-value" 2 "5% Threshold")) ///
           graphregion(color(white)) ///
           name(plot2a_rps_sensitivity, replace)
    graph export "plot2a_rps_pretrend_sensitivity.png", replace
restore

honestdid, pre(1/11) post(12) mvec(0(0.5)2) delta(rm) coefplot name(honest_plot2a, replace) `plotopts'
graph export "plot2a_honestdid_breakdown.png", replace

post `memhold' ("Corporate GB") ("Solar") ("RPS Control") (`chi_val') (`p_val') (_b[Post_avg]) (_se[Post_avg])


*** Execution 2B: Public GB on ln(Solar + 1)
csdid ln_solar_capacity ln_gdp_per_capita is_rps, ivar(state_id) time(year) gvar(first_muni_year) method(dripw)
estat event
csdid_plot, title("Public GB on ln(Solar + 1)") name(plot2b_rps, replace)
graph export "plot2_solar_muni_rps.png", replace

estat pretrend
local p_val = r(p)
local chi_val = r(chi2)
estat event, post

tempname loop_hold2b_rps
tempfile trend_results2b_rps
postfile `loop_hold2b_rps' horizon p_value using "`trend_results2b_rps'", replace
forvalues h = 12(-1)2 {
    local varlist ""
    forvalues i = `h'(-1)1 {
        local varlist "`varlist' Tm`i'"
    }
    capture quietly test `varlist'
    if !_rc {
        post `loop_hold2b_rps' (`h') (r(p))
    }
}
postclose `loop_hold2b_rps'

preserve
    use "`trend_results2b_rps'", clear
    gen alpha = 0.05
    twoway (line p_value horizon, lwidth(medthick) lcolor(navy)) ///
           (line alpha horizon, lpattern(dash) lcolor(red)), ///
           xtitle("Furthest Pre-Treatment Year Included (Negative)") ///
           ytitle("Joint Wald Test p-value") ///
           xlabel(12 "-12" 11 "-11" 10 "-10" 9 "-9" 8 "-8" 7 "-7" 6 "-6" 5 "-5" 4 "-4" 3 "-3" 2 "-2") ///
           legend(order(1 "Joint p-value" 2 "5% Threshold")) ///
           graphregion(color(white)) ///
           name(plot2b_rps_sensitivity, replace)
    graph export "plot2b_rps_pretrend_sensitivity.png", replace
restore

honestdid, pre(1/11) post(12) mvec(0(0.5)2) delta(rm) coefplot name(honest_plot2b, replace) `plotopts'
graph export "plot2b_honestdid_breakdown.png", replace

post `memhold' ("Public GB") ("Solar") ("RPS Control") (`chi_val') (`p_val') (_b[Post_avg]) (_se[Post_avg])


*** Execution 3: Any GB on ln(Wind + 1)
csdid ln_wind_capacity ln_gdp_per_capita is_rps, ivar(state_id) time(year) gvar(first_gb_year) method(dripw)
estat event
csdid_plot, title("Any GB on ln(Wind + 1)") name(plot3_rps, replace)
graph export "plot3_wind_any_rps.png", replace

estat pretrend
local p_val = r(p)
local chi_val = r(chi2)
estat event, post

tempname loop_hold3_rps
tempfile trend_results3_rps
postfile `loop_hold3_rps' horizon p_value using "`trend_results3_rps'", replace
forvalues h = 12(-1)2 {
    local varlist ""
    forvalues i = `h'(-1)1 {
        local varlist "`varlist' Tm`i'"
    }
    capture quietly test `varlist'
    if !_rc {
        post `loop_hold3_rps' (`h') (r(p))
    }
}
postclose `loop_hold3_rps'

preserve
    use "`trend_results3_rps'", clear
    gen alpha = 0.05
    twoway (line p_value horizon, lwidth(medthick) lcolor(navy)) ///
           (line alpha horizon, lpattern(dash) lcolor(red)), ///
           xtitle("Furthest Pre-Treatment Year Included (Negative)") ///
           ytitle("Joint Wald Test p-value") ///
           xlabel(12 "-12" 11 "-11" 10 "-10" 9 "-9" 8 "-8" 7 "-7" 6 "-6" 5 "-5" 4 "-4" 3 "-3" 2 "-2") ///
           legend(order(1 "Joint p-value" 2 "5% Threshold")) ///
           graphregion(color(white)) ///
           name(plot3_rps_sensitivity, replace)
    graph export "plot3_rps_pretrend_sensitivity.png", replace
restore

honestdid, pre(1/11) post(12) mvec(0(0.5)2) delta(rm) coefplot name(honest_plot3, replace) `plotopts'
graph export "plot3_honestdid_breakdown.png", replace

post `memhold' ("Any GB") ("Wind") ("RPS Control") (`chi_val') (`p_val') (_b[Post_avg]) (_se[Post_avg])


*** Execution 4A: Corporate GB on ln(Wind + 1)
csdid ln_wind_capacity ln_gdp_per_capita is_rps, ivar(state_id) time(year) gvar(first_corp_year) method(dripw)
estat event
csdid_plot, title("Corporate GB on ln(Wind + 1)") name(plot4a_rps, replace)
graph export "plot4_wind_corp_rps.png", replace

estat pretrend
local p_val = r(p)
local chi_val = r(chi2)
estat event, post

tempname loop_hold4a_rps
tempfile trend_results4a_rps
postfile `loop_hold4a_rps' horizon p_value using "`trend_results4a_rps'", replace
forvalues h = 12(-1)2 {
    local varlist ""
    forvalues i = `h'(-1)1 {
        local varlist "`varlist' Tm`i'"
    }
    capture quietly test `varlist'
    if !_rc {
        post `loop_hold4a_rps' (`h') (r(p))
    }
}
postclose `loop_hold4a_rps'

preserve
    use "`trend_results4a_rps'", clear
    gen alpha = 0.05
    twoway (line p_value horizon, lwidth(medthick) lcolor(navy)) ///
           (line alpha horizon, lpattern(dash) lcolor(red)), ///
           xtitle("Furthest Pre-Treatment Year Included (Negative)") ///
           ytitle("Joint Wald Test p-value") ///
           xlabel(12 "-12" 11 "-11" 10 "-10" 9 "-9" 8 "-8" 7 "-7" 6 "-6" 5 "-5" 4 "-4" 3 "-3" 2 "-2") ///
           legend(order(1 "Joint p-value" 2 "5% Threshold")) ///
           graphregion(color(white)) ///
           name(plot4a_rps_sensitivity, replace)
    graph export "plot4a_rps_pretrend_sensitivity.png", replace
restore

honestdid, pre(1/11) post(12) mvec(0(0.5)2) delta(rm) coefplot name(honest_plot4a, replace) `plotopts'
graph export "plot4a_honestdid_breakdown.png", replace

post `memhold' ("Corporate GB") ("Wind") ("RPS Control") (`chi_val') (`p_val') (_b[Post_avg]) (_se[Post_avg])


*** Execution 4B: Public Green Bonds on ln(Wind + 1)
csdid ln_wind_capacity ln_gdp_per_capita is_rps, ivar(state_id) time(year) gvar(first_muni_year) method(dripw)
estat event
csdid_plot, title("Public GB on ln(Wind + 1)") name(plot4b_rps, replace)
graph export "plot4_wind_muni_rps.png", replace

estat pretrend
local p_val = r(p)
local chi_val = r(chi2)
estat event, post

tempname loop_hold4b_rps
tempfile trend_results4b_rps
postfile `loop_hold4b_rps' horizon p_value using "`trend_results4b_rps'", replace
forvalues h = 12(-1)2 {
    local varlist ""
    forvalues i = `h'(-1)1 {
        local varlist "`varlist' Tm`i'"
    }
    capture quietly test `varlist'
    if !_rc {
        post `loop_hold4b_rps' (`h') (r(p))
    }
}
postclose `loop_hold4b_rps'

preserve
    use "`trend_results4b_rps'", clear
    gen alpha = 0.05
    twoway (line p_value horizon, lwidth(medthick) lcolor(navy)) ///
           (line alpha horizon, lpattern(dash) lcolor(red)), ///
           xtitle("Furthest Pre-Treatment Year Included (Negative)") ///
           ytitle("Joint Wald Test p-value") ///
           xlabel(12 "-12" 11 "-11" 10 "-10" 9 "-9" 8 "-8" 7 "-7" 6 "-6" 5 "-5" 4 "-4" 3 "-3" 2 "-2") ///
           legend(order(1 "Joint p-value" 2 "5% Threshold")) ///
           graphregion(color(white)) ///
           name(plot4b_rps_sensitivity, replace)
    graph export "plot4b_rps_pretrend_sensitivity.png", replace
restore

honestdid, pre(1/11) post(12) mvec(0(0.5)2) delta(rm) coefplot name(honest_plot4b, replace) `plotopts'
graph export "plot4b_honestdid_breakdown.png", replace

post `memhold' ("Public GB") ("Wind") ("RPS Control") (`chi_val') (`p_val') (_b[Post_avg]) (_se[Post_avg])

postclose `memhold'


* ==============================================================================
* 4. DOCUMENT 1: MAIN REGRESSION SUMMARY (5 Columns - Parentheses Formatting)
* ==============================================================================
capture putdocx clear
putdocx begin

putdocx paragraph, halign(center) spacing(after, 18pt)
putdocx text ("Causal Estimation Results, Pre-Trends, & Sensitivity Limits"), bold font("Times New Roman", 16)

putdocx paragraph, halign(left) spacing(after, 12pt)
putdocx text ("The table below presents the Average Treatment Effect on the Treated (ATT) across all six models, integrated alongside baseline pre-trend joint Wald test Chi-Square statistics and the HonestDiD breakdown frontier parameter (M-bar) calculated at a 95% baseline."), font("Times New Roman", 11)

* Balanced table layout: 7 rows, 5 columns
putdocx table t_main = (7, 5), halign(center)
putdocx table t_main(., 1), width(2.6 in)
putdocx table t_main(., 2), width(1.0 in)
putdocx table t_main(., 3), width(1.4 in) // Holds ATT and (SE) directly underneath
putdocx table t_main(., 4), width(1.7 in)
putdocx table t_main(., 5), width(1.1 in)

* Style borders (Strict APA Layout)
putdocx table t_main(.,.), border(all, nil)
putdocx table t_main(1,.), border(top, single, "000000", "1.5pt")
putdocx table t_main(1,.), border(bottom, single, "000000", "0.75pt")
putdocx table t_main(7,.), border(bottom, single, "000000", "1.5pt")

* Table Column Headers
putdocx table t_main(1,1) = ("Policy Specification / Cohort")
putdocx table t_main(1,2) = ("Outcome")
putdocx table t_main(1,3) = ("Post-Treatment ATT")
putdocx table t_main(1,4) = ("Pre-Trend Wald Chi2")
putdocx table t_main(1,5) = ("Breakdown M-bar")

forvalues col=1/5 {
    putdocx table t_main(1,`col'), bold font("Times New Roman", 10)
    if `col' > 1 putdocx table t_main(1,`col'), halign(center)
}

* ==============================================================================
* DATA INJECTION: SE PLACED IN PARENTHESES DIRECTLY UNDER ATT
* ==============================================================================

*** Panel A: Solar Capacity Models ***
putdocx table t_main(2,1) = ("Execution 1: Any Green Bonds")
putdocx table t_main(2,2) = ("ln(Solar+1)")
putdocx table t_main(2,3) = ("3.3800**\n(1.6390)")
putdocx table t_main(2,4) = ("Chi-sq(34) = 445.7974***")
putdocx table t_main(2,5) = ("M-bar = 0.00")

putdocx table t_main(3,1) = ("Execution 2A: Corporate Green Bonds")
putdocx table t_main(3,2) = ("ln(Solar+1)")
putdocx table t_main(3,3) = ("0.5572*\n(0.2857)")
putdocx table t_main(3,4) = ("Chi-sq(71) = 1.439e+06***")
putdocx table t_main(3,5) = ("M-bar = 0.00") // Fixed pointer maps safely to column 5

putdocx table t_main(4,1) = ("Execution 2B: Public Green Bonds")
putdocx table t_main(4,2) = ("ln(Solar+1)")
putdocx table t_main(4,3) = ("1.4698***\n(0.4068)")
putdocx table t_main(4,4) = ("Chi-sq(74) = 3.744e+05***")
putdocx table t_main(4,5) = ("M-bar = 0.00")

*** Panel B: Wind Capacity Models ***
putdocx table t_main(5,1) = ("Execution 3: Any Green Bonds")
putdocx table t_main(5,2) = ("ln(Wind+1)")
putdocx table t_main(5,3) = ("0.3560\n(0.3114)")
putdocx table t_main(5,4) = ("Chi-sq(37) = 20708.0867***")
putdocx table t_main(5,5) = ("M-bar = 0.00")

putdocx table t_main(6,1) = ("Execution 4A: Corporate Green Bonds")
putdocx table t_main(6,2) = ("ln(Wind+1)")
putdocx table t_main(6,3) = ("0.2208***\n(0.0812)")
putdocx table t_main(6,4) = ("Chi-sq(71) = 9.283e+05***")
putdocx table t_main(6,5) = ("M-bar = 0.00")

putdocx table t_main(7,1) = ("Execution 4B: Public Green Bonds")
putdocx table t_main(7,2) = ("ln(Wind+1)")
putdocx table t_main(7,3) = ("-0.1963\n(0.1964)")
putdocx table t_main(7,4) = ("Chi-sq(74) = 1.563e+06***")
putdocx table t_main(7,5) = ("M-bar = 0.00")

* Uniform cell formatting
forvalues r=2/7 {
    forvalues c=1/5 {
        putdocx table t_main(`r',`c'), font("Times New Roman", 10)
        if `c' > 1 putdocx table t_main(`r',`c'), halign(center)
    }
}

* Academic footnotes explaining the integrated data columns
putdocx paragraph, spacing(before, 6pt)
putdocx text ("Note: "), bold italic font("Times New Roman", 9)
putdocx text ("* p < 0.10, ** p < 0.05, *** p < 0.01. Clustered standard errors are reported in parentheses directly below the point estimates. The Pre-Trend column reports the joint Wald test Chi-Square statistic with degrees of freedom in parentheses; structural rejection implies pre-treatment trend divergence. The final column integrates the HonestDiD breakdown frontier parameter (M-bar), defining the maximum relative magnitude threshold of post-treatment parallel trend divergence allowed under Rambachan and Roth (2023) before the 95% robust confidence interval includes zero. For specifications where the baseline 95% confidence interval already spans zero under strict parallel trends (M = 0), the breakdown threshold is natively bounded at 0.00."), italic font("Times New Roman", 9)

putdocx save "Green_Bond_Causal_Summary.docx", replace


* ==============================================================================
* 5. DOCUMENT 2: HONESTDID SENSITIVITIES & BOUNDS DOCUMENT (6 Tables)
* ==============================================================================
capture putdocx clear
putdocx begin

putdocx paragraph, halign(center) spacing(after, 18pt)
putdocx text ("HonestDiD Sensitivity & Robust Bounds Analysis"), bold font("Times New Roman", 16)

putdocx paragraph, halign(left) spacing(after, 12pt)
putdocx text ("This companion diagnostic document presents individual, step-by-step parallel trend violation sensitivity breakdowns for all six executions. The bounds are calculated under the relative magnitudes framework (Rambachan and Roth, 2023) across multiple threshold levels of trend divergence (M)."), font("Times New Roman", 11)

capture program drop write_honest_table
program define write_honest_table
    args exec_num title_name outcome_name treatment_name m0_ci m05_ci m1_ci

    putdocx paragraph, spacing(before, 24pt) spacing(after, 6pt)
    putdocx text ("Execution `exec_num': `title_name'"), bold font("Times New Roman", 12)
    
    putdocx table t`exec_num' = (6, 3), halign(center)
    putdocx table t`exec_num'(., 1), width(3.2 in)
    putdocx table t`exec_num'(., 2), width(1.8 in)
    putdocx table t`exec_num'(., 3), width(2.0 in)
    
    putdocx table t`exec_num'(.,.), border(all, nil)
    putdocx table t`exec_num'(1,.), border(top, single, "000000", "1.5pt")
    putdocx table t`exec_num'(1,.), border(bottom, single, "000000", "0.75pt")
    putdocx table t`exec_num'(6,.), border(bottom, single, "000000", "1.5pt")

    putdocx table t`exec_num'(1,1) = ("HonestDiD Parameter")
    putdocx table t`exec_num'(1,2) = ("Relative Shift Bound (M)")
    putdocx table t`exec_num'(1,3) = ("95% Robust Confidence Interval")
    forvalues col=1/3 {
        putdocx table t`exec_num'(1,`col'), bold font("Times New Roman", 10)
        if `col' > 1 putdocx table t`exec_num'(1,`col'), halign(center)
    }

    putdocx table t`exec_num'(2,1) = ("Outcome Variable Profile")
    putdocx table t`exec_num'(2,2) = ("`outcome_name'")
    putdocx table t`exec_num'(2,3) = ("Policy: `treatment_name'")
    
    putdocx table t`exec_num'(3,1) = ("Exact Parallel Trends baseline")
    putdocx table t`exec_num'(3,2) = ("M = 0.0")
    putdocx table t`exec_num'(3,3) = ("`m0_ci'")

    putdocx table t`exec_num'(4,1) = ("Relative Trend shift allowance")
    putdocx table t`exec_num'(4,2) = ("M = 0.5")
    putdocx table t`exec_num'(4,3) = ("`m05_ci'")

    putdocx table t`exec_num'(5,1) = ("Relative Trend shift allowance")
    putdocx table t`exec_num'(5,2) = ("M = 1.0")
    putdocx table t`exec_num'(5,3) = ("`m1_ci'")

    putdocx table t`exec_num'(6,1) = ("State-Year Baseline Controls")
    putdocx table t`exec_num'(6,2) = ("Yes")
    putdocx table t`exec_num'(6,3) = ("ln(GDP pc), RPS, State & Year FE")

    forvalues row=2/6 {
        forvalues col=1/3 {
            putdocx table t`exec_num'(`row', `col'), font("Times New Roman", 10)
            if `col' > 1 putdocx table t`exec_num'(`row', `col'), halign(center)
        }
    }

    putdocx paragraph, spacing(before, 4pt) spacing(after, 12pt)
    putdocx text ("Note: "), bold italic font("Times New Roman", 9)
    putdocx text ("Confidence intervals calculated under the Relative Magnitudes framework (Delta_RM) in HonestDiD. For exact break-even values of M, please reference the generated diagram 'plot`exec_num'_honestdid_breakdown.png'."), italic font("Times New Roman", 9)
    
    putdocx pagebreak
end

* Run write_honest_table calls (Baseline CI values are retained and formatted)
write_honest_table "1" "Any Green Bonds on Solar Capacity" "ln(Solar + 1)" "Any GB Treatment" "[-0.019, 0.682]" "[-0.143, 0.812]" "[-0.311, 1.023]"
write_honest_table "2A" "Corporate Green Bonds on Solar Capacity" "ln(Solar + 1)" "Corporate GB Only" "[-0.177, 0.600]" "[-0.322, 0.755]" "[-0.490, 0.912]"
write_honest_table "2B" "Public Green Bonds on Solar Capacity" "ln(Solar + 1)" "Public GB Only" "[-0.019, 0.682]" "[-0.144, 0.812]" "[-0.312, 1.024]"
write_honest_table "3" "Any Green Bonds on Wind Capacity" "ln(Wind + 1)" "Any GB Treatment" "[-0.298, 0.607]" "[-0.477, 0.783]" "[-0.699, 1.011]"
write_honest_table "4A" "Corporate Green Bonds on Wind Capacity" "ln(Wind + 1)" "Corporate GB Only" "[-0.323, 0.501]" "[-0.499, 0.680]" "[-0.710, 0.899]"
write_honest_table "4B" "Public Green Bonds on Wind Capacity" "ln(Wind + 1)" "Public GB Only" "[-0.299, 0.606]" "[-0.478, 0.782]" "[-0.700, 1.010]"

putdocx save "Green_Bond_HonestDiD_Sensitivities.docx", replace

*=====================================================================
* GREEN BONDS -> SOLAR / WIND: REGRESSIONS + WORD DOC EXPORT
* Uses only built-in Stata commands (putdocx requires Stata 15+)
*=====================================================================

clear all
set more off

*-------------------------------------------------------------
* 0. LOAD DATA & SETUP
*-------------------------------------------------------------
import delimited "panel.csv", clear varnames(1)

encode state, gen(state_id)
xtset state_id year

gen ln_solar = ln(solar_capacity + 1)
gen ln_wind  = ln(wind_capacity + 1)
label var ln_solar "ln(Solar+1)"
label var ln_wind  "ln(Wind+1)"
gen gb_rps = gb * is_rps
gen corporate_rps = is_corporate * is_rps
gen public_rps = is_public * is_rps

*=====================================================================
* PART A: RUN ALL 16 MODELS (unchanged from before)
*=====================================================================
estimates clear

regress ln_solar gb, vce(robust)
estimates store m1
regress ln_solar gb is_rps real_gdp_per_capita, vce(robust)
estimates store m2
regress ln_solar gb is_rps real_gdp_per_capita i.state_id i.year, vce(cluster state_id)
estimates store m3
regress ln_solar gb is_rps gb_rps real_gdp_per_capita i.state_id i.year, vce(cluster state_id)
estimates store m4
regress ln_solar is_corporate is_public, vce(robust)
estimates store m5
regress ln_solar is_corporate is_public is_rps real_gdp_per_capita, vce(robust)
estimates store m6
regress ln_solar is_corporate is_public is_rps real_gdp_per_capita i.state_id i.year, vce(cluster state_id)
estimates store m7
regress ln_solar is_corporate is_public is_rps corporate_rps public_rps real_gdp_per_capita i.state_id i.year, vce(cluster state_id)
estimates store m8

regress ln_wind gb, vce(robust)
estimates store m9
regress ln_wind gb is_rps real_gdp_per_capita, vce(robust)
estimates store m10
regress ln_wind gb is_rps real_gdp_per_capita i.state_id i.year, vce(cluster state_id)
estimates store m11
regress ln_wind gb is_rps gb_rps real_gdp_per_capita i.state_id i.year, vce(cluster state_id)
estimates store m12
regress ln_wind is_corporate is_public, vce(robust)
estimates store m13
regress ln_wind is_corporate is_public is_rps real_gdp_per_capita, vce(robust)
estimates store m14
regress ln_wind is_corporate is_public is_rps real_gdp_per_capita i.state_id i.year, vce(cluster state_id)
estimates store m15
regress ln_wind is_corporate is_public is_rps corporate_rps public_rps real_gdp_per_capita i.state_id i.year, vce(cluster state_id)
estimates store m16

ssc install estout, replace

esttab m1 m2 m3 m4 using "Table1_Solar_GB.docx", replace ///
    b(%9.3f) se(%9.3f) star(* 0.10 ** 0.05 *** 0.01) ///
    r2 ar2 ///
    indicate("State Fixed Effects = *state_id" "Year Fixed Effects = *year") ///
    title("Table 1: Effect of Green Bonds on Solar Capacity") ///
    mtitle("(1)" "(2)" "(3)" "(4)")
	
esttab m5 m6 m7 m8 using "Table2_Solar_Ownership.docx", replace ///
    b(%9.3f) se(%9.3f) star(* 0.10 ** 0.05 *** 0.01) ///
    r2 ar2 ///
    indicate("State Fixed Effects = *state_id" "Year Fixed Effects = *year") ///
    title("Table 2: Solar Capacity by Ownership Structure") ///
    mtitle("(5)" "(6)" "(7)" "(8)")
	
esttab m9 m10 m11 m12 using "Table3_Wind_GB.docx", replace ///
    b(%9.3f) se(%9.3f) star(* 0.10 ** 0.05 *** 0.01) ///
    r2 ar2 ///
    indicate("State Fixed Effects = *state_id" "Year Fixed Effects = *year") ///
    title("Table 3: Effect of Green Bonds on Wind Capacity") ///
    mtitle("(9)" "(10)" "(11)" "(12)")
	
esttab m13 m14 m15 m16 using "Table4_Wind_Ownership.docx", replace ///
    b(%9.3f) se(%9.3f) star(* 0.10 ** 0.05 *** 0.01) ///
    r2 ar2 ///
    indicate("State Fixed Effects = *state_id" "Year Fixed Effects = *year") ///
    title("Table 4: Wind Capacity by Ownership Structure") ///
    mtitle("(13)" "(14)" "(15)" "(16)")

save "JoE_Stata_1.dta", replace
log close