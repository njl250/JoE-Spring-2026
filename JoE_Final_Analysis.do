* ==============================================================================
* GREEN BOND FINANCING AND STATE RENEWABLE CAPACITY (2010-2024)
* Author: Noah Leibowitz
*
* Data:       panel_improved.csv (on the Desktop)
* Outcomes:   ln_solar_pc, ln_wind_pc  = ln(W per person + 1)
* Treatment:  gb_solar = green bond proceeds funded in-state solar
*             gb_wind  = green bond proceeds funded in-state wind
* Controls:   ln_gdp_pc (ln real GDP per capita), is_rps (RPS in force)
*
* All output is saved to the Desktop:
*   Table_TWFE_Solar.rtf, Table_TWFE_Wind.rtf
*   cs_event_solar.png, cs_event_wind.png
*   cs_honestdid_solar.png, cs_honestdid_wind.png
*   CSDiD_Results.docx, CSDiD_Results.csv
*   JoE_Final_Analysis.log
* ==============================================================================

clear all
set more off
capture log close

* Work from the Desktop and save everything there
global DESKTOP "C:/Users/njl119/OneDrive - Case Western Reserve University/Desktop"
cd "$DESKTOP"

log using "$DESKTOP/JoE_Final_Analysis.log", replace

* Install packages if needed
foreach pkg in estout drdid csdid honestdid coefplot {
    capture which `pkg'
    if _rc ssc install `pkg', replace
}

* ------------------------------------------------------------------------------
* LOAD DATA
* ------------------------------------------------------------------------------
import delimited "$DESKTOP/panel_improved.csv", clear varnames(1) case(lower)
xtset state_id year

* Treatment x RPS interactions
gen solar_rps = gb_solar * is_rps
gen wind_rps  = gb_wind  * is_rps

* First year each state was treated (0 = never treated), for csdid
foreach t in gb_solar gb_wind {
    bysort state_id: egen g_`t' = min(cond(`t' == 1, year, .))
    replace g_`t' = 0 if missing(g_`t')
}

* Labels for the tables
label var gb_solar    "Green Bond Financing of In-State Solar"
label var gb_wind     "Green Bond Financing of In-State Wind"
label var is_rps      "Renewable Portfolio Standard"
label var ln_gdp_pc   "ln(Real GDP per Capita)"
label var solar_rps   "Solar Financing x RPS"
label var wind_rps    "Wind Financing x RPS"
label var ln_solar_pc "ln(Solar W per Capita + 1)"
label var ln_wind_pc  "ln(Wind W per Capita + 1)"


* ==============================================================================
* PART 1: TWFE REGRESSIONS
*   (1) treatment only
*   (2) + RPS and ln GDP per capita
*   (3) + state and year fixed effects
*   (4) + treatment x RPS interaction
* ==============================================================================
eststo clear

* Solar
eststo s1: regress ln_solar_pc gb_solar, vce(robust)
eststo s2: regress ln_solar_pc gb_solar is_rps ln_gdp_pc, vce(robust)
eststo s3: regress ln_solar_pc gb_solar is_rps ln_gdp_pc i.state_id i.year, vce(cluster state_id)
eststo s4: regress ln_solar_pc gb_solar is_rps solar_rps ln_gdp_pc i.state_id i.year, vce(cluster state_id)

* Wind
eststo w1: regress ln_wind_pc gb_wind, vce(robust)
eststo w2: regress ln_wind_pc gb_wind is_rps ln_gdp_pc, vce(robust)
eststo w3: regress ln_wind_pc gb_wind is_rps ln_gdp_pc i.state_id i.year, vce(cluster state_id)
eststo w4: regress ln_wind_pc gb_wind is_rps wind_rps ln_gdp_pc i.state_id i.year, vce(cluster state_id)

local tabopts b(%9.3f) se(%9.3f) star(* 0.10 ** 0.05 *** 0.01) label nonumbers ///
    mtitles("(1)" "(2)" "(3)" "(4)") ///
    stats(N r2 r2_a, fmt(%9.0f %9.3f %9.3f) labels("Observations" "R-squared" "Adjusted R-squared")) ///
    indicate("State Fixed Effects = *.state_id" "Year Fixed Effects = *.year", labels("Yes" "No")) ///
    nonotes addnotes("Standard errors in parentheses: robust in (1)-(2), clustered by state in (3)-(4)." ///
                     "* p<0.10, ** p<0.05, *** p<0.01")

esttab s1 s2 s3 s4 using "$DESKTOP/Table_TWFE_Solar.rtf", replace `tabopts' ///
    keep(gb_solar is_rps solar_rps ln_gdp_pc _cons) ///
    order(gb_solar is_rps solar_rps ln_gdp_pc _cons) ///
    title("TWFE: Green Bond Financing of Solar on ln(Solar W per Capita + 1)")

esttab w1 w2 w3 w4 using "$DESKTOP/Table_TWFE_Wind.rtf", replace `tabopts' ///
    keep(gb_wind is_rps wind_rps ln_gdp_pc _cons) ///
    order(gb_wind is_rps wind_rps ln_gdp_pc _cons) ///
    title("TWFE: Green Bond Financing of Wind on ln(Wind W per Capita + 1)")

* Same tables in the Results window / log
esttab s1 s2 s3 s4, `tabopts' keep(gb_solar is_rps solar_rps ln_gdp_pc _cons) title("TWFE Solar")
esttab w1 w2 w3 w4, `tabopts' keep(gb_wind is_rps wind_rps ln_gdp_pc _cons) title("TWFE Wind")


* ==============================================================================
* PART 2: CALLAWAY & SANT'ANNA + HONESTDID
* ==============================================================================

* --- Mata helper: breakdown value M-bar from the HonestDiD confidence sets.
*     M-bar = largest M at which the robust 95% CI still excludes zero.
*     Returns -1 if the CI already includes zero at M = 0.
capture mata: mata drop mbar_calc()
mata:
real scalar mbar_calc(real matrix CI, real scalar sgn)
{
    real scalar j, mcol, i, lo, hi, best, k
    real colvector others
    // The M column is the one with a missing value (the "Original" row)
    mcol = 0
    for (j = 1; j <= cols(CI); j++) {
        if (missing(CI[., j]) > 0) mcol = j
    }
    if (mcol == 0) mcol = 1
    best = -1
    for (i = 1; i <= rows(CI); i++) {
        if (missing(CI[i, mcol])) continue
        others = J(0, 1, .)
        for (k = 1; k <= cols(CI); k++) {
            if (k != mcol) others = others \ CI[i, k]
        }
        lo = min(others)
        hi = max(others)
        if ((sgn > 0 & lo > 0) | (sgn < 0 & hi < 0)) {
            if (CI[i, mcol] > best) best = CI[i, mcol]
        }
        else break
    }
    return(best)
}
end

* --- Results file
capture postclose csres
tempfile csresults
postfile csres str10 outcome double n_treated att se pval post_avg ///
    pre_chi2 pre_df pre_p mbar using "`csresults'", replace

* --- One Callaway & Sant'Anna run: ATT, event study, pre-trend test, HonestDiD
capture program drop run_cs
program define run_cs
    args y g stub label

    display _newline as result "=================================================="
    display as result "Callaway & Sant'Anna: `label'"
    display as result "=================================================="

    quietly levelsof state_id if `g' > 0, local(tr)
    local ntr : word count `tr'

    csdid `y' ln_gdp_pc is_rps, ivar(state_id) time(year) gvar(`g') ///
        notyet long2 method(dripw)
    estimates store cs_`stub'

    * Overall ATT
    estat simple, post
    local att = _b[ATT]
    local se  = _se[ATT]
    local pv  = 2 * normal(-abs(`att' / `se'))

    * Event study, full window -> plot
    estimates restore cs_`stub'
    estat event
    csdid_plot, title("Green Bonds Allocated to `label'") ///
        ytitle("Effect on ln(`label' W per Person + 1)") ///
        xtitle("Years Relative to Treatment")
    graph export "$DESKTOP/cs_event_`stub'.png", replace width(2000)

    * Event study used for the pre-trend test and HonestDiD.
    * Leads limited to t = -6 to -2: longer leads rest on only one or two states.
    estimates restore cs_`stub'
    estat event, window(-5 8) estore(ev_`stub')
    estimates restore ev_`stub'
    local post_avg = _b[Post_avg]

    local names : colnames e(b)
    local leads ""
    local pre ""
    local post ""
    local i = 0
    foreach n of local names {
        local ++i
        if substr("`n'", 1, 2) == "Tm" & "`n'" != "Tm1" {
            local leads "`leads' `n'"
            local pre "`pre' `i'"
        }
        if substr("`n'", 1, 2) == "Tp" local post "`post' `i'"
    }

    * Pre-trend test: joint Wald test that all leads equal zero
    test `leads'
    local chi = r(chi2)
    local df  = r(df)
    local pp  = r(p)

    * HonestDiD (relative magnitudes) for the average post-treatment effect
    local pre_first  : word 1 of `pre'
    local pre_last   : word `: word count `pre'' of `pre'
    local post_first : word 1 of `post'
    local post_last  : word `: word count `post'' of `post'
    local npost : word count `post'
    matrix lvec = J(1, `npost', 1 / `npost')

    * Plot (M = 0 to 2 in steps of 0.25)
    honestdid, pre(`pre_first'/`pre_last') post(`post_first'/`post_last') ///
        l_vec(lvec) mvec(0(0.25)2) delta(rm) coefplot ///
        title("HonestDiD Sensitivity: `label'") ///
        xtitle("M (Relative Magnitude of Violation)") ///
        ytitle("95% Robust CI, Average Post-Treatment Effect")
    graph export "$DESKTOP/cs_honestdid_`stub'.png", replace width(2000)

    * Breakdown value M-bar on a finer grid (steps of 0.05)
    local mbar = .
    capture noisily {
        honestdid, pre(`pre_first'/`pre_last') post(`post_first'/`post_last') ///
            l_vec(lvec) mvec(0(0.05)2) delta(rm)
        local hes `s(HonestEventStudy)'
        mata: st_local("mbar", strofreal(mbar_calc(`hes'.CI, `post_avg')))
    }
    if missing(`mbar') {
        display as error "M-bar could not be read from honestdid; use the table above."
    }
    else if `mbar' < 0 {
        display as result "M-bar (`label'): CI includes zero even at M = 0"
    }
    else {
        display as result "M-bar (`label') = `mbar'"
    }

    post csres ("`label'") (`ntr') (`att') (`se') (`pv') (`post_avg') ///
        (`chi') (`df') (`pp') (`mbar')
end

* --- Runs
run_cs ln_solar_pc g_gb_solar solar "Solar"
run_cs ln_wind_pc  g_gb_wind  wind  "Wind"

postclose csres

* --- Results table
preserve
use "`csresults'", clear

gen stars = cond(pval < 0.01, "***", cond(pval < 0.05, "**", cond(pval < 0.10, "*", "")))
gen str20 ATT_str     = string(att, "%9.3f") + stars
gen str20 SE_str      = "(" + string(se, "%9.3f") + ")"
gen str20 Pct_str     = string(100 * (exp(att) - 1), "%9.1f") + "%"
gen str20 P_str       = string(pval, "%9.3f")
gen str30 Pre_str     = string(pre_p, "%9.3f") + " (chi2(" + string(pre_df, "%9.0f") + ") = " + ///
                        string(pre_chi2, "%9.2f") + ")"
gen str30 Mbar_str    = string(mbar, "%9.2f")
replace Mbar_str      = "None (CI includes 0 at M = 0)" if mbar < 0
replace Mbar_str      = "> 2" if mbar >= 2 & !missing(mbar)
replace Mbar_str      = "See log" if missing(mbar)
gen str5 Treated_str  = string(n_treated, "%9.0f")

export delimited outcome n_treated att se pval pre_chi2 pre_df pre_p mbar ///
    using "$DESKTOP/CSDiD_Results.csv", replace

list outcome ATT_str SE_str Pct_str P_str Pre_str Mbar_str Treated_str, noobs abbrev(12)

capture putdocx clear
putdocx begin, font("Times New Roman", 10)
putdocx paragraph, halign(center)
putdocx text ("Callaway & Sant'Anna (2021) Estimates of Green Bond Financing on Renewable Capacity"), ///
    bold font("Times New Roman", 13)

putdocx table tbl = data(outcome ATT_str SE_str Pct_str P_str Pre_str Mbar_str Treated_str), ///
    varnames border(all, nil)
putdocx table tbl(1, 1) = ("Outcome")
putdocx table tbl(1, 2) = ("ATT")
putdocx table tbl(1, 3) = ("SE")
putdocx table tbl(1, 4) = ("% Effect")
putdocx table tbl(1, 5) = ("p-value")
putdocx table tbl(1, 6) = ("Pre-Trend p-value")
putdocx table tbl(1, 7) = ("Breakdown M-bar")
putdocx table tbl(1, 8) = ("Treated States")
putdocx table tbl(1, .), bold border(top, single)
putdocx table tbl(1, .), border(bottom, single)
putdocx table tbl(3, .), border(bottom, single)

putdocx paragraph
putdocx text ("Notes: Outcome is ln(installed W per person + 1). Controls: ln real GDP per capita and RPS status (doubly robust estimator); not-yet-treated states as controls. "), ///
    font("Times New Roman", 9)
putdocx text ("ATT is the simple aggregation across cohorts and post-treatment years; % effect = 100 x (exp(ATT) - 1). Standard errors clustered by state. "), ///
    font("Times New Roman", 9)
putdocx text ("Pre-trend p-value: joint Wald test that event-time leads t = -6 to -2 equal zero (base t = -1). "), ///
    font("Times New Roman", 9)
putdocx text ("Breakdown M-bar: largest relative-magnitude violation of parallel trends (Rambachan & Roth, 2023) at which the 95% robust CI for the average post-treatment effect excludes zero. "), ///
    font("Times New Roman", 9)
putdocx text ("* p<0.10, ** p<0.05, *** p<0.01."), font("Times New Roman", 9)

putdocx save "$DESKTOP/CSDiD_Results.docx", replace
restore

log close