* ==============================================================================
* CALLAWAY & SANT'ANNA: GREEN BONDS ALLOCATED TO WIND AND SOLAR
* Outcome:   ln(W per person + 1)
* Controls:  ln(real GDP per capita), RPS
* Plots saved to: "CSDiD & Sensitivity" folder in the repository
* ==============================================================================
clear all
set more off
capture log close
log using "JoE_CSDiD_improved.log", replace

* Install packages if needed
foreach pkg in drdid csdid honestdid {
    capture which `pkg'
    if _rc ssc install `pkg', replace
}

* Folder where all graphs are saved (created if it doesn't exist)
local outdir "CSDiD & Sensitivity"
capture mkdir "`outdir'"

* Load data
import delimited "panel_improved.csv", clear varnames(1) case(lower)
xtset state_id year

* First year each state was treated (0 = never treated)
foreach t in gb_wind gb_solar {
    bysort state_id: egen g_`t' = min(cond(`t' == 1, year, .))
    replace g_`t' = 0 if missing(g_`t')
}

* ------------------------------------------------------------------------------
* HonestDiD helper: finds the pre- and post-period coefficients in a stored
* event study and plots robust confidence intervals for the AVERAGE
* post-treatment effect as parallel-trends violations grow (M = 0 to 2)
* ------------------------------------------------------------------------------
capture program drop hd_plot
program define hd_plot
    args est title file
    estimates restore `est'
    local names : colnames e(b)
    local pre ""
    local post ""
    local i = 0
    foreach n of local names {
        local ++i
        if substr("`n'", 1, 2) == "Tm" & "`n'" != "Tm1" local pre  "`pre' `i'"
        if substr("`n'", 1, 2) == "Tp" local post "`post' `i'"
    }
    local pre_first  : word 1 of `pre'
    local pre_last   : word `: word count `pre'' of `pre'
    local post_first : word 1 of `post'
    local post_last  : word `: word count `post'' of `post'
    local npost : word count `post'
    matrix lvec = J(1, `npost', 1 / `npost')

    honestdid, pre(`pre_first'/`pre_last') post(`post_first'/`post_last') ///
        l_vec(lvec) mvec(0(0.25)2) delta(rm) coefplot ///
        title("`title'") xtitle("M (Relative Magnitude of Violation)") ///
        ytitle("95% Robust CI, Average Post-Treatment Effect")
    graph export "`file'", replace width(2000)
end

* ------------------------------------------------------------------------------
* WIND: green bonds allocated to wind -> wind capacity per person
* ------------------------------------------------------------------------------
csdid ln_wind_pc ln_gdp_pc is_rps, ivar(state_id) time(year) gvar(g_gb_wind) ///
    notyet long2 method(dripw)
estimates store cs_wind

estat simple                      // <- THE WIND EFFECT (ATT)

* Event study, full window
estimates restore cs_wind
estat event
csdid_plot, title("Green Bonds Allocated to Wind") ///
    ytitle("Effect on ln(Wind W per Person + 1)") xtitle("Years Relative to Treatment")
graph export "`outdir'/cs_event_wind.png", replace width(2000)

* HonestDiD (uses leads -5 to -2: the far leads rest on 1-2 states and would
* make the sensitivity bounds reflect noise rather than real pre-trends)
estimates restore cs_wind
estat event, window(-5 8) estore(ev_wind)
hd_plot ev_wind "HonestDiD Sensitivity: Wind" "`outdir'/cs_honestdid_wind.png"

* ------------------------------------------------------------------------------
* SOLAR: green bonds allocated to solar -> solar capacity per person
* ------------------------------------------------------------------------------
csdid ln_solar_pc ln_gdp_pc is_rps, ivar(state_id) time(year) gvar(g_gb_solar) ///
    notyet long2 method(dripw)
estimates store cs_solar

estat simple                      // <- THE SOLAR EFFECT (ATT)

* Event study, full window
estimates restore cs_solar
estat event
csdid_plot, title("Green Bonds Allocated to Solar") ///
    ytitle("Effect on ln(Solar W per Person + 1)") xtitle("Years Relative to Treatment")
graph export "`outdir'/cs_event_solar.png", replace width(2000)

* HonestDiD
estimates restore cs_solar
estat event, window(-5 8) estore(ev_solar)
hd_plot ev_solar "HonestDiD Sensitivity: Solar" "`outdir'/cs_honestdid_solar.png"

* ==============================================================================
* ROBUSTNESS: TOTAL CAPACITY (NOT PER CAPITA)
* Outcome: ln(installed MW + 1)
* ==============================================================================
gen ln_wind_mw  = ln(wind_capacity + 1)
gen ln_solar_mw = ln(solar_capacity + 1)
label var ln_wind_mw  "ln(Wind MW + 1)"
label var ln_solar_mw "ln(Solar MW + 1)"

* ------------------------------------------------------------------------------
* WIND: green bonds allocated to wind -> total wind capacity
* ------------------------------------------------------------------------------
csdid ln_wind_mw ln_gdp_pc is_rps, ivar(state_id) time(year) gvar(g_gb_wind) ///
    notyet long2 method(dripw)
estimates store cs_wind_mw

estat simple                      // <- WIND EFFECT ON TOTAL MW (ATT)

estimates restore cs_wind_mw
estat event
csdid_plot, title("Green Bonds Allocated to Wind: Total Capacity") ///
    ytitle("Effect on ln(Wind MW + 1)") xtitle("Years Relative to Treatment")
graph export "`outdir'/cs_event_wind_mw.png", replace width(2000)

estimates restore cs_wind_mw
estat event, window(-5 8) estore(ev_wind_mw)
hd_plot ev_wind_mw "HonestDiD Sensitivity: Wind (Total MW)" "`outdir'/cs_honestdid_wind_mw.png"

* ------------------------------------------------------------------------------
* SOLAR: green bonds allocated to solar -> total solar capacity
* ------------------------------------------------------------------------------
csdid ln_solar_mw ln_gdp_pc is_rps, ivar(state_id) time(year) gvar(g_gb_solar) ///
    notyet long2 method(dripw)
estimates store cs_solar_mw

estat simple                      // <- SOLAR EFFECT ON TOTAL MW (ATT)

estimates restore cs_solar_mw
estat event
csdid_plot, title("Green Bonds Allocated to Solar: Total Capacity") ///
    ytitle("Effect on ln(Solar MW + 1)") xtitle("Years Relative to Treatment")
graph export "`outdir'/cs_event_solar_mw.png", replace width(2000)

estimates restore cs_solar_mw
estat event, window(-5 8) estore(ev_solar_mw)
hd_plot ev_solar_mw "HonestDiD Sensitivity: Solar (Total MW)" "`outdir'/cs_honestdid_solar_mw.png"

log close