* ==============================================================================
* TWFE REGRESSIONS ON panel_improved.csv
* Green Bond Financing and State Renewable Energy Capacity (2010-2024)
* Author: Noah Leibowitz
*
* Outcomes:  ln(W per person + 1) for solar and wind  (ln_solar_pc, ln_wind_pc)
* Control:   ln(real GDP per capita)                    (ln_gdp_pc), is_rps
*
* Tables (written as .rtf, open in Word):
*   Table1_Solar_GB          corrected market-entry gb on solar
*   Table2_Solar_Ownership   corrected corporate / public on solar
*   Table3_Wind_GB           corrected market-entry gb on wind
*   Table4_Wind_Ownership    corrected corporate / public on wind
*   Table5_Tech_Specific     gb_solar on solar, gb_wind on wind (main spec)
*   Table6_Placebo           gb_other (non-solar/wind bonds) placebo tests
* ==============================================================================

clear all
set more off
capture log close
log using "JoE_TWFE_improved.log", replace

capture which esttab
if _rc ssc install estout, replace

* ------------------------------------------------------------------------------
* 1. LOAD DATA & SETUP
* ------------------------------------------------------------------------------
import delimited "panel_improved.csv", clear varnames(1) case(lower)

xtset state_id year

* Interaction terms (treatment x active RPS)
gen gb_rps        = gb * is_rps
gen corporate_rps = is_corporate * is_rps
gen public_rps    = is_public * is_rps
gen solar_rps     = gb_solar * is_rps
gen wind_rps      = gb_wind * is_rps

* Labels used in the output tables
label var gb            "Green Bond Issuer (gb)"
label var is_corporate  "Corporate Green Bond Issuer"
label var is_public     "Public Green Bond Issuer"
label var gb_solar      "Green Bond Financing of In-State Solar"
label var gb_wind       "Green Bond Financing of In-State Wind"
label var gb_other      "Placebo: Non-Solar/Wind Green Bond"
label var is_rps        "Renewable Portfolio Standard"
label var ln_gdp_pc     "ln(Real GDP per Capita)"
label var gb_rps        "Green Bond x RPS"
label var corporate_rps "Corporate x RPS"
label var public_rps    "Public x RPS"
label var solar_rps     "Solar Financing x RPS"
label var wind_rps      "Wind Financing x RPS"
label var ln_solar_pc   "ln(Solar W per Capita + 1)"
label var ln_wind_pc    "ln(Wind W per Capita + 1)"

* Quick checks
assert _N == 720
foreach v in gb is_corporate is_public gb_solar gb_wind gb_other is_rps placebo_sample {
    assert inlist(`v', 0, 1)
}
summarize ln_solar_pc ln_wind_pc ln_gdp_pc gb is_corporate is_public gb_solar gb_wind gb_other is_rps

* Common esttab options
local tabopts b(%9.3f) se(%9.3f) star(* 0.10 ** 0.05 *** 0.01) label ///
    stats(N r2 r2_a, fmt(%9.0f %9.3f %9.3f) labels("Observations" "R-squared" "Adjusted R-squared")) ///
    indicate("State Fixed Effects = *.state_id" "Year Fixed Effects = *.year", labels("Yes" "No")) ///
    nonotes addnotes("Standard errors in parentheses: robust in pooled models, clustered by state in fixed-effects models." ///
                     "* p<0.10, ** p<0.05, *** p<0.01")

* ------------------------------------------------------------------------------
* 2. TABLES 1-4: CORRECTED MARKET-ENTRY TREATMENT (comparable to the draft)
* ------------------------------------------------------------------------------
eststo clear

* Solar: any green bond
eststo m1:  regress ln_solar_pc gb, vce(robust)
eststo m2:  regress ln_solar_pc gb is_rps ln_gdp_pc, vce(robust)
eststo m3:  regress ln_solar_pc gb is_rps ln_gdp_pc i.state_id i.year, vce(cluster state_id)
eststo m4:  regress ln_solar_pc gb is_rps gb_rps ln_gdp_pc i.state_id i.year, vce(cluster state_id)

* Solar: by issuer type
eststo m5:  regress ln_solar_pc is_corporate is_public, vce(robust)
eststo m6:  regress ln_solar_pc is_corporate is_public is_rps ln_gdp_pc, vce(robust)
eststo m7:  regress ln_solar_pc is_corporate is_public is_rps ln_gdp_pc i.state_id i.year, vce(cluster state_id)
eststo m8:  regress ln_solar_pc is_corporate is_public is_rps corporate_rps public_rps ln_gdp_pc i.state_id i.year, vce(cluster state_id)

* Wind: any green bond
eststo m9:  regress ln_wind_pc gb, vce(robust)
eststo m10: regress ln_wind_pc gb is_rps ln_gdp_pc, vce(robust)
eststo m11: regress ln_wind_pc gb is_rps ln_gdp_pc i.state_id i.year, vce(cluster state_id)
eststo m12: regress ln_wind_pc gb is_rps gb_rps ln_gdp_pc i.state_id i.year, vce(cluster state_id)

* Wind: by issuer type
eststo m13: regress ln_wind_pc is_corporate is_public, vce(robust)
eststo m14: regress ln_wind_pc is_corporate is_public is_rps ln_gdp_pc, vce(robust)
eststo m15: regress ln_wind_pc is_corporate is_public is_rps ln_gdp_pc i.state_id i.year, vce(cluster state_id)
eststo m16: regress ln_wind_pc is_corporate is_public is_rps corporate_rps public_rps ln_gdp_pc i.state_id i.year, vce(cluster state_id)

esttab m1 m2 m3 m4 using "Table1_Solar_GB.rtf", replace `tabopts' ///
    keep(gb is_rps gb_rps ln_gdp_pc _cons) order(gb is_rps gb_rps ln_gdp_pc _cons) ///
    title("Table 1: Green Bonds on Solar Capacity per Capita (Corrected Market Entry)") ///
    mtitles("(1)" "(2)" "(3)" "(4)") nonumbers

esttab m5 m6 m7 m8 using "Table2_Solar_Ownership.rtf", replace `tabopts' ///
    keep(is_corporate is_public is_rps corporate_rps public_rps ln_gdp_pc _cons) ///
    order(is_corporate is_public is_rps corporate_rps public_rps ln_gdp_pc _cons) ///
    title("Table 2: Solar Capacity per Capita by Issuer Type (Corrected)") ///
    mtitles("(5)" "(6)" "(7)" "(8)") nonumbers

esttab m9 m10 m11 m12 using "Table3_Wind_GB.rtf", replace `tabopts' ///
    keep(gb is_rps gb_rps ln_gdp_pc _cons) order(gb is_rps gb_rps ln_gdp_pc _cons) ///
    title("Table 3: Green Bonds on Wind Capacity per Capita (Corrected Market Entry)") ///
    mtitles("(9)" "(10)" "(11)" "(12)") nonumbers

esttab m13 m14 m15 m16 using "Table4_Wind_Ownership.rtf", replace `tabopts' ///
    keep(is_corporate is_public is_rps corporate_rps public_rps ln_gdp_pc _cons) ///
    order(is_corporate is_public is_rps corporate_rps public_rps ln_gdp_pc _cons) ///
    title("Table 4: Wind Capacity per Capita by Issuer Type (Corrected)") ///
    mtitles("(13)" "(14)" "(15)" "(16)") nonumbers

* ------------------------------------------------------------------------------
* 3. TABLE 5: TECHNOLOGY-SPECIFIC TREATMENT (main specification)
*    gb_solar enters only the solar equation, gb_wind only the wind equation
* ------------------------------------------------------------------------------
eststo t1: regress ln_solar_pc gb_solar is_rps ln_gdp_pc, vce(robust)
eststo t2: regress ln_solar_pc gb_solar is_rps ln_gdp_pc i.state_id i.year, vce(cluster state_id)
eststo t3: regress ln_solar_pc gb_solar is_rps solar_rps ln_gdp_pc i.state_id i.year, vce(cluster state_id)

eststo t4: regress ln_wind_pc gb_wind is_rps ln_gdp_pc, vce(robust)
eststo t5: regress ln_wind_pc gb_wind is_rps ln_gdp_pc i.state_id i.year, vce(cluster state_id)
eststo t6: regress ln_wind_pc gb_wind is_rps wind_rps ln_gdp_pc i.state_id i.year, vce(cluster state_id)

esttab t1 t2 t3 t4 t5 t6 using "Table5_Tech_Specific.rtf", replace `tabopts' ///
    keep(gb_solar solar_rps gb_wind wind_rps is_rps ln_gdp_pc _cons) ///
    order(gb_solar solar_rps gb_wind wind_rps is_rps ln_gdp_pc _cons) ///
    title("Table 5: Green Bond Financing of In-State Solar and Wind Capacity") ///
    mgroups("ln(Solar W per Capita + 1)" "ln(Wind W per Capita + 1)", pattern(1 0 0 1 0 0)) ///
    mtitles("(1)" "(2)" "(3)" "(4)" "(5)" "(6)") nonumbers

* ------------------------------------------------------------------------------
* 4. TABLE 6: PLACEBO
*    (1)-(2): gb_other among states never treated for solar or wind
*    (3)-(4): full sample, technology treatment and placebo together
* ------------------------------------------------------------------------------
eststo p1: regress ln_solar_pc gb_other is_rps ln_gdp_pc i.state_id i.year ///
    if placebo_sample == 1, vce(cluster state_id)
eststo p2: regress ln_wind_pc gb_other is_rps ln_gdp_pc i.state_id i.year ///
    if placebo_sample == 1, vce(cluster state_id)
eststo p3: regress ln_solar_pc gb_solar gb_other is_rps ln_gdp_pc i.state_id i.year, vce(cluster state_id)
eststo p4: regress ln_wind_pc gb_wind gb_other is_rps ln_gdp_pc i.state_id i.year, vce(cluster state_id)

esttab p1 p2 p3 p4 using "Table6_Placebo.rtf", replace `tabopts' ///
    keep(gb_other gb_solar gb_wind is_rps ln_gdp_pc _cons) ///
    order(gb_other gb_solar gb_wind is_rps ln_gdp_pc _cons) ///
    title("Table 6: Placebo Tests Using Green Bonds That Funded Non-Solar, Non-Wind Uses") ///
    mgroups("Placebo sample" "Full sample", pattern(1 0 1 0)) ///
    mtitles("Solar" "Wind" "Solar" "Wind")

* Also print every table to the Results window / log
esttab m1 m2 m3 m4, `tabopts' keep(gb is_rps gb_rps ln_gdp_pc _cons) title("Table 1")
esttab m5 m6 m7 m8, `tabopts' keep(is_corporate is_public is_rps corporate_rps public_rps ln_gdp_pc _cons) title("Table 2")
esttab m9 m10 m11 m12, `tabopts' keep(gb is_rps gb_rps ln_gdp_pc _cons) title("Table 3")
esttab m13 m14 m15 m16, `tabopts' keep(is_corporate is_public is_rps corporate_rps public_rps ln_gdp_pc _cons) title("Table 4")
esttab t1 t2 t3 t4 t5 t6, `tabopts' keep(gb_solar solar_rps gb_wind wind_rps is_rps ln_gdp_pc _cons) title("Table 5")
esttab p1 p2 p3 p4, `tabopts' keep(gb_other gb_solar gb_wind is_rps ln_gdp_pc _cons) title("Table 6")

log close
