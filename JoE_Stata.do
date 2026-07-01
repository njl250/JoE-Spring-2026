clear
set more off
capture log close
log using JoE_Stata.log, replace
import delimited "panel.csv", clear


*** Begin Stata Coding Here
summarize
describe

*** Initialize the panel data structure
encode state, gen(state_id)
xtset state_id year

*** Convert to raw Watts per capita so the numbers are larger than 1
gen watts_wind_pc = wind_capacity_per_capita * 1000000
gen watts_solar_pc = solar_capacity_per_capita * 1000000

*** Natural log our wind capacity per capita, solar capacity per capita, and GDP per capita variables
*** To avoid generating missing values with ln(0), add 1 to the capacity (ln(x + 1))
gen ln_wind_capacity = ln(watts_wind_pc + 1)
gen ln_solar_capacity = ln(watts_solar_pc + 1)

*** GDP per capita is always > 0
gen ln_gdp_per_capita = ln(real_gdp_per_capita)

*** RPS x Corporate Interaction Term
gen rps_x_corporate = is_rps * is_corporate

*** RPS x Public Interaction Term
gen rps_x_public = is_rps * is_public

*** RPS x Any Green Bond Interaction Term
gen rps_x_anyGB = is_rps * gb

summarize
describe


save "JoE_Stata_1.dta", replace