clear
set more off
capture log close
log using JoE_Stata.log, replace
import delimited "panel.csv", clear


*** Begin Stata Coding Here
summarize

save "JoE_Stata.dta", replace