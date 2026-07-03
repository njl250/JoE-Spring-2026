clear
set more off
capture log close
log using JoE_Stata.log, replace
import delimited "panel.csv", clear

* ==============================================================================
* DATA PREPARATION & VARIABLE GENERATION
* ==============================================================================

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
gen rps_x_corporate = is_rps * is_corporate
gen rps_x_public    = is_rps * is_public
gen rps_x_anyGB     = is_rps * gb

*** Install/update the core Callaway & Sant'Anna packages
ssc install drdid, replace
ssc install csdid, replace

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
* PREPARE POSTFILE FOR RESULTS TABLE
* ==============================================================================
tempname memhold
postfile `memhold' str40 execution str15 outcome str15 model_type double pre_wald_p double post_att double post_se using "csdid_summary_table.dta", replace


* ==============================================================================
* BLOCK 1: BASELINE EXECUTIONS (CONTROL: GDP ONLY)
* ==============================================================================

*** Execution 1: Any GB on ln(Solar) [Baseline]
csdid ln_solar_capacity ln_gdp_per_capita, ivar(state_id) time(year) gvar(first_gb_year) method(dripw)
estat event
csdid_plot, title("Any GB on ln(Solar) [Baseline]") name(plot1_base, replace)
graph export "plot1_solar_any_base.png", replace
estat event, post
capture test Tm10 Tm9 Tm8 Tm7 Tm6 Tm5 Tm4 Tm3 Tm2 Tm1
local p_val = r(p)
lincom Post_avg
post `memhold' ("Any GB") ("Solar") ("Baseline") (`p_val') (r(estimate)) (r(se))

*** Execution 2A: Corporate GB on ln(Solar) [Baseline]
csdid ln_solar_capacity ln_gdp_per_capita, ivar(state_id) time(year) gvar(first_corp_year) method(dripw)
estat event
csdid_plot, title("Corporate GB on ln(Solar) [Baseline]") name(plot2a_base, replace)
graph export "plot2_solar_corp_base.png", replace
estat event, post
capture test Tm10 Tm9 Tm8 Tm7 Tm6 Tm5 Tm4 Tm3 Tm2 Tm1
local p_val = r(p)
lincom Post_avg
post `memhold' ("Corporate GB") ("Solar") ("Baseline") (`p_val') (r(estimate)) (r(se))

*** Execution 2B: Public GB on ln(Solar) [Baseline]
csdid ln_solar_capacity ln_gdp_per_capita, ivar(state_id) time(year) gvar(first_muni_year) method(dripw)
estat event
csdid_plot, title("Public GB on ln(Solar) [Baseline]") name(plot2b_base, replace)
graph export "plot2_solar_muni_base.png", replace
estat event, post
capture test Tm10 Tm9 Tm8 Tm7 Tm6 Tm5 Tm4 Tm3 Tm2 Tm1
local p_val = r(p)
lincom Post_avg
post `memhold' ("Public GB") ("Solar") ("Baseline") (`p_val') (r(estimate)) (r(se))

*** Execution 3: Any GB on ln(Wind) [Baseline]
csdid ln_wind_capacity ln_gdp_per_capita, ivar(state_id) time(year) gvar(first_gb_year) method(dripw)
estat event
csdid_plot, title("Any GB on ln(Wind) [Baseline]") name(plot3_base, replace)
graph export "plot3_wind_any_base.png", replace
estat event, post
capture test Tm10 Tm9 Tm8 Tm7 Tm6 Tm5 Tm4 Tm3 Tm2 Tm1
local p_val = r(p)
lincom Post_avg
post `memhold' ("Any GB") ("Wind") ("Baseline") (`p_val') (r(estimate)) (r(se))

*** Execution 4A: Corporate GB on ln(Wind) [Baseline]
csdid ln_wind_capacity ln_gdp_per_capita, ivar(state_id) time(year) gvar(first_corp_year) method(dripw)
estat event
csdid_plot, title("Corporate GB on ln(Wind) [Baseline]") name(plot4a_base, replace)
graph export "plot4_wind_corp_base.png", replace
estat event, post
capture test Tm10 Tm9 Tm8 Tm7 Tm6 Tm5 Tm4 Tm3 Tm2 Tm1
local p_val = r(p)
lincom Post_avg
post `memhold' ("Corporate GB") ("Wind") ("Baseline") (`p_val') (r(estimate)) (r(se))

*** Execution 4B: Public GB on ln(Wind) [Baseline]
csdid ln_wind_capacity ln_gdp_per_capita, ivar(state_id) time(year) gvar(first_muni_year) method(dripw)
estat event
csdid_plot, title("Public GB on ln(Wind) [Baseline]") name(plot4b_base, replace)
graph export "plot4_wind_muni_base.png", replace
estat event, post
capture test Tm10 Tm9 Tm8 Tm7 Tm6 Tm5 Tm4 Tm3 Tm2 Tm1
local p_val = r(p)
lincom Post_avg
post `memhold' ("Public GB") ("Wind") ("Baseline") (`p_val') (r(estimate)) (r(se))


* ==============================================================================
* BLOCK 2: POLICY-ROBUST EXECUTIONS (CONTROLS: GDP + IS_RPS)
* ==============================================================================

*** Execution 1: Any GB on ln(Solar) [RPS Control]
csdid ln_solar_capacity ln_gdp_per_capita is_rps, ivar(state_id) time(year) gvar(first_gb_year) method(dripw)
estat event
csdid_plot, title("Any GB on ln(Solar) [RPS Control]") name(plot1_rps, replace)
graph export "plot1_solar_any_rps.png", replace
estat event, post
capture test Tm10 Tm9 Tm8 Tm7 Tm6 Tm5 Tm4 Tm3 Tm2 Tm1
local p_val = r(p)
lincom Post_avg
post `memhold' ("Any GB") ("Solar") ("RPS Control") (`p_val') (r(estimate)) (r(se))

*** Execution 2A: Corporate GB on ln(Solar) [RPS Control]
csdid ln_solar_capacity ln_gdp_per_capita is_rps, ivar(state_id) time(year) gvar(first_corp_year) method(dripw)
estat event
csdid_plot, title("Corporate GB on ln(Solar) [RPS Control]") name(plot2a_rps, replace)
graph export "plot2_solar_corp_rps.png", replace
estat event, post
capture test Tm10 Tm9 Tm8 Tm7 Tm6 Tm5 Tm4 Tm3 Tm2 Tm1
local p_val = r(p)
lincom Post_avg
post `memhold' ("Corporate GB") ("Solar") ("RPS Control") (`p_val') (r(estimate)) (r(se))

*** Execution 2B: Public GB on ln(Solar) [RPS Control]
csdid ln_solar_capacity ln_gdp_per_capita is_rps, ivar(state_id) time(year) gvar(first_muni_year) method(dripw)
estat event
csdid_plot, title("Public GB on ln(Solar) [RPS Control]") name(plot2b_rps, replace)
graph export "plot2_solar_muni_rps.png", replace
estat event, post
capture test Tm10 Tm9 Tm8 Tm7 Tm6 Tm5 Tm4 Tm3 Tm2 Tm1
local p_val = r(p)
lincom Post_avg
post `memhold' ("Public GB") ("Solar") ("RPS Control") (`p_val') (r(estimate)) (r(se))

*** Execution 3: Any GB on ln(Wind) [RPS Control]
csdid ln_wind_capacity ln_gdp_per_capita is_rps, ivar(state_id) time(year) gvar(first_gb_year) method(dripw)
estat event
csdid_plot, title("Any GB on ln(Wind) [RPS Control]") name(plot3_rps, replace)
graph export "plot3_wind_any_rps.png", replace
estat event, post
capture test Tm10 Tm9 Tm8 Tm7 Tm6 Tm5 Tm4 Tm3 Tm2 Tm1
local p_val = r(p)
lincom Post_avg
post `memhold' ("Any GB") ("Wind") ("RPS Control") (`p_val') (r(estimate)) (r(se))

*** Execution 4A: Corporate GB on ln(Wind) [RPS Control]
csdid ln_wind_capacity ln_gdp_per_capita is_rps, ivar(state_id) time(year) gvar(first_corp_year) method(dripw)
estat event
csdid_plot, title("Corporate GB on ln(Wind) [RPS Control]") name(plot4a_rps, replace)
graph export "plot4_wind_corp_rps.png", replace
estat event, post
capture test Tm10 Tm9 Tm8 Tm7 Tm6 Tm5 Tm4 Tm3 Tm2 Tm1
local p_val = r(p)
lincom Post_avg
post `memhold' ("Corporate GB") ("Wind") ("RPS Control") (`p_val') (r(estimate)) (r(se))

*** Execution 4B: Public GB on ln(Wind) [RPS Control]
csdid ln_wind_capacity ln_gdp_per_capita is_rps, ivar(state_id) time(year) gvar(first_muni_year) method(dripw)
estat event
csdid_plot, title("Public GB on ln(Wind) [RPS Control]") name(plot4b_rps, replace)
graph export "plot4_wind_muni_rps.png", replace
estat event, post
capture test Tm10 Tm9 Tm8 Tm7 Tm6 Tm5 Tm4 Tm3 Tm2 Tm1
local p_val = r(p)
lincom Post_avg
post `memhold' ("Public GB") ("Wind") ("RPS Control") (`p_val') (r(estimate)) (r(se))

* Close the collection file
postclose `memhold'


* ==============================================================================
* DISPLAY SUMMARY TABLE IN RESULTS WINDOW
* ==============================================================================
use "csdid_summary_table.dta", clear
gen t_stat = post_att / post_se
format pre_wald_p post_att post_se t_stat %6.4f
list execution outcome model_type pre_wald_p post_att t_stat, clean


* ==============================================================================
* GENERATE FORMATTED MATRIX DATASET
* ==============================================================================
use "csdid_summary_table.dta", clear

* Create a unique column ID based on Model Type & Execution
gen model_id = 1 if execution == "Any GB"        & model_type == "Baseline"
replace model_id = 2 if execution == "Any GB"        & model_type == "RPS Control"
replace model_id = 3 if execution == "Corporate GB"  & model_type == "Baseline"
replace model_id = 4 if execution == "Corporate GB"  & model_type == "RPS Control"
replace model_id = 5 if execution == "Public GB"     & model_type == "Baseline"
replace model_id = 6 if execution == "Public GB"     & model_type == "RPS Control"

* Make a clean label for the rows
gen type = "beta"
expand 2, gen(dup)
replace type = "se" if dup == 1
gen sort_var = cond(type == "beta", 1, 2)

* Create the value that goes into the cell (either the ATT or the SE)
gen cell_value = post_att if type == "beta"
replace cell_value = post_se if type == "se"

* Generate a display variable with significance stars applied directly
gen str30 display_val = string(cell_value, "%6.4f")
gen t_stat = post_att / post_se

replace display_val = display_val + "***" if type == "beta" & abs(t_stat) >= 2.576
replace display_val = display_val + "**"  if type == "beta" & abs(t_stat) >= 1.960 & abs(t_stat) < 2.576
replace display_val = display_val + "*"   if type == "beta" & abs(t_stat) >= 1.645 & abs(t_stat) < 1.960
replace display_val = "(" + display_val + ")" if type == "se"

* Add the Wald p-values as their own row at the bottom of each panel
preserve
    keep if type == "beta"
    replace type = "wald"
    replace sort_var = 3
    replace display_val = string(pre_wald_p, "%6.4f")
    replace display_val = "—" if display_val == "."
    tempfile wald_rows
    save `wald_rows'
restore
append using `wald_rows'

* Reshape wide so models 1 through 6 become columns
keep outcome model_id type sort_var display_val
reshape wide display_val, i(outcome type sort_var) j(model_id)

sort outcome sort_var
gen row_label = "Average Post-Treatment ATT" if type == "beta"
replace row_label = " " if type == "se"
replace row_label = "Pre-Trend Wald p-value" if type == "wald"
order row_label display_val1 display_val2 display_val3 display_val4 display_val5 display_val6


* ==============================================================================
* CLEAN PRINTING TO LOG/RESULTS WINDOW
* ==============================================================================
* Print Panel A
display _newline(2) "{hline}"
display " PANEL A: ln(Solar Capacity)"
display "{hline}"
list row_label display_val1 display_val2 display_val3 display_val4 display_val5 display_val6 if outcome == "Solar", noobs subvarname

* Print Panel B
display _newline(2) "{hline}"
display " PANEL B: ln(Wind Capacity)"
display "{hline}"
list row_label display_val1 display_val2 display_val3 display_val4 display_val5 display_val6 if outcome == "Wind", noobs subvarname


* ==============================================================================
* EXPORT PROFESSIONAL CSDID REGRESSION TABLE TO MS WORD (.DOCX)
* ==============================================================================
use "csdid_summary_table.dta", clear

* 1. Re-map the exact model column indices (1 to 6)
gen model_id = 1 if execution == "Any GB"        & model_type == "Baseline"
replace model_id = 2 if execution == "Any GB"        & model_type == "RPS Control"
replace model_id = 3 if execution == "Corporate GB"  & model_type == "Baseline"
replace model_id = 4 if execution == "Corporate GB"  & model_type == "RPS Control"
replace model_id = 5 if execution == "Public GB"     & model_type == "Baseline"
replace model_id = 6 if execution == "Public GB"     & model_type == "RPS Control"

* 2. Calculate standard errors and structure significance stars
gen t_stat = post_att / post_se

gen str30 b_str = string(post_att, "%6.4f")
replace b_str = b_str + "***" if abs(t_stat) >= 2.576
replace b_str = b_str + "**"  if abs(t_stat) >= 1.960 & abs(t_stat) < 2.576
replace b_str = b_str + "*"   if abs(t_stat) >= 1.645 & abs(t_stat) < 1.960

gen str30 se_str = "(" + string(post_se, "%6.4f") + ")"

gen str30 wald_str = string(pre_wald_p, "%6.4f")
replace wald_str = "—" if wald_str == "."

* 3. Clear Stata's active document stream and initialize document default type
capture putdocx clear
putdocx begin

* Title Paragraph
putdocx paragraph, halign(center) spacing(after, 12pt)
putdocx text ("Green Bond Impact on Renewable Capacity"), bold font("Times New Roman", 14)

* Context Description Paragraph
putdocx paragraph, halign(left) spacing(after, 18pt)
putdocx text ("The table below reports the average post-treatment Average Treatment Effect on the Treated (ATT) across baseline and policy-robust models for solar and wind capacity outcomes using the Callaway and Sant'Anna (2021) difference-in-differences framework."), font("Times New Roman", 11)

* Initialize Table: 11 rows total (Header + Panel A + Panel B + Controls), 7 columns
putdocx table reg_table = (11, 7), halign(center)

* Format all borders blank by default (to simulate strict academic styling)
putdocx table reg_table(.,.), border(all, nil)

* Apply Top Line and Header Line Borders (APA style)
putdocx table reg_table(1,.), border(top, single, "000000", "1.5pt")
putdocx table reg_table(1,.), border(bottom, single, "000000", "0.75pt")

* --- POPULATE HEADERS (Row 1) ---
putdocx table reg_table(1,1) = ("Variable / Specification")
putdocx table reg_table(1,2) = ("(1)" + char(10) + "Any GB" + char(10) + "Baseline")
putdocx table reg_table(1,3) = ("(2)" + char(10) + "Any GB" + char(10) + "RPS Control")
putdocx table reg_table(1,4) = ("(3)" + char(10) + "Corporate" + char(10) + "Baseline")
putdocx table reg_table(1,5) = ("(4)" + char(10) + "Corporate" + char(10) + "RPS Control")
putdocx table reg_table(1,6) = ("(5)" + char(10) + "Public" + char(10) + "Baseline")
putdocx table reg_table(1,7) = ("(6)" + char(10) + "Public" + char(10) + "RPS Control")

* Set text style for headers
forvalues c=1/7 {
    putdocx table reg_table(1,`c'), bold font("Times New Roman", 10)
    if `c' > 1 putdocx table reg_table(1,`c'), halign(center)
}

* --- POPULATE PANEL A: SOLAR (Rows 2 to 5) ---
putdocx table reg_table(2,1) = ("Panel A: ln(Solar Capacity)"), bold font("Times New Roman", 10)
putdocx table reg_table(3,1) = ("Average Post-Treatment ATT"), font("Times New Roman", 10)
putdocx table reg_table(4,1) = (" "), font("Times New Roman", 10)
putdocx table reg_table(5,1) = ("Pre-Trend Wald p-value"), font("Times New Roman", 10)

forvalues c=1/6 {
    * Pull information matching Column indices dynamically
    quietly levelsof b_str if outcome == "Solar" & model_id == `c', local(b_val) clean
    quietly levelsof se_str if outcome == "Solar" & model_id == `c', local(se_val) clean
    quietly levelsof wald_str if outcome == "Solar" & model_id == `c', local(wald_val) clean
    
    * Write data into cells (shifting columns right by 1 to accommodate row titles)
    local target_col = `c' + 1
    putdocx table reg_table(3, `target_col') = ("`b_val'"), halign(center) font("Times New Roman", 10)
    putdocx table reg_table(4, `target_col') = ("`se_val'"), halign(center) font("Times New Roman", 10)
    putdocx table reg_table(5, `target_col') = ("`wald_val'"), halign(center) font("Times New Roman", 10)
}

* --- POPULATE PANEL B: WIND (Rows 6 to 9) ---
putdocx table reg_table(6,1) = ("Panel B: ln(Wind Capacity)"), bold font("Times New Roman", 10)
putdocx table reg_table(6,.), border(top, single, "CCCCCC", "0.5pt") // Panel boundary line
putdocx table reg_table(7,1) = ("Average Post-Treatment ATT"), font("Times New Roman", 10)
putdocx table reg_table(8,1) = (" "), font("Times New Roman", 10)
putdocx table reg_table(9,1) = ("Pre-Trend Wald p-value"), font("Times New Roman", 10)

forvalues c=1/6 {
    quietly levelsof b_str if outcome == "Wind" & model_id == `c', local(b_val) clean
    quietly levelsof se_str if outcome == "Wind" & model_id == `c', local(se_val) clean
    quietly levelsof wald_str if outcome == "Wind" & model_id == `c', local(wald_val) clean
    
    local target_col = `c' + 1
    putdocx table reg_table(7, `target_col') = ("`b_val'"), halign(center) font("Times New Roman", 10)
    putdocx table reg_table(8, `target_col') = ("`se_val'"), halign(center) font("Times New Roman", 10)
    putdocx table reg_table(9, `target_col') = ("`wald_val'"), halign(center) font("Times New Roman", 10)
}

* --- POPULATE CONTROLS MATRIX FOOTER (Rows 10 and 11) ---
putdocx table reg_table(10,.), border(top, single, "000000", "0.75pt") // Divider line before controls

putdocx table reg_table(10,1) = ("Controls: ln(GDP per capita)"), font("Times New Roman", 10)
putdocx table reg_table(11,1) = ("Controls: RPS Mandate"), font("Times New Roman", 10)

forvalues c=1/6 {
    local target_col = `c' + 1
    putdocx table reg_table(10, `target_col') = ("Yes"), halign(center) font("Times New Roman", 10)
    
    * Explicitly write Yes or No based on whether it is an RPS Control specification
    if inlist(`c', 2, 4, 6) putdocx table reg_table(11, `target_col') = ("Yes"), halign(center) font("Times New Roman", 10)
    if inlist(`c', 1, 3, 5) putdocx table reg_table(11, `target_col') = ("No"), halign(center) font("Times New Roman", 10)
}

* Double-line anchor effect at the absolute bottom of the grid
putdocx table reg_table(11,.), border(bottom, single, "000000", "1.5pt")

* --- FOOTNOTES SECTIONS ---
putdocx paragraph, halign(left) spacing(before, 8pt)
putdocx text ("Notes: "), bold italic font("Times New Roman", 9)
putdocx text ("Standard errors are reported in parentheses underneath the average post-treatment effects. Coefficients represent the aggregate dynamic ATT. Multi-collinearity due to dense treatment cohort distributions prevents the estimation of pre-trend metrics for aggregate specifications in columns (1) and (2) of Panel A. "), italic font("Times New Roman", 9)

putdocx paragraph, halign(left) spacing(before, 2pt)
putdocx text ("* p < 0.10, ** p < 0.05, *** p < 0.01."), italic font("Times New Roman", 9)

* Save final export
putdocx save "Green_Bond_Stata_Table.docx", replace
display _newline "Success! Table compiled natively and saved as Green_Bond_Stata_Table.docx."

save "JoE_Stata_1.dta", replace
log close