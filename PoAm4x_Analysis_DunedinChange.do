
//global db3275 "db3275"
global db3275 "danielbelsky"

//****************************************************************************************//
//****************************************************************************************//
// CHANGE SCORE ANALYSIS 38-45
//****************************************************************************************//
//****************************************************************************************//

use "/Users/$db3275/Box/Belsky/mPoA/mPoA45/Dunedin/PoAm4x_Dunedin.dta", clear 
destring wbc45-ntrphlsp45, replace force 
putexcel set "/Users/$db3275/OneDrive - cumc.columbia.edu/Projects/mPoA/DunedinPoAm45/Code/PoAm4xTables.xlsx", sheet(dunedin_change) modify

count if poam4x!=. & (balance38!=. & balance45!=. | grip38!=. & grip45!=. | limitations38!=. & limitations 45!=. | srh38!=. & srh45!=. | face38!=. & face45!=.)
putexcel A3 = "DunedinPoAm4x N="
putexcel B3 = matrix(r(N))
count if poam4x!=.	& (balance38!=. & balance45!=. | grip38!=. & grip45!=. | limitations38!=. & limitations 45!=. | srh38!=. & srh45!=. | face38!=. & face45!=.)
putexcel A4 = "Pace of Aging N="
putexcel B4 = matrix(r(N))
count if paceofaging!=.	& poam4x!=. & (balance38!=. & balance45!=. | grip38!=. & grip45!=. | limitations38!=. & limitations 45!=. | srh38!=. & srh45!=. | face38!=. & face45!=.)
putexcel A5 = "DunedinPoAm4x & Pace of Aging N="
putexcel B5 = matrix(r(N))

//****************************************************************************************//
//Summary Statistics
//****************************************************************************************//
global DPHENO "DB DG d_gpdom DL Diq iPoor DF"
	//Change Score Summary Statistics
mean $DPHENO if poam4x !=.	
tabstat $DPHENO, s(mean sd n) by(sex) save
matrix A = r(StatTotal)' , r(Stat1)', r(Stat2)'
putexcel E8 = "Women"
putexcel H8 = "Men"
putexcel A9 = "Change Score Summary Statistics - All Available Data"
putexcel A10 = matrix(A), names
	//Summary Statistics for SMs w DunedinPoAm4x
tabstat $DPHENO if poam4x!=., s(mean sd n) by(sex) save
matrix A = r(StatTotal)' , r(Stat1)', r(Stat2)'
putexcel A19 = "Change Score Summary Statistics - SMs w/ DunedinPoAm4x Data"
putexcel A20 = matrix(A), names
	//Summary Statistics for SMs w Pace of Aging
tabstat $DPHENO if paceofaging!=., s(mean sd n) by(sex) save
matrix A = r(StatTotal)' , r(Stat1)', r(Stat2)'
putexcel A29 = "Change Score Summary Statistics - SMs w/ Pace of Aging Data"
putexcel A30 = matrix(A), names

	//Means & 95% CIs of change in terms of baseline SDs for in-text reporting
matrix Fx = J(1,4,999)
foreach y in zDB zDG zDL Diq DF {
	mean `y' if poam4x!=.
	matrix A = _b[`y'] , _b[`y'] - invttail(e(df_r),0.025)*_se[`y'], _b[`y'] + invttail(e(df_r),0.025)*_se[`y'], e(N)
	matrix rownames A = `y'
	matrix Fx = Fx \ A
	}
matrix Fx =Fx[2...,1...]
matrix colnames Fx = M lb ub N
matrix list Fx
putexcel A68 = "Change in Terms of Baseline SD"
putexcel A69 = matrix(Fx), names

	//% with Fair or Poor SRH at 38 and 45 for in-text reporting
tab srh45 if poam4x!=. & DSRH!=.
tab srh38 if poam4x!=. & DSRH!=.


use "/Users/$db3275/Box/Belsky/mPoA/mPoA45/Dunedin/PoAm4x_Dunedin.dta", clear 
destring wbc45-ntrphlsp45, replace force 
putexcel set "/Users/$db3275/OneDrive - cumc.columbia.edu/Projects/mPoA/DunedinPoAm45/Code/PoAm4xTables.xlsx", sheet(dunedin_change) modify
	//Regression Analysis
foreach y in br_poam4x br_poam3x_45 br_horvath br_hannum br_phenoage br_grimage paceofaging{
	capture drop z`y'
	egen z`y' =std(`y')
	matrix DFx_`y' = J(1,5,999)
	matrix DFx_`y'_cells = J(1,5,999)
	matrix DFx_`y'_smk = J(1,5,999)
	foreach x in $DPHENO{
		if `"`x'"'=="iPoor"{
				//Base Model
			quietly poisson `x' z`y' sex, robust 
			matrix A = exp(_b[z`y']) ,  exp(_b[z`y'] - invnormal(0.975)*_se[z`y']), exp(_b[z`y'] + invnormal(0.975)*_se[z`y']), 2*normal(-abs(_b[z`y']/_se[z`y'])), e(N)
			matrix rownames A = `x'
			matrix DFx_`y' = DFx_`y' \ A			
				//Adjusted for Cells
			quietly poisson `x' z`y' sex wbc45-ntrphlsp45, robust 
			matrix A = exp(_b[z`y']) ,  exp(_b[z`y'] - invnormal(0.975)*_se[z`y']), exp(_b[z`y'] + invnormal(0.975)*_se[z`y']), 2*normal(-abs(_b[z`y']/_se[z`y'])), e(N)
			matrix rownames A = `x'
			matrix DFx_`y'_cells = DFx_`y'_cells \ A
				//Adjusted for Smoking
			quietly poisson `x' z`y' sex PackYrLifTm45, robust 
			matrix A = exp(_b[z`y']) ,  exp(_b[z`y'] - invnormal(0.975)*_se[z`y']), exp(_b[z`y'] + invnormal(0.975)*_se[z`y']), 2*normal(-abs(_b[z`y']/_se[z`y'])), e(N)
			matrix rownames A = `x'
			matrix DFx_`y'_smk = DFx_`y'_smk \ A
			}
		else{
		capture drop Z 
		gen Z =.
		foreach s in 1 2{
			quietly sum `x' if sex==`s'
			replace Z = (`x'-r(mean))/r(sd) if sex==`s'
			}
		//Base Model
		quietly reg Z z`y' sex, robust 
		matrix A = _b[z`y'] , _b[z`y'] - invttail(e(df_r),0.025)*_se[z`y'], _b[z`y'] + invttail(e(df_r),0.025)*_se[z`y'], 2*ttail(e(df_r),abs(_b[z`y']/_se[z`y'])), e(N)
		matrix rownames A = `x'
		matrix DFx_`y' = DFx_`y' \ A
		//Adjusted for Cells
		quietly reg Z z`y' sex wbc45-ntrphlsp45, robust 		
		matrix A = _b[z`y'] , _b[z`y'] - invttail(e(df_r),0.025)*_se[z`y'], _b[z`y'] + invttail(e(df_r),0.025)*_se[z`y'], 2*ttail(e(df_r),abs(_b[z`y']/_se[z`y'])), e(N)	
		matrix rownames A = `x'
		matrix DFx_`y'_cells = DFx_`y'_cells \ A
		//Adjusted for Smoking (pack years)
		quietly reg Z z`y' sex PackYrLifTm45, robust 
		matrix A = _b[z`y'] , _b[z`y'] - invttail(e(df_r),0.025)*_se[z`y'], _b[z`y'] + invttail(e(df_r),0.025)*_se[z`y'], 2*ttail(e(df_r),abs(_b[z`y']/_se[z`y'])), e(N)
		matrix rownames A = `x'
		matrix DFx_`y'_smk = DFx_`y'_smk \ A
		drop Z
		}
		}
		drop z`y'
	foreach q in "" _cells _smk{ 
		matrix DFx_`y'`q' = DFx_`y'`q'[2...,1...]
		matrix colnames DFx_`y'`q' = r lb ub p N
		matrix list DFx_`y'`q' 
		}
	}

foreach y in br_poam4x br_poam3x_45 br_horvath br_hannum br_phenoage br_grimage paceofaging {
	foreach v in 1 2 3{
		matrix A`v' = J(1,5,`v')
		}
	matrix rownames A1 = `y'
	matrix rownames A2 = `y'_cells
	matrix rownames A3 = `y'_smk
	matrix `y' = A1 \ DFx_`y' \ A2 \ DFx_`y'_cells \ A3 \ DFx_`y'_smk
	}
matrix A = br_poam4x \ br_poam3x_45 \ br_horvath \ br_hannum \ br_phenoage \ br_grimage \ paceofaging
matrix colnames A = b lb ub p N 
putexcel B40 = matrix(A), names 
		
	
//Sensitivity Analysis - in-text reporting of residualized change analysis of IQ
preserve
	capture drop zbr_poam4x
	egen zbr_poam4x=std(br_poam4x)
	egen ziq45=std(fsiq45a)
	egen zDiq=std(Diq)
	reg zDiq zbr_poam4x sex, robust  
	reg zDiq zbr_poam4x wfsiq713 sex, robust 
restore



