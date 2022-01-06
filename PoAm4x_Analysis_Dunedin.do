//global db3275 "db3275"
global db3275 "danielbelsky"


//*********************************************************************************//
//Data to R for Matrix Plot
use "/Users/$db3275/Box/Belsky/mPoA/mPoA45/Dunedin/PoAm4x_Dunedin.dta", clear 
destring wbc45-ntrphlsp45, replace force
	export delimited snum poam4x poam3x_38 paceofaging paceofaging_38 using "/Users/$db3275/OneDrive - cumc.columbia.edu/Projects/MK/DWB/PoAforMatrixPlot_210529.csv", delim(,) replace  
//*********************************************************************************//

//*********************************************************************************//
//Summary Statistics 
use "/Users/$db3275/Box/Belsky/mPoA/mPoA45/Dunedin/PoAm4x_Dunedin.dta", clear 
destring wbc45-ntrphlsp45, replace force
putexcel set "/Users/$db3275/OneDrive - cumc.columbia.edu/Projects/mPoA/DunedinPoAm45/Code/PoAm4xTables.xlsx", sheet(dunedin) modify

	//PoA Cross-wave Correlation
corr paceofaging paceofaging_38 
matrix A = r(rho), r(N)
matrix colnames A = r N
putexcel B2 = "r(PoA38,PoA45)"
putexcel C1 = matrix(A), colnames
	//PoA Cross-wave Correlation for those with DNAm at both time points
corr paceofaging paceofaging_38 if poam3x_38!=. & poam4x!=. 
matrix A = r(rho), r(N)
putexcel B3 = "DNAm sample r(PoA38,PoA45)"
putexcel C3 = matrix(A)
	//mPoA Cross-wave Correlation
corr poam3x_38 poam4x  
matrix A = r(rho), r(N)
putexcel B4 = "r(PoAm3x38,PoAm4x45)"
putexcel C4 = matrix(A)
	//mPoA within-wave Correlation
corr poam3x_45 poam4x 
matrix A = r(rho), r(N)
putexcel B5 = "r(PoAm3x45,PoAm4x45)"
putexcel C5 = matrix(A)

	//Summary Stats for PoA variables
sum paceofaging
matrix A = (r(mean), r(sd), r(N))
sum paceofaging if poam4x!=.
matrix A = A \ (r(mean), r(sd), r(N))
sum poam4x 
matrix A = A \ (r(mean), r(sd), r(N))
sum poam4x if sex==1
matrix A = A \ (r(mean), r(sd), r(N))
sum poam4x if sex==2
matrix A = A \ (r(mean), r(sd), r(N))
matrix rownames A = PoA_all PoA_DNAm PoAm4x PoAm4x_f PoAm4x_m
matrix colnames A = M SD N 
matrix list A
putexcel B7 = matrix(A), names

	//Corr of PoAm w/ PoA 
corr paceofaging poam4x
matrix A = r(rho), r(N)
corr paceofaging poam3x_45
matrix A = A \ (r(rho), r(N))
corr paceofaging_38 poam3x_38
matrix A = A \ (r(rho), r(N))
corr paceofaging_38 poam4x
matrix A = A \ (r(rho), r(N))
matrix rownames A = PoA_4x PoA_3x45 PoA38_3x38 PoA38_4x45
matrix colnames A = r N 
putexcel B15 = matrix(A), names
	
	//Sample Sizes for Supp Table with outcome measures 
count if poam4x!=. &  (balClsMax45!=. | Velocity_m45!=. | StepPlace45!=. | ChairStands45!=. | GripMax45!=. | PhyLimts45!=. )
matrix A = r(N)
count if poam4x!=. &  (pri45!=. | wmi45!=. | psi45!=. )
matrix A = A \ r(N)
count if poam4x!=. &  (Health45!=. | ZFacialAge45!=. )
matrix A = A \ r(N)
matrix rownames A = Phys Cog Subjective
matrix colnames A = N
putexcel B22=matrix(A), names

	//Intercorr of aging measures 
corr br_poam4x br_poam3x_45 br_horvath br_hannum br_phenoage br_grimage 
putexcel B30=matrix(r(C)), names
//*********************************************************************************//


//*********************************************************************************//
//Analysis of Validation Metrics
use "/Users/$db3275/Box/Belsky/mPoA/mPoA45/Dunedin/PoAm4x_Dunedin.dta", clear 
destring wbc45-ntrphlsp45, replace force
putexcel set "/Users/$db3275/OneDrive - cumc.columbia.edu/Projects/mPoA/DunedinPoAm45/Code/PoAm4xTables.xlsx", sheet(dunedin) modify

	//Associations with Age 45 Validation Metrics 
rename ZFacialAge45 zFacialAge45
#delimit ;
global PHENO "
	balClsMax45 
	Velocity_m45 
	StepPlace45 
	ChairStands45 
	GripMax45
	gpdom45
	PhyLimts45 
	pri45
	wmi45
	psi45
	Health45  
	zFacialAge45 
	" ; #delimit cr		
	//Regression Analysis
capture drop br_paceofaging	
gen br_paceofaging = paceofaging 
foreach y in br_poam4x br_poam3x_45 br_horvath br_hannum br_phenoage br_grimage br_paceofaging {
	capture drop z`y'
	egen z`y' =std(`y') if br_poam4x!=.
	matrix RFx_`y' = J(1,5,999)
	matrix RFx_`y'_cells = J(1,5,999)
	matrix RFx_`y'_smk = J(1,5,999)
	foreach x in $PHENO{
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
		matrix RFx_`y' = RFx_`y' \ A
		//Adjusted for Cells
		quietly reg Z z`y' sex wbc45-ntrphlsp45, robust 		
		matrix A = _b[z`y'] , _b[z`y'] - invttail(e(df_r),0.025)*_se[z`y'], _b[z`y'] + invttail(e(df_r),0.025)*_se[z`y'], 2*ttail(e(df_r),abs(_b[z`y']/_se[z`y'])), e(N)	
		matrix rownames A = `x'
		matrix RFx_`y'_cells = RFx_`y'_cells \ A
		//Adjusted for Smoking (pack years)
		quietly reg Z z`y' sex PackYrLifTm45, robust 
		matrix A = _b[z`y'] , _b[z`y'] - invttail(e(df_r),0.025)*_se[z`y'], _b[z`y'] + invttail(e(df_r),0.025)*_se[z`y'], 2*ttail(e(df_r),abs(_b[z`y']/_se[z`y'])), e(N)
		matrix rownames A = `x'
		matrix RFx_`y'_smk = RFx_`y'_smk \ A
		drop Z
		}
		drop z`y'
	foreach q in "" _cells _smk{ 
		matrix RFx_`y'`q' = RFx_`y'`q'[2...,1...]
		matrix colnames RFx_`y'`q' = r lb ub p N
		matrix list RFx_`y'`q' 
		}
	}	



	
