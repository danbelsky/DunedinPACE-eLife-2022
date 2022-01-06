
//global db3275 "db3275"
global db3275 "danielbelsky"

putexcel set "/Users/$db3275/OneDrive - cumc.columbia.edu/Projects/mPoA/DunedinPoAm45/Code/PoAm4xTables.xlsx", sheet(lehne) modify

use "/Users/$db3275/Box/Belsky/mPoA/mPoA45/UndSoc/PoAm4x_LehneReliability.dta", clear 

//DuneidnPoAm Reliability 
use "/Users/$db3275/OneDrive - cumc.columbia.edu/Projects/mPoA/DunedinPoAm45/Lehne/Lehne_Reliability.dta", clear
	//DunedinPoAm Test/Re-test reliability plot
preserve 
keep id batch dunedinpoam* 
reshape wide dunedinpoam_38 dunedinpoam_45, i(id) j(batch)
#delimit ;
twoway scatter dunedinpoam_452 dunedinpoam_451, msymbol(O) mcolor(gold) mlcolor(orange_red)
	|| lfit dunedinpoam_452 dunedinpoam_451 , lpattern(dash) lwidth(medthick) lcolor(orange_red)
	
	scheme(plotplain)
	ylabel(,labsize(medlarge))
	xlabel(,labsize(medlarge))
	ytitle(Replicate 2, size(medlarge))
	xtitle(Replicate 1, size(medlarge))
	legend(off)
	name(trtscatter_4x, replace)
	xsize(4) ysize(4)
	title(DunedinPACE, size(medium))
	; #delimit cr 
restore
graph export "/Users/$db3275/OneDrive - cumc.columbia.edu/Projects/mPoA/DunedinPoAm45/Code/Figures/Lehne_PoAm4x_Reliability.pdf" ,replace

	
//ICC
mixed dunedinpoam_45 || id : batch 
	estat icc  
matrix A = r(icc2), r(ci2), e(N), e(N_g)
mixed dunedinpoam_38 || id : batch 
	estat icc  
matrix A = A \ (r(icc2), r(ci2)), e(N), e(N_g)
matrix rownames A = PoAm4x PoAm3x
matrix colnames A = ICC lb ub N Ng
matrix list A 
putexcel B10 = matrix(A), names 


//Cross-measure comparison of reliability 
use "/Users/$db3275/OneDrive - cumc.columbia.edu/Projects/mPoA/DunedinPoAm45/Lehne/Lehne_Reliability.dta", clear

preserve 
keep id batch ar_dnamage ar_dnamagehannum ar_dnamphenoage ar_dnamgrimage ar_dunedinpoam_38 ar_dunedinpoam_45
reshape wide ar_dnamage ar_dnamagehannum ar_dnamphenoage ar_dnamgrimage ar_dunedinpoam_38 ar_dunedinpoam_45, i(id) j(batch)
foreach y in 1 2 { 
	foreach x in dnamage dnamagehannum dnamphenoage dnamgrimage dunedinpoam_38 dunedinpoam_45 { 
		rename ar_`x'`y' ar`y'`x'
			}
	}
reshape long ar1 ar2, i(id) j(clock) str
gen X = 3 if clock == "dnamage"
replace X = 4 if clock== "dnamagehannum" 
replace X = 5 if clock== "dnamphenoage"
replace X = 6 if clock== "dnamgrimage"
replace X = 2 if clock== "dunedinpoam_38" 
replace X = 1 if clock== "dunedinpoam_45"
capture label drop clock
label define clock 1 "Horvath" 2 "Hannum" 3 "PhenoAge" 4 "GrimAge" 5 "DunedinPoAm" 6 "DunedinPACE"
label values X clock
bys X: corr ar1 ar2 
#delimit ;
twoway scatter ar2 ar1, by(X, legend(off) note("") xrescale yrescale)
	msymbol(O) mcolor(dknavy)
	scheme(plotplain)
	ylabel(,labsize(medlarge))
	xlabel(,labsize(medlarge))
	ytitle(Replicate 2, size(medlarge))
	xtitle(Replicate 1, size(medlarge))
	|| lfit ar2 ar1 , lpattern(dash) lwidth(medthick) lcolor(red)
	name(trtscatter, replace)
	; #delimit cr 
restore 
graph export "/Users/$db3275/OneDrive - cumc.columbia.edu/Projects/mPoA/DunedinPoAm45/Code/Figures/Lehne_Comp_Reliability.pdf" ,replace



//Calculate ICCs 
matrix Fx = (999,999,999,999)	
foreach x of varlist dnamage dnamagehannum dnamphenoage dnamgrimage dunedinpoam_38 dunedinpoam_45 { 
	quietly mixed `x'  || id : batch 
	di `"***	`x'		***"' 
	estat icc 
	matrix A = r(icc2), r(ci2), e(N)	
	matrix rownames A = `x'
	matrix Fx = Fx \ A 
	}
matrix Fx_raw = Fx[2...,1...]
matrix colnames Fx_raw = ICC lb ub n 
matrix list Fx_raw, format(%9.2f)
putexcel B15 = "ICCs for Unadjusted Measures"
putexcel B16 = matrix (Fx_raw), names

matrix Fx = (999,999,999,999)	
foreach x of varlist dnamage dnamagehannum dnamphenoage dnamgrimage dunedinpoam_38 dunedinpoam_45 { 
	quietly mixed `x' agech1 sex   || id : batch 
	di `"***	`x'		***"' 
	estat icc 
	matrix A = r(icc2), r(ci2), e(N)	
	matrix rownames A = `x'
	matrix Fx = Fx \ A 
	}
matrix Fx_adj = Fx[2...,1...]
matrix colnames Fx_adj = ICC lb ub n 
matrix list Fx_adj, format(%9.2f)
putexcel B25 = "ICCs Adjusted for Age and Sex"
putexcel B26 = matrix (Fx_adj), names

//Calculate ICCs 
matrix Fx = (999,999,999,999)	
foreach x of varlist ar_* { 
	quietly mixed `x'  || id : batch 
	di `"***	`x'		***"' 
	estat icc 
	matrix A = r(icc2), r(ci2), e(N)	
	matrix rownames A = `x'
	matrix Fx = Fx \ A 
	}
matrix Fx_resid = Fx[2...,1...]
matrix colnames Fx_resid = ICC lb ub n 
matrix list Fx_resid, format(%9.2f)
putexcel B35 = "ICCs for Age-Residuals"
putexcel B36 = matrix (Fx_adj), names














	
