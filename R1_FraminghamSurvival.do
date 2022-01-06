//global db3275 "db3275"
global db3275 "danielbelsky"
//***************************************************************************//
//Survival Analysis 
//***************************************************************************//
use "/Users/$db3275/OneDrive - cumc.columbia.edu/Projects/mPoA/DunedinPoAm45/eLife/R1/R1_Framingham.dta", clear


foreach Y in cvd stroke_tia death { 
stset T_`Y', failure(`Y'=1)
}


//Clock Var macro
global clocks "zpace zpoam zar_horvath zar_hannum zar_phenoage zar_grimage" 
//Cell control macro
global cells "cd4t cd8t nk bcell mono gran cd8pcd28ncd45ran cd8naive cd4naive plasmablast"

foreach Y in cvd chf chd afix dem death cvddeath chddeath stroke stroke_tia { 
	
stset T_`Y', failure(`Y'=1)

matrix M1 = J(1,5,999) 
matrix M2 = J(1,5,999) 
matrix M3 = J(1,5,999) 
foreach y in $clocks{ 
	//M1
	stcox `y' age8 sex B1 B2, cluster(familyid) robust 
		matrix A = exp(_b[`y']) , exp(_b[`y'] - invnormal(0.975)*_se[`y']), exp(_b[`y'] + invnormal(0.975)*_se[`y']), 2*normal(-abs(_b[`y']/_se[`y'])), e(N) 
		matrix rownames A = `y'
		matrix M1 = M1 \ A
	//M2
	stcox `y' age8 sex B1 B2 $cells, cluster(familyid) robust
		matrix A = exp(_b[`y']) , exp(_b[`y'] - invnormal(0.975)*_se[`y']), exp(_b[`y'] + invnormal(0.975)*_se[`y']), 2*normal(-abs(_b[`y']/_se[`y'])), e(N)  
		matrix rownames A = `y'
		matrix M2 = M2 \ A	
	//M3
	stcox `y' age8 sex B1 B2 smk cpd, cluster(familyid) robust
		matrix A = exp(_b[`y']) , exp(_b[`y'] - invnormal(0.975)*_se[`y']), exp(_b[`y'] + invnormal(0.975)*_se[`y']), 2*normal(-abs(_b[`y']/_se[`y'])), e(N) 
		matrix rownames A = `y'
		matrix M3 = M3 \ A	
	}
foreach x in 1 2 3{
	matrix `Y'M`x'=M`x'[2...,1...]
	matrix colnames `Y'M`x'=HR lb ub p N
	matrix list `Y'M`x'
	}

//Clock-adjusted PACE Effects 	
matrix X = J(1,5,999) 	
foreach x in zpoam zar_horvath zar_hannum zar_phenoage zar_grimage { 
	stcox zpace `x' age8 sex B1 B2, cluster(familyid)
		matrix A = (exp(_b[zpace]) , exp(_b[zpace] - invnormal(0.975)*_se[zpace]), exp(_b[zpace] + invnormal(0.975)*_se[zpace]), 2*normal(-abs(_b[zpace]/_se[zpace])), e(N) )  \ (  exp(_b[`x']) , exp(_b[`x'] - invnormal(0.975)*_se[`x']), exp(_b[`x'] + invnormal(0.975)*_se[`x']), 2*normal(-abs(_b[`x']/_se[`x'])), e(N) )
		matrix rownames A = pace `x'
		matrix X = X \ A
		}
	matrix `Y'X=X[2...,1...]
	matrix colnames `Y'X=HR lb ub p N
	matrix list `Y'X
		
putexcel set "/Users/$db3275/OneDrive - cumc.columbia.edu/Projects/mPoA/DunedinPoAm45/eLife/R1/R1_Framingham.xlsx", sheet(`Y') modify

count if `Y'==1 & zpace<. & T_`Y'>0
local N = r(N)
capture drop temp 
gen temp = age8+T_`Y'
sum temp if `Y'==1 & T_`Y'>0 & zpace<.
local M = round(r(mean), .01)
local SD = round(r(sd), .01)
sum T_`Y'
local X = round(r(max),.01)
drop temp 

putexcel A1 = `"`Y': N=`N' cases at mean age = `M' (SD=`SD') over up to `X' years of follow-up"'
putexcel B2 = matrix(`Y'M1), names 
putexcel B11 = matrix(`Y'M2), names 
putexcel B21 = matrix(`Y'M3), names 
putexcel B31 = matrix(`Y'X), names

}

