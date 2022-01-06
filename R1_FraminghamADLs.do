global db3275 "db3275"
//global db3275 "danielbelsky"
//***************************************************************************//
//ADL Analysis 
//***************************************************************************//
use "/Users/$db3275/OneDrive - cumc.columbia.edu/Projects/mPoA/DunedinPoAm45/eLife/R1/R1_Framingham.dta", clear
	
foreach x in dnagi katz rosbres {
tab i`x' if zpace!=. & `x'8 == 0 
recode i`x' (2=1)
tab i`x' if zpace!=. & `x'8 == 0  
}	
	

//Clock Var macro
global clocks "zpace zpoam zar_horvath zar_hannum zar_phenoage zar_grimage" 
//Cell control macro
global cells "cd4t cd8t nk bcell mono gran cd8pcd28ncd45ran cd8naive cd4naive plasmablast"
	
foreach Y in rosbres katz dnagi { 			
	
//Univariate Fx 	
matrix M1 = J(1,5,999) 
matrix M2 = J(1,5,999) 
matrix M3 = J(1,5,999) 
foreach y in $clocks{ 
	//M1
	poisson i`Y' `y' `Y'8 sex age8 B1 B2 if `Y'8==0, robust nolog irr cluster(familyid )
		matrix A = exp(_b[`y']) , exp(_b[`y'] - invnormal(0.975)*_se[`y']), exp(_b[`y'] + invnormal(0.975)*_se[`y']), 2*normal(-abs(_b[`y']/_se[`y'])), e(N) 
		matrix rownames A = `y'
		matrix M1 = M1 \ A
	//M2
	poisson i`Y' `y' `Y'8 sex age8 B1 B2 $cells if `Y'8==0, robust nolog irr cluster(familyid )
		matrix A = exp(_b[`y']) , exp(_b[`y'] - invnormal(0.975)*_se[`y']), exp(_b[`y'] + invnormal(0.975)*_se[`y']), 2*normal(-abs(_b[`y']/_se[`y'])), e(N)  
		matrix rownames A = `y'
		matrix M2 = M2 \ A	
	//M3
	poisson i`Y' `y' `Y'8 sex age8 B1 B2 i.smk cpd if `Y'8==0, robust nolog irr cluster(familyid )
		matrix A = exp(_b[`y']) , exp(_b[`y'] - invnormal(0.975)*_se[`y']), exp(_b[`y'] + invnormal(0.975)*_se[`y']), 2*normal(-abs(_b[`y']/_se[`y'])), e(N) 
		matrix rownames A = `y'
		matrix M3 = M3 \ A	
	}
foreach x in 1 2 3{
	matrix `Y'M`x'=M`x'[2...,1...]
	matrix colnames `Y'M`x'=HR lb ub p N
	matrix list `Y'M`x'
	}
	
//Multivariate Fx 	
matrix X = J(1,5,999) 	
foreach x in zpoam zar_horvath zar_hannum zar_phenoage zar_grimage { 
		poisson i`Y' zpace `Y'8 `x' sex age8 B1 B2 if `Y'8==0, robust nolog irr cluster(familyid )
		matrix A = (exp(_b[zpace]) , exp(_b[zpace] - invnormal(0.975)*_se[zpace]), exp(_b[zpace] + invnormal(0.975)*_se[zpace]), 2*normal(-abs(_b[zpace]/_se[zpace])), e(N) )  \ (  exp(_b[`x']) , exp(_b[`x'] - invnormal(0.975)*_se[`x']), exp(_b[`x'] + invnormal(0.975)*_se[`x']), 2*normal(-abs(_b[`x']/_se[`x'])), e(N) )
		matrix rownames A = pace `x'
		matrix X = X \ A
		}
	matrix `Y'X=X[2...,1...]
	matrix colnames `Y'X=HR lb ub p N
	matrix list `Y'X	

putexcel set "/Users/$db3275/OneDrive - cumc.columbia.edu/Projects/mPoA/DunedinPoAm45/eLife/R1/R1_Framingham.xlsx", sheet(`Y') modify
putexcel B2 = matrix(`Y'M1), names 
putexcel B11 = matrix(`Y'M2), names 
putexcel B21 = matrix(`Y'M3), names 
putexcel B31 = matrix(`Y'X), names	

} 
//***************************************************************************//

	


