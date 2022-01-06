//global db3275 "danielbelsky"
global db3275 "db3275"


//******************************************************************************//
//LEHNE DATASET 
//******************************************************************************//
import delimited using "/Users/$db3275/OneDrive - cumc.columbia.edu/Projects/MK/DunedinPACE eLife/SugdenAlternateScores/LehneData_POAScores_11012021.csv", clear varn(1) delim(comma)
split sample_num, parse("sample ") gen(id)
destring id2, force gen(id)
split group, parse("Technical replicate group ") gen(batch)
destring batch2, force gen(batch)
drop id1 id2 batch1 batch2 
order id batch 

keep id batch dunedinage38_poamreliable dunedinage45_allepic450probes 

merge 1:1 id batch using "/Users/$db3275/OneDrive - cumc.columbia.edu/Projects/mPoA/DunedinPoAm45/Lehne/Lehne_Reliability.dta", nogen

rename dunedinage38_poamreliable poam_reliable 
rename dunedinage45_allepic450probes pace_all 
rename dunedinpoam_38 poam  
rename dunedinpoam_45 pace 

keep dnamage dnamagehannum dnamphenoage dnamgrimage poam* pace* agech1 sex id batch 

global clocks "pace pace_all poam poam_reliable dnamage dnamagehannum dnamphenoage dnamgrimage"

foreach x in $clocks{ 
	capture drop ar_`x'
	gen ar_`x'=.
	foreach b in 1 2 { 
		reg `x' agech1 if batch == `b'
		capture drop Z 
		predict Z if e(sample), r 
		replace ar_`x' = Z if batch ==`b'
		drop Z 
		}
	}
	//ICCs 
capture drop N 
egen N = count(dnamage), by(id)
	//Raw versions of measures 
matrix Fx = (999,999,999,999)	
foreach x of varlist $clocks { 
	quietly mixed ar_`x' if N==2 || id : batch 
	di `"***	`x'		***"' 
	estat icc 
	matrix A = r(icc2), r(ci2), e(N)	
	matrix rownames A = `x'
	matrix Fx = Fx \ A 
	}
matrix Fx_adj = Fx[2...,1...]
matrix colnames Fx_adj = ICC lb ub n 
matrix list Fx_adj, format(%9.2f)
	//Correlations 
preserve 
	matrix X =999
	keep id batch ar_* 
	reshape wide ar_* , i(id) j(batch)
	foreach x in $clocks {
		corr ar_`x'1 ar_`x'2
		matrix X = X \ r(rho)
		}
restore
matrix X = X[2...,1]
matrix rownames X = $clocks 
matrix colnames X = r
matrix list X 
matrix LADJ = X, Fx_adj
matrix list LADJ

//******************************************************************************//
putexcel set "/Users/$db3275/OneDrive - cumc.columbia.edu/Projects/mPoA/DunedinPoAm45/eLife/R1/R1_Reliability.xlsx", sheet(LEHNE) modify	
putexcel B2 = matrix(LADJ) ,names  
//******************************************************************************//


//******************************************************************************//
//******************************************************************************//
//SUGDEN  
//******************************************************************************//
//******************************************************************************//

	// "Sugden EPIC-EPIC n=28"  (DUNEDIN)
clear 
input ICC lb ub 
0.97 0.94 0.98
.89 .8 .94
0.93 0.87 0.96
.97 .94 .98
0.64 0.31 0.81
0.87 0.74 0.93
0.81 0.64 0.9
0.96 0.92 0.98
end 
mkmat ICC lb ub , mat(S1ADJ)
matrix A =(0.94,0.81,0.87,0.94,0.46,0.79,0.67, 0.92)'
matrix colnames A = r 
matrix B=J(8,1,28)
matrix colnames B = n 
matrix S1ADJ = A ,S1ADJ, B 
//******************************************************************************//
	// "Sugden 450k-EPIC n=350" (ERISK)
clear 
input ICC lb ub 
0.87 0.82 0.9
.6 .2 .79
0.80 0.76 0.83
.86 .83 .88
0.38 -0.03 0.6
0.71 0.66 0.76
0.75 0.52 0.84
0.78 0.48 0.88
end 
mkmat ICC lb ub , mat(S2ADJ)
matrix A =(0.79,0.65, 0.67, 0.76, 0.4, 0.56,0.68, 0.74)'
matrix colnames A = r 
matrix B=J(8,1,350)
matrix colnames B = n 
matrix S2ADJ = A ,S2ADJ, B 

foreach x in S1ADJ S2ADJ { 
	matrix rownames `x' = pace pace_all poam poam_reliable dnamage dnamagehannum dnamphenoage dnamgrimage 
	matrix list `x'
	}
	
//******************************************************************************//
//******************************************************************************//
