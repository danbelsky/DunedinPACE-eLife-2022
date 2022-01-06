
//global db3275 "db3275"
global db3275 "danielbelsky"

//Fx
use "/Users/$db3275/Box/Belsky/mPoA/mPoA45/ERisk/PoAm4x_ERisk.dta", clear 

//******************************************************************************//
//Regression Analysis
//******************************************************************************//
preserve
foreach x in horvath hannum phenoage grimage{ 
	replace br_`x' = br_`x'-age 
	}	
foreach x in poam4x poam3x horvath hannum phenoage grimage {
	egen zbr_`x' =  std(br_`x')
foreach z in "" _cells _smk _nsmk { 
	matrix Fx`z'=J(1,5,999)
	matrix colnames Fx`z' = b lb ub p N
	if `"`z'"'=="" {
		global C "sampsex"
		}
	if `"`z'"'=="_cells" {
		global C "sampsex plasmablast cd8pcd28ncd45ran cd8naive cd4naive cd8t cd4t nk bcell mono gran"
		}
	if `"`z'"'=="_smk" {
		global C "sampsex smkcure18 smkcnume18" //smkpkyre18
		}
	if `"`z'"'=="_nsmk" {
		global C "sampsex if smkcure==0"
		}		
	foreach y in zses zvic{
		quietly reg zbr_`x' `y' sampsex $C, cluster(familyid) robust 
		matrix A = _b[`y'] , _b[`y'] - invttail(e(df_r),0.025)*_se[`y'], _b[`y'] + invttail(e(df_r),0.025)*_se[`y'], 2*ttail(e(df_r),abs(_b[`y']/_se[`y'])), e(N)
		matrix rownames A = `y'
		matrix Fx`z' = Fx`z' \ A
		}
	quietly reg zbr_`x' zses zvic sampsex $C, cluster(familyid) robust 
		matrix B = ( _b[zses] , _b[zses] - invttail(e(df_r),0.025)*_se[zses], _b[zses] + invttail(e(df_r),0.025)*_se[zses], 2*ttail(e(df_r),abs(_b[zses]/_se[zses])), e(N) \ _b[zvic] , _b[zvic] - invttail(e(df_r),0.025)*_se[zvic], _b[zvic] + invttail(e(df_r),0.025)*_se[zvic], 2*ttail(e(df_r),abs(_b[zvic]/_se[zvic])), e(N) )
		matrix rownames B = Zses_mv Zvic_mv
		matrix Fx`z' = Fx`z' \ B	
	foreach y in "seswq35"  { 
		quietly reg zbr_`x' ib3.`y' $C, cluster(familyid) robust 
		#delimit ;
		matrix A = (_b[2.`y'] , _b[2.`y'] - invttail(e(df_r),0.025)*_se[2.`y'], _b[2.`y'] + invttail(e(df_r),0.025)*_se[2.`y'], 2*ttail(e(df_r),abs(_b[2.`y']/_se[2.`y'])), e(N) \ 
		_b[1.`y'] , _b[1.`y'] - invttail(e(df_r),0.025)*_se[1.`y'], _b[1.`y'] + invttail(e(df_r),0.025)*_se[1.`y'], 2*ttail(e(df_r),abs(_b[1.`y']/_se[1.`y'])), e(N) ); #delimit cr
		matrix rownames A = `y'_2 `y'_3
		matrix Fx`z' = Fx`z' \ A			
		}
	foreach y in "polyve512c" {
		quietly reg zbr_`x' i.`y' $C, cluster(familyid) robust 
		#delimit ;
		matrix A = (_b[1.`y'] , _b[1.`y'] - invttail(e(df_r),0.025)*_se[1.`y'], _b[1.`y'] + invttail(e(df_r),0.025)*_se[1.`y'], 2*ttail(e(df_r),abs(_b[1.`y']/_se[1.`y'])), e(N) \
		_b[2.`y'] , _b[2.`y'] - invttail(e(df_r),0.025)*_se[2.`y'], _b[2.`y'] + invttail(e(df_r),0.025)*_se[2.`y'], 2*ttail(e(df_r),abs(_b[2.`y']/_se[2.`y'])), e(N) \ 
		_b[3.`y'] , _b[3.`y'] - invttail(e(df_r),0.025)*_se[3.`y'], _b[3.`y'] + invttail(e(df_r),0.025)*_se[3.`y'], 2*ttail(e(df_r),abs(_b[3.`y']/_se[3.`y'])), e(N) ); #delimit cr
		matrix rownames A = `y'_1 `y'_2 `y'_3
		matrix Fx`z' = Fx`z' \ A	
		}
matrix `x'`z' = Fx`z'[2...,1...]
matrix list `x'`z'
		}
		}
//******************************************************************************//		
//Matrix of Correlation Effects (Exposure as continuous)
//******************************************************************************//
foreach x in poam4x poam3x horvath hannum phenoage grimage {		
	matrix SES_`x' = (`x'[1,1..5] \ `x'_cells[1,1..5] \ `x'_smk[1,1..5] \ `x'_nsmk[1,1..5]), (1\2\3\4) 
		matrix rownames SES_`x'= Base Cells Smk NonSmk
	matrix VIC_`x' = (`x'[2,1..5] \ `x'_cells[2,1..5] \ `x'_smk[2,1..5] \ `x'_nsmk[2,1..5]), (1\2\3\4)
		matrix rownames VIC_`x'= Base Cells Smk NonSmk
	matrix R`x'=J(1,6,999)
	matrix rownames R`x'=`x'
	} 
matrix SES= Rpoam4x \ SES_poam4x 
foreach x in 	poam3x horvath hannum phenoage grimage {
	matrix SES = SES \ R`x' \ SES_`x'
	}
matrix colnames SES = b lb ub p N 
matrix list SES

matrix VIC= Rpoam4x \ VIC_poam4x 
foreach x in 	poam3x horvath hannum phenoage grimage {
	matrix VIC = VIC \ R`x' \ VIC_`x'
	}
matrix colnames VIC = b lb ub p N 
matrix list VIC

putexcel set "/Users/$db3275/OneDrive - cumc.columbia.edu/Projects/mPoA/DunedinPoAm45/Code/PoAm4xTables_Checked.xlsx", sheet(Erisk) modify
putexcel B15 = matrix(SES), names
putexcel B50 = matrix(VIC), names

//******************************************************************************//
//Matrix of Cohen's d Effects (exposure as factor)
//******************************************************************************//
foreach x in poam4x poam3x horvath hannum phenoage grimage {		
	matrix SES_`x' = (`x'[5..6,1..5] \ `x'_cells[5..6,1..5] \ `x'_smk[5..6,1..5] \ `x'_nsmk[5..6,1..5]), (1\1\2\2\3\3\4\4) 
		matrix rownames SES_`x'= Base_Mid Base_High Cells_Mid Cells_High Smk_Mid Smk_High NonSmk_Mid NonSmk_High
	matrix VIC_`x' = (`x'[7..9,1..5] \ `x'_cells[7..9,1..5] \ `x'_smk[7..9,1..5] \ `x'_nsmk[7..9,1..5]), (1\1\1\2\2\2\3\3\3\4\4\4)
		matrix rownames VIC_`x'= Base_1 Base_2 Base_3 Cells_1 Cells_2 Cells_3 Smk_1 Smk_2 Smk_3 NonSmk_1 NonSmk_2 NonSmk_3
	matrix R`x'=J(1,6,999)
	matrix rownames R`x'=`x'
	} 
matrix SES= Rpoam4x \ SES_poam4x 
foreach x in 	poam3x horvath hannum phenoage grimage {
	matrix SES = SES \ R`x' \ SES_`x'
	}
matrix colnames SES = b lb ub p N 
matrix list SES
matrix VIC= Rpoam4x \ VIC_poam4x 
foreach x in 	poam3x horvath hannum phenoage grimage {
	matrix VIC = VIC \ R`x' \ VIC_`x'
	}
matrix colnames VIC = b lb ub p N 
matrix list VIC

putexcel set "/Users/$db3275/OneDrive - cumc.columbia.edu/Projects/mPoA/DunedinPoAm45/Code/PoAm4xTables_Checked.xlsx", sheet(Erisk_Cat) modify
putexcel B15 = matrix(SES), names
putexcel B75 = matrix(VIC), names
//******************************************************************************//
//******************************************************************************//

//******************************************************************************//
//MAIN TEXT FIGURE 5. E-Risk Cohen's D Effect Sizes vs. reference 
//******************************************************************************//
	//Regression Results 
foreach x in poam4x poam3x horvath hannum phenoage grimage {		
	matrix SES_`x' = (`x'[5..6,1..5] \ `x'_cells[5..6,1..5] \ `x'_smk[5..6,1..5] \ `x'_nsmk[5..6,1..5]), (1\1\2\2\3\3\4\4) 
		matrix rownames SES_`x'= Base2 Base1 Cells2 Cells1 Smk2 Smk1 NonSmk2 NonSmk1
	matrix VIC_`x' = (`x'[7..9,1..5] \ `x'_cells[7..9,1..5] \ `x'_smk[7..9,1..5] \ `x'_nsmk[7..9,1..5]), (1\1\1\2\2\2\3\3\3\4\4\4)
		matrix rownames VIC_`x'= Base 2 3 Cells 2 3 Smk 2 3  NonSmk 2 3
	matrix R`x'=J(1,6,999)
	matrix rownames R`x'=`x'
	} 
matrix SES =( SES_poam4x , J(8,1,1) \ SES_poam3x, J(8,1,2) \ SES_grimage, J(8,1,3) ), J(24,1,1)
matrix VIC =( VIC_poam4x , J(12,1,1) \ VIC_poam3x, J(12,1,2) \ VIC_grimage, J(12,1,3) ), J(36,1,2)
matrix list SES 
matrix FX = SES \ VIC
matrix colnames FX = b lb ub p N M BA EXP
matrix list FX


preserve 
use "/Users/$db3275/Box/Belsky/mPoA/mPoA45/ERisk/PoAm4x_ERisk.dta", clear 
egen zbr_poam4x = std(br_poam4x)
reg zbr_poam4x i.sesw sampsex 
margins, over(seswq35)
matrix A = r(table)
matrix B = (A[1,1...]\A[5,1...]\A[6,1...]\(1,2,3)\(1,1,1)\(1,1,1))'
matrix colnames B = b ll ul level exposure model 
reg zbr_poam4x i.sesw sampsex smkpkyre18
margins, over(seswq35) at(smkpkyre18==0)
matrix A = r(table)
matrix B = B \ (A[1,1...]\A[5,1...]\A[6,1...]\(1,2,3)\(1,1,1)\(2,2,2))'
reg zbr_poam4x i.sesw sampsex if smkcnume18==0
margins, over(seswq35)
matrix A = r(table)
matrix B = B \ (A[1,1...]\A[5,1...]\A[6,1...]\(1,2,3)\(1,1,1)\(3,3,3))'

reg zbr_poam4x i.polyve512c sampsex smkpkyre18
margins, over(polyve512c)
matrix A = r(table)
matrix B = B \ (A[1,1...]\A[5,1...]\A[6,1...]\(1,2,3,4)\J(1,4,2)\J(1,4,1))'
matrix colnames B = b ll ul level exposure model 
margins, over(polyve512c) at(smkpkyre18==0)
matrix A = r(table)
matrix B = B \ (A[1,1...]\A[5,1...]\A[6,1...]\(1,2,3,4)\J(1,4,2)\J(1,4,2))'
reg zbr_poam4x i.polyve512c sampsex if smkcnume18==0
margins, over(polyve512c)
matrix A = r(table)
matrix B = B \ (A[1,1...]\A[5,1...]\A[6,1...]\(1,2,3,4)\J(1,4,2)\J(1,4,3))'
restore


