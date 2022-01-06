
//global db3275 "db3275"
global db3275 "danielbelsky"

use "/Users/$db3275/Box/Belsky/mPoA/mPoA45/UndSoc/PoAm4x_UndSoc.dta", clear

//Data to R for Matrix Plot
export delimited id br_poam4x br_poam3x ar_horvath ar_hannum ar_phenoage ar_grimage using "/Users/$db3275/OneDrive - cumc.columbia.edu/Projects/MK/DWB/PoAforMatrixPlot_US210607.csv", delim(,) replace  


putexcel set "/Users/$db3275/OneDrive - cumc.columbia.edu/Projects/mPoA/DunedinPoAm45/Code/PoAm4xTables.xlsx", sheet(undsoc) modify

//Correlation with Chronological age [r=0.26 for Women and r=0.41 for men; r=0.32 overall] 
sum poam4x 
matrix A = r(mean) , r(sd)
corr poam4x age 
matrix A = A, r(rho), r(N)
sum poam4x if sex==1
matrix B = r(mean) , r(sd)
corr poam4x age if sex==1
matrix B = B, r(rho), r(N)
matrix A = A \ B
sum poam4x if sex==2
matrix B = r(mean) , r(sd)
corr poam4x age if sex==2
matrix B = B, r(rho), r(N)
matrix A = A \ B
matrix colnames A = M SD r_Age N
matrix rownames A = All F M 
matrix list A 
putexcel B2 = matrix(A), names


//Regression tables 
egen zbr_poam4x = std(br_poam4x)
egen zbr_poam3x = std(br_poam3x)
foreach x in poam4x poam3x horvath hannum phenoage grimage{
	egen zar_`x' = std(ar_`x')
	}

foreach y in srh hd_log kdm_advance PAA mmgsdval {
capture drop z`y'
	gen z`y'=.
	foreach S in 1 2 { 	
		quietly sum `y' if sex ==`S'
		replace z`y' = (`y'-r(mean))/r(sd) if sex ==`S'
		}
	}

global base "sex age"	
global cells "sex age cd8t-gran plasmablast cd8pcd28ncd45ran cd8naive cd4naive"
global smk "sex age i.smk"	
foreach x in poam4x poam3x horvath hannum phenoage grimage { 
foreach M in base cells smk{ 
	matrix Fx = J(1,5,999)
	foreach y in srh hd_log kdm_advance PAA mmgsdval {
		quietly reg z`y' zar_`x' $`M' , robust 
		#delimit ;
		matrix A = _b[zar_`x'], 
			_b[zar_`x'] - invttail(e(df_r),0.025)*_se[zar_`x'], 
			_b[zar_`x'] + invttail(e(df_r),0.025)*_se[zar_`x'], 
			2*ttail(e(df_r),abs(_b[zar_`x']/_se[zar_`x'])), e(N) ; #delimit cr 
		matrix rownames A = `y'
		matrix Fx = Fx \ A
		}
	matrix `x'_`M' = Fx[2...,1...]
	matrix colnames `x'_`M' = b lb ub p N
	matrix list `x'_`M' 
	}
}
foreach n in 1 2 3 {
	matrix A`n' = J(1,5,`n')
	}
foreach x in poam4x poam3x horvath hannum phenoage grimage { 
	matrix rownames A1 = `x'
	matrix rownames A2 = `x'_cells
	matrix rownames A3 = `x'_smk
	matrix `x' = A1 \ `x'_base \ A2 \ `x'_cells \ A3 \ `x'_smk
	matrix colnames `x' = b lb ub p N 
	}
matrix A = poam4x \poam3x \horvath \hannum \phenoage \grimage	

putexcel B10 = matrix(A), names 

//Effect-size comparison
matrix Fx = J(5,1,1) \ J(5,1,2) \ J(5,1,3) \ J(5,1,4) \ J(5,1,5) \ J(5,1,6)
matrix X = (1\2\3\4\5)
matrix Fx = Fx , (X\X\X\X\X\X)
matrix colnames Fx = X Y 
matrix B = poam4x_base 
foreach x in  poam3x horvath hannum phenoage grimage { 
matrix B = B \ `x'_base
} 
matrix Fx = B, Fx 
matrix list Fx

//Figure for effect-size comparison
preserve 
clear 
svmat2 Fx, names(col) rnames(ystring)
capture label drop X 
label define X 1 "DunedinPACE" 2 "DunedinPoAm" 3 "Horvath Clock" 4 "Hannum Clock" 5 "PhenoAge Clock" 6 "GrimAge Clock"
label values X X
recode Y (4=1) (3=2) (2=3) (1=4)
drop if Y ==5
capture label drop Y 
label define Y 1 "Phenotypic Age Advancement" 2 "KDM BA Advancement"  3 "Homeostatic Dysregulation" 4 "Self-rated Health"
label values Y Y
list 
#delimit ; 
twoway rcap lb ub X if X==1, lcolor(orange_red) by(Y, legend(off) note("")) scheme(plotplain)
	|| rcap lb ub X if X==2, lcolor(orange)  
	|| rcap lb ub X if X==3, lcolor(maroon)  
	|| rcap lb ub X if X==4, lcolor(purple) 
	|| rcap lb ub X if X==5, lcolor(navy) 
	|| rcap lb ub X if X==6, lcolor(dknavy) 
	
	|| scatter b X if X==1, mcolor(gold) mlcolor(orange_red) msize(large) msymbol(O)  
	|| scatter b X if X==2, mcolor(orange) msize(large) msymbol(O)  
	|| scatter b X if X==3, mcolor(maroon) msize(large) msymbol(O)  
	|| scatter b X if X==4, mcolor(purple) msize(large) msymbol(O)  
	|| scatter b X if X==5, mcolor(navy) msize(large) msymbol(O)  
	|| scatter b X if X==6, mcolor(dknavy) msize(large) msymbol(O)  
	xscale(range(.75 6.25))
	xlabel(1(1)6, valuelabels angle(50) labsize(medlarge))
	xtitle("")
	ylabel(,labsize(medlarge))
	ytitle(Effect-size - Pearson r, size(medlarge))
	name(compfx,replace)
; #delimit cr 
restore
graph export "/Users/$db3275/OneDrive - cumc.columbia.edu/Projects/mPoA/DunedinPoAm45/Code/Figures/UndSoc_CompFx.pdf" ,replace



reg zbr_poam4x i.srh sex age, robust 

	//PoAm4x by Age 
#delimit ; 
twoway scatter br_poam4x age , mcolor(gold) msymbol(O)
	|| lfit br_poam4x age, lcolor(blue) lwidth(medthick) lpattern(dash)
	by(sex, note("") legend(off)) 
	xtitle(Chronological Age, size(medlarge))
	ytitle(DunedinPACE, size(medlarge))
	ylabel(.5(.25)2, labsize(medlarge) format(%9.2f))
	xlabel(20(20)100, labsize(medlarge) format(%9.0f))
	yline(1, lcolor(gs13))	
	scheme(plotplain)
name(DunedinPoAm4x, replace)
; #delimit cr
graph export "/Users/$db3275/OneDrive - cumc.columbia.edu/Projects/mPoA/DunedinPoAm45/Code/Figures/UndSoc_PoAm4x_Age_bysex.pdf" ,replace



	//PoAm4x by Age -- matched to elife
#delimit ; 
twoway scatter br_poam4x age if sex==1 , mcolor(gold%70) mlcolor(gold%1) msymbol(O)
	|| scatter br_poam4x age if sex==2 , mcolor(blue%50) msymbol(+)
	|| lfitci br_poam4x age if sex==1, lcolor(orange) lwidth(medthick) lpattern(dash)
		ciplot(rarea) acolor(orange%30) alpattern(solid) alcolor(orange%5)
		
	|| lfitci br_poam4x age if sex==2, lcolor(blue) lwidth(medthick) lpattern(dash)
		ciplot(rarea) acolor(blue%20) alpattern(solid) alcolor(blue%5)
	
	legend(cols(1) pos(11) ring(0) lab(1 "Women") lab(2 "Men") order (1 2))
	xtitle(Chronological Age, size(medlarge))
	ytitle(DunedinPACE, size(medlarge))
	ylabel(.5(.25)2, labsize(medlarge) format(%9.2f))
	xlabel(20(20)100, labsize(medlarge) format(%9.0f))
	yline(1, lcolor(gs13))	
	scheme(plotplain)
	xsize(4) ysize(4)
name(DunedinPoAm4x_elife, replace)
; #delimit cr
graph export "/Users/$db3275/OneDrive - cumc.columbia.edu/Projects/mPoA/DunedinPoAm45/Code/Figures/UndSoc_PoAm4x_Age.pdf" ,replace

//PAA
#delimit ;
	twoway scatter  zbr_poam4x zPAA, mcolor(purple%30) mlcolor(purple%1) msymbol(O) 
		|| lfit  zbr_poam4x zPAA, lcolor(pink) lwidth(medthick) lpattern(dash)
		scheme(plotplain)
		legend(off)
		xlabel(,labsize(medlarge))
		ylabel(,labsize(medlarge))
		ytitle(DunedinPACE (z-score), size(medlarge))
		xtitle("Phenotypic-Age Advancement (z-score)", size(medlarge)) 
		yline(0,lcolor(gs10) )
		name(paa, replace)
	; #delimit cr
//KDM
#delimit ;
	twoway scatter  zbr_poam4x zkdm_advance, mcolor(blue%30) mlcolor(blue%1) msymbol(O) 
		|| lfit  zbr_poam4x zkdm_advance, lcolor(midblue) lwidth(medthick) lpattern(dash)
		scheme(plotplain)
		legend(off)
		xlabel(,labsize(medlarge))
		ylabel(,labsize(medlarge))
		ytitle(DunedinPACE (z-score), size(medlarge))
		xtitle("KDM-BA Advancement (z-score)", size(medlarge)) 
		yline(0,lcolor(gs10) 	)
		name(kdmba,replace)
	; #delimit cr	
//HD
#delimit ;
	twoway scatter  zbr_poam4x zhd_log, mcolor(orange%30) mlcolor(orange%1) msymbol(O) 
		|| lfit zbr_poam4x zhd_log , lcolor(orange_red) lwidth(medthick) lpattern(dash)
		scheme(plotplain)
		legend(off)
		xlabel(,labsize(medlarge))
		ylabel(,labsize(medlarge))
		ytitle(DunedinPACE (z-score), size(medlarge))
		xtitle("Homeostatic Dysregulation (z-score)", size(medlarge)) 
		yline(0,lcolor(gs10) 	)
		name(hd,replace)
	; #delimit cr	
	//PoAm4x & SRH 
foreach x in poam4x  {
	quietly reg zbr_`x' i.srh c.age##c.age##sex , robust 
	quietly margins, over(srh)
	#delimit ;
	if `"`x'"'=="poam4x" { ; local T "DunedinPACE"; local C "gold"; local Cl "orange_red" ; }; 
	if `"`x'"'=="grimage" { ; local T "GrimAge"; local C "navy"; local Cl "navy"; }; 
	marginsplot, scheme(plotplain) 
		plot1opts(connect(none) msymbol(O) msize(vlarge) mcolor(`C') mlcolor(`Cl'))
		ci1opts(lcolor(`C') lwidth(medthick))
		title("")
		ytitle(`T' Z-score, size(medlarge))
		xtitle(Self-rated Health, size(medlarge))
		ylabel(-.5(.25)1, labsize(medlarge) format(%9.2f))
		xlabel( ,labsize(medlarge) angle(40))
		yline(0)
		xscale(range(.75 5.25))
		name(`x', replace)
	; #delimit cr 
	}
graph combine kdmba paa hd poam4x, scheme(s1mono)
graph export "/Users/$db3275/OneDrive - cumc.columbia.edu/Projects/mPoA/DunedinPoAm45/Code/Figures/UndSoc_PoAm4x_BA.pdf" ,replace






/*


//Sex-standardized aging measures and outcomes 
capture drop cage 
gen cage = age-60 
foreach y in  mmgsdval srh htfev ar_poam4x ar_grimage { 
	capture drop z`y'
	gen z`y'=.
	foreach S in 1 2 { 	
		quietly sum `y' if sex ==`S'
		replace z`y' = (`y'-r(mean))/r(sd) if sex ==`S'
		}
	}

	//Dichotomous fair/poor self-rated health
tab srh 
capture drop lsrh 
recode srh (1/3=0) (4/5=1), gen(lsrh)	




//Grip by Age 
#delimit ; 
twoway scatter mmgsdval age , mcolor(orange_red%40) mlcolor(orange_red%1) msymbol(O)
	|| qfit mmgsdval age, lcolor(blue) lwidth(medthick) lpattern(dash)
	by(sex, note("") legend(off)) 
	xtitle(Chronological Age, size(medlarge))
	ytitle(Grip Strength, size(medlarge))
	ylabel(, labsize(medlarge) format(%9.0f))
	xlabel(20(20)100, labsize(medlarge) format(%9.0f))
	scheme(plotplain)
name(Grip, replace)
; #delimit cr

//FEV1
#delimit ; 
twoway scatter htfev age , mcolor(midblue%40) mlcolor(midblue%1) msymbol(O)
	|| qfit htfev age, lcolor(blue) lwidth(medthick) lpattern(dash)
	by(sex, note("") legend(off)) 
	xtitle(Chronological Age, size(medlarge))
	ytitle(Lung Function (FEV1 - liters), size(medlarge))
	ylabel(, labsize(medlarge) format(%9.0f))
	xlabel(20(20)100, labsize(medlarge) format(%9.0f))
	scheme(plotplain)
name(FEV, replace)
; #delimit cr

//HbA1C
capture drop hba1cp 
gen hba1cp = (hba1c/10.929) + 2.15
#delimit ; 
twoway scatter hba1cp age , mcolor(cranberry%30) mlcolor(cranberry%1) msymbol(O)
	|| qfit hba1cp age if hba1c<50, lcolor(blue) lwidth(medthick) lpattern(dash)
	by(sex, note("") legend(off)) 
	xtitle(Chronological Age, size(medlarge))
	ytitle("HbA1C (%)", size(medlarge))
	ylabel(, labsize(medlarge) format(%9.0f))
	xlabel(20(20)100, labsize(medlarge) format(%9.0f))
	scheme(plotplain)
name(hba1c, replace)
; #delimit cr


//Self-rated Health by Age 
quietly logit lsrh c.age##sex 
margins, at(age=(20(5)100)) over(sex)
#delimit ; 
marginsplot, legend(ring(0) pos(11) cols(1) symxsize(5))
recastci(rarea) scheme(plotplain)
plot1opts(msymbol(none) lcolor(gold) lwidth(thick))
plot2opts(msymbol(none) lcolor(orange) lwidth(thick))
ci1opts(fcolor(gold%30) lcolor(none))
ci2opts(fcolor(orange%30) lcolor(none)) 
xlabel(20(20)100, labsize(medlarge))
ylabel(0(.2)1, labsize(medlarge))
title("")
xtitle(Chronological Age, size(medlarge))
ytitle(Probability of Fair or Poor Self-rated Health, size(medlarge))
xsize(4) ysize(4)
; #delimit cr 



//Sick
quietly logit sick c.age##sex 
margins, at(age=(20(5)100)) over(sex)
#delimit ; 
marginsplot, legend(ring(0) pos(11) cols(1) symxsize(5))
recastci(rarea) scheme(plotplain)
plot1opts(msymbol(none) lcolor(rose) lwidth(thick))
plot2opts(msymbol(none) lcolor(cranberry) lwidth(thick))
ci1opts(fcolor(gold%30) lcolor(none))
ci2opts(fcolor(orange%30) lcolor(none)) 
xlabel(20(20)100, labsize(medlarge))
ylabel(0(.2)1, labsize(medlarge))
title("")
xtitle(Chronological Age, size(medlarge))
ytitle(Prob. of Prevalent Chronic Disease/ Disability, size(medlarge))
xsize(4) ysize(4)
; #delimit cr 





	//PoAm4x & SRH 
capture drop zar_poam4x 
gen zar_poam4x = zbr_poam4x
foreach x in poam4x grimage {
	quietly reg zar_`x' i.srh c.cage##c.cage##sex , robust 
	quietly margins, over(srh)
	#delimit ;
	if `"`x'"'=="poam4x" { ; local T "DunedinPoAm4x"; local C "gold"; local Cl "orange_red" ; }; 
	if `"`x'"'=="grimage" { ; local T "GrimAge"; local C "navy"; local Cl "navy"; }; 
	marginsplot, scheme(plotplain) 
		plot1opts(connect(none) msymbol(O) msize(vlarge) mcolor(`C') mlcolor(`Cl'))
		ci1opts(lcolor(`C') lwidth(medthick))
		title(`T', size(medlarge))
		ytitle(`T' Z-score, size(medlarge))
		xtitle("")
		ylabel(-.5(.25)1, labsize(medlarge) format(%9.2f))
		xlabel( ,labsize(medlarge) angle(40))
		yline(0)
		xscale(range(.75 5.25))
		name(`x', replace)
	; #delimit cr 
	}
graph combine poam4x grimage, cols(2) scheme(s1mono) ycommon name(srh, replace)

	//PoAm4x & Grip 
foreach x in poam4x grimage {
	 reg  mmgsdval c.zar_`x' c.cage##c.cage##sex c.zar_`x'#c.cage, robust 
	 margins, at(cage=(-30(10)30) zar_`x'=(-1 1))
	#delimit ;
	if `"`x'"'=="poam4x" { ; local T "DunedinPoAm4x"; local C1 "gold"; local C2 "orange_red"; local l1 "Slow"; local l2 "Fast"; }; 
	if `"`x'"'=="grimage" { ; local T "GrimAge"; local C1 "midblue"; local C2 "navy"; local l1 "Young"; local l2 "Old"; }; 
	marginsplot, scheme(plotplain) 
		plot1opts( msymbol(none) lwidth(thick) lcolor(`C1'))
		plot2opts( msymbol(none) lwidth(thick) lcolor(`C2'))
		ci1opts(lcolor(`C1') lwidth(medthick))
		ci2opts(lcolor(`C2') lwidth(medthick) )
		title(`T', size(medlarge))
		ytitle("Grip Strength (kg)", size(medlarge))
		xtitle("")
		ylabel( , labsize(medlarge) format(%9.0f))
		xlabel(-30 "30" -20 "40" -10 "50" 0 "60" 10 "70" 20 "80" 30 "90" ,labsize(medlarge))
		yline(0)
		legend(ring(0) pos(7) cols(1) symxsize(5) order(3 4) lab(3 `"`l1'"') lab(4 `"`l2'"') title(`T'))
		name(`x', replace)
	; #delimit cr 
	}
graph combine poam4x grimage, cols(2) scheme(s1mono) ycommon name(grip, replace)



	//PoAm4x & Lung 
foreach x in poam4x grimage {
	 reg  htfev c.zar_`x' c.cage##sex c.zar_`x'#c.cage, robust 
	 margins, at(cage=(-30(10)30) zar_`x'=(-1 1))
	#delimit ;
	if `"`x'"'=="poam4x" { ; local T "DunedinPoAm4x"; local C1 "gold"; local C2 "orange_red"; local l1 "Slow"; local l2 "Fast"; }; 
	if `"`x'"'=="grimage" { ; local T "GrimAge"; local C1 "midblue"; local C2 "navy"; local l1 "Young"; local l2 "Old"; }; 
	marginsplot, scheme(plotplain) 
		plot1opts( msymbol(none) lwidth(thick) lcolor(`C1'))
		plot2opts( msymbol(none) lwidth(thick) lcolor(`C2'))
		ci1opts(lcolor(`C1') lwidth(medthick))
		ci2opts(lcolor(`C2') lwidth(medthick) )
		title(`T', size(medlarge))
		ytitle("Lung Function (FEV1 - Liters)", size(medlarge))
		xtitle("")
		ylabel( , labsize(medlarge) format(%9.0f))
		xlabel(-30 "30" -20 "40" -10 "50" 0 "60" 10 "70" 20 "80" 30 "90" ,labsize(medlarge))
		yline(0)
		legend(ring(0) pos(7) cols(1) symxsize(5) order(3 4) lab(3 `"`l1'"') lab(4 `"`l2'"') title(`T'))
		name(`x', replace)
	; #delimit cr 
	}
graph combine poam4x grimage, cols(2) scheme(s1mono) ycommon name(fev, replace)


	//PoAm4x & HbA1C 
foreach x in poam4x grimage {
	 reg  hba1cp c.zar_`x' c.cage##sex c.zar_`x'#c.cage#sex, robust 
	 margins, at(cage=(-30(10)30) zar_`x'=(-1 1))
	#delimit ;
	if `"`x'"'=="poam4x" { ; local T "DunedinPoAm4x"; local C1 "gold"; local C2 "orange_red"; local l1 "Slow"; local l2 "Fast"; }; 
	if `"`x'"'=="grimage" { ; local T "GrimAge"; local C1 "midblue"; local C2 "navy"; local l1 "Young"; local l2 "Old"; }; 
	marginsplot, scheme(plotplain) 
		plot1opts( msymbol(none) lwidth(thick) lcolor(`C1'))
		plot2opts( msymbol(none) lwidth(thick) lcolor(`C2'))
		ci1opts(lcolor(`C1') lwidth(medthick))
		ci2opts(lcolor(`C2') lwidth(medthick) )
		title(`T', size(medlarge))
		ytitle("HbA1C (mmol/mol)", size(medlarge))
		xtitle("")
		ylabel( , labsize(medlarge) format(%9.0f))
		xlabel(-30 "30" -20 "40" -10 "50" 0 "60" 10 "70" 20 "80" 30 "90" ,labsize(medlarge))
		yline(0)
		legend(ring(0) pos(5) cols(1) symxsize(5) order(3 4) lab(3 `"`l1'"') lab(4 `"`l2'"') title(`T'))
		name(`x', replace)
	; #delimit cr 
	}
graph combine poam4x grimage, cols(2) scheme(s1mono) ycommon name(hba1c, replace)

//PoAm4x & Sick 
foreach x in poam4x grimage {
	 logit  sick c.zar_`x' c.cage##sex c.zar_`x'#c.cage#sex, robust 
	 margins, at(cage=(-30(10)30) zar_`x'=(-1 1))
	#delimit ;
	if `"`x'"'=="poam4x" { ; local T "DunedinPoAm4x"; local C1 "gold"; local C2 "orange_red"; local l1 "Slow"; local l2 "Fast"; }; 
	if `"`x'"'=="grimage" { ; local T "GrimAge"; local C1 "midblue"; local C2 "navy"; local l1 "Young"; local l2 "Old"; }; 
	marginsplot, scheme(plotplain) 
		plot1opts( msymbol(none) lwidth(thick) lcolor(`C1'))
		plot2opts( msymbol(none) lwidth(thick) lcolor(`C2'))
		ci1opts(lcolor(`C1') lwidth(medthick))
		ci2opts(lcolor(`C2') lwidth(medthick) )
		title(`T', size(medlarge))
		ytitle("Prob. Prevalent Chronic Disease/Disability", size(medlarge))
		xtitle("")
		ylabel( , labsize(medlarge) format(%9.1f))
		xlabel(-30 "30" -20 "40" -10 "50" 0 "60" 10 "70" 20 "80" 30 "90" ,labsize(medlarge))
		legend(ring(0) pos(5) cols(1) symxsize(5) order(3 4) lab(3 `"`l1'"') lab(4 `"`l2'"') title(`T'))
		name(`x', replace)
	; #delimit cr 
	}
graph combine poam4x grimage, cols(2) scheme(s1mono) ycommon name(sick, replace)


//PoAm4x & Poor Health 
foreach x in poam4x grimage {
	 logit  lsrh c.zar_`x' c.cage##sex c.zar_`x'#c.cage#sex, robust 
	 margins, at(cage=(-30(10)30) zar_`x'=(-1 1))
	#delimit ;
	if `"`x'"'=="poam4x" { ; local T "DunedinPoAm4x"; local C1 "gold"; local C2 "orange_red"; local l1 "Slow"; local l2 "Fast"; }; 
	if `"`x'"'=="grimage" { ; local T "GrimAge"; local C1 "midblue"; local C2 "navy"; local l1 "Young"; local l2 "Old"; }; 
	marginsplot, scheme(plotplain) 
		plot1opts( msymbol(none) lwidth(thick) lcolor(`C1'))
		plot2opts( msymbol(none) lwidth(thick) lcolor(`C2'))
		ci1opts(lcolor(`C1') lwidth(medthick))
		ci2opts(lcolor(`C2') lwidth(medthick) )
		title(`T', size(medlarge))
		ytitle("Prob. Fair/Poor Self-rated Health", size(medlarge))
		xtitle("")
		ylabel( , labsize(medlarge) format(%9.1f))
		xlabel(-30 "30" -20 "40" -10 "50" 0 "60" 10 "70" 20 "80" 30 "90" ,labsize(medlarge))
		legend(ring(0) pos(5) cols(1) symxsize(5) order(3 4) lab(3 `"`l1'"') lab(4 `"`l2'"') title(`T'))
		name(`x', replace)
	; #delimit cr 
	}
graph combine poam4x grimage, cols(2) scheme(s1mono) ycommon name(lsrh, replace)




global cells "cd8t cd4t nk bcell mono gran plasmablast cd8pcd28ncd45ran cd8naive cd4naive"
global smoking "i.smk cpd1"	
foreach x in poam4x grimage {
	foreach y in mmgsdval srh htfev {
		//Base Model
		quietly reg z`y' zar_`x' c.cage##c.cage##sex , robust 
		matrix A = _b[zar_`x'] , _b[zar_`x'] - invttail(e(df_r),0.025)*_se[zar_`x'], _b[zar_`x'] + invttail(e(df_r),0.025)*_se[zar_`x'], 2*ttail(e(df_r),abs(_b[zar_`x']/_se[zar_`x'])), e(N)
		//Adjusted for Cell Counts 
		quietly reg z`y' zar_`x' c.cage##c.cage##sex $cells, robust 
		matrix A = A \ _b[zar_`x'] , _b[zar_`x'] - invttail(e(df_r),0.025)*_se[zar_`x'], _b[zar_`x'] + invttail(e(df_r),0.025)*_se[zar_`x'], 2*ttail(e(df_r),abs(_b[zar_`x']/_se[zar_`x'])), e(N)
		//Adjusted for Smoking 
		quietly reg z`y' zar_`x' c.cage##c.cage##sex $smoking , robust 
		matrix A = A \ _b[zar_`x'] , _b[zar_`x'] - invttail(e(df_r),0.025)*_se[zar_`x'], _b[zar_`x'] + invttail(e(df_r),0.025)*_se[zar_`x'], 2*ttail(e(df_r),abs(_b[zar_`x']/_se[zar_`x'])), e(N)	
		//Non-Smokers 
		quietly reg z`y' zar_`x' c.cage##c.cage##sex i.smk if smk<2 , robust 
		matrix A = A \ _b[zar_`x'] , _b[zar_`x'] - invttail(e(df_r),0.025)*_se[zar_`x'], _b[zar_`x'] + invttail(e(df_r),0.025)*_se[zar_`x'], 2*ttail(e(df_r),abs(_b[zar_`x']/_se[zar_`x'])), e(N)			
		matrix `y'=A , (1\2\3\4)
		matrix rownames `y' = `y' adj_cells adj_smoking nonsmokers 
		}	
	matrix `x' = (mmgsdval \htfev \ srh)  , (J(4,1,1) \ J(4,1,2) \ J(4,1,2))	
	}
matrix Fx = (poam4x, J(12,1,1)) \ (grimage,J(12,1,2))
matrix colnames Fx = r ll ul p N M Y X 
matrix list Fx 	
preserve 
	clear 
	svmat2 Fx, names(col)
	capture label drop M 
	label define M 1 "Base Model" 2 "Cell-count Adjusted" 3 "Smoking-Adjusted" 4 "Non-Smokers"
	label values M M 
	capture label drop Y 
	label define Y 1 "Grip Strength" 2 "Self-rated Health"
	label values Y Y
	capture label drop X 
	label define X 1 "DunedinPoAm4x" 2 "GrimAge" 	
	label values X X
	#delimit ; 
	twoway bar r X if X==1 & Y==2, bcolor(gold) barwidth(.8) by(M, legend(off) note(""))
		|| bar r X if X==2 & Y==2, bcolor(navy) barwidth(.8) 
		|| rcap ll ul X  if Y==2, lcolor(gs11) lwidth(medthick)
		xscale(range(.75 2.25))
		xlabel(1 "PoAm4x" 2 "GrimAge", labsize(medlarge) )
		ylabel(0(.1).3,labsize(medlarge) format(%9.2f))
		ytitle(Effect-size (Pearson r))
		xtitle("")
		scheme(plotplain)
		xsize(4) ysize(4)
		name(UndSocFx_SRH, replace) ; #delimit cr 	
	#delimit ; 
	twoway bar r X if X==1 & Y==1, bcolor(gold) barwidth(.8) by(M, legend(off) note(""))
		|| bar r X if X==2 & Y==1, bcolor(navy) barwidth(.8) 
		|| rcap ll ul X  if Y==1, lcolor(gs11) lwidth(medthick)
		xscale(range(.75 2.25))
		yscale(reverse)
		xlabel(1 "PoAm4x" 2 "GrimAge", labsize(medlarge) )
		ylabel(-.2(.1).1,labsize(medlarge) format(%9.2f))
		ytitle(Effect-size (Pearson r))
		xtitle("")
		scheme(plotplain)
		xsize(4) ysize(4)
		name(UndSocFx_Grip, replace) ; #delimit cr 		
restore 

foreach x in poam4x grimage {
	poisson lsrh zar_`x' c.cage##c.cage##sex , robust irr 
	}


//Grip Strength
capture drop cage 
gen cage = age-60
foreach x in poam4x grimage {
	capture drop zar_`x'
	egen zar_`x' = std(ar_`x')
reg mmgsdval zar_`x' c.zar_`x'#sex c.cage##c.cage##sex if age>50, robust 
margins, dydx(zar_`x') over(sex)
#delimit ;
if `"`x'"'=="poam4x" { ; local T "DunedinPoAm4x"; local C "gold"; }; 
if `"`x'"'=="grimage" { ; local T "GrimAge"; local C "navy"; }; 
marginsplot, scheme(plotplain) 
	plot1opts(connect(none) msymbol(O) msize(vlarge) mcolor(`C'))
	ci1opts(lcolor(`C') lwidth(medthick))
	title(`T', size(medlarge))
	ytitle(Effect-size, size(medlarge))
	xtitle("")
	ylabel(, labsize(medlarge))
	xlabel( ,labsize(medlarge))
	yline(0)
	xscale(range(.75 2.25))
	name(`x', replace)
; #delimit cr 
}
graph combine poam4x grimage, cols(2) scheme(s1mono) ycommon 




//Effect-sizes 
capture drop cage 
gen cage = age-60
capture drop zsrh 
egen zsrh = std(srh)
foreach x in poam4x grimage {
	capture drop zar_`x'
	egen zar_`x' = std(ar_`x')
reg zsrh zar_`x' c.zar_`x'#sex c.cage##c.cage##sex, robust 
margins, dydx(zar_`x') over(sex)
#delimit ;
if `"`x'"'=="poam4x" { ; local T "DunedinPoAm4x"; local C "gold"; }; 
if `"`x'"'=="grimage" { ; local T "GrimAge"; local C "navy"; }; 
marginsplot, scheme(plotplain) 
	plot1opts(connect(none) msymbol(O) msize(vlarge) mcolor(`C'))
	ci1opts(lcolor(`C') lwidth(medthick))
	title(`T', size(medlarge))
	ytitle(Effect-size, size(medlarge))
	xtitle("")
	ylabel(0(.1).4, labsize(medlarge))
	xlabel( ,labsize(medlarge))
	yline(0)
	xscale(range(.75 2.25))
	name(`x', replace)
; #delimit cr 
}
graph combine poam4x grimage, cols(2) scheme(s1mono) ycommon name(srh, replace)





/*
//Grip Strength
foreach x in poam4x grimage {
	capture drop zar_`x'
	egen zar_`x' = std(ar_`x')
#delimit ; 
reg mmgsdval zar_`x' c.zar_`x'#sex 
					c.zar_`x'#c.cage
					c.zar_`x'#c.cage#sex 
				/*	c.zar_`x'#c.cage#c.cage 
					c.zar_`x'#c.cage#c.cage#sex 
				*/	c.cage##c.cage##sex , robust ; #delimit cr 
margins, dydx(zar_`x') over(sex) at(cage=(-20(10)30))
#delimit ;
if `"`x'"'=="poam4x" { ; local T "DunedinPoAm4x"; local C "gold"; }; 
if `"`x'"'=="grimage" { ; local T "GrimAge"; local C "navy"; }; 
marginsplot, scheme(plotplain) 
	plot1opts( msymbol(O) msize(large) mcolor(`C') lcolor(`C'))
	plot2opts( msymbol(T) msize(large) mcolor(`C') lcolor(`C'))
	ci1opts(lcolor(`C') lwidth(medthick))
	ci2opts(lcolor(`C') lwidth(medthick))
	title(`T', size(medlarge))
	ytitle(Effect-size, size(medlarge))
	xtitle("")
	ylabel(, labsize(medlarge))
	xlabel(-20 "40" -10 "50" 0 "60" 10 "70" 20 "80 "30 "90" ,labsize(medlarge))
	yline(0)
	xscale(range(.75 2.25))
	name(`x', replace)
; #delimit cr 
}
graph combine poam4x grimage, cols(2) scheme(s1mono) ycommon 

*/


preserve 
quietly sum mmgsdval if sex == 1 
gen Z = (mmgsdval-r(mean))/r(sd) if sex==1 
quietly sum mmgsdval if sex == 2 
replace Z = (mmgsdval-r(mean))/r(sd) if sex == 2
foreach x in poam4x grimage {
	capture drop zar_`x'
	egen zar_`x' = std(ar_`x')
quietly reg Z zar_`x' c.zar_`x'#sex c.cage##c.cage##sex, robust 
margins, dydx(zar_`x') over(sex)
}

  
