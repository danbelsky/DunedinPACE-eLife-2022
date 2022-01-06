
//NAS 
preserve
clear 
input X hr ll ul M Y  
1 1.26 1.14 1.40 1 1
1 1.24 1.11 1.38 2 1
1 1.24 1.11 1.38 3 1
1 1.16 1.12 1.20 1 2  
1 1.16 1.12 1.20 2 2
1 1.14 1.04 1.18 3 2
1 1.23 1.07 1.42 1 3  
1 1.21 1.05 1.40 2 3
1 1.20 1.04 1.36 3 3
2 1.29 1.16 1.45 1 1
2 1.28 1.13 1.44 2 1
2 1.26 1.12 1.42 3 1
2 1.15 1.11 1.20 1 2  
2 1.15 1.10 1.19 2 2
2 1.16 1.12 1.21 3 2
2 1.19 1.03 1.38 1 3  
2 1.18 1.01 1.38 2 3
2 1.15 0.99 1.34 3 3
end 
capture label drop X 
label define X 1 "DunedinPoAm4x" 2 "DunedinPoAm3x"
label values X X
capture label drop Y 
label define Y 1 "Mortality" 2 "Prevalent Chronic Disease" 3 "Incident Chronic Disease"
label values Y Y 
capture label drop M 
label define M 1 `""Base" "Model""' 2 `""Cell-count" "Adjusted""' 3 `""Smoking" "Adjusted""' 
label values M M 

drop if Y == 2 
drop if M ==2 
bys X: gen N = _n
replace N = N + 1 if Y ==3 
 
#delimit ;
twoway rcap ll ul N if X==1 , lcolor(gold)  lwidth(medthick)
	|| rcap ll ul N if X==2 , lcolor(orange_red) lwidth(medthick)
	|| scatter hr N if X==1 & M==1, msymbol(O) mfcolor(gold) mlcolor(orange_red) msize(vlarge)
	|| scatter hr N if X==1 & M==2, msymbol(T) mfcolor(gold) mlcolor(orange_red) msize(vlarge)
	|| scatter hr N if X==1 & M==3, msymbol(D) mfcolor(gold) mlcolor(orange_red) msize(vlarge)
	|| scatter hr N if X==2 & M==1, msymbol(O) mcolor(orange_red) msize(vlarge)
	|| scatter hr N if X==2 & M==2, msymbol(T) mcolor(orange_red) msize(vlarge)	
	|| scatter hr N if X==2 & M==3, msymbol(D) mcolor(orange_red) msize(vlarge)	
	yline(1)
	xlabel(1.5 "Mortality" 4.5 `""Incident" "Chronic Disease""', labsize(medlarge))
	ylabel(1(.1)1.5,labsize(medlarge) format(%9.2f))
	xtitle("")
	ytitle(Hazard Ratio)
	xtitle("")
	xscale(range(.7 3.25))
	by(X, legend(off) note("")) 
	scheme(plotplain)
	name(NAS, replace )
; #delimit cr 
restore 
