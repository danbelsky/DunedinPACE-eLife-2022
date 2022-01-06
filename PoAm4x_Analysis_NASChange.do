
//Code to confirm change over time in DunedinPoAm variables 
import delimited using "/Users/db3275/OneDrive - cumc.columbia.edu/NAS/DNAm_XuGaoAugust2020/fixed effects200911.csv", delim(comma) varn(1) clear

gen D = date(date,"MDY")
egen B = min(D), by(id)
gen T = D-B 
gen baseline = T==0

gen T5=T/(365*5)
gen T1=T/(365)

tab T5

capture drop bage
capture drop cage5
egen bage =min(age), by(id)
gen cage5=(bage-65)/5

quietly sum mpoa if baseline==1 
gen zmpoa = (mpoa-r(mean))/r(sd)
quietly sum mpoa45_new if baseline==1 
gen zpoam4x = (mpoa45_new-r(mean))/r(sd)

//Test change per 5y
foreach x in mpoa mpoa45_new { 
	xtreg `x' T5  , fe i(id)	
	}
//Quantify in SD units	
foreach x in zmpoa zpoam4x { 
	xtreg `x' T5  , fe i(id)	
	}	
	
//Test if change varies by age at baseline (no)	
foreach x in mpoa mpoa45_new { 
	xtreg `x' c.T5##c.cage5  , fe i(id)	
	}
		
xtreg mpoa45 c.T1  , fe i(id)
margins, at(T=(0(2)12) )
	
#delimit ; 	
twoway scatter mpoa45 age if baseline==1, msymbol(O) mfcolor(navy%30) mlcolor(navy%1)	
	|| lfit mpoa45 age if baseline==1, lcolor(dknavy) lwidth(thick) lpattern(dash) 
legend(off)
scheme(plotplain)
name(nasage, replace)
; #delimit cr

spagplot mpoa45 age if cage5==1, id(id) nofit

#delimit ; 	
twoway scatter mpoa45 age if baseline==1, msymbol(O) mfcolor(navy%30) mlcolor(navy%1)	
	|| scatter mpoa45 age if baseline==1 & age==70, msymbol(O) mcolor(red)  
	|| scatter mpoa45 age if cage5==1 , msymbol(D) mcolor(red)  
	|| scatter mpoa45 age if baseline==1 & age==70, msymbol(O) mcolor(red)  
	|| lfit mpoa45 age, lcolor(dknavy) lwidth(thick) lpattern(dash) 
scheme(plotplain)
name(nasage_long, replace)
; #delimit cr

	
	
clear 
import delimited using "/Users/db3275/OneDrive - cumc.columbia.edu/NAS/DNAm_XuGaoAugust2020/mortality200911.csv", delim(comma) clear 	

capture drop zpoam4x3
egen poam4x3 = cut(mpoa45_new), group(3)
recode poam4x3 (0=1) (1=2) (2=3) 

egen zpoam4x= std(mpoa45_new)
capture drop zpoam4xhal
recode zpoam4x (-10/-1=1) (-1/1=2) (1/10=3), gen(zpoam4xhal)

stset time_fu, failure(dead==1) scale(365)
//local Z "zpoam4x3"
local Z "zpoam4xhal"
#delimit ;
sts graph, survival by(`Z') censored(single) censopts(lcolor(gs6))
risktable xlabel(0(3)15) 
risktable(, order(1 "Slow DunedinPoAm4x" 2 "Average DunedinPoAm4x" 3 "Fast DunedinPoAm4x") failevents title(Number At Risk (Deaths)))
title("")
xtitle(Analysis Time (years))
ytitle(Survival)
legend(pos(3) cols(1) symxsize(5) lab(1 "Slow") lab(2 "Average") lab(3 "Fast")
	title(DunedinPoAm4x, size(medsmall)) region(lcolor(white)) )
plot1opts(lwidth(medthick) lcolor(blue))
plot2opts(lwidth(medthick) lcolor(gs10))
plot3opts(lwidth(medthick) lcolor(cranberry))
ylabel(,angle(horiz) nogrid labsize(medlarge) format(%9.2f))
xlabel(,labsize(medlarge))
graphregion(color(white)) plotregion(color(white))
xsize(7) ysize(4)
name(naswcens4x,replace)
; #delimit cr

graph export "/Users/db3275/OneDrive - cumc.columbia.edu/Projects/mPoA/DunedinPoAm45/Code/Figures/NAS_PoAm4x_Survival.pdf", replace	


