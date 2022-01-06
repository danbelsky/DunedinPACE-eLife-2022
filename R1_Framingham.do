global db3275 "db3275"
//global db3275 "danielbelsky"


//***************************************************************************//
//***************************************************************************//
//Survival
//***************************************************************************//
//***************************************************************************//
//Death 
import delimited "/Users/$db3275/OneDrive - cumc.columbia.edu/Framingham/Data/Mortality/DataPull_DunedinPACE/vr_survdth_2018_a_1268s.csv", delim(comma) varn(1) clear
drop shareid 
save death, replace 
//CVD
import delimited "/Users/$db3275/OneDrive - cumc.columbia.edu/Framingham/Data/Mortality/DataPull_DunedinPACE/vr_survcvd_2018_a_1267s.csv", delim(comma) varn(1) clear
drop shareid 
save cvd, replace 
//Stroke
	//V1
import delimited using "/Users/db3275/OneDrive - cumc.columbia.edu/Framingham/Data/Mortality/DataPull_DunedinPACE/phs000007.v32.pht006024.v4.p13.c1.vr_svstktia_2018_a_1270s.HMB-IRB-MDS.txt", clear varn(11) delim(tab)
drop shareid 
keep if idtype==1 // keep only Offspring Cohort 
save stroke1 , replace 
	//V2
import delimited using "/Users/db3275/OneDrive - cumc.columbia.edu/Framingham/Data/Mortality/DataPull_DunedinPACE/phs000007.v32.pht006023.v4.p13.c1.vr_svstk_2018_a_1269s.HMB-IRB-MDS.txt", clear varn(11) delim(tab)
drop shareid 
keep if idtype==1 // keep only Offspring Cohort 
save stroke2 , replace 

//Assemble Survival Dataset 
use death, clear 
merge 1:1 dbgap_subject_id using cvd, nogen
merge 1:1 dbgap_subject_id using stroke1, nogen
merge 1:1 dbgap_subject_id using stroke2, nogen
foreach x in cvddeath chddeath datedth dthrvwd lastcon lastsoe { 
	destring `x', replace force 
}
preserve 
	import delimited using "/Users/$db3275/OneDrive - cumc.columbia.edu/Framingham/Data/Mortality/framingham_mortality_MK07052021.csv", clear varn(1) delim(comma)
	foreach x in chddeath cvddeath datedth dthrvwd lastatt lastcon lastsoe {
	destring `x', replace force 
	}
	drop if cohort == "GEN3"
	save mkframinghamoffspring, replace
restore 
merge 1:1 dbgap_subject_id using mkframinghamoffspring,
keep if _merge ==3 | _merge==2 
drop _merge 

//Residualize clocks for age 
gen horvath = dnamage 
gen hannum = dnamagehannum
gen phenoage = dnamphenoage 
gen grimage = dnamgrimage 
destring age8, replace force 
foreach x in horvath hannum phenoage grimage{ 
	reg `x' age8 
	predict ar_`x', r 
	egen zar_`x' = std(ar_`x')
	}
/*
//Residualize GrimAge for age and sex 
reg grimage age8 sex
predict ar_grimage, r 
egen zar_grimage = std(ar_grimage)
*/
//Standardize PoA measures 
egen zpoam = std(dunedin_poam38)
egen zpace = std(dunedin_poam45)
//Smoking covariates 
gen smk = 0 if smoking_status == "Never Smoker"
replace smk = 1 if smoking_status == "Former Smoker"
replace smk = 2 if smoking_status == "Current Smoker"
label var smk "Smoking Status"
capture label drop smk 
label define smk 0 "Never Smoker" 1 "Former Smoker" 2 "Current Smoker"
label values smk smk 
gen cpd = smoking_quantity
destring cpd, replace force
//Clock Var macro
global clocks "zpace zpoam zar_horvath zar_hannum zar_phenoage zar_grimage" 
//Cell control macro
global cells "cd4t cd8t nk bcell mono gran cd8pcd28ncd45ran cd8naive cd4naive plasmablast"

//Survival Variables 
destring date8, replace force 
	//Time from Visit 8 baseline to CVD 
capture drop T_cvd* 
gen T_cvd = (cvddate - date8) / 365
	//Time from Visit 8 baseline to CHF 
capture drop T_chf* 
gen T_chf = (chfdate - date8) / 365
	//Time from Visit 8 baseline to CHD 
capture drop T_chd* 
gen T_chd = (chddate - date8) / 365
	//Time from Visit 8 baseline to Stroke
gen T_stroke = (strokedate - date8) / 365
	//Time from Visit 8 baseline to Stroke / TIA
gen T_stroke_tia = (stroketiadate - date8) / 365
	//Time from Visit 8 baseline to DEATH 
capture drop T_death* 
gen T_death = (datedth - date8) / 365
replace T_death = (lastcon - date8) / 365 if T_death ==.
capture drop death 
gen death=(dthrvwd==1)

tab packs_set , gen(B)
	
save "/Users/$db3275/OneDrive - cumc.columbia.edu/Projects/mPoA/DunedinPoAm45/eLife/R1/R1_Framingham_Survival.dta", replace 

//***************************************************************************//
//***************************************************************************//
//DISABILITY 
//***************************************************************************//
//***************************************************************************//
import delimited using "/Users/$db3275/OneDrive - cumc.columbia.edu/Framingham/Data/Mortality/DataPull_DunedinPACE/ex1_8s_full.txt", clear delim(tab) varn(11) 

//Rosow-Breslau ADLs
//ARE YOU ABLE TO DO HEAVY WORK AROUND THE HOUSE, LIKE SHOVELING SNOW OR WASHING WINDOWS, WALLS, OR FLOORS WITHOUT HELP?
tab h468
//ARE YOU ABLE TO WALK HALF A MILE WITHOUT HELP? (ABOUT 4-6 BLOCKS)
tab h469
//ARE YOU ABLE TO WALK UP AND DOWN ONE FLIGHT OF STAIRS WITHOUT HELP?
tab h470
capture drop rosbres 
egen rosbres=rowtotal(h468 h469 h470)
recode rosbres (0=3) (1=2) (2=1) (3=0)
rename rosbres rosbres8 


//KATZ ADLs
// No help needed / uses device / human assistance / dependent / do not do 
//CAN YOU GET DRESSED (UNDRESSING REDRESSING) INDEPENDENTLY OR DO YOU NEED HUMAN ASSISTANCE OR THE USE OF A DEVICE (SUCH AS VELCRO, ELASTIC LACES)?
tab h474
//CAN YOU BATHE (INCLUDING GETTING IN AND OUT OF THE TUB/SHOWER) INDEPENDENTLY OR DO YOU NEED HUMAN ASSISTANCE OR USE OF A DEVICE (BATH CHAIR, LONG HANDLED SPONGE, HAND HELD SHOWER, SAFETY BARS)?
tab h475
// CAN YOU EAT INDEPENDENTLY OR DO YOU NEED HUMAN ASSISTANCE OR THE USE OF A DEVICE (SUCH AS ROCKING KNIFE, SPORK, LONG STRAW, PLATE GUARD)?
tab h476
//CAN YOU TRANSFER (GETTING IN AND OUT OF A CHAIR) INDEPENDENTLY OR DO YOU NEED HUMAN ASSISTANCE OR THE USE OF A DEVICE (SUCH AS, SLIDING BOARD, GRAB BARS, SPECIAL SEAT)?
tab h477
//CAN YOU DO TOILETING ACTIVITIES (USING BATHROOM FACILITIES AND HANDLE CLOTHING) INDEPENDENTLY OR DO YOU NEED HUMAN ASSISTANCE OR THE USE OF A DEVICE (SUCH AS, SPECIAL TOILET SEAT, COMMODE)?
tab h478

capture drop katz
egen katz = rowtotal(h474-h478)
tab katz 
rename katz katz8 

//NAGI ADLs
// No difficulty / little difficulty/ some difficulty / lot of difficulty / unable to do + MD orders not to do
//PULLING OR PUSHING LARGE OBJECTS LIKE A LIVING ROOM CHAIR
tab h569
//EITHER STOOPING, CROUCHING, OR KNEELING
tab h570
//REACHING OR EXTENDING ARMS BELOW SHOULDER LEVEL
tab h571
//REACHING OR EXTENDING ARMS ABOVE SHOULDER LEVEL
tab h572
//EITHER WRITING, OR HANDLING, OR FINGERING SMALL OBJECTS
tab h573
//STANDING IN ONE PLACE FOR LONG PERIODS, SAY 15 MINUTES
tab h574 
//SITTING FOR LONG PERIODS, SAY 1 HOUR
tab h575 
//LIFTING OR CARRYING WEIGHTS UNDER 10 POUNDS (LIKE A BAG OF POTATOES)
tab h576
//LIFTING OR CARRYING WEIGHTS OVER 10 POUNDS (LIKE A VERY HEAVY BAG OF GROCERIES)
forvalues v=69(1)75{
	recode h5`v' (6=.) (5=4)
	}
capture drop nagi 
egen nagi = rowtotal(h569-h575)
tab nagi 
forvalues v=69(1)75{
	capture drop dh5`v'
	recode h5`v' (1/3=0) (4=1), gen(dh5`v')
	}
capture drop dnagi 
egen dnagi = rowtotal(dh569-dh575)
tab dnagi 
rename nagi nagi8 
rename dnagi dnagi8 

keep dbgap_subject_id rosbres8 katz8 nagi8 dnagi8 	
save temp, replace 


//WAVE 9 ADLs 
import delimited using "/Users/$db3275/OneDrive - cumc.columbia.edu/Framingham/Data/Mortality/DataPull_DunedinPACE/e_exam_ex09_1b_0844s_full.txt", clear delim(tab) varn(11) 

//ROSOW-BRESLAU
capture drop rosbres 
egen rosbres=rowtotal(j609 j610 j611)
recode rosbres (0=3) (1=2) (2=1) (3=0)
tab rosbres 
capture drop rosbres9 
rename rosbres rosbres9

//KATZ ADLs 
foreach var of varlist j612-j616{
	tab `var'
	}
capture drop katz
egen katz = rowtotal(j612-j616)
tab katz 
rename katz katz9 

//NAGI (no/little/some/lot/unable)
forvalues v=597(1)605{
	tab j`v'
	}
capture drop nagi
egen nagi = rowtotal(j597-j605)
tab nagi 
forvalues v=597(1)605{
	capture drop dj5`v'
	recode j`v' (1/3=0) (4=1), gen(dj`v')
	}
capture drop dnagi 
egen dnagi = rowtotal(dj597-dj605)
tab dnagi 
rename nagi nagi9
rename dnagi dnagi9

keep dbgap_subject_id rosbres9 katz9 nagi9 dnagi9 	
merge 1:1 dbgap_subject_id using temp , nogen 
save temp, replace 

use "/Users/$db3275/OneDrive - cumc.columbia.edu/Projects/mPoA/DunedinPoAm45/eLife/R1/R1_Framingham_Survival.dta", clear 
keep dbgap_subject_id sex age1 age2 age8 zpace idtype cohort

merge 1:1 dbgap_subject_id using temp, nogen 

keep if zpace!=. 

//ROSOW-BRESLAU
tab rosbres8 
tab rosbres9 
foreach x in 8 9 {
	recode rosbres`x' (3=2)
	}
tab rosbres*
//KATZ
tab katz8 
tab katz9
foreach x in 8 9 {
	recode katz`x' (3/100=2)
	}
tab katz*
//NAGI
tab dnagi8 
tab dnagi9 
foreach x in 8 9 {
	recode dnagi`x' (3/100=2)
	}
tab dnagi* 

foreach y in rosbres katz dnagi { 
	capture drop i`y' 
	gen i`y' = `y'9-`y'8
	replace i`y'=0 if i`y'<0
	}
	
save "/Users/$db3275/OneDrive - cumc.columbia.edu/Projects/mPoA/DunedinPoAm45/eLife/R1/R1_Framingham_ADL.dta", replace 

use "/Users/$db3275/OneDrive - cumc.columbia.edu/Projects/mPoA/DunedinPoAm45/eLife/R1/R1_Framingham_Survival.dta", clear
merge 1:1 dbgap_subject_id using "/Users/$db3275/OneDrive - cumc.columbia.edu/Projects/mPoA/DunedinPoAm45/eLife/R1/R1_Framingham_ADL.dta", nogen
save "/Users/$db3275/OneDrive - cumc.columbia.edu/Projects/mPoA/DunedinPoAm45/eLife/R1/R1_Framingham.dta", replace 



