********************************************************************************************************
/* THIS CODE SELECTS THE US INPUT OUTPUT MATRIX AND CONSTRUCTS ITS TECHNICAL COEFFICIENT MATRIX */
********************************************************************************************************
use "WIOT2013_October16_ROW.dta", clear
keep IndustryCode IndustryDescription Country RNr Year vUSA* TOT
keep if Country=="USA"
drop if IndustryCode=="U"
*drop vUSA57 vUSA58 vUSA59 vUSA60 vUSA61
* We stop at 55 because 56 sums to 0
mkmat vUSA1-vUSA55, matrix(wiot)
mkmat TOT, matrix(x)
mat li x
matrix A = J(55,55,0)
matrix list A 
forvalues i = 1/55 {
	forvalues j = 1/55 {
		 matrix A[`i',`j']= wiot[`i',`j']/x[`j',1]
	}
}
        
matrix list A
drop vUSA*
rename TOT Output
svmat A, names(SIC)
saveold "wiot_US_2013.dta",version(12) replace
********************************************************************************************************
