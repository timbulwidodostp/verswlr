* verswlr.ado

*! verswlr v1.0.0 Theodore Karrison 24dec2014

** Versatile Weighted Logrank Test

** Syntax:  verswlr group [if] [in] [, options]
** by is allowed

program verswlr, rclass byable(recall)
 version 13.1
 syntax varlist(max=1 numeric) [if] [in] [,rho1(real 0) gamma1(real 0) ///
  rho2(real 1) gamma2(real 0) rho3(real 0) gamma3(real 1)]
 
   local res fh11 fh22 fh33 maxz pval sampsize
   tempname `res' 
   tempvar group 
   gen `group'=`1'
   
   tempname sigma11 sigma22 sigma33 cov12 corr12 cov13 corr13 cov23 corr23 ///
    fhmax tv1 tv2 tv3 tv4 tv5 tv6 tv7 tv8 ulr vlr uearly vearly ulate vlate u v 
   
   marksample touse
   quietly count if `touse'
   if `r(N)'==0 {
            error 2000
   }
   scalar `sampsize'=`r(N)'
   
   * logrank (Fleming-Harrington 0,0) or non-default specification 
   quietly sts test `group' if `touse', fh(`=`rho1'' `=`gamma1'') mat(`ulr' `vlr')
   scalar `fh11'=r(chi2)
   scalar `sigma11'=sqrt(`vlr'[1,1])
   
   * early (Fleming-Harrington 1,0) or non-default specification
   quietly sts test `group' if `touse', fh(`=`rho2'' `=`gamma2'') mat(`uearly' `vearly')
   scalar `fh22'=r(chi2)
   scalar `sigma22'=sqrt(`vearly'[1,1])
   
   * late (Fleming-Harrington 0,1) or non-default specification
   quietly sts test `group' if `touse', fh(`=`rho3'' `=`gamma3'') mat(`ulate' `vlate')
   scalar `fh33'=r(chi2)
   scalar `sigma33'=sqrt(`vlate'[1,1])
   
   * Max statistic
   scalar `fhmax'=max(`fh11',`fh22',`fh33')
   scalar `maxz'=sqrt(`fhmax')
   
   * obtain pval for maxz
   ** covariance term for fh11 and fh22 statistics
   quietly sts test `group' if `touse', ///
    fh(`=(`rho1'+`rho2')/2' `=(`gamma1'+`gamma2')/2') mat(`u' `v')
   scalar `cov12'=`v'[1,1]
   scalar `corr12'=`cov12'/(`sigma11'*`sigma22')
   
   ** covariance term for fh11 and fh33 statistics
   quietly sts test `group' if `touse', /// 
    fh(`=(`rho1'+`rho3')/2' `=(`gamma1'+`gamma3')/2') mat(`u' `v')
   scalar `cov13'=`v'[1,1]
   scalar `corr13'=`cov13'/(`sigma11'*`sigma33')
   
   ** covariance term for fh22 and fh33 statistics
   quietly sts test `group' if `touse', ///
    fh(`=(`rho2'+`rho3')/2' `=(`gamma2'+`gamma3')/2') mat(`u' `v')
   scalar `cov23'=`v'[1,1]
   scalar `corr23'=`cov23'/(`sigma22'*`sigma33')
   
   trivarn -`maxz' -`maxz' -`maxz' `corr12' `corr13' `corr23'
   scalar `tv1'=r(tvr)
   trivarn -`maxz' -`maxz' `maxz' `corr12' `corr13' `corr23'
   scalar `tv2'=r(tvr)
   trivarn -`maxz' `maxz' -`maxz' `corr12' `corr13' `corr23'
   scalar `tv3'=r(tvr)
   trivarn -`maxz'  `maxz'  `maxz' `corr12' `corr13' `corr23'
   scalar `tv4'=r(tvr)
   trivarn `maxz' -`maxz' -`maxz' `corr12' `corr13' `corr23'
   scalar `tv5'=r(tvr)
   trivarn `maxz' -`maxz' `maxz' `corr12' `corr13' `corr23'
   scalar `tv6'=r(tvr)
   trivarn `maxz' `maxz' -`maxz' `corr12' `corr13' `corr23'
   scalar `tv7'=r(tvr)
   trivarn `maxz' `maxz' `maxz' `corr12' `corr13' `corr23'
   scalar `tv8'=r(tvr)
   scalar `pval'=1-(`tv1'-`tv2'-`tv3'+`tv4'-`tv5'+`tv6'+`tv7'-`tv8')
 
 display as txt "  "
 display as result "Versatile Logrank Test"
 display as txt "  "
 display as txt "Total Sample Size: n = " `sampsize'
 display as txt "  "
 display as txt "Rho-Gamma tests:"
 display as txt "fh(rho1,gamma1) = fh(" `rho1' "," `gamma1' "): Chi-square = " `fh11'
 display as txt "fh(rho2,gamma2) = fh(" `rho2' "," `gamma2' "): Chi-square = " `fh22'
 display as txt "fh(rho3,gamma3) = fh(" `rho3' "," `gamma3' "): Chi-square = " `fh33'
 display as txt "  "
 display as txt "Maximum Test:"
 display as txt "Max abs z: " `maxz'
 display as txt "Two sided p-value: " `pval'
 
 foreach r of local res {
  return scalar `r' = ``r''
 }
 
 end 
 

program trivarn, rclass
 version 13.1
 tempname ier h1 h2 h3 r12 r13 r23 hh rr h12 h13 h122 h132 tv rr12 rr13 ///
  del fac rr122 rr133 f1 f2 f3 hp1 hp2 
 tempname x w
 
 matrix `x' = (0.04691008, 0.23076534, 0.5, 0.76923466, 0.95308992)
 matrix `w' = (0.018854042, 0.038088059, 0.0452707394, 0.038088059, 0.018854042)
 
 scalar `ier'=1
 
 scalar `h1'=`1'
 scalar `h2'=`2'
 scalar `h3'=`3'
 scalar `r12'=`4'
 scalar `r13'=`5'
 scalar `r23'=`6'


 forvalues k =1/3 {
  if abs(`r23')>=abs(`r12') & abs(`r23')>=abs(`r13') {
   continue, break
  }
  scalar `hh'=`h1'
  scalar `h1'=`h2'
  scalar `h2'=`h3'
  scalar `h3'=`hh'
  scalar `rr'=`r12'
  scalar `r12'=`r23'
  scalar `r23'=`r13'
  scalar `r13'=`rr'
 }
 
 scalar `h12'=`h1'*`h2'
 scalar `h13'=`h1'*`h3'
 scalar `h122'=(`h1'*`h1'+`h2'*`h2')*0.5
 scalar `h132'=(`h1'*`h1'+`h3'*`h3')*0.5
 scalar `tv'=0

 forvalues i=1/5 {
  scalar `rr12'=`r12'*`x'[1,`i']
  scalar `rr13'=`r13'*`x'[1,`i']
  scalar `del'=1.0-`rr12'*`rr12'-`rr13'*`rr13'-`r23'*`r23'+2*`rr12'*`rr13'*`r23'
  if `del'<=0.0 {
   continue, break
  }
  scalar `fac'=sqrt(`del')
  scalar `rr122'=1.0-`rr12'*`rr12'
  scalar `rr133'=1.0-`rr13'*`rr13'
  scalar `f1'=`rr13'-`r23'*`rr12'
  scalar `f2'=`r23'-`rr12'*`rr13'
  scalar `f3'=`rr12'-`r23'*`rr13'
  
  scalar `hp1'=(`h3'*`rr122'-`h1'*`f1'-`h2'*`f2')/(`fac'*sqrt(`rr122'))
  scalar `hp2'=(`h2'*`rr133'-`h1'*`f3'-`h3'*`f2')/(`fac'*sqrt(`rr133'))
  
  scalar `tv'=`tv'+`w'[1,`i']* ///
  (exp((`rr12'*`h12'-`h122')/`rr122')/sqrt(`rr122'))*(1.0-normal(`hp1'))*`r12'
   
  scalar `tv'=`tv'+`w'[1,`i']* ///
   (exp((`rr13'*`h13'-`h132')/`rr133')/sqrt(`rr133'))*(1.0-normal(`hp2'))*`r13'
  return scalar tvr=`tv'
 }
 
 if `del'>0.0 {
  scalar `ier'=0
  scalar `tv'=(1.0-normal(`h1'))*binormal(-`h2',-`h3',`r23')+`tv'
  return scalar tvr=`tv'
 }
 
 end
