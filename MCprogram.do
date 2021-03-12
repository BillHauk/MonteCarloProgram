version 7.0
clear
set more off
set mem 50m

set matsize 200
*true value for alpha
global alpha 0.27
*true value for beta
global beta 0.27

*number of countries in Monte Carlo sample draw
global N 69
*Number of Monte Carlo draws
global drawnum 1000

*measurement error for yini variable
global shockyini 0.1
*measurement error for sk variable
global shocksk 0.1
*measurement error for sh variable
global shocksh 0.1
*measurement error for ndg variable
global shockndg 0.1
*measurement error autocorrelation term
global errrho 0.1
*correlation scalar for fixed effects and sk (data has 0.6395 average)
global skfecorr 1
*correlation scalar for fixed effects and sh (data has 0.8694 average)
global shfecorr 1
*correlation scalar for fixed effects and ndg (data has -0.5842 average)
global ndgfecorr 1
*correlation scalar for fixed effects and yini (data has 0.9377 average)
global yinifecorr 1
*Choose uniform (dataqual=0) or country specific (dataqual=1) measurement errors
global dataqual 1
*Error term correlation with sk
global skerrcorr 0.4
*Error term correlation with sh
global sherrcorr 0.4
*Error term correlation with ndg
global ndgerrcorr -0.4

use MCdatabal

* Take the raw data and generate residuals and country fixed effects

global lambda=(1-$alpha-$beta)*(0.08)
global multiple=1-exp(-$lambda*5)

sort tavawac period5y
generate rgdpch_end=rgdpch_ini[_n+1] if tavawac==tavawac[_n+1]
drop if period5y<3 | period5y==11
generate lki=log(ki2)
generate lenroll=log(senroll)
generate lgpop=log(gpop+0.07)
generate lini=log(rgdpch_ini)
generate lend=log(rgdpch_end)


xi: xtreg lend lki lenroll lgpop lini i.period5y, fe
predict resids, e
predict fixedeffect, u
replace fixedeffect=fixedeffect+_b[_cons]


keep  tavawac period5y lini lenroll lki lgpop dqual fixedeffect resids
sort tavawac period5y
reshape wide lini lenroll lki lgpop resids, i(tavawac) j(period5y)
save baselinedatawide.dta, replace

global t=1

*Each one of these loops generates coefficient values using fixed effects (i.e. the true model)
*shocked fixed effects, and shocked ols and saves them as data
*Note that the coefficients generated are of the form alpha/(1-alpha-beta), etc.

while $t==1{
	do MCgenerate
	keep  country period yini sk sh ndg yend sndg ssh ssk syini syend
	* Define variables and tsset command for Arellano-Bond
	sort period country
	by period: gen meanssk=sum(ssk)
	by period: gen meanssh=sum(ssh)
	by period: gen meansndg=sum(sndg)
	by period: gen meansyini=sum(syini)
	by period: gen meansyend=sum(syend)
	by period: replace meanssk=meanssk[_N]/$N
	by period: replace meanssh=meanssh[_N]/$N
	by period: replace meansndg=meansndg[_N]/$N
	by period: replace meansyini=meansyini[_N]/$N
	by period: replace meansyend=meansyend[_N]/$N
	gen sskd=ssk-meanssk
	gen sshd=ssh-meanssh
	gen sndgd=sndg-meansndg
	gen syinid=syini-meansyini
	gen syendd=syend-meansyend
	drop mean*
	sort country period
	tsset country period
	xi: xtreg yend sk sh ndg yini i.period, fe
	matrix coefs=e(b)
	svmat coefs
	rename coefs1 skfe 	
	rename coefs2 shfe
	rename coefs3 ndgfe
	rename coefs4 yinife
	rename coefs5 consfe
	drop coefs*
	drop _Iper*
	xi: xtreg syend ssk ssh sndg syini i.period, fe
	matrix coefs=e(b)
	svmat coefs
	rename coefs1 skfeshock
	rename coefs2 shfeshock
	rename coefs3 ndgfeshock
	rename coefs4 yinifeshock
	rename coefs5 consfeshock
	drop coefs*
	drop _Iper*
	xi: xtreg syend ssk ssh sndg syini i.period, be
	matrix coefs=e(b)
	svmat coefs
	rename coefs1 skbeshock
	rename coefs2 shbeshock
	rename coefs3 ndgbeshock
	rename coefs4 yinibeshock
	rename coefs5 consbeshock
	drop coefs*
	drop _Iper*
	xi: xtreg syend ssk ssh sndg syini i.period, re
	matrix coefs=e(b)
	svmat coefs
	rename coefs1 skreshock
	rename coefs2 shreshock
	rename coefs3 ndgreshock
	rename coefs4 yinireshock
	rename coefs5 consreshock
	drop coefs*
	drop _Iper*
	xtabond syendd, pre(sskd sshd sndgd)
	matrix coefs=e(b)
	svmat coefs
	rename coefs1 yiniabondshock
	rename coefs2 skabondshock
	rename coefs3 shabondshock
	rename coefs4 ndgabondshock
	rename coefs5 consabondshock
	xtabond2 syendd sskd sshd sndgd L.syendd, gmmstyle(sskd sshd sndgd L.syendd)
	matrix coefs=e(b)
	svmat coefs
	rename coefs1 skbbondshock
	rename coefs2 shbbondshock
	rename coefs3 ndgbbondshock
	rename coefs4 yinibbondshock
	rename coefs5 consbbondshock
	#delimit ;
	keep skfe shfe ndgfe yinife consfe skfeshock shfeshock ndgfeshock yinifeshock consfeshock
	skreshock shreshock ndgreshock yinireshock consreshock
	skabondshock shabondshock ndgabondshock yiniabondshock consabondshock
	skbbondshock shbbondshock ndgbbondshock yinibbondshock consbbondshock
	skbeshock shbeshock ndgbeshock yinibeshock consbeshock;
	#delimit cr
	keep if _n==1
	save montecarloresults, replace
	noisily disp $t
	global t=$t+1
}

quietly while $t<=$drawnum{
	do MCgenerate
	keep  country period yini sk sh ndg yend sndg ssh ssk syini syend
	* Define variables and tsset command for Arellano-Bond
	sort period country
	by period: gen meanssk=sum(ssk)
	by period: gen meanssh=sum(ssh)
	by period: gen meansndg=sum(sndg)
	by period: gen meansyini=sum(syini)
	by period: gen meansyend=sum(syend)
	by period: replace meanssk=meanssk[_N]/$N
	by period: replace meanssh=meanssh[_N]/$N
	by period: replace meansndg=meansndg[_N]/$N
	by period: replace meansyini=meansyini[_N]/$N
	by period: replace meansyend=meansyend[_N]/$N
	gen sskd=ssk-meanssk
	gen sshd=ssh-meanssh
	gen sndgd=sndg-meansndg
	gen syinid=syini-meansyini
	gen syendd=syend-meansyend
	drop mean*
	sort country period
	tsset country period
	xi: xtreg yend sk sh ndg yini i.period, fe
	matrix coefs=e(b)
	svmat coefs
	rename coefs1 skfe 	
	rename coefs2 shfe
	rename coefs3 ndgfe
	rename coefs4 yinife
	rename coefs5 consfe
	drop coefs*
	drop _Iper*
	xi: xtreg syend ssk ssh sndg syini i.period, fe
	matrix coefs=e(b)
	svmat coefs
	rename coefs1 skfeshock
	rename coefs2 shfeshock
	rename coefs3 ndgfeshock
	rename coefs4 yinifeshock
	rename coefs5 consfeshock
	drop coefs*
	drop _Iper*
	xi: xtreg syend ssk ssh sndg syini i.period, be
	matrix coefs=e(b)
	svmat coefs
	rename coefs1 skbeshock
	rename coefs2 shbeshock
	rename coefs3 ndgbeshock
	rename coefs4 yinibeshock
	rename coefs5 consbeshock
	drop coefs*
	drop _Iper*
	xi: xtreg syend ssk ssh sndg syini i.period, re
	matrix coefs=e(b)
	svmat coefs
	rename coefs1 skreshock
	rename coefs2 shreshock
	rename coefs3 ndgreshock
	rename coefs4 yinireshock
	rename coefs5 consreshock
	drop coefs*
	drop _Iper*
	xtabond syendd, pre(sskd sshd sndgd)
	matrix coefs=e(b)
	svmat coefs
	rename coefs1 yiniabondshock
	rename coefs2 skabondshock
	rename coefs3 shabondshock
	rename coefs4 ndgabondshock
	rename coefs5 consabondshock
	xtabond2 syendd sskd sshd sndgd L.syendd, gmmstyle(sskd sshd sndgd L.syendd)
	matrix coefs=e(b)
	svmat coefs
	rename coefs1 skbbondshock
	rename coefs2 shbbondshock
	rename coefs3 ndgbbondshock
	rename coefs4 yinibbondshock
	rename coefs5 consbbondshock
	#delimit ;
	keep skfe shfe ndgfe yinife consfe skfeshock shfeshock ndgfeshock yinifeshock consfeshock
	skreshock shreshock ndgreshock yinireshock consreshock
	skabondshock shabondshock ndgabondshock yiniabondshock consabondshock
	skbbondshock shbbondshock ndgbbondshock yinibbondshock consbbondshock
	skbeshock shbeshock ndgbeshock yinibeshock consbeshock;
	#delimit cr
	keep if _n==1
	append using montecarloresults
	save montecarloresults, replace
	noisily disp $t
	global t=$t+1
}


#delimit ;
order skfe shfe ndgfe yinife consfe skfeshock shfeshock ndgfeshock yinifeshock consfeshock
skbeshock shbeshock ndgbeshock yinibeshock consbeshock
skreshock shreshock ndgreshock yinireshock consreshock
skabondshock shabondshock ndgabondshock yiniabondshock consabondshock
skbbondshock shbbondshock ndgbbondshock yinibbondshock consbbondshock;
#delimit cr

