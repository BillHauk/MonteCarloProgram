*This routine takes the modified raw data generated in datasetgenerate and
*calculates the mean, standard deviations, and correlations of those variables
*Monte Carlo draws are then generated using these moments

clear
do covariance

#delimit ;
drawnorm yini3 sk3 sh3 ndg3 yini4 sk4 sh4 ndg4 yini5 sk5 sh5 ndg5 yini6 sk6 sh6
ndg6 yini7 sk7 sh7 ndg7 yini8 sk8 sh8 ndg8 yini9 sk9 sh9 ndg9 yini10 sk10 sh10
ndg10 fes dq , n($N) means(means) cov(covariance);
#delimit cr

replace dq =1 if dq <1
replace dq =4 if dq >4

generate country=_n
order country

save generateddata, replace

* These commands add white noise to the observations
use baselinedatawide
global i 3
while $i==3{
	sum resids$i
	matrix N=(r(mean))
	matrix T=(r(sd))
	matrix M=N
	matrix S=T
	global i=$i+1
}
while $i<=10{
	sum resids$i
	matrix N=(r(mean))
	matrix T=(r(sd))
	matrix M=(M, N)
	matrix S=(S, T)
	global i=$i+1
}
matrix accum X=resids3 resids4 resids5 resids6 resids7 resids8 resids9 resids10, nocons dev
matrix R=corr(X)
drawnorm e3 e4 e5 e6 e7 e8 e9 e10, n($N) corr(R) sds(S)
generate country=_n
sort country
keep country e*
save errors, replace
clear

use generateddata
sort country
merge country using errors
drop _merge
save generateddata, replace

*These commands make the residual terms correlated with the other RHS variables

*Note: The version of this program used to generate Table 5 in the published
*version of the paper had a minor typo that caused the residual correlations to
*not vary as described in Table 5.  The corrected version follows below.  This error
*did not substantively alter the results found in Table 5 or described in 
*Section 3.2.3 of the published paper.

clear
use generateddata

keep country sk* sh* ndg* e*

global i 3

while $i<=10 {
	sum e$i
	global errsum=r(mean)
	global errsd=r(sd)
	sum sk$i
	global sksd=r(sd)
	sum sh$i
	global shsd=r(sd)
	sum ndg$i
	global ndgsd=r(sd)
	matrix accum X=(sk$i sh$i ndg$i), nocons dev
	matrix cov=X/($N-1)
	matrix covinv=inv(cov)
	matrix corrs=($skerrcorr*$sksd*$errsd, $sherrcorr*$shsd*$errsd, $ndgerrcorr*$ndgsd*$errsd)
	matrix rhos=covinv*corrs'
	gen err$i=rhos[1,1]*sk$i+rhos[2,1]*sh$i+rhos[3,1]*ndg$i
	sum err$i
	global epsmean=$errsum-r(mean)
	global epssd=(($errsd)^2-(r(sd))^2)^.5
	drawnorm eps$i, means($epsmean) sds($epssd)
	gen newe$i=err$i+eps$i
	global i=$i+1
}

keep country newe*
sort country
save errors, replace

use generateddata
sort country
merge country using errors

global i 3
while $i<=10 {
	replace e$i=newe$i
	global i=$i+1
}

drop newe* _merge
save generateddata, replace

*Now dependent variables generated using the "true" values of alpha and beta and the independent
*variables previously drawn

#delimit ;

* yend are end of period values of generated income, yini is the same for beginning of period;

generate yend3=$multiple*($alpha/(1-$alpha-$beta))*sk3+$multiple*($beta/(1-$alpha-$beta))*sh3
-$multiple*(($alpha+$beta)/(1-$alpha-$beta))*ndg3+(1-$multiple)*yini3+fes+e3+0.02*5*(4-(1-$multiple)*3);

generate yend4=$multiple*($alpha/(1-$alpha-$beta))*sk4+$multiple*($beta/(1-$alpha-$beta))*sh4
-$multiple*(($alpha+$beta)/(1-$alpha-$beta))*ndg4+(1-$multiple)*yend3+fes+e4+0.02*5*(5-(1-$multiple)*4);

generate yend5=$multiple*($alpha/(1-$alpha-$beta))*sk5+$multiple*($beta/(1-$alpha-$beta))*sh5
-$multiple*(($alpha+$beta)/(1-$alpha-$beta))*ndg5+(1-$multiple)*yend4+fes+e5+0.02*5*(6-(1-$multiple)*5);

generate yend6=$multiple*($alpha/(1-$alpha-$beta))*sk6+$multiple*($beta/(1-$alpha-$beta))*sh6
-$multiple*(($alpha+$beta)/(1-$alpha-$beta))*ndg6+(1-$multiple)*yend5+fes+e6+0.02*5*(7-(1-$multiple)*6);

generate yend7=$multiple*($alpha/(1-$alpha-$beta))*sk7+$multiple*($beta/(1-$alpha-$beta))*sh7
-$multiple*(($alpha+$beta)/(1-$alpha-$beta))*ndg7+(1-$multiple)*yend6+fes+e7+0.02*5*(8-(1-$multiple)*7);

generate yend8=$multiple*($alpha/(1-$alpha-$beta))*sk8+$multiple*($beta/(1-$alpha-$beta))*sh8
-$multiple*(($alpha+$beta)/(1-$alpha-$beta))*ndg8+(1-$multiple)*yend7+fes+e8+0.02*5*(9-(1-$multiple)*8);

generate yend9=$multiple*($alpha/(1-$alpha-$beta))*sk9+$multiple*($beta/(1-$alpha-$beta))*sh9
-$multiple*(($alpha+$beta)/(1-$alpha-$beta))*ndg9+(1-$multiple)*yend8+fes+e9+0.02*5*(10-(1-$multiple)*9);

generate yend10=$multiple*($alpha/(1-$alpha-$beta))*sk10+$multiple*($beta/(1-$alpha-$beta))*sh10
-$multiple*(($alpha+$beta)/(1-$alpha-$beta))*ndg10+(1-$multiple)*yend9+fes+e10+0.02*5*(11-(1-$multiple)*10);

#delimit cr

drop yini4 yini5 yini6 yini7 yini8 yini9 yini10

*generate yini3=yini3
generate yini4=yend3
generate yini5=yend4
generate yini6=yend5
generate yini7=yend6
generate yini8=yend7
generate yini9=yend8
generate yini10=yend9

keep country yend* sk* sh* ndg* yini* dq fes

iis country
save generateddata, replace

global Ntot=$N*8

global i 3
while $i==3{
	sum yini$i
	matrix U=(((r(sd))^2*$shockyini)^.5)
	matrix N=(-((r(sd))^2*$shockyini)/2)
	sum sk$i
	matrix U=(U, ((r(sd))^2*$shocksk)^.5)
	matrix N=(N, ((r(sd))^2*$shocksk)/2)
	sum sh$i
	matrix U=(U, ((r(sd))^2*$shocksh)^.5)
	matrix N=(N, ((r(sd))^2*$shocksh)/2)
	sum ndg$i
	matrix U=(U, ((r(sd))^2*$shockndg)^.5)
	matrix N=(N, ((r(sd))^2*$shockndg)/2)
	matrix V=U
	matrix T=N
	global i=$i+1
}

while $i<=10{
	sum yini$i
	matrix U=(((r(sd))^2*$shockyini)^.5)
	matrix N=(-((r(sd))^2*$shockyini)/2)
	sum sk$i
	matrix U=(U, ((r(sd))^2*$shocksk)^.5)
	matrix N=(N, ((r(sd))^2*$shocksk)/2)
	sum sh$i
	matrix U=(U, ((r(sd))^2*$shocksh)^.5)
	matrix N=(N, ((r(sd))^2*$shocksh)/2)
	sum ndg$i
	matrix U=(U, ((r(sd))^2*$shockndg)^.5)
	matrix N=(N, ((r(sd))^2*$shockndg)/2)
	matrix V=(V, U)
	matrix T=(T, N)
	global i=$i+1
}
sum yini10
matrix V=(V, ((r(sd))^2*$shockyini)^.5)
matrix T=(T, ((r(sd))^2*$shockyini)/2)
matrix V=(V, 0)
matrix T=(T, 0)

#delimit;
matrix W=(1,$errrho,$errrho^2,$errrho^3,$errrho^4,$errrho^5,$errrho^6,$errrho^7\
$errrho,1,$errrho,$errrho^2,$errrho^3,$errrho^4,$errrho^5,$errrho^6\
$errrho^2,$errrho,1,$errrho,$errrho^2,$errrho^3,$errrho^4,$errrho^5\
$errrho^3,$errrho^2,$errrho,1,$errrho,$errrho^2,$errrho^3,$errrho^4\
$errrho^4,$errrho^3,$errrho^2,$errrho,1,$errrho,$errrho^2,$errrho^3\
$errrho^5,$errrho^4,$errrho^3,$errrho^2,$errrho,1,$errrho,$errrho^2\
$errrho^6,$errrho^5,$errrho^4,$errrho^3,$errrho^2,$errrho,1,$errrho\
$errrho^7,$errrho^6,$errrho^5,$errrho^4,$errrho^3,$errrho^2,$errrho,1);
#delimit cr

matrix Y=(1,0,0,0\0,1,0,0\0,0,1,0\0,0,0,1)
matrix Z=W#Y
matrix A=($errrho^8,0,0,0,$errrho^7,0,0,0,$errrho^6,0,0,0,$errrho^5,0,0,0,$errrho^4,0,0,0,$errrho^3,0,0,0,$errrho^2,0,0,0,$errrho,0,0,0)
matrix C=(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
matrix A=(A\C)
matrix B=(1,0\0,1)
matrix Z=(Z,A'\A,B)

use baselinedatawide
#delimit ;
drawnorm shockyini3 shocksk3 shocksh3 shockndg3 shockyini4 shocksk4
shocksh4 shockndg4 shockyini5 shocksk5 shocksh5 shockndg5 shockyini6 
shocksk6 shocksh6 shockndg6 shockyini7 shocksk7 shocksh7 shockndg7 shockyini8 shocksk8
shocksh8 shockndg8 shockyini9 shocksk9 shocksh9 shockndg9 shockyini10 shocksk10
shocksh10 shockndg10 shockyend10 shockfes, n($N) sds(V) corr(Z);
#delimit cr
generate country=_n
drop shockfes
reshape long shockndg shocksh shocksk shockyini shockyend shockfes, i(country) j(period)
sort country period
save shocks, replace

use generateddata
reshape long yend sk sh ndg yini, i(country) j(period)
sort country period
merge country period using shocks
keep country period sk sh ndg yend yini shockyini shocksk shocksh shockndg shockyend shockfes dq
sort country period
save generateddata, replace

sum dq
generate dqmultiple=(dq/r(mean))^.5
replace dqmultiple=1 if $dataqual==0

generate sndg=ndg+shockndg*dqmultiple
generate ssh=sh+shocksh*dqmultiple
generate ssk=sk+shocksk*dqmultiple
generate syini=yini+shockyini*dqmultiple
by country: generate syend=yend+shockyini[_n+1]*dqmultiple
by country: replace syend=yend+shockyini[_n-7]*dqmultiple if period==10
