*These loops get the means and standard deviations from the dataset generated
*in datasetgenerate

use baselinedatawide.dta
global i 3
while $i==3{
	sum lini$i
	matrix N=(r(mean))
	sum lki$i
	matrix N=(N, r(mean))
	sum lenroll$i
	matrix N=(N, r(mean))
	sum lgpop$i
	matrix N=(N, r(mean))
	matrix M=N
	global i=$i+1
}

while $i<=10{
	sum lini$i
	matrix N=(r(mean))
	sum lki$i
	matrix N=(N, r(mean))
	sum lenroll$i
	matrix N=(N, r(mean))
	sum lgpop$i
	matrix N=(N, r(mean))
	matrix M=(M, N)
	global i=$i+1
}
sum fixedeffect
matrix M=(M, r(mean))
sum dqual 
matrix M=(M, r(mean))
matrix A=($yinifecorr,0,0,0\0,$skfecorr,0,0\0,0,$shfecorr,0\0,0,0,$ndgfecorr)
matrix B=(1,0,0,0,0,0,0,0\0,1,0,0,0,0,0,0\0,0,1,0,0,0,0,0\0,0,0,1,0,0,0,0\0,0,0,0,1,0,0,0\0,0,0,0,0,1,0,0\0,0,0,0,0,0,1,0\0,0,0,0,0,0,0,1)
matrix C=B#A


* And then puts them into matricies to be used as parameters for a drawnorm command

#delimit ;
matrix accum X=lini3 lki3 lenroll3 lgpop3 lini4 lki4 lenroll4 lgpop4 lini5 lki5
lenroll5 lgpop5 lini6 lki6 lenroll6 lgpop6 lini7 lki7 lenroll7 lgpop7 lini8 lki8
lenroll8 lgpop8 lini9 lki9 lenroll9 lgpop9 lini10 lki10 lenroll10 lgpop10 fixedeffect dqual, nocons dev;
#delimit cr
matrix R=X/($N-1)
matrix D=C*R[1..32,33]
matrix R[33,1]=D'
matrix R[1,33]=D
matrix F=$yinifecorr*R[33,34]
matrix R[33,34]=F
matrix R[34,33]=F'

global i=1
matrix R[$i,1]=(1/(1+$shockyini))^.5*R[$i,1..34]
matrix R[1,$i]=(1/(1+$shockyini))^.5*R[1..34,$i]
matrix R[$i+4,1]=(1/(1+$shockyini))^.5*R[$i+4,1..34]
matrix R[1,$i+4]=(1/(1+$shockyini))^.5*R[1..34,$i+4]
matrix R[$i+8,1]=(1/(1+$shockyini))^.5*R[$i+8,1..34]
matrix R[1,$i+8]=(1/(1+$shockyini))^.5*R[1..34,$i+8]
matrix R[$i+12,1]=(1/(1+$shockyini))^.5*R[$i+12,1..34]
matrix R[1,$i+12]=(1/(1+$shockyini))^.5*R[1..34,$i+12]
matrix R[$i+16,1]=(1/(1+$shockyini))^.5*R[$i+16,1..34]
matrix R[1,$i+16]=(1/(1+$shockyini))^.5*R[1..34,$i+16]
matrix R[$i+20,1]=(1/(1+$shockyini))^.5*R[$i+20,1..34]
matrix R[1,$i+20]=(1/(1+$shockyini))^.5*R[1..34,$i+20]
matrix R[$i+24,1]=(1/(1+$shockyini))^.5*R[$i+24,1..34]
matrix R[1,$i+24]=(1/(1+$shockyini))^.5*R[1..34,$i+24]
matrix R[$i+28,1]=(1/(1+$shockyini))^.5*R[$i+28,1..34]
matrix R[1,$i+28]=(1/(1+$shockyini))^.5*R[1..34,$i+28]

global i=$i+1
matrix R[$i,1]=(1/(1+$shocksk))^.5*R[$i,1..34]
matrix R[1,$i]=(1/(1+$shocksk))^.5*R[1..34,$i]
matrix R[$i+4,1]=(1/(1+$shocksk))^.5*R[$i+4,1..34]
matrix R[1,$i+4]=(1/(1+$shocksk))^.5*R[1..34,$i+4]
matrix R[$i+8,1]=(1/(1+$shocksk))^.5*R[$i+8,1..34]
matrix R[1,$i+8]=(1/(1+$shocksk))^.5*R[1..34,$i+8]
matrix R[$i+12,1]=(1/(1+$shocksk))^.5*R[$i+12,1..34]
matrix R[1,$i+12]=(1/(1+$shocksk))^.5*R[1..34,$i+12]
matrix R[$i+16,1]=(1/(1+$shocksk))^.5*R[$i+16,1..34]
matrix R[1,$i+16]=(1/(1+$shocksk))^.5*R[1..34,$i+16]
matrix R[$i+20,1]=(1/(1+$shocksk))^.5*R[$i+20,1..34]
matrix R[1,$i+20]=(1/(1+$shocksk))^.5*R[1..34,$i+20]
matrix R[$i+24,1]=(1/(1+$shocksk))^.5*R[$i+24,1..34]
matrix R[1,$i+24]=(1/(1+$shocksk))^.5*R[1..34,$i+24]
matrix R[$i+28,1]=(1/(1+$shocksk))^.5*R[$i+28,1..34]
matrix R[1,$i+28]=(1/(1+$shocksk))^.5*R[1..34,$i+28]

global i=$i+1
matrix R[$i,1]=(1/(1+$shocksh))^.5*R[$i,1..34]
matrix R[1,$i]=(1/(1+$shocksh))^.5*R[1..34,$i]
matrix R[$i+4,1]=(1/(1+$shocksh))^.5*R[$i+4,1..34]
matrix R[1,$i+4]=(1/(1+$shocksh))^.5*R[1..34,$i+4]
matrix R[$i+8,1]=(1/(1+$shocksh))^.5*R[$i+8,1..34]
matrix R[1,$i+8]=(1/(1+$shocksh))^.5*R[1..34,$i+8]
matrix R[$i+12,1]=(1/(1+$shocksh))^.5*R[$i+12,1..34]
matrix R[1,$i+12]=(1/(1+$shocksh))^.5*R[1..34,$i+12]
matrix R[$i+16,1]=(1/(1+$shocksh))^.5*R[$i+16,1..34]
matrix R[1,$i+16]=(1/(1+$shocksh))^.5*R[1..34,$i+16]
matrix R[$i+20,1]=(1/(1+$shocksh))^.5*R[$i+20,1..34]
matrix R[1,$i+20]=(1/(1+$shocksh))^.5*R[1..34,$i+20]
matrix R[$i+24,1]=(1/(1+$shocksh))^.5*R[$i+24,1..34]
matrix R[1,$i+24]=(1/(1+$shocksh))^.5*R[1..34,$i+24]
matrix R[$i+28,1]=(1/(1+$shocksh))^.5*R[$i+28,1..34]
matrix R[1,$i+28]=(1/(1+$shocksh))^.5*R[1..34,$i+28]

global i=$i+1
matrix R[$i,1]=(1/(1+$shockndg))^.5*R[$i,1..34]
matrix R[1,$i]=(1/(1+$shockndg))^.5*R[1..34,$i]
matrix R[$i+4,1]=(1/(1+$shockndg))^.5*R[$i+4,1..34]
matrix R[1,$i+4]=(1/(1+$shockndg))^.5*R[1..34,$i+4]
matrix R[$i+8,1]=(1/(1+$shockndg))^.5*R[$i+8,1..34]
matrix R[1,$i+8]=(1/(1+$shockndg))^.5*R[1..34,$i+8]
matrix R[$i+12,1]=(1/(1+$shockndg))^.5*R[$i+12,1..34]
matrix R[1,$i+12]=(1/(1+$shockndg))^.5*R[1..34,$i+12]
matrix R[$i+16,1]=(1/(1+$shockndg))^.5*R[$i+16,1..34]
matrix R[1,$i+16]=(1/(1+$shockndg))^.5*R[1..34,$i+16]
matrix R[$i+20,1]=(1/(1+$shockndg))^.5*R[$i+20,1..34]
matrix R[1,$i+20]=(1/(1+$shockndg))^.5*R[1..34,$i+20]
matrix R[$i+24,1]=(1/(1+$shockndg))^.5*R[$i+24,1..34]
matrix R[1,$i+24]=(1/(1+$shockndg))^.5*R[1..34,$i+24]
matrix R[$i+28,1]=(1/(1+$shockndg))^.5*R[$i+28,1..34]
matrix R[1,$i+28]=(1/(1+$shockndg))^.5*R[1..34,$i+28]

matrix means=M
matrix covariance=R


