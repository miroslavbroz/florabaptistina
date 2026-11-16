#!/usr/bin/env gnuplot

H(a,a_c,C)= log10(abs(a-a_c)/C)/ 0.2

set samples 1000
set xr [2.18:2.34]
set yr [6:20]

set table "aH__.tmp"
p H(x,2.20,2.1e-4), H(x,2.26,0.2e-4)
unset table

set xl "a [au]"
set yl "H [mag]"
set cbl "p(a|H) [1]"

set zr [:]
set cbr [:10]
set palette rgbformulae 33,13,10
set pm3d map
set key above

sp \
   "<paste ../lnlikea.tmp ../lnlikeh.tmp" u 1:5:(exp($2)) w d lw 2 lc 'gray' t "1",\
   "<paste ../lnlikea.tmp ../lnlikeh.tmp" u 1:5:(exp($3)) w d lw 2 lc palette z t "2",\
   "aH__.tmp" u 1:2:(5.0) w l lw 3 dt (3,1.5) lc 'black' t "H(a)"

pa -1

set term post eps enh color solid "Helvetica,10" size 6.5cm,6cm
set out "aH__.eps"
rep

system("epstopdf aH__.eps")


q

  H(x,a_c3,C3)     not sprintf("C = %.2f * 1e-4", C3/1.e-4) lc 'orange' lw 3 dt 2,\

   "<paste ../lnlikea.tmp ../lnlikeh.tmp" u 1:5:(exp($4)) w d t "3",\
   "<paste ../lnlikea.tmp ../lnlikeh.tmp" u 1:5:(0.0) w d lc 'gray' t "observed",\


