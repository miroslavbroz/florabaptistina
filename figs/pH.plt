#!/usr/bin/env gnuplot

set xl "p_V [1]"
set yl "H [mag]"
set zl "p(a|H) p(H) [1]"

tmp=2
#set zr [-tmp:tmp]
set xyplane 0.0
set view 0,0

sp \
   "<paste ../lnlikep.tmp ../lnlikeh.tmp" u 1:5:(0.0) w d lc 'black',\

pa -1

set term png small
set out "pH.png"
rep

q

