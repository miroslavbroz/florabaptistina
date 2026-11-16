#!/usr/bin/env gnuplot

set xl "{/Helvetica-Oblique a} [au]"
set yl "{/Helvetica-Oblique p}{/=2{ }}({/Helvetica-Oblique a}{/=2{ }}|{/Helvetica-Oblique H}{/=4{ }}) [1]"

set palette rgbformulae 33,13,10

set yr [:]
set yr [-0.1:15]
#set logscale y
set zeroaxis
set key above samplen 1.5 width -2

p \
  "<gawk '{ print $1,$2; }' hist1.out  | ./norm.awk" u 1:2 w histeps lw 3 lc 'gray' t "observed",\
  "<gawk '(NF>=4){ print $1,exp($4); }' ../lnlikea.tmp | sort -g" u 1:2 w p pt 7 ps 0.2 lw 1 lc 'black' t "F. plus B.",\
  "<gawk '(NF>=2){ print $1,exp($2); }' ../lnlikea.tmp | sort -g" u 1:2 w p pt 7 ps 0.2 lw 1 dt 1 lc '#ff6600' t "Flora",\
  "<gawk '(NF>=3){ print $1,exp($3); }' ../lnlikea.tmp | sort -g" u 1:2 w p pt 7 ps 0.2 lw 1 dt (3,1.5) lc '#6600ff' t "Baptistina",\

pa -1

set term post eps enh color solid "Helvetica,10" size 6.5cm,6cm
set out "hist1.eps"
rep

system("epstopdf hist1.eps")


q

  "<gawk '(NF>=4){ print $1,exp($4); }' ../lnlikea.tmp | sort -g" u 1:2 w p pt 6 t "3",\
  "<gawk '(NF>=5){ print $1,exp($5); }' ../lnlikea.tmp | sort -g" u 1:2 w p t "4",\
  "<gawk '(NF>=6){ print $1,exp($6); }' ../lnlikea.tmp | sort -g" u 1:2 w p t "all",\
  "<paste ../lnlikea.tmp ../lnlikeh.tmp | sort -g" u 1:(exp($6+$12)) w d,\

