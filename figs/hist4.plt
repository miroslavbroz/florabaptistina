#!/usr/bin/env gnuplot

set xl "{/Helvetica-Oblique H} [1]"
set yl "{/Helvetica-Oblique p}{/=2{ }}({/Helvetica-Oblique H}{/=2{ }}) [1]"

set palette rgbformulae 33,13,10

set yr [:1.5]
#set logscale y
set zeroaxis
set key above samplen 1.5 width -2

set arrow from 18.5,graph 0 rto 0,graph 1 nohead lt 0

p \
  "<gawk '{ print $1,$2; }' hist4.out  | ./norm.awk" u 1:2 w histeps lw 3 lc 'gray' t "observed",\
  "<gawk '(NF>=4){ print $1,exp($4); }' ../lnlikeh.tmp | sort -g" u 1:2 w l lw 4 lc 'black' t "F. plus B.",\
  "<gawk '(NF>=2){ print $1,exp($2); }' ../lnlikeh.tmp | sort -g" u 1:2 w l lw 3 dt 1 lc '#ff6600' t "Flora",\
  "<gawk '(NF>=3){ print $1,exp($3); }' ../lnlikeh.tmp | sort -g" u 1:2 w l lw 2 dt (3,1.5) lc '#6600ff' t "Baptistina",\

pa -1

set term post eps enh color solid "Helvetica,10" size 6.5cm,6cm
set out "hist4.eps"
rep

system("epstopdf hist4.eps")


q

  "<gawk '($1<18.5){ print $1,$2; }' hist4.out  | ./norm.awk" u 1:2 w histeps t "observed"

