#!/usr/bin/env gnuplot

set xl "sin {/Helvetica-Oblique i} [1]"
set yl "{/Helvetica-Oblique p}{/=2{ }}(sin {/Helvetica-Oblique i}{/=4{ }}) [1]"

set palette rgbformulae 33,13,10

set yr [:]
set zeroaxis
set key above samplen 1.5 width -2

set arrow from 0.0825,graph 0 rto 0,graph 1 nohead lt 0
set arrow from 0.1175,graph 0 rto 0,graph 1 nohead lt 0

p \
  "<gawk '{ print $1,$2; }' hist3.out  | ./norm.awk" u 1:2 w histeps lw 3 lc 'gray' t "observed",\
  "<gawk '(NF>=4){ print $1,exp($4); }' ../lnlikei.tmp | sort -g" u 1:2 w l lw 4 lc 'black' t "F. plus B.",\
  "<gawk '(NF>=2){ print $1,exp($2); }' ../lnlikei.tmp | sort -g" u 1:2 w l lw 3 dt 1 lc '#ff6600' t "Flora",\
  "<gawk '(NF>=3){ print $1,exp($3); }' ../lnlikei.tmp | sort -g" u 1:2 w l lw 2 dt (3,1.5) lc '#6600ff' t "Baptistina",\

pa -1

set term post eps enh color solid "Helvetica,10" size 6.5cm,6cm
set out "hist3.eps"
rep

system("epstopdf hist3.eps")

q

  "<gawk '{ print $5,exp($6); }' lnlike1.tmp | sort -g | ./norm.awk" u 1:2 w lp,\
  "<gawk '{ print $5,exp($6); }' lnlike2.tmp | sort -g | ./norm.awk" u 1:2 w lp,\
