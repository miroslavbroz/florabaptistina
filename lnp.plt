#!/usr/bin/env gnuplot

set palette rgbformulae 33,13,10

#set yr [-10000:]
set zeroaxis
set key bottom

p "lnp.tmp" u 1:(column($#)):2 lc palette z

pa -1

set term png small
set out "lnp.png"
rep

q

p "lnp.tmp" u 1:6:2 lc palette z

