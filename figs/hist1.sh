#!/bin/sh

gawk '{ print $1; }' < ../lnlikea.tmp > hist1.in

hist6 100 0 -1 hist1.in > hist1.out

./hist1.plt

