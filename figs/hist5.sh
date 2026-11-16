#!/bin/sh

gawk '!/\?/{ print $5; }' < elmts.txt > hist5.in

hist6 100 0 -1 hist5.in > hist5.out

./hist5.plt

