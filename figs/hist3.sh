#!/bin/sh

gawk '{ print $3; }' < elmts.txt > hist3.in

hist6 100 0 -1 hist3.in > hist3.out

./hist3.plt

