#!/bin/sh

gawk '{ print $2; }' < elmts.txt > hist2.in

hist6 100 0 -1 hist2.in > hist2.out

./hist2.plt

