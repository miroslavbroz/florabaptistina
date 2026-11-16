#!/bin/sh

gawk '{ print $4; }' < elmts.txt | sort -g > hist4.in

hist6 100 0 -1 hist4.in > hist4.out

./hist4.plt

