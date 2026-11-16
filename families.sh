#!/bin/sh

#J=2  # j_nonrandom <- JUST FOR TESTING PURPOSES!!!
J=3  # k_random

for I in 0 1 ; do

  gawk -vI=$I -vJ=$J '
(ARGIND==1) && !/^#/ && ($J==I){
  s[$1+0]=1;
}
(ARGIND==1) && !/^#/ && ($1<=2){
  if (I==0) { s[1]=1; s[2]=0; }  # (8) Flora
  if (I==1) { s[1]=0; s[2]=1; }  # (298) Baptistina
}
(ARGIND>1) && !/^#/ && (NF>0){
  k++;
  if (s[k]==1){ print; }
}
' families.tmp 8_Flora_halo.list_110_WO_INTERLOPERS_298 > family$I.dat

done

