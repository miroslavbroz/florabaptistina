#!/bin/sh

gawk '!/^#/ && (NF>0){ i++; if ($57>0.12){ print $36,$37,$38,$4,$57,0.0,i; } }' < 8_Flora_halo.list_110_WO_INTERLOPERS_298 > elmts.txt

