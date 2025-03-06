#! /bin/sh

module load r/4.4.1
Rscript 03_Stats_inversions.r
module unload r/4.4.1
