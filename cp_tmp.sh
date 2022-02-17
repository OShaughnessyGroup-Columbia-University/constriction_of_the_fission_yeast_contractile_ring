#!/bin/bash
FIRST=1;
LAST=750;
for i in $(eval echo "{$FIRST..$LAST}"); 
do 
    cp wtfc4_1_10min.mat gm3_$i\_0sec.mat
    cp wtfc4_1_10min.mat gmdrop_$i\_0sec.mat
done
