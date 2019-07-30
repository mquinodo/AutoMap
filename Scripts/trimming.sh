#!/bin/bash

# trimming terminal heterozygous variants on both sides of regions 

awk -F"\t" 'BEGIN{OFS="\t"} {if(state=="No" && $13=="Yes" && $6=="het") $13="No"; print $0} {state=$13} {zyg=$6}' $1 \
| tac | awk -F"\t" 'BEGIN{OFS="\t"} {if(state=="No" && $13=="Yes" && $6=="het") $13="No"; print $0} {state=$13} {zyg=$6}' | tac > $2

