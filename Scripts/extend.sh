#!/bin/bash

# extending homozygous regions to closest SNP

maxsize=$3
awk -v maxsize="$maxsize" -F"\t" 'function abs(v) {return v < 0 ? -v : v} BEGIN{OFS="\t"} \
{if(state=="Yes" && $13=="No" && extend=="No") {$13="Yes"; extend="Yes"; if(abs(pos-$2)>maxsize*1000000) {$2=pos+maxsize*1000000}; $6="ext"} else {extend="No"} print $0} {state=$13; pos=$2}' $1 \
| tac | awk -v maxsize="$maxsize" -F"\t" 'function abs(v) {return v < 0 ? -v : v} BEGIN{OFS="\t"} \
{if(state=="Yes" && $13=="No" && extend=="No") {$13="Yes"; extend="Yes"; if(abs(pos-$2)>maxsize*1000000) {$2=pos-maxsize*1000000}; $6="ext"} else {extend="No"} print $0} {state=$13; pos=$2}' | tac > $2



