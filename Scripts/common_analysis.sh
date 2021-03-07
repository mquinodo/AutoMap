#!/bin/bash

usage() { echo "Usage: $0 [--res <string>] [--name <string>] [--ids <string>] [--panel <string>] [--panelname <string>]" 1>&2; exit 1; }


while getopts ":-:" o; do
    case "${o}" in
    	-)  
        	case $OPTARG in
                res)
                    res="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                name)
                    out="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                ids)
                    ids="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                panel)
                    panel="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                panelname)
                    panelname="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                
				*)
            		usage
            	;;
         esac ;;        
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))

if [ -z "${panel}" ]; then
    panel=""
fi
if [ -z "${panelname}" ]; then
    panelname=""
fi

IFS=',' read -ra pats <<< "$ids"

here="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

mkdir -p $res/$out 2> $here.log

#### processing ids
count=0
for pat in "${pats[@]}"
do
    if [ ! -f $res/$pat/$pat.HomRegions.tsv ]; then
        echo "$res/$pat/$pat.HomRegions.tsv is not found. -> Exit"
        exit 1
    else
        echo "# $pat used."
    fi
    grep -v "#" $res/$pat/$pat.HomRegions.tsv > $res/$pat/$pat.HomRegions.noH.tsv
    if [ "$count" -gt "0" ]; then
        bedtools intersect -a $res/$pat/$pat.HomRegions.noH.tsv -b $res/$out/$out.temp.$(($count-1)).tsv | cut -f1-3 > $res/$out/$out.temp.$count.tsv
    else
        cat $res/$pat/$pat.HomRegions.noH.tsv > $res/$out/$out.temp.$count.tsv
    fi
    count=$(($count+1))
    rm $res/$pat/$pat.HomRegions.noH.tsv
done
cat $res/$out/$out.temp.$(($count-1)).tsv > $res/$out/$out.temp.ids.tsv
cat $res/$out/$out.temp.ids.tsv > $res/$out/$out.temp.0.tsv
count=1

echo -e "#Chr\tBegin\tEnd" > $res/$out/$out.HomRegions.temp1.tsv
echo -e "Size(Mb)" > $res/$out/$out.HomRegions.temp2.tsv
cat $res/$out/$out.temp.$(($count-1)).tsv >> $res/$out/$out.HomRegions.temp1.tsv
cat $res/$out/$out.temp.$(($count-1)).tsv | awk -F"\t" '{printf("%.2f\n",($3-$2)/1000000)}' >> $res/$out/$out.HomRegions.temp2.tsv
paste -d"\t" $res/$out/$out.HomRegions.temp1.tsv $res/$out/$out.HomRegions.temp2.tsv > $res/$out/$out.HomRegions.tsv
rm $res/$out/$out.temp*.tsv $res/$out/$out.HomRegions.temp*.tsv

# produce output text file
file=$res/$out/$out.HomRegions.tsv
tot="$(grep -v "#" $file | grep -v "chrX" | awk '{s+=($3-$2)/1000000} END {print s}')"
LC_NUMERIC="en_US.UTF-8" tot2="$(printf "%.2f" $tot)"
if [ "$tot2" == "" ]; then
tot2=0
fi
echo "## INFO: $tot2 Mb are in Homozygous Regions (autosomal chromosomes)" >> $file;

# annotate outpu file with panel
if [ -s $panel ]; then
    perl $here/gene_names.pl $res/$out/$out.HomRegions.tsv $res/$out/$out.HomRegions.$panelname.tsv $panel $panelname 2> $here.log
fi 

# produce output pdf file
size=$(more $file | grep INFO | awk -F " " '{print $3}' )
outputR=$res/$out/$out.HomRegions.pdf
Rscript $here/make_graph_common.R $out $file $outputR $size 2> $here.log

echo "## There is $tot2 Mb are in common Homozygous Regions (autosomal chromosomes) after analysis. ##"

