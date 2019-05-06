#!/bin/bash

usage() { echo "Usage: $0 [--res <string>] [--name <string>] [--affected <string>] [--healthy <string>] [--panel <string>] [--panelname <string>]" 1>&2; exit 1; }


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
                affected)
                    affected="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                healthy)
                    healthy="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
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

if [ -z "${res}" ]; then
    echo "## ERROR: You need to provide the location of previous results used by main script through --res option"
    exit 1
fi
if [ -z "${out}" ]; then
    echo "## ERROR: You need the name of the analysis through --name option"
    exit 1
fi
if [ -z "${affected}" ]; then
    echo "## ERROR: You need to provide at least one affected individual through --affected option"
    exit 1
fi
if [ -z "${healthy}" ]; then
    healthy=""
fi
if [ -z "${panel}" ]; then
    panel=""
fi
if [ -z "${panelname}" ]; then
    panelname=""
fi

IFS=',' read -ra pats <<< "$affected"
IFS=',' read -ra cons <<< "$healthy"

here="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

mkdir -p $res/$out 2> $here.log

#### processing affected
count=0
for pat in "${pats[@]}"
do
    if [ ! -f $res/$pat/$pat.HomRegions.tsv ]; then
        echo "$res/$pat/$pat.HomRegions.tsv is not found. -> Exit"
        exit 1
    else
        echo "# $pat used as affected."
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
cat $res/$out/$out.temp.$(($count-1)).tsv > $res/$out/$out.temp.affected.tsv
cat $res/$out/$out.temp.affected.tsv > $res/$out/$out.temp.0.tsv

#### processing controls
count=1
for con in "${cons[@]}"
do
    if [ ! -f $res/$con/$con.HomRegions.tsv ]; then
        echo "$res/$con/$con.HomRegions.tsv is not found. -> Exit"
        exit 1
    else
        echo "# $con used as control."
    fi
    grep -v "#" $res/$con/$con.HomRegions.tsv > $res/$con/$con.HomRegions.noH.tsv

    bedtools subtract -a $res/$out/$out.temp.$(($count-1)).tsv  -b $res/$con/$con.HomRegions.noH.tsv | cut -f1-3 > $res/$out/$out.temp.$count.tsv

    count=$(($count+1))
    rm $res/$con/$con.HomRegions.noH.tsv
done
echo -e "#Chr\tBegin\tEnd" > $res/$out/$out.HomRegions.tsv
cat $res/$out/$out.temp.$(($count-1)).tsv >> $res/$out/$out.HomRegions.tsv
rm $res/$out/$out.temp*.tsv

# produce output text file
file=$res/$out/$out.HomRegions.tsv
tot="$(grep -v "#" $file | grep -v "chrX" | awk '{s+=($3-$2)/1000000} END {print s}')"
tot2="$(printf "%.2f" $tot)"
echo "## INFO: $tot2 Mb are in Homozygous Regions (autosomal chromosomes)" >> $file;

# annotate outpu file with panel
if [ -s $panel ]; then
    perl $here/gene_names.pl $res/$out/$out.HomRegions.tsv $res/$out/$out.HomRegions.$panelname.tsv $panel $panelname 2> $here.log
fi 

# produce output pdf file
size=$(more $file | grep INFO | awk -F " " '{print $3}' )
outputR=$res/$out/$out.HomRegions.pdf
Rscript $here/script_to_make_graph_HomRegions_simple.R $out $file $outputR $size 2> $here.log

echo "## There is $tot2 Mb are in Homozygous Regions (autosomal chromosomes) after analysis. ##"

