#!/bin/bash

usage() { echo "Usage: $0 [--vcf <string>] [--genome <hg19|hg38>] [--out <string>] [--pat <string>] [--panel <string>] [--panelname <string>] [--DP <0-99>] [--binomial <0-1>] [--percaltlow <0-1.0>] [--percalthigh <0-1.0>] [--window <3-99>] [--windowthres <1-999>] [--minsize <0-99>] [--minvar <1-999>] [--minperc <0-100>] [--maxgap <0-1000Mb>] [--chrX] [--extend <0-100>]" 1>&2; exit 1; }
numbervar() { echo "Less than 10'000 variants with AD and DP available. Exit." 1>&2; exit 1; }

while getopts ":-:" o; do
    case "${o}" in
    	-)  
        	case $OPTARG in
                vcf)
                    vcf="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    nb=$(grep -v "#" $vcf | grep -P "AD|DP4" | grep GT | wc -l)
                    (($(echo "$nb >= 10000" | bc -l))) || numbervar
                    ;;
                genome)
                    genome="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    [ "$genome" == "hg19" ] || [ "$genome" == "hg38" ] || usage
                    ;;
                out)
                    out="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                pat)
                    pat="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                panel)
                    panel="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                panelname)
                    panelname="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                DP)
                    DP=$(echo "${!OPTIND}" | bc); OPTIND=$(( $OPTIND + 1 ))
                    (($(echo "$DP >= 0" | bc -l) && $(echo "$DP<=99" | bc -l))) || usage
                    ;;
                binomial)
                    binomial=$(echo "${!OPTIND}" | bc); OPTIND=$(( $OPTIND + 1 ))
                    (($(echo "$binomial >= 0" | bc -l) && $(echo "$binomial<=1" | bc -l))) || usage
                    ;;
                percaltlow)
                    percaltlow=$(echo "${!OPTIND}" | bc); OPTIND=$(( $OPTIND + 1 ))
                    (($(echo "$percaltlow >= 0" | bc -l) && $(echo "$percaltlow<=1" | bc -l))) || usage
                    ;;
                percalthigh)
                    percalthigh=$(echo "${!OPTIND}" | bc); OPTIND=$(( $OPTIND + 1 ))
                    (($(echo "$percalthigh >= 0" | bc -l) && $(echo "$percalthigh<=1" | bc -l))) || usage
                    ;;
                window)
                    window=$(echo "${!OPTIND}" | bc); OPTIND=$(( $OPTIND + 1 ))
                    (($(echo "$window >= 3" | bc -l) && $(echo "$window<=999" | bc -l))) || usage
                    ;;
                windowthres)
                    windowthres=$(echo "${!OPTIND}" | bc); OPTIND=$(( $OPTIND + 1 ))
                    (($(echo "$windowthres >= 1" | bc -l) && $(echo "$windowthres<=999" | bc -l))) || usage
                    ;;
                minsize)
                    minsize=$(echo "${!OPTIND}" | bc); OPTIND=$(( $OPTIND + 1 ))
                    (($(echo "$minsize >= 0" | bc -l) && $(echo "$minsize<=99" | bc -l))) || usage
                    ;;
                minvar)
                    minvar=$(echo "${!OPTIND}" | bc); OPTIND=$(( $OPTIND + 1 ))
                    (($(echo "$minvar >= 1" | bc -l) && $(echo "$minvar<=999" | bc -l))) || usage
                    ;;
                minperc)
                    minperc=$(echo "${!OPTIND}" | bc); OPTIND=$(( $OPTIND + 1 ))
                    (($(echo "$minperc >= 0" | bc -l) && $(echo "$minperc<=100" | bc -l))) || usage
                    ;;
                maxgap)
                    maxgap=$(echo "${!OPTIND}" | bc); OPTIND=$(( $OPTIND + 1 ))
                    (($(echo "$maxgap >= 0" | bc -l) && $(echo "$maxgap<=1000" | bc -l))) || usage
                    ;;
                chrX)
                    chrx="Yes"
                    ;;
                extend)
                    extend=$(echo "${!OPTIND}" | bc); OPTIND=$(( $OPTIND + 1 ))
                    (($(echo "$extend >= 0" | bc -l) && $(echo "$extend<=100" | bc -l))) || usage
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

echo "## Parameters used by default:"
if [ -z "${vcf}" ]; then
    echo "## ERROR: You need to provide an input vcf through --vcf option"
    exit 1
fi
if [ -z "${genome}" ]; then
    echo "## ERROR: You need to provide the genome version through --genome option (hg19 or hg38)"
    exit 1
fi
if [ -z "${out}" ]; then
    echo "## ERROR: You need to provide an output directory through --out option"
    exit 1
fi
if [ -z "${panel}" ] && [ -z "${panelname}" ]; then
    echo "No panel used, you need to provide both a panel name and panel through --panelname and --panel options"
fi
if [ -z "${panel}" ]; then
    panelname=""
fi
if [ -z "${panelname}" ]; then
    panel=""
fi
if [ -z "${DP}" ]; then
    DP=15
    echo " -> No use of --DP option, value set as default: 15"
fi
if [ -z "${binomial}" ]; then
    binomial=0.000001
    echo " -> No use of --binomial option, value set as default: 0.000001"
fi
if [ -z "${percaltlow}" ]; then
    percaltlow=0.25
    echo " -> No use of --percaltlow option, value set as default: 0.25"
fi
if [ -z "${percalthigh}" ]; then
    percalthigh=0.75
    echo " -> No use of --percalthigh option, value set as default: 0.75"
fi
if [ -z "${window}" ]; then
    window=7
    echo " -> No use of --window option, value set as default: 7"
fi
if [ -z "${windowthres}" ]; then
    windowthres=5
    echo " -> No use of --windowthres option, value set as default: 5"
fi
if [ -z "${minsize}" ]; then
    minsize=2
    echo " -> No use of --minsize option, value set as default: 2"
fi
if [ -z "${minvar}" ]; then
    minvar=25
    echo " -> No use of --minvar option, value set as default: 25"
fi
if [ -z "${minperc}" ]; then
    minperc=88
    echo " -> No use of --minperc option, value set as default:88"
fi
if [ -z "${maxgap}" ]; then
    maxgap=10
    echo " -> No use of --maxgap option, value set as default:10"
fi
if [ -n "${chrx}" ]; then
    echo " -> chrX will be included in the graphics."
fi
if [ -z "${chrx}" ]; then
    echo " -> chrX will NOT be included in the graphics."
    chrx="No"
fi
if [ -n "${extend}" ]; then
    echo " -> Homozygosity regions will be extended to nearest variant with maximum of $extend Mb."
fi
if [ -z "${extend}" ]; then
    echo " -> Homozygosity regions will be extended to nearest variant with maximum of 1.5 Mb."
    extend=1.5
fi

here="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

nb="$(bcftools query -l $vcf 2> $here.log | wc -l | cut -d" " -f1 )"
if [ "$nb" == "1" ]; then
    if [ -z "${pat}" ]; then
        pat="$(bcftools query -l $vcf 2> $here.log)"
        echo "## WARNING: No sample name provided through --pat option, name will be taken from the VCF: $pat"
    fi
fi
if [ "$nb" -gt "1" ]; then
    echo "## EROOR: Mutli-sample VCF file, please run AutoMap only on individual VCF files"
    exit 1
fi
if [ "$nb" == "0" ]; then
    echo "## ERROR: The input VCF format is incorrect."
    exit 1
fi

mkdir -p $out/$pat


#### conversion VCF

echo
echo "1) Parsing of VCF file and variant filtering"
numb="$(grep -v "#" $vcf | wc -l)"
echo " * $numb variants before filtering"

# removing variants with multiple additional alleles and variants in repeats

if [ "$genome" == "hg19" ]; then
   rep=$here/Ressources/repeats.bed
fi
if [ "$genome" == "hg38" ]; then
   rep=$here/Ressources/repeats_hg38.bed
fi
if [ -s $out/$pat/$pat.tsv ] || [ -s $out/$pat/$pat.clean.tsv ]; then
    :
else
    #cat $vcf > $vcf.norepeats.vcf
    bedtools subtract -a $vcf -b $rep -header > $vcf.norepeats.vcf
    grep -v "##" $vcf.norepeats.vcf | egrep -v "1/2"  > $out/$pat/$pat.tsv
fi

# parsing of the vcf file
if [ -s $out/$pat/$pat.clean.tsv ]; then
    :
else 
    perl $here/Scripts/parse_vcf.pl $out/$pat/$pat.tsv $out/$pat/$pat.clean.tsv 2> $here.log
    rm $out/$pat/$pat.tsv
fi

# filtering of variants on quality
grep -v "#" $out/$pat/$pat.clean.tsv  | awk -v percalthigh="$percalthigh" -v binomial="$binomial" -v percaltlow="$percaltlow" -F"\t" '{if($6 == "hom") print $0; if($6 == "het" && $11<=percalthigh && $12>=binomial && $11 >= percaltlow) print $0;}' | awk -v DP="$DP"  -F"\t" '{if($9 >= DP) print $0}' > $out/$pat/$pat.clean.qual.tsv

sort  -k1,1V -k2,2n -t $'\t' $out/$pat/$pat.clean.qual.tsv  > $out/$pat/$pat.clean.qual.sort.tsv
numb="$(grep -v "#" $out/$pat/$pat.clean.qual.sort.tsv | wc -l)"
echo " * $numb variants after filtering"

echo
echo "2) Detection of ROHs with sliding window, trimming and extension"
input=$out/$pat/$pat.clean.qual.sort.tsv
output_path=$out/$pat
output=$output_path/$pat.HomRegions
perl $here/Scripts/homo_regions_exome.pl $input $output $panel $panelname $window $windowthres $here/Scripts/trimming.sh $maxgap $here/Scripts/extend.sh $extend 2> $here.log

echo 
echo "3) Filtering of regions found and output to text file"
numb="$(grep -v "#" $output.tsv | wc -l)"
echo " * $numb regions before filtering"
if [ "$chrx" == "Yes" ]; then
    awk -v minsize="$minsize" -v minvar="$minvar" -v minperc="$minperc" -F "\t" '{if(($4>minsize && $5>minvar && $6>minperc) || $1 ~ /^#/) print $0}' $output.$panelname.tsv grep -v "chrY" > $output.strict.$panelname.tsv
    awk -v minsize="$minsize" -v minvar="$minvar" -v minperc="$minperc" -F "\t" '{if(($4>minsize && $5>minvar && $6>minperc) || $1 ~ /^#/) print $0}' $output.tsv grep -v "chrY" > $output.strict.tsv
else
    awk -v minsize="$minsize" -v minvar="$minvar" -v minperc="$minperc" -F "\t" '{if(($4>minsize && $5>minvar && $6>minperc) || $1 ~ /^#/) print $0}' $output.$panelname.tsv | grep -P -v "chrX|chrY" > $output.strict.$panelname.tsv
    awk -v minsize="$minsize" -v minvar="$minvar" -v minperc="$minperc" -F "\t" '{if(($4>minsize && $5>minvar && $6>minperc) || $1 ~ /^#/) print $0}' $output.tsv | grep -P -v "chrX|chrY" > $output.strict.tsv
fi
mv $output.strict.$panelname.tsv $output.$panelname.tsv
mv $output.strict.tsv $output.tsv
numb="$(grep -v "#" $output.tsv | wc -l)"
tot="$(grep -v "#" $output.tsv | grep -P -v "chrX|chrY" | cut -f4 | awk '{s+=$1} END {print s}')"
echo " * $numb regions after filtering with $tot Mb in total"

for file in $output.tsv $output.$panelname.tsv
do
    tot="$(grep -v "#" $file | grep -P -v "chrX|chrY" | cut -f4 | awk '{s+=$1} END {print s}')"
    echo "## INFO: $tot Mb are in Homozygous Regions (autosomal chromosomes)" >> $file;
    echo "## Variant filtering parameters used: DP=$DP, percaltlow=$percaltlow, percalthigh=$percalthigh, binomial=$binomial, maxgap=$maxgap" >> $file;
    echo "## Other parameters used: window=$window, windowthres=$windowthres, minsize=$minsize, minvar=$minvar, minperc=$minperc, chrX=$chrx, extend=$extend" >> $file;
done


echo
echo "4) Generating PDF"
size=$(more $output.tsv | grep INFO | awk -F " " '{print $3}' )
outputR=$output
if [ "$chrx" == "Yes" ]; then
    
    Rscript $here/Scripts/make_graph_chrX.R $pat $output.tsv $outputR.chrX.pdf $size 2> $here.log
else
    Rscript $here/Scripts/make_graph.R $pat $output.tsv $outputR.pdf $size 2> $here.log
fi

rm -f $out/$pat/$pat.clean.qual* $out/$pat/$pat.HomRegions.homozygosity* $here.log
