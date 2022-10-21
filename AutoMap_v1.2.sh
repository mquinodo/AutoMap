#!/bin/bash

usage() { echo "## ERROR: Usage: $0 [--vcf <string>] [--vcflist <string>] [--genome <hg19|hg38>] [--out <string>] [--common] [--id <string>] [--panel <string>] [--panelname <string>] [--DP <0-99>] [--binomial <0-1.0>] [--percaltlow <0-1.0>] [--percalthigh <0-1.0>] [--window <3-999>] [--windowthres <1-999>] [--minsize <0-99>] [--minvar <1-999>] [--minperc <0-100>] [--maxgap <0-1000Mb>] [--chrX] [--extend <0-100Mb>]. Exit." 1>&2; exit 1; }
numbervar() { echo "## ERROR: Less than 10,000 variants ($nbvar detected variants) with AD (or AO) and DP available. Exit." 1>&2; exit 1; }
numbervar2() { echo "## ERROR: Less than 10,000 variants with good quality detected (only $nbvar passing QC). Exit." 1>&2; exit 1; }
multivcf() { echo "## ERROR: Multi-sample VCF file, please run AutoMap with --multivcf option. Exit." 1>&2; exit 1; }
multivcf2() { echo "## ERROR: Please analyze only one mutli-sample VCF file at a time. Exit." 1>&2; exit 1; }
multivcf3() { echo "## ERROR: Please do not use --id or --common option with mutli-sample VCF file. Exit." 1>&2; exit 1; }
nobcftools() { echo "## ERROR: bcftools lower than v1.9 -> Please Update! Exit." 1>&2; exit 1; }
nobedtools() { echo "## ERROR: bedtools lower than v2.24.0 -> Please Update! Exit." 1>&2; exit 1; }
noperl() { echo "## ERROR: perl lower than v5.22.0 -> Please Update! Exit." 1>&2; exit 1; }
noRversion() { echo "## ERROR: R lower than v3.2.0 -> Please Update! Exit." 1>&2; exit 1; }
novcf() { echo "## ERROR: You need to provide an input vcf through --vcf option. Exit." 1>&2; exit 1; }
nogenome() { echo "## ERROR: You need to provide the genome version through --genome option (hg19 or hg38). Exit." 1>&2; exit 1; }
nooutput() { echo "## ERROR: You need to provide an output directory through --out option. Exit." 1>&2; exit 1; }
nobothpanel() { echo "## ERROR: No panel used, you need to provide both a panel name and panel through --panelname and --panel options. Exit." 1>&2; exit 1; }
nosameidvcf() { echo "## ERROR: Not the same number of ids and vcf files. Exit." 1>&2; exit 1; }
emptypanel() { echo "## ERROR: Panel file is empty. Exit." 1>&2; exit 1; }
nopanelfile() { echo "## ERROR: Panel file does not exist. Exit." 1>&2; exit 1; }
vcfbadformat() { echo "## ERROR: The input VCF format is incorrect ('bcftools query -l' was unsuccessful). Exit." 1>&2; exit 1; }
vcflistempty() { echo "## ERROR: The list of VCF files is empty. Exit." 1>&2; exit 1; }

currentver="$(bcftools -v | head -n1 | cut -d" " -f2)"
requiredver="1.9"
 if [ "$(printf '%s\n' "$requiredver" "$currentver" | sort -V | head -n1)" = "$requiredver" ]; then 
        echo "# bcftools higher or equal to v1.9"
 else
        nobcftools
 fi

currentver="$(bedtools --version | cut -d" " -f2)"
requiredver="v2.24.0"
 if [ "$(printf '%s\n' "$requiredver" "$currentver" | sort -V | head -n1)" = "$requiredver" ]; then 
        echo "# bedtools higher or equal to v2.24.0"
 else
        nobedtools
 fi

currentver="$(perl -v | grep "This is perl" | cut -d"(" -f2 | cut -d")" -f1)"
requiredver="v5.22.0"
 if [ "$(printf '%s\n' "$requiredver" "$currentver" | sort -V | head -n1)" = "$requiredver" ]; then 
        echo "# perl higher or equal to v5.22.0"
 else
        noperl
 fi

currentver="$(R --version | grep "R version" | cut -d" " -f3)"
requiredver="3.2.0"
 if [ "$(printf '%s\n' "$requiredver" "$currentver" | sort -V | head -n1)" = "$requiredver" ]; then 
        echo "# R higher or equal to v3.2.0"
 else
        noRversion
 fi

while getopts ":-:" o; do
    case "${o}" in
    	-)  
        	case $OPTARG in
                vcf)
                    vcf="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    vcfs=(${vcf//,/ })
                    numbervcf=${#vcfs[@]}
                    ;;
                vcflist)
                    vcflist="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    readarray -t vcfs < $vcflist
                    numbervcf=${#vcfs[@]}
                    ;;
                genome)
                    genome="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    [ "$genome" == "hg19" ] || [ "$genome" == "hg38" ] || usage
                    ;;
                out)
                    out="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                id)
                    id="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ids=(${id//,/ })
                    numberid=${#ids[@]}
                    allid=$id
                    identer="Yes"
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
                multivcf)
                    multivcf="Yes"
                    ;;
                common)
                    common="Yes"
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
if [ -z "${vcf}" ] && [ "${numbervcf}" == "0" ]; then
    novcf
fi
if [ "${numbervcf}" == "0" ]; then
    vcflistempty
fi
if [ -z "${genome}" ]; then
    nogenome
fi
if [ -z "${out}" ]; then
    nooutput
fi
if [ ! -z "${panel}" ] && [ -z "${panelname}" ]; then
    nobothpanel
fi
if [ -z "${panel}" ] && [ ! -z "${panelname}" ]; then
    nobothpanel
fi
if [ "${numbervcf}" != "${numberid}" ] && [ ! -z "${id}" ]; then
    nosameidvcf
fi
if [ ! -z "${panel}" ]; then
    if [ -f "${panel}" ]
    then
        if [ -s "${panel}" ]
        then
            echo ""
        else
            emptypanel
        fi
    else
        nopanelfile
    fi
fi
if [ -z "${panel}" ]; then
    panel="NA"
fi
if [ -z "${panelname}" ]; then
    panelname="NA"
fi
if [ -z "${DP}" ]; then
    DP=8
    echo " -> No use of --DP option, value set as default: 8"
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
    minsize=1
    echo " -> No use of --minsize option, value set as default: 1"
fi
if [ -z "${minvar}" ]; then
    minvar=25
    echo " -> No use of --minvar option, value set as default: 25"
fi
if [ -z "${minperc}" ]; then
    minperc=88
    echo " -> No use of --minperc option, value set as default: 88"
fi
if [ -z "${maxgap}" ]; then
    maxgap=10
    echo " -> No use of --maxgap option, value set as default: 10"
fi
if [ -n "${chrx}" ]; then
    echo " -> chrX will be included in the analysis and in the graphics."
fi
if [ -n "${common}" ]; then
    echo " -> ROHs common to all samples will be computed."
fi
if [ -z "${chrx}" ]; then
    echo " -> chrX will NOT be included in the analysis and in the graphics."
    chrx="No"
fi
if [ -z "${multivcf}" ]; then
    multivcf="No"
fi
if [ -n "${extend}" ]; then
    echo " -> Homozygosity regions will be extended to nearest variant with maximum of $extend Mb."
fi
if [ -z "${extend}" ]; then
    echo " -> Homozygosity regions will be extended to nearest variant with maximum of 1 Mb."
    extend=1
fi


###### multi-sample VCFs #######

if [ "$multivcf" == "Yes" ]; then
    if [ "$numbervcf" -gt "1" ]; then
        multivcf2
    fi
    if [ "$identer" == "Yes" ]; then
        multivcf3
    fi
    if [ "$common" == "Yes" ]; then
        multivcf3
    fi

    nb="$(bcftools query -l $vcf 2> $here/.log | wc -l | cut -d" " -f1 )"

    for (( k=0; k<$nb; k++ ))
    do
        numb=$(($k+1))
        pat=$(bcftools query -l $vcf | head -n $numb | tail -1)
        mkdir -p $out/$pat
        bcftools view -c1 -Ov -s $pat -o $out/$pat/$pat.individual.vcf $vcf

        here="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

        echo "## Launching analyis for sample: $pat"

        if [[ "$chrx" == "Yes" ]] && [[ "$panelname" != "NA" ]]; then
            bash $here/AutoMap_v1.2.sh --id $pat --chrX --vcf $out/$pat/$pat.individual.vcf --out $out --genome $genome --panel $panel --panelname $panelname --DP $DP --binomial $binomial --percaltlow $percaltlow --percalthigh $percalthigh --window $window --windowthres $windowthres --minsize $minsize --minvar $minvar --minperc $minperc --maxgap $maxgap --extend $extend
        fi
        if [[ "$chrx" == "Yes" ]] && [[ "$panelname" == "NA" ]]; then
            bash $here/AutoMap_v1.2.sh --id $pat --chrX --vcf $out/$pat/$pat.individual.vcf --out $out --genome $genome --DP $DP --binomial $binomial --percaltlow $percaltlow --percalthigh $percalthigh --window $window --windowthres $windowthres --minsize $minsize --minvar $minvar --minperc $minperc --maxgap $maxgap --extend $extend
        fi
        if [[ "$chrx" == "No" ]] && [[ "$panelname" != "NA" ]]; then
            bash $here/AutoMap_v1.2.sh --id $pat --vcf $out/$pat/$pat.individual.vcf --out $out --genome $genome --panel $panel --panelname $panelname --DP $DP --binomial $binomial --percaltlow $percaltlow --percalthigh $percalthigh --window $window --windowthres $windowthres --minsize $minsize --minvar $minvar --minperc $minperc --maxgap $maxgap --extend $extend
        fi
        if [[ "$chrx" == "No" ]] && [[ "$panelname" == "NA" ]]; then
            bash $here/AutoMap_v1.2.sh --id $pat --vcf $out/$pat/$pat.individual.vcf --out $out --genome $genome --DP $DP --binomial $binomial --percaltlow $percaltlow --percalthigh $percalthigh --window $window --windowthres $windowthres --minsize $minsize --minvar $minvar --minperc $minperc --maxgap $maxgap --extend $extend
        fi

    done
fi

###### LOOP on VCFs ######
if [ "$multivcf" == "No" ]; then
    for (( k=0; k<$numbervcf; k++ ))
    do
        vcf=${vcfs[$k]}
        id=${ids[$k]}

        here="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

        nbvar=$(grep -v "#" $vcf | grep -P "AD|DP4|AO" | grep GT | wc -l)

        nb="$(bcftools query -l $vcf 2> $here/.log | wc -l | cut -d" " -f1 )"
        if [ "$nb" == "1" ]; then
            if [ -z "${id}" ]; then
                id="$(bcftools query -l $vcf 2> $here/.log)"
                if [ "$k" == "0" ]; then
                    allid=$id
                fi
                if [ "$k" -gt "0" ]; then
                    allid="$allid,$id"
                fi
                echo "## WARNING: No sample name provided through --id option, name will be taken from the VCF: $id"
            fi

            if [ "$nbvar" -lt "10000" ]; then
                numbervar
            fi
        fi
        if [ "$nb" -gt "1" ]; then
            multivcf
        fi
        if [ "$nb" == "0" ]; then
            vcfbadformat
        fi

        mkdir -p $out/$id


        #### conversion VCF

        echo
        echo "1) Parsing of VCF file and variant filtering"
        numb="$(grep -v "#" $vcf | wc -l)"
        echo " * $numb variants before filtering"

        # removing variants with multiple additional alleles and variants in repeats
        if [ "$genome" == "hg19" ]; then
            rep=$here/Resources/repeats.bed
            if [ ! -f "$rep" ]; then
                bash $here/Resources/merge_repeats.sh $here/Resources/repeats.part1.bed.gz $here/Resources/repeats.part2.bed.gz $rep
            fi
        fi
        if [ "$genome" == "hg38" ]; then
            rep=$here/Resources/repeats_hg38.bed
            if [ ! -f "$rep" ]; then
                bash $here/Resources/merge_repeats.sh $here/Resources/repeats_hg38.part1.bed.gz $here/Resources/repeats_hg38.part2.bed.gz $rep
            fi
        fi
        if [ -s $out/$id/$id.tsv ] || [ -s $out/$id/$id.clean.tsv ]; then
            :
        else
            awk '{if($0 !~ /^#/ && $0 !~ /^chr/) print "chr"$0; else print $0}' $vcf > $vcf.chr
            bedtools subtract -a $vcf.chr -b $rep -header > $vcf.norepeats.vcf
            grep -v "##" $vcf.norepeats.vcf | egrep -v "1/2"  > $out/$id/$id.tsv
            rm $vcf.norepeats.vcf
        fi

        # parsing of the vcf file
        if [ -s $out/$id/$id.clean.tsv ]; then
            :
        else 
            perl $here/Scripts/parse_vcf_v1.1.pl $out/$id/$id.tsv $out/$id/$id.clean.tsv 2> $here/.log
            rm $out/$id/$id.tsv
        fi

        # filtering of variants on quality
        grep -v "#" $out/$id/$id.clean.tsv  | awk -v percalthigh="$percalthigh" -v binomial="$binomial" -v percaltlow="$percaltlow" -F"\t" '{if($6 == "hom") print $0; if($6 == "het" && $11<=percalthigh && $12>=binomial && $11 >= percaltlow) print $0;}' | awk -v DP="$DP"  -F"\t" '{if($9 >= DP) print $0}' > $out/$id/$id.clean.qual.tsv

        nbvar=$(cat $out/$id/$id.clean.qual.tsv | wc -l)

        if [ "$nbvar" -lt "10000" ]; then
                numbervar2
        fi


        sort  -k1,1V -k2,2n -t $'\t' $out/$id/$id.clean.qual.tsv  > $out/$id/$id.clean.qual.sort.tsv
        numb="$(grep -v "#" $out/$id/$id.clean.qual.sort.tsv | wc -l)"
        echo " * $numb variants after filtering"

        echo
        echo "2) Detection of ROHs with sliding window, trimming and extension"
        input=$out/$id/$id.clean.qual.sort.tsv
        output_path=$out/$id
        output=$output_path/$id.HomRegions
        perl $here/Scripts/homo_regions.pl $input $output $panel $panelname $window $windowthres $here/Scripts/trimming.sh $maxgap $here/Scripts/extend.sh $extend 2> $here/.log

        echo 
        echo "3) Filtering of regions found and output to text file"
        numb="$(grep -v "#" $output.tsv | wc -l)"
        echo " * $numb regions before filtering"
        if [ "$chrx" == "Yes" ]; then

            if [ "$panelname" != "NA" ]; then
                head -n1 $output.$panelname.tsv > $output.$panelname.head.tsv
                awk -v minsize="$minsize" -v minvar="$minvar" -v minperc="$minperc" -F "\t" '{if(($4>minsize && $5>minvar && $6>minperc) || $1 ~ /^#/) print $0}' $output.$panelname.tsv | grep -v "chrY" | tail -n+2 | sort -k1,1V -k2,2n > $output.strict.$panelname.temp.tsv
                cat $output.$panelname.head.tsv $output.strict.$panelname.temp.tsv > $output.strict.$panelname.tsv
                rm $output.$panelname.head.tsv $output.strict.$panelname.temp.tsv
            fi
            head -n1 $output.tsv > $output.head.tsv
            awk -v minsize="$minsize" -v minvar="$minvar" -v minperc="$minperc" -F "\t" '{if(($4>minsize && $5>minvar && $6>minperc) || $1 ~ /^#/) print $0}' $output.tsv | grep -v "chrY" | tail -n+2 | sort -k1,1V -k2,2n > $output.strict.temp.tsv
            cat $output.head.tsv $output.strict.temp.tsv > $output.strict.tsv
            rm $output.head.tsv $output.strict.temp.tsv
        else
            if [ "$panelname" != "NA" ]; then
                head -n1 $output.$panelname.tsv > $output.$panelname.head.tsv
                awk -v minsize="$minsize" -v minvar="$minvar" -v minperc="$minperc" -F "\t" '{if(($4>minsize && $5>minvar && $6>minperc) || $1 ~ /^#/) print $0}' $output.$panelname.tsv | grep -P -v "chrX|chrY" | tail -n+2 | sort -k1,1V -k2,2n > $output.strict.$panelname.temp.tsv
                cat $output.$panelname.head.tsv $output.strict.$panelname.temp.tsv > $output.strict.$panelname.tsv
                rm $output.$panelname.head.tsv $output.strict.$panelname.temp.tsv
            fi
            head -n1 $output.tsv > $output.head.tsv
            awk -v minsize="$minsize" -v minvar="$minvar" -v minperc="$minperc" -F "\t" '{if(($4>minsize && $5>minvar && $6>minperc) || $1 ~ /^#/) print $0}' $output.tsv | grep -P -v "chrX|chrY" | tail -n+2 | sort -k1,1V -k2,2n > $output.strict.temp.tsv
            cat $output.head.tsv $output.strict.temp.tsv > $output.strict.tsv
            rm $output.head.tsv $output.strict.temp.tsv
        fi
        mv $output.strict.tsv $output.tsv
        numb="$(grep -v "#" $output.tsv | wc -l)"
        tot="$(grep -v "#" $output.tsv | grep -P -v "chrX|chrY" | cut -f4 | awk '{s+=$1} END {print s}')"
        if [ "$tot" == "" ]; then
            tot=0
        fi
        echo " * $numb regions after filtering with $tot Mb in total"

        file=$output.tsv
        tot="$(grep -v "#" $file | grep -P -v "chrX|chrY" | cut -f4 | awk '{s+=$1} END {print s}')"
        if [ "$tot" == "" ]; then
            tot=0
        fi
        echo "## INFO: $tot Mb are in Homozygous Regions (autosomal chromosomes)" >> $file;
        echo "## AutoMap v1.2 used for analysis" >> $file;
        echo "## Variant filtering parameters used: DP=$DP, percaltlow=$percaltlow, percalthigh=$percalthigh, binomial=$binomial, maxgap=$maxgap" >> $file;
        echo "## Other parameters used: window=$window, windowthres=$windowthres, minsize=$minsize, minvar=$minvar, minperc=$minperc, chrX=$chrx, extend=$extend" >> $file;

        if [ "$panelname" != "NA" ]; then
            file=$output.$panelname.tsv
            tot="$(grep -v "#" $file | grep -P -v "chrX|chrY" | cut -f4 | awk '{s+=$1} END {print s}')"
            if [ "$tot" == "" ]; then
            tot=0
            fi
            echo "## INFO: $tot Mb are in Homozygous Regions (autosomal chromosomes)" >> $file;
            echo "## AutoMap v1.2 used for analysis" >> $file;
            echo "## Variant filtering parameters used: DP=$DP, percaltlow=$percaltlow, percalthigh=$percalthigh, binomial=$binomial, maxgap=$maxgap" >> $file;
            echo "## Other parameters used: window=$window, windowthres=$windowthres, minsize=$minsize, minvar=$minvar, minperc=$minperc, chrX=$chrx, extend=$extend" >> $file;
        fi

        echo
        echo "4) Generating PDF"
        size=$(cat $output.tsv | grep INFO | awk -F " " '{print $3}' )
        outputR=$output
        if [ "$chrx" == "Yes" ]; then
            
            Rscript $here/Scripts/make_graph_chrX.R $id $output.tsv $outputR.pdf $size 2> $here/.log
        else
            Rscript $here/Scripts/make_graph.R $id $output.tsv $outputR.pdf $size 2> $here/.log
        fi

        rm -f $out/$id/$id.clean* $out/$id/$id.HomRegions.homozygosity* $here/.log $vcf.chr
    done

    # Regions common to all

    if [ "$numbervcf" -gt "1" ] && [ "$common" == "Yes" ] ; then
        echo
        echo "5) Computing common ROHs"
        output=${allid//,/_}
        if [ "$panelname" != "NA" ]; then
            bash $here/Scripts/common_analysis.sh --res $out --name $output --ids $allid --panelname $panelname --panel $panel
            echo "## AutoMap v1.2 used for analysis" >> $out/$output/$output.HomRegions.tsv;
            echo "## Variant filtering parameters used: DP=$DP, percaltlow=$percaltlow, percalthigh=$percalthigh, binomial=$binomial, maxgap=$maxgap" >> $out/$output/$output.HomRegions.tsv;
            echo "## Other parameters used: window=$window, windowthres=$windowthres, minsize=$minsize, minvar=$minvar, minperc=$minperc, chrX=$chrx, extend=$extend" >> $out/$output/$output.HomRegions.tsv;
            echo "## AutoMap v1.2 used for analysis" >> $out/$output/$output.HomRegions.$panelname.tsv;
            echo "## Variant filtering parameters used: DP=$DP, percaltlow=$percaltlow, percalthigh=$percalthigh, binomial=$binomial, maxgap=$maxgap" >> $out/$output/$output.HomRegions.$panelname.tsv;
            echo "## Other parameters used: window=$window, windowthres=$windowthres, minsize=$minsize, minvar=$minvar, minperc=$minperc, chrX=$chrx, extend=$extend" >> $out/$output/$output.HomRegions.$panelname.tsv;
        fi 
        if [ "$panelname" == "NA" ]; then
            bash $here/Scripts/common_analysis.sh --res $out --name $output --ids $allid
            echo "## AutoMap v1.2 used for analysis" >> $out/$output/$output.HomRegions.tsv;
            echo "## Variant filtering parameters used: DP=$DP, percaltlow=$percaltlow, percalthigh=$percalthigh, binomial=$binomial, maxgap=$maxgap" >> $out/$output/$output.HomRegions.tsv;
            echo "## Other parameters used: window=$window, windowthres=$windowthres, minsize=$minsize, minvar=$minvar, minperc=$minperc, chrX=$chrx, extend=$extend" >> $out/$output/$output.HomRegions.tsv;
        fi 
    fi
fi

rm -f $here/.log
