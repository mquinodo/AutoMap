# AutoMap
Tool to find regions of homozygosity (ROHs) from exome sequencing data.

This software was written by Quinodoz *et al.* from the University of Lausanne, Switzerland. It was presented at XX and is published in XXX.
The online version can be found at: 
```
www.XXX.ch
```

## Requirements
This code need BCFTools, BEDTools, perl and R installed for proper function.

## Usage
This tool contain two main features: detection of ROHs and family analyis to find common ROHs between individuals.
### AutoMap
The main script AutoMap_v1.0.sh is taking as input a VCF file which need contain GT (genotype) and AD (allelic depths for the ref and alt alleles) or DP4 (number of high-quality ref-forward bases, ref-reverse, alt-forward and alt-reverse bases) fields for variants.
It is outputing a text file containing the detected ROHs and a pdf file for graphical representation.

It is called with bash:
```
bash AutoMap_v1.0.sh --vcf VCF_file --out output_directory --genome [hg19|hg38] [other options]
```

#### Mandatory options
+ --vcf
+ --genome
+ --out

#### Other options
+ --pat
+ --panel
+ --panelname
+ --DP
+ --binomial
+ --percaltlow
+ --percalthigh
+ --window
+ --windowthres
+ --minsize
+ --minvar
+ --miperc + --maxgap
+ --extend
+ --chrX

#### Output

### Family analysis
The script is extracting ROHs common to affected individuals and not present in healthy family members. It is outputing a text file containing the detected ROHs and a pdf file for graphical representation.

It is called with bash:
```
bash family_analysis.sh --res results_directory --name output_name --affected list_of_affected [other options]
```

#### Mandatory options
+ --res
+ --name
+ --affected
#### Other options
+ --healthy
+ --panel
+ --panelname
#### Output

