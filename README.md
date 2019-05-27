# AutoMap
Tool to find regions of homozygosity (ROHs) from exome sequencing data.

This software was written by Quinodoz *et al.* from the University of Lausanne, Switzerland. It was presented at XX and is published in XXX.
The online version can be found at: 
```
www.XXX.ch
```

## Prerequisites
+ BCFTools [[Link](https://samtools.github.io/bcftools/howtos/install.html)] (>= v1.9)
+ BEDTools [[Link](https://bedtools.readthedocs.io/en/latest/content/installation.html)] (>= v2.24.0)
+ Perl with Text::CSV and List::Util modules [[Link](https://www.perl.org/get.html)] (>= v5.22.0)
+ R [[Link](https://cran.r-project.org/mirrors.html)] (>= v3.4.1)

## Usage
This tool contains two main features: detection of ROHs and family analysis to find common ROHs between individuals.
### AutoMap
The main script AutoMap_v1.0.sh is taking as input a VCF file which need contain GT (genotype) and AD (allelic depths for the ref and alt alleles) or DP4 (number of high-quality ref-forward bases, ref-reverse, alt-forward and alt-reverse bases) fields for variants.
It is outputing a text file containing the detected ROHs and a pdf file for graphical representation.

It is called with bash:
```
bash AutoMap_v1.0.sh --vcf VCF_file --out output_directory --genome [hg19|hg38] [other options]
```

#### Mandatory options
Option | Value | Description
--- | --- | ---
--vcf | STRING | VCF file of the individual to analyze
--genome | [hg19/hg38] | Genome build used in the VCF file
--out | STRING | Output directory in which the directory with the individual name will be created

#### Other options
Option | Default | Value | Description
--- | --- | --- | ---
--pat | from VCF | STRING | Name of the individual analyzed
--panel | None | STRING | File containg a gene or region panel (see panel format)
--panelname | None | STRING | Name of the panel file for output
--DP | 15 | 0-99 | Minimal depth for variants
--binomial | 0.000001 | 0-1 | Minimal p-value for Binomial test for reference and alternative alleles counts
--percaltlow | 0.25 | 0-1 | Minimal allternative reads ratio for heterozygous variants
--percalthigh | 0.75 | 0-1 | Maximal allternative reads ratio for heterozygous variants
--window | 7 | 3-999 | Size of the sliding window
--windowthres | 5 | 1-999 | Threshold of homozygous variants in the window
--minsize | 2.0 | 0-99 | Minimal size of detected ROH [Mb]
--minvar | 25 | 1-999 | Minimal number of variant in detected ROH
--miperc | 88 | 0-100 | Minimal percentage of homozygous variants in detected ROH
--maxgap | 10 | 0-1000 | Maximal gap allowed between two variants in one ROH [Mb]
--extend | 1.5 | 0-100 | Maximal extension at both ROH boundaries (if no heterozygous SNPs closer)
--chrX   | - | - | Outputs will contain chromosome X  

#### Panel format
The file provided throught --panel option should contain four tab-separted fields with:

1° Feature name (gene name for example)

2° Chrmosome (chr1, chr2,...)

3° Beginning position

4° End position

#### Outputs
The output text file (.HomRegions.tsv) and pdf file (.HomRegions.pdf) will be place in a folder with the name of the individual analyzed in the output directory.

#### Test
AutoMap comes with a test VCF file to check for correct software installations.
The command to use inside AutoMap directory is:
```
bash AutoMap_v0.1.sh --vcf Test/TestSample.vcf --out Test/ --genome hg19
```
It will produce pdf and text results in a new directory: TestSample. The pdf and text results then can be compared to the expected one which are in Test directory (TestSample.HomRegions_correct.pdf and TestSample.HomRegions_correct.tsv).


### Family analysis
The script is extracting ROHs common to affected individuals and not present in healthy family members. It is outputing a text file containing the detected ROHs and a pdf file for graphical representation.

It is called with bash:
```
bash family_analysis.sh --res results_directory --name output_name --affected list_of_affected [other options]
```

#### Mandatory options
Option | Value | Description
--- | --- | ---
--res | STRING | Directory containing the individual results (same as --out for AutoMap.sh)
--name | STRING | Output name of the analysis
--affected | STRING | Comma-separated string of affected individuals (Ex. "Pat1,Pat2,Pat3")
#### Other options
Option | Default | Value | Description
--- | --- | --- | ---
--healthy | None | STRING | Comma-separated string of affected individuals (Ex. "Pat1,Pat2,Pat3")
--panel | None | STRING | File containg a gene or region panel (see panel format)
--panelname | None | STRING | Name of the panel file for output
#### Output
Same as AutoMap.sh
