## Update Version 1.2
AutoMap v1.2 accepts list of VCF files in a text file with the --vcflist option.

## Update Version 1.1
AutoMap v1.1 is now compatible with VCF made by freebayes variant caller (with AO and not AD field).
It is now compatible with GATK, samtools mpileup, strelka and freebayes.

# AutoMap
Tool to find regions of homozygosity (ROHs) from sequencing data of human samples.

AutoMap for non-human data is available here: [[Link](https://github.com/mquinodo/AutoMap-nonhuman)] 

This software was written by Mathieu Quinodoz in the group of Prof. Rivolta from the IOB in Basel, Switzerland. It is presented at ASHG 2019. It was developped on Ubuntu 16.04.6 LTS (GNU/Linux 4.4.0-101-generic x86_64).

The online version can be found at: 
https://automap.iob.ch

## Prerequisites
+ BCFTools [[Link](https://samtools.github.io/bcftools/howtos/install.html)] (>= v1.9)
+ BEDTools [[Link](https://bedtools.readthedocs.io/en/latest/content/installation.html)] (>= v2.25.0)
+ Perl [[Link](https://www.perl.org/get.html)] (>= v5.22.0)
+ R [[Link](https://cran.r-project.org/mirrors.html)] (>= v3.2.0)

## Installation
The tool does not require compilation.

## Usage
This tool contains two main features: detection of ROHs for one or multiple individuals and additional analysis to find common ROHs.
### AutoMap for single individual
The main script AutoMap_v1.0.sh takes as input a VCF file which needs to contain GT (genotype) and AD (or AO) (allelic depths for the ref and alt alleles) or DP4 (number of high-quality ref-forward bases, ref-reverse, alt-forward and alt-reverse bases) fields for variants. The output is a text file containing the detected ROHs and a pdf file with the ROHs graphical representation.
It is called with bash:
```
bash AutoMap_v1.0.sh --vcf VCF_file --out output_directory --genome [hg19|hg38] [other options]
```
The approximate computation time per sample is 30 seconds for exome and few minutes for genome sequencing data.

### AutoMap for multiple individuals
The same script AutoMap_v1.0.sh can be used to compute ROHs for a list of individuals. VCFs files must be specified in --vcf option separated with commas.
```
bash AutoMap_v1.0.sh --vcf VCF1,VCF2,VCF3 --out output_directory --genome [hg19|hg38] [other options]
```
Ids can be specified through --id option also separated with commas. If IDs are not specified, they will be taken from VCF files directly.
```
bash AutoMap_v1.0.sh --vcf VCF1,VCF2,VCF3 --id ID1,ID2,ID3 --out output_directory --genome [hg19|hg38] [other options]
```

### Common ROHs to multiple individuals 
Autosomal ROHs common to multiple individuals can be computed with the --common option when multiple samples are analyzed simultaneously:
```
bash AutoMap_v1.0.sh --vcf VCF1,VCF2,VCF3 --out output_directory --genome [hg19|hg38] --common [other options]
```

#### Required arguments
Option | Value | Description
--- | --- | ---
--vcf | STRING | Single VCF file or list of VCF files of the individual(s) to analyse, separated with commas (--vcflist can be used instead)
--genome | [hg19/hg38] | Genome build used in the VCF file
--out | STRING | Output directory in which the directory with the individual name will be created

#### Optional arguments
Option | Default | Value | Description
--- | --- | --- | ---
--common | - | - | If multiple samples are analysed, common ROHs will be computed.
--id | from VCF | STRING | ID of the individual analysed or list of individuals analysed separated by commas (same order as in the list of VCF files)
--panel | None | STRING | File containing a gene or region panel (see panel format)
--panelname | None | STRING | Name of the panel file for output
--DP | 8 | 0-99 | Minimal depth for variants
--binomial | 0.000001 | 0-1 | Minimal p-value for Binomial test for reference and alternative alleles counts
--percaltlow | 0.25 | 0-1 | Minimal allternative reads ratio for heterozygous variants
--percalthigh | 0.75 | 0-1 | Maximal allternative reads ratio for heterozygous variants
--window | 7 | 3-999 | Size of the sliding window
--windowthres | 5 | 1-999 | Threshold of homozygous variants in the window
--minsize | 1.0 | 0-99 | Minimal size of detected ROH [Mb]
--minvar | 25 | 1-999 | Minimal number of variant in detected ROH
--minperc | 88 | 0-100 | Minimal percentage of homozygous variants in detected ROH
--maxgap | 10 | 0-1000 | Maximal gap allowed between two variants in one ROH [Mb]
--extend | 1.0 | 0-100 | Maximal extension at both ROH boundaries (if no heterozygous SNPs closer)
--chrX   | - | - | Outputs will contain chromosome X
--vcflist | - | - | list of vcf files to analyze (compatible with the --common option)
--multivcf | - | - | To use for analysis of a multi-sample VCF (cannot be used with --id and --common options). It will analyze each sample separately and generate individual VCFs that can be further used for identification of common ROHs.

#### Panel format
The file provided through --panel option should contain four tab-separated fields with:

1° Feature name (gene name for example)

2° Chromosome (with “chr” notation: chr1, chr2, ...)

3° Beginning of the gene

4° End of the gene

An example can be found in the Resources (RetNet_hg19.txt). The genome reference must be the same that the analysed VCF files.

#### Outputs for single individual
The output text file (.HomRegions.tsv) and pdf file (.HomRegions.pdf) will be place in a folder with the name of the individual analysed in the output directory.

#### Outputs for multiple individuals
The output will be the same as individual output with an additional one containing the common regions of all analysed individuals if option --common has been used.

#### Test sample
A testing VCF with random variants can be found in the Test folder. The following command can be done to test AutoMap:
```
AUTOMAP_HOME=<path/to/automap/script>
bash $AUTOMAP_HOME/AutoMap_v1.0.sh
  --vcf $AUTOMAP_HOME/Test/TestSample.vcf
  --genome hg19
  --out $AUTOMAP_HOME/Results-test
```
And the output can be compared with the expected output files in the Test folder.
