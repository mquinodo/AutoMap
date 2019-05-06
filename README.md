# AutoMap
Suite of tools for analysing off-target reads to find CNVs, homozygous regions, and shared haplotypes.

This software was written by Matthew Wakeling at the University of Exeter, and was presented at the 2017 ASHG meeting in Orlando, Florida.

## Running Java
This code requires the htsjdk library and the JAMA matrix maths library. The easiest way to get everything required is to download the GATK Jar. All operations require this GATK jar and the SavvySuite java to be in the Java classpath, in order for Java to find it.

This can be done in two ways. The first option is to set the CLASSPATH environment variable:
```
export CLASSPATH=/path/to/GenomeAnalysisTK.jar:/path/to/SavvySuite/directory
```
The ":" character separates the two parts of this path, to specify that code can be found in the two places. The second option is to specify the "-cp" option every time you run java, like this:
```
java -cp /path/to/GenomeAnalysisTK.jar:/path/to/SavvySuite/directory blah blah blah
```
For all subsequent code fragments, where "java" or "javac" is specified, it is assumed that the classpath is correctly configured as specified above, either by adding the "-cp" option or using the CLASSPATH environment variable.

If you have a large server, it is also sensible to add the -XX:ParallelGCThreads=2 -XX:ConcGCThreads=2 options to java, to prevent it creating too many garbage collection threads. This is a small performance enhancement, and is done like this:
```
java -XX:ParallelGCThreads=2 -XX:ConcGCThreads=2 blah blah blah
```

## Compiling
Compiling the code is then done by:
```
javac *.java
```
in the SavvySuite directory.

## Usage
This suite contains three separate tools for analysing off-target reads.
### SavvyCNV
This software analyses the read depth of off-target reads to detect CNVs. It requires a reasonable number of samples sequenced using the same method. (Don't mix samples sequenced using different methods - it won't work well.) The sample data must be provided in aligned BAM files. First, each BAM file must be converted to a coverage summary file, with the following command:
```
java -Xmx1g CoverageBinner sample.bam >sample.coverageBinner
```
This step can be performed in parallel on each sample, and will produce a file approximately 7MB in size. The operation requires very little RAM.

To perform the analysis, the following command should be used:
```
java -Xmx30g SavvyCNV -d (size) *.coverageBinner >cnv_list.csv 2>log_messages.txt
```
The <size> parameter is the size of the chunks that the genome is split into. If you have targeted sequencing with three million reads and about 50% off target reads, then a chunk size of 200,000 is appropriate. It is sensible to process male and female samples separately if CNVs in the X/Y chromosomes are to be detected.
In addition, the following arguments can be provided:
+ -trans (transition probability) - This is the transition probability to use for the Viterbi algorithm. The default is 0.00001. To increase the sensitivity and false positive rate, increase this parameter.
+ -cutoff (noise cutoff) - This is the noise threshold above which a chunk of the genome will be excluded from analysis. The default of 0.25 is probably best in most situations.
+ -g - Switches on the generation of graphs for all samples that have a detected CNV. The generated graphs will be placed in the same directory as the *.coverageBinner file.
+ -a - The same as ```-g```, but produces a graph for every sample.
+ -cytoBands (file) - The location of a cytoBands file for the genome reference, for example downloadable from http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz - this will be plotted in the graph behind the data.
+ -mosaic - Switches the software into mosaic mode. Normally, the state probability calculations assume that the relative dosage is either <=0.5, 1, or >=1.5. Dosage levels must cross the mid-point between 1 and 0.5 or 1.5 before they become evidence of a CNV. This increases sensitivity and specificity at the cost of being able to detect mosaic CNVs. With this switch, mosaics can be detected. The size parameter will need to be increased, and small CNVs will not be detected as effectively.
+ -sv (number) - This changes the number of singular vectors that are removed for noise reduction. The default is 5. This must be less than the number of samples.
+ -minReads (number) - This sets the minimum number of reads that a genome chunk must have on average across the samples in order to be analysed. The default is 20.
+ -minProb (number) - This sets the greatest (lowest) probability that a single chunk can contribute to a CNV. The number is a phred score. This is the largest quality score that a very small CNV can have.

The output cnv_list.csv contains a tab-separated list of detected CNVs. The columns in the output are:
1. Chromosome
2. CNV start position
3. CNV end position
4. Deletion/duplication
5. Number of genome chunks used as evidence for the CNV
6. Width of CNV in chunks (including gaps containing noisy chunks, such as centromeres)
7. The phred score of the CNV
8. The phred score divided by the width of the CNV in chunks. Most valid CNVs have a value more than ten in this column
9. The relative dosage in the CNV, so 0.5 for a heterozygous deletion etc.
10. The filename containing the sample summary.

The log_messages.txt file contains log messages, and also a summary of each sample. Each sample has a line containing the following columns:
1. Noisyness of the sample before excluding known CNVs
2. Noisyness of the sample after excluding known CNVs (the CNV calling is performed twice, as the calling depends on the noisyness of the sample). This figure should be below 0.2 for good results. If too many of the samples have a high noisyness, then increase the chunk size parameter (-d <size>).
3. Number of deletions found in the sample
4. Number of duplications found in the sample
5. The number of reads in the sample.
6. The filename containing the sample summary.

### PrepareLinkageData
This software is used to pre-process linkage disequilibrium data, to get it into a format that can be read quickly by SavvyHomozygosity and SavvySharedHaplotypes. It requires a vcf file containing whole genome genotype data for many samples (a few hundred would be appropriate). It can be run as follows:
```
java -Xmx5g PrepareLinkageData whole_genome.vcf >linkage_data
```

### SavvyHomozygosity
This software analyses the off-target reads to determine homozygous regions of the genome for a single sample. It can be run as follows:
```
java -Xmx5g SavvyHomozygosity linkage_data sample.bam >sample.bed
```
The sample.bam contains the sample targeted sequencing data. The software produces a BED file containing homozygous regions. The analysis requires at least a million off-target reads, and preferably more, so a sample with three million reads and about 50% off-target is about right. The software searches for read pairs, so the amount of detail that can be extracted (and consequently the time taken) is proportional to the square of the number of off-target reads. Samples with too few reads do not produce a valid result.

### SavvySharedHaplotypes
This software analyses the off-target reads from two samples to determine areas with homozygous shared haplotypes, using the same methodology as SavvyHomozygosity. It can be run as follows:
```
java -Xmx5g SavvySharedHomozygosity linkage_data sample1.bam sample2.bam >shared.bed
```

### SavvyVcfHomozygosity
This software analyses a VCF file to determine homozygous regions of the genome. It can be run as follows:
```
java -Xmx5g SavvyVcfHomozygosity linkage_data input.vcf (options) (samplenames) >sample.bed
```
If a single samplename is used with no options, then the output will be a BED file containing homozygous regions for the sample. If multiple samplenames are used, then only regions that are homozygous in all samples are output. Adding the "-p" option alters this behaviour, so that regions where all the named samples share a haplotype are output.

The linkage_data file is used to filter the variants in the VCF file to only allow variants that are close to a Hardy-Weinberg equilibrium, to reduce the number of false heterozygous variants in homozygous regions. The variants are then analysed using a multi-resolution hidden markov model to find regions of at least 16 variants where the proportion of heterozygous variants is no more than 10%.
