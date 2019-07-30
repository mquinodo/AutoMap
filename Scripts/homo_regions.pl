#!/usr/bin/perl

use strict;

##################
# INPUT DATA
my ( $input, $output, $panel, $panelname, $window, $windowthres, $trimming, $maxgap, $extend, $maxsize) = @ARGV;

# look if there is a gene panel provided
my $pan=1;
if ( !defined($panelname) || $panelname eq 'NA'){$pan=0;}

# Extract genomic position of genes from the panel and put in "paneldata"
our @paneldata = ();
if ($pan==1){
    print " * Import the genes information\n";
    open my $fj, '<', $panel or die "Could not open $panel: $!";
    my $j=0;
    while(<$fj>) { 
        chomp($_);
        my @row = split("\t", $_);
        for (my $i=0; $i<scalar(@row);$i++){
            $paneldata[$j][$i] = $row[$i];
        }
        $j++;
    }
    close $fj;
}

# Extract variants from input file and put in "inputdata"
our @inputdata = (); 
open my $fh, '<', $input or die "Could not open $input: $!";
my $j=0;
while(<$fh>) { 
    chomp($_);
    my @row = split("\t", $_);
    for (my $i=0; $i<scalar(@row);$i++){
        $inputdata[$j][$i] = $row[$i];
    }
    $j++;
}
close $fh;

#opening input file and creating output file
open my $in, "<", $input or die $!;
my $outfile = join("", $output,".homozygosity");
open my $out, ">", $outfile or die $!;

print " * Treatment of the data\n";
#going throught data to extract variants in homozygous regions
my $k = 0; # iterations in the loop
while (<$in>){
    my $score = 0;
    my $score_max = 0;
    my $homoz = "No";
    my $no = 0;
    chomp($_);
    my @info = split("\t", $_); # info contain all the information of the variant
    if ($info[0] =~ /#/){
        print $out join("\t",@info) . "\t" . "Homo_region" . "\n"; #adding title in header
    }
    else {  
        my @array3 = sort { $a <=> $b } ($k-($window-1),0);
        my @array4 = sort { $a <=> $b } ($k+($window-1),$j);

        my @zygosity = map $_->[ 5 ], @inputdata[$array3[-1]..$array4[0]]; # taking het/hom info for variants $window before to $window after the one analyzed
        my @chromosome = map $_->[ 0 ], @inputdata[$array3[-1]..$array4[0]]; # same but for chromosome

        # in score we put number of these $window variants which are homozygous and on same chromosome
        for (my $z = 0; $z<$window*2; $z++){
            if ($zygosity[$z] eq "hom" && $chromosome[$z] eq $info[0]){
                $score++;
            }
        }
        # if less than $windowthres variants before and after the one analyzed are homozygous we will stop analysis for this one
        if ($score<=$windowthres-1){$no=1;} 
        # taking a $window-variant window , sliding it around the variant and taking best number of homozygous variant in score_max
        for (my $p =0; $p<$window; $p++){

            my @array1 = sort { $a <=> $b } ($k-($window-1)+$p,0);
            my @array2 = sort { $a <=> $b } ($k+$p,$j);

            my @zygosity = map $_->[ 5 ], @inputdata[$array1[-1]..$array2[0]];
            my @chromosome = map $_->[ 0 ], @inputdata[$array1[-1]..$array2[0]];

            $score = 0;
            for (my $z = 0; $z<$window; $z++){
                if ($zygosity[$z] eq "hom" && $chromosome[$z] eq $info[0]){
                    $score++;
                }
            }
            if ($score > $score_max){$score_max = $score}
            if ($no == 1){$p=80} # stop if previous criteria is wrong to go faster
        }
        # printing Yes in output if the highest score if score bigger than threshold
        if($score_max>=$windowthres){$homoz = "Yes"}
        print $out join("\t",@info) . "\t" . $homoz . "\n";
    }
    $k++;
}
close $in;
close $out;

# trimming of ends of regions
my $command = join("", "bash ", $trimming, " ", $outfile," ",$outfile,".trim");
system($command);

# extension of the ends of the regions
my $command = join("", "bash ", $extend, " ", $outfile,".trim ",$outfile,".trim.extend ", $maxsize);
system($command);

######################## 
print " * Printing of the homozygous regions\n";
my $outfile = join("", $output);
open my $out3, ">", join("", $outfile,".tsv") or die $!;

my $state = "";
my $nb_var = 0;
my $nb_hom = 0;
my $nb_het = 0;
my $chromosome = "chr1";
my $last_position = 0;
my $begin =0;
my $chromo = "";
my @listG =();

if ($pan==1){
    open my $out2, ">", join("", $outfile,".",${panelname},".tsv") or die $!;
    open my $in1, "<:encoding(UTF-8)", join("", $output,".homozygosity.trim.extend") or die $!;
    # print header
    print $out2 join("\t","#Chr","Begin", "End", "Size(Mb)","Nb_variants", "Percentage_homozygosity", "${panelname}") . "\n";
    #going throught variants in homozygosity regions to define them
    while (<$in1>){
        chomp($_);
        my @info = split("\t", $_);
        #$info[-1] is the status (with Yes or No)
        # if we were in a homo. region but now the next variant is not inside or if we change chromosome
        if (($info[-1] eq "No" && $state eq "Yes") || ($chromosome ne $info[0] && $info[-1] eq "Yes" && $state eq "Yes") || $info[5] eq "ext"){
            my $percent = sprintf "%.2f", $nb_hom/($nb_hom+$nb_het+0.000000000001)*100;
            # fetching genes which are in the terminated region
            for (my $n =0; $n<scalar(@paneldata); $n++){
                if ($paneldata[$n][1] eq $chromo){
                    if (($paneldata[$n][2]<$begin && $paneldata[$n][2]>$last_position) ||
                    ($paneldata[$n][2]>$begin && $paneldata[$n][2]<$last_position) ||
                    ($paneldata[$n][3]<$begin && $paneldata[$n][3]>$last_position) ||
                    ($paneldata[$n][3]>$begin && $paneldata[$n][3]<$last_position)){push @listG, $paneldata[$n][0]}
                }
            }
            my $listG2 = join(";", @listG); 
            #in case there are no genes of interest in the homozygous region
            if($#listG <0){
                $listG2="none";
            }
            if ($info[5] eq "ext" && $chromo eq $info[0]){
                ${last_position}=$info[1];
            }
            my $size=sprintf "%.2f",(${last_position}-$begin)/1000000;#size of the homozygous region in Mb
            # printing in output file
            print $out2 $last_position . "\t" . $size . "\t" . $nb_var . "\t" . $percent . "\t" . $listG2 . "\n";
            # reset values to default
            $state = "";
            $nb_var = 0;
            $nb_hom = 0;
            $nb_het = 0;
            $begin = 0;
            $chromo = "";
            $size=0;
            @listG = ();
        }
        # if we begin a new region
        if ($info[-1] eq "Yes" && $state eq ""){
            print $out2 $info[0] . "\t" . $info[1] . "\t"; # print chr and begin in output file
            $state = "Yes";
            $begin = $info[1];
            $chromo = $info[0];
            if ($info[5] eq "hom"){$nb_hom++; $nb_var++;}
            if ($info[5] eq "het"){$nb_het++; $nb_var++;}
        }
        # if we are staying in the homozygous region
        if ($info[-1] eq "Yes" && $state eq "Yes" && $info[1]<(${last_position}+1000000*$maxgap)){     
            if ($info[5] eq "hom"){$nb_hom++; $nb_var++;}
            if ($info[5] eq "het"){$nb_het++; $nb_var++;}
        }
        # if we are staying in the homozygous region but the gap is higher than accepted by --maxgap option
        if ($info[-1] eq "Yes" && $state eq "Yes" && $info[1]>=(${last_position}+1000000*$maxgap)){
            if ($info[5] eq "hom"){$nb_hom++}
            if ($info[5] eq "het"){$nb_het++}
            my $percent = sprintf "%.2f", $nb_hom/($nb_hom+$nb_het+0.000000000001)*100;
            # fetching genes which are in the terminated region
            for (my $n =0; $n<scalar(@paneldata); $n++){
                if ($paneldata[$n][1] eq $chromo){
                    if (($paneldata[$n][2]<$begin && $paneldata[$n][2]>$last_position) ||
                    ($paneldata[$n][2]>$begin && $paneldata[$n][2]<$last_position) ||
                    ($paneldata[$n][3]<$begin && $paneldata[$n][3]>$last_position) ||
                    ($paneldata[$n][3]>$begin && $paneldata[$n][3]<$last_position)){push @listG, $paneldata[$n][0]}
                }
            }
            my $listG2 = join(";", @listG); 
            #in case there are no genes of interest in the homozygous region
            if($#listG <0){
                $listG2="none";
            }
            my $size=sprintf "%.2f",(${last_position}-$begin)/1000000;#size of the homozygous region in Mb
            # printing in output file
            print $out2 $last_position . "\t" . $size . "\t" . $nb_var . "\t" . $percent . "\t" . $listG2 . "\n";
            # reset values to default
            $state = "";
            $nb_var = 0;
            $nb_hom = 0;
            $nb_het = 0;
            $begin = 0;
            $chromo = "";
            $size=0;
            @listG = ();
            print $out2 $info[0] . "\t" . $info[1] . "\t"; # print chr and begin in output file
            $state = "Yes";
            $begin = $info[1];
            $chromo = $info[0];
            if ($info[5] eq "hom"){$nb_hom++; $nb_var++;}
            if ($info[5] eq "het"){$nb_het++; $nb_var++;}
        }
        $chromosome = $info[0];
        $last_position = $info[1];
    }
    close $out2;
    close $in1;
}

my $state = "";
my $nb_var = 0;
my $nb_hom = 0;
my $nb_het = 0;
my $chromosome = "chr1";
my $last_position = 0;
my $begin =0;
my $chromo = "";
open my $in, "<", join("", $output,".homozygosity.trim.extend") or die $!;

# print header
print $out3 join("\t","#Chr","Begin", "End", "Size(Mb)","Nb_variants", "Percentage_homozygosity") . "\n";

#going throught variants in homozygosity regions to define them
while (<$in>){
    chomp($_);
    my @info = split("\t", $_);
    #$info[-1] is the status (with Yes or No)
    # if we were in a homo. region but now the next variant is not inside or if we change chromosome
    if (($info[-1] eq "No" && $state eq "Yes") || ($chromosome ne $info[0] && $info[-1] eq "Yes" && $state eq "Yes") || $info[5] eq "ext"){
        my $percent = sprintf "%.2f", $nb_hom/($nb_hom+$nb_het+0.000000000001)*100;
        if ($info[5] eq "ext" && $chromo eq $info[0]){
            ${last_position}=$info[1];
        }
        my $size=sprintf "%.2f",(${last_position}-$begin)/1000000;#size of the homozygous region in Mb
        # printing in output file
        print $out3 $last_position . "\t" . $size . "\t" . $nb_var . "\t" . $percent . "\n";
        # reset values to default
        $state = "";
        $nb_var = 0;
        $nb_hom = 0;
        $nb_het = 0;
        $begin = 0;
        $chromo = "";
        $size=0;
    }
    # if we begin a new region
    if ($info[-1] eq "Yes" && $state eq ""){
        print $out3 $info[0] . "\t" . $info[1] . "\t"; # print chr and begin in output file
        $state = "Yes";
        $begin = $info[1];
        $chromo = $info[0];
        if ($info[5] eq "hom"){$nb_hom++; $nb_var++;}
        if ($info[5] eq "het"){$nb_het++; $nb_var++;}
    }
    # if we are staying in the homozygous region
    if ($info[-1] eq "Yes" && $state eq "Yes" && $info[1]<(${last_position}+1000000*$maxgap)){
        if ($info[5] eq "hom"){$nb_hom++; $nb_var++;}
        if ($info[5] eq "het"){$nb_het++; $nb_var++;}
    }
    # if we are staying in the homozygous region but the gap is higher than accepted by --maxgap option
    if ($info[-1] eq "Yes" && $state eq "Yes" && $info[1]>=(${last_position}+1000000*$maxgap)){
        if ($info[5] eq "hom"){$nb_hom++}
        if ($info[5] eq "het"){$nb_het++}
        my $percent = sprintf "%.2f", $nb_hom/($nb_hom+$nb_het+0.000000000001)*100;
        my $size=sprintf "%.2f",(${last_position}-$begin)/1000000;#size of the homozygous region in Mb
        # printing in output file
        print $out3 $last_position . "\t" . $size . "\t" . $nb_var . "\t" . $percent . "\n";
        # reset values to default
        $state = "";
        $nb_var = 0;
        $nb_hom = 0;
        $nb_het = 0;
        $begin = 0;
        $chromo = "";
        $size=0;
        print $out3 $info[0] . "\t" . $info[1] . "\t"; # print chr and begin in output file
        $state = "Yes";
        $begin = $info[1];
        $chromo = $info[0];
        if ($info[5] eq "hom"){$nb_hom++; $nb_var++;}
        if ($info[5] eq "het"){$nb_het++; $nb_var++;}
    }
    $chromosome = $info[0];
    $last_position = $info[1];
}

# remove temporary file
my $command = join("", "rm ", $output, ".homozygosity");
#system($command);
my $command = join("", "rm ", $output, ".homozygosity.trim");
#system($command);
my $command = join("", "rm ", $output, ".homozygosity.trim.extend");
#system($command);

close $out3;
close $in;
exit;