#!/usr/bin/perl

use strict;
use Text::CSV;
use List::Util qw[min max];

##################
# INPUT DATA
my ( $input, $output, $panel, $panelname) = @ARGV;

# Extract genomic position of genes from the panel and put in "paneldata"
our @paneldata = ();

my $csv = Text::CSV->new;
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

open my $out2, ">", $output or die $!;
open my $in1, "<:encoding(UTF-8)", $input or die $!;
# print header
print $out2 join("\t","#Chr","Begin", "End", "${panelname}") . "\n";
#going throught variants in homozygosity regions to define them
while (<$in1>){
    my @listG =();
    chomp($_);
    my @info = split("\t", $_);
    if ($info[0] =~m/chr/){
        my $chrom=$info[0];
        my $begin=$info[1];
        my $end=$info[2];

        for (my $n =0; $n<scalar(@paneldata); $n++){
            if ($paneldata[$n][1] eq $chrom){
                if (($paneldata[$n][2]<$begin && $paneldata[$n][2]>$end) ||
                ($paneldata[$n][2]>$begin && $paneldata[$n][2]<$end) ||
                ($paneldata[$n][3]<$begin && $paneldata[$n][3]>$end) ||
                ($paneldata[$n][3]>$begin && $paneldata[$n][3]<$end)){push @listG, $paneldata[$n][0]}
            }
        }
        my $listG2 = join(";", @listG); 
        #in case there are no genes of interest in the homozygous region
        if($#listG <0){
            $listG2="none";
        }
        # printing in output file
        print $out2 $chrom . "\t" . $begin . "\t" . $end . "\t" . $listG2 . "\n";  
    }
}
close $out2;
close $in1;


exit;
