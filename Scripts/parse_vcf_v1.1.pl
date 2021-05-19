#!/usr/bin/perl
#!/bin/bash
use strict;                # Terminate when error occurs
use warnings;              # Display all warning messages

##### FUNCTIONS FROM: Copyleft, T. J. Finney, 2002-11-17. Free for non-commercial use.
##### Available at: https://www.halotype.com/RKM/figures/TJF/binomial.txt

sub binomial_approx {
    my $k = shift(@_);
    my $n = shift(@_);
    my $p = shift(@_);
    my($sigma, $mu, $pi, $const, $exponent, $prob);

    $mu = $n * $p;
    $sigma = sqrt($mu * (1 - $p));
    $pi = atan2(1, 1) * 4;
    $const = 1 / ($sigma * sqrt(2 * $pi));
    $exponent = -0.5 * (($k - $mu) / $sigma)**2;
    $prob = $const * exp($exponent);

    return $prob;
}

sub binomial {
    my $k = shift(@_);
    my $n = shift(@_);
    my $p = shift(@_);
    my $prob;

    $prob = ($p**$k) * ((1 - $p)**($n - $k)) * &factorial($n) / (&factorial($k) * &factorial($n - $k));

    return $prob;
}

sub factorial {
    my $n = shift(@_);
    my $fact = 1;

    if (($n < 0) or (170 < $n)) {
	die "Factorial out of range";
    }

    for(my $i = 1; $i <= $n; $i++) {
	$fact *= $i;
    }

    return $fact;
}

###################
### EXPLANATION ###
###################
### script that cleans the input file into a more readable file (arrange the format, delete some doublet column, ...)


##################
# get "input file" and "output file" from command line
my ( $input, $output) = @ARGV;

#####################################
### Optional but recommended
#####################################

if (!(-f $input)){	die "Invalid input file";	}
if ( !defined($output) || $output eq ''){	die "Missing output file";	}

#####################################
### END OF OPTIONAL PART
#####################################

open( IN, $input ) or die "Unable to open file $input :$!"; # open the input file (in reading mode)
open( OUT, ">".$output ) or die "Unable to open file $output :$!"; # open (or create) the output file in writing mode

while (<IN>){
	if(/^#CHROM/){
		chomp;
		my @record = split(/\t/); #list of elements in the record
		print OUT "$record[0]\t$record[1]\t$record[2]\t$record[3]\t$record[4]\tZYG\t$record[5]\t$record[6]\t";
		print OUT "read_depth\tread_ref_alt\t%_alt\tbinomial_prob\n";
	}
	if (/^\w/){
		chomp;
		my @record = split(/\t/); #list of elements in the record

		print OUT "$record[0]\t$record[1]\t$record[2]\t$record[3]\t$record[4]\t";

		my @cov_info=split(/:/,$record[9]);
		my $zyg="";
		if($cov_info[0] eq "1/1"){
			$zyg="hom";	
		}
		if($cov_info[0] eq "0/1"){
			$zyg="het";	
		}
		print OUT "$zyg\t$record[5]\t$record[6]\t";

		my $read_ref_alt="0,0";
		my $done = 0;
		my $read_alt = 0;
		my $read_ref = 0;
		my $read_tot = 0;
		my $i=0;

		for my $number (split(/;/,$record[7])){ # split by transcript				
			if($number =~ /^DP4/){
				$read_ref_alt=(split(/=/,$number))[1];
				$done=1;
				$read_alt=(split(/,/,$read_ref_alt))[2]+(split(/,/,$read_ref_alt))[3];
				$read_ref=(split(/,/,$read_ref_alt))[0]+(split(/,/,$read_ref_alt))[1];
			}
		}

		if ($done == 0){
			for my $number (split(/:/,$record[8])){ # split by transcript				
				if($number eq "AD"){
					$read_ref_alt=(split(/:/,$record[9]))[$i];
					if($read_ref_alt ne "."){
						$read_alt=(split(/,/,$read_ref_alt))[1];
						$read_ref=(split(/,/,$read_ref_alt))[0];
						$done=1;
					}
				}
				$i=$i+1;
			}
		}

		$i=0;

		if ($done == 0){
			for my $number (split(/:/,$record[8])){ # split by transcript		
				if($number eq "DP"){
					$read_tot=(split(/:/,$record[9]))[$i];
				}
				if($number eq "AO"){
					$read_alt=(split(/:/,$record[9]))[$i];
				}
				$i=$i+1;
			}
			$read_ref=$read_tot-$read_alt;
		}

		my $read_depth=$read_alt+$read_ref;
		my $perc_alt=0;
		if($read_alt > 0){
			$perc_alt=sprintf "%.2f", $read_alt/$read_depth;
		}

		my $n = $read_depth;
		my $k = $read_alt;
		my $m = $read_depth-$read_alt; # ref reads
		my $p = 0.5;
		my $proba = 0;
		my $prob = 0;

		if ($n > 100) {
		    if ($m >= $k){
		    	for (my $i = 0; $i <= $k; $i++) {
					$prob = &binomial_approx($i, $n, $p);
					$proba += $prob;
			    }
		    }
		    else{
		    	for (my $i = $k; $i <= $n; $i++) {
					$prob = &binomial_approx($i, $n, $p);
					$proba += $prob;
			    }
		    }
		}
		else {
			if ($m >= $k){
		    	for (my $i = 0; $i <= $k; $i++) {
					$prob = &binomial($i, $n, $p);
					$proba += $prob;
			    }
		    }
		    else{
		    	for (my $i = $k; $i <= $n; $i++) {
					$prob = &binomial($i, $n, $p);
					$proba += $prob;
			    }
		    }
		}

		print OUT "$read_depth\t$read_ref,$read_alt\t$perc_alt\t$proba\n";
	}
}
close(IN);
close(OUT);

exit;
