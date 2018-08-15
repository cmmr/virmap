#!/usr/bin/env perl
#
#
#
use warnings;
use strict;

#converts all letters in a FASTA to upper case

open IN, "-";
my $all = 0;
if (defined($ARGV[0]) and $ARGV[0] eq "all") {
	$all = 1;
}

while (my $line = <IN>) {
	if ($line =~ m/^>/) {
		print "$line";
	} else {
		if ($all) {
			$line = uc($line);
		} else {
			$line =~ tr/acgt/ACGT/;;
		}
		print "$line";
	}

}
