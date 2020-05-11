#!/usr/bin/env perl

use warnings;
use strict;


#convert fastq to virmap-style fasta

my $in = $ARGV[0];
my $countPrefix = $ARGV[1];

my $cmd = "cat";

my $fileType = `file -k $in`;

if ($fileType =~ m/bzip2/) {
	$cmd = "lbzip2 -d -c -n6";
} elsif ($fileType =~ m/gzip/) {
	$cmd = "pigz -d -c";
}

open R1, "$cmd $in |";
while (my $head1 = <R1>) {
	chomp $head1;
	my $seq1 = <R1>;
	if (length($seq1) < 50) {
		next;
	}
	my $space1 = <R1>;
	my $qual1 = <R1>;
	unless ($seq1 and $head1) {
		last;
	}
	$head1 =~ s/^@/>$countPrefix./;
	$head1 =~ s/\s/./g;
	print "$head1\n$seq1";

}
