#!/usr/bin/env perl

use warnings;
use strict;
use FAlite;

#convert fasta input to virmap-style fasta 

my $in = $ARGV[0];
my $countPrefix = $ARGV[1];

my $cmd = "cat";

my $fileType = `file -k $in`;

if ($fileType =~ m/bzip2/) {
	$cmd = "lbzip2 -d -c -n4";
} elsif ($fileType =~ m/gzip/) {
	$cmd = "pigz -d -c -p4";
}


open IN, "$cmd $in |";
my $fasta_file = new FAlite(\*IN); # or any other filehandle
while(my $entry = $fasta_file->nextEntry) {
	my $head = $entry->def;
	my $seq = $entry->seq;
	if (length($seq) < 50) {
		next;
	}
	$head =~ s/^>/>$countPrefix./g;
	$head =~ s/\s/./g;
	unless ($head and $seq) {
		last;
	}
	print "$head\n$seq\n";

}
