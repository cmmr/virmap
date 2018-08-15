#!/usr/bin/env perl
#
#
use warnings;
use strict;
use FAlite;

#add a pad of Ns to 5' and 3' ends of every FASTA entry

open IN, "$ARGV[0]";
my $padLength = $ARGV[1];
unless ($padLength =~ m/^\d+$/) {
	die "pad is not a number\n";
}
my $pad = "N" x $padLength;
my $fasta_file = new FAlite(\*IN); # or any other filehandle
while (my $entry = $fasta_file->nextEntry) {
	my $def = $entry->def;
	my $seq = $entry->seq;
	print "$def\n";
	$seq = $pad . $seq . $pad;
	print "$seq\n";
}
