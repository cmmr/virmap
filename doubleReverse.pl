#!/usr/bin/env perl


use warnings;
use strict;
use FAlite;

#doubles every sequence with its reverse complement


open IN, "$ARGV[0]";


my $fasta_file = new FAlite(\*IN); # or any other filehandle
while (my $entry = $fasta_file->nextEntry) {
	my $def = $entry->def;
	my $seq = $entry->seq;
	print "$def\n";
	print "$seq\n";
	my $revSeq = reverse($seq);
	$revSeq =~ tr/ACGT/TGCA/;
	$revSeq =~ tr/acgt/tgca/;	
	print "$def;rev=1\n";
	print "$revSeq\n";
}

