#!/usr/bin/env perl
#
#
use warnings;
use strict;
use FAlite;


open IN, "blastdbcmd -db $ARGV[0] -line_length 1000 -entry all |";


my $fasta_file = new FAlite(\*IN); # or any other filehandle
while (my $entry = $fasta_file->nextEntry) {
	my $def = $entry->def;
	my @parts = split /\|/, $def;
	$parts[4] =~ s/.*=//g;
	print ">$parts[3]|kraken:taxid|$parts[4]\n";
	my $seq = $entry->seq;
	print "$seq\n";
}
