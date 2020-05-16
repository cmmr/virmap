#!/usr/bin/env perl
#
#
use warnings;
use strict;
use FAlite;


my $seen = {};

open IN, "-";
my $fasta_file = new FAlite(\*IN); # or any other filehandle
while (my $entry = $fasta_file->nextEntry) {
	my $def = $entry->def;
	my $pos = index($def, "\t");
	my $sum = substr($def, $pos + 1);
	if (exists $seen->{$sum}) {
		next;
	}
	else {
		my $seq = $entry->seq;
		my $name = substr($def, 0, $pos);
		$seen->{$sum} = 1;
		print "$name\n";
		print "$seq\n";
	}
}
