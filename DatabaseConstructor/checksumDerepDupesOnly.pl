#!/usr/bin/env perl
#
#


use warnings;
use strict;


open IN, "-";

my $seen = {};

while (my $line = <IN>) {
	my @parts = split /\t/, $line;
	if (exists $seen->{$parts[1]}) {
		next;
	} else {
		$seen->{$parts[1]} = 1;
		print "$parts[0]\t$parts[2]";
	}
}
