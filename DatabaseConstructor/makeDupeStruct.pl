#!/usr/bin/env perl
#
#
use warnings;
use strict;
use Sereal;
use Compress::Zstd qw(compress);;


open IN, "-";
my $seen = {};
while (my $line = <IN>) {
	chomp $line;
	$seen->{$line} = 1;
}
my $encoder = Sereal::Encoder->new();
my $encoded = $encoder->encode($seen);
my $compressOut = compress($encoded);
print "$compressOut";
