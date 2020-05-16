#!/usr/bin/env perl
#
#
use warnings;
use strict;
use Sereal;
use Compress::Zstd qw(decompress);;



open IN, "-";
open PRE, "$ARGV[0]";
my $serealData;
binmode(PRE);
{
	local $/;
	undef $/;
	$serealData = <PRE>;
}
my $decoder = Sereal::Decoder->new();
my $seen = $decoder->decode(decompress($serealData));
close PRE;
open OUT, ">$ARGV[1]";
open UNIQ, "| shuf >$ARGV[2]";
while (my $line = <IN>) {
	my @parts = split /\t/, $line;
	if (exists $seen->{$parts[1]}) {
		print OUT "$line";
		next;
	}
	else {
		print UNIQ "$parts[0]\t$parts[2]";
	}
}
close OUT;
close UNIQ;
