#!/usr/bin/env perl
#
#
use warnings;
use strict;
#use DataBrowser qw(browse);
use Sereal;
use Compress::Zstd qw(compress decompress compress_mt);


my $decoder = Sereal::Decoder->new();
my $in = $ARGV[0];
open NODE, "$in";
my $serealData;
binmode(NODE);
{
	local $/;
	undef $/;
	$serealData = <NODE>;
}
close NODE;
my $data = $decoder->decode(decompress($serealData));
#browse($data);
