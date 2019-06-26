#!/usr/bin/env perl
#
#
use warnings;
use strict;
use Sereal;
use Compress::Zstd qw(compress decompress compress_mt);
use FAlite;






my $decoder = Sereal::Decoder->new();

my ($parents, $names, $ranks, $children) = makeTaxonomy($ARGV[2]);
my $groups = {
        "131567" => "cellularOrganisms",
        "12908" => "unclassified",
        "12884" => "viroids",
        "28384" => "others",
        "10239" => "virus"
};


open IN, "lbzcat -n4 $ARGV[1] |";

my $filter = {};

while (my $line = <IN>) {
	chomp $line;
	my @parts = split /\t/, $line;
	my $contig = $parts[0];
	my $taxaId = $parts[1];
	if ($taxaId <= 1) {
		next;
	}
	my $group = traceToGroup($groups, $parents, $taxaId);
	if ($group == 10239 or $group == 12284) {
		next;
	}
	$filter->{$contig} = 1;
}
close IN;

open IN, "$ARGV[0]";
my $fasta_file = new FAlite(\*IN); # or any other filehandle
while (my $entry = $fasta_file->nextEntry) {
        my $head = $entry->def;
        $head =~ s/>//g;
	if (exists $filter->{$head}) {
		print STDERR "KRAKEN FILTERED: $head\n";
		next;
	}
	
	my $seq = $entry->seq;
	print ">$head\n";
	print "$seq\n";
}











sub makeTaxonomy {
	my $in = $_[0];
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

	my $parents = $data->{'parents'};
	my $ranks = $data->{'ranks'};
	my $names = $data->{'names'};
	my $children = $data->{'children'};
	return($parents, $names, $ranks, $children);

}

sub traceToGroup {
	my $groups = $_[0];
	my $parents = $_[1];
	my $node = $_[2];
	my $group = 0;

	my $escape = 0;
	while (1) {
		$escape++;
		if ($escape == 100) {
			last;
		}
		if (exists $groups->{$node}) {
			$group = $node;
			last
		} elsif (exists $parents->{$node}) {
			$node = $parents->{$node};
		} else {
			last;
		}
	}
	return($group);

}

