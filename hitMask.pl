#!/usr/bin/env perl
#
#
use warnings;
use strict;
use FAlite;

#mask fasta per base where an alignment event happened

my $in = "-";


open IN, "$in";
my $bitScores = {};



while (my $line = <IN>) {
	chomp $line;
	my @parts = split /\t/, $line;
	my $bitScore =  $parts[11];
	my $length = $parts[3];
	my $origDef = $parts[0];
	$origDef =~ s/;coords.*//g;
	my $bitPerLength = $bitScore / $length;
	if (not exists $bitScores->{$origDef}->{$parts[0]} or $bitScores->{$origDef}->{$parts[0]} < $bitPerLength) {
		$bitScores->{$origDef}->{$parts[0]} = $bitPerLength;
	}
}
close IN;
open IN, "$ARGV[0]";
my $fasta_file = new FAlite(\*IN); # or any other filehandle
my $seqs = {};
while (my $entry = $fasta_file->nextEntry) {
	my $def = $entry->def;
	$def =~ s/^>//g;
	$seqs->{$def} = $entry->seq;

}
close IN;
open IN, "$ARGV[1]";

$fasta_file = new FAlite(\*IN);
my $lengths = {};
while (my $entry = $fasta_file->nextEntry) {
	my $def = $entry->def;
	$def =~ s/^>//g;
	$lengths->{$def} = length($entry->seq);
}
foreach my $def (keys %$seqs) {
	unless (exists $bitScores->{$def}) {
		print ">$def\n";
		print "$seqs->{$def}\n";
		next;
	}
	my @letters = split //, $seqs->{$def};
	foreach my $hit (keys %{$bitScores->{$def}}) {
		$hit =~ m/coords=(\d+)/;
		my $start = $1;
		my $length = $lengths->{$hit};
		for (my $i = $start; $i < ($start + $length); $i++) {
			$letters[$i] = "N";
		}
	}
	my $outStr = join "", @letters;
	print ">$def\n";
	print "$outStr\n";
	
}
