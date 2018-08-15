#!/usr/bin/env perl

use warnings;
use strict;
use FAlite;


#split scaffolds into N free segments (not n free segments)

open IN, "$ARGV[0]";
my $defaultChunk = "inf";
if ($ARGV[1]) {
	$defaultChunk = $ARGV[1];;
}
if ($defaultChunk < 100 and $defaultChunk != 0) {
	$defaultChunk = 100;
}

my $fasta_file = new FAlite(\*IN); # or any other filehandle
my $count = 0;
while (my $entry = $fasta_file->nextEntry) {
	my $def = $entry->def;
	my $seq = $entry->seq;
	my $length = length($seq);
	my @letters = split //, $seq;
	my $on = 0;
	my $startCoord = 0;
	my @localChain;
	for (my $i = 0; $i < scalar(@letters); $i++) {
		if ($letters[$i] eq "N" and $on) {
			$on = 0;
			printNow(\@localChain, $startCoord, $def, $defaultChunk);
			@localChain = ();
		} elsif (uc($letters[$i]) ne "N" and not $on) {
			$startCoord = $i; 
			$on = 1;
			push @localChain, $letters[$i];
		} elsif ($letters[$i] ne "N" and $on) {
			push @localChain, $letters[$i];
		}
	}
	if (scalar(@localChain)) {
		printNow(\@localChain, $startCoord, $def, $defaultChunk);
	}
	$count++;
}
sub printNow {
	
	my @localChain = @{$_[0]};
	my $startCoord = $_[1];
	my $def = $_[2];
	my $chunk = $_[3];
	my $outStr = join "", @localChain;
	my $length = length($outStr);
	if ($length < $chunk) {
		$chunk = $length;
	}
	if ($length < 50) {
		return;
	}
	my $denom = int($length / $chunk);
	if ($denom == 0) {
		#shouldn't happen
		$denom = 1;
	}
	if (scalar(@localChain) < 17) {
		return;
	}
	for (my $j = 0; $j < $denom; $j++) {
		my $start = $j * $chunk;
		my $end = $start + $chunk;
		if ($j == ($denom - 1)) {
			$end = $length;
		}
		my $outCoords = $startCoord + $start;
		my $localStr = uc(substr($outStr, $start, ($end - $start)));
		print "$def;coords=$outCoords\n";
		print "$localStr\n";    
		$start = $start + int($chunk / 2);
		$end = $end + int($chunk / 2);
		if ($end > $length) {
			next;
		}
		$outCoords = $startCoord + $start;
		$localStr = uc(substr($outStr, $start, ($end - $start)));
		print "$def;coords=$outCoords\n";
		print "$localStr\n";
	}
	



}
