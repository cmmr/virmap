#!/usr/bin/env perl
#
#
use warnings;

use strict;
use FAlite;

#find contigs that have big kmer pool overlaps, and small length differences
#naive attempt to remove multiple nearly identical circular contigs with different origin offsets

my $kmers = {};
my $seqs = {};
my $sizes = {};
my $alignments = {};
my $lengths = {};
my $count = 0;
my $kmerLength = 101;
my $radius = .99;
my $self = {};
open IN, "$ARGV[0]";
my $segments = {};
my $fasta_file = new FAlite(\*IN); 
while (my $entry = $fasta_file->nextEntry) {
	my $def = $entry->def;
	$def =~ s/;coords=\d+//g;
	$def =~ s/^>//g;
	$segments->{$def}++;
	my $seq = $entry->seq;
	my $length = length($seq);
	if ($length < $kmerLength) {
		next;
	}
	$seq = uc($seq);
	for (my $i = 0; $i < ($length - $kmerLength); $i += 3) {
		my $kmer = substr($seq, $i, $kmerLength);
		my $revKmer = reverse($kmer);
		$revKmer =~ tr/ACGT/TGCA/;
		$kmers->{$kmer} = $revKmer;
		$kmers->{$revKmer} = $kmer;
		$self->{$def}++;
	}
}
close IN;
my $orig = {};
open IN, "$ARGV[1]";
$fasta_file = new FAlite(\*IN);
while (my $entry = $fasta_file->nextEntry) {
	my $def = $entry->def;
	my $seq = $entry->seq;
	$def =~ s/^>//g;
	unless (exists $self->{$def} and $self->{$def} > .98 * $segments->{$def}) {
		print STDERR "skipped $def\n";
		print ">$def\n";
		print "$seq\n";
		if (exists $self->{$def}) {
			delete $self->{$def};
		}
	} else {
		$orig->{$def} = $seq;
	}
	
}

print STDERR "Clustering\n";
my $acceptedKmers = {};
my $clusters = {};
foreach my $def (sort {$self->{$b} <=> $self->{$a}} keys %$self) {
	my $seq = $orig->{$def};
	$seq = uc($seq);
	my $overlap = 0;
	my $overallHits = 0;
	my $seenkmers = {};
	my $denomOverlap = 0;
	for (my $j = 0; $j < length($seq) - $kmerLength; $j++) {
		my $kmer = substr($seq, $j, $kmerLength);
		if (exists $kmers->{$kmer}) {
			$seenkmers->{$kmer} = 1;
			$seenkmers->{$kmers->{$kmer}} = 1;
			$overallHits++
		}
		if (exists $acceptedKmers->{$kmer}) {
			$overlap++;
		}
	}
	my $avgOverlap = 0;
	my $randomKmerOverlap = 0;
	if ($overallHits) {
		$randomKmerOverlap = $overlap / $overallHits;
	}
	
	if ($randomKmerOverlap <= $radius) {
		$clusters->{$def}->{$def} = $seq;
		foreach my $kmer (keys %$seenkmers) {
			$acceptedKmers->{$kmer}->{$def} = 1;
		}
		
	} else {
		my $hits = {};
		foreach my $kmer (keys %$seenkmers) {
			foreach my $other (keys %{$acceptedKmers->{$kmer}}) {
				if ($other eq $def) {
					next;
				}
				$hits->{$other}++;
			}
		}
		my $topHit;
		my @hits = sort {$hits->{$b} <=> $hits->{$a} || $self->{$b} <=> $self->{$a}} keys %$hits;
		if (scalar(@hits)) {
			$topHit = $hits[0];
		}
		if ($topHit and (abs(($self->{$topHit}/$self->{$def}) - 1) < .1 and abs(($self->{$def}/$self->{$topHit}) - 1) < .1 )) {
			$clusters->{$topHit}->{$def} = $seq;
			print STDERR "Added $def to cluster $hits[0] with $hits->{$hits[0]} common kmers\n";
		} else {
			$clusters->{$def}->{$def} = $seq;
			foreach my $kmer (keys %$seenkmers) {
				$acceptedKmers->{$kmer}->{$def} = 1;
			}	
		}
	}

}
my $clusterCount = 0;
foreach my $cluster (keys %$clusters) {
	my $def = $cluster;
	unless ($def =~ m/^contig/) {
		print ">$def\n";
	} else {
		my $backup = "";
		foreach my $check (sort keys %{$clusters->{$cluster}}) {
			unless ($check =~ m/^contig/) {
				$backup = $check;
			}
		}
		if ($backup) {
			print ">$backup\n";
		} else {
			print ">$def\n";
		}
	}
	print "$orig->{$def}\n";
}

