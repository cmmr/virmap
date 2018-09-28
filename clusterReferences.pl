#!/usr/bin/env perl


use warnings;
use strict;


my $in = $ARGV[0];
my $data = {};


my $prot2NucMap = {};
my $nuc2FullName = {};
my $nucHits = {};
open IN, "lbzip2 -d -c -n3 $ARGV[1] | cut -f1,3,4 |";

while (my $line = <IN>) {
	my @split = split /\t/, $line;
	$data->{$split[1]}->{$split[0]} = 1;
	$nucHits->{$split[1]}++;
	my @parts = split /\|/, $split[1], 5;
	$nuc2FullName->{$parts[3]} = $split[1];
	
	
}
close IN;
my $full2parts = {};
open IN, "lbzip2 -d -c -n3 $in | cut -f1,3,4 |";


while (my $line = <IN>) {
	my @splitParts = split /\t/, $line;
	my @lineparts = split /\|/, $splitParts[1];
	

	unless (exists $nuc2FullName->{$lineparts[2]}) {
		$lineparts[5] =~ s/;.*//g;
		$splitParts[1] =~ m/;taxId=(\d+)/;
		$nuc2FullName->{$lineparts[2]} = "gi|NoNucHits|gb|$lineparts[2]|$lineparts[6];taxId=$1";
		$nucHits->{$nuc2FullName->{$lineparts[2]}} = 0;
	}
	unless (exists $full2parts->{$nuc2FullName->{$lineparts[2]}}->{$splitParts[1]}) {
		$full2parts->{$nuc2FullName->{$lineparts[2]}}->{$splitParts[1]} = 1;
	}
	$data->{$nuc2FullName->{$lineparts[2]}}->{$splitParts[0]} = 1;
}
close IN;
my @references = sort {scalar(keys %{$data->{$b}}) <=> scalar(keys %{$data->{$a}}) || $nucHits->{$b} <=> $nucHits->{$a} || $a cmp $b} keys %{$data};
my %good;
my $clusterRadius = .97;
my $blackList = {};
my $chimeras = {};
my $whiteList = {};
my $e = 2.7182818284590452353602874713527;
my $gr = 1.6180339887498948482045868343656;
my $pi = 3.1415926535897932384626433832795;
my $ratio = $e;
foreach my $reference (@references) {
	my $nonInto = {};
	my $blackListCount = 0;
	my $count = 0;
	my $intersect = 0;
	foreach my $read (keys %{$data->{$reference}}) {
		if (exists $blackList->{$read}) {
			$blackListCount++;
			next;
		}
		$count++;
		if (exists $whiteList->{$read}) {
			$intersect++;
		} else {
			$nonInto->{$read} = 1;
		}
	}
	unless ($count) {
		next;
	}
	my $goodRefReadAmount = scalar(keys %$whiteList);
	unless ($goodRefReadAmount) {
		$goodRefReadAmount = $count;
	}
	my $adjustment = ((log($goodRefReadAmount/$count)/log($ratio)) / 100) * $ratio;
	my $intersectionPerc = ($intersect/$count);
	my $overallReadCount = scalar(keys %{$data->{$reference}});
	my $blackListIntersection = ($blackListCount/$overallReadCount);
	if ($blackListIntersection > $intersectionPerc) {
		$intersectionPerc = $blackListIntersection;
	}
	my $explainPerc = $intersect / $overallReadCount;
	my $radius = $clusterRadius - $adjustment;
	if ($intersectionPerc > $radius) {
		print STDERR "Absorbed because ($intersectionPerc > ($radius) [adjustment: $adjustment] or ($blackListIntersection) ($blackListCount/$overallReadCount) > ($radius))\n";
		print STDERR "Rejected: $reference\n$intersect explained by more than one previously accepted centroid; OverallReads: $overallReadCount; perc: $explainPerc\n\n";
		foreach my $read (keys %$nonInto) {
			$blackList->{$read} = 1;
		}
		next;
	} else {
		print STDERR "Not Absorbed because ($intersectionPerc > ($radius) [adjustment: $adjustment] or ($blackListIntersection) ($blackListCount/$overallReadCount) > ($radius))\n";
		print STDERR "Accepted: $reference\n$intersect explained by more than one previously accepted centroid; OverallReads: $overallReadCount; perc: $explainPerc\n\n";
		print STDERR "Added $reference as a new centroid\n";
		$good{$reference} = 1;
		foreach my $read (keys %{$data->{$reference}}) {
			$whiteList->{$read} = 1;
		}
	}

}
foreach my $reference (sort keys %good) {
	print "Parent:$reference\n";
	if (exists $full2parts->{$reference}) {
		foreach my $acc (sort keys %{$full2parts->{$reference}}) {
			print "Child:$acc\n";
		}	
	}
}	
