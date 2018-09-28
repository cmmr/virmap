#!/usr/bin/env perl

use warnings;
use strict;
use FAlite;
use File::Temp qw/tempdir/;

my $in = $ARGV[0];
my $tmpdir = tempdir( CLEANUP => 1 );
my $basename = `basename $ARGV[0]`;
chomp $basename;
my $iterateCount = 0;
my $THREADS = $ARGV[1];
my $outStrict = "";
if ($ARGV[2]) {
	$outStrict = $ARGV[2];
}
my $outCount = 0;

while (1) {
	my $outBefore = {};
	my $outAfter = {};
	print STDERR "Outer iteration $outCount.$iterateCount of merge\n";
	open IN, "$in";
	my $fasta_file = new FAlite(\*IN);
	while (my $entry = $fasta_file->nextEntry) {
		my $def = $entry->def;
		$outBefore->{$def} = 1;
	}
	close IN;
	my $strict = "";
	if ($outCount == 0 and $outStrict eq "init") {
		$strict = "init";
	} 
	if ($outStrict eq "init" and $outCount > 0) {
		$strict = "strict";
	}
	if ($outStrict eq "strict") {
		$strict = "strict";
	}	
	while (1) {
		my $before = {};
		my $after = {};
		if ($outStrict eq "strict" and $outCount > 0) {
			system("cat $in > $tmpdir/$basename.final.fa");
			$in = "$tmpdir/$basename.final.fa";
		}
		print STDERR "Iteration $outCount.$iterateCount of merge\n";
		open IN, "$in";
		$fasta_file = new FAlite(\*IN); # or any other filehandle
		while (my $entry = $fasta_file->nextEntry) {
			my $def = $entry->def;
			$before->{$def} = 1;
		}
	
		my $oldIn = $in;
		close IN;
		system("doubleReverse.pl $in > $tmpdir/$basename.prepped.iterationMerge.$outCount.$iterateCount.fa");
		print STDERR "Called merge with $strict\n";
		system("easyMerge.pl $tmpdir/$basename.prepped.iterationMerge.$outCount.$iterateCount.fa $THREADS $strict > $tmpdir/$basename.iterationMerge.$outCount.$iterateCount.fa");
		$in = "$tmpdir/$basename.iterationMerge.$outCount.$iterateCount.fa";
		open IN, "$in";
		$fasta_file = new FAlite(\*IN);
		while (my $entry = $fasta_file->nextEntry) {
			my $def = $entry->def;
			$after->{$def} = 1;
		}
		close IN;
		print STDERR "END Inner iteration $outCount.$iterateCount of merge\n";
		my $beforeCount = scalar(keys %$before);
		my $afterCount = scalar(keys %$after);
		print STDERR "inner before after merge: $beforeCount\n";
		print STDERR "inner after after mege count: $afterCount\n";
		if ($beforeCount  == $afterCount) {
			print STDERR "Finished inner loop after $outCount.$iterateCount iterations\n";
			print STDERR "in was $in\n";
			last;
		} elsif ($afterCount > $beforeCount) {
			$iterateCount++;
			print STDERR "ERROR: After cannot be greater than before\n";
			foreach my $def (keys %$after) {
				$def =~ s/;noChange=1//g;
				unless (exists $before->{$def}) {
					print STDERR "$def is not in before\n";
				}
				
			}
			die;
			
		}
		$iterateCount++;
	}
	print STDERR "\n";
	print STDERR "BACK to outer loop\n";
	open IN, "$in";
	$fasta_file = new FAlite(\*IN);
	while (my $entry = $fasta_file->nextEntry) {
		my $def = $entry->def;
		$outAfter->{$def} = 1;
	}
	close IN;
	print STDERR "END Outer iteration $outCount.$iterateCount of merge\n";
	my $outBeforeCount = scalar(keys %$outBefore);
	my $outAfterCount = scalar(keys %$outAfter);
	print STDERR "before count $outBeforeCount\n";
	print STDERR "out after count $outAfterCount\n";
	if (scalar(keys %$outBefore) == scalar(keys %$outAfter)) {
		print STDERR "Finished merging after $iterateCount iterations\n";
		print STDERR "printing infile $in\n";
		open IN, "$in";
		$fasta_file = new FAlite(\*IN);
		while (my $entry = $fasta_file->nextEntry) {
			my $def = $entry->def;
			my $seq = $entry->seq;
			$def =~ s/;noChange=1//g;
			print "$def\n";
			print "$seq\n";
		}		

		last;
	} else {
		$outCount++;
		open IN, "$in";
		$fasta_file = new FAlite(\*IN);
		open OUT, ">$tmpdir/$basename.iterationMerge.$outCount.$iterateCount.fa";
		while (my $entry = $fasta_file->nextEntry) {
			my $def = $entry->def;
			my $seq = $entry->seq;
			$def =~ s/;noChange=1//g;
			print OUT "$def\n";
			print OUT "$seq\n";
		}
		close OUT;
		close IN;
		$in = "$tmpdir/$basename.iterationMerge.$outCount.$iterateCount.fa";
	}
}

