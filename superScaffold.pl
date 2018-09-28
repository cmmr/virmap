#!/usr/bin/env perl



use warnings;
use strict;
use FAlite;
my $parents = {};
my $nucAcc2Parent = {};

open IN,"$ARGV[0]";

while (my $line = <IN>) {
	chomp $line;
	if ($line =~ m/^Parent:/) {
		$line =~ s/^Parent://g;
		my @parts = split /\|/, $line;
		$nucAcc2Parent->{$parts[3]} = $line;
		$parents->{$line} = {};
	}
	if ($line =~ m/^Child:/) {
		$line =~ s/^Child://g;
		my @parts = split /\|/, $line;
		unless (exists $nucAcc2Parent->{$parts[2]}) {
			die "$parts[2] has no parent, line = $line\n";
		}
		$parents->{$nucAcc2Parent->{$parts[2]}}->{$line} = 1;
	}

}
close IN;
open IN, "$ARGV[1]";

my $seqs = {};
my $offsets = {};
my $fasta_file = new FAlite(\*IN); # or any other filehandle
while (my $entry = $fasta_file->nextEntry) {

	my $def = $entry->def;
	$def =~ s/^>.*Scaffold.//g;
	if ($def =~ s/;startPos=(\d+)//g) {
		my $offset = $1;
		$offsets->{$def} = $offset;
	}
	my $seq = $entry->seq;
	$seqs->{$def} = $seq;
}
foreach my $parent (sort keys %$parents) {
	my @seqArray;
	if (exists $seqs->{$parent}) {
		my $seq = $seqs->{$parent};
		@seqArray = split //, $seq;
		unshift @seqArray, "N";
	}
	my $offset = "NA";
	if (exists ($offsets->{$parent})) {
		$offset = $offsets->{$parent};
	} else {
		$offset = "inf";
		foreach my $child (keys %{$parents->{$parent}}) {
			unless ($child =~ s/.*;pos=//g) {
				die "no position information in $child\n";
			} 
			unless ($child =~ s/;codonStart=.*//g) {
				die "can't strip end info off $child, no ;codonStart=\n";
			}
			my @matches = $child =~ m/(\d+)/g;
			unless (scalar @matches) {
				die "$child has no coordinates\n";
			}
			foreach my $match (@matches) {
				if ($match < $offset) {
					$offset = $match;
				}
			}
			
		}

	}
	foreach my $child (keys %{$parents->{$parent}}) {
		$child =~ m/;pos=(.*);codonStart/;
		my $posInfo = $1;
		my @matches = $posInfo =~ m/(\d+)/g;
		foreach my $match (@matches) {
			for (my $i = 0; $i < ($offset - $match); $i++) {
				unshift @seqArray , "N";
			}	
			if ($match < $offset) {
				$offset = $offset - ($offset - $match);
			}
		}		
		unless ($posInfo) {
			die "can't grab position information on $child\n";
		}	
		my $complement = 0;
		if ($posInfo =~ m/complement\((.*)\)$/) {
			$complement = 1;
			$posInfo = $1;
		}
		my @parts;
		if ($posInfo =~ m/join\((.*)\)$/) {
			$posInfo = $1;
			@parts = split /,/, $posInfo;
		} else {
			push @parts, $posInfo;
		}
		unless (scalar @parts) {
			die "no valid position information\n";
		}
		my $childSeq = $seqs->{$child};
		unless ($childSeq) {
			next;
		}
		if ($complement) {
			$childSeq =~ tr/ACGT/TGCA/;
			@parts = reverse(@parts);
		}
		my @childSeqArray = split //, $childSeq;
		my $length = scalar(@childSeqArray);
		my $i = 0;
		$child =~ m/codonStart=(\d+)/;
		my $codonOffset = 1 - $1;
		foreach my $segment (@parts) {
			$segment =~ m/<?(\d+)\.\.>?(\d+)/;
			unless (defined($1) and defined($2)) {
				next;
			}
			if ($complement) {
				my $start = $1 - ($offset - 1) + $codonOffset;
				my $stop = $2 - ($offset - 1) + $codonOffset;
				for (my $j = $stop; $j >=$start; $j--) {
					if ($i == scalar(@childSeqArray)) {
						last;
					}
					if (not defined $seqArray[$j]) {
						$seqArray[$j] = $childSeqArray[$i];
					} elsif ($seqArray[$j] eq "N" or $seqArray[$j] eq "n") {
						$seqArray[$j] = $childSeqArray[$i];
					} elsif ($childSeqArray[$i] eq "N" or $childSeqArray[$i] eq "n") {
						$i++;
						next;
					} elsif ($seqArray[$j] ne $childSeqArray[$i]) {
						$seqArray[$j] = lc($seqArray[$j]);;
					}
					$i++;
				}
			} else {
				my $start = $1 - ($offset - 1) - $codonOffset;
				my $stop = $2 - ($offset - 1) - $codonOffset;
				for (my $j = $start; $j <= $stop; $j++) {
					if ($i == scalar(@childSeqArray)) {
						last;
					}
					if (not defined $seqArray[$j]) {
						$seqArray[$j] = $childSeqArray[$i];
					} elsif ($seqArray[$j] eq "N" or $seqArray[$j] eq "n") {
						$seqArray[$j] = $childSeqArray[$i];
					} elsif ($childSeqArray[$i] eq "N" or $childSeqArray[$i] eq "n") {
						$i++;
						next;
					} elsif ($seqArray[$j] ne $childSeqArray[$i]) {
						$seqArray[$j] = lc($seqArray[$j]);;
					}
					$i++;
				}
			}
		}
	}
	for (my $i = 0; $i < scalar(@seqArray); $i++) {
		unless (defined($seqArray[$i])) {
			$seqArray[$i] = "N";
		}
	}
	shift(@seqArray);
	my $newSeq = join "", @seqArray;
	$newSeq =~ s/^[N|n]+//g;
	$newSeq =~ s/[N|n]+$//g;
	unless ($newSeq) {
		next;
	}
	print ">$parent\n";
	print "$newSeq\n";
	
}
