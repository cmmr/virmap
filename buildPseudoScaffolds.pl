#!/usr/bin/env perl


use warnings;
use strict;
use FAlite;


open REF, "cat $ARGV[0] |";
my $goodReferences = {};
while (my $line = <REF>) {
	chomp $line;
	if ($line =~ m/^Parent:/) {
		$line =~ s/^Parent://g;
	}
	if ($line =~ m/^Child:/) {
		$line =~ s/^Child://g;
	}
	$goodReferences->{$line} = 1;
}
close REF;

my $referenceHits = {};

open AASAM, "lbzip2 -d -c -n3 $ARGV[2] |";
my $nuc3 = {};
my $bestHitPerRead = {};
my $lowestMismatchPerRead = {};
my $lowestScorePerRead = {};
my $reference2Reads = {};
my $readInfo = {};
my $sizeQR = qr/size=(\d+)/;
while (my $line = <AASAM>) {

	chomp $line;
	if ($line =~ m/^@/) {
		next;
	}
	my @parts = split /\t/, $line, 4;
	my $reference = $parts[2];
	unless (exists $goodReferences->{$reference}) {
		next;
	}
	my @subparts = split /\t/, $parts[3];
	unless ($subparts[0] and $subparts[6] and $subparts[9] and $subparts[14] and $subparts[15]) {
		print STDERR "line broke = $line\n";
	}
	my $readName = $parts[0];
	$subparts[15] =~ s/ZS:i://g;
	$subparts[14] =~ s/ZF:i://g;
	$subparts[11] =~ s/ZR:i://g;
	$subparts[13] =~ s/ZI:i://g;

	my $mismatch = (100 - $subparts[13]);
	my $score = $subparts[11];
	if (exists $readInfo->{$readName}->{$reference}) {
		if ($readInfo->{$readName}->{$reference}->{'score'} >= $score) {
			next;
		}
	}
	$readInfo->{$readName}->{$reference}->{'score'} = $score;
	$readInfo->{$readName}->{$reference}->{'mismatch'} = $mismatch;
	$readInfo->{$readName}->{$reference}->{'start'} = $subparts[0];
	$readInfo->{$readName}->{$reference}->{'cigar'} = $subparts[2];
	$readInfo->{$readName}->{$reference}->{'frame'} = $subparts[14];
	$readInfo->{$readName}->{$reference}->{'dnaStart'} = $subparts[15];
	$reference2Reads->{$reference}->{$readName} = 1;
	unless (exists $lowestMismatchPerRead->{$readName}) {
		$lowestMismatchPerRead->{$readName} = $mismatch;
		$lowestScorePerRead->{$readName} = $score;
		$bestHitPerRead->{$readName}->{$reference} = 1;
	} elsif ($mismatch < $lowestMismatchPerRead->{$readName}) {
		delete $bestHitPerRead->{$readName};
		$lowestMismatchPerRead->{$readName} = $mismatch;
		$lowestScorePerRead->{$readName} = $score;
		$bestHitPerRead->{$readName}->{$reference} = 1;
	} elsif ($score == $lowestMismatchPerRead->{$readName}) {
		if ($score > $lowestScorePerRead->{$readName}) {
			delete $bestHitPerRead->{$readName};
			$lowestMismatchPerRead->{$readName} = $mismatch;
			$lowestScorePerRead->{$readName} = $score;
			$bestHitPerRead->{$readName}->{$reference} = 1;
		} elsif ($score == $lowestScorePerRead->{$readName}) {
			$bestHitPerRead->{$readName}->{$reference} = 1;
		}
	}
}
close AASAM;
open NUCSAM, "lbzip2 -d -c -n3 $ARGV[3] |";
my $nucReadInfo = {};
my $nucLowestMismatchPerRead = {};
my $nucBestPerRead = {};
my $nucReferenceToReads = {};
while (my $line = <NUCSAM>) {
	chomp $line;
	if ($line =~ m/^@/) {
		chomp $line;
	}
	my @parts = split /\t/, $line, 7;
	unless (exists $goodReferences->{$parts[2]}) {
		next;
	}
	my $readName = $parts[0];
	my $referenceHit = $parts[2];
	$line =~ m/NM:i:(\d+)/;
	my $mismatch = $1;
	$nucReadInfo->{$readName}->{$referenceHit}->{'start'} = $parts[3];
	$nucReadInfo->{$readName}->{$referenceHit}->{'mismatch'} = $mismatch;
	$nucReadInfo->{$readName}->{$referenceHit}->{'samcode'} = $parts[1];
	$nucReadInfo->{$readName}->{$referenceHit}->{'cigar'} = $parts[5];
	unless ($nucReadInfo->{$readName}->{$referenceHit}->{'cigar'}) {
		print STDERR "no cigar on $readName on $referenceHit\n";
	}	
	$nucReferenceToReads->{$referenceHit}->{$readName} = 1;
	unless (exists $nucBestPerRead->{$readName}) {
		$nucBestPerRead->{$readName}->{$referenceHit} = 1;
		$nucLowestMismatchPerRead->{$readName} = $mismatch;
	} elsif ($nucLowestMismatchPerRead->{$readName} > $mismatch) {
		delete $nucBestPerRead->{$readName};
		$nucBestPerRead->{$readName}->{$referenceHit} = 1;
		$nucLowestMismatchPerRead->{$readName} = $mismatch;
	} elsif ($nucLowestMismatchPerRead->{$readName} == $mismatch) {
		$nucBestPerRead->{$readName}->{$referenceHit} = 1;
	}
}
close NUCSAM;
open READS, "lbzip2 -d -c -n3 $ARGV[1] |";
my $fasta_file = new FAlite(\*READS);
my $reads = {};

while(my $entry = $fasta_file->nextEntry) {
	my $header = $entry->def;
	$header =~ s/^>//g;
	unless (exists $readInfo->{$header} or exists $nucReadInfo->{$header}) {
		next;
	}
	my $seq = $entry->seq;
	$reads->{$header} = $seq;
}
close READS;


my @references = sort keys %$reference2Reads;
foreach my $reference (@references) {
	my @reads = keys %{$reference2Reads->{$reference}};
	my $startCount = scalar(@reads);
	my $soloCount = 0;
	foreach my $read (@reads) {
		unless (exists $bestHitPerRead->{$read}->{$reference}) {
			delete $reference2Reads->{$reference}->{$read};
		} else {
			if (scalar(keys %{$bestHitPerRead->{$read}}) == 1) {
				$soloCount++;
			} else {
				delete $reference2Reads->{$reference}->{$read};
			}
		}

	}
	@reads = keys %{$reference2Reads->{$reference}};
	my $endCount = scalar(@reads);
	if ($endCount == 0 or $soloCount/$endCount <= .20) {
		delete $reference2Reads->{$reference};
	}
}
my @nucReferences = sort keys %$nucReferenceToReads;
foreach my $reference (@nucReferences) {
	my @reads = keys %{$nucReferenceToReads->{$reference}};
	my $startCount = scalar(@reads); 
	my $soloCount = 0;
	foreach my $read (@reads) {
		unless (exists $nucBestPerRead->{$read}->{$reference}) {
			delete $nucReferenceToReads->{$reference}->{$read};
		} else {
			if (scalar keys (%{$nucBestPerRead->{$read}}) == 1) {
				$soloCount++;
			} else {
				delete $nucReferenceToReads->{$reference}->{$read};
			}
		}
	}
	@reads = keys %{$nucReferenceToReads->{$reference}};
	my $endCount = scalar(@reads);
	if ($endCount == 0 or $soloCount/$endCount <= .20) {
		delete $nucReferenceToReads->{$reference};
		delete $goodReferences->{$reference};
	} 
}


my $alternateLoci = {};
my $aaStops = {};
foreach my $ref (sort keys %$reference2Reads) {
	my $highest = 0;
	foreach my $readName (sort {$readInfo->{$a}->{$ref}->{'start'} <=> $readInfo->{$b}->{$ref}->{'start'}} keys %{$reference2Reads->{$ref}}) {
		$readName =~ m/$sizeQR/;
		my $countAdd = $1;
		my $cigar = $readInfo->{$readName}->{$ref}->{'cigar'};
		my $start = $readInfo->{$readName}->{$ref}->{'start'};
		my $frame = $readInfo->{$readName}->{$ref}->{'frame'};
		my $dnaStart = $readInfo->{$readName}->{$ref}->{'dnaStart'};
		my $codonLength = 0;
		my @directives;
		my $deleteCount = 0;
		if ($cigar =~ m/^(\d+)M$/) {
			$codonLength = $1;
			for (my $i = 0; $i < $codonLength; $i++) {
				push @directives, "M";
			}
		} else {
			my $return = parseCigar($cigar);
			@directives = @{$return};
			foreach my $directive (@directives) {
				if ($directive eq "D") {
					$deleteCount++;
				}
			}
			$codonLength = scalar(@directives) - $deleteCount;
		}
		$start = $start - 1;
		$referenceHits->{$ref} += $countAdd;
		my $read = $reads->{$readName};
		my $origSeq = $read;
		my $readLength = length($read);
		if ($frame < 0) {
			$read = reverse($read);
			$read =~ tr/ACGT/TGCA/;
			$dnaStart = $readLength - $dnaStart + 1;
		}
		$dnaStart = $dnaStart - 1;
		$read = substr($read, $dnaStart, ($readLength - $dnaStart));
		my $nucStart = ($start * 3);
		my @nucLetters = split //, $read;
		if ($nucStart < 0) {
			my $adjust = $nucStart * -1;
			for (my $i = 0; $i < $adjust; $i++) {
				shift @nucLetters;
			}
			$nucStart = 0;
		}
		my $adjustedRead = join "", @nucLetters;
		my @nucLetters3 = ( $adjustedRead =~ m/.../g );
		my $Dcount = $countAdd;
		if ($deleteCount > (.3 * scalar(@nucLetters3))) {
			$Dcount = 0;
		}
		unless (scalar(@nucLetters3)) {
			next;
		}
		if ($start + $codonLength > $highest) {
			$highest = $start + $codonLength;
		}
		my $directivePos = 0;
		my $dist = 0;
		for (my $i = 0; $i < $codonLength; $i++) {
			unless (defined($nucLetters3[$i])) {
				print STDERR "broken codon on pos $i on $ref with read $readName\n";
			}
			if ($directives[$directivePos] eq "M" or $directives[$directivePos] eq "S") {
				$nuc3->{$ref}->{$dist + $start}->{$nucLetters3[$i]} += $countAdd;
				$dist++;
				$directivePos++;
				next;
			}
			if ($directives[$directivePos] eq "D") {
				if ($Dcount) {
					$nuc3->{$ref}->{$dist + $start}->{"DDD"} += $Dcount;
				}
				$directivePos++;
				$dist++;
				redo;
			}
			if ($directives[$directivePos] eq "I") {
				$alternateLoci->{$ref}->{$dist + $start}->{$nucLetters3[$i]} += $countAdd;
				$directivePos++;
				next;
			}
			print STDERR "Should never reach here, cumulative cigar length exceeds expectations\n";
			die;
		}
	}
	$aaStops->{$ref} = $highest;
}

foreach my $reference (sort keys %{$nuc3}) {
	my @scaffold;
	for (my $pos = 0; $pos < $aaStops->{$reference}; $pos++) {
		my @top;
		if (exists $nuc3->{$reference}->{$pos} or exists $alternateLoci->{$reference}->{$pos}) {
			my @top;
			my $noNuc = 0;
			if ($nuc3->{$reference}->{$pos}) {
				@top = sort {$nuc3->{$reference}->{$pos}->{$b} <=> $nuc3->{$reference}->{$pos}->{$a}} keys %{$nuc3->{$reference}->{$pos}};
			} else {
				$noNuc = 1;
				$top[0] = "X";
				$nuc3->{$reference}->{$pos}->{$top[0]} = 0;
			}
			if (exists $alternateLoci->{$reference}->{$pos}) {
				my @altTop = sort {$alternateLoci->{$reference}->{$pos}->{$b} <=> $alternateLoci->{$reference}->{$pos}->{$a}} keys %{$alternateLoci->{$reference}->{$pos}};
				my $altTop = $altTop[0];
				if ($alternateLoci->{$reference}->{$pos}->{$altTop[0]} >= $nuc3->{$reference}->{$pos}->{$top[0]}) {
					if (defined $altTop[1] and $alternateLoci->{$reference}->{$pos}->{$altTop[0]} == $alternateLoci->{$reference}->{$pos}->{$altTop[1]}) {
						$altTop = "nnn";
					}
				}
				if ($noNuc) {
					delete $nuc3->{$reference}->{$pos};
				}
			}
			if (exists $nuc3->{$reference}->{$pos}) {
				my $letter = $top[0];
				if (defined $top[1] and ($nuc3->{$reference}->{$pos}->{$top[0]} == $nuc3->{$reference}->{$pos}->{$top[1]})) {
					$letter = "nnn";
				}
				if ($letter eq "DDD") {
					push @scaffold, "NNN";
					next;
				} else {
					push @scaffold, $letter;
				}
			}
		} else {
			push @scaffold, "NNN";
		}
	}
	my $outFrag = join "", @scaffold;
	print ">Scaffold.$reference\n";
	print "$outFrag\n";
}
my $nucArrays = {};
my $lowestPerNuc = {};
my $highestPerNuc = {};
my $nucReferenceHits = {};
my $alternateNuc = {};
foreach my $reference (keys %$nucReferenceToReads) {	
	my $lowest = "inf";
	my $highest = 0;
	foreach my $read (keys %{$nucReferenceToReads->{$reference}}) {
		$read =~ m/$sizeQR/;
		my $sizeAdd = $1;
		$nucReferenceHits->{$reference} += $sizeAdd;	
		my $start = $nucReadInfo->{$read}->{$reference}->{'start'};
		my $cigar = $nucReadInfo->{$read}->{$reference}->{'cigar'};
		unless ($cigar) {
			print STDERR "Broken nucleotide cigar on $read on $reference\n";
			next;
		}
		my @directives;
		my $matchLength = 0;
		my $deleteCount = 0;
		if ($cigar =~ m/^(\d+)M$/) {
			$matchLength = $1;
			my $temp = "M" x $matchLength;
			@directives = split //, $temp;
		} else {
			my $return = parseCigar($cigar);
			@directives = @{$return};
			foreach my $directive (@directives) {
				if ($directive eq "D") {
					$deleteCount++;
				}
			}
			$matchLength = scalar(@directives) - $deleteCount;
		}
		if ($start < $lowest) {
			$lowest = $start;
		}
		my $status = $nucReadInfo->{$read}->{$reference}->{'samcode'};
		if ($lowest > $start) {
			$lowest = $start;
		}
		unless (exists $reads->{$read}) {
			print STDERR "$read is not in read set\n";
			next;
		}
		my $localLetters = $reads->{$read};
		if ($status == 16 or $status == 272) {
			$localLetters = reverse($localLetters);
			$localLetters =~ tr/ACGT/TGCA/;
		}
		my @readLetters = split //, $localLetters;
		my $Dadd = $sizeAdd;
		if ($deleteCount > (.3 * scalar(@readLetters))) {
			$Dadd = 0;
		}
		if ($start + $matchLength > $highest) {
			$highest = $start + scalar(@readLetters);
		}
		my $directivePos = 0;
		my $dist = 0;
		for (my $i = 0; $i < $matchLength; $i++) {
			if ($directives[$directivePos] eq "M" or $directives[$directivePos] eq "S") {
				$nucArrays->{$reference}->{$start + $dist}->{$readLetters[$i]} += $sizeAdd;
				$dist++;
				$directivePos++;
				next;
			}
			if ($directives[$directivePos] eq "D") {
				if ($Dadd) {
					$nucArrays->{$reference}->{$start + $dist}->{"D"} += $Dadd;
				}
				$directivePos++;
				$dist++;
				redo;
			}
			if ($directives[$directivePos] eq "I") {
				$alternateNuc->{$reference}->{$start + $dist}->{$readLetters[$i]} += $sizeAdd;
				$directivePos++;
				next;
			}
			print STDERR "i= $i; directivePos = $directivePos($directives[$directivePos]); dist = $dist, cigar = $cigar code = $status\n";
			print STDERR "This should never be reached, too many nuc letters\n";
			die;
		}
	}
	$lowestPerNuc->{$reference} = $lowest;
	$highestPerNuc->{$reference} = $highest;
		
}

foreach my $reference (sort keys %$nucReferenceToReads) {
	my @newRef;
	for (my $i = $lowestPerNuc->{$reference}; $i < $highestPerNuc->{$reference}; $i++) {
		if (exists $nucArrays->{$reference}->{$i} or exists $alternateNuc->{$reference}->{$i}) {
			my @letters;
			my $noNuc = 0;
			if (exists $nucArrays->{$reference}->{$i}) {
				@letters = sort {$nucArrays->{$reference}->{$i}->{$b} <=> $nucArrays->{$reference}->{$i}->{$a}} keys %{$nucArrays->{$reference}->{$i}};
			} else {
				$noNuc = 1;
				@letters = qw(X);
				$nucArrays->{$reference}->{$i}->{$letters[0]} = 0;
			}
			if (exists $alternateNuc->{$reference}->{$i}) {
				my @altTop = sort {$alternateNuc->{$reference}->{$i}->{$b} <=> $alternateNuc->{$reference}->{$i}->{$a}} keys %{$alternateNuc->{$reference}->{$i}};
				my $altTop = $altTop[0];
				if ($alternateNuc->{$reference}->{$i}->{$altTop} >= $nucArrays->{$reference}->{$i}->{$letters[0]}) {
					if (defined $altTop[1] and $alternateNuc->{$reference}->{$i}->{$altTop[1]} == $alternateNuc->{$reference}->{$i}->{$altTop}) {
						$altTop = "n"
					}
				}
				if ($noNuc) {
					delete $nucArrays->{$reference}->{$i};
				}
			}
			if (exists $nucArrays->{$reference}->{$i}) {
				if (scalar(@letters) > 1 and $nucArrays->{$reference}->{$i}->{$letters[0]} == $nucArrays->{$reference}->{$i}->{$letters[1]}) {
					push @newRef, "n";
				} elsif ($letters[0] eq "D") {
					push @newRef, "N";
					next;
				} else {
					push @newRef, $letters[0];
				}
			}
		} else {
			push @newRef, "N";
		}
	}
	my $newRef = join "", @newRef;
	print ">NUC.Scaffold.$reference;startPos=$lowestPerNuc->{$reference}\n";
	print "$newRef\n";

}





sub parseCigar {
	my $cigar = $_[0];
	my @directives;
	my @cigNum = split /[A-Z]/, $cigar;
	my @cigLet = split /\d+/, $cigar;
	shift @cigLet;
	for (my $i = 0; $i < scalar(@cigNum); $i++) {
		for (my $j = 0; $j < $cigNum[$i]; $j++) {
			push @directives, $cigLet[$i];
		}
	}
	return(\@directives);

}












