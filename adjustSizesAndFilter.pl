#!/usr/bin/env perl


use warnings;
use strict;
use Time::HiRes qw(time);
use threads;
use Thread::Queue qw( );
use threads::shared;
use FAlite;
use File::Temp qw/tempdir/;
use RocksDB;
use Compress::Zstd qw(compress decompress compress_mt);
use Sereal;


my $q = Thread::Queue->new();
my $goQ = Thread::Queue->new();
my $keysQ = Thread::Queue->new();
my $finished :shared;
$finished = 0;

my $function = shift @ARGV;
my $quant = 0;
unless ($function eq "build" or $function eq "quant") {
	die;
}
if ($function eq "quant") {
	$quant = 1;
}
my $pad = shift @ARGV;
unless ($pad =~ m/^\d+$/) {
	die "pad is not a number\n";
}

my $t1 = time();
my $startTime = $t1;
my $THREADS = $ARGV[3];
my $tmpdir = tempdir( DIR => "/dev/shm", CLEANUP => 1);
my @threads;
for (my $i = 0; $i < $THREADS; $i++) {
	push @threads, threads->create(\&worker, $quant, $tmpdir, $pad);
}

my $encoder = Sereal::Encoder->new();
my $decoder = Sereal::Decoder->new();
my $dbMASK = RocksDB->new("$tmpdir/MASK", { create_if_missing => 1, db_log_dir => "/dev/null", keep_log_file_num => 1, max_log_file_size => 1, db_stats_log_interval => 9999999, max_bytes_for_level_base => 100000000, target_file_size_base => 20000000, allow_mmap_reads => 'true'});
my $dbREAD = RocksDB->new("$tmpdir/READ", { create_if_missing => 1, db_log_dir => "/dev/null", keep_log_file_num => 1, max_log_file_size => 1, db_stats_log_interval => 9999999, max_bytes_for_level_base => 100000000, target_file_size_base => 20000000, allow_mmap_reads => 'true'});

my $skip = {};
open IN, "$ARGV[2]";
while (my $line = <IN>) {
	chomp $line;
	$skip->{$line} = 1;
}
close IN;

open IN, "$ARGV[1]";
my $fasta_file = new FAlite(\*IN); # or any other filehandle
while (my $entry = $fasta_file->nextEntry) {
	my $def = $entry->def;
	$def =~ s/^>//g;
	if (exists $skip->{$def}) {
		next;
	}
	my $compressed = compress($entry->seq);
	$dbMASK->put("$def", $compressed);

}
close IN;
my $sortThreads = $THREADS;
if ($sortThreads > 8) {
	$sortThreads = 8;
}
open IN, "lbzip2 -dc -n$THREADS $ARGV[0] | sort -k3 --buffer-size=10G --parallel=$sortThreads --compress-program=lz4 | cut -f1,3,4,6,10 |";



my $noHits = qr/noHits=1/;
my $lines = {};
my $lineAmount = {};
my $keysAmount :shared;
my $printed = {};
while (my $line = <IN>) {
	chomp $line;
	my @parts = split /\t/, $line, 3;
	my $reference = $parts[1];
	if (exists $skip->{$reference}) {
		next;
	}
	$lines->{$reference}->{$line} = 1;	
	$lineAmount->{$reference}++;
}
$keysAmount = scalar(keys %$lines);
close IN;
for (my $i = 0; $i < $THREADS; $i++) {
	$keysQ->enqueue($keysAmount);
}
close IN;
my $t2 = time();
my $elapsed = $t2 - $t1;
print STDERR "finished reading SAM after $elapsed seconds\n";
$t2 = $t1;


foreach my $hit (sort {$lineAmount->{$b} <=> $lineAmount->{$a}} keys %$lines) {
	my $compressed = compress($encoder->encode($lines->{$hit}));
	$dbREAD->put("$hit", $compressed);
	$q->enqueue($hit);
}
undef($lines);
undef($dbMASK);
undef($dbREAD);
sub worker {
	select(STDERR);
	$| = 1;
	my $tid = threads->tid();	
	my $quat = $_[0];
	my $tmpdir = $_[1];
	my $padded = $_[2];
	my $encoder = Sereal::Encoder->new();
	my $decoder = Sereal::Decoder->new();
	my $start = $goQ->dequeue();
	my $keysAmount = $keysQ->dequeue();
	my $dbMASK = RocksDB->new("$tmpdir/MASK", { read_only => 1, db_log_dir => '/dev/null', keep_log_file_num => 1, allow_mmap_reads => 'true'}) or die "can't open DB\n";
	my $dbREAD = RocksDB->new("$tmpdir/READ", { read_only => 1, db_log_dir => '/dev/null', keep_log_file_num => 1, allow_mmap_reads => 'true'}) or die "can't open DB\n";
	my $beginStrip = qr/^(\d+)S/;
	my $sizeQr = qr/size=(\d+);/;
	my $stripQr = qr/S\.(.)/;
	my $return = {};
	while (1) {
		my $isContig = 0;
		my $reference = $q->dequeue();
		unless (defined $reference) {
			my $peek = $q->peek();
			if (not defined $peek) {
				last;
			} else {
				redo;
			}
		}
		my $t1 = time;
		my $lineCont = $decoder->decode(decompress($dbREAD->get($reference)));
		if ($reference =~ m/^contig/) {
			$isContig = 1;
		}
		my $pileups = {};
		my $count;
		my $alternate = {};
		my $alternateAdjust = {};
		unless (scalar(keys %$lineCont)) {
			lockPrintStderr("$reference has no reads");
			next;
		}
		my $printed = 0;
		my $origScaffold = decompress($dbMASK->get($reference));
		my $scaffold = "a" . $origScaffold;
		my @scafLetters = split //, $scaffold;
		my $caseMask = [];
		foreach my $letter (@scafLetters) {
			if ($isContig) {
				push @$caseMask, 0;
				next;
			}
			if ($letter =~ m/[acgt]/) {
				push @$caseMask, 1;
			} elsif ($letter eq "n") {
				push @$caseMask, 2;
			} else {
				push @$caseMask, 0;
			}
		}
		unless(scalar(@{$caseMask})) {
			die "casemask is empty\n";
		}
		my $lastPos = scalar(@{$caseMask}) - $padded;
		my $firstPos = $padded + 1;
		my $time1 = time;
		my $parseCigarTime = 0;
		foreach my $line (keys %$lineCont) {
			my @parts = split /\t/, $line;
			my $cigar = $parts[3];
			my $seq = $parts[4];
			my $start = $parts[2];
			$parts[0] =~ m/$sizeQr/;
			my $size = $1;
			unless ($size) {
				$size = 1;
			}
			$count += $size;
			if ($quant) {
				next;
			}
			my $beforeParse = time;
			my ($return, $Mtotal, $Dtotal) = parseCigar($cigar);
			my $afterParse = time;
			$parseCigarTime += ($afterParse - $beforeParse);
			my $length = length($seq);
			if ($Mtotal < $length * .5) {
				next;
			}
			my $negative = 0;
			if ($cigar =~ m/$beginStrip/) {
				$negative = $1;
			}
			my @letters = split //, $seq;
			my @directives = @{$return};
			my $directivePos = 0;
			my $dist = 0;
			my $Dsize = $size;
			if ($Dtotal > .2 * $length) {
				$Dsize = 0;
			}
			my $consect = 0;
			my $insertAddPos = 0;
			my $insertStr = "";
			for (my $i = 0; $i < scalar(@letters); $i++) {
				my $insertPos = $start + $dist - $negative;
				if (not $isContig and $insertPos > $lastPos) {
					$letters[$i] = lc($letters[$i]);
				} 
				if (not $isContig and $insertPos < $firstPos) {
					$letters[$i] = lc($letters[$i]);
				}
				if ($consect and $directives[$directivePos] eq "I") {
					$insertStr .= $letters[$i];
					$consect++;
					$directivePos++;
					next;
				}
				if ($consect and $directives[$directivePos] ne "I") {
					$alternate->{$insertAddPos}->{$insertStr} += $size;
					$alternateAdjust->{$insertAddPos} += $size;
					#$dist += length($insertStr);		
					#$insertPos += length($insertStr);
					$insertStr = "";
					$consect = 0;
					
				}
				if ($directives[$directivePos] eq "M" or ($directives[$directivePos] eq "S" and $insertPos < $firstPos)) {
					$pileups->{$insertPos}->{$letters[$i]} += $size;
					$directivePos++;
					$dist++;
					next;
				}
				if ($directives[$directivePos] eq "D") {
					if ($Dsize) {
						$pileups->{$insertPos}->{"D"} += $Dsize;
					}
					$directivePos++;
					$dist++;
					redo;
				}
				if (not $consect and $directives[$directivePos] eq "I") {
					$directivePos++;
					$insertAddPos = $insertPos;
					$consect = 1;
					$insertStr = "$letters[$i]";
					next;
				}
				if ($directives[$directivePos] eq "S" and $insertPos >= 1) {
					$pileups->{$insertPos}->{"S.$letters[$i]"} += $size; 
					$dist++;
					$directivePos++;
					next;
				}
				if ($directives[$directivePos] eq "X.M") {
					$dist++;
					$directivePos++;
					next;
				}
				if ($directives[$directivePos] eq "X.D") {
					$dist++;
					$directivePos++;
					redo;
				}
				if ($directives[$directivePos] eq "X.I") {
					$directivePos++;
					next;
				}
			
			}
		}
		my $time2 = time;
		my $pileTime = $time2 - $time1;
		$time1 = $time2;
		$lineCont = {};
		my @sequence;
		if (not scalar(keys %$pileups) and not scalar(keys %{$alternate}) and not $quant) {
			lockPrintStderr("$reference has no valid hits");
			next;
		}
		if ($quant) {
			$pileups->{0} = 1;
		}
		my @positions = sort {$a <=> $b} keys %{$pileups};
		my $start = $positions[0];
		my $stop = $positions[$#positions];
		if (scalar(keys %$alternate) and not $quant) {
			my @altPos = sort {$a <=> $b} keys %{$alternate};
			if ($altPos[0] < $start) {
				$start = $altPos[0];
			}
			if ($altPos[$#altPos] > $stop) {
				$stop = $altPos[$#altPos]
			}
		}
		for (my $i = $start; $i <= $stop; $i++) {
			if ($quant) {
				last
			}
			if (not exists $pileups->{$i} and not exists $alternate->{$i}) {
				if ($caseMask->[$i] == 2) {
					push @sequence, "n";
				} else {
					push @sequence, "N";
				}
			} else {
				my @order;
				my $noOrig = 0;
				if (exists $pileups->{$i}) {
					my @keys = sort {$pileups->{$i}->{$b} <=> $pileups->{$i}->{$a}} keys %{$pileups->{$i}};
					if (scalar(@keys) > 1) {
						my @toDelete;
						foreach my $letter (@keys) {
							if ($letter =~ m/S\./) {
								push @toDelete, $letter;
							}
						}
						if (scalar(@toDelete) >= 1 and scalar(@toDelete) != scalar(@keys)) {
							foreach my $key (@toDelete) {
								delete $pileups->{$i}->{$key};
							}
							@keys = keys %{$pileups->{$i}};
						} 
					}
					@order = sort {$pileups->{$i}->{$b} <=> $pileups->{$i}->{$a}} @keys;
				} else {
					$noOrig = 1;
					@order = qw(X);
					$pileups->{$i}->{$order[0]} = -1000000;
				}
				if (exists $alternate->{$i}) {
					my @altOrder = sort {$alternate->{$i}->{$b} <=> $alternate->{$i}->{$a} || $b cmp $a} keys %{$alternate->{$i}};
					my $altInsert = $altOrder[0];
					if ($alternate->{$i}->{$altInsert} >= ($pileups->{$i}->{$order[0]} - $alternateAdjust->{$i})) {
						if (defined $altOrder[1] and $alternate->{$i}->{$altInsert} == $alternate->{$i}->{$altOrder[1]}) {
							$altInsert = "n" x length($altInsert);;
						}
						if ($caseMask->[$i] == 1) {
							$altInsert = lc($altInsert);
						}
						push @sequence, $altInsert;
					}
					if ($noOrig) {
						delete $pileups->{$i};
					}
				}
				if (exists $pileups->{$i}) {
					if (defined $order[1] and $pileups->{$i}->{$order[0]} == $pileups->{$i}->{$order[1]}) {
						push @sequence, "n";
					} else {
						unless ($order[0] eq "D") {
							if ($order[0] =~ m/$stripQr/) {
								$order[0] = $1;
							}
							unless($order[0]) {
								lockPrintStderr("WARNING: putting nothing: $order[0] in sequence for $reference at pos $i");
							}
							if (not defined $caseMask->[$i] or $caseMask->[$i] == 1) {
								$order[0] = lc($order[0]);
							}
							push @sequence, $order[0];
						}
					}
				}
			}	
		}
		$time2 = time;
		my $rebuildTime = $time2 - $time1;
		my $outSeq = "";
		unless ($quant) {
			$outSeq = join "", @sequence;
		} else {
			$outSeq = $origScaffold;
		}
		$pileups = {};
		$alternate = {};
		my $sizeDiff = 0;
		my $lengthDiff = 0;
		my $sizeErrStr = "";
		my $beforeLength;
		my $newLength = length($outSeq);
		if ($reference =~ m/length=(\d+)/) {
			$beforeLength = $1;
		} else {
			$beforeLength = $newLength;
		}
		if ($reference =~ m/length=\d+/) {
			$reference =~ s/length=\d+/length=$newLength/g;
		} else {
			$reference .= ";length=$newLength";
		}
		$lengthDiff = $newLength / $beforeLength;
		$sizeErrStr = " (lengthDiff = $lengthDiff)";
		if ($quant) {
			my $beforeSize;
			if ($reference =~ m/;size=(\d+)/g) {
				$beforeSize = $1;
			}
			unless ($beforeSize) {
				$beforeSize = $count;
			}
			$sizeDiff = $count / $beforeSize;
			if ($reference =~ m/size=\d+/) {
				$reference =~ s/size=\d+/size=$count/g;
			} else {
				$reference .= ";size=$count";
			}
			$sizeErrStr = " (sizeDiff = $sizeDiff; lengthDiff = $lengthDiff)";
		}
		$return->{$reference} = $outSeq;
		my $t2 = time;
		my $elapsed = $t2 - $t1;
		{lock($finished);
		$finished++;
		if ($elapsed > 5) {
			lockPrintStderr("finished on thread $tid ($finished / $keysAmount) time: $elapsed seconds $reference$sizeErrStr, pileTime = $pileTime, rebuildTime = $rebuildTime, parseCigarTime = $parseCigarTime");;
			select->flush();
		}
		}
	}
	return($return);
}
for (my $i = 0; $i < $THREADS; $i++) {
	$q->enqueue(undef);
}
for (my $i = 0; $i < $THREADS; $i++) {	
	$goQ->enqueue("GO");
}


$t2 = time();
$elapsed = $t2 - $t1;
$t1 = $t2;
my $seqs = {};
foreach my $thr (@threads) {
	my $tid = $thr->tid();
	my $return = $thr->join();
	foreach my $def (keys %$return) {
		$seqs->{$def} = $return->{$def};
	}
}
foreach my $def (sort keys %$seqs) {
	unless (exists $seqs->{$def} and defined $seqs->{$def}) {
		next;
	}
	print ">$def\n";
	print "$seqs->{$def}\n";
}


my $endTime = time();
my $overallTime = $endTime - $startTime;
print STDERR "Overall pileup time $overallTime seconds\n";

sub parseCigar {
	my $cigar = $_[0];
	my @directives;
	my @cigNum = split /[A-Z]/, $cigar;
	my @cigLet = split /\d+/, $cigar;
	shift @cigLet;
	my $Mtotal = 0;
	my $Dtotal = 0;
	my $largestD = 0;
	my $largestDi = 0;
	for (my $i = 0; $i < scalar(@cigNum); $i++) {
		if ($cigLet[$i] eq "M") {
			$Mtotal += $cigNum[$i];
		}
		if ($cigLet[$i] eq "D") {
			$Dtotal += $cigNum[$i];
			if ($cigNum[$i] > $largestD) {
				$largestD = $cigNum[$i];
				$largestDi = $i;
			}
		}
		for (my $j = 0; $j < $cigNum[$i]; $j++) {
			push @directives, $cigLet[$i];
		}
	} 
	if ($largestD > .5 * $Mtotal) {
		@directives = ();
		$Mtotal = 0;
		$Dtotal = 0;
		my $leftM = 0;
		my $rightM = 0;
		for (my $i = 0; $i < scalar(@cigNum); $i++) {
			if ($i < $largestDi and $cigLet[$i] eq "M") {
				$leftM += $cigNum[$i];
			}
			if ($i > $largestDi and $cigLet[$i] eq "M") {
				$rightM += $cigNum[$i];
			}
		}
		for (my $i = 0; $i < scalar(@cigNum); $i++) {
			my $letter = $cigLet[$i];
			if ($leftM > $rightM and $i > $largestDi) {
				$letter = "X.$letter";
			} elsif ($rightM > $leftM and $i < $largestDi) {
				$letter = "X.$letter";
			} elsif ($leftM == $rightM or $i == $largestDi) {
				$letter = "X.$letter";
			}
			for (my $j = 0; $j < $cigNum[$i]; $j++) {
				if ($letter eq "M") {
					$Mtotal++;
				}
				if ($letter eq "D") {
					$Dtotal++;
				}
				push @directives, $letter;
			}
		}
	}
	my $directive = \@directives;
	return($directive, $Mtotal, $Dtotal);


}
sub lockPrintStderr {
	lock($finished);
	my $tid = threads->tid();
	my $printout = $_[0];
	print STDERR "THREAD $tid - $printout\n";
}
