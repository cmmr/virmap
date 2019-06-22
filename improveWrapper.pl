#!/usr/bin/env perl


use warnings;
use strict;
use threads;
use threads::shared;
use Thread::Queue qw( );
use FAlite;
use File::Temp qw/tempdir/;
use Digest::MD5 qw(md5_hex);
use RocksDB;
use Compress::Zstd qw(compress decompress compress_mt);
use Sereal;
use Time::HiRes qw(usleep);



my $t1 = time;
my $startTime = $t1;
my $in = $ARGV[0];
my $reads = $ARGV[1];
my $threads = $ARGV[2];
my $gbBlastn = $ARGV[3];
my $gbBlastx = $ARGV[4];
my $taxaJson = $ARGV[5];
my $timeLimit = $ARGV[6];
my $strict = $ARGV[7];
my @threads;
my $lock :shared;
my $print :shared;
my $totalJobs :shared;
my $finishedJobs :shared;
my $devShmTmp = tempdir( DIR => "/dev/shm", CLEANUP => 1);
my $returnQ = Thread::Queue->new();
my $goQ = Thread::Queue->new();
my $encoder = Sereal::Encoder->new();
my $decoder = Sereal::Decoder->new();
my $blastThreads = $threads;
if ($blastThreads > 1) {
	$blastThreads = 1;
}
for (my $i = 0; $i < $blastThreads; $i++) {
	$threads[$i] = threads->create(\&worker, $devShmTmp);
}


if (not defined ($strict) or $strict ne "strict") {
	$strict = "";
}
my $whiteList = "";
if (-e $ARGV[$#ARGV]) {
	$whiteList = $ARGV[$#ARGV];
}
my $basename = `basename $in`;
chomp $basename;
my $iterateCount = 0;
my $tmpdir = tempdir( CLEANUP => 1 );
my $accQR = qr/^.*\|.*\|.*\|([^\|]+)\|/;
my $seen = {};
my $pass = {};
my $max = {};
my $finalFilter = 0;
my $beforeSizes = {};
my $afterSizes = {};
my $seenSizes = {};
my $saved = {};
my $savedDef = {};
my $sizeQr = qr/size=(\d+)/;
my $onePerc = {};
my $status = {};
my $radius = .80;
my $maxCycles = 40;
my $maxRadius = .97;
my $radGrowth = ($maxRadius - $radius) / ($maxCycles / 2);
my $lastAdded = {};
my $padLength = 300;
my $inStrict = "";
my $local = "local=t ignorefrequentkmers=f greedy=f excludefraction=0";
my $params = "vslow=t";
while (1) {
	print STDERR "Cycle $iterateCount started\n";
	if ($radius > $maxRadius) {
		$radius = $maxRadius;
		$padLength = 0;
		$params = "ignorefrequentkmers=f greedy=f excludefraction=0";
		$inStrict = "strict";
		$local = "";
	}
	if ($strict) {
		$inStrict = $strict;
	}
	print STDERR "Current radius = $radius\n";
	my $before = {};
	my $after = {};
	if ($iterateCount == $maxCycles) {
		print STDERR "cycles exhausted $iterateCount\n";
		last;
	}
	open IN, "$in";
	my $fasta_file = new FAlite(\*IN); # or any other filehandle
	my $beforeDef = {};
	my $beforeSize = {};
	my $sortLengths = {};
	while (my $entry = $fasta_file->nextEntry) {
		my $def = $entry->def;
		$def =~ s/^>//g;
		my $seq = $entry->seq;
		my $acc = getAcc($def);;
		unless (defined $acc) {
			print STDERR "getAcc failed $def\n";
			next;
		}
		if ($acc eq "NA") {
			die "no acc on $def\n";
		}
		$before->{$acc} = $seq;
		$beforeDef->{$acc} = $def;
		my $temp = $seq;
		$temp =~ s/N//g;
		$sortLengths->{$acc} = length($temp);
		if ($iterateCount > 2) {
			my $md5 = md5_hex($seq);
			$seen->{$acc}->{$md5} = 1;
		}
	}
	close IN;
	print STDERR "Merge Wrapper Before\n";
	system("mergeWrapper.pl $in $threads $inStrict > $tmpdir/$basename.iterationImprove.$iterateCount.merged.fa");
	print STDERR "END Merge Wrapper Before\n";
	$in = "$tmpdir/$basename.iterationImprove.$iterateCount.merged.fa";
	system("addPad.pl $in $padLength > $tmpdir/$basename.iterationImprove.$iterateCount.merged.padded.fa");
	$in = "$tmpdir/$basename.iterationImprove.$iterateCount.merged.padded.fa";
	system("cat $in | toUC.pl all > $devShmTmp/ForBbmap.fa");
	print STDERR "bbMap\n";
	system("bbmap.sh -Xmx32g usejni=t nodisk=t deterministic=t ref=$devShmTmp/ForBbmap.fa ambiguous=random in=$reads $params 32bit=t sam=1.3 threads=$threads $local outm=stdout minid=$radius noheader=t | lbzip2 -c -n$threads  > $tmpdir/$basename.sam.bz2");
	print STDERR "Adjust Sizes\n";
	open OUT, ">$tmpdir/$basename.good";
	close OUT;
	system("adjustSizesAndFilter.pl build $padLength $tmpdir/$basename.sam.bz2 $in $tmpdir/$basename.good $threads > $tmpdir/$basename.iterationImprove.$iterateCount.fa");
	$in = "$tmpdir/$basename.iterationImprove.$iterateCount.fa";
	open IN, "$in";
	my $afterDef = {};
	$fasta_file = new FAlite(\*IN);	
	while (my $entry = $fasta_file->nextEntry) {
		my $def = $entry->def;
		my $seq = $entry->seq;
		$def =~ s/^>//g;
		my $acc = getAcc($def);
		if (exists $pass->{$def}) {
			next;
		}
		if ($acc eq "NA") {
			die "no acc on after $def\n";
		}
		unless (length($seq)) {
			print STDERR "$acc has no seq. but has this def: $def\n";
			next;
		}

		#acc is short name
		$after->{$acc} = $seq;
		#real header names after
		$afterDef->{$acc} = $def;
	}
	close IN;
	my $bad = 0;
	print STDERR "\n";
	my $hitTimeLimit = 0;
	my $nowTime = time;
	my $runTime = $nowTime - $startTime;
	if ($runTime > $timeLimit) {
		print STDERR "Time limit breached: $runTime seconds elasped, $timeLimit seconds limit\n";
		$hitTimeLimit = 1;
	}
	if (-d "$devShmTmp/SEND") {
		system("rm -rf $devShmTmp/SEND");
	}
	if (-d "$devShmTmp/REC") {
		system("rm -rf $devShmTmp/REC");	
	}
	my $dbSEND = RocksDB->new("$devShmTmp/SEND", { create_if_missing => 1, db_log_dir => "/dev/null", keep_log_file_num => 1, max_log_file_size => 1, db_stats_log_interval => 9999999, max_bytes_for_level_base => 100000000, target_file_size_base => 20000000});
	my @toBeQueued;
	my $threadConts = {};
	my $threadPut = 1;
	$totalJobs = 0;
	$finishedJobs = 0;
	for (my $i = 1; $i <= $blastThreads; $i++) {
		$threadConts->{$i} = {};
	}	
	foreach my $def (sort {$sortLengths->{$b} <=> $sortLengths->{$a}} keys %$before) {
		unless (exists $after->{$def}) {
			$bad = 1;
			next;
		}
		my $beforePrint = $before->{$def};
		$beforePrint =~ s/N//gi;
		my $afterPrint = $after->{$def};
		$afterPrint =~ s/N//gi;
		my $container = {};
		$container->{'beforePrint'} = $beforePrint;
		$container->{'afterPrint'} = $afterPrint;
		$threadConts->{$threadPut}->{$def} = $container;
		$threadPut++;
		$totalJobs++;
		if ($threadPut > $blastThreads) {
			$threadPut = 1;
		}
		push @toBeQueued, $def;
	}
	foreach my $thrPut (keys %$threadConts) {
		$dbSEND->put("$thrPut", compress($encoder->encode($threadConts->{$thrPut})));
	}
	$dbSEND->compact_range();
	$dbSEND = "";
	print STDERR "THREADING START\n";
	my $threadTime1 = time;
	print STDERR "DONE ENQUEUE JOBS\n";
	for (my $i = 0; $i < $blastThreads; $i++) {
		$goQ->enqueue("GO");
	}
	print STDERR "DONE ENQUEUE GOS\n";
	for (my $i = 0; $i < $blastThreads; $i++) {
		$returnQ->dequeue();
	}
	print STDERR "GOT ALL RETURNS\n";
	sub saveDef {
		my $def = $_[0];
		$saved->{$def} = $after->{$def};
		$savedDef->{$def} = $afterDef->{$def};
	}
	print STDERR "THREADING END\n";


	my $dbREC = RocksDB->new("$devShmTmp/REC", { read_only => 1, db_log_dir => '/dev/null', keep_log_file_num => 1}) or die "can't open DB\n";
	my $denom = scalar(@toBeQueued);
	my $totalDistance = 0;
	my $maxDist = 0;
	my $maxAcc = "";
	foreach my $def (@toBeQueued) {
		my $distance = $dbREC->get("$def");
		if ($distance > $maxDist) {
			$maxDist = $distance;
			$maxAcc = $def;
		}
		$totalDistance += $distance;
		if ($hitTimeLimit) {
			print STDERR "accepted $def because time limit hit\n";
			saveDef($def);
			next;
		}
		my $md5 = md5_hex($after->{$def});
		if (exists $seen->{$def}->{$md5} and $iterateCount > 2) {
			saveDef($def);
			next;
		} elsif ((defined($distance) and $distance < .005 and $iterateCount > 2)) {
			saveDef($def);
			next;
		} elsif (defined($distance) and $distance == 0) {
			saveDef($def);
			next;
		} elsif (defined($distance) and $distance < .01  and $iterateCount > 5) {
			saveDef($def);
			next
		} elsif (defined($distance)){
			$bad = 1;
			next;
		} else {
			print STDERR "$def failed to calculate distance (should never hit this)\n";
			$bad = 1;
			next;
		}
	}
	my $avgDistance = $totalDistance / $denom;
	print STDERR "Average distance to previous $avgDistance, max distance: $maxDist on $maxAcc\n";
	my $savedCount = 0;
	my $savedSeq = 0;
	foreach my $def (keys %$saved) {
		$savedCount++;
		$savedSeq += length($saved->{$def});
	}
	my $threadTime2 = time;
	my $elapsedThread = $threadTime2 - $threadTime1;
	print STDERR "THREADING took $elapsedThread seconds at $blastThreads threads\n";
	print STDERR "Cycle $iterateCount finished\n";
	unless ($bad) {
		print STDERR "Converged on cycle $iterateCount\n";
		last;
	} else {
		print STDERR "Didn't converge on cycle $iterateCount\n";
	}
	my $t2 = time;
	my $elapsed = $t2 - $t1;
	$t1 = $t2;
	print STDERR "$elapsed seconds for cycle $iterateCount\n\n\n";
	$iterateCount++;
	$radius += $radGrowth;
	
}
my $usedAcc = {};
foreach my $def (sort keys %$saved) {
	$usedAcc->{$def} = 1;
	my $seq = $saved->{$def};
	print ">$savedDef->{$def}\n";
	print "$seq\n";
}
open IN, "$in";
my $fasta_file = new FAlite(\*IN); # or any other filehandle
while (my $entry = $fasta_file->nextEntry) {
	my $def = $entry->def;
	my $acc = getAcc($def);
	if (exists $usedAcc->{$acc}) {
		next;
	}
	my $seq = $entry->seq;
	print "$def\n";
	print "$seq\n";
}
close IN;
for (my $i = 0; $i < $blastThreads; $i++) {
	$goQ->enqueue(undef);
}
for (my $i = 0; $i < $blastThreads; $i++) {
	$threads[$i]->join();
}
my $finishTime = time;
my $elapsed = $finishTime - $startTime;
print STDERR "Overall improve time: $elapsed seconds\n";

sub getAcc {
	my $def = $_[0];
	my $acc = "NA";
	if ($def =~ m/$accQR/) {
		$acc = $1;
	} else {
		$def =~ m/contig_(\d+)/;
		$acc = $1;
	}
	return($acc);

}


sub worker {
	use RocksDB;
	use Text::Levenshtein::XS qw(distance);
	use Text::Levenshtein::Damerau::XS qw/xs_edistance/;
	select(STDERR);
	$| = 1;
	my $tmpdir = $_[0];
	my $tid = threads->tid();
	my $encoder = Sereal::Encoder->new();
	my $decoder = Sereal::Decoder->new();
	sub hd {
		return ($_[0] ^ $_[1]) =~ tr/\001-\255//;
	}

	while (1) {
		my $go = $goQ->dequeue();
		lockPrintStderr("got go");
		if (not defined ($go)) {
			my $peek = $goQ->peek();
			if (not defined($peek)) {
				last;
			}
		}
		my $usleep = int(rand(500));
		usleep($usleep);
		my $thrCont;
		{lock($lock);
		my $db = RocksDB->new("$tmpdir/SEND", { read_only => 1, db_log_dir => '/dev/null', keep_log_file_num => 1}) or die "can't open DB\n";
		$thrCont = $decoder->decode(decompress($db->get("$tid")));
		}
		my $localFinish = 0;
		my $localJobsCount = scalar(keys %$thrCont);
		my $distances = {};
		my $t1 = time;
		foreach my $def (keys %$thrCont) {
			my $t2 = time;
			if ($t2 - $t1 > 10) {
				$t1 = $t2;
				lock($finishedJobs);
				lockPrintStderr("finished $finishedJobs/$totalJobs, $localFinish/$localJobsCount done");
			}
			my $container = $thrCont->{$def};
			my $beforePrint = uc($container->{'beforePrint'});
			my $afterPrint = uc($container->{'afterPrint'});
			my $lengthB = length($beforePrint);
			my $lengthA = length($afterPrint);
			my $afterDenom = length($afterPrint);
			my $totalLength = $afterDenom + length($beforePrint);
			my $distance = "";
			if ($beforePrint eq $afterPrint) {
				$distances->{$def} = 0;
			} elsif ($lengthA == $lengthB) {
				my $diff = hd($beforePrint, $afterPrint);
				$distances->{$def} = $diff / $afterDenom;
			} elsif ($lengthA <= 5000 and $lengthB <= 5000) {
				my $diff = xs_edistance($beforePrint, $afterPrint);
				if (not defined $diff) {
					$diff = $afterDenom;
				}
				$distances->{$def} = $diff / $afterDenom;
			} else {
				my $beforeStr = ">Before\n$beforePrint\n";
				my $afterStr = ">After\n$afterPrint\n";
				{lockPrintStderr("long blast $def");
					lock($lock);
					open BEF,">$tmpdir/$tid.Query.fa";
					print BEF "$beforeStr";
					close BEF;
					open AFT, ">$tmpdir/$tid.Subject.fa";
					print AFT "$afterStr";
					close AFT;				
				}
				my $blastReturn = `blastn -task megablast -word_size 50 -query $tmpdir/$tid.Query.fa -subject $tmpdir/$tid.Subject.fa -outfmt "6 std gaps" -max_target_seqs 1 -culling_limit 1 -qcov_hsp_perc 1 | grep Before | grep After`;    
				chomp $blastReturn;
				my $distance;
				sub calcBlastDistance {
					my $distance = $_[0];
					my $blastReturn = $_[1];
					unless ($blastReturn =~ m/Before.*After/) {
						return undef;
					}
					my @lines = split /\n/, $blastReturn;
					foreach my $line (@lines) {
						my @parts = split /\t/, $line;
						$distance -= int((($parts[2] / 100) * ($parts[3] * 2)));
						$distance += $parts[12];
					}
					return $distance;
				}
				$distance = calcBlastDistance($totalLength, $blastReturn);
				if (not defined $distance or $distance < 0) {
					$blastReturn = `blastn -task megablast -word_size 50 -query $tmpdir/$tid.Query.fa -subject $tmpdir/$tid.Subject.fa -outfmt "6 std gaps" -max_hsps 1 | grep Before | grep After`;
					chomp $blastReturn;
					$distance = calcBlastDistance($totalLength, $blastReturn);
					if ($distance < 0) {
						$distance = $afterDenom;
						lock($print);
						print STDERR "negative distance for $def\n";
					} else {
						lock($print);
						print STDERR "Levenshtein for $def\n";
						$distance = distance($beforePrint, $afterPrint);
					}
				}
				unless (defined ($distance)) {
					$distance = $afterDenom;
				}
				$distances->{$def} = $distance / $afterDenom;
			}
			{lock($finishedJobs);
			$localFinish++;
			$finishedJobs++;
			}
		}
		{lock($lock);
		my $dbREC = RocksDB->new("$tmpdir/REC", { create_if_missing => 1, db_log_dir => "/dev/null", keep_log_file_num => 1, max_log_file_size => 1, db_stats_log_interval => 9999999, max_bytes_for_level_base => 100000000, target_file_size_base => 20000000});
		foreach my $def (keys %$distances) {
			$dbREC->put("$def", "$distances->{$def}");
		}
		$dbREC->compact_range();
		$returnQ->enqueue($tid);
		}
		
	}
	return();
}
sub lockPrintStderr {
	lock($print);
	my $tid = threads->tid();
	my $printout = $_[0];
	print STDERR "THREAD $tid - $printout\n";
}
