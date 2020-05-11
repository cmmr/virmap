#!/usr/bin/env perl


use warnings;
use strict;
use threads;
use threads::shared;
use Thread::Queue;
use RocksDB;
use Compress::Zstd qw(compress decompress compress_mt);
use FAlite;
use File::Temp qw/tempdir/;
use Sereal;
use FAlite;
use IPC::Open2;
use Time::HiRes qw(time usleep);


my $t1 = time();
my $exactRatio = .03;
my $zipperRatio = .03;
my $maxOverlapError = .04;
my $minOverlap = 100;
my $evalQ = Thread::Queue->new();
my $returnQ = Thread::Queue->new();
my $goQ = Thread::Queue->new();
my $lock :shared;
my $THREADS = $ARGV[1];
my $printedContig = 0;
my $seed = 27;
my $encoder = Sereal::Encoder->new();
my $decoder = Sereal::Decoder->new();
my $mer :shared;


if ($THREADS > 12) {
	$THREADS = 12;
}


my $sortThreads = $THREADS;
if ($sortThreads > 8) {
	$sortThreads = 8;
}
my $tmpdir = tempdir( DIR => "/dev/shm", CLEANUP => 1);
if (not $THREADS =~ m/^\d+$/ or $THREADS < 4) {
	$THREADS = 4;
}

if (defined $ARGV[2] and $ARGV[2] eq "strict") {
	$seed = 41;
	$exactRatio = .02;
	$zipperRatio = .02;
	$maxOverlapError = .03;
}
if (defined $ARGV[2] and $ARGV[2] eq "init") {
	$seed = 51;
	$exactRatio = .04;
	$zipperRatio = .04;
	$maxOverlapError = .05;
}
open IN, "cat $ARGV[0] |";
open OUT, ">$tmpdir/easyMerge.query.fa";
my $fasta_file = new FAlite(\*IN); # or any other filehandle
my $noChange = 0;
my $different = 0;
my $total = 0;
my $initialNoRev = 0;
while (my $entry = $fasta_file->nextEntry) {
	$total++;
	my $def = $entry->def;
	if ($def =~ m/;rev=1.*;rev=1/) {
		die "double double reverse reverse bug\n";
	}
	unless ($def =~ m/;rev=1/) {
		$initialNoRev++;
	}
	if ($def =~ m/;noChange=1/) {
		$noChange++;
		next;
	}
	$different++;
	my $seq = $entry->seq;
	print OUT "$def\n";
	print OUT "$seq\n";
}
close IN;
close OUT;
unless ($different) {
	open IN, "cat $ARGV[0] |";
		
	$fasta_file = new FAlite(\*IN);
	while (my $entry = $fasta_file->nextEntry) {
		my $def = $entry->def;
		if ($def =~ m/;rev=1/) {
			next;
		}
		$printedContig++;
		my $seq = $entry->seq;
		print "$def\n";
		print "$seq\n";
	}
	close IN;
	exit();
}
my @threads;
for (my $i = 0; $i < $THREADS; $i++) {
	$threads[$i] = threads->create(\&worker, $tmpdir);
}
system("makeblastdb -in $ARGV[0] -out $tmpdir/easyMergeDb -dbtype nucl 1>&2");
my $blastout = `blastn -word_size $seed -query $tmpdir/easyMerge.query.fa -db $tmpdir/easyMergeDb -outfmt 6 -strand plus -ungapped -max_target_seqs 10000000`;
my $dbHL = RocksDB->new("$tmpdir/HL", { create_if_missing => 1, db_log_dir => "/dev/null", keep_log_file_num => 1, max_log_file_size => 1, db_stats_log_interval => 9999999, max_bytes_for_level_base => 100000000, target_file_size_base => 20000000, allow_mmap_reads => 'true'});
my $dbRO = RocksDB->new("$tmpdir/RO", { create_if_missing => 1, db_log_dir => "/dev/null", keep_log_file_num => 1, max_log_file_size => 1, db_stats_log_interval => 9999999, max_bytes_for_level_base => 100000000, target_file_size_base => 20000000, allow_mmap_reads => 'true'});
my $dbAUX = RocksDB->new("$tmpdir/AUX", { create_if_missing => 1, db_log_dir => "/dev/null", keep_log_file_num => 1, max_log_file_size => 1, db_stats_log_interval => 9999999, max_bytes_for_level_base => 100000000, target_file_size_base => 20000000, allow_mmap_reads => 'true'});
my @lines = split /\n/, $blastout;
$blastout = "";
my $results = {};
my $highestLeft = {};
my $rangeOrders = {};
my $hits = {};
my $sizes = {};
my $identities = {};
my $hspLength = {};
my $instant = 0;
my $seqs = {};
my $merged = {};

my $newLines = ();
my $reinsert = 0;
foreach my $rawline (@lines) {
	my $line = [split /\t/, $rawline];
	my @oldParts = @$line;
	push @{$newLines}, \@oldParts;
	unless ($line->[1] =~ m/noChange=1/) {
		next;
	}
	my $query = $line->[0];
	my $subject = $line->[1];
	$line->[0] = $subject;
	$line->[1] = $query;
	my $Qstart = $line->[6];
	my $Qend = $line->[7];
	my $Sstart = $line->[8];
	my $Send = $line->[9];
	$line->[6] = $Sstart;
	$line->[7] = $Send;
	$line->[8] = $Qstart;
	$line->[9] = $Qend;
	$reinsert++;
	push @{$newLines}, $line;
}
@lines = sort {$a->[0] cmp $b->[0] || $a->[6] <=> $b->[6]} @{$newLines};
undef($newLines);	
my $noInstant = {};
foreach my $line (@lines) {
	if ($line->[0] eq $line->[1]) {
		next;
	}
	if ($line->[0] =~ m/(.*);rev=1/) {
		if ($line->[1] eq $1) {
			next;
		}
	}
	if ($line->[1] =~ m/(.*);rev=1/) {
		if ($line->[0] eq $1) {
			next;
		}
	}
	if ($line->[9] < $line->[8]) {
		next;
	}
	my $alignLength = $line->[3];
	my $mismatches = $line->[4];
	if ($line->[0] =~ m/^contig.*length=(\d+)/ and not $line->[1] =~ m/rev=1/ and not exists $noInstant->{$line->[1]}) {
		my $length = $1;
		if ($line->[2] > 99 and $line->[3] > .99 * $length) {
			$instant++;
			my $put = "INSTANT $line->[1]";
			if ($line->[0] =~ m/(.*);rev=1/) {
				$merged->{$1} = $put;
			} else {
				$merged->{"$line->[0];rev=1"} = $put;
			}
			$merged->{$line->[0]} = $put;
			if ($line->[1] =~ m/^contig/) {
				$noInstant->{$line->[1]} = 1;
				if ($line->[1] =~ m/(.*);rev=1/) {
					$noInstant->{$1} = 1;
				} else {
					$noInstant->{$line->[1]} = 1
				}
			}
			next;
		}
	} elsif ($line->[1] =~ m/^contig.*length=(\d+)/ and not $line->[0] =~ m/rev=1/ and not exists $noInstant->{$line->[0]}) {
		my $length = $1;
		if ($line->[2] > 99 and $line->[3] > .99 * $length) {
			$instant++;
			my $put = "INSTANT $line->[0]";
			if ($line->[1] =~ m/(.*);rev=1/) {
				$merged->{$1} = $put;
			} else {
				$merged->{"$line->[1];rev=1"} = $put;
			}
			$merged->{$line->[1]} = $put;
			if ($line->[0] =~ m/^contig/) {
				$noInstant->{$line->[0]} = 1;
				if ($line->[0] =~ m/(.*);rev=1/) {
					$noInstant->{$1} = 1;
				} else {
					$noInstant->{$line->[0]} = 1
				}
			}
			next;
		}
	}
	$hits->{$line->[0]}->{$line->[1]}->{"$line->[6]..$line->[7]"}->{"$line->[8]..$line->[9]"} = $line->[11];
	$identities->{$line->[0]}->{$line->[1]} += ($line->[11]);
}
print STDERR "$instant contigs instantly absorbed\n";
foreach my $line (@lines) {
	if ($line->[0] eq $line->[1]) {
		next;
	}
	if ($line->[9] < $line->[8]) {
		next;
	}
	if (exists $merged->{$line->[0]} or exists $merged->{$line->[1]}) {
		next;
	}
	unless (exists $hits->{$line->[1]}->{$line->[0]}->{"$line->[8]..$line->[9]"}->{"$line->[6]..$line->[7]"} and exists ($hits->{$line->[0]}->{$line->[1]}->{"$line->[6]..$line->[7]"}->{"$line->[8]..$line->[9]"})) {
		next;
	}
	$results->{$line->[0]}->{$line->[1]}->{"$line->[6]..$line->[7]"}->{"$line->[8]..$line->[9]"} = $line->[11];
	unless (exists $rangeOrders->{$line->[0]}->{$line->[1]}) {
		$rangeOrders->{$line->[0]}->{$line->[1]} = ();
		push @{$rangeOrders->{$line->[0]}->{$line->[1]}}, "$line->[6]..$line->[7]";
	} else {
		push @{$rangeOrders->{$line->[0]}->{$line->[1]}}, "$line->[6]..$line->[7]";
	}
	unless (exists $highestLeft->{$line->[0]}->{$line->[1]}) {
		$highestLeft->{$line->[0]}->{$line->[1]} = $line->[6];
	} elsif ($line->[6] < $highestLeft->{$line->[0]}->{$line->[1]}) {
		$highestLeft->{$line->[0]}->{$line->[1]} = $line->[6];
	}
}
undef($hits);
undef(@lines);
open IN, "cat $ARGV[0] |";
$fasta_file = new FAlite(\*IN); # or any other filehandle
while (my $entry = $fasta_file->nextEntry) {
	my $def = $entry->def;
	$def =~ s/^>//g;
	my $seq = $entry->seq;
	if (exists $merged->{$def}) {
		next;
	}
	unless (exists $results->{$def}) {
		if ($def =~ m/rev=1/) {
			next;
		}
		$merged->{$def} = "NOHITS";
		$merged->{"$def;rev=1"} = "NOHITS";
		my $printDef = $def;
		unless ($printDef =~ m/;noChange=1/) {
			$printDef = "$printDef;noChange=1";
		}
		$printedContig++;
		print ">$printDef\n";
		print "$seq\n";
		next;
	}	
	$seqs->{$def} = $seq;
	$sizes->{$def} = () = $seq =~ /[^Nn]/gi;

}
close IN;
my $evaluated = {};
foreach my $key (keys %$results) {
	if (exists $merged->{$key}) {
		delete $results->{$key};
		delete $sizes->{$key};
		next;
	}
	unless (exists $sizes->{$key}) {
		die "$key isn't in sizes\n";
	}
}
$dbAUX->put("seqs", compress($encoder->encode($seqs)));
$dbAUX->put("sizes", compress($encoder->encode($sizes)));
my @order = sort {$sizes->{$b} <=> $sizes->{$a}} keys %$results;


my $threadIdx = 0;
my $jobs = {};
my $centroidOrderRev = {};
my $centroidOrderFwd = {};
my $centroidIdx = 0;
my $t2 = time();
my $elapsed = $t2 - $t1;
print STDERR "blast search and contig processing: $elapsed seconds\n";
$t1 = $t2;
foreach my $centroid (@order) {
	foreach my $hit (keys %{$results->{$centroid}}) {
		unless (exists $results->{$hit}) {
			delete $results->{$centroid}->{$hit};
		}
	}
	my $isRev = 0;
	if ($centroid =~ m/;rev=1/) {
		$isRev = 1;
	}
	if ($isRev and not exists $centroidOrderRev->{$centroid}) {
		$centroidOrderRev->{$centroid} = $centroidIdx;
	} elsif (not $isRev and not exists $centroidOrderFwd->{$centroid}) {
		$centroidOrderFwd->{$centroid} = $centroidIdx
	}
	$centroidIdx++;
	foreach my $hit (sort {$identities->{$centroid}->{$b} <=> $identities->{$centroid}->{$a} || $sizes->{$b} <=> $sizes->{$a}} keys %{$results->{$centroid}}) {
		if (exists $jobs->{$hit} and exists $jobs->{$hit}->{$centroid}) {
			next;
		}
		$jobs->{$centroid}->{$hit} = $threadIdx;
		$threadIdx++;
	}
}

my $indexPerCentroid = {};
foreach my $centroid (@order) {
	if (not exists $jobs->{$centroid}) {
		if (exists $centroidOrderFwd->{$centroid}) {
			delete $centroidOrderFwd->{$centroid};
		}
		if (exists $centroidOrderRev->{$centroid}) {
			delete $centroidOrderRev->{$centroid};
		}
		next;
	}
	my @centroidOrder = sort {$jobs->{$centroid}->{$a} <=> $jobs->{$centroid}->{$b}} keys %{$jobs->{$centroid}};
	$indexPerCentroid->{$centroid} = \@centroidOrder;
}
my $queuedJobs = 0;
my $queue = {};
while ($queuedJobs < $threadIdx) {
	foreach my $centroid (sort {$centroidOrderFwd->{$a} <=> $centroidOrderFwd->{$b}} keys %$centroidOrderFwd) {
		if (scalar @{$indexPerCentroid->{$centroid}}) {
			my $hit = shift @{$indexPerCentroid->{$centroid}};
			$evalQ->enqueue($queuedJobs);
			$queue->{$queuedJobs} = [$centroid, $hit];
			$queuedJobs++;
		} else {
			delete $centroidOrderFwd->{$centroid};
		}
	}
	foreach my $centroid (sort {$centroidOrderRev->{$a} <=> $centroidOrderRev->{$b}} keys %$centroidOrderRev) {
		if (scalar @{$indexPerCentroid->{$centroid}}) {
			my $hit = shift @{$indexPerCentroid->{$centroid}};
			$evalQ->enqueue($queuedJobs);
			$queue->{$queuedJobs} = [$centroid, $hit];
			$queuedJobs++;
		} else {
			delete $centroidOrderRev->{$centroid};
		}
	}
}
$dbAUX->put("queue", compress($encoder->encode($queue)));
foreach my $head (keys %$highestLeft) {
	my $compressed = compress($encoder->encode($highestLeft->{$head}));
	$dbHL->put("$head", $compressed);
}
undef($highestLeft);
foreach my $head (keys %$rangeOrders) {
	my $compressed = compress($encoder->encode($rangeOrders->{$head}));
	$dbRO->put("$head", $compressed);
}
undef($rangeOrders);
$dbHL->compact_range;
$dbRO->compact_range;
$dbAUX->compact_range;
undef($dbHL);
undef($dbRO);
undef($dbAUX);
$t2 = time();
$elapsed = $t2 - $t1;
print STDERR "setting up threads and database structures: $elapsed seconds\n";
$t1 = $t2;
print STDERR "total jobs = $threadIdx\n";
undef($results);
for (my $i = 0; $i < $THREADS; $i++) {
	$evalQ->enqueue(undef);
}
for (my $i = 0; $i < $THREADS; $i++) {
	$goQ->enqueue("GO");
}
my $undefCount = 0;
my $returnedJobs = {};
my $buffer = {};
my $lowestSeen = 0;
my $wastedTime = 0;
my $goodTime = 0;
my $mergedCount = 0;
my $deleted = {};
my $timeDelete = time;
while (1) {
	my $evalQCount = $evalQ->pending();
	my $returnQCount = $returnQ->pending();
	if ($undefCount == $THREADS and $evalQCount == 0 and $returnQCount == 0 and scalar(keys %$buffer) == 0) {
		last;
	}
	my $lowest = "inf";
	if ($returnQCount) {	
		my @returnJobs = $returnQ->extract(0, $returnQCount);
		foreach my $jobReturn (@returnJobs) {
			unless (defined $jobReturn and $undefCount < $THREADS) {
				$undefCount++;
				next;
			}
			$jobReturn = $decoder->decode($jobReturn);
			$buffer->{$jobReturn->[0]} = $jobReturn;
			if ($jobReturn->[0] < $lowest) {
				
				$lowest = $jobReturn->[0];
			}
		}
	} else {
		usleep(10);
		next;
	}
	foreach my $jobIdx (sort {$a <=> $b} keys %$buffer) {
		unless ($jobIdx == $lowestSeen) {
			last;
		}
		$lowestSeen++;
		my $return = $buffer->{$jobIdx};
		delete $buffer->{$jobIdx};
		my $fail = $return->[1];
		my $centroid = $return->[2];
		my $hit = $return->[3];
		my $messages = $return->[4];
		if (exists $merged->{$centroid} or exists $merged->{$hit}) {
			$wastedTime += $return->[9];
			next;
		}
		if ($fail) {
			$goodTime += $return->[9];
			next;
		}
		my $head = $return->[5];
		my $killList = {};
		$killList->{$centroid} = 1;
		$killList->{$hit} = 1;
		if ($centroid =~ m/(.*);rev=1/) {
			$killList->{$1} = 1;
		} else {
			$killList->{"$centroid;rev=1"} = 1;
		}
		if ($hit =~ m/(.*);rev=1/) {
			$killList->{$1} = 1;
		} else {
			$killList->{"$hit;rev=1"} = 1;
		}
		foreach my $kill (keys %$killList) {
			$merged->{$kill} = 1;
		}
		my $seq = $return->[6];
		$head =~ s/;noChange=1//g;
		print ">$head\n";
		print "$seq\n";
		foreach my $message (@$messages) {
			print STDERR "MERGE: $message\n";
		}
		print STDERR "MERGE: printed as $head\n";
		print STDERR "GOOD MERGE: used $return->[8]\n";
		$goodTime += $return->[9];
		$printedContig++;
		$mergedCount++;
	}	
	my $timeAfter = time;
	if ($timeAfter - $timeDelete > 60) {
		lock($evalQ);
		my $pending = $evalQ->pending();
		if ($pending) {
			my @remain = $evalQ->extract(0, $pending);
			foreach my $jobRemain (@remain) {
				if (defined $jobRemain and (exists $merged->{$queue->{$jobRemain}->[0]} or exists $merged->{$queue->{$jobRemain}->[1]})) {
					my $failResult = [$jobRemain, 1, $queue->{$jobRemain}->[0], $queue->{$jobRemain}->[1], ["DELETED"], "", "", 0, "0 seconds", 0];
					$buffer->{$jobRemain} = $failResult;
				} else {
					$evalQ->enqueue($jobRemain);
				}
			}
		}
		$timeDelete = time;
	}
}
$t2 = time();
$elapsed = $t2 - $t1;
print STDERR "$mergedCount contigs merged; threading portion took $elapsed seconds\n";

sub worker {
	select(STDERR);
	$| = 1;	
	my $encoder = Sereal::Encoder->new();
	my $decoder = Sereal::Decoder->new();
	my $tmpdir = $_[0];
	my $go = $goQ->dequeue();
	if ($go eq "NO") {
		return(0);
	}
	select(STDERR);
	$| = 1;
	my $dbHL = RocksDB->new("$tmpdir/HL", { read_only => 1, db_log_dir => '/dev/null', keep_log_file_num => 1, allow_mmap_reads => 'true'}) or die "can't open DB\n";
	my $dbRO = RocksDB->new("$tmpdir/RO", { read_only => 1, db_log_dir => '/dev/null', keep_log_file_num => 1, allow_mmap_reads => 'true'}) or die "can't open DB\n";
	my $dbAUX;
	{lock($lock);
	$dbAUX = RocksDB->new("$tmpdir/AUX", { read_only => 1, db_log_dir => '/dev/null', keep_log_file_num => 1, allow_mmap_reads => 'true'}) or die "can't open DB\n";
	}
	my $queue = $decoder->decode(decompress($dbAUX->get("queue")));
	my $highestLeft = {};
	my $rangeOrders = {};
	my $localSizes = $decoder->decode(decompress($dbAUX->get("sizes")));
	my $localSeqs = $decoder->decode(decompress($dbAUX->get("seqs")));
	my $tid = threads->tid();
	while (1) {
		my $jobIdx = $evalQ->dequeue();
		my $t1 = time();
		unless (defined($jobIdx)) {
			lock($evalQ);
			my $peek = $evalQ->peek();
			if (not defined $peek) {
				last;
			} else {
				redo;
			}	
		}
		my $job = $queue->{$jobIdx};
		my $centroid = $job->[0];
		my $hit = $job->[1];
		my $messages = [];
		my $seconds = 0;
		my $timeRaw = 0;
		my $failResult = [$jobIdx, 1, $centroid, $hit, $messages, "", "", $tid];
		my $mergeRatio = $exactRatio; 
		my $mergeZipper = $zipperRatio;
		
		if ($centroid =~ m/^contig/ and $hit =~ m/^contig/) {
			if ($mergeRatio > .01) {
				$mergeRatio = .01;
			}
			if ($mergeZipper > .01) {
				$mergeZipper = .01;
			}
		}
		my $left = $centroid;
		my $right = $hit;
		unless (exists $highestLeft->{$centroid}) {
			$highestLeft->{$centroid} = $decoder->decode(decompress($dbHL->get("$centroid")));
			$rangeOrders->{$centroid} = $decoder->decode(decompress($dbRO->get("$centroid")));
		}
		unless (exists $highestLeft->{$hit}) {
			$highestLeft->{$hit} = $decoder->decode(decompress($dbHL->get("$hit")));
			$rangeOrders->{$hit} = $decoder->decode(decompress($dbRO->get("$hit")));
		}
		if ($highestLeft->{$centroid}->{$hit} < $highestLeft->{$hit}->{$centroid}) {
			$left = $hit;
			$right = $centroid;
		}
		my $leftIsContig = 0;
		my $rightIsContig = 0;
		if ($left =~ m/^contig/) {
			$leftIsContig = 1;
		}
		if ($right =~ m/^contig/) {
			$rightIsContig = 1;
		}
		my $leftSequence = $localSeqs->{$left};
		my $origLeftLength = length($leftSequence);
		my $rightSequence = $localSeqs->{$right};
		my $origRightLength = length($rightSequence);
		my $initialOffset = $highestLeft->{$left}->{$right} - $highestLeft->{$right}->{$left};
		my $leftOffset = 0;
		my $rightOffset = $initialOffset;
		my @leftOrder = @{$rangeOrders->{$left}->{$right}};
		my @rightOrder = @{$rangeOrders->{$right}->{$left}};
		my $addN = "N" x $initialOffset;
		$rightSequence = $addN . $rightSequence;
		shift @leftOrder;
		shift @rightOrder;
		my $addedNs = 0;
		unless (scalar(@leftOrder) == scalar(@rightOrder)) {
			push @{$messages}, "hits for $left and $right are not symmetrical";
			my $t2 = time();
			my $elapsed = $t2 - $t1;
			$timeRaw = $elapsed;
			$seconds = "$elapsed seconds";
			push @{$failResult}, ($seconds, $timeRaw);
			$returnQ->enqueue($encoder->encode($failResult));
			next;
		}
		my $failIn = 0;
		for (my $i = 0; $i < scalar(@leftOrder); $i++) {
			$leftOrder[$i] =~ m/(\d+)\.\.\d+/;
			my $leftStart = $1;
			$rightOrder[$i] =~ m/(\d+)\.\.\d+/;
			my $rightStart = $1;
			unless ($leftStart and $rightStart) {
				push @{$messages}, "can't get hits for $left vs $right at pos:\n$leftOrder[$i] vs $rightOrder[$i]";
				$failIn = 2;
				last;
			}
			my $diffRightLeft = ($rightStart + $rightOffset) - ($leftStart + $leftOffset);
			my $absAdd = abs($diffRightLeft);
			$addedNs += $absAdd;
			my $negativeLeft = 0;
			my $leftCheckStop = $leftStart + $leftOffset;
			my $leftCheckStart = $leftCheckStop - $absAdd;
			if ($leftCheckStart < 0) {
				$negativeLeft = 1;
				$leftCheckStart = 0;
			}
			my $negativeRight = 0;
			my $rightCheckStop = $rightStart + $rightOffset;
			my $rightCheckStart = $rightCheckStop - $absAdd;
			if ($rightCheckStart < 0) {
				$negativeRight = 1;
				$rightCheckStart = 0;
			}
			my $tt1 = time;
			my $leftCheck = uc(substr($leftSequence, $leftCheckStart, ($leftCheckStop - $leftCheckStart)));
			my $rightCheck = uc(substr($rightSequence, $rightCheckStart, ($rightCheckStop - $rightCheckStart)));
			my $tt2 = time;
			my $leftN = $leftCheck =~ tr/N//;
			my $rightN = $rightCheck =~ tr/N//;
			my $tt3 = time;
			my $subStrTime = $tt2 - $tt1;
			my $checkNTime = $tt3 - $tt2;
			if ($negativeLeft) {
				$leftN = "inf";
			}
			if ($negativeRight) {
				$rightN = "inf";
			}
			if ($diffRightLeft < 0) {
				if (($leftN != 0 or $rightN != 0) and $leftN >= $rightN) {
					if ($leftStart + $leftOffset + $absAdd > length($leftSequence)) {
						$failIn = 1;
						last;
					}
					my $leftBegin = substr($leftSequence, 0, $leftStart + $leftOffset);
					my $leftEnd = substr($leftSequence, ($leftStart + $leftOffset + $absAdd));
					$leftSequence = $leftBegin . $leftEnd;
					$leftOffset -= $absAdd;
				} else  {
					my $tmpStr = "N" x $absAdd;
					if (($rightStart + $rightOffset) > length($rightSequence)) {
						$failIn = 1;
						last;
					}
					my $rightBegin = substr($rightSequence, 0, $rightStart + $rightOffset);
					my $rightEnd = substr($rightSequence, ($rightStart + $rightOffset));
					$rightSequence = $rightBegin . $tmpStr . $rightEnd;
					$rightOffset += $absAdd;
				}
			} elsif ($diffRightLeft > 0) {
				if (($leftN != 0 or $rightN != 0) and $rightN >= $leftN) {
					if (($rightStart + $rightOffset + $absAdd) > length($rightSequence)) {
						$failIn = 1;
						last
					}
					my $rightBegin = substr($rightSequence, 0, $rightStart + $rightOffset);
					my $rightEnd = substr($rightSequence, ($rightStart + $rightOffset + $absAdd));
					$rightSequence = $rightBegin . $rightEnd;
					$rightOffset -= $absAdd;
				} else {
					if (($leftStart + $leftOffset) > length($leftSequence)) {
						$failIn = 1;
						last;
					}
					my $tmpStr = "N" x $absAdd;
					my $leftBegin = substr($leftSequence, 0, $leftStart + $leftOffset);
					my $leftEnd = substr($leftSequence, ($leftStart + $leftOffset));
					$leftSequence = $leftBegin . $tmpStr . $leftEnd;
					$leftOffset += $absAdd;
				}
			}
		}
		if ($addedNs > $origRightLength or $addedNs > $origLeftLength) {
			push @{$messages}, "added N's exceeds length $left vs $right";
			my $t2 = time();
			my $elapsed = $t2 - $t1;
			$seconds = "$elapsed seconds";
			$timeRaw = $elapsed;
			push @{$failResult}, ($seconds, $timeRaw);
			$returnQ->enqueue($encoder->encode($failResult));
			next;
		}
		if ($failIn) {
			if ($failIn == 1) {
				push @{$messages}, "coordinate out of reach $left vs $right";
			}
			my $t2 = time();
			my $elapsed = $t2 - $t1;
			$seconds = "$elapsed seconds";
			$timeRaw = $elapsed;
			push @{$failResult}, ($seconds, $timeRaw);
			$returnQ->enqueue($encoder->encode($failResult));
			next;
		}
		my $lengthRight = length($rightSequence);
		my $lengthLeft = length($leftSequence);
		my $leftPad = $initialOffset;
		my $last = $lengthLeft;
		my $before = substr($leftSequence, 0, $leftPad);
		if ($left =~ m/^contig/ and not $right =~ m/^contig/) {
			$before = lc($before);
		}
		my $after;
		my $rightPad;
		if ($lengthRight > $last) {
			$last = $lengthRight;
			$rightPad = $lengthLeft;
			$after = substr($rightSequence, $rightPad);
			if ($right =~ m/^contig/ and not $left =~ m/^contig/) {
				$after = lc($after);
			}
		} else {
			$rightPad = $lengthRight;
			$after = substr($leftSequence, $rightPad);;
			if ($left =~ m/^contig/ and not $right =~ m/^contig/) {
				$after = lc($after);
			}
		}
		my $leftMiddle = substr($leftSequence, $leftPad, ($rightPad - $leftPad));
		my $rightMiddle = substr($rightSequence, $leftPad, ($rightPad - $leftPad));
		my $limit = length($leftMiddle);
		my $good = 0;
		my $bad = 0;
		my @merged;
		my $fit = 0;
		#Try to determine if sticky ends on contigs actually matter, maybe pad ends with n instead of N?
		my $fromRight = 0;
		my $fromLeft = 0;
		my $lasti = 0;
		my $printed = 0;
		my @leftSequence = split //, $leftMiddle;
		my @rightSequence = split //, $rightMiddle;
		$addedNs = 0;
		for (my $i = 0; $i < $limit; $i++) {	
			$lasti = $i;
			if ($bad / $limit > $maxOverlapError) {
				last;
			}
			if ($addedNs > $origRightLength or $addedNs > $origLeftLength) {
				$bad = $limit;
				last;
			}
			my $leftLetter = "N";
			my $rightLetter = "N";
			if (defined($leftSequence[$i])) {
				$leftLetter = $leftSequence[$i];
			}
			if (defined($rightSequence[$i])) {
				$rightLetter = $rightSequence[$i];
			}
			if ($rightLetter eq "N" and $leftLetter eq "N") {
				push @merged, "N";
				$addedNs++;
			} elsif ($rightLetter eq "n" and $leftLetter eq "n") {
				push @merged, "n";
				$addedNs++;
			} elsif (uc($rightLetter) eq "N" and uc($leftLetter) eq "N") {
				push @merged, "n";
				$addedNs++;
			} elsif (uc($rightLetter) ne "N" and uc($leftLetter) eq "N") {
				$fit++;
				if ($rightIsContig and not $leftIsContig) {
					$fromRight++;
					$rightLetter = lc($rightLetter);
				}
				push @merged, $rightLetter;
				$addedNs++;
			} elsif (uc($rightLetter) eq "N" and uc($leftLetter) ne "N") {
				$fit++;
				if ($leftIsContig and not $rightIsContig) {
					$fromLeft++;
					$leftLetter = lc($leftLetter);
				}
				push @merged, $leftLetter;
				$addedNs++;
			} elsif (uc($rightLetter) ne "N" and uc($rightLetter) eq uc($leftLetter)) {
				if (not $rightIsContig and not $leftIsContig) {
					if ($leftLetter =~ m/[ACGT]/ or $rightLetter =~ m/[ACGT]/) {
						$rightLetter = uc($rightLetter);
					}
				} elsif ($rightIsContig and $leftLetter =~ m/[acgt]/) {
					$rightLetter = lc($rightLetter);
				}
				$good++;
				push @merged, $rightLetter;
			} elsif (uc($rightLetter) ne "N" and uc($rightLetter) ne uc($leftLetter)) {
				$bad++;
				push @merged, "n";
			}
		}
		my $totalOverlap = $bad + $good;
		unless ($totalOverlap) {
			print STDERR "shouldn't have no overlap at this stage $centroid vs $hit\n";
			die;
		}
		my $overlapToFit = "";
		my $weightedOverlap = $totalOverlap + (.5 * $fit);
		my $ratio = $bad / ($weightedOverlap);
		my $ratio2 = $bad / $totalOverlap;
		my $contigLetters = 0;
		my $noFit = 0;
		if ($totalOverlap > .9 * $weightedOverlap) {
			$noFit = 1;
		}
		if (($totalOverlap < $minOverlap) or ($bad and ((($bad / $totalOverlap) > $mergeRatio) and ($noFit or ((($bad / $weightedOverlap) > $mergeZipper) or (($bad / $weightedOverlap) < $mergeZipper and ($bad / $totalOverlap) > $maxOverlapError))))))  {
			push @{$messages}, "failed merge conditions";
			push @{$messages}, "$left vs $right";
			push @{$messages}, "($bad and ((($bad / $totalOverlap) > $mergeRatio) and ((($bad / $weightedOverlap) > $mergeZipper) or (($bad / $weightedOverlap) < $mergeZipper and ($bad / $totalOverlap) > $maxOverlapError)))) ($ratio)";;
			my $t2 = time();
			my $elapsed = $t2 - $t1;
			$seconds = "$elapsed seconds";
			$timeRaw = $elapsed;
			$failResult->[1] = 2;
			push @{$failResult}, ($seconds, $timeRaw);
			$returnQ->enqueue($encoder->encode($failResult));
			#don't merge
		} else {
			push @{$messages}, "($bad and ((($bad / $totalOverlap) > $mergeRatio) ($ratio2) and noFit = $noFit or ((($bad / $weightedOverlap) > $mergeZipper) or (($bad / $weightedOverlap) < $mergeZipper and ($bad / $totalOverlap) > $maxOverlapError))))  ($ratio) [ratio of good events to bad events]";
			push @{$messages}, "Merged $left and $right";
			if ($fromLeft or $fromRight) {
				push @{$messages}, "Previous has contig into N";
			}
			my $newSeq = join "", @merged;
			$newSeq = $before . $newSeq . $after;
			my $print = "";
			my $length = length($newSeq);
			if ($left =~ m/^contig/ and $right =~ m/^contig/) {
				$left =~ m/(contig_\d+);src=([a-z]+);length=(\d+);cov=(\d+\.?\d?)/;
				my $contigPrefix = $1;
				my $leftSrc = $2;
				my $leftLen = $3;
				my $leftCov = $4;
				$right =~ m/src=([a-z]+);length=(\d+);cov=(\d+\.?\d?)/;
				my $rightSrc = $1;
				my $rightLen = $2;
				my $rightCov = $3;
				my $source = "";
				if ($leftSrc eq $rightSrc) {
					$source = $leftSrc
				} else {
					$source = "mixed";
				}
				
				my $newCov = int((($leftLen * $leftCov) + ($rightLen * $rightCov)) / $length);
				
				$print = "$contigPrefix;src=$source;length=$length;cov=$newCov;merged=1;taxId=10239";
			} elsif ((not $left =~ m/rev=1/ and not $right =~ m/rev=1/) or ($left =~ m/rev=1/ and $right =~ m/rev=1/)) {
				if ($left =~ m/^contig/) { 
					$contigLetters = $fromLeft;	
					$print = $right;
				} elsif ($right =~ m/^contig/) {
					$contigLetters = $fromRight;
					$print = $left;
				} elsif ($localSizes->{$left} < $localSizes->{$right}) {
					$print = $right;
				} else {
					$print = $left;
				}
			} elsif ($left =~ m/^contig/ and not $right =~ m/^contig/) {
				$contigLetters = $fromLeft;
				$print = $right;
			} elsif (not $left =~ m/^contig/ and $right =~ m/^contig/) {
				$contigLetters = $fromRight;
				$print = $left;
			} elsif ($left =~ m/rev=1/) {
				$print = $right;
			} elsif ($right =~ m/rev=1/) {
				$print = $left;
			} else {
				lock($lock);
				print STDERR "FAILED all conditions\n";
				print STDERR "left = $left\n";
				print STDERR "right = $right\n";
				die;
			}
			if ($print =~ m/;rev=1/) {
				$print =~ s/;rev=1//g;
				$newSeq = reverse($newSeq);
				$newSeq =~ tr/ACGT/TGCA/;
				$newSeq =~ tr/acgt/tgca/;
			}
			if ($print =~ m/length=(\d+)/) {
				$print =~ s/length=\d+/length=$length/;
			}
			my $t2 = time();
			my $elapsed = $t2 - $t1;
			$seconds = "$elapsed seconds";
			$timeRaw = $elapsed;
			$returnQ->enqueue($encoder->encode([$jobIdx, 0, $centroid, $hit, $messages, $print, $newSeq, $tid, $seconds, $timeRaw]));
			$print =~ s/;noChange=1//g;
			
		}
	}		
	lock($lock);
	$returnQ->enqueue(undef);
	return();
}

foreach my $thr (@threads) {
	$thr->join();
}

print STDERR "Time wasted = $wastedTime\n";
print STDERR "Good time = $goodTime\n";


foreach my $def (sort keys %$seqs) {
	unless (exists $merged->{$def}) {
		if ($def =~ m/rev=1/) {
			next;
		}
		my $seq = $seqs->{$def};;
		if (not $def =~ m/;noChange=1/) {
			$def .= ";noChange=1";
		}
		$printedContig++;
		print ">$def\n";
		print "$seq\n";
	} 
}

print STDERR "started with $initialNoRev non reversed contigs, printed $printedContig contigs, $instant instantly removed\n";













