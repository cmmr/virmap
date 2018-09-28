#!/usr/bin/env perl


use warnings;
use strict;
use threads;
use threads::shared;
use Thread::Queue qw( );
use FAlite;
use English;
use Cpanel::JSON::XS;
use Compress::Zstd qw(compress decompress compress_mt);
use Scalar::Util qw(reftype);
use File::Temp qw/tempdir/;
use RocksDB;
use Sereal;
use Statistics::Basic qw(:all);



my $tmpdir = tempdir( DIR => "/dev/shm", CLEANUP => 1);








my $verb = shift @ARGV;
my $infoFloor = 0;
unless ($verb eq "filter" or $verb eq "classify" or $verb eq "entropy") {
	die;
}
my $filter = 0;
if ($verb eq "filter") {
	$filter = 1;
}
my $entropy = 0;
if ($verb eq "entropy") {
	$entropy = 1;
}
if ($verb eq "classify") {
	$infoFloor = shift @ARGV;
}
unless ($infoFloor =~ m/^\d+$/) {
	print STDERR "infoFloor must be an integer\n";
	die;
}


my $finished :shared;
$finished = 0;
my $total :shared;
my $lock :shared = 0;
my $THREADS = $ARGV[4];
my $q = Thread::Queue->new();
my $startQ = Thread::Queue->new();
my $groupQ = Thread::Queue->new();
my $encoder = Sereal::Encoder->new();
my $decoder = Sereal::Decoder->new();
my $lookUpTree = {
	"11320" => ["1"],
	"95342" => ["2", "1"],
	"11983" => ["2", "1"]
};
my $groups = {
	"virus" => "10239",
	"others" => "28384",
	"viroids" => "12884",
	"unclassified" => "12908",
	"cellularOrganisms" => "131567"
};

my $whiteListFile;

if ($entropy) {
	$THREADS = 0;
	$whiteListFile = $ARGV[2];
} 

my $whiteList = {};
if ($whiteListFile) {
	if ($entropy) {
		my ($parents, $names, $ranks, $children) = makeTaxonomy($ARGV[0]);
	}
	open IN, "$whiteListFile";
	while (my $line = <IN>) {
		chomp $line;
		$whiteList->{$line} = 1;
	}
	close IN;

}
$whiteList = shared_clone($whiteList);

if (not $entropy and $THREADS < 4) {
	$THREADS = 4;
} elsif ($entropy) {
	$THREADS = 0;
}
my @threads;
for (my $i = 0; $i < $THREADS; $i++) {
	$threads[$i] = threads->create(\&worker, $tmpdir, $whiteList, $filter, $infoFloor)
}
my $RocksDBthreads = $THREADS;
if ($RocksDBthreads > 8) {
	$RocksDBthreads = 8;
}
my $sortThreads = $THREADS;
if ($sortThreads > 8) {
	$sortThreads = 8;
}
open IN, "$ARGV[1]";
my $fasta_file = new FAlite(\*IN);
my $seqCount = 0;

my $dbTAX = RocksDB->new("$tmpdir/TAX", { create_if_missing => 1, db_log_dir => "/dev/null", keep_log_file_num => 1, max_log_file_size => 1, db_stats_log_interval => 9999999, max_bytes_for_level_base => 100000000, target_file_size_base => 20000000});
$dbTAX->put("groups", compress($encoder->encode($groups)));


my $scaffolds = {};
while(my $entry = $fasta_file->nextEntry) {
	my $head = $entry->def;
	$head =~ s/^>//g;
	$seqCount++;
	my $seq = $entry->seq;
	if ($entropy) {
		$head =~ m/taxId=(\d+)/;
		my $taxId = $1;
		my $print = 1;
		my $printErr = "";
		if (not defined $taxId) {
			$printErr = "WARNING: $head has no taxaId";
		} else {
			my $seqEnt = entropy($seq);
			if ($seqEnt <= 5) {
				$print = 0;
				$printErr = "FILTERED: $head has entropy of $seqEnt";
			}
		}
		if ($print) {
			print ">$head\n";
			print "$seq\n";
		}
		if ($printErr) {
			print STDERR "$printErr\n";
		}
	} else {
		$dbTAX->put("scaf.$head", $seq);
		$scaffolds->{$head} = 1;
	}
}
close IN;
undef($fasta_file);
if ($entropy) {
	exit();
}

my $dbAA = RocksDB->new("$tmpdir/AA", { create_if_missing => 1, db_log_dir => "/dev/null", keep_log_file_num => 1, max_log_file_size => 1, db_stats_log_interval => 9999999, max_bytes_for_level_base => 100000000, target_file_size_base => 20000000});
my $dbNUC = RocksDB->new("$tmpdir/NUC", { create_if_missing => 1, db_log_dir => "/dev/null", keep_log_file_num => 1, max_log_file_size => 1, db_stats_log_interval => 9999999, max_bytes_for_level_base => 100000000, target_file_size_base => 20000000});
my ($parents, $names, $ranks, $children) = makeTaxonomy($ARGV[0]);


sub readBlastLines {

	my $reference = $_[0];
	my $db = $_[1];
	my $topPerPos = $_[2];
	my $topPerPosPerTaxa = $_[3];
	my $parents = $_[4];
	my $ranks = $_[5];
	my $names = $_[6];
	my $corrections = $_[7];
	my $ok = $_[8];
	my $localScore = $_[9];
	my $overallScores = $_[10];
	my $maxLocal = $_[11];
	my $maxPoss = $_[12];
	my $groups = $_[13];
	my $selfMax = $_[14];
	my $localMax = $_[15];
	my $scale = $_[16];
	my $codons = $_[17];
	my $repPerPosPerTaxa = $_[18];
	my $decoder = $_[19];
	my $otherCont = $_[20];
	my $tmpdir = $_[21];
	my $scaffold = $_[22];
	my $blosum62 = $_[23];
	my $singleMode = $_[24];
	my $lambda = $_[25];
	my $k = $_[26];
	my $expectedBitPerNuc = $_[27];
	my $localDepth = {};
	my $tid = threads->tid();
	my $hit = 0;
	my $goodAfterCycle = {};
	my $lines = $decoder->decode(decompress($db->get($reference)));
	my $maxPerAlignPos = {};
	my $isProt = 0;
	my $totalHitPerPos = {};
	my $taxaPerPos = {};
	my $prevTopPerPos = {};
	if (scalar(keys %$codons)) {
		$isProt = 1;
	}
	for (my $cycle = 0; $cycle < 2; $cycle++) {
		foreach my $line (@$lines) {
			my @parts = split /\t/, $line;
			my $head = $parts[0];
			my $coordAdjust = 0;
			if ($head =~ m/;coords=(\d+)/) {
				$coordAdjust = $1;
				$head =~ s/;coords.*//g;
			}		
			my $evalue = $parts[10];
			if ($evalue > .1) {
				next;
			}
			my $refHit = $parts[1];
			$refHit =~ m/;taxId=(\d+)/;
			my $refTaxa = $1;
			unless (defined $refTaxa) {
				next;
			} elsif ($refTaxa == 0) {
				next;
			}
			if (exists $corrections->{$refTaxa}) {
				$refTaxa = $corrections->{$refTaxa};
			}
			$hit = 1;
			if ($cycle > 0) {
				unless (exists $goodAfterCycle->{$refTaxa}) {
					next;
				}
			}	
			my $queryStart = $parts[6] + $coordAdjust;
			my $queryStop = $parts[7] + $coordAdjust;
			my $frame = 1;
			if ($queryStop < $queryStart) {
				my $temp = $queryStart;
				$frame = -1;
				$queryStart = $queryStop;
				$queryStop = $temp;
			}
			my $bitScore = $parts[11];
			my $bitScorePerAlign = $bitScore / (($queryStop - $queryStart) + 1);
			if (not $isProt and $bitScorePerAlign > $expectedBitPerNuc) {
				$bitScorePerAlign = $expectedBitPerNuc;
			}
			my $gotTop = 0;
			for (my $i = $queryStart; $i <= $queryStop; $i++) {
				unless (exists $totalHitPerPos->{$i}) {
					$totalHitPerPos->{$i} = 0;
				}
				if ($cycle == 0) {
					unless (exists $topPerPos->{$i}) {
						$topPerPosPerTaxa->{$i}->{$refTaxa} = $bitScorePerAlign;
						$topPerPos->{$i} = $bitScorePerAlign;
						$gotTop = 1;
					} elsif ($bitScorePerAlign > $topPerPos->{$i}) {
						$topPerPosPerTaxa->{$i}->{$refTaxa} = $bitScorePerAlign;
						$topPerPos->{$i} = $bitScorePerAlign;
						$gotTop = 1;
						foreach my $checkTaxa (keys %{$topPerPosPerTaxa->{$i}}) {
							if ($topPerPosPerTaxa->{$i}->{$checkTaxa} < (.90 * $bitScorePerAlign) ) {
								delete $topPerPosPerTaxa->{$i}->{$checkTaxa};
							}
						}
					} elsif ((not exists $topPerPosPerTaxa->{$i}->{$refTaxa} and $bitScorePerAlign > (.90 * $topPerPos->{$i})) or (exists $topPerPosPerTaxa->{$i}->{$refTaxa} and $bitScorePerAlign > $topPerPosPerTaxa->{$i}->{$refTaxa})) {
						$topPerPosPerTaxa->{$i}->{$refTaxa} = $bitScorePerAlign;
					}
				} else {
					if (exists $topPerPos->{$i}) {
						my $ratio = ($bitScorePerAlign / $topPerPos->{$i});	
						if ($ratio > 1) {
							lock($lock);
							lockPrintStderr("top per pos $i = $topPerPos->{$i}, bitscoreperslign = $bitScorePerAlign (scaled), ratio = $ratio, should never happen, taxa = $refTaxa");
							
						}
						$taxaPerPos->{$i}->{$refTaxa}++;
						if (not exists $topPerPosPerTaxa->{$i}->{$refTaxa} or $topPerPosPerTaxa->{$i}->{$refTaxa} < $bitScorePerAlign) {
							$topPerPosPerTaxa->{$i}->{$refTaxa} = $bitScorePerAlign;
						}
						my $localScale = $scale;
						if (not $isProt and not exists $otherCont->{$i} and not $singleMode) {
							$localScale = 3;
						}
						my $delta = (($bitScorePerAlign / $topPerPos->{$i}) ** (3 ** $localScale)) / $localScale;
						if ($delta > 1) {
							lockPrintStderr("delta = $delta, should never happen, on taxa $refTaxa");
						}
						#### DB Rep per base per accepted alignmet per taxa ####
						$repPerPosPerTaxa->{$i}->{$refTaxa} += (($bitScorePerAlign / $topPerPos->{$i}) ** (3 ** $localScale)) / $localScale;
						$totalHitPerPos->{$i}++;
					} else {
						lockPrintStderr("Top hit not in second round: $i");
						$topPerPos->{$i} = 0;
						$totalHitPerPos->{$i} = 0;
					}
				}
			}
			if ($gotTop) {
				my $segment = uc(substr($scaffold, $queryStart - 1, ($queryStop - $queryStart) + 1));
				my $score = 0;
				my $segLength = length($segment);
				if ($isProt) {
					if ($frame == -1) {
						$segment = reverse($segment);
						$segment =~ tr/ACGT/TGCA/;
					}
					my @codons = $segment =~ /(...)/g;
					foreach my $codon (@codons) {
						if (exists $codons->{$codon} and exists $blosum62->{$codons->{$codon}}) {
							$score += $blosum62->{$codons->{$codon}};
						} else {
							$score--;
						}
					}
				} else {
					my $nCount = $segment =~ tr/Nn//;
					$score = (2 * ($segLength - $nCount)) - (2 * $nCount);
				}
				my $bestBitPerBase = raw2bit($score, $lambda, $k) / $segLength;
				if (not $isProt and $bestBitPerBase > $expectedBitPerNuc) {
					$bestBitPerBase = $expectedBitPerNuc;
				} elsif ($isProt and $bestBitPerBase < $bitScorePerAlign) {
					#rounding errors
					$bestBitPerBase = $bitScorePerAlign;
				} 
				for (my $i = $queryStart; $i <= $queryStop; $i++) {
					if (not exists $maxPerAlignPos->{$i} or $maxPerAlignPos->{$i} < $bestBitPerBase) {
						$maxPerAlignPos->{$i} = $bestBitPerBase;
					}
				}	
			}
		}
		if ($cycle == 0) {
			foreach my $pos (keys %$topPerPosPerTaxa) {
				foreach my $putTaxa (keys %{$topPerPosPerTaxa->{$pos}}) {
					$goodAfterCycle->{$putTaxa} = 1;	
				}
			}
		}
	}
	my @reps = values %$totalHitPerPos;
	my $minHitProt = "inf";;
	my $avgHit = 1;
	if (scalar @reps) {
		$avgHit = mean(\@reps);
		$minHitProt = .1 * $avgHit;
	}
	my $totalBitScore = 0;
	my $totalDenom = 0;
	my $kept = {};
	foreach my $pos (sort {$a <=> $b} keys %$topPerPos) {
		my $add = $maxPerAlignPos->{$pos};
		unless ($totalHitPerPos->{$pos}) {
			next;
		}
		$kept->{$pos} = 0;
		if (not $isProt or ($avgHit < 10) or ($totalHitPerPos->{$pos} > $minHitProt)) {
			$totalBitScore += $add;
			$totalDenom++;
			$kept->{$pos} = 1;
			${$selfMax} += $add;
			${$localMax} += $add;
		} else {
			if (exists $topPerPosPerTaxa->{$pos}) {
				delete($topPerPosPerTaxa->{$pos});
			}
			if (exists $topPerPos->{$pos}) {
				delete($topPerPos->{$pos});
			}
			if (exists $totalHitPerPos->{$pos}) {
				delete($totalHitPerPos->{$pos});
			}
			if (exists $repPerPosPerTaxa->{$pos}) {
				delete($repPerPosPerTaxa->{$pos});
			}
			next;
		}
		foreach my $taxa (keys %{$topPerPosPerTaxa->{$pos}}) {
			$maxLocal->{$taxa} += $add;
			$maxPoss->{$taxa} += $add;
			$localScore->{$taxa} += $topPerPosPerTaxa->{$pos}->{$taxa};
			$overallScores->{$taxa} += $topPerPosPerTaxa->{$pos}->{$taxa};
		}
	}
	my $avgBestBit = 0;
	if ($totalDenom) {
		$avgBestBit = $totalBitScore / $totalDenom;
	}
	return($hit, $avgBestBit);
}

sub checkWhiteList {
	my $whiteList = $_[0];
	my $checkTaxa = $_[1];
	my $parents = $_[2];
	my $escape = 0;
	while ($checkTaxa != 1) {
		$escape++;
		if ($escape == 100) {
			return(0);
		}
		if (exists $whiteList->{$checkTaxa}) {
			return(1);
		}
		if (exists $parents->{$checkTaxa}) {
			$checkTaxa= $parents->{$checkTaxa};
		}
	}
	return(0);
	
}




sub worker {
	use RocksDB;
	select(STDERR);
	$| = 1;
	my $tid = threads->tid();
	my $tmpdir = $_[0];
	my $whiteList = $_[1];
	my $filter = $_[2];
	my $infoFloor = $_[3];
	my $return = {};
	my $ok = {};
	my $nucLambda = 0.625;
	my $nucK =  0.410;
	my $protLambda = 0.267;
	my $protK = 0.041;
	my $maxAddAA = 0;
	my $maxAddNucl = 0;
	my $aaScale = 3;
	my $nucScale = 1;
	my $expectedBitPerNuc = (2 * ($nucLambda)) / log(2);
	my $e = 2.7182818284590452353602874713527;
	my $pi = 3.1415926535897932384626433832795;
	my $sqrt2 = sqrt(2);
	my $cubrt2 = 2 ** (1/3);
	my $taxaAccuracyBase = 10;
	my $entropyLimit = 5.5;
	my $filterRatio = .99;
	my $rescueRatio = .98;
	my $startNow = $startQ->dequeue();
	my $encoder = Sereal::Encoder->new();
	my $decoder = Sereal::Decoder->new();
	my $dbAA = RocksDB->new("$tmpdir/AA", { read_only => 1, db_log_dir => '/dev/null', keep_log_file_num => 1}) or die "can't open DB\n";
	my $dbNUC = RocksDB->new("$tmpdir/NUC", { read_only => 1, db_log_dir => '/dev/null', keep_log_file_num => 1}) or die "can't open DB\n";
	my $dbTAX = RocksDB->new("$tmpdir/TAX", { read_only => 1, db_log_dir => '/dev/null', keep_log_file_num => 1}) or die "can't open DB\n";
	my $groups = $decoder->decode(decompress($dbTAX->get("groups")));
	my $singleMode = $dbTAX->get("singleMode");
	my $suppressors = {
		"others" => 1,
		"unclassified" => 1
	};
	my $codons = {
		TTT => "F", TTC => "F", TTA => "L", TTG => "L",
		TCT => "S", TCC => "S", TCA => "S", TCG => "S",
		TAT => "Y", TAC => "Y", TAA => "*", TAG => "*",
		TGT => "C", TGC => "C", TGA => "*", TGG => "W",
		CTT => "L", CTC => "L", CTA => "L", CTG => "L",
		CCT => "P", CCC => "P", CCA => "P", CCG => "P",
		CAT => "H", CAC => "H", CAA => "Q", CAG => "Q",
		CGT => "R", CGC => "R", CGA => "R", CGG => "R",
		ATT => "I", ATC => "I", ATA => "I", ATG => "M",
		ACT => "T", ACC => "T", ACA => "T", ACG => "T",
		AAT => "N", AAC => "N", AAA => "K", AAG => "K",
		AGT => "S", AGC => "S", AGA => "R", AGG => "R",
		GTT => "V", GTC => "V", GTA => "V", GTG => "V",
		GCT => "A", GCC => "A", GCA => "A", GCG => "A",
		GAT => "D", GAC => "D", GAA => "E", GAG => "E",
		GGT => "G", GGC => "G", GGA => "G", GGG => "G",
	};
	my $blosum62 = {
		A => 4,
		R => 5,
		N => 6,
		D => 6,
		C => 9,
		Q => 5,
		E => 5,
		G => 6,
		H => 8,
		I => 4,
		L => 4,
		K => 5,
		M => 5,
		F => 6,
		P => 7,
		S => 4,
		T => 5,
		W => 11,
		Y => 7,
		V => 4,
		X => -1,
		"*" => 1
	};



	
	while (1) {
		my $reference = $q->dequeue();
		unless (defined $reference) {
			my $peek = $q->peek();
			if (not defined $peek) {
				last;
			} else {
				redo;
			}
		}
		lockPrintStderr("GOT $reference");
		my $isContig = 0;
		if ($reference =~ m/^contig/) {
			$isContig = 1;
		}
		my $fromProt = 0;
		if ($reference =~ m/NoNucHits/) {
			$fromProt = 1;
		}
		my $mixedSources = 0;
		$reference =~ m/taxId=(\d+)/;
		my $refTaxId = $1;
		unless ($refTaxId) {
			lockPrintStderr("$reference has no detectable taxId");
		}
		my $parents = $decoder->decode(decompress($dbTAX->get("par.$reference")));;
		my $ranks = $decoder->decode(decompress($dbTAX->get("ran.$reference")));;
		my $names = $decoder->decode(decompress($dbTAX->get("nam.$reference")));
		my $taxa2Group = $decoder->decode(decompress($dbTAX->get("t2g.$reference")));
		my $corrections = $decoder->decode(decompress($dbTAX->get("cor.$reference")));
		my $whiteListed = checkWhiteList($whiteList, $refTaxId, $parents);
		my $scaffold = $dbTAX->get("scaf.$reference");
		my $scafEnt = entropy($scaffold);
		if ($filter and $scafEnt <= 5) {
			lockPrintStderr("FILTERED: $reference is low complexity $scafEnt bit per triplet");
			next;
		}
		my $lcCount = () = $scaffold =~ /[acgt]/g;	
		if ($lcCount > .2 * length($scaffold)) {
			$mixedSources = 1;
		}
		my $topPerPosAA = {};
		my $topPerPosNucl = {};
		my $topPerPosPerTaxaAA = {};
		my $topPerPosPerTaxaNucl = {};
		my $repPerPosPerTaxaAA = {};
		my $repPerPosPerTaxaNucl = {};
		my $repPerPosPerTaxa = {};
		my $overallScoresAA = {};
		my $overallScoresNucl = {};
		my $overallScores = {};
		my $maxPossAA = {};
		my $maxPossNucl = {};
		my $maxPoss = {};
		my $selfMaxAA = 0;
		my $selfMaxNucl = 0;
		my $selfMax = 0;
		my $nonFilterAddMax = 0;
		my $protHits = 0;
		my $nucHits = 0;
		($protHits, $maxAddAA) = readBlastLines($reference, $dbAA, $topPerPosAA, $topPerPosPerTaxaAA, $parents, $ranks, $names, $corrections, $ok, $overallScoresAA, $overallScores, $maxPossAA, $maxPoss, $groups, \$selfMax, \$selfMaxAA, $aaScale, $codons, $repPerPosPerTaxaAA, $decoder, {}, $tmpdir, $scaffold, $blosum62, $singleMode, $protLambda, $protK, 0);
		($nucHits, $maxAddNucl) = readBlastLines($reference, $dbNUC, $topPerPosNucl, $topPerPosPerTaxaNucl, $parents, $ranks, $names, $corrections, $ok, $overallScoresNucl, $overallScores, $maxPossNucl, $maxPoss, $groups, \$selfMax, \$selfMaxNucl, $nucScale, {}, $repPerPosPerTaxaNucl, $decoder, $topPerPosAA, $tmpdir, $scaffold, {}, $singleMode, $nucLambda, $nucK, $expectedBitPerNuc);
		my @letters = split //, $scaffold;
		unshift @letters, "N";
		my $total = 0;
		my $unal = 0;
		my $badTop = 0;
		my $virusDist = 0;
		my $badTopAA = 0;
		my $badTopNucl = 0;
		my $virusLost = 0;
		my $suppressRescue = 0;
		my $filterUnal = 0;
		my $rescueDist = 0;
		my $groupScorePerPosAA = {};	
		my $groupScorePerPosNucl = {};
		my $groupScorePerPos = {};
		my $avgFail = {};
		my $avgFailOverlap = {};
		my $suppressCount = {};
		my $repPerGroup = {};
		my $avgRepPerPos = 0;
		my $unalVirusRepAdd = 0;
		my $taxaRepAA = {};
		my $taxaRepNucl = {};
		my $taxaRep = {};
		my $posHit = 0;
		my $highUnknown = 0;
		my $denovoUnal = 0;
		my $nonFilterAddAA = 0;
		my $nonFilterAddNucl = 0;
		my $avgBitDensity = 0;
		my $highHitLowDensity = 0;
		if (not $protHits and not $nucHits and $filter) {
			lockPrintStderr("KEPT: NO hits bypass on filter for $reference");
			$return->{$reference} = $scaffold;
			next;
		}
		for (my $pos = 0; $pos < scalar(@letters); $pos++) {
			if (not defined $letters[$pos]) {
				next;
			}
			if ($letters[$pos] eq "N" or $letters[$pos] eq "n") {
				next;
			}
			$total++;
			if (not exists $topPerPosAA->{$pos} and not exists $topPerPosNucl->{$pos}) {
				unless ($letters[$pos] =~ m/[acgt]/) {
					$unal++;
				} else {
					$denovoUnal++;
				}
				if ($filter) {
					$filterUnal++;
				}
				$nonFilterAddNucl += $maxAddNucl;
				$nonFilterAddAA += $maxAddAA;
				next;
			}
			$posHit++;
			my $localGroupRep = {};
			if (exists $repPerPosPerTaxaAA->{$pos}) {
				foreach my $taxa (keys %{$repPerPosPerTaxaAA->{$pos}}) {
					$repPerGroup->{$taxa2Group->{$taxa}} += $repPerPosPerTaxaAA->{$pos}->{$taxa};
					$localGroupRep->{$taxa2Group->{$taxa}} += $repPerPosPerTaxaAA->{$pos}->{$taxa};
					$taxaRepAA->{$taxa} += $repPerPosPerTaxaAA->{$pos}->{$taxa};
					$taxaRep->{$taxa} += $repPerPosPerTaxaAA->{$pos}->{$taxa};
					$avgRepPerPos += $repPerPosPerTaxaAA->{$pos}->{$taxa};	
				}
			}
			if (exists $repPerPosPerTaxaNucl->{$pos}) {
				foreach my $taxa (keys %{$repPerPosPerTaxaNucl->{$pos}}) {
					$repPerGroup->{$taxa2Group->{$taxa}} += $repPerPosPerTaxaNucl->{$pos}->{$taxa};
					$localGroupRep->{$taxa2Group->{$taxa}} += $repPerPosPerTaxaNucl->{$pos}->{$taxa};
					$taxaRep->{$taxa} += $repPerPosPerTaxaNucl->{$pos}->{$taxa};
					$taxaRepNucl->{$taxa} += $repPerPosPerTaxaNucl->{$pos}->{$taxa};
					$avgRepPerPos += $repPerPosPerTaxaNucl->{$pos}->{$taxa};
				}
			}
			if (exists $topPerPosPerTaxaAA->{$pos}) {
				foreach my $taxa (keys %{$topPerPosPerTaxaAA->{$pos}}) {
					if (not exists $groupScorePerPosAA->{$pos}->{$taxa2Group->{$taxa}} or $groupScorePerPosAA->{$pos}->{$taxa2Group->{$taxa}} < $topPerPosPerTaxaAA->{$pos}->{$taxa}) {
						$groupScorePerPosAA->{$pos}->{$taxa2Group->{$taxa}} = $topPerPosPerTaxaAA->{$pos}->{$taxa};
					}
				}
			}
			if (exists $topPerPosPerTaxaNucl->{$pos}) {
				foreach my $taxa (keys %{$topPerPosPerTaxaNucl->{$pos}}) {
					if (not exists $groupScorePerPosNucl->{$pos}->{$taxa2Group->{$taxa}} or $groupScorePerPosNucl->{$pos}->{$taxa2Group->{$taxa}} < $topPerPosPerTaxaNucl->{$pos}->{$taxa}) {
						$groupScorePerPosNucl->{$pos}->{$taxa2Group->{$taxa}} = $topPerPosPerTaxaNucl->{$pos}->{$taxa};
					}
				}	
			}
			foreach my $group (keys %$groups) {
				if (exists $groupScorePerPosNucl->{$pos}->{$group}) {
					$groupScorePerPos->{$pos}->{$group} += $groupScorePerPosNucl->{$pos}->{$group};
				}
				if (exists $groupScorePerPosAA->{$pos}->{$group}) {
					$groupScorePerPos->{$pos}->{$group} += $groupScorePerPosAA->{$pos}->{$group};
				}
			}
			my $topSuppress = {};
			my $topSuppressAA = {};
			my $topSuppressNucl = {};
			my $suppress = 0;
			foreach my $suppressor (keys %$suppressors) {
				if (exists $groupScorePerPosAA->{$pos}->{$suppressor}) {
					$topSuppressAA->{$suppressor} = $groupScorePerPosAA->{$pos}->{$suppressor};
					$topSuppress->{$suppressor} += $groupScorePerPosAA->{$pos}->{$suppressor};
					$suppress = 1;
				}
				if (exists $groupScorePerPosNucl->{$pos}->{$suppressor}) {
					$topSuppressNucl->{$suppressor} = $groupScorePerPosNucl->{$pos}->{$suppressor};
					$topSuppress->{$suppressor} += $groupScorePerPosNucl->{$pos}->{$suppressor};
					$suppress = 1;
				}
			}
			my @bestScores = sort {$groupScorePerPos->{$pos}->{$b} <=> $groupScorePerPos->{$pos}->{$a}} keys %{$groupScorePerPos->{$pos}};
			if (exists $groupScorePerPos->{$pos}->{'virus'}) {
				my $ratios = {};
				foreach my $group (@bestScores) {
					$ratios->{$group} = $groupScorePerPos->{$pos}->{'virus'} / $groupScorePerPos->{$pos}->{$group};
				}
				if ($ratios->{$bestScores[0]} < $filterRatio) {
					if ($suppress) {
						my @bestSuppress = sort {$topSuppress->{$b} <=> $topSuppress->{$a}} keys %$topSuppress;
						if ($ratios->{$bestSuppress[0]} >= $rescueRatio and $ratios->{$bestSuppress[0]} < 1) {
							foreach my $suppressor (@bestSuppress) {
								$repPerGroup->{$suppressor} -= $localGroupRep->{$suppressor};
							}
							$suppressRescue += $ratios->{$bestSuppress[0]};
							$rescueDist++;
							$suppressCount->{$bestSuppress[0]}++;
						} else {
							$suppress = 0;
						}
					} 
					foreach my $group (keys %$ratios) {
						if ($group eq 'virus') {
							next;
						}
						$avgFail->{$group} += $ratios->{$group};
						$avgFailOverlap->{$group}++;
					}
					unless ($suppress) {
						$virusLost++;
						$badTop++;
					}	
				}
			} else {
				$badTop++;
			}
		}
		my $totalScore = 0;
		foreach my $val (values %$topPerPosAA) {
			$totalScore += $val
		}
		foreach my $val (values %$topPerPosNucl) {
			$totalScore += $val
		}
		if ($posHit) {
			$avgBitDensity = $totalScore / $posHit;
			$avgRepPerPos = $avgRepPerPos / $posHit;
		}
		
		$unalVirusRepAdd = $avgRepPerPos * $filterUnal;
		my $unalLimit = .8;
		$highUnknown = (($isContig or $mixedSources) and (($unal + $denovoUnal) / $total > $unalLimit) and $total > 500);
		my $freePass = ($highUnknown or $whiteListed);
		if (not $freePass and not $filter and not $posHit and not $isContig) {
			lockPrintStderr("FILTERED: $reference no hits on pseudo constructed contig");
			next;
		}		
		if (not $freePass and not $filter and ($unal / ($total - $denovoUnal)) > $unalLimit and not ($isContig or $mixedSources)) {
			lockPrintStderr("FILTERED: Contig has <20% of pseudoconstructed segments hitting the database: $reference");
			next;
		}

		if ($unalVirusRepAdd) {
			$repPerGroup->{'virus'} += $unalVirusRepAdd;
		}
		foreach my $group (keys %$repPerGroup) {
			$repPerGroup->{$group} = $repPerGroup->{$group} / $total;
		}
		unless (exists $repPerGroup->{'virus'}) {
			$repPerGroup->{'virus'} = 0;
		}
		my $avgRescue = 0;
		my $rescueStr = "";
		if ($rescueDist) {	
			$avgRescue = $suppressRescue / $rescueDist;
			my @bestSuppress = sort {$suppressCount->{$b} <=> $suppressCount->{$a}} keys %$suppressCount;
			$rescueStr = "best suppressor: $bestSuppress[0] with $suppressCount->{$bestSuppress[0]} based recovered";
		}
		foreach my $group (keys %$avgFail) {
			$avgFail->{$group} = $avgFail->{$group} / $avgFailOverlap->{$group};
		}
		my $remove = 0;
		my $status = "";
		my @repOrder = sort {$repPerGroup->{$b} <=> $repPerGroup->{$a}} keys %$repPerGroup;
		my $bestRepGroup = $repPerGroup->{$repOrder[0]};
		if (not $freePass and ($badTop / $total > .5 or $bestRepGroup == 0 or $repPerGroup->{'virus'} < $bestRepGroup)) {
			if (exists $suppressors->{$repOrder[0]} and $repPerGroup->{'virus'} > (.5 * $bestRepGroup)) {
				$status = "KEPT";
				lockPrintStderr("RECOVERED: $repOrder[0] is a suppressed category, virus is > 50% is it's score: $repPerGroup->{'virus'} vs $bestRepGroup");
			} else {
				$remove = 1;
				$status = "FILTERED";
			}
			if ($repPerGroup->{'virus'} < $bestRepGroup and $status ne "KEPT") {
				lockPrintStderr("FILTERED: $repOrder[0] had a higher score than virus");
			}
		} else {
			$status = "KEPT";
		}	
		unless ($filter) {
			$filterUnal = $unal;
		}
		{lock($lock);
		lockPrintStderr("$status: $reference ($badTop bp/$total bp) bad bases; $virusLost bp virus lost; $filterUnal bp unaligned; $rescueDist bp rescued; avg rescue: $avgRescue; entropy: $scafEnt");
		lockPrintStderr("$status: avgRepPerPos = $avgRepPerPos");
		}
		if ($filter) {
			unless ($remove) {
				$return->{$reference} = $scaffold;
			}
			next;
			
		} elsif ($remove) {
			next;
		}
		#### END FILTER MODE ####		



		#### TAXA ENGINE ####

		my $groupScores = {};		
		my $groupScoresAA = {};
		my $groupSciresNucl = {};
		my $maxAA = 0;
		my $maxNucl = 0;
		my $maxNuclTaxa = 0;
		my $maxAATaxa = 0;
		my $max = 0;
		foreach my $taxa (keys %$overallScores) {
			$groupScores->{$taxa2Group->{$taxa}}->{$taxa} = $overallScores->{$taxa};
			if ($overallScores->{$taxa} > $max) {
				$max = $overallScores->{$taxa};
			}
			if (exists $overallScoresAA->{$taxa} and $overallScoresAA->{$taxa} > $maxAA) {
				$maxAA = $overallScoresAA->{$taxa};
				$maxAATaxa = $taxa;
			}	
			if (exists $overallScoresNucl->{$taxa} and $overallScoresNucl->{$taxa} > $maxNucl) {
				$maxNucl = $overallScoresNucl->{$taxa};
				$maxNuclTaxa = $taxa;
			}
		}
		my $groupScoresSortedOrder = {};
		my $maxes = {};
		
		foreach my $group (keys %$groupScores) {
			my @sorted = sort {$groupScores->{$group}->{$b} <=> $groupScores->{$group}->{$a}} keys %{$groupScores->{$group}};
			$groupScoresSortedOrder->{$group} = \@sorted;
			$maxes->{$group} = $groupScores->{$group}->{$sorted[0]};
			if ($maxes->{$group} > $selfMax) {
				$selfMax = $maxes->{$group};
			}
		}
		unless (exists $maxes->{'virus'}) {
			$maxes->{'virus'} = 0;
		}
		lockPrintStderr("BITDENSITY: $reference = $avgBitDensity");
		if (not $freePass and $maxes->{'virus'} < $infoFloor and ($selfMax < 1.25 * $infoFloor) and not $highUnknown) {
			lockPrintStderr("FILTERED: $reference has a max virus score of $maxes->{'virus'} and a selfMax of $selfMax which is less than the limit of $infoFloor and 125% of $infoFloor");
			next;
		}
		foreach my $taxa (keys %$taxaRep) {
			unless (exists $taxaRep->{$taxa}) {
				lockPrintStderr("WARNING: taxa id: $taxa not present in $reference db calculations");
				$taxaRep->{$taxa} = 0;
			}
			if (exists $taxaRepAA->{$taxa}) {
				$taxaRepAA->{$taxa} = $taxaRepAA->{$taxa} / $total;
			} else {
				$taxaRepAA->{$taxa} = 0;
			}
			if (exists $taxaRepNucl->{$taxa}) {			
				$taxaRepNucl->{$taxa} = $taxaRepNucl->{$taxa} / $total;
			} else {
				$taxaRepNucl->{$taxa} = 0;
			}
			$taxaRep->{$taxa} = $taxaRep->{$taxa} / $total;
		}
		my $minPresence = (.90 * $max);
		if (not $freePass and $maxes->{'virus'} < $minPresence) { 
			lockPrintStderr("lower limit of scoring: $minPresence");
			lockPrintStderr("FILTERED: virus below $minPresence threshold $reference");
			lockBrowseErr("Max scores of groups:", $maxes);
			next;
		}
		{lock($lock);
		foreach my $group (keys %$groupScores) {
			my $checkLimit = 10;
			my @keys = @{$groupScoresSortedOrder->{$group}};
			my $groupMax = $overallScores->{$keys[0]};
			my $groupMemberCount = scalar(@keys);
			if ($groupMemberCount < $checkLimit) {
				$checkLimit = $groupMemberCount;
			}
			unless ($checkLimit) {
				next;
			}
			lockPrintStderr("Top $checkLimit scoring taxa IDs in group $group");
			for (my $i = 0; $i < $checkLimit; $i++) {
				my $printNucl = 0;
				my $printAA = 0;
				unless (exists $overallScoresAA->{$keys[$i]}) {
					$overallScoresAA->{$keys[$i]} = 0;	
				} else {
					$overallScoresAA->{$keys[$i]} = $overallScoresAA->{$keys[$i]};
				}
				unless (exists $overallScoresNucl->{$keys[$i]}) {
					$overallScoresNucl->{$keys[$i]} = 0;
				} else {
					$overallScoresNucl->{$keys[$i]} = $overallScoresNucl->{$keys[$i]};
				}
				if ($i < $checkLimit) {
					my $printScore = int($overallScores->{$keys[$i]});
					my $printRep = int($taxaRep->{$keys[$i]});
					my $printTaxaRepNucl = int($taxaRepNucl->{$keys[$i]});
					my $printTaxaRepAA = int($taxaRepAA->{$keys[$i]});
					my $printScoreAA = int($overallScoresAA->{$keys[$i]});
					my $printScoreNucl = int($overallScoresNucl->{$keys[$i]});
					lockPrintStderr("SCORE: $printScore\tprot: $printScoreAA\tnuc: $printScoreNucl\tdbRep: $printRep\tAA: $printTaxaRepAA\tNucl: $printTaxaRepNucl\t\tTaxID: $keys[$i] - $names->{$keys[$i]}");
				}
			}
			lockPrintStderr("");
		}
		}
		if ($maxAA > $selfMaxAA) {
			$selfMaxAA = $maxAA;
		}
		if ($maxNucl > $selfMaxNucl) {
			$selfMaxNucl = $maxNucl;
		}
		
		my $baseNucl = ($e / ($pi * $sqrt2)) / 10;
		my $baseAA = ($cubrt2 / 4) / 10;
		my $baseBlurNucl = 0;
		my $baseBlurAA = 0;
		my $baseRadiusNucl = 0;
		my $accuracyNucl = 0;
		my $accuracyAdjustNucl = 0;
		my $maxNuclFrac = 0;
		my $approxIdNucl = 0;
		#### LCA radius calculated here ####
		if ($maxNucl) {
			$accuracyNucl = $maxNucl / $selfMaxNucl;
			$approxIdNucl = (((($accuracyNucl * $expectedBitPerNuc) * log(2)) / $nucLambda) + 3) / 5;
			my $denom = log($maxNucl) * log10($maxNucl);
			$baseBlurNucl = (-1 * log($approxIdNucl) / $denom) + $baseNucl;
			$baseRadiusNucl = ($e/log10($maxNucl)) * $baseBlurNucl;
			$accuracyAdjustNucl = (-1 * log10($approxIdNucl) / log10(100 * $baseBlurNucl)) * $baseBlurNucl;
			$maxNuclFrac = $maxNucl / ($maxNucl + $maxAA);
		}
		my $baseRadiusAA = 0;
		my $accuracyAA = 0;
		my $accuracyAdjustAA = 0;
		my $maxAAFrac = 0;
		my $approxIdAA = 0;
		if ($maxAA) {
			$accuracyAA = $maxAA / $selfMaxAA;
			my $denom = log($maxAA) * log10($maxAA);
			$baseBlurAA = (-1 * log($accuracyAA) / $denom) + $baseAA;
			$baseRadiusAA = ($e/log10($maxAA)) * $baseBlurAA;
			$accuracyAdjustAA = (-1 * log10($accuracyAA) / log10(100 * $baseBlurAA)) * $baseBlurAA;
			$maxAAFrac = $maxAA / ($maxNucl + $maxAA);
		}
		my $minUnit = 50;
		my $bitRadiusNucl = ($baseRadiusNucl + $accuracyAdjustNucl) * $maxNucl;
		my $bitRadiusAA = ($baseRadiusAA + $accuracyAdjustAA) * $maxAA;
		my $radNucl = ($baseRadiusNucl + $accuracyAdjustNucl);
		my $radAA = ($baseRadiusAA + $accuracyAdjustAA);
		my $bitRadius = $bitRadiusNucl + $bitRadiusAA;
		if ($bitRadius < $minUnit) {
			$bitRadius = $minUnit;
			$bitRadiusNucl = $maxNuclFrac * $minUnit;
			$bitRadiusAA = $maxAAFrac * $minUnit;
		}
		my $lowerLimit = $max - $bitRadius;
		my $lowerLimitNucl = $maxNucl - $bitRadiusNucl;
		my $lowerLimitAA = $maxAA - $bitRadiusAA;
		my $printLimit = int($lowerLimit);
		my $printLimitNucl = int($lowerLimitNucl);
		my $printLimitAA = int($lowerLimitAA);
		my $worryLimit = $max - (.14 * $max);
		my $printWorry = int($worryLimit);
		lockPrintStderr("OrigTaxId:$refTaxId");
		my $bestAA = {};
		my $bestNucl = {};
		my $pool = {};
		my $bestOverall = {};
		my $poolScores = {};
		my $poolRep = {};	
		sub scaleForAccuracy {
			my $rep = $_[0];
			my $accuracy = $_[1];
			my $base = $_[2];
			my $minRep = $_[3];
			if ($accuracy == 1 or $rep == 0) {
				return($rep);
			}
			return(((1/$minRep) * $rep) ** ($accuracy ** $base));
		}
		if (not $freePass) {
			lock($lock);
			lockPrintStderr("bitRadius = $bitRadius");
			lockPrintStderr("max possible Nucleotide Score = $selfMaxNucl, max Poss AA = $selfMaxAA");
			if ($maxNucl) {
				lockPrintStderr("lower limit of scoring nucl: $printLimitNucl max nucl = $maxNucl set by $maxNuclTaxa - $names->{$maxNuclTaxa}");
			}
			if ($maxAA) {
				lockPrintStderr("lower limit of scoring prot: $printLimitAA max AA = $maxAA set by $maxAATaxa - $names->{$maxAATaxa}");
			}
			lockPrintStderr("baseRadiusNucl: $baseRadiusNucl, accuracyAdjustNucl: $accuracyAdjustNucl, bitRadiusNucl: $bitRadiusNucl");
			lockPrintStderr("baseRadiusAA: $baseRadiusAA, accuracyAdjustAA: $accuracyAdjustAA, bitRadiusAA: $bitRadiusAA");
			lockPrintStderr("lower limit of scoring overall if needed: $printLimit");
			lockPrintStderr("worry limit of scoring: $printWorry");
			lockPrintStderr("");
			lockPrintStderr("Nucl Accuracy: $accuracyNucl, nuclRadius: $radNucl, baseBlurNucl = $baseBlurNucl, approximate Nucl identity: $approxIdNucl");
			lockPrintStderr("Prot Accuracy $accuracyAA, protRadius: $radAA, baseBlurAA = $baseBlurAA");
		}
		foreach my $group (keys %$groupScores) {
			if (exists $suppressors->{$group}) {
				next;
			}
			foreach my $taxa (@{$groupScoresSortedOrder->{$group}}) {
				if ($taxa == $maxNuclTaxa or $taxa == $maxAATaxa) {
					$pool->{$group} = 1;
					$poolScores->{$taxa} = $overallScores->{$taxa};
				} elsif  ($overallScores->{$taxa} > $lowerLimit and exists $overallScoresNucl->{$taxa} and $overallScoresNucl->{$taxa} > $lowerLimitNucl) {
					if (exists $overallScoresAA->{$taxa} and $overallScoresAA->{$taxa} < $lowerLimitAA and $overallScoresAA->{$taxa} > 200) {
						next;
					}
					$pool->{$group} = 1;
					$poolScores->{$taxa} = $overallScores->{$taxa};
				}
			}
			unless (exists $pool->{$group}) {
				foreach my $taxa (@{$groupScoresSortedOrder->{$group}}) {
					if ($overallScores->{$taxa} > $lowerLimit) {
						if (exists $overallScoresNucl->{$taxa} and $overallScoresNucl->{$taxa} < $lowerLimitNucl) {
							next;
						}
						$pool->{$group} = 1;
						$poolScores->{$taxa} = $overallScores->{$taxa};
					}
				}
			}
			unless (exists $pool->{$group}) {
				foreach my $taxa (@{$groupScoresSortedOrder->{$group}}) {
					if ($overallScores->{$taxa} > $lowerLimit) {
						$pool->{$group} = 1;
						$poolScores->{$taxa} = $overallScores->{$taxa};
					}
				}
			}
		}
		my $minRepAA = 1;
		my $minRepNucl = 1;
		my $totalRawRep = 0;
		foreach my $taxa (keys %$poolScores) {
			if ($taxaRepAA->{$taxa} < $minRepAA and $taxaRepAA->{$taxa} > 0) {
				$minRepAA = $taxaRepAA->{$taxa};
			}
			if ($taxaRepNucl->{$taxa} < $minRepNucl and $taxaRepNucl->{$taxa} > 0) {
				$minRepNucl = $taxaRepNucl->{$taxa};
			}
			$totalRawRep += $taxaRepAA->{$taxa} + $taxaRepNucl->{$taxa};
		}
		foreach my $taxa (keys %$poolScores) {
			$poolRep->{$taxa} = scaleForAccuracy($taxaRepAA->{$taxa}, $accuracyAA, $taxaAccuracyBase, $minRepAA) + scaleForAccuracy($taxaRepNucl->{$taxa}, $approxIdNucl, $taxaAccuracyBase, $minRepNucl);
		}
		if (not $freePass) {
			lock($lock);
			lockPrintStderr("Scaled DB Rep");
			lock($lock);
			select->flush()
		}
		if (not $freePass and not exists $pool->{'virus'}) {
				lockPrintStderr("FILTERED: $reference has no virus in the top 10% of scoring taxa. Max score = $max, virus max = $maxes->{'virus'}");
				next;
		}
		foreach my $pos (keys %$groupScorePerPosAA) {
			foreach my $group (keys %{$groupScorePerPosAA->{$pos}}) {
				$bestAA->{$group} += $groupScorePerPosAA->{$pos}->{$group}
			}
		}
		foreach my $pos (keys %$groupScorePerPosNucl) {
			foreach my $group (keys %{$groupScorePerPosNucl->{$pos}}) {
				$bestNucl->{$group} += $groupScorePerPosNucl->{$pos}->{$group};
			}
		}

		foreach my $group (keys %$groupScores) {
			unless (exists $bestAA->{$group}) {
				$bestAA->{$group} = 0;
			}
			unless (exists $bestNucl->{$group}) {
				$bestNucl->{$group} = 0;
			}
			$bestAA->{$group} = int($bestAA->{$group});
			$bestNucl->{$group} = int($bestNucl->{$group});
			$bestOverall->{$group} = $bestNucl->{$group} + $bestAA->{$group};
			if ($bestOverall->{$group} > $selfMax) {
				$selfMax = $bestOverall->{$group};
			}
		}
		my $newTaxaId = "";
		my $misannotations = {};
		my @poolMembers = keys %$poolRep;
		if (scalar(@poolMembers) and not $freePass) {
			if (scalar(@poolMembers) > 1) {		
				($newTaxaId, $misannotations) = findBestLCA($parents, $ranks, $poolRep, $poolScores);
				lockPrintStderr("AFTER LCA = $newTaxaId");
			} else {
				$newTaxaId = $poolMembers[0];
			}
		} else {
			$newTaxaId = $refTaxId;
		}
		if (not $freePass and not isDescendant($parents, $newTaxaId, $groups->{'virus'})) {
			lockPrintStderr("FILTERED: $reference resolves to $newTaxaId - $names->{$newTaxaId}");
			next;
		}
		
		my @misannotations = keys %$misannotations;
		my $misannoStr = "";
		if (scalar(@misannotations)) {
			$misannoStr = join ",", @misannotations;
		}
		if ($newTaxaId != $refTaxId) {
			lockPrintStderr("Changed Taxonomy of $reference from $refTaxId to $newTaxaId - $names->{$newTaxaId}");
			$reference =~ s/taxId=\d+/taxId=$newTaxaId/;
			$reference = "$reference" . ";origTaxId=$refTaxId";
		}
		$selfMax = int($selfMax);
		my $totalMax = int($selfMax + $nonFilterAddNucl + $nonFilterAddAA);
		unless (exists $bestOverall->{'virus'}) {
			$bestOverall->{'virus'} = 0;
		}
		if ($totalMax < $bestOverall->{'virus'}) {
			$totalMax = $bestOverall->{'virus'};
		}
		unless (exists $bestNucl->{'virus'}) {
			$bestNucl->{'virus'} = 0;
		}
		unless (exists $bestAA->{'virus'}) {
			$bestAA->{'virus'} = 0;
		}
		$selfMaxAA = int($selfMaxAA);
		$selfMaxNucl = int($selfMaxNucl);
		$reference = "$reference;maxScore=$bestOverall->{'virus'};maxAlign=$selfMax;maxPossibleScore=$totalMax;protScore=$bestAA->{'virus'};maxProtAlign=$selfMaxAA;nuclScore=$bestNucl->{'virus'};maxNuclAlign=$selfMaxNucl";
		my $nUsage = () = $scaffold =~ /[Nn]/gi;
		my $totalLength = length($scaffold);
		my $letterOccupancy = (($totalLength - $nUsage) / $totalLength) * 100;
		$letterOccupancy = sprintf("%.2f", $letterOccupancy);
		$reference = "$reference;letterOccupancy=$letterOccupancy";
		if (exists $bestOverall->{'virus'} and $bestOverall->{'virus'} < 1000) {
			lockPrintStderr("WEAK:$reference");
			$reference .= ";weak=1";
		}
		if (not $highUnknown and exists $bestOverall->{'virus'} and $bestOverall->{'virus'} and $bestOverall->{'virus'}/$totalMax < .8) {
			lockPrintStderr("HIGH DIV:$reference");
			$reference .= ";highDivergence=1";
		}
		if ($highUnknown) {
			$reference .= ";highUnknownInformation=1";
		} else {
			foreach my $group (keys %$bestOverall) {
				if ($group eq 'virus' or $bestOverall->{$group} < $worryLimit) {
					next;
				} else {
					my $groupPrint = ucfirst($group);
					$reference .= ";high$groupPrint=$bestOverall->{$group}";
				}
			}
		}
		if ($misannoStr) {
			$reference .= ";potentialMisannotations=$misannoStr";
			lockPrintStderr("POTENTIL MISANNOTATIONS: $misannoStr");
		}
		if ($whiteListed) {
			$reference .= ";whiteListed=1";
			lockPrintStderr("WHITELISTED: $reference");
		}
		lockPrintStderr("END $reference on thread $tid");
		$return->{$reference} = $scaffold;

	}
	lockPrintStderr("$tid finished processing");
	return($return);
}

my $mapped = {};
foreach my $head (keys %$scaffolds) {
	$mapped->{$head} = 0;
}
my $aaHit = {};
my $ref2Taxa = {};
my $singleMode = 0;
if (-e $ARGV[2]) {
	my $addCount = packBlastDbIntoRocksDb($aaHit, $scaffolds, $dbAA, $mapped, $ARGV[2], $ref2Taxa, $THREADS, $verb);
	print STDERR "finished parsing aa\n";
	unless ($addCount) {
		$singleMode = 1;
	}
} else {
	$singleMode = 1;
}
my $nucHit = {};
if (-e $ARGV[3]) {
	my $addCount = packBlastDbIntoRocksDb($nucHit, $scaffolds, $dbNUC, $mapped, $ARGV[3], $ref2Taxa, $THREADS, $verb);
	print STDERR "Finished parsing nuc\n";
	unless ($addCount) {
		$singleMode = 1;
	}
} else {
	$singleMode = 1;
}
if ($singleMode) {
	print STDERR "working in single database mode\n";
}
$dbTAX->put("singleMode", "$singleMode");
foreach my $head (keys %$scaffolds) {
	my $taxId;
	$head =~ m/taxId=(\d+)/;
	$taxId = $1;
	unless ($taxId) {
		die "can't get taxId on reference $head\n";
	}
	$ref2Taxa->{$head}->{$taxId} = 1;
	unless (exists $aaHit->{$head}) {
		my $temp = [];
		my $compressed = compress($encoder->encode($temp));
		$dbAA->put("$head" => $compressed);
	}
	unless (exists $nucHit->{$head}) {
		my $temp = [];
		my $compressed = compress($encoder->encode($temp));
		$dbNUC->put("$head" => $compressed);
	}

	
}

print STDERR "prebuild taxa structures\n";
my $corrections = {};
my $knownGood = {};
my $taxa2Group = {};
foreach my $reference (keys %$ref2Taxa) {
	my $newParents = {};
	my $newRanks = {};
	my $newNames = {};
	my $newCorrections = {};
	my $newTaxa2Group = {};
	foreach my $taxa (keys %{$ref2Taxa->{$reference}}) {
		unless (defined $taxa) {
			die;
		}
		my $good = 1;
		my $check = 0;
		foreach my $collapseTaxa (keys %$lookUpTree) {
			if (exists $corrections->{$taxa}) {
				$newCorrections->{$taxa} = $corrections->{$taxa};
				$taxa = $corrections->{$taxa};
				last;
			} elsif (exists $knownGood->{$taxa}) {
				last;
			} elsif (exists $children->{$taxa}) {
				last;
			}
			$check = 1;
			foreach my $layers (@{$lookUpTree->{$collapseTaxa}}) {
				my $lookUp = checkUpToSpecies($parents, $ranks, $names, $taxa, $collapseTaxa, $layers);
				if ($lookUp and $lookUp != $taxa) {
					$good = 0;
					$corrections->{$taxa} = $lookUp;
					$newCorrections->{$taxa} = $lookUp;
					$taxa = $lookUp;
					last;
				}
			}
		}
		if ($check and $good) {
			$knownGood->{$taxa} = 1;
		}
		foreach my $group (keys %$groups) {
			if (exists $taxa2Group->{$taxa}) {
				$newTaxa2Group->{$taxa} = $taxa2Group->{$taxa};
				last;
			}
			if (isDescendant($parents, $taxa, $groups->{$group})) {
				$newTaxa2Group->{$taxa} = $group;
				$taxa2Group->{$taxa} = $group;
				last;
			}
		}
		unless (exists $newTaxa2Group->{$taxa}) {
			$newTaxa2Group->{$taxa} = "unknown";
			$taxa2Group->{$taxa} = "unknown";
		}
		$newParents->{$taxa} = $parents->{$taxa};
		$newRanks->{$taxa} = $ranks->{$taxa};
		$newNames->{$taxa} = $names->{$taxa};
		my $escape = 0;
		while (1) {
			$escape++;
			if ($escape == 100) {
				print STDERR "escaped trace\n";
				last;
			}
			unless (exists $parents->{$taxa}) {
				last;
			}
			$taxa = $parents->{$taxa};
			if ($taxa == 1) {
				last;
			}
			if (exists $newParents->{$taxa}) {
				last;
			}
			$newParents->{$taxa} = $parents->{$taxa};
			$newRanks->{$taxa} = $ranks->{$taxa};
			$newNames->{$taxa} = $names->{$taxa};
			
		
		}
	}
	$newParents->{"1"} = $parents->{"1"};
	$newRanks->{"1"} = $ranks->{"1"};
	$newNames->{"1"} = $names->{"1"};
	my $compressedParents = compress($encoder->encode($newParents));
	my $compressedNames = compress($encoder->encode($newNames));
	my $compressedRanks = compress($encoder->encode($newRanks));
	my $compressedTaxa2Group = compress($encoder->encode($newTaxa2Group));
	my $compressecCorrections = compress($encoder->encode($newCorrections));
	
	$dbTAX->put("par.$reference", $compressedParents);
	$dbTAX->put("nam.$reference", $compressedNames);
	$dbTAX->put("ran.$reference", $compressedRanks);
	$dbTAX->put("t2g.$reference", $compressedTaxa2Group);
	$dbTAX->put("cor.$reference", $compressecCorrections);
		
}
print STDERR "done prebuild\n";
$dbAA->compact_range;
$dbNUC->compact_range;
$dbTAX->compact_range;
print STDERR "done compact\n";
undef($dbAA);
undef($dbNUC);
undef($dbTAX);
undef($parents);
undef($ranks);
undef($names);
undef($ref2Taxa);
my $goCount = 0;
print STDERR "copying structures\n";
print STDERR "done copying\n";

my $scaffoldCount = scalar(keys %$scaffolds);






print STDERR "before enqueue scaffolds\n";
foreach my $head (sort {$mapped->{$b} <=> $mapped->{$a}} keys %$scaffolds) {
	$q->enqueue($head);
}
print STDERR "after enqueue scaffolds\n";

for (my $i = 0; $i < $THREADS; $i++) {
	$startQ->enqueue("GO");
}

print STDERR "after enqueue go\n";

undef($mapped);

sleep(10);

for (my $i = 0; $i < $THREADS; $i++) {
	$q->enqueue(undef);
}
print STDERR "After enqueue undef\n";

my $seqs = {};
my @running = threads->list(threads::running);
while (1) {
	if (scalar(@running)) {
		sleep(1);
		@running = threads->list(threads::running);
		foreach my $thr (@running) {
			my $tid = $thr->tid();
		}
	} else {
		last;
	}
}
foreach my $thr (@threads) {
	my $tid = $thr->tid();
	my $temp = $thr->join();
	foreach my $ref (keys %{$temp}) {
		$seqs->{$ref} = $temp->{$ref};
	}
}

foreach my $ref (sort keys %$seqs) {
	print ">$ref\n";
	print "$seqs->{$ref}\n";
}

sub findBestLCA {
	my $parents = $_[0];
	my $ranks = $_[1];
	my $rep = $_[2];
	my $scores = $_[3];
	my $tid = threads->tid();
	my $misannotations = {};
	my $scoreAmount = {
		"superkingdom" => 1,
		"phylum" => 2,
		"class" => 3,
		"order" => 6,
		"family" => 13,
		"genus" => 25,
		"species" => 50
	};
	my $wordsBest = {};
	my $best = {};
	my $e = 2.7182818284590452353602874713527;
	my $lineagesPerRank = [];
	my $scoresPerNode = {};
	my $totalRep = 0;
	foreach my $taxaId (keys %$rep) {
		my $node = $taxaId;
		$totalRep += $rep->{$taxaId};
		my $escape = 0;
		my @localRoute;
		while (1) {
			$escape++;
			if ($escape == 100) {
				last;
			}
			push @localRoute, $node;
			if (not exists $scoresPerNode->{$node} or $scoresPerNode->{$node} < $scores->{$taxaId}) {
				$scoresPerNode->{$node} = $scores->{$taxaId};
			}
			unless (exists $parents->{$node}) {
				last;
			}
			$node = $parents->{$node};
			if ($node == 1) {
				last;
			}	
		}
		@localRoute = reverse(@localRoute);
		for (my $i = 0; $i < scalar(@localRoute); $i++) {
			my $parent = 1;
			if ($i > 0) {
				my $parentIdx = $i - 1;
				$parent = $localRoute[$parentIdx]
			}
			$lineagesPerRank->[$i]->{$localRoute[$i]} = $parent;
		}
	}
	my @reverseTree = reverse(@{$lineagesPerRank});
	my @bestScoresPerNode = sort {$scoresPerNode->{$b} <=> $scoresPerNode->{$a}} keys %$scoresPerNode;
	my $maxNodeScore = $scoresPerNode->{$bestScoresPerNode[0]};
	my $limitNodeScore = .99 * $maxNodeScore;
	lockBrowseErr("Lineage:", $lineagesPerRank);
	lockPrintStderr("TotalRep: $totalRep");
	my $outMessage = "noMajority";
	my $denom = $totalRep;
	if ($denom < $e) {
		$denom = $e;
	}
	my $superMajBase = (2/3);
	my $adjustment = ($superMajBase - .5) / log2($denom);
	my $superMajRepLimit = ($superMajBase * $totalRep);
	#### Majority Threshold calulated here ####
	my $majRepLimit = (.5 + $adjustment) * $totalRep;
	lockPrintStderr("Supermajority Limit = $superMajRepLimit");
	lockPrintStderr("Slight majority Limit = $majRepLimit");
	lockPrintStderr("Adjustment to make majority limit = $adjustment");
	my $tieRanks = {};	
	for (my $i = 0; $i < scalar(@reverseTree); $i++) {
		foreach my $taxaId (sort {$rep->{$b} <=> $rep->{$a}} keys %{$reverseTree[$i]}) {
			$tieRanks->{$taxaId} = $i;
			if ($taxaId == 10239) {
				last;
			}
			if ($rep->{$taxaId} > $majRepLimit) {
				lockBrowseErr("DbScores at super majority:", $rep);
				if ($scoresPerNode->{$taxaId} < $limitNodeScore) {
					lockPrintStderr("WARNING: Score for node $taxaId is less than 99% of the top possible score: $maxNodeScore");
				}
				if ($rep->{$taxaId} < $superMajRepLimit) {
					$misannotations->{'noSupermajority'} = 1;
				}
				return($taxaId, $misannotations);
			} 
			if (exists $parents->{$taxaId}) {
				$rep->{$parents->{$taxaId}} += $rep->{$taxaId};
			}
		}
	}
	my @sortedPoolScore = sort {$rep->{$b} <=> $rep->{$a} || $tieRanks->{$a} <=> $tieRanks->{$b} || $a <=> $b} keys %$rep;
	shift(@sortedPoolScore);
	lockBrowseErr("DbScores end:", $rep);
	lockPrintStderr("WARNING: Score for node $sortedPoolScore[0] did not achieve required representation");
	if ($sortedPoolScore[0] > .55 * $totalRep and $sortedPoolScore[0] < $superMajRepLimit) {
		$outMessage = "noSupermajority";
	} elsif ($sortedPoolScore[0] < .55 * $totalRep and $sortedPoolScore[0] > .5 * $totalRep) {
		$outMessage = "noSlightMajority";
	} else {
		$outMessage = "noMajority";
	}
	$misannotations->{$outMessage} = 1;
	return($sortedPoolScore[0], $misannotations);
}



sub isDescendant {
	my $parents = $_[0];
	my $node = $_[1];
	my $check = $_[2];
	my $escape = 0;
	while (1) {
		$escape++;
		if ($escape == 100) {
			warn "escaped isVirus\n";
			return 0;
		}
		if ($node == $check) {
			return 1;
		}
		if ($node == 1 or $node == 0) {
			return 0;
		}
		unless (exists $parents->{$node}) {
			return 0;
		}
		$node = $parents->{$node};
	}

}

sub entropy {
	my $calcStr = $_[0]; 
	$calcStr =~ s/[Nn]//g;
	$calcStr = uc($calcStr);
	my $k = 3;
	my $end = length($calcStr);
	my $total = 0;
	my $kmers = {};
	my @kmers;
	for (my $i = 0; $i < $end - $k; $i++) {
		$total++;
		my $kmer = substr($calcStr, $i, $k);
		$kmers->{$kmer}++;
		push @kmers, $kmer;
	}
	my $entropy = 0;
	foreach my $kmer (@kmers) {
		my $probability = $kmers->{$kmer} / $total;
		$entropy += log2($probability);
	}
	$entropy *= -1;
	$entropy = $entropy/$total;
	return($entropy);

}


sub checkUpToSpecies {
	my $parents = $_[0];
	my $ranks = $_[1];
	my $names = $_[2];
	my $refTaxa = $_[3];
	my $checkTaxa = $_[4];
	my $layersUp = $_[5];
	if ($refTaxa == 0) {
		return(0);
	}
	my $escape = 0;
	my $node = $refTaxa;
	my $hasIt = 0;
	my $killRanks = {"superkingdom" => 1, "kingdom" => 1, "phylum" => 1, "class" => 1, "order" => 1, "family" => 1, "genus" => 1, "species" => 1};
	while (1) {
		$escape++;
		if ($escape == 100) {
			last;
		}
		my $localClimb = $node;
		for (my $i = 0; $i < $layersUp; $i++) {
			unless (exists $parents->{$localClimb}) {
				last;
			}
			$localClimb = $parents->{$localClimb};
		}
		unless (exists $ranks->{$node}) {
			return(0);
		}
		if ($localClimb == $checkTaxa) {
			my $refName = $names->{$refTaxa};
			my $nodeName = $names->{$node};
			$hasIt = $node;
			last;
		} elsif (exists $killRanks->{$ranks->{$node}} or $parents->{$node} == 1) {
			last;
		} else {
			$node = $parents->{$node};
		}
	}
	return ($hasIt);
}

sub packBlastDbIntoRocksDb {

	my $hits = $_[0];
	my $scaffolds = $_[1];
	my $db = $_[2];
	my $mapped = $_[3];
	my $file = $_[4];
	my $ref2Taxa = $_[5];
	my $THREADS = $_[6];
	my $mode = $_[7];
	my $skipCount = 100000;;
	if ($mode eq "filter") {
		$skipCount = 1000;;
	}
	my $sortStr = "--stable -k1,1 -k12,12nr -k3,3nr --buffer-size=25G --parallel=$sortThreads --compress-program=lz4";
	my $lineCount = 0;
	open IN, "lbzip2 -d -c -n$THREADS $file | sort $sortStr |";
	my $lines = [];
	my $currentHead;
	my $line;
	my $nextLine;
	my $seenLine = 0;
	my $addCount = 0;
	my $countCont = {};
	while (1) {
		my $eof = 0;
		my $head = "";
		$nextLine = <IN>;
		unless (defined($nextLine)) {
			$eof = 1;
		}
		if ($eof and not $seenLine) {
			last;
		}
		$seenLine = 1;
		if ($line) {
			chomp $line;
			my @parts = split /\t/, $line;
			$head = $parts[0];
			$countCont->{$head}++;
			if ($countCont->{$head} > $skipCount and not $eof) {
				$line = $nextLine;
				next;
			}
			$head =~ s/;coords=.*//g;
			my $target = $parts[1];
			if ($target =~ m/;taxId=(\d+)/) {
				$ref2Taxa->{$head}->{"$1"} = 1;
			} else {
				$line = $nextLine;
				next;
			}
			$mapped->{$head} += $parts[3];
		} else {
			my @parts = split /\t/, $nextLine;
			$currentHead = $parts[0];
			$currentHead =~ s/;coords=.*//g;
		}
		if (not exists $scaffolds->{$head} and not $eof) {
			$line = $nextLine;
			next;
		}
		if ($currentHead ne $head or $eof) {
			if ($eof and $currentHead ne $head) {
				packLines($db, $currentHead, $lines);
				$addCount++;
				$hits->{$currentHead} = 1;
				$lines = [];
				$currentHead = $head;
			} 
			if ($currentHead eq $head) {
				push @{$lines}, $line;
			}
			$hits->{$currentHead} = 1;
			packLines($db, $currentHead, $lines);
			$addCount++;
			if ($eof) {
				last;
			}
			$lines = [];
			$currentHead = $head;
		}
		if ($line) {
			push @{$lines}, $line;
		}
		$line = $nextLine;
	}
	close IN;
	return($addCount);
}
sub packLines {
	my $db = $_[0];
	my $head = $_[1];
	my $lines = $_[2];
	my $compressed = compress($encoder->encode($lines));
	$db->put("$head" => $compressed);

}

sub lockPrintStderr {
	lock($lock); 
	my $tid = threads->tid();
	my $printout = $_[0];
	print STDERR "THREAD $tid - $printout\n";
}
sub lockBrowseErr {
	lock($lock);
	my $message = $_[0];
	my $toBrowse = $_[1];
	my $tid = threads->tid();
	lockPrintStderr("$message");

}
		

sub makeTaxonomy {
	my $in = $_[0];
	open NODE, "$in";
	my $serealData;
	binmode(NODE);
	{
		local $/;
		undef $/;
		
		$serealData = <NODE>;
	}
	close NODE;
	my $data = $decoder->decode(decompress($serealData));

	my $parents = $data->{'parents'};
	my $ranks = $data->{'ranks'};
	my $names = $data->{'names'};
	my $children = $data->{'children'};
	return($parents, $names, $ranks, $children);

}

sub log10 {
	return(log($_[0]) / log(10));
}

sub log2 {
	return(log($_[0]) / log(2))
}

sub bit2raw {
	my $bitScore = $_[0];
	my $lambda = $_[1];
	my $k = $_[2];
	return((($bitScore * log(2)) + log($k)) / $lambda);

}

sub raw2bit {
	my $rawScore = $_[0];
	my $lambda = $_[1];
	my $k = $_[2];
	return((($lambda * $rawScore) - log($k)) / log(2))


}











