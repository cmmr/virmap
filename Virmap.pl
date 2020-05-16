#!/usr/bin/env perl

use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;
use FAlite;
use File::Temp qw/tempdir/;
use Time::HiRes qw(time);
use Fcntl; 
use OpenSourceOrg::API;


my $client = OpenSourceOrg::API->new();

unless($client) {
	die;
}

my $startTime = time();
my $cpuStart = 0;
my @files;
my @readF;
my @readR;
my @interleaved;
my $gbBlastx;
my $gbBlastn;
my $virBbmap;
my $virDmnd;
my $taxaJson;
my $outputDir;
my $fastaInput;
my $sampleName;
my $bigRam;
my $hugeRam;
my $krakenFilter;
my $krakenDb;
my $loose;
my $bbMapLimit = 4;
my $threads = 4;
my $tmpdir;
my $sensitive;
my $whiteList = "";
my $useMegaHit = 0;
my $improveTimeLimit = 36000;
my $commandLineArgs = join " ", @ARGV;
my $noNucMap;
my $noAaMap;
my $noIterative;
my $skipTaxonomy;
my $noNormalize;
my $bothTadpoleAndMegahit = 0;
my $useTadpole = 1;
my $keepTemp = 0;
my $noCorrection;
my $noFilter;
my $noAssm;
my $strict;
my $noMerge;
my $noEntropy;
my $useBbnorm;
my $infoFloor = 300;


my $startMsg = "Virmap called with: $commandLineArgs\n";

unless (scalar(@ARGV)) {
	my $help = <<'HELPSTRING';
usage: Virmap.pl [options] <databases> --readF ReadSet1_R1.fastq.bz2 ReadSet2_R1.fastq.gz --readR ReadSet1_R2.fastq.bz2 ReadSet2_R2.fastq.bz2 --readUnpaired ReadsUnpaired.fastq.bz2

HELPSTRING

print "$help\n";
exit();
}
print STDERR "$startMsg\n";


GetOptions ("interleaved:s{,}" => \@interleaved, "readF:s{,}" => \@readF, "readR:s{,}" => \@readR, "readUnpaired:s{,}" => \@files, "gbBlastx=s" => \$gbBlastx, "gbBlastn=s" => \$gbBlastn, "virBbmap=s" => \$virBbmap, "virDmnd=s" => \$virDmnd, "taxaJson=s" => \$taxaJson, "outputDir=s" => \$outputDir, "fasta" => \$fastaInput, "sampleName=s" => \$sampleName, "threads=i" => \$threads, "tmpdir=s" => \$tmpdir, "bigRam" => \$bigRam, "hugeRam" => \$hugeRam, "whiteList=s" => \$whiteList, "bbMapLimit=s" => \$bbMapLimit, "loose" => \$loose, "improveTimeLimit=s" => \$improveTimeLimit, "useMegahit" => \$useMegaHit, "sensitive" => \$sensitive, "noNucMap" => \$noNucMap, "noAaMap" => \$noAaMap, "noIterImp" => \$noIterative, "skipTaxonomy" => \$skipTaxonomy, "both" => \$bothTadpoleAndMegahit, "noCorrection" => \$noCorrection, "noFilter" => \$noFilter, "noNormalize" => \$noNormalize, "noAssembly" => \$noAssm, "strict" => \$strict, "noMerge" => \$noMerge, "noEntropy" => \$noEntropy, "infoFloor" => \$infoFloor, "keepTemp" => \$keepTemp, "useBbnorm" => \$useBbnorm, "krakenFilter" => \$krakenFilter, "krakenDb=s" => \$krakenDb);

if ($useMegaHit) {
	$useTadpole = 0;
}

if ($useMegaHit and not $bothTadpoleAndMegahit) {
	$useTadpole = 0;
}
if ($useTadpole and not $bothTadpoleAndMegahit) {
	$useMegaHit = 0;
}

if ($bothTadpoleAndMegahit) {
	$useTadpole = 1;
	$useMegaHit = 1;
}
if ($noAssm) {
	$useTadpole = 0;
	$useMegaHit = 0;
}
if ($noAaMap and $noNucMap) {
	$noNormalize = 1;
}

unless ($infoFloor =~ m/^\d+$/) {
	$infoFloor = 300;
}
my $strictMerge = "";
if ($sensitive) {
	$strictMerge = "strict";
}
if ($noFilter and not $noIterative) {
	$noFilter = 0;
}

#hard coded to >=4 threads, probably should fix/be fine.
if ($threads < 4 or $threads =~ m/[^\d]/) {
	$threads = 4;
}
my $fail = 0;
my $noTaxaDep = 0;
if ($skipTaxonomy and $noFilter and $noIterative) {
	$noTaxaDep = 1;
}
if (scalar(@readF) != scalar(@readR)) {
	print STDERR "Amount of forward reads must match reverse reads\n";
	$fail = 1;
}
push @files, @readF;
push @files, @readR;
push @files, @interleaved;
unless ($files[0]) {
	print STDERR "No input files in --files\n";
	$fail = 1;
}
unless ($gbBlastx or $noTaxaDep) {
	print STDERR "Amino acid Genbank database not defined in --gbBlastx\n";
	$fail = 1;
}
unless ($gbBlastn or $noTaxaDep) {
	print STDERR "Nucleotide Genbank database not defined in --gbBlastn\n";
	$fail = 1;
}
unless ($virBbmap or $noNucMap) {
	print STDERR "Virus Bowtie2 database not defined in --virBbmap\n";
	$fail = 1;
}
unless ($virDmnd or $noAaMap) {
	print STDERR "Virus Diamond database not defined in --virDmnd\n";
	$fail = 1;
}
unless ($taxaJson or $noTaxaDep) {
	print STDERR "Taxonomy JSON not defined in --taxaJson\n";
	$fail = 1;
}
unless ($outputDir) {
	print STDERR "Output directory not defined in --outputDir\n";
	$fail = 1;
}
unless ($sampleName) {
	print STDERR "Sample name not defined in --sampleName\n";
	$fail = 1;
}
if ($strict and $loose) {
	print STDERR "Cannot enable both strict and loose\n";
	$fail = 1;
}
if ($noAaMap and $noNucMap and $noAssm) {
	print STDERR "All mapping and assembly options are disabled; nothing to do\n";
	$fail = 1;
}
my $cleanUp = 1;
if ($keepTemp) {
	$cleanUp = 0;
}
unless ($tmpdir) {
	unless ($ENV{'TMPDIR'}) {
		print STDERR "TMPDIR not set\n";
		$fail = 1;
	}
	$tmpdir = tempdir( CLEANUP => $cleanUp );
}
if ($fail) {
	exit()
}


#setup

unless ($improveTimeLimit =~ m/^\d+$/) {
	$improveTimeLimit = "36000";	
}
my $ram = "50g";
my $RAM = "50G";
my $cdHitRam = "48000";
my $megaHitRam = 53687091200;
my $diamondParams = "";
if ($bigRam) {
	print STDERR "Enabling big RAM\n";
	$ram = "100g";
	$RAM = "100G";
	$cdHitRam = 100000;
	$megaHitRam = 107374182400;
	$diamondParams = "-b 6"
}
if ($hugeRam) {
	print STDERR "Enabling huge RAM\n";
	$ram = "1200g";
	$RAM = "120G";
	$cdHitRam = 1000000;
	$megaHitRam = 1073741824000;
	$diamondParams = "-b 100";
	
}
my $diamondPct = "80";
my $diamondRad = "8";
my $bbmapPct = ".90";
my $bbmapRad = ".95";
if ($loose) {
	$diamondPct = "60";
	$bbmapPct = ".80";
}
if ($strict) {
	$diamondPct = "85";
	$bbmapPct = ".93";
}
my $filePrep = "renameFastq2FastaVirome.pl";

if ($fastaInput) {
	$filePrep = "renameFasta2FastaVirome.pl";
}
if (-d $outputDir) {
	print STDERR "$outputDir already exists, output directory is not clean\n";
	exit();
} else {
	system("mkdir -p $outputDir");
}
my $prefix = "$outputDir/$sampleName";
my $tmpPrefix = "$tmpdir/$sampleName";
my $compressStdout = "lbzip2 -c -n$threads";
my $compress = "lbzip2 -n$threads";
if (-e "$prefix.fa" or -e "$prefix.fa.bz2") {
	print STDERR "FATAL: $prefix.fa and/or $prefix.fa.bz2 already exists, output directory is dirty\n";
	exit();
}
system("echo '$startMsg' >> $prefix.log");


#Compact all files into one fasta
my $seenFiles = {};
my $fileCount = 1;
foreach my $file (@files) {
	if (exists $seenFiles->{$file}) {
		print STDERR "$file is present more than once in --files\n";
		next;
	} else {
		$seenFiles->{$file} = 1;
	}
	system("$filePrep $file $fileCount >> $tmpPrefix.fa");
	$fileCount++;
}
system("$compressStdout $tmpPrefix.fa > $prefix.fa.bz2");
my (($timeDiff, $cpuDiff)) = timeDiff($startTime, "0", "decompress");

#dereplicate
system("vsearch --fasta_width 0 --derep_fulllength $tmpPrefix.fa --sizeout --strand both --threads $threads -output $tmpPrefix.derep.fa 1>$prefix.derep.err 2>&1");
system("$compressStdout $tmpPrefix.derep.fa > $prefix.derep.fa.bz2");
($timeDiff, $cpuDiff) = timeDiff($timeDiff, $cpuDiff, "dereplicate");

my $fileToMap = "$tmpPrefix.derep.fa"; 
my $filePseudoScaf = "$prefix.derep.fa.bz2";

#normalize
unless ($noNormalize) {
	unless ($useBbnorm) {
		system("normalize-by-median.py -M $RAM -C 5 --ksize 31 -o $tmpPrefix.normalized.fa $tmpPrefix.derep.fa 2>$prefix.normalize.err");
	} else {
		system("bbnorm.sh -Xmx$ram in=$tmpPrefix.derep.fa out=$tmpPrefix.normalized.fa bits=16 passes=1 tmpdir=$tmpdir fastawrap=0 k=31 target=5 threads=$threads percentile=50 min=0 minprob=0 2>$prefix.normalize.err");
	}
	system("$compressStdout $tmpPrefix.normalized.fa > $prefix.normalized.fa.bz2");
	$fileToMap = "$tmpPrefix.normalized.fa";
	$filePseudoScaf = "$prefix.normalized.fa.bz2";
	($timeDiff, $cpuDiff) = timeDiff($timeDiff, $cpuDiff, "normalize");
}

#align to virus bbmap
unless ($noNucMap) {
	system("bbmap.sh -Xmx$ram noheader=t threads=$threads in=$fileToMap path=$virBbmap noheader=t ambiguous=all ignorefrequentkmers=f slow=t excludefraction=0 greedy=f usejni=t maxsites2=10000000 minid=$bbmapPct outm=stdout.sam secondary=t sssr=$bbmapRad maxsites=100000000 sam=1.3 nmtag=t 32bit=t statsfile=$prefix.bbmap.err 2>$tmpPrefix.bbmap.err | zstd -q -c -T$threads > $tmpPrefix.nuc.sam.zst");
	system("zstd -q -dc $tmpPrefix.nuc.sam.zst | sort -k3,3 -k4,4n --buffer-size=$RAM --parallel=$threads --compress-program=lz4 | lbzip2 -c -n$threads > $prefix.nuc.sam.bz2");
	system("cat $tmpPrefix.bbmap.err >> $prefix.bbmap.err");
	($timeDiff, $cpuDiff) = timeDiff($timeDiff, $cpuDiff, "bbmap to virus");
} else {
	system("touch $prefix.nuc.sam; bzip2 $prefix.nuc.sam");
}

#align to virus diamond
unless ($noAaMap) {
	system("diamond blastx --algo 0 -p $threads -f 101 -q $fileToMap --db $virDmnd --unal 0 --masking 0 --top $diamondRad --id $diamondPct --query-cover 80 -c 1 -l 16 --evalue .1 --freq-sd 250 --index-mode 4 --shapes 1 --verbose -t $tmpdir 2>>$tmpPrefix.buildSuperScaffolds.err | grep -v \"^@\" | zstd -q -c -T$threads > $tmpPrefix.aa.sam.zst"); 
	system("zstd -q -dc $tmpPrefix.aa.sam.zst | sort -k3,3 -k4,4n --buffer-size=$RAM --parallel=$threads --compress-program=lz4 | lbzip2 -c -n$threads > $prefix.aa.sam.bz2");
	($timeDiff, $cpuDiff) = timeDiff($timeDiff, $cpuDiff, "diamond to virus");
} else {
	system("touch $prefix.aa.sam; bzip2 $prefix.aa.sam");
}

#cluster references, build pseudoscaffolds, scaffold into super scaffolds
my $superScaffolds = "";
unless ($noNucMap and $noAaMap) {
	system("clusterReferences.pl $prefix.aa.sam.bz2 $prefix.nuc.sam.bz2 1>$prefix.both.centroids.txt 2>>$tmpPrefix.buildSuperScaffolds.err");
	
	#make pseudo-scaffolds 
	system("buildPseudoScaffolds.pl $prefix.both.centroids.txt $filePseudoScaf $prefix.aa.sam.bz2 $prefix.nuc.sam.bz2 1>$prefix.pseudoScaffolds.fa 2>>$tmpPrefix.buildSuperScaffolds.err");

	#super-scaffold prot onto nuc
	system("superScaffold.pl $prefix.both.centroids.txt $prefix.pseudoScaffolds.fa > $prefix.superScaffolds.init.fa 2>>$tmpPrefix.buildSuperScaffolds.err");
	

	#self-improve superscaffolds 
	system("cat $prefix.superScaffolds.init.fa | toUC.pl > $tmpPrefix.superScaffolds.fa");
	system("bbmap.sh -Xmx$ram in=$tmpPrefix.derep.fa ref=$tmpPrefix.superScaffolds.fa nodisk=t usejni=t vslow=t minid=.8 local=t ignorefrequentkmers=f greedy=f excludefraction=0 ambiguous=random deterministic=t sam=1.3 32bit=t noheader=t outm=$tmpPrefix.superScaffolds.sam threads=$threads 2>>$tmpPrefix.buildSuperScaffolds.err");
	system("$compress $tmpPrefix.superScaffolds.sam");
	system("adjustSizesAndFilter.pl build 0 $tmpPrefix.superScaffolds.sam.bz2 $tmpPrefix.superScaffolds.fa /dev/null $threads > $prefix.superScaffolds.fa 2>>$tmpPrefix.buildSuperScaffolds.err");
	if (-s "$tmpPrefix.buildSuperScaffolds.err") {
		system("$compressStdout $tmpPrefix.buildSuperScaffolds.err > $prefix.superScaffolds.err.bz2");
	}
	$superScaffolds = "$prefix.superScaffolds.fa";
	($timeDiff, $cpuDiff) = timeDiff($timeDiff, $cpuDiff, "construct superscaffolds");
}


sub dumpPair {
	my @fileCont = @{$_[0]};
	my $outfile = $_[1];
	my $errorFile = $_[2];
	foreach my $file (@fileCont) {
		system("reformat.sh in=$file out=stdout.fq >> $outfile 2>>$errorFile");
	}
}

sub dumpInterleaved {
	my @fileCont = @{$_[0]};
	my $outfileF = $_[1];
	my $outfileR = $_[2];
	my $errorFile = $_[3];
	my $count = 0;
	my @forwards;
	my @reverses;
	foreach my $file (@fileCont) {
		system("reformat.sh in=$file verifyinterleaved=t out1=$outfileF.$count.tmp.1.fq out2=$outfileR.$count.tmp.2.fq 2>>$errorFile");
		push @forwards, "$outfileF.$count.tmp.1.fq";
		push @reverses, "$outfileR.$count.tmp.2.fq";
		$count++;
	}
	my $fwdTmp = join " ", @forwards;
	my $revTmp = join " ", @reverses;
	system("cat $fwdTmp >> $outfileF");
	system("rm $fwdTmp");
	system("cat $revTmp >> $outfileR");
	system("rm $revTmp");
}

#multi-k assemble
my $contigs = "";
if ($useMegaHit or $useTadpole) {
	my $paired = 0;
	if (scalar(@readF) and scalar(@readR)) {
		$paired = 1;
		dumpPair(\@readF, "$tmpPrefix.1.fq", "$tmpPrefix.assembly.err");
		dumpPair(\@readR, "$tmpPrefix.2.fq", "$tmpPrefix.assembly.err");
	}
	if (scalar(@interleaved)) {
		$paired = 1;
		dumpInterleaved(\@interleaved, "$tmpPrefix.1.fq", "$tmpPrefix.2.fq", "$tmpPrefix.assembly.err");
	}
	if ($useMegaHit) {
		my $megahitInputStr = "";
		if ($paired) {
			$megahitInputStr = "-1 $tmpPrefix.1.fq -2 $tmpPrefix.2.fq";
		} else {
			$megahitInputStr = "-r $tmpPrefix.fa";
		}
		unless ($sensitive) {
			system("megahit $megahitInputStr --min-contig-len 500 -t $threads -m $megaHitRam --tmp-dir $tmpdir -o $tmpPrefix.MegaHitAssm --out-prefix MegaHit 1>>$tmpPrefix.assembly.err 2>&1");
		} else {
			system("megahit $megahitInputStr --min-contig-len 500 -t $threads -m $megaHitRam --tmp-dir $tmpdir --presets meta-sensitive -o $tmpPrefix.MegaHitAssm --out-prefix MegaHit 1>>$tmpPrefix.assembly.err 2>&1");
		}
		system("cat $tmpPrefix.MegaHitAssm/MegaHit.contigs.fa >> $tmpPrefix.contigContainer.fa");
		($timeDiff, $cpuDiff) = timeDiff($timeDiff, $cpuDiff, "megahit assembly");
	}
	if ($useTadpole) {
		my $tadPoleInputStr = "";
		if ($paired) {
			$tadPoleInputStr = "in1=$tmpPrefix.1.fq in2=$tmpPrefix.2.fq";
		} else {
			$tadPoleInputStr = "in=$tmpPrefix.fa";
		}
		my $superParam = "";
		if ($superScaffolds) {
			$superParam = "extra=$superScaffolds";
		}
		if (not $noCorrection) {
			system("tadpole.sh -Xmx$ram usejni=t mode=correct ecc=t k=31 reassemble=t pincer=t tail=t eccfull=t conservative=t $tadPoleInputStr fastawrap=0 out=$tmpPrefix.tadpole.fq 1>>$tmpPrefix.assembly.err 2>&1");
			$tadPoleInputStr = "in=$tmpPrefix.tadpole.fq interleaved=t";
		}
		for (my $i = 29; $i <= 85; $i += 4) {
				system("tadpole.sh -Xmx$ram $tadPoleInputStr $superParam out=$tmpPrefix.assembly.$i.fa k=$i mincontig=500 minprobmain=f threads=$threads mincountseed=3 mincountextend=2 branchmult1=2 branchmult2=2 branchlower=2 shave=t rinse=t contigpasses=256 contigpassmult=1.05 minprob=0 usejni=t 1>>$tmpPrefix.assembly.err 2>&1");
			unless ($sensitive) {
				$i += 4;
			}
		}
		system("cat $tmpPrefix.assembly.*.fa >> $tmpPrefix.contigContainer.fa");
		($timeDiff, $cpuDiff) = timeDiff($timeDiff, $cpuDiff, "tadpole assembly");
	}
	if (-s "$tmpPrefix.contigContainer.fa") {
		system("renameContigs.pl $tmpPrefix.contigContainer.fa > $tmpPrefix.forDedupe.fa");
		system("dedupe.sh in=$tmpPrefix.forDedupe.fa out=$tmpPrefix.contigs.fa fastawrap=0 threads=1 ordered=t 2>>$tmpPrefix.assembly.err");
	} else {
		system("touch $prefix.contigs.fa");
	}
	if ($krakenFilter) {
		system("kraken2 --db $krakenDb --threads=$threads $tmpPrefix.contigs.fa 2>>$tmpPrefix.assembly.err | cut -f2,3 | lbzip2 -c -n$threads > $tmpPrefix.kraken.out.bz2");
		system("krakenFilter.pl $tmpPrefix.contigs.fa $tmpPrefix.kraken.out.bz2 $taxaJson > $prefix.contigs.fa 2>>$tmpPrefix.assembly.err");	
	
	} else {
		system("cp $tmpPrefix.contigs.fa $prefix.contigs.fa");
	}
	system("$compressStdout $tmpPrefix.assembly.err > $prefix.assembly.err.bz2");
	$contigs = "$prefix.contigs.fa";
	($timeDiff, $cpuDiff) = timeDiff($timeDiff, $cpuDiff, "dedupe assembly");
}

#combine de novo and pseudo assemblies or not
my $inFile = "";
if ($contigs and $superScaffolds) {
	system("cat $superScaffolds $contigs > $tmpPrefix.combined.fa");
	system("mergeWrapper.pl $tmpPrefix.combined.fa $threads strict > $prefix.combined.fa 2>>$prefix.combine.err");
	$inFile = "$prefix.combined.fa";
	($timeDiff, $cpuDiff) = timeDiff($timeDiff, $cpuDiff, "merge assembly");
} elsif ($contigs) {
	$inFile = $contigs;
} else {
	$inFile = $superScaffolds;
}

#filter contigs
unless ($noFilter) {
	unless ($noEntropy) {
		system("determineTaxonomy.pl entropy $taxaJson $inFile $whiteList > $tmpPrefix.entropy.fa 2>>$tmpPrefix.filter.err");
		$inFile = "$tmpPrefix.entropy.fa"; 
	}
	##diamond map
	system("cat $inFile | toUC.pl > $tmpPrefix.forFilter.UC.fa");
	scaffoldSplit("$tmpPrefix.forFilter.UC.fa", "$tmpPrefix.forDiamondFilter.fa", "diamond filter", 300);
	system("diamond blastx -q $tmpPrefix.forDiamondFilter.fa --db $gbBlastx -o $tmpPrefix.diamondFilter.out --algo 1 -f 6 --top 5 --tmpdir $tmpdir -e .1 --unal 0 -c 1 -p $threads --verbose 1>>$tmpPrefix.filter.err 2>>$tmpPrefix.filter.err");
	system("$compressStdout $tmpPrefix.diamondFilter.out > $prefix.diamondFilter.out.bz2");
	($timeDiff, $cpuDiff) = timeDiff($timeDiff, $cpuDiff, "diamond filter map");

	##megablast map
	scaffoldSplit("$tmpPrefix.forFilter.UC.fa", "$tmpPrefix.forBlastnFilter.split.fa", "blastn filter", 0);
	system("blastn -task blastn -word_size 39 -query $tmpPrefix.forBlastnFilter.split.fa -db $gbBlastn -max_target_seqs 1000 -num_threads $threads -outfmt 6 -out $tmpPrefix.blastnFilter.out 2>>$tmpPrefix.filter.err");
	system("$compressStdout $tmpPrefix.blastnFilter.out > $prefix.blastnFilter.out.bz2");
	($timeDiff, $cpuDiff) = timeDiff($timeDiff, $cpuDiff, "blastn filter map");

	##filter seqs
	system("determineTaxonomy.pl filter $taxaJson $inFile $prefix.diamondFilter.out.bz2 $prefix.blastnFilter.out.bz2 $threads $whiteList 1>$prefix.filtered.fa 2>>$tmpPrefix.filter.err");
	system("$compressStdout $tmpPrefix.filter.err > $prefix.filter.err.bz2");
	($timeDiff, $cpuDiff) = timeDiff($timeDiff, $cpuDiff, "filter contigs");	

	$inFile = "$prefix.filtered.fa";
}


#improve/merge superscaffolds if possible
unless ($noIterative) {
	if ($noFilter and not $noEntropy) {
		system("determineTaxonomy.pl entropy $taxaJson $inFile $whiteList > $tmpPrefix.entropy.fa 2>>$tmpPrefix.iterateImprove.err");
		$inFile = "$tmpPrefix.entropy.fa";
	}
	system("mergeWrapper.pl $inFile $threads init > $tmpPrefix.forCluster.fa 2>>$tmpPrefix.iterateImprove.err");
	scaffoldSplit("$tmpPrefix.forCluster.fa", "$tmpPrefix.forCluster.split.fa", "cluster", 0);
	system("clusterByKmer.pl $tmpPrefix.forCluster.split.fa $tmpPrefix.forCluster.fa > $tmpPrefix.forImprove.fa 2>>$tmpPrefix.iterateImprove.err");
	system("improveWrapper.pl $tmpPrefix.forImprove.fa $tmpPrefix.derep.fa $threads $gbBlastn $gbBlastx $taxaJson $improveTimeLimit $strictMerge $whiteList 1>$prefix.improved.fa 2>>$tmpPrefix.iterateImprove.err");
	system("$compressStdout $tmpPrefix.iterateImprove.err > $prefix.iterateImprove.err.bz2");
	$inFile = "$prefix.improved.fa";
	($timeDiff, $cpuDiff) = timeDiff($timeDiff, $cpuDiff, "iterative improvement");
} else {
	if (not $noMerge) {
		system("mergeWrapper.pl $inFile $threads $strictMerge > $prefix.selfMerged.fa 2>$tmpPrefix.merge.err");
		system("$compressStdout $tmpPrefix.merge.err > $prefix.merge.err.bz2");
		$inFile = "$prefix.selfMerged.fa";
	}
}

system("cat $inFile | toUC.pl all > $tmpPrefix.selfAlign.fa");
system("bbmap.sh -Xmx$ram usejni=t in=$tmpPrefix.derep.fa ref=$tmpPrefix.selfAlign.fa nodisk=t vslow=t minid=.8 noheader=t 32bit=t ignorefrequentkmers=f greedy=f excludefraction=0 deterministic=t sam=1.3 threads=$threads local=t outm=stdout ambiguous=random 2>$tmpPrefix.selfAlign.err | lbzip2 -c -n$threads > $tmpPrefix.selfAlign.sam.bz2");
system("adjustSizesAndFilter.pl quant 0 $tmpPrefix.selfAlign.sam.bz2 $inFile /dev/null $threads > $prefix.selfAlign.fa 2>>$tmpPrefix.selfAlign.err");
system("$compressStdout $tmpPrefix.selfAlign.err > $prefix.selfAlign.err.bz2");
$inFile = "$prefix.selfAlign.fa";
($timeDiff, $cpuDiff) = timeDiff($timeDiff, $cpuDiff, "self align and quantify");


#final blastn and blastx for determining taxonomy
unless ($skipTaxonomy) {
	if ($noIterative and $noFilter and not $noEntropy) {
		system("determineTaxonomy.pl entropy $taxaJson $inFile $whiteList > $tmpPrefix.entropy.fa 2>>$tmpPrefix.taxonomy.err");
		$inFile = "$tmpPrefix.entropy.fa";
	}
	system("cat $inFile | toUC.pl > $tmpPrefix.ucBlastInput.fa");
	scaffoldSplit("$tmpPrefix.ucBlastInput.fa", "$tmpPrefix.blastnInput.fa", "blastn", 0);
	#if you change blast scoring or method, you need to adjust the lambda and kappa in determineTaxonomy.pl
	system("blastn -task blastn -query $tmpPrefix.blastnInput.fa -db $gbBlastn -max_target_seqs 100000 -outfmt 6 -word_size 17 -num_threads $threads -out $tmpPrefix.blastn.out 2>$prefix.blastn.err");
	system("$compressStdout $tmpPrefix.blastn.out > $prefix.blastn.out.bz2");
	($timeDiff, $cpuDiff) = timeDiff($timeDiff, $cpuDiff, "blastn full");


	scaffoldSplit("$tmpPrefix.ucBlastInput.fa", "$tmpPrefix.diamondInput.fa", "diamond", 100);
	system("diamond blastx --algo 1 -q $tmpPrefix.diamondInput.fa --db $gbBlastx --max-target-seqs 100000 --max-hsps 100000 -p $threads --outfmt 6 --verbose -o $tmpPrefix.diamondBlastx.out -e 1e-5 --unal 0 -c 1 --masking 0 --comp-based-stats 0 -t $tmpdir --verbose 1>>$prefix.diamondBlastx.err 2>>$prefix.diamondBlastx.err");
	system("cat $tmpPrefix.diamondBlastx.out | hitMask.pl $tmpPrefix.ucBlastInput.fa $tmpPrefix.diamondInput.fa > $tmpPrefix.diamondInput.remain.fa");
	system("cat $tmpPrefix.diamondInput.remain.fa | toUC.pl > $tmpPrefix.diamondInput.remain.UC.fa");
	scaffoldSplit("$tmpPrefix.diamondInput.remain.UC.fa", "$tmpPrefix.diamondInput.remain.split.fa", "uc.remain", 0);
	system("diamond blastx --algo 1 -q $tmpPrefix.diamondInput.remain.split.fa --db $gbBlastx --max-target-seqs 100000 --max-hsps 100000 -p $threads --verbose --outfmt 6 -o $tmpPrefix.diamondBlastx.remain.out -e .1 --unal 0 -c 1 --masking 0 --comp-based-stats 0 -t $tmpdir --verbose 1>>$prefix.diamondBlastx.err 2>>$prefix.diamondBlastx.err");
	system("$compressStdout $tmpPrefix.diamondBlastx.out $tmpPrefix.diamondBlastx.remain.out > $prefix.diamondBlastx.out.bz2");
	($timeDiff, $cpuDiff) = timeDiff($timeDiff, $cpuDiff, "diamond full");
	
	#actually determine taxonomy
	system("determineTaxonomy.pl classify $infoFloor $taxaJson $inFile $prefix.diamondBlastx.out.bz2 $prefix.blastn.out.bz2 $threads $whiteList 1>$prefix.final.fa 2>>$tmpPrefix.taxonomy.err");
	system("$compressStdout $tmpPrefix.taxonomy.err > $prefix.taxonomy.err.bz2");
	($timeDiff, $cpuDiff) = timeDiff($timeDiff, $cpuDiff, "determine taxonomy");
} else {
	system("cp $inFile $prefix.final.fa");
}
timeDiff($startTime, "0", "Overall Virmap time");

sub timeDiff {
	my $t1 = $_[0];
	my $c1 = $_[1];
	my $message = $_[2];
	my $t2 = time;
#	my $ticks = sysconf('_SC_CLK_TCK');
#	my $cpuLine = `cat /proc/$$/stat`;
#	chomp $cpuLine;
#	my @cpuParts = split / /, $cpuLine;
#	my $c2 = ($cpuParts[13] + $cpuParts[14] + $cpuParts[15] + $cpuParts[16]) / $ticks;
	my $c2 = 0;
	my $timeDiff = sprintf("%.2f", ($t2 - $t1));
#	my $cpuDiff = sprintf("%.2f", ($c2 - $c1));
#	my $cpuRatio = sprintf("%.2f", ($cpuDiff / $timeDiff));
#	my $msg = "TIME $sampleName $message: $timeDiff seconds, $cpuDiff CPU seconds, $cpuRatio CPU ratio";
	my $msg = "TIME $sampleName $message: $timeDiff seconds";
	print STDERR "$msg\n";
	system("echo '$msg' >> $prefix.log");
	return($t2, $c2);
}



sub scaffoldSplit {
	my $source = $_[0];
	my $target = $_[1];
	my $message = $_[2];
	my $default = $_[3];
	system("scaffoldSplit.pl $source $default > $target");

}
