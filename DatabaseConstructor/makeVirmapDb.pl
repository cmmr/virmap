#!/usr/bin/env perl
#
#
#
use warnings;
use strict;
use Getopt::Long qw/GetOptionsFromString GetOptions/;
use File::Temp qw/tempdir/;
use Cwd;
#use DataBrowser qw(browse);
use Sys::MemInfo qw(totalmem freemem totalswap);
use Array::Shuffle qw(shuffle_array);
use Digest::SHA qw(sha256_hex);
use Time::HiRes qw(time);
use threads;



my $outputDir = "";
my $divisionsFile = "";
my $saveFasta = 0;
GetOptions("outputDir=s" => \$outputDir, "divisionsFile=s" => \$divisionsFile, "saveFasta" => \$saveFasta);



my $procs = `cat /proc/cpuinfo | grep "^processor" | tail -1 | sed "s/.* //g"`;
chomp $procs;
$procs = $procs + 1;


my $totalMem = totalmem();
if ($totalMem / (1024 ** 3) < 170) {
	die "Machine doesn't have enough memory to safely build the database\n";
}


unless ($ENV{'TMPDIR'}) {
	die "TMPDIR not set\n";
}
unless (-e $ENV{'TMPDIR'}) {
	die "TMPDIR set, but doesn't exist: $ENV{'TMPDIR'}\n";
}
my $tmpdir = tempdir( CLEANUP => 1 );
unless (testFreeSpaceGb($tmpdir) > 700) {
	die "$tmpdir has less than 700GB of free space\n";
}
if ($saveFasta and testFreeSpaceGb($tmpdir) < 2000) {
	die "$tmpdir has less than 2000GB of free space (saveFastaMode)\n";
}

my $startDir = cwd();
 

my $ncbiGenbankFtp = "ftp://ftp.ncbi.nlm.nih.gov/genbank/";


my @wantedDivisions;
my @defaultDivisions = qw(gbbct gbcon gbenv gbhtg gbinv gbmam gbpat gbphg gbpln gbpri gbrod gbsyn gbuna gbvrl gbvrt);
my @viralDivisions = qw(gbvrl gbphg);
my $viralDivisions = {};
foreach my $div (@viralDivisions) {
	$viralDivisions->{$div} = 1;
}



unless ($outputDir) {
	die "--outputDir not defined\n";
}
unless (-e $outputDir) {
	die "$outputDir does not exist\n";
}
unless (testFreeSpaceGb($outputDir) > 300) {
	die" $outputDir does not have 300GB of free space\n";
}
my $absOutputDir = `readlink -e $outputDir`;
chomp $absOutputDir;
$outputDir = $absOutputDir;


if (-f $divisionsFile) {
	my $temp = {};
	open IN, "$divisionsFile";
	while (my $line = <IN>) {
		chomp $line;
		if (testValidDivision($line)) {
			$temp->{$line} = 1;
		} else {
			print STDERR "No valid .seq.gz files for wanted division $line\n";
		}
	}
	@wantedDivisions = keys %$temp;
	my $printOutDiv = join " ", @wantedDivisions;
	print STDERR "Using divisions from file $divisionsFile: $printOutDiv\n";
} else {
	print STDERR "Using default divisions\n";
	foreach my $div (@defaultDivisions) {
		if (testValidDivision($div)) {
			push @wantedDivisions, $div;
		} else {
			print STDERR "Default division $div no longer has any files, please update the build script\n";
		}
		
	}
}
my $hasVirus = 0; 
foreach my $div (@wantedDivisions) {
	if (exists $viralDivisions->{$div}) {
		$hasVirus = 1;
	}
}

unless ($hasVirus) {
	die "no viral divisions in wanted divisions\n";
}


chdir($tmpdir);

#### Get raw files from Genbank ####

my $files = {};
foreach my $div (@wantedDivisions) {
	my $return = `lftp -e "rels $div*.seq.gz; exit" ftp://ftp.ncbi.nlm.nih.gov/genbank/ 2>/dev/null | sed "s/.* //g"`;
	chomp $return;
	my @files = split /\n/, $return;
	$files->{$div} = \@files;

}

my @includeString;
foreach my $div (@wantedDivisions) {
	push @includeString, "--include-glob=$div*.seq.gz";
}
my $includeString = join " ", @includeString;
system("lftp -e \"mirror -r --parallel=8 $includeString; exit\" $ncbiGenbankFtp 2>/dev/null");
print STDERR "Finished grabbing division dumps from $ncbiGenbankFtp\n";

my $fileCheckFail = 0;
foreach my $div (keys %$files) {
	foreach my $file (@{$files->{$div}}) {
		unless (-e "$tmpdir/$file") {
			print STDERR "$file missing from division $div\n";
			$fileCheckFail = 1;
		}
	}
}
if ($fileCheckFail) {
	die "not all expected files present\n";
}


#### Make Viral databases ####
open PROT, ">$tmpdir/vrl.prot.commands";
open PROTGB, ">$tmpdir/vrl.gb.prot.commands";
open NUC, ">$tmpdir/vrl.nuc.commands";
foreach my $div (@viralDivisions) {
	foreach my $file (@{$files->{$div}}) {
		print NUC "pigz -dc $tmpdir/$file | genbankToNuc.pl noSum\n";
		print PROTGB "pigz -dc $tmpdir/$file | genbankToProt.pl noPos noSum | paste - -\n";
		print PROT "pigz -dc $tmpdir/$file | genbankToProt.pl noSum\n";
	}
}
close PROTGB;
close PROT;
close NUC;


#my $vrlProtThr = threads->create(\&makeVrlFaa, $tmpdir);
#my $vrlGbProtThr = threads->create(\&makeGbVrlFaa, $tmpdir);
#my $vrlNucThr = threads->create(\&makeVrlFna, $tmpdir);


my $halfprocs = int($procs / 2);

#### make generic genbank databases ####
open PROT, ">$tmpdir/prot.commands";
open NUC, ">$tmpdir/nuc.commands";
open PROTS, ">$tmpdir/prot.speed.commands";
open PROTS2, ">$tmpdir/prot.speed.commands.2";
open NUCS, ">$tmpdir/nuc.speed.commands";
open NUCS2, ">$tmpdir/nuc.speed.commands.2";
open NUCS3, ">$tmpdir/nuc.speed.commands.3";
my @speedFiles;
my @toDelete;
my $devShmTmp = tempdir( DIR => "/dev/shm", CLEANUP => 1);
my $fileOrder = {};
foreach my $div (@wantedDivisions) {
	if (exists $viralDivisions->{$div}) {
		next;
	}
	foreach my $file (@{$files->{$div}}) {
		my $fileSize = -s "$tmpdir/$file";
		$fileOrder->{$file} = $fileSize;
	}
}
foreach my $file (sort {$fileOrder->{$b} <=> $fileOrder->{$a}} keys %$fileOrder) {
	push @speedFiles, "$tmpdir/$file";
	print PROT "pigz -dc $tmpdir/$file | genbankToProt.pl noPos\n";
	print NUC "pigz -dc $tmpdir/$file | genbankToNuc.pl\n";
	print PROTS "pigz -dc $tmpdir/$file | genbankToProt.pl noPos pasteMode | tee >(zstd -cq > $tmpdir/$file.aa.faa.zst) | cut -f2\n";
	print PROTS2 "zstd -dcq $tmpdir/$file.aa.faa.zst | checksumDerepPreLoad.pl $devShmTmp/aa.dupes $tmpdir/$file.dupes.faa.zst $tmpdir/$file.uniq.faa.zst\n";
	push @toDelete, "$tmpdir/$file.aa.faa.zst";
	print NUCS "pigz -dc $tmpdir/$file | genbankToNuc.pl pasteMode | tee >(zstd -cq > $tmpdir/$file.nuc.fna.zst) | cut -f2\n";
	print NUCS2 "zstd -dcq $tmpdir/$file.nuc.fna.zst | checksumDerepPreLoad.pl $devShmTmp/nuc.dupes $tmpdir/$file.dupes.fna.zst $tmpdir/$file.uniq.fna.zst\n";
	print NUCS3 "zstd -dcq $tmpdir/$file.nuc.fna.zst | cut -f1 | sed 's/.*taxId=//g'\n";
	push @toDelete, "$tmpdir/$file.nuc.fna.zst";
}

close PROT;
close NUC;
close PROTS;
close NUCS;
close PROTS2;
close NUCS2;
#my $makeGbBlastxRawThr = threads->create(\&makeGbBlastxFaa, $tmpdir, $halfprocs);
#my $makeGbBlastnRawThr = threads->create(\&makeGbBlastnFna, $tmpdir, $halfprocs, $outputDir);
my $toDelete = join " ", @toDelete;
system("bash", "-c", "cat $tmpdir/nuc.speed.commands | parallel | sort --buffer-size=30G --parallel=$halfprocs | uniq -c | sed -re 's/^\\s+//g' | grep -v '^1 ' | cut -f2 -d ' ' | makeDupeStruct.pl > $devShmTmp/nuc.dupes");


my $makeGbBlastnRawThr = threads->create(\&makeGbBlastnFnaSpeed, $tmpdir, $procs, \@speedFiles, $devShmTmp, $outputDir);
#my $makeGbBlastxRawThr = threads->create(\&makeGbBlastxFaaSpeed, $tmpdir, $halfprocs, \@speedFiles, $devShmTmp);


my $makeGbBlastxRawThr;
my $vrlProtThr;
my $vrlGbProtThr;
my $vrlNucThr;





#my $vrlProtThr = threads->create(\&makeVrlFaa, $tmpdir);
#my $vrlGbProtThr = threads->create(\&makeGbVrlFaa, $tmpdir);
#my $vrlNucThr = threads->create(\&makeVrlFna, $tmpdir, $outputDir);
my $virProtReady = 0;
my $virGbProtReady = 0;
my $virNucReady = 0;
my $gbBlastxReady = 0;
my $gbBlastnReady = 0;
my $combinedGbBlastnReady = 0;

my $gbBlastxDone = 0;
my $gbBlastnDone = 0;
my $virDmndDone = 0;
my $virBbmapDone = 0;
my $taxonomyDone = 0;
my $deleteDone = 0;

my $virDmndThr;
my $virDmndRunning = 0;
my $virBbmapThr;
my $virBbmapRunning = 0;
my $gbBlastxThr;
my $gbBlastxRunning = 0;
my $combineGbBlastnThr;
my $combineGbBlastnRunning = 0;
my $gbBlastnThr;
my $gbBlastnRunning = 0;
my $taxonomyThr;
my $taxonomyRunning = 0;
my $deleteThr;
my $deleteRunning = 0;
my $stage2 = 0;
my $stage2Running = 0;

my $extraFiles;
while (1) {
	sleep 1;
	if ($stage2 and not $stage2Running) {
		$makeGbBlastxRawThr = threads->create(\&makeGbBlastxFaaSpeed, $tmpdir, $halfprocs, \@speedFiles, $devShmTmp);
		$vrlProtThr = threads->create(\&makeVrlFaa, $tmpdir);
		$vrlGbProtThr = threads->create(\&makeGbVrlFaa, $tmpdir);
		$vrlNucThr = threads->create(\&makeVrlFna, $tmpdir, $outputDir);
		$stage2Running = 1;
		next;
	}
	if ($stage2 and not $gbBlastxReady and $makeGbBlastxRawThr->is_joinable()) {
		$makeGbBlastxRawThr->join();
		$gbBlastxReady = 1;
		print STDERR "Finished parsing for gbBlastx\n";
	}
	if (not $gbBlastnReady and $makeGbBlastnRawThr->is_joinable()) {
		$extraFiles = $makeGbBlastnRawThr->join();
		$gbBlastnReady = 1;
		print STDERR "Finished parsing for gbBlastn\n";
	}
	if ($stage2 and $gbBlastnReady and $gbBlastxReady and not $deleteRunning) {
		$deleteThr = threads->create(\&deleteFiles, $toDelete);
		$deleteRunning = 1;
	}
	if ($stage2 and not $deleteDone and $deleteRunning and $deleteThr->is_joinable()) {
		$deleteThr->join();
		$deleteDone = 1;
	}
	if ($stage2 and not $virProtReady and $vrlProtThr->is_joinable()) {
		$vrlProtThr->join();
		$virProtReady = 1;
		print STDERR "Finished parsing virus for virDmnd.dmnd\n";
	}
	if ($stage2 and not $virGbProtReady and $vrlGbProtThr->is_joinable()) {
		$vrlGbProtThr->join();
		$virGbProtReady = 1;
		print STDERR "Finished parsing virus for virDmnd.dmnd\n";
	}
	if ($stage2 and not $virNucReady and $vrlNucThr->is_joinable) {
		$vrlNucThr->join();
		$virNucReady = 1;
		print STDERR "Finished parsing virus for virBbmap\n";
	}

	if ($stage2 and $virProtReady and not $virDmndRunning) {
		$virDmndThr = threads->create(\&makeVirdmnd, $tmpdir, $outputDir, $halfprocs, $saveFasta);
		$virDmndRunning = 1;
	}
	if ($stage2 and not $virDmndDone and $virDmndRunning and $virDmndThr->is_joinable()) {
		$virDmndThr->join();
		$virDmndDone = 1;
		print STDERR "Finished building virDmnd.dmnd\n";
	}

	if ($stage2 and $virNucReady and not $virBbmapRunning) {
		$virBbmapThr = threads->create(\&makeVirBbmap, $tmpdir, $outputDir, $halfprocs, $saveFasta);
		$virBbmapRunning = 1;
	}
	if ($stage2 and not $virBbmapDone and $virBbmapRunning and $virBbmapThr->is_joinable()) {
		$virBbmapThr->join();
		$virBbmapDone = 1;
		print STDERR "Finished building virBbmap\n";
	}

	if ($stage2 and $gbBlastxReady and $virGbProtReady and not $gbBlastxRunning) {
		$gbBlastxThr = threads->create(\&makeGbBlastx, $tmpdir, $outputDir, $halfprocs, $saveFasta);
		$gbBlastxRunning = 1;
	}
	if ($stage2 and not $gbBlastxDone and $gbBlastxRunning and $gbBlastxThr->is_joinable()) {
		$gbBlastxThr->join();
		$gbBlastxDone = 1;
		print STDERR "Finished building gbBlastx.dmnd\n";
	}
	if ($gbBlastnReady and not $gbBlastnRunning) {
		$gbBlastnThr = threads->create(\&makeGbBlastn, $tmpdir, $extraFiles, $outputDir, $procs, $saveFasta);
		$stage2 = 1;
		$gbBlastnRunning = 1;
	}
	if ($stage2 and $gbBlastnReady and $virNucReady and not $taxonomyRunning) {
		$taxonomyThr = threads->create(\&makeTaxonomy, $tmpdir, $outputDir);
		$taxonomyRunning = 1;
	}
	if ($stage2 and not $gbBlastnDone and $gbBlastnRunning and $gbBlastnThr->is_joinable()) {
		$gbBlastnThr->join();
		$gbBlastnDone = 1;
		print STDERR "Finished building gbBlastn\n";
	}
	if ($stage2 and not $taxonomyDone and $taxonomyRunning and $taxonomyThr->is_joinable()) {
		$taxonomyThr->join();
		$taxonomyDone = 1;
		print STDERR "Finished building taxonomy structure\n";
	}
	if ($stage2 and $gbBlastnDone and $gbBlastxDone and $virDmndDone and $virBbmapDone and $taxonomyDone) {
		my $checksum = sha256_hex(time());
		open CS, ">$outputDir/VirmapDatabases.checksum";
		print CS "$checksum\n";
		close CS;
		last;
	}


}



print STDERR "Finished building all virmap databases\n";


chdir($startDir);

sub deleteFiles {
	my $toDelete = $_[0];
	system("rm $toDelete");
	return(0);
}

sub makeVirBbmap {
	my $tmpdir = $_[0];
	my $outdir = $_[1];
	my $halfprocs = $_[2];
	my $saveFasta = $_[3];
	system("bbmap.sh -Xmx70G ref=$tmpdir/viral.nuc.fna path=$outdir/virBbmap build=1 2>$outdir/makeVirBbmap.err");
	if ($saveFasta) {
		system("cat $tmpdir/viral.nuc.fna | lbzip2 -c -n $halfprocs > $outdir/viral.nuc.fna.bz2");
	}
	return(0);
}

sub makeVirdmnd {
	my $tmpdir = $_[0];
	my $outdir = $_[1];
	my $halfprocs = $_[2];
	my $saveFasta = $_[3];
	system("diamond makedb --in $tmpdir/viral.aa.faa --db $outdir/virDmnd.dmnd 2>$outdir/makeVirDmnd.err");
	if ($saveFasta) {
		system("cat $tmpdir/viral.aa.faa | lbzip2 -c -n $halfprocs > $outdir/viral.aa.faa.bz2");
	}
	return(0);
}

sub makeVrlFaa {
	my $tmpdir = $_[0];
	system("cat $tmpdir/vrl.prot.commands | parallel | paste - - | shuf | tr '\\t' '\\n' > $tmpdir/viral.aa.faa");
	return(0);
}
sub makeGbVrlFaa {
	my $tmpdir = $_[0];
	system("cat $tmpdir/vrl.gb.prot.commands | parallel | zstd -cq -T$halfprocs> $tmpdir/viral.gb.aa.pre.zst");
	return(0);
}
sub makeVrlFna {
	my $tmpdir = $_[0];
	my $outdir = $_[1];
	system("bash", "-c", "cat $tmpdir/vrl.nuc.commands | parallel | paste - - | shuf | tee >(zstd -cq > $tmpdir/viral.nuc.pre.zst) | tr '\\t' '\\n' > $tmpdir/viral.nuc.fna");
	system("cat $tmpdir/viral.nuc.fna | grep '^>' | sed 's/.*taxId=//g' | sort --buffer-size=10G --parallel=8 | uniq > $tmpdir/viral.taxids");
	system("cp $tmpdir/viral.taxids $outdir/");
	return(0);
}

sub makeGbBlastxFaa {
	my $tmpdir = $_[0];
	my $halfprocs = $_[1];
	system("cat $tmpdir/prot.commands | parallel -j $halfprocs --lb | checksumDerep.pl | paste - - > $tmpdir/gb.aa.pre");
	return(0);
}
sub makeGbBlastxFaaSpeed {
	my $tmpdir = $_[0];
	my $halfprocs = $_[1];
	my @inFiles = @{$_[2]};
	my $devShmTmp = $_[3];
	system("bash", "-c", "cat $tmpdir/prot.speed.commands | parallel -j $halfprocs | sort --buffer-size=30G --parallel=$halfprocs | uniq -c | sed -re 's/^\\s+//g' | grep -v '^1 ' | cut -f2 -d ' ' | makeDupeStruct.pl > $devShmTmp/aa.dupes");
	system("cat $tmpdir/prot.speed.commands.2 | parallel -j $halfprocs");	
	my @dupeFiles;
	my @uniqFiles;
	foreach my $file (@inFiles) {
		push @dupeFiles, "$file.dupes.faa.zst";
		push @uniqFiles, "$file.uniq.faa.zst";
	}
	my $dupesLine = join " ", @dupeFiles;
	my $uniqLine = join " ", @uniqFiles;
	system("zstd -dcq $dupesLine | checksumDerepDupesOnly.pl | zstd -cq -T$halfprocs > $tmpdir/gb.aa.pre.zst");
	system("rm $dupesLine");
	system("zstd -dcq $uniqLine | zstd -cq -T$halfprocs >> $tmpdir/gb.aa.pre.zst");
	system("rm $uniqLine");
	return(0);
}
sub makeGbBlastx {
	my $tmpdir = $_[0];
	my $outdir = $_[1];
	my $halfprocs = $_[2];
	my $saveFasta = $_[3];
	system("bash", "-c", "zstd -dcq $tmpdir/gb.aa.pre.zst $tmpdir/viral.gb.aa.pre.zst | shuf | tr '\\t' '\\n' | zstd -cq -T$halfprocs > $tmpdir/gb.aa.faa.zst; rm $tmpdir/gb.aa.pre.zst $tmpdir/viral.gb.aa.pre.zst");
	system("zstd -cq $tmpdir/gb.aa.faa.zst | diamond makedb --in /dev/stdin --db $outdir/gbBlastx.dmnd 2>$outdir/makeGbBlastx.err");
	if ($saveFasta) {
		system("zstd -cq $tmpdir/gb.aa.faa.zst | lbzip2 -c -n $halfprocs > $outdir/gb.aa.faa.bz2");
	}
	system("rm $tmpdir/gb.aa.faa.zst");
	return(0);
}

sub makeGbBlastnFna {
	my $tmpdir = $_[0];
	my $halfprocs = $_[1];
	my $outdir = $_[2];
	system("cat $tmpdir/nuc.commands | parallel -j $halfprocs --lb | checksumDerep.pl | paste - - | tee $tmpdir/gb.nuc.pre | cut -f1 | sed 's/.*taxIx=//g' | sort --buffer-size=30G --parallel=8 | uniq > $tmpdir/gb.taxids");
	system("cp $tmpdir/gb.taxids $outdir/");
	return(0);
}
sub makeGbBlastnFnaSpeed {
	my $tmpdir = $_[0];
	my $subThreads = $_[1];
	my @inFiles = @{$_[2]};
	my $outdir = $_[3];
	system("cat $tmpdir/nuc.speed.commands.2 | shuf | parallel -j $subThreads");
	system("cat $tmpdir/nuc.speed.commands.3 | parallel -j $subThreads | sort --buffer-size=30G --parallel=$halfprocs | uniq > $tmpdir/gb.taxids");
	system("cp $tmpdir/gb.taxids $outdir/");
	my @dupeFiles;
	my @uniqFiles;
	foreach my $file (@inFiles) {
		push @dupeFiles, "$file.dupes.fna.zst";
		push @uniqFiles, "$file.uniq.fna.zst";
	}
	shuffle_array(@uniqFiles);
	my $dupesLine = join " ", @dupeFiles;
	my $uniqLine = join " ", @uniqFiles;
	system("zstd -dcq $dupesLine | checksumDerepDupesOnly.pl | zstd -cq -T$subThreads > $tmpdir/gb.nuc.pre.zst");
	system("rm $dupesLine");
#	system("cat $uniqLine >> $tmpdir/gb.nuc.pre");
#	system("rm $uniqLine");
	return(\@uniqFiles);
}


#sub combineGbBlastn {
sub makeGbBlastn {
	my $tmpdir = $_[0];
	my $extra = $_[1];
	my $outdir = $_[2];
	my $procs = $_[3];
	my $saveFasta = $_[4];
	my $extraFiles = "";
	if ($extra) {
		$extraFiles = join " ", @{$extra};
	}
	my $saveCmd = "";
	if ($saveFasta) {
		$saveCmd = "| tee >(lbzip2 -c -n $halfprocs > $outdir/gb.nuc.fna.bz2)";
	}
	system("bash", "-c", "zstd -dqc $tmpdir/gb.nuc.pre.zst $extraFiles $tmpdir/viral.nuc.pre.zst | tr '\\t' '\\n' $saveCmd | makeblastdb -title gbBlastn -out $outdir/gbBlastn -dbtype nucl -max_file_sz 2000000000 1>$outdir/makeGbBlastn.err 2>&1; rm $tmpdir/gb.nuc.pre.zst $tmpdir/viral.nuc.pre.zst");
	if ($extra) {
		system("rm $extraFiles");
	}
#	return(0);
#}
#sub makeGbBlastn {
#	my $tmpdir = $_[0];
#	my $outdir = $_[1];
#	system("makeblastdb -out $outdir/gbBlastn -in $tmpdir/gb.nuc.fna -dbtype nucl -max_file_sz 2000000000 1>$outdir/makeGbBlastn.err 2>&1");
	return(0);
}
sub makeTaxonomy {
	my $tmpdir = $_[0];
	my $outdir = $_[1];
	system("wget -q --no-check-certificate https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz -O $tmpdir/taxdump.tar.gz");
	unless (-e "$tmpdir/taxdump.tar.gz" and -s "$tmpdir/taxdump.tar.gz") {
		print STDERR "WARNING: unable to grab taxdump.tar.gz, please manually make taxonomy structure\n";
		return(0);
	}
	system("tar zxf $tmpdir/taxdump.tar.gz -C $tmpdir/");
	system("cat $tmpdir/gb.taxids $tmpdir/viral.taxids > $tmpdir/taxaIds");
	system("makeStruct.pl $tmpdir/nodes.dmp $tmpdir/names.dmp $tmpdir/merged.dmp $tmpdir/delnodes.dmp $tmpdir/taxaIds > $outdir/Taxonomy.virmap 2>$outdir/makeTaxaStruct.err");
	return(0);
}

sub testValidDivision {
	my $div = $_[0];
	my $valid = `lftp -e "rels $div*.seq.gz; exit" ftp://ftp.ncbi.nlm.nih.gov/genbank/ 2>/dev/null`;
	chomp $valid;
	if (length($valid)) {
		return 1;
	} else {
		return 0;
	}

}

sub testFreeSpaceGb {
	my $testDir = $_[0];
	my $statLines = `stat -f $tmpdir`;
	my $blockSize;
	my $freeBlocks;
	my @statLines = split /\n/, $statLines;
	foreach my $line (@statLines) {
		chomp $line; #might be unneccessary
		if ($line =~ m/Block size: (\d+)/) {
			$blockSize = $1;
		}
		if ($line =~ m/Blocks: Total:.*Available: (\d+)/) {
			$freeBlocks = $1;
		}
	}
	my $gb = (1024 ** 3);
	my $gbFree = ($blockSize * $freeBlocks) / $gb;
	return($gbFree);
	


}
