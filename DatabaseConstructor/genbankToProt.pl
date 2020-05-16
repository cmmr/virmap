#!/usr/bin/env perl

use warnings;
use strict;
use Digest::MD5 qw(md5_hex);
use Digest::SHA qw(sha1_hex);


open IN, "-";
my $yearLimit = "inf";
my $noPos = 0;
my $noSum = 0;
my $pasteMode = 0;
my $args = {};
foreach my $elem (@ARGV) {
	if ($elem =~ m/\d\d\d\d/) {
		$yearLimit = $elem;
	} else {
		$args->{$elem} = 1
	}
}
if (exists $args->{"noPos"}) {
	$noPos = 1;
}
if (exists $args->{"noSum"}) {
	$noSum = 1;
}
if (exists $args->{"pasteMode"}) {
	$pasteMode = 1;
}
my $min = 30;
my $grabdesc = 0;
my $grabaa = 0;
my $grabpos = 0;
my $grabdef = 0;
my $accession = "";
my $desc = "";
my $translation = "";
my $protId = "";
my $dbref = "";
my $prod = "";
my $def = "";
my $grabProd = 0;
my $pos = "";
my $taxaId = 1;
my $soloTaxaId = 0;
my $reject = 0;
my $codonStart = 0;
my $badTime = 0;
my $date = "";
my $defLine = "";
while (my $line = <IN>) {
	chomp $line;
	if ($line =~ m/^LOCUS.*\s([^\s]+)$/ and $yearLimit) {
		$date = $1;
		my @parts = split /-/, $date;
		
		if (defined $parts[2] and $parts[2] > $yearLimit) {
			$badTime = 1;
		}
		if (not defined $parts[2]) {
			print STDERR "$line has no year\n";
		}
		next;
	}
	if ($line =~ m/^\/\/$/) {
		$grabdesc = 0;
		$grabaa = 0;
		$grabpos = 0;
		$grabdef = 0;
		$accession = "";
		$desc = "";
		$translation = "";
		$protId = "";
		$dbref = "";
		$prod = "";
		$def = "";
		$pos = "";
		$taxaId = 1;
		$soloTaxaId = 0;
		$grabProd = 0;
		$reject = 0;
		$codonStart = 0;
		$badTime = 0;	
		next;
	}
	if ($line =~ m/^SOURCE\s+(.*)/) {
		$desc = $1;
		next
	}
	if ($line =~ m/^\s+\/codon_start=(\d+)/) {
		$codonStart = $1;
		next;
	}
	if ($line =~ m/^DEFINITION\s+(.*)/) {
		$soloTaxaId = 0;
		$reject = 0;
		$def = $1;
		$grabdef = 1;
		next;
	}
	if ($line =~ m/^ACCESSION/) {
		$grabdef = 0;
		next;
	}
	if ($grabdef) {
		$line =~ m/\s+(.*)/;
		unless ($1) {
			next;
		}
		$def = $def . " $1";
		next;
	}
	if ($line =~ m/^     CDS\s+(\S+)/) {
		$pos = $1;
		$grabpos = 1;
		next;
	}
	if ($grabpos) {
		if ($line =~ m/^\s+\//) {
			$grabpos = 0;
			redo;
		}
		$line =~ m/^\s+(.*)/;
		$pos .= $1;
		next;
	} 
	if ($line =~ m/^VERSION\s+(\S+)/) {
		$accession = "NEWPROTDB|$1";
		next;
	}
	if ($line =~ m/^\s+\/protein_id="(\S+)"/){
		$protId = $1;
		next;
	}
	if ($line =~ m/^\s+\/db_xref="taxon:(\d+)/) {
		$taxaId = $1;
		next;
	}
	if ($line =~ m/product/) {
	}
	if ($line =~ m/^\s+\/product="(.*)/) {
		$prod = $1;
		if ($prod =~ m/"$/) {
			$prod =~ s/"$//g;
		} else {
			$grabProd = 1;
		}
		next;
	}
	if ($grabProd) {
		$line =~ m/^\s+(.*)/;
		$prod .= " $1";
		if ($prod =~ m/"$/) {
			$grabProd = 0;
			$prod =~ s/"$//g;
		}
		next;
	}
	if ($line =~ m/^\s+\/translation="(\S+)/) {
		my $aa = $1;
		$desc =~ s/\W/./g;
		$desc =~ s/\.+/./g;
		$prod =~ s/\W/./g;
		$prod =~ s/\.+/./g;
		my $posStr = "pos=$pos;codonStart=$codonStart;";
		if ($noPos) {
			$posStr = "";
		}
		$defLine = ">GI|$accession|$protId|$prod|$desc;${posStr}taxId=$taxaId";
		if ($aa =~ m/"/) {
			$aa =~ s/"//g;
			if (not $reject and not $badTime and length($aa) >= $min) {
				unless ($noSum) {
					$defLine = addSuffix($defLine, $aa);
				}
				unless ($pasteMode) {
					print "$defLine\n";
					print "$aa\n";
				} else {
					print "$defLine\t$aa\n";
				}
			}
			next;
		} else {
			$translation = $aa;
			$grabaa = 1;
			next
		}
	}
	if ($grabaa) {
		$line =~ s/^\s+//g;
		if ($line =~ m/"/) {
			$line =~ s/"//g;
			$translation = $translation . $line;
			$grabaa = 0;
			if (not $reject and not $badTime and length($translation) >= $min) {
				unless ($noSum) {
					$defLine = addSuffix($defLine, $translation);
				}
				unless ($pasteMode) {
					print "$defLine\n";
					print "$translation\n";
				} else {
					print "$defLine\t$translation\n";
				}
			}
			next;
		} else {
			$translation = $translation . $line;
			next;
		}
	}

}
sub addSuffix {
	my $def = $_[0];
	my $aa = $_[1];
	my $md5 = md5_hex($aa);
	my $sha = sha1_hex($aa);
	my $tempSum = "$md5$sha";
	my $sum = md5_hex($tempSum);
	$def = "$def\t$sum";
	return $def;

}
