#!/usr/bin/env perl
#
use warnings;
use strict;
use Digest::MD5 qw(md5_hex);
use Digest::SHA qw(sha1_hex);

my $in = "-";
my $yearLimit = "inf";
my $args = {};
my $noSum = 0;
my $pasteMode = 0;
foreach my $elem (@ARGV) {
	if ($elem =~ m/^\d\d\d\d$/) {
		$yearLimit = $elem;
	} else {
		$args->{$elem} = 1;
	}
}
if (exists $args->{"noSum"}) {
	$noSum = 1;
}
if (exists $args->{"pasteMode"}) {
	$pasteMode = 1;
}
open IN, "$in";
my $print = 0;
my $canprint = 0;
my $qr = qr/^ORIGIN/;
my $qr1 = qr/^\/\/$/;
my $sqr = qr/[\s+|\d+]/;
my $version = qr/^VERSION\s+(\S+)/;
my $definition = qr/^DEFINITION\s+(.*)/;
my $definitionLine = "";
my $accQR = qr/^ACCESSION/;
my $inSourceMetadata = 0;
my $note = "";
my $grabNote = 0;
my $min = 80;
my $defStr = "";
my $taxIdQr = qr/db_xref="taxon:(\d+)/;
my $taxId = 1;
my $indef = 0;
my $seq = "";
my $badTime = 0;
while (my $line = <IN>) {	
	chomp $line;
	if ($line =~ m/^LOCUS.*\s([^\s]+)$/ and $yearLimit) {
		my $date = $1;
		my @parts = split /-/, $date;
		if (defined $parts[2] and $parts[2] > $yearLimit) {
			$badTime = 1;
		}
	}
	if ($line =~ m/$accQR/) {
		$indef = 0;
		$canprint = 1;
		next;
	}
	if ($indef) {
		$line =~ s/^\s+//g;
		$definitionLine .= " $line";
		next;
	}
	if ($line =~ m/$definition/i) {
		$definitionLine = $1;
		$indef = 1;
		next
	}
	if ($line =~ m/$taxIdQr/) {
		$taxId = $1;
	}
	if ($line =~ m/$version/) {
		my $acc = $1;
		if ($canprint) {
			$definitionLine =~ s/\W/./g;
			$definitionLine =~ s/\.+/./g;
			$defStr = ">gi|NEWDB|gb|$acc|$definitionLine";
		}
		next;
	}
	if ($line =~ m/$qr/ and $canprint) {
		if ($taxId) {
			$defStr = $defStr . ";taxId=$taxId";
		} else {
			$defStr = $defStr . ";taxId=1";
		}
		$defStr =~ s/\s+/./g;
		$defStr =~ s/\.+/./g;
		$print = 1;
		next;
	}
	if ($line =~ m/$qr1/) {
		if (length($seq) >= $min and not $badTime) {
			my $suffix = "";
			unless ($noSum) {
				my $md5 = md5_hex($seq);
				my $sha = sha1_hex($seq);
				my $tempSum = "$md5$sha";
				my $sum = md5_hex($tempSum);
				$suffix = "\t$sum";
			}
			unless ($pasteMode) {
				print "$defStr$suffix\n";
				print "$seq\n";
			} else {
				print "$defStr$suffix\t$seq\n";
			}
		}
		$badTime = 0;
		$print = 0;
		$canprint = 0;
		$taxId = 1;
		$defStr = "";
		$seq = "";
		next;
	}
	if ($print and $canprint) {
		$line =~ s/$sqr//g;
		$line = uc($line);
		$line =~ s/[^ACGTN]/N/g;
		$seq .= $line;
	}
	


}
