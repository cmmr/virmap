#!/usr/bin/env perl


use warnings;
use strict;
use FAlite;



#preps contigs from tadpole and megahit to have a tadpole-esque name
#renames contigs to prevent collisions

open IN, "$ARGV[0]";


my $fasta_file = new FAlite(\*IN); # or any other filehandle
my $contigCount = 0;
while (my $entry = $fasta_file->nextEntry) {
	my $def = $entry->def;
	my $seq = $entry->seq;
	my $gc = () = $seq =~ m/[GC]/g;
	if ($def =~ m/^>contig_/) {
		$def =~ s/,/;/g;
		$def =~ s/contig_\d+/contig_$contigCount;src=tadpole/;
		unless ($def =~ m/;taxId=/) {
			$def =~ s/$/;taxId=10239/g;
		}
	}
	if ($def =~ m/^>k\d+_.*multi=(\d+\.?\d?).*len=(\d+)/) {
		my $cov = $1;
		my $length = $2;
		my $gcContent = $gc / $length;
		$gcContent = sprintf("%.3f", $gcContent);
		$def = ">contig_$contigCount;src=megahit;length=$length;cov=$cov;gc=$gcContent;taxId=10239";
	}
	$contigCount++;
	print "$def\n";
	print "$seq\n";
}

