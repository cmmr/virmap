#!/usr/bin/env perl
#
#
use warnings;
use strict;
use Sereal;
use Compress::Zstd qw(compress);;
use XML::Hash::XS qw();
#use DataBrowser qw(browse);
use Array::Split qw( split_by split_into );


my $data = {};
my $encoder = Sereal::Encoder->new();

open IN, "$ARGV[0]";
while (my $line = <IN>) {
	chomp $line;
	my @parts = split /\t/, $line;
	my $self = $parts[0];
	my $parent = $parts[2];
	my $rank = $parts[4];
	$data->{'parents'}->{$self} = $parent;
	$data->{'ranks'}->{$self} = $rank;
	$data->{'chilren'}->{$parent}->{$self} = 1;
}
close IN;


open IN, "$ARGV[1]";
while (my $line = <IN>) {
	chomp $line;
	my @parts = split /\t/, $line;
	unless ($parts[6] eq "scientific name") {
		next;
	}	
	$data->{'names'}->{$parts[0]} = $parts[2];
}
close IN;
foreach my $node (keys %{$data->{'parents'}}) {
	unless (exists $data->{'names'}->{$node}) {
		print STDERR "$node has no name\n";
	}
}


my $merged = {};
open IN, "$ARGV[2]";
while (my $line = <IN>) {
	chomp $line;
	my @parts = split /\t/, $line;
	my $orig = $parts[0];
	my $new = $parts[2];
	$merged->{$orig} = $new;
}
close IN;

my $deleted = {};
open IN, "$ARGV[3]";
while (my $line = <IN>) {
	chomp $line;
	my @parts = split /\t/, $line;
	$deleted->{$parts[0]} = 1;
}

my $missing = {};



open IN, "$ARGV[4]";
my $conv = XML::Hash::XS->new(utf8 => 0, encoding => 'utf-8');
while (my $line = <IN>) {
	chomp $line;
	unless (defined($line) and $line) {
		next;
	}
	unless ($line > 0) {
		next;
	}
	if (exists $merged->{$line} and exists $data->{'parents'}->{$merged->{$line}}) {
		$data->{'parents'}->{$line} = $data->{'parents'}->{$merged->{$line}};
		$data->{'names'}->{$line} = $data->{'names'}->{$merged->{$line}};
		$data->{'ranks'}->{$line} = $data->{'ranks'}->{$merged->{$line}};
		if (exists $data->{'children'}->{$merged->{$line}}) {
			$data->{'children'}->{$line} = $data->{'children'}->{$merged->{$line}};
		}
		next;
	}
	if (exists $deleted->{$line}) {
		makeDeleted($line);
	}
	unless (exists $data->{'parents'}->{$line}) {
		$missing->{$line} = 1;
	}
}
my @missingOverall = keys %$missing;
my @missing;
if (scalar(@missingOverall) and scalar @missingOverall < 200) {
	push @missing, \@missingOverall;
} else {
	@missing = split_by(200, @missingOverall);
}
if (scalar(@missing)) {
	fillInMissing($data, \@missing);
}
foreach my $node (@missingOverall) {
	unless (exists $data->{'parents'}->{$node}) {
		print STDERR "$node still missing\n";
	}
	makeDeleted($node);
}
sub makeDeleted {
	my $node = $_[0];
	$data->{'parents'}->{$node} = 1;
	$data->{'names'}->{$node} = "DELETED_ENTRY";
	$data->{'ranks'}->{$node} = "DELETED_RANK";
}	
sub fillInMissing {
	my $data = $_[0];
	my @missing = @{$_[1]};
	my $neededParents = {};
	foreach my $arrayRef (@missing) {
		my $queryString = join ",", @$arrayRef;
		#print STDERR "queryString = $queryString\n";
		my $xmlData = `efetch -db taxonomy -mode xml -id $queryString`;
	#	print STDERR "$xmlData\n";
		my $hash = $conv->xml2hash($xmlData);
		my $refType = ref($hash);
		#browse($hash);
		unless ($refType eq "HASH") {
			print STDERR "efetch failed\n";
			die;
		}
		my $refTypeTaxon = ref($hash->{'Taxon'});
		if ($refTypeTaxon eq "HASH") {
			my $temp = {};
			$temp->{'Taxon'} = [];
			push @{$temp->{'Taxon'}}, $hash->{'Taxon'};
			$hash = $temp;
		}
		foreach my $subHash (@{$hash->{'Taxon'}}) {
			unless (exists $subHash->{'AkaTaxIds'}->{'TaxId'} and exists $subHash->{'AkaTaxIds'}->{'TaxId'} and exists $subHash->{'ParentTaxId'} and exists $subHash->{'ScientificName'} and exists $subHash->{'Rank'}) {
				#print "subHash Missing info\n";	
				#browse($subHash);
				#print "\n";
				next;
			}
			my $self = $subHash->{'AkaTaxIds'}->{'TaxId'};
			my $parent = $subHash->{'ParentTaxId'};
			my $name = $subHash->{'ScientificName'};
			my $rank = $subHash->{'Rank'};
			unless (exists $missing->{$self}) {
				#print "AkaTaxIds not missing for $self\n";
				#browse($subHash);
				#print "\n";
				next;
			}
			unless (exists $data->{'parents'}->{$parent}) {
				#print "missing parent for $self, parent = $parent\n";
				$neededParents->{$parent} = 1;
			}
			$data->{'parents'}->{$self} = $parent;
			$data->{'names'}->{$self} = $name;
			$data->{'ranks'}->{$self} = $rank;
			$data->{'children'}->{$parent}->{$self} = 1;
		}
	}
}
my $encoded = $encoder->encode($data);
my $compressOut = compress($encoded);
print "$compressOut";















__END__
