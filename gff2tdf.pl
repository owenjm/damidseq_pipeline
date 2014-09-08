#!/usr/bin/perl -w
use strict;

die "gff2tdf.pl [filenames to process]" unless @ARGV;

# Path to IGVtools
# leave blank if igvtools executable (and igvtools.jar) are in your path
my $IGV_path = "";
my $igvtools = $IGV_path."igvtools";

foreach my $fn (@ARGV) {
	print STDERR "Now working on $fn ...\n";
	print STDERR "  converting to wig ...\n";
	
	my ($fhead) = $fn =~ m/(.*).gff/;
	
	my @data;
	
	open (IN, "<$fn") || die "Cannot open $fn for reading: $!\n";
	open (OUT, ">$fhead.wig") || die "Cannot open output wig file for writing: $!\n";
	
	print OUT qq(track type=bedGraph name="$fn" description="$fn"\n);

	while (<IN>) {
		chomp;

		my ($ref,$start,$end,$score) = (split)[0,3,4,5];
		next if $ref =~ m/^#/;

		$ref =~ s/CHROMOSOME_|chr//;
		$ref = "chr".$ref;

		$start--; # zero-based, half-open
		push @data, [$ref,$start,$end,$score];
	}
	close IN;
	
	@data = sort {$a->[1] <=> $b->[1]} @data;
	@data = sort {$a->[0] cmp $b->[0]} @data;

	foreach my $l (@data) {
		my ($ref,$start,$end,$score) = @{$l};
		print OUT join("\t", $ref,$start,$end,$score), "\n";
	}
	
	close OUT;
	
	print STDERR "  converting to TDF ...\n";
	
	`$igvtools toTDF $fhead.wig $fhead.tdf dmel_r5.33`;
	
	unlink("$fhead.wig")
}