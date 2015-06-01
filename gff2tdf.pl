#!/usr/bin/perl -w
use strict;

unless (@ARGV) {
	help();
}

my %vars = (
	'genome' => 'dmel_r5.33',
	'igvpath' => '',
	'format' => 'gff',
);

my %vars_details = (
	'genome' => 'IGVTools genome file to use',
	'igvpath' => "path to IGVtools\n\r(leave blank if igvtools and igvtools.jar are in your path)",
	'format' => 'input format to use (gff or bedgraph)',
);

my @in_files;

process_cli();

my $igvtools = $vars{'igvpath'}."igvtools";
print STDERR "Creating TDF using $vars{'genome'} genome\n  (Use --genome=[genome] to change)\n\n";

convert();
exit 0;

sub process_cli {
	foreach (@ARGV) {
		if (/--(.*)=(.*)/) {
			unless (defined($vars{$1})) {
				print STDERR "Did not understand $_ ...\n";
				help();
			}
			my ($v, $opt) = ($1,$2);
			$vars{$v} = $opt;
			next;
		} elsif (/--h[elp]*/) {
			help();
		} elsif (/--(.*)/) {
			print STDERR "Please add a parameter to $_ ...\n\n";
			exit 1;
		}
		push @in_files, $_;
	}
}
	

sub convert {
	foreach my $fn (@in_files) {
		next if $fn =~ /^-/;
		my $in = $fn;
		print STDERR "Now working on $fn ...\n";
		
		my ($fhead) = $fn =~ m/(.*)\./;
		
		if ($vars{'format'} eq 'gff') {
			print STDERR "  converting to bedgraph ...\n";
			my @data;
			
			open (IN, "<$fn") || die "Cannot open $fn for reading: $!\n";
			open (OUT, ">$fhead.bedgraph") || die "Cannot open output wig file for writing: $!\n";
			
			print OUT qq(track type=bedGraph name="$fn" description="$fn"\n);
		
			while (<IN>) {
				chomp;
		
				my ($ref,$start,$end,$score) = (split)[0,3,4,5];
				next if $ref =~ m/^#/;
				next if $score =~ m/inf/i;
				next if $score =~ m/NA/i;
		
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
			$in = "$fhead.bedgraph";
		}
		
		print STDERR "  converting to TDF ...\n";
		
		`$igvtools toTDF $in $fhead.tdf $vars{'genome'}`;
		
		unlink("$fhead.bedgraph") if $vars{'format'} eq 'gff'; 
	}
}


sub help {
	print STDOUT <<EOT;

gff2tdf.pl -- converts gff or bedgraph files to IGV TDF format
	
Options:
EOT
	
	my $opt_len = 0;
	foreach (keys %vars) {
		my $l = length($_);
		#print "--> $_: $l\n";
		$opt_len = $l if $l > $opt_len;
	}
	
	$opt_len+=2;
	
	my $cols= `tput cols` || 80;
	
	my ($v, $val, $def, $def_format);
	my $help_format = "format STDOUT =\n"
		.' '.'^'.'<'x$opt_len . ' '. '^' . '<'x($cols-$opt_len-4) . "\n"
		.'$v, $def_format'."\n"
		.' '.'^'.'<'x$opt_len . '   '. '^' . '<'x($cols-$opt_len-6) . "~~\n"
		.'$v, $def_format'."\n"
		.".\n";
		
	eval $help_format;
	die $@ if $@;
	
	foreach my $k (sort (keys %vars)) {
		($v, $val, $def) = ($k, $vars{$k}, $vars_details{$k});
		$def||="";
		$def_format = $val ? "$def\n\r[Current value: $val]" : $def;
		$v = "--$v";
#		format =
# ^<<<<<<<<<<<<<<<<<<<< ^<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#$v, $def_format
# ^<<<<<<<<<<<<<<<<<<<<   ^<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ~~
#$v, $def_format
#.

		write();
		
	}
	print STDOUT "\n";
	exit 1;
}
