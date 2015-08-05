#!/usr/bin/perl -w
use strict;
$|++;

my @in_files;

my %vars = (
	'name' => '',
	'scaffolds' => 0,
	'mito' => 0,
);

my %vars_details = (
	'name' => 'Name of organism (for output file)',
	'scaffolds' => 'Process scaffold assemblies (not recommended)',
	'mito' => 'Process mitochondrial chromosome (not recommended)',
);

process_cli();

# Globals
my $motif = "GATC";
my $motif_len = 4;

my $file = shift(@in_files);
my ($fhead) = $file =~ m/(.*?)./; 

$vars{'name'} ||= $fhead;

generate_track($file);

print STDOUT "All done.                \n\n";

sub generate_track {
	my $fn = shift;
	my $chr;
	
	# open output files
	my $track_name = $vars{'name'}.".GATC.gff";
	print STDOUT "Writing data to $track_name ...\n\n";
	open (TRACK, ">$track_name") || die "Cannot open $track_name for writing: $!\n";
	
	# load sequence data
	print STDOUT "Opening $fn ...\n";
	my $in="";
	
	if ($fn =~ m/\.gz$/) {
		# gzipped fasta file
		open (FN, "gunzip -c $fn |") || die "Error: cannot open $fn: $!\n\n";
	} else {
		open (FN, "<$fn") || die "Error: cannot open $fn: $!\n\n";
	}
	
	while (<FN>) {
		if (m/^\>/) {
			# New chromosome header
			if (m/scaffold/i) {
				unless ($vars{'scaffolds'}) {
					process($chr, $in);
					$in = "";
					$chr = "";
					next;
				}
			}
			
			if (m/mito/i) {
				unless ($vars{'mito'}) {
					process($chr, $in);
					$in = "";
					$chr = "";
					next;
				}
			}
			
			process($chr, $in);
			
			# get new chromosome name
			($chr) = m/^>(.*?)\s/;
			
			# Reset $in, and don't add the header line!
			$in = "";
			next;
		}
		chomp;
		s/\s//g;
		$in .= $_;
	}
	
	# process the last chromosome
	process($chr, $in);
	
	close FN;	
	close TRACK;
}

sub process {
	my ($chr, $in) = @_;
	return unless $chr && $in;

	# process what we've got ...
	print STDERR "Processing $chr ...                \n";
	motif_hash($in, $chr);
}

sub motif_hash {
	my ($in, $chr) = @_;
		
	my $l = length($in);
	for my $base (0 .. ($l-($motif_len+1))) {
		
		if ($base%100000==0) {
			my $pc = sprintf("%0.0f",($base*100)/$l);
			print STDERR "  $pc% done ...\r";
		}
		
		my $w = substr($in,$base,$motif_len);
		
		if ($w =~ m/$motif/gi) {
			print TRACK join("\t",$chr, ".", ".", $base, ($base+$motif_len), 1, '+', '.', '.'), "\n";
		} 
	}
}

sub process_cli {
	# CLI processing
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
			# if no parameter is specified we assume it's a switch ...
			# (could be a bit nicer and check this is ok with a hash representing data type ...)
			if (defined($vars{$1})) {
				$vars{$1} = 1;
			} else {
				print STDERR "Did not understand $_ ...\n";
				help();
			}
			next;
		}
		push @in_files, $_;
	}
	help() unless @in_files;
}

sub help {
	print STDERR <<EOT;

gatc.track.maker.pl -- generates a GFF file containing the locations of all GATC sites in the genome sequence
  (Use on a single genome FASTA file)

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