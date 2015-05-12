#!/usr/bin/perl -w

# Copyright � 2013-15, Owen Marshall

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or (at
# your option) any later version. 
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details. 
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 
# USA

use strict;
use File::Copy;

$|++;

my $version = "1.2pre2";
print STDERR "\nDamID-seq pipeline v$version\nCopyright � 2013-15, Owen Marshall\n\n";

# Global parameters
my %vars = (
	'bowtie' => 1,
	'bamfiles' => 0,
	'extend_reads' => 1,
	'len' => 300,
	'q' => 30,
	'bowtie2_path' => '',
	'samtools_path' => '',
	'bedtools_path' => '',
	'gatc_frag_file' => '',
	'kde_plot' => 0,
	'bowtie2_genome_dir' => '',
	'threads' => 7,
	'full_data_files' => 1,
	'qscore1min' => 0.4,
	'qscore1max' => 1.0,
	'qscore2max' => 0.9,
	'norm_override'=> 0,
	'output_format'=>'gff',
	'method_subtract' => 0,
	'pseudocounts' => 0,
	'ps_factor' => 10,
	'min_norm_value' => -5,
	'max_norm_value' => 5,
	'norm_steps' => 300,
	'just_align' => 0,
	'save_defaults' => 0,
	'load_defaults' => '',
	'reset_defaults' => 0,
	'bins' => '75',
);

my %vars_details = (
	'bowtie' => 'Perform bowtie2 alignment [1/0]',
	'bamfiles' => 'Only process BAM files',
	'extend_reads' => 'Perform read extension [1/0]',
	'len' => 'Length to extend reads to',
	'q' => 'Cutoff average Q score for aligned reads',
	'bowtie2_path' => 'path to bowtie2 executable (leave blank if in path)',
	'samtools_path' => 'path to samtools executable (leave blank if in path)',
	'bedtools_path' => 'path to BEDTools executable (leave blank if in path)',
	'gatc_frag_file' => 'GFF file containing all instances of the sequence GATC',
	'kde_plot' => 'create an Rplot of the kernel density fit for normalisation (requires R)',
	'bowtie2_genome_dir' => 'Directory and basename for bowtie2 .bt2 indices',
	'threads' => 'threads for bowtie2 to use',
	'full_data_files' => 'Output full binned ratio files (not only GATC array)',
	'qscore1min' => 'min decile for normalising from Dam array',
	'qscore1max' => 'max decile for normalising from Dam array',
	'qscore2max' => 'max decile for normalising from fusion-protein array',
	'norm_override'=> 'Normalise by this amount instead',
	'output_format' => 'Output tracks in this format [gff/bedgraph]',
	'method_subtract' => 'Output values are (Dam_fusion - Dam) instead of log2(Dam_fusion/Dam) (not recommended)',
	'pseudocounts' => 'Add this value of psuedocounts instead (default: optimal number of pseudocounts determined algorithmically)',
	'ps_factor' => 'Value of c in c*(reads/bins) formula for calculating pseudocounts (default = 10)',
	'min_norm_value' => 'Minimum log2 value to limit normalisation search at (default = -5)',
	'max_norm_value' => 'Maximum log2 value to limit normalisation search at (default = +5)',
	'norm_steps' => 'Number of points in normalisation routine (default = 300)',
	'just_align' => 'Just align the FASTQ files, generate BAM files, and exit',
	'save_defaults' => "Save runtime parameters as default\n\r(provide a name to differentiate different genomes -- these can be loaded with 'load_defaults')",
	'load_defaults' => "Load this saved set of defaults\n\r(use 'list' to list current saved options)",
	'reset_defaults' => 'Delete user-defined parameters',
	'bins' => 'Width of bins to use for mapping reads',
);

# Global variables
my @gatc_simple;
my @cli;
my @in_files;

my %ampnorm;
my %files;
my %index;
my %ar;
my %array;
my %bowtie_output;
my %gatc;
my %gatc_reverse_hash;
my %gatc_chr;
my %gatc_frags;
my %gatc_frag_score;
my %prot_hash;
my %norm_factors;
my %seg;
my %full_tracks;
my %counts;
my %bins;

my $HOME = (getpwuid($<))[7];
my $pi = 3.1415926536;
my $gatc_fragments;
my $denom;
my $damname;
my $frags;
my $no_paths;
my $windows_file;

# Read parameters if exist
read_defaults();
parameter_check();
#process_cli(0);

process_cli(1);


# CLI processing
check_paths();
save_defaults();

# Log file
init_log_file();

# Read input files
process_input_files();


#############################
### Pipeline workflow
###

align_sequences();
extend_reads();
load_gatc_frags() unless $vars{'just_align'};
calc_bins();

exit 0 if $vars{'just_align'};

for my $i (keys %files) {
	find_quants($i);
}

find_norm_factor();
normalize();
generate_ratio();

printout("All done.\n\n");

exit 0;


#############################
### Subroutines start here
###

sub check_paths {
	# Check paths
	if ($vars{'gatc_frag_file'} && !(-e "$vars{'gatc_frag_file'}")) {
		print STDERR "WARTNING: --gatc_frag_file option specified but file does not appear to exist.\n";
		if ($no_paths) {
			print STDERR "*** not saving file paths automatically (use --save_defaults=1 to override)\n\n";
			$no_paths=0;
		}
	}
	
	if ($vars{'bowtie2_genome_dir'} && !(-e "$vars{'bowtie2_genome_dir'}.1.bt2")) {
		print STDERR "WARTNING: --bowtie2_genome_dir option specified but files with basename does not appear to exist.\n\n(Please ensure you specify both the directory and the basename of the .bt2 index files -- i.e. if files are dmel_r5.57.1.bt2 dmel_r5.57.2.bt2 etc, located inside the directory Dm_r5.57, use '[path to]/Dm_r5.57/dmel_r5.57' as the option value ...\n\n";
		if ($no_paths) {
			print STDERR "*** not saving file paths automatically (use --save_defaults=1 to override)\n\n";
			$no_paths=0;
		}
	}
	
	# Parameter checks:
	for my $p ('gatc_frag_file','bowtie2_genome_dir') {
		unless ($vars{$p}) {
			die("Please use the --$p option to specifiy the $vars_details{$p} ...\n\n");
		}
	}
	
	# Save file paths if they haven't been set before
	if ($no_paths) {
		print STDOUT "All external files correctly located.  Paths will be saved as defaults and be used from now on.\n";
		$vars{'save_defaults'}=1;
	}
}


sub process_cli {
	# CLI processing
	my $order = shift;
	foreach (@ARGV) {
		if (/--(.*)=(.*)/) {
			unless (defined($vars{$1})) {
				print STDERR "Did not understand $_ ...\n";
				help();
			}
			my ($v, $opt) = ($1,$2);
			$opt =~ s/~/$HOME/;
			$vars{$v} = $opt;
			push @cli, "$v=$opt" if $order;
			next;
		} elsif (/--h[elp]*/) {
			help() if $order;
		} elsif (/--(.*)/) {
			# if no parameter is specified we assume it's a switch ...
			# (could be a bit nicer and check this is ok with a hash representing data type ...)
			if (defined($vars{$1})) {
				$vars{$1} = 1;
			} else {
				print STDERR "Did not understand $_ ...\n";
				help();
			}
			push @cli, "$1" if $order;
			next;
		}
		push @in_files, $_ if $order;
	}
}


sub process_input_files {
	
	# Index file
	my $index_file = "index.txt";
	
	# CLI files 
	# if no files specified, process all .gz files in directory.
	unless (@in_files) {
		printout("Searching for files ...\n");
		@in_files = glob("*.fastq.gz") || glob("*.gz") || glob("*.fastq") || glob("*.bam");
		if ($vars{'bamfiles'}) {
			@in_files = glob("*.bam");
		}
	}
	
	unless (@in_files) {
		printout("ERROR: no FASTQ or BAM files found (and no files specified on the command-line).\n(Use --help to see command-line options ...)\n\n")
	}
	
	printout("\n*** Reading data files ...\n");
	
	my $index_txt_error = <<EOT;
The pipeline script requires a single tab-delimited file named index.txt that lists sequencing adaptors with sample names -- for e.g.:

A6	Dam
A12	polII

Adaptors are taken from the .fastq file name (e.g. the file name for the Dam sample above should contain 'A006') and do not need to match the actual adaptors used.
EOT

	if ($vars{'bamfiles'}) {
		# Bamfiles
		foreach my $l (@in_files) {
			my ($name) = $l =~ m/(.*?)[-|\.].*bam/i;
			print "$name\t$l\n";
			$files{$name}[0]=$l;
			if ($name =~ m/^dam/i) {
				die("\nERROR: more than one Dam sample detected.  Please only use one Dam control per run -- specify the files to process on the command-line\n\n") if $damname;
				$damname = $name;
			}
		}
	} elsif (@in_files && !(-e $index_file)) {
		foreach my $l (@in_files) {
			my ($name) = $l =~ m/(.*?)(_|-ext|\.)/i;
			my ($ext) = $l =~ m/.*\.(.*)/;
			$vars{'bamfiles'} = 1 if $ext =~ /bam/;
			
			print "$name\t$l\n";
			$files{$name}[0]=$l;
			if ($name =~ m/^dam/i) {
				die("\nERROR: more than one Dam sample detected.  Please only use one Dam control per run -- specify the files to process on the command-line\n\n") if $damname;
				$damname = $name;
			}
		}
	} elsif (-e $index_file) {
		open (NORM, "<$index_file") || die "Unable to open index file for reading: $!\n";
		my @norm = <NORM>;
		close NORM;
		chomp(@norm);
		printout("\nIndex\tName\n");
		foreach my $l (@norm) {
			next if $l =~ m/^$/; # skip new lines
			next if $l =~ m/^#/; # skip comments
			
			my ($i, $name) = split(/\s+/,$l);
			
			unless (($name) && ($i =~ m/.*?(?:0*)+\d+/i)) {
				die("\nERROR: misformatted index.txt file.\n\n$index_txt_error\n")
			}
			
			printout("$i\t$name\n");
			my ($i_num) = ($i =~ m/.*?(?:0*)+(\d+)/i);
			$index{$i_num}=$name;
		}
		
		printout("\n*** Matching adaptors ...\n");
		foreach my $l (@in_files) {
			print "$l\n";
			
			## Change this next line's regexp to match your sequencing format (currently matches, e.g. "Index6" or "A006")
			my ($i) = $l =~ m/.*?(?:index|a)(?:0*)(\d+)/i;
			
			if ($index{$i}) {
				printout("$index{$i}: Index $i\n");
				$files{$index{$i}}[0]=$l;
				
				if ($index{$i} =~ m/^dam/i) {
					die("\nERROR: more than one Dam sample detected.  Please only use one Dam control per run -- specify the files to process on the command-line\n\n") if $damname;
					$damname = $index{$i}
				}
			} else {
				die ("\nERROR: cannot find adaptor reference in fastq files.\n\n$index_txt_error\n");
			}
		}
	} else {
		help();
	}
	
	# Check that there's a Dam control sample ...
	unless ($vars{'just_align'}) {
		die ("Error: no Dam control sample detected!\n\n") unless $damname;
	}
	
	# If we're only processing .bam files, we won't need to do an alignment or extend reads...
	if ($vars{'bamfiles'}) {
		$vars{'bowtie'} = 0;
		$vars{'extend_reads'} = 0;
	}
}


sub read_defaults {
	# Create config directory if it doesn't exist
	unless (-d "$HOME/.config") {
		mkdir("$HOME/.config");
	}
	
	unless (-d "$HOME/.config/damid_pipeline") {
		mkdir("$HOME/.config/damid_pipeline");
	}
	
	# Migration for version 1.0 
	if (-e "$HOME/.config/damid_pipeline_defaults") {
		move("$HOME/.config/damid_pipeline_defaults","$HOME/.config/damid_pipeline/defaults")
	}
	
	if ($vars{'load_defaults'}) {
		if ($vars{'load_defaults'} eq 'list') {
			print_defaults_files();
		} else {
			unless (-e "$HOME/.config/damid_pipeline/defaults.$vars{'load_defaults'}") {
				print STDERR "Error: cannot find defaults file $vars{'load_defaults'}\n\n";
				print_defaults_files();
			}
			# load the defaults
			read_defaults_file("$HOME/.config/damid_pipeline/defaults.$vars{'load_defaults'}");
		}
	} elsif (-e "$HOME/.config/damid_pipeline/defaults") {
		# read default parameters if exist
		read_defaults_file("$HOME/.config/damid_pipeline/defaults");
	}
}

sub print_defaults_files {
	my @defaults_files = glob("$HOME/.config/damid_pipeline/defaults.*");
	print STDOUT "Available defaults:\n";
	if (@defaults_files) {
		foreach (@defaults_files) {
			my ($name) = m#$HOME/\.config/damid_pipeline/defaults\.(.*)$#;
			print STDOUT "  $name\n";
		}
	} else {
		print STDOUT "  No genome-specific defaults files found\n  (Use --save_defaults=[name] to save them)\n";
	}
	print STDOUT "\n";
	exit 0;
}

sub read_defaults_file {
	my $file = shift;
	open(DEFAULTS, "<$file") || die("Cannot open defaults file for writing: $!\n\n");
	while (<DEFAULTS>) {
		chomp;
		my ($p,$v) = split;
		$vars{$p} = $v;
	}
	close DEFAULTS;
}


sub save_defaults {
	# Save parameters if requested
	if ($vars{'save_defaults'}) {
		my $name;
		if ($vars{'save_defaults'} eq '1') {
			$name = '';
		} else {
			$name = ".$vars{'save_defaults'}";
		}
		
		print STDERR "Writing defaults to file ...\n";
		open(DEFAULTS,">$HOME/.config/damid_pipeline/defaults$name") || die("Cannot open defaults file for writing: $!\n\n");
		for my $p (keys %vars) {
			next if $p eq 'save_defaults';
			next if $p eq 'reset_defaults';
			next unless $vars{$p};
			print DEFAULTS "$p\t$vars{$p}\n";
		}
		print STDERR "Done.\n\n";
		close DEFAULTS;
		exit 0;
	}
	
	# reset defaults
	if ($vars{'reset_defaults'}) {
		my $name;
		if ($vars{'reset_defaults'} eq '1') {
			$name = '';
		} else {
			$name = ".$vars{'reset_defaults'}";
		}
		
		unlink("$HOME/.config/damid_pipeline/defaults$name") if -e "$HOME/.config/damid_pipeline/defaults$name";
		print STDERR "Defaults reset ... please restart.\n\n";
		exit 0;
	}
}

sub parameter_check {
	# Parameter check -- first run:
	for my $p ('gatc_frag_file','bowtie2_genome_dir') {
		unless ($vars{$p}) {
			$no_paths = 1;
		}
	}
}

sub align_sequences {
	if ($vars{'bowtie'}) {
		
		printout("\n\n*** Aligning files with bowtie2 ...\n");
		foreach my $fn (keys %files) {
			my $pair1 = $files{$fn}[0];
			
			printout("\nNow working on $fn ...\n");
			
			if ($vars{'bowtie'}) {
				$bowtie_output{$fn} = `$vars{'bowtie2_path'}bowtie2 -p $vars{'threads'} -x $vars{'bowtie2_genome_dir'} -U $pair1 -S $fn.sam 2>&1`;
				printout("$bowtie_output{$fn}\n");
			}
		}
	}
}

sub extend_reads {
	# Extend reads
	if ($vars{'extend_reads'}) {
		printout("\n*** Extending reads to $vars{'len'} bases ...\n\n");
		foreach my $fn (keys %files) {
			printout("Reading input file: $fn ...\n");
			open (IN, "<$fn.sam") || die "Unable to read $fn: $!\n";
			
			printout("Processing data ...\n");
			open (OUT, ">$fn-ext$vars{'len'}.sam");
			
			my $c=0;
			my $seqs;
			while (<IN>) {
				chomp;
				my ($qname, $flag, $chr, $pos, $mapq, $cigar, $rnext, $pnext, $tlen, $seq, $qual) = split('\s+');
				
				$c++;
				unless ($c%100000) {
					print "$c lines processed ...\r";
				}
				
				unless ($seq) {
					print OUT "$_\n";
					next;
				}
				
				next unless $mapq>=$vars{'q'};
				
				$seqs++;
				
				$cigar="$vars{'len'}M";
				$seq = "*"; # We're extending reads and we have no sequence information for the extension ...
				$qual= "*"; #  ... ditto for quality.  Thankfully the SAMfile spec allows for no sequence or quality information.
				
				print OUT join("\t", $qname, $flag, $chr, $pos, $mapq, $cigar, $rnext, $pnext, $tlen, $seq, $qual), "\n";
				
			}
			close IN;
			close OUT;
			
			# the original SAM files are huge, so we nuke them
			unlink("$fn.sam");
			
			printout("Seqs extended (>q30) = $seqs\n\n");
		}
	}
}

sub calc_bins {
	printout("\n*** Calculating bins ...");
	foreach my $n (keys %files) {
		my $fn = $vars{'bamfiles'} ? $files{$n}[0] : "$n-ext$vars{'len'}";
		$fn =~ s/\.bam$//;
		
		printout("\nNow working on $fn ...\n");
		
		unless ($vars{'bamfiles'}) {
			printout("Generating .bam file ...\n");
			`$vars{'samtools_path'}samtools view -Sb $fn.sam > $fn.bam`;
			
			unlink("$fn.sam");
		}
		
		# Obtain readcounts via samtools ...
		my $seqs = `$vars{'samtools_path'}samtools view -c -F 4 $fn.bam`;
		chomp($seqs);
		
		$counts{$n}=$seqs;
		printout("  $seqs reads\n");
		if ($n =~ m/^dam$/i) {
			$denom = $counts{$n};
		}
		
		# short-circuit for just_align option
		next if $vars{'just_align'};
				
		printout("  Generating bins from $fn.bam ...\n");
		
		
		my %cov;
		my %gff;
		my %bases;
		my $lines;
		# We now use manually calculate coverage, rather than using bedtools' -coverage option, for both memory considerations and ease of use.
		# (bedtools' -coverage memory consumption is pretty crazy, getting up to 8Gb for ~25M reads ... this implementation uses about 1/10th of that memory.)
		# Speed considerations are similar (this is seems to be slightly faster ...)
		open(COV, "samtools view -h $fn.bam |") || die ("Error: unable to open bamfile $!\n");
		while (<COV>) {
			chomp();
			$lines++;
			my ($ref, $bit, $chr, $pos, $q, $len) = split(/\t/);
			unless ($len) {
				if (m/\@SQ\s+SN:([\d|\w]+)\s+LN:(\d+)/) {
					# Grab chromosome sizes from BAM file
					my ($chr, $size) = ($1, $2);
					$bases{$chr} = $size;
				} 
				next;
			}
			
			$len =~ s/[^\d]//g;
			
			if ($lines % 20000 == 0) {
				my $pc = sprintf("%0.1f", ($lines*100 / $counts{$n}));
				print STDERR "  $pc% processed ...\r"
			}
						
			# break into bins
			my $start =  int($pos/$vars{'bins'});
			my $end =  int(($pos + $len)/$vars{'bins'});
						
			foreach my $bin ($start .. $end) {
				$cov{$chr}{$bin*$vars{'bins'}}++;
			}	
		}
		close COV;
		
		my $read_bins=0;
		foreach my $chr (keys %cov) {
			for (my $bin = 0; $bin < $bases{$chr}+$vars{'bins'}; $bin += $vars{'bins'}) {
				# need to fill every bin ...
				$read_bins++;
				my $end = $bin+($vars{'bins'}-1);
				my $score = $cov{$chr}{$bin} || 0;
				
				push @{$gff{$chr}}, [$bin, $end, $score];
				push @{$full_tracks{$n}}, [$chr, $bin, $end, $score] if $vars{'full_data_files'};
			}
		}
		
		printout("  Sorting ...                 \n");
		my %tmp_gff;
		foreach my $k (keys %gff) {
			@{$tmp_gff{$k}} = sort { $a->[0] <=> $b->[0] } @{$gff{$k}};
		}
		%gff = %tmp_gff;
	
		
		$bins{$n} = $read_bins;
		
		my $gatc_nonhits = 0;
					
		my $index = 0;
		
		foreach my $chr (keys %gff) {
			unless ($gatc{$chr}) {
				print STDERR "  Warning: GFF file contains chromosome identity ($chr) not found in GATC file.\n";
				next;
			}
			
			my $last_input=0;
			foreach my $i (0 .. @{$gatc{$chr}}-2) {
				$index++;
				# Get two adjacent GATC sites
				my ($mida) = @{$gatc{$chr}}[$i];
				my ($midb) = @{$gatc{$chr}}[$i+1];
		
				if ($index%100 == 0) {
					my $pc = sprintf("%0.0f",($index*100)/$#gatc_simple);
					print STDOUT "  $pc\% processed ...\r";
				}
				
				my $sum;
				my @data;
				my $count;
				# Run through the probe array to find the probes that lie within the fragment
				foreach my $l ($last_input .. @{$gff{$chr}}-2) {
					# Takes way too long to scan through the entire array for each GATC fragment, and is in theory unnecessary
					# This proceedure uses $last_input to store the last probe found within an array, and start the next search from there
					# (relies on an ordered array of data, hence the sorting above)
					
					my ($start, $end, $score) = @{$gff{$chr}[$l]};
					my ($startn, $endn, $scoren) = @{$gff{$chr}[$l+1]};
					
					$last_input= ($l-1 > 0 ? $l-1 : 0);
					
					if ($start > $midb) {
						last;
					}
					next if $end < $mida;
					
					if (($end > $mida) && ($start < $midb)) {
						push @data, $score;
					}
				}
				
				
				unless (@data) {
					$gatc_nonhits++;
					next;
				}
				my ($mean, $sd, $se) = stdev(@data);
				
				push(@{$gatc_frags{"$chr-$mida-$midb"}}, $n);
				$gatc_frag_score{"$chr-$mida-$midb"}{$n}=$mean;
				push @{$array{$n}}, [ $chr, $mida, $midb, $mean];
			}
		} 
		printout("\n");
	}
}

sub normalize {
	printout("\n*** Normalising counts and adding background ...\n");
	foreach my $n (keys %files) {
		printout("\nProcessing input: $n ...\n");
		
		my $norm;
		unless ($n =~ m/^dam$/i) {
			if ($vars{'norm_override'}) {
				printout("Normalisation override!\n  Would have normalised by $norm_factors{$n}\n  ... will instead normalise by $vars{'norm_override'}\n");
				$norm = $vars{'norm_override'};
			} elsif ($vars{'rpm_norm'}) {
				# RPM normalisation goes here ...
				
			} else {
				printout("  ... normalising by $norm_factors{$n}\n");
				$norm = $norm_factors{$n};
			}
		}
		
		# move pseudocounts into ratio generation routine ...
		foreach (@{$array{$n}}) {
			@{$_}[3] *= $norm unless $n =~ m/^dam$/i;
			# @{$_}[3] += $pseudocounts;
		}
		
		if ($vars{'full_data_files'}) {
			foreach (@{$full_tracks{$n}}) {
				@{$_}[3] *= $norm unless $n =~ m/^dam$/i;
				# @{$_}[3] += $pseudocounts;
			}
		}
	}
}

sub printout {
	my $s = shift;
	print STDERR $s;
	print STATS $s;
}


sub generate_ratio {
	printout("\n\n*** Generating ratios ...\n");
	foreach my $n (keys %files) {
		next if $n =~ m/^dam$/i;
		
		printout("\nNow working on $n ...\n"); 

		write_file($n, 1);
		write_file($n, 0) if $vars{'full_data_files'};
	}
}


sub write_file {
	my ($n, $write_gatcs) = @_;

	my %data;
	
	# find the lowest number of compared counts -- this will now use a different number of pseudocounts per sample
	my $psc_min = ($counts{$n}, $denom)[$counts{$n} > $denom];
	my $pseudocounts = ($vars{'pseudocounts'} ? $vars{'pseudocounts'} : $vars{'ps_factor'}*$psc_min/$bins{$n}); # pseudocounts value is related to total reads/number of bins
		
	my $dref = ($write_gatcs ? \%array : \%full_tracks);
	
	printout("  Reading Dam ...\n");
	foreach (@{ $dref -> {$damname}}) {
		my ($chr, $start, $end, $score) = @{$_};
		push @{ $data{$chr}{$start}}, $end;
		push @{ $data{$chr}{$start}}, $score;
	}

	printout("  Reading $n ...\n");
	foreach (@{ $dref-> {$n}}) {
		my ($chr, $start, $end, $score) = @{$_};
		push @{ $data{$chr}{$start}}, $score;
	}
	
	my $track_type = $write_gatcs ? 'gatc.' : '';
	my $name_spacer = $vars{'method_subtract'} ? ".sub." : ".";
	my $fname = "$n-vs-Dam$name_spacer$track_type";
	my $fout=$fname.($vars{'output_format'} eq 'bedgraph' ? 'bedgraph' : 'gff');
	
	my $print_psc = sprintf("%0.2f",$pseudocounts);
	printout("  ... adding $print_psc pseudocounts to each sample\n");
	
	open (OUT, ">$fout") || die "Unable to open $fout for writing: $!\n";
	print OUT qq(track type=bedGraph name="$fname" description="$fname"\n) if $vars{'output_format'} eq 'bedgraph';
	
	foreach my $chr (sort keys %data) {
		foreach my $start (sort {$a <=> $b} keys %{ $data{$chr}}) {
			my ($end, $score1, $score2) = @{ $data{$chr}{$start}};
			
			my $score = generate_score($score1,$score2, $pseudocounts);
			
			if ($vars{'output_format'} eq 'bedgraph') {
				print OUT join("\t", $chr, $start, $end, $score), "\n";
			} else {
				print OUT join("\t", $chr, '.', '.', $start, $end, $score, '.', '.', '.'), "\n";
			}
		}
	}
	close OUT;
}


sub generate_score {
	my ($score1, $score2, $pseudocounts) = @_;
	my $score;
	
	if ($vars{'method_subtract'}) {
		$score = $score2 - $score1;
	} else {
		
		$score1 += $pseudocounts;
		$score2 += $pseudocounts;
		
		unless (($score2) && ($score1)) {
			$score = 0;
		} else {
			$score = log($score2 /$score1)/log(2);
		}
	}
	
	return $score;
}

sub find_norm_factor {
	printout("\n\n*** Calculating normalisation factor ...\n");
		
	foreach my $n (keys %files) {
		printout("Now working on $n ...\n");
		
		next if $n =~ m/^dam$/i;
		my @ratios;
		my $total_frags;
		foreach my $key (keys %gatc_frag_score) {
			$total_frags++;
			my %score;
			my %qscore;
			$score{1} = $gatc_frag_score{$key}{$damname};
			$score{2} = $gatc_frag_score{$key}{$n};
			
			next unless ($score{1} && $score{2});
			
			$qscore{1} = qsc($score{1},$damname);
			$qscore{2} = qsc($score{2}, $n);
			
			if (($qscore{1}>=$vars{'qscore1min'}) && ($qscore{1}<=$vars{'qscore1max'}) && ($qscore{2}<=$vars{'qscore2max'})) {
				if ($vars{'old_norm_method'}) {
					# old approximation method provided for compatibility -- not recommended
					push @ratios, $score{2}/$score{1};
				} else {
					push @ratios, log($score{2}/$score{1})/log(2);
				}
			}
			
		}
		my $total_measurements = @ratios;
		
		my $avg;
		if ($vars{'old_norm_method'}) {
			$avg = stdev(@ratios);
		} else {
			# new normalisation method finds the maximum of the gaussian kernel density estimate of log2 ratio data
			$avg = 2**kdenmax(@ratios);
		}
		
		# Just in case ...
		if ($avg == 0) {
			$avg = 1;
			printout("\n*** WARNING: average ratio was zero ... no normalisation applied.\n(Unless you are processing the provided example files or *very* small read numbers, this shouldn't happen ...)\n\n")
		}
		
		my $norm = 1/$avg;
		printout("  Norm factor = $norm based off $total_measurements frags (total $total_frags)\n");
		$norm_factors{$n}=$norm;
	}
}

sub qsc {
	my ($score,$prot) = @_;
	my $qscore;
	foreach my $c (sort {$a <=> $b} keys %seg) {
		if ($score < $seg{$c}{$prot}) {
			$qscore = $c ;
			last;
		}
	}
	$qscore||=1;
	return $qscore;
}


sub find_quants {
	my $prot = shift;
	my @frags;
	
	my $total_coverage;
	
	printout("\nNow working on $prot ...\n");

	foreach my $k (keys %gatc_frag_score) {
		my $score = $gatc_frag_score{$k}{$prot} || 0;
		
		#sprint "$k $prot $score\n";
		
		next unless $score > 0; # use only non-zero elements of the array
		push (@frags, $score);
	}
	
	printout("Sorting frags ...\n");
	my @sorted_frags = sort {$a <=> $b} @frags;
	
	my @quants;
	for (my $q=0.1;$q<=1;$q+=0.1) {
		push @quants, [$q, int($q * @sorted_frags)];
	}
	
	printout("Finding quants ...\n");
	foreach  (@quants) {
		my $cut_off = @{$_}[0];
		my $score = $sorted_frags[@{$_}[1]];
		printout("   Quant $cut_off:".sprintf("%0.2f",$score)."\n");
		$seg{$cut_off}{$prot} = $score;
	}
}

sub stdev {
	my @a = @_;
	my $n = @a;
	return ($a[0], 0, 0, $a[0]) if $n<2;
	
	my $total;
	foreach (@a) {
		$total+=$_;
	}
	my $mean = $total/$n;
	
	my $sum;
	foreach my $x (@a) {
		$sum+=abs($x-$mean);
	}
	my $sd= $sum/($n-1);
	my $se = $sd/sqrt($n);
	
	my $half = int($n/2);
	my @b = sort {$a <=> $b} @a;
	my $median = $b[$half];
	
	return ($mean, $sd, $se, $median);
}

sub kden {
	# gaussian kernel density estimator
	my @x = @_;
	
	my $s = 0;
	my $n = $#x;
	
	# KDE is measured over an equidistant grid:
	my $xmin = max((min(@x))[1],$vars{'min_norm_value'});
	my $xmax = min((max(@x))[1],$vars{'max_norm_value'});
	my $steps = $vars{'norm_steps'};
		
	## bandwidth via Silverman's estimator (eq 3.31 in Silverman 1986)
	my ($mean, $sd) = stdev(@x);
	my $h = min($sd,iqr(@x))/1.34 * $n**(-1/5);
	
	## kernel density
	my @sample;
	foreach my $st (0 .. $steps) {
		my $y = $xmin + $st*($xmax-$xmin)/$steps;
		my $pc = sprintf("%0.0f",100*$st/$steps);
		print STDERR "  Fitting kernel density $pc% complete ...\r";
		
		my $s = 0;
		for my $i (0..$n) {
			$s += exp(-1/2*(($x[$i]-$y)/$h)**2)/(sqrt(2*$pi));
		}
		my $fh = $s/($n*$h);
		
		push @sample, [$y, $fh];
	}
	
	return (@sample);
}

sub iqr {
	my @x = @_;
	@x = sort {$a <=> $b} @x;
	
	my $qmin_in = int($#x/4);
	my $qmax_in = int(3*$#x/4);
	
	my $qmin = $x[$qmin_in];
	my $qmax = $x[$qmax_in];
	
	my $iqr = $qmax-$qmin;
	return $iqr;
}

sub kdenmax {
	my @x = @_;
	
	my @kd = kden(@x);
	
	my @col1 = map $_ -> [0], @kd;
	my @col2 = map $_ -> [1], @kd;
	
	my $col1 = join(",",@col1);
	my $col2 = join(",",@col2);
	
	my ($index, $peak_max) = max(@col2);
	
	if ($vars{'kde_plot'}) {
		# plot the kernel density fit if the user desires ...
		open(DAT, ">dat") || die ("Cannot open file for writing: $!\n");
		foreach (@x) {
			print DAT "$_\n";
		}
		close DAT;
		
		`r -e 'dev.new(file="kden.pdf"); d <- read.table("dat",header=F); hist(d\$V1,breaks=500,prob=T,main="log2 ratio of data used for normalisation"); points(c($col1),c($col2),xlab="x",ylab="Density",col="red");  lines(density(d\$V1),col="green"); abline(v=$col1[$index],col="blue"); dev.off()'`;
		unlink("dat");
	}
	
	return $col1[$index];
}

sub max {
    my ($max, @vars) = @_;
	my $index=0;
	$max||=0;
    for my $i (0..$#vars) {
        ($max, $index) = ($vars[$i], $i+1) if $vars[$i] > $max;
    }
    return ($index, $max);
}

sub min {
    my ($min, @vars) = @_;
	my $index=0;
	$min||=0;
    for my $i (0..$#vars) {
        ($min, $index) = ($vars[$i],$i+1) if $vars[$i] < $min;
    }
    return ($index, $min);
}

sub help {
	print STDERR "Options:\n";
	
	my $opt_len = 0;
	foreach (keys %vars) {
		my $l = length($_);
		#print "--> $_: $l\n";
		$opt_len = $l if $l > $opt_len;
	}
	
	$opt_len+=2;
	
	my $cols= `tput cols` || 80;
	
	my ($v, $val, $def, $def_format, $current);
	my $curr_len = 40;
	my $help_format = "format STDOUT =\n"
		.' '.'^'.'<'x$opt_len . ' '. '^' . '<'x($cols-$opt_len-4) . "\n"
		.'$v, $def_format'."\n"
		.' '.'^'.'<'x$opt_len . '   '. '^' . '<'x($cols-$opt_len-6) . "~~\n"
		.'$v, $def_format'."\n"
		.".\n";
		
	#my $help_format = "format STDOUT =\n"
	#.' '.'^'.'<'x$opt_len . ' '. '^'.'<'x($cols-$opt_len-4-2-$curr_len) . '^'.'<'x($curr_len-2)."\n"
	#.'$v, $def_format, $val'."\n"
	#.' '.'^'.'<'x$opt_len . '   '. '^'.'<'x($cols-$opt_len-6-2-$curr_len) . '^'.'<'x($curr_len-2). "~~\n"
	#.'$v, $def_format, $val'."\n"
	#.".\n";
		
		
	eval $help_format;
	die $@ if $@;
	
	foreach my $k (sort (keys %vars)) {
		($v, $val, $def) = ($k, $vars{$k}, $vars_details{$k});
		$def||="";
		$def_format = $val ? "$def\n\r[Current value: $val]" : $def;
		$v = "--$v";
		write();		
	}
	print STDOUT "\n";
	exit 1;
}

sub load_gatc_frags {
	printout("\n*** Reading GATC file ...\n");
	
	if ($vars{'gatc_frag_file'} =~ m/\.gz$/) {
		# gzipped file
		open (GATC, "gunzip -c $vars{'gatc_frag_file'} |") || die "Error: cannot open GATC file $vars{'gatc_frag_file'}: $!\n\n";
	} else {
		open (GATC, "<$vars{'gatc_frag_file'}") || die "Error: cannot open GATC file $vars{'gatc_frag_file'}: $!\n\n";
	}
	
	#open (GATC, "<","$vars{'gatc_frag_file'}") || die "Unable to read GATC file: $!\n";
	
	foreach (<GATC>) {
		my ($chr, $source, $type, $start, $end, $score, $b, $c, $name) = split('\t');
		my $mid = ($start+$end)/2;
		push (@gatc_simple, [$chr, $mid]);
		push @{$gatc{$chr}}, $mid;
	}
	close GATC;
	
	print STDERR "Sorting ...\n";
	my %tmp;
	foreach my $k (keys %gatc) {
		@{$tmp{$k}} = sort { $a <=> $b } @{$gatc{$k}};
	}
	%gatc = %tmp;
}

sub init_log_file {
	my $date = get_date();
	open (STATS, ">pipeline-$date.log") || die "Could not open bowtie output file for writing: $!\n";
	my $args = join("|",@cli);
	
	printout("Version $version\n\n");
	print STATS "Command-line options: @ARGV\n\n";
}

sub get_date {
	my $date = localtime();
	$date =~ s/:\d\d\s/ /;
	$date =~ s/\s+/_/g;
	$date =~ s/:/-/g;
	return($date);
}