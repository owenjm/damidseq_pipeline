#!/usr/bin/perl -w

# Copyright © 2013-14, Owen Marshall

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
$|++;

my $version = "1.0";
print STDERR "\nBrand lab DamID-seq pipeline v$version\nCopyright © 2013-14, Owen Marshall\n\n";

# Global parameters
my %vars = (
	'bowtie' => 1,
	'extend_reads' => 1,
	'len' => 300,
	'q' => 30,
	'bins' => 75,
	'bowtie2_path' => '',
	'samtools_path' => '',
	'coveragebed_path' => '',
	'window_file_dir' => '',
	'gatc_frag_file' => '',
	'bowtie2_genome_dir' => '',
	'threads' => 7,
	'full_data_files' => 1,
	'qscore1min' => 0.6,
	'qscore1max' => 0.9,
	'qscore2max' => 0.9,
	'norm_override'=> 0,
	'pseudocounts' => 0,
	'save_defaults' => 0,
	'reset_defaults' => 0,
);

my %vars_details = (
	'bowtie' => 'Perform bowtie2 alignment [1/0]',
	'extend_reads' => 'Perform read extension [1/0]',
	'len' => 'Length to extend reads to',
	'q' => 'Cutoff Q score for aligned reads',
	'bins' => 'Genome bin size for coverageBed (requires premade file)',
	'bowtie2_path' => 'path to bowtie2 executable (leave blank if in path)',
	'samtools_path' => 'path to samtools executable (leave blank if in path)',
	'coveragebed_path' => 'path to coverageBed executable (leave blank if in path)',
	'window_file_dir' => 'directory for window files for coverageBed',
	'gatc_frag_file' => 'GFF file containing all instances of the sequence GATC',
	'bowtie2_genome_dir' => 'Directory and basename for bowtie2 .bt2 indices',
	'threads' => 'threads for bowtie2 to use',
	'full_data_files' => 'Output full binned ratio files (not only GATC array)',
	'qscore1min' => 'min decile for normalising from Dam array',
	'qscore1max' => 'max decile for normalising from Dam array',
	'qscore2max' => 'max decile for normalising from fusion-protein array',
	'norm_override'=> 'Normalise by this amount instead',
	'pseudocounts' => 'Add this value of psuedocounts instead (default: optimal number of pseudocounts determined algorithmically)',
	'save_defaults' => 'Save runtime parameters as default',
	'reset_defaults' => 'Delete user-defined parameters'
);

# Home directory
my $HOME = (getpwuid($<))[7];

# Global variables
my %ampnorm;
my %files;
my %index;
my %ar;
my %array;
my %bowtie_output;
my %gatc_reverse_hash;
my %gatc_chr;
my @gatc;
my @gatc_simple;
my %gatc_frags;
my %gatc_frag_size;
my %gatc_frag_score;
my $gatc_fragments;
my %prot_hash;
my %norm_factors;
my %seg;
my %full_tracks;
my %counts;
my $denom;
my $damname;

# Read parameters if exist
if (-e "$HOME/.config/damid_pipeline_defaults") {
	open(DEFAULTS, "<$HOME/.config/damid_pipeline_defaults") || die("Cannot open defaults file for writing: $!\n\n");
	while (<DEFAULTS>) {
		chomp;
		my ($p,$v) = split;
		$vars{$p} = $v;
	}
	close DEFAULTS;
}

# CLI processing
my @cli;
my @in_files;
foreach (@ARGV) {
	if (/--(.*)=(.*)/) {
		unless (defined($vars{$1})) {
			print STDERR "Did not understand $_ ...\n";
			help();
		}
		my ($v, $opt) = ($1,$2);
		$opt =~ s/~/$HOME/;
		$vars{$v} = $opt;
		push @cli, "$v=$opt";
		next;
	} elsif (/--h[elp]*/) {
		help();
	} elsif (/--(.*)/) {
		print STDERR "Please add a parameter to $_ ...\n\n";
		exit 1;
	}
	push @in_files, $_;
}

# Save parameters if requested
if ($vars{'save_defaults'}) {
	unless (-d "$HOME/.config/") {
		mkdir("$HOME/.config/") || die("Cannot create $HOME/.config directory: $!\n\n")
	}
	
	print STDERR "Writing defaults to file ...\n";
	open(DEFAULTS,">$HOME/.config/damid_pipeline_defaults") || die("Cannot open defaults file for writing: $!\n\n");
	for my $p (keys %vars) {
		next if $p eq 'save_defaults';
		next if $p eq 'reset_defaults';
		next unless $vars{$p};
		print DEFAULTS "$p\t$vars{$p}\n";
	}
	print STDERR "Done.\n\n";
	close DEFAULTS;
}

# reset defaults
if ($vars{'reset_defaults'}) {
	unlink("$HOME/.config/damid_pipeline_defaults") if -e "$HOME/.config/damid_pipeline_defaults";
	print STDERR "Defaults reset ... please restart.\n\n";
	exit 0;
}

# Parameter checks:
for my $p ('window_file_dir','gatc_frag_file','bowtie2_genome_dir') {
	unless ($vars{$p}) {
		die("Please use the --$p option to specifiy the $vars_details{$p} ...\n\n")
	}
}

# Global input files
my $windows_file = "$vars{'window_file_dir'}/dm3_windows_$vars{'bins'}";
my $bowtie_genome_file = $vars{'bowtie2_genome_dir'};

# Log file
my $date = localtime();
$date =~ s/:\d\d\s/ /;
$date =~ s/\s+/_/g;
$date =~ s/:/-/g;

open (STATS, ">pipeline-$date.log") || die "Could not open bowtie output file for writing: $!\n";
my $args = join("|",@cli);
printout("Version $version\n\n");
print STATS "Command-line options: @ARGV\n\n";

# Index file
my $index_file = "index.txt";

# CLI files 
# if no files specified, process all .gz files in directory.  If no .gz files, show a help message
unless (@in_files) {
	printout("Searching for files ...\n");
	@in_files = glob("*.gz");
}
help() unless @in_files;

# Read Index file
printout("\n\n*** Reading index file ...\n");
if (-e $index_file) {
	open (NORM, "<$index_file") || die "Unable to open normalisation file for reading: $!\n";
	my @norm = <NORM>;
	chomp(@norm);
	printout("\nIndex\tName\n");
	foreach my $l (@norm) {
		my ($i, $name) = split(/\s+/,$l);
		printout("$i\t$name\n");
		$i =~ s/A//;
		$index{$i}=$name;
	}
}

printout("\n\n*** Matching adaptors ...\n");
foreach my $l (@in_files) {
	print "$l\n";
	
	## Change this next line's regexp to match your sequencing format (currently matches eg "Index6" or "A6")
	my ($i) = $l =~ m/.*?(?:index|a0*)+(\d+)/i;
	if ($index{$i}) {
		printout("$index{$i}: Index $i\n\n");
		$files{$index{$i}}[0]=$l;
		
		if ($index{$i} =~ m/^dam$/i) {
			die("Error: more than one Dam sample detected.  Please only use one Dam control per run.\n\n") if $damname;
			$damname = $index{$i}
		}
	}
}

# Check that there's a Dam control sample ...
die ("Error: no Dam control sample detected!\n\n") unless $damname;

align_sequences();
extend_reads();
load_gatc_frags();
calc_bins();

for my $i (keys %files) {
	find_quants($i);
}

quantize_data();
normalize();
generate_ratio();

printout("All done.\n\n");




sub align_sequences {
	if ($vars{'bowtie'}) {
		
		printout("\n\n*** Aligning files with bowtie2 ...\n");
		foreach my $fn (keys %files) {
			my $pair1 = $files{$fn}[0];
			
			printout("\nNow working on $fn ...\n");
			
			if ($vars{'bowtie'}) {
				$bowtie_output{$fn} = `$vars{'bowtie2_path'}bowtie2 -p $vars{'threads'} -x $bowtie_genome_file -U $pair1 -S $fn.sam 2>&1`;
				printout("$bowtie_output{$fn}\n");
			}
		}
	}
}

sub extend_reads {
	# Extend reads
	if ($vars{'extend_reads'}) {
		printout("\n\n*** Extending reads to $vars{'len'} bases ...\n\n");
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
			
			printout("\nSeqs extended (>q30) = $seqs\n");
			
			# Store numbers of reads for normalisation later
			$counts{$fn}=$seqs;
			if ($fn =~ m/^dam$/i) {
				$denom = $counts{$fn};
			}
			
			close OUT;
		}
	}
}

sub calc_bins {
	printout("\n\n*** Calculating bins ...\n");
	foreach my $n (keys %files) {
		my $fn = "$n-ext$vars{'len'}";
		
		printout("\nNow working on $fn ...\n\n");
		
		if ($vars{'extend_reads'}) {
			printout("Generating .bam file ...\n");
			`$vars{'samtools_path'}samtools view -Sb $fn.sam > $fn.bam`;
			
			printout("Sorting ...\n");
			`$vars{'samtools_path'}samtools sort $fn.bam $fn-sorted`;
		}
		
		my $fout = "$fn-sorted.$vars{'bins'}.gff";
		
		printout("Generating bins from $fn.bam ...\n");
		
		my @a = `$vars{'coveragebed_path'}coverageBed -abam $fn-sorted.bam -b $windows_file`;
		
		die "Unable to process data!" unless @a;
		
		my @ar;
		foreach my $l (@a) {
			my ($chr, $start, $end, $score) =split("\t", $l);
			$chr=~s/chr//g;
			
			unless (defined($start)) {
				printout("Warning -- misread line:\n$l\n\n");
			}
			
			# move everything into an array from here on in ...
			push @ar, [ $chr, $start, $end, $score];

		}
		# Clear some memory ...
		@a=();
		
		printout("Sorting ...                 \n");
		@ar = sort { $a->[1] <=> $b->[1] } @ar; # sort start
		@ar = sort { $a->[0] cmp $b->[0] } @ar; # sort chr
		
		@{$full_tracks{$n}} = @ar if $vars{'full_data_files'};
		
		my $gatc_nonhits = 0;
					
		my $last_input=0;
		my @old_data;
		foreach my $i (0 .. ($#gatc_simple-1)) {
			# Get two adjacent GATC sites
			my ($chra, $mida) = @{$gatc_simple[$i]};
			my ($chrb, $midb) = @{$gatc_simple[$i+1]};

			if ($i%100 == 0) {
				my $pc = sprintf("%0.2f",($i*100)/$#gatc_simple);
				my $non_pc = sprintf("%0.2f",($gatc_nonhits*100)/($i+1));
				print STDERR "$pc\% processed ... [$non_pc% missing frags]\r";
			}
			next unless $chra eq $chrb;
						
			my $sum;
			my @data;
			my $count;
			# Run through the probe array to find the probes that lie within the fragment
			foreach my $l ($last_input .. $#ar-1) {
				# Takes way too long to scan through the entire array for each GATC fragment, and is in theory unnecessary
				# This proceedure uses $last_input to store the last probe found within an array, and start the next search from there
				# (relies on an ordered array of data, hence the sorting above)
				
				my ($chr, $start, $end, $score) = @{$ar[$l]};
				my ($chrn, $startn, $endn, $scoren) = @{$ar[$l+1]};
				
				next unless $chr eq $chra;
				last unless $chrn eq $chra; # Short circuit if we've reached the end of the array for this chromosome
							
				$last_input=$l-1;

				last if $start > $midb;
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
			my $mean_rnd = sprintf("%0.3f", $mean);
			
			push(@{$gatc_frags{"$chra-$mida-$midb"}}, $n);
			$gatc_frag_size{"$chra-$mida-$midb"}=$midb-$mida;
			$gatc_frag_score{"$chra-$mida-$midb"}{$n}=$mean_rnd;
			push @{$array{$n}}, [ $chra, $mida, $midb, $mean_rnd];
		}
	}
}

sub normalize {
	printout("\n\n*** Normalising counts and adding background ...\n");
	foreach my $n (keys %files) {
		printout("\nProcessing input: $n ...\n");
		
		my $norm;
		unless ($n =~ m/^dam$/i) {
			if ($vars{'norm_override'}) {
				printout("Normalisation override!\n  Would have normalised by $norm_factors{$n}\n  ... will instead normalise by $vars{'norm_override'}\n");
				$norm = $vars{'norm_override'};
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
		
		printout("\nNow working on $n ...\n\n"); 

		my %data;
		
		# find the lowest number of compared counts -- this will now use a different number of pseudocounts per sample
		my $psc_min = ($counts{$n}, $denom)[$counts{$n} > $denom];
		my $pseudocounts = ($vars{'pseudocounts'} ? $vars{'pseudocounts'} : 10*$psc_min/@{$full_tracks{$n}}); # pseudocounts value is related to total reads/number of bins
		printout("  ... adding $pseudocounts pseudocounts to each sample\n");
		
		printout("Reading Dam ...\n");
		foreach (@{$array{$damname}}) {
			my ($chr, $start, $end, $score) = @{$_};
			push @{ $data{$chr}{$start}}, $end;
			push @{ $data{$chr}{$start}}, $score+$pseudocounts;
		}

		printout("Reading $n ...\n");
		foreach (@{$array{$n}}) {
			my ($chr, $start, $end, $score) = @{$_};
			push @{ $data{$chr}{$start}}, $score+$pseudocounts;
		}
		
		my $fout="$n-vs-$damname.gatc.gff";
		open (OUT, ">$fout") || die "Unable to open $fout for writing: $!\n";
		foreach my $chr (sort keys %data) {
			foreach my $start (sort {$a <=> $b} keys $data{$chr}) {
				my ($end, $score1, $score2) = @{ $data{$chr}{$start}};
				next unless defined($score2);
				my $score = log($score2 /$score1)/log(2);
			
				print OUT join("\t",$chr, '.', '.', $start, $end, $score, '.', '.', '.'), "\n";
			}
		}
		close OUT;
	}
	
	if ($vars{'full_data_files'}) {
		foreach my $n (keys %files) {
			next if $n =~ m/^dam$/i;
			
			printout("\nNow working on $n (full track) ...\n\n"); 

			my %data;
			
			# find the lowest number of compared counts -- this will now use a different number of pseudocounts per sample
				my $psc_min = ($counts{$n}, $denom)[$counts{$n} > $denom];
				my $pseudocounts = ($vars{'pseudocounts'} ? $vars{'pseudocounts'} : 10*$psc_min/@{$full_tracks{$n}}); # pseudocounts value is related to total reads/number of bins
				printout("  ... adding $pseudocounts pseudocounts to each sample\n");
			
			printout("Reading Dam ...\n");
			foreach (@{$full_tracks{$damname}}) {
				my ($chr, $start, $end, $score) = @{$_};
				push @{ $data{$chr}{$start}}, $end;
				push @{ $data{$chr}{$start}}, $score+$pseudocounts;
			}

			printout("Reading $n ...\n");
			foreach (@{$full_tracks{$n}}) {
				my ($chr, $start, $end, $score) = @{$_};
				push @{ $data{$chr}{$start}}, $score+$pseudocounts;
			}
			
			my $fout="$n-vs-$damname.gff";
			open (OUT, ">$fout") || die "Unable to open $fout for writing: $!\n";
			foreach my $chr (sort keys %data) {
				foreach my $start (sort {$a <=> $b} keys $data{$chr}) {
					my ($end, $score1, $score2) = @{ $data{$chr}{$start}};
					next unless defined($score2);
					my $score = log($score2 /$score1)/log(2);
				
					print OUT join("\t",$chr, '.', '.', $start, $end, $score, '.', '.', '.'), "\n";
				}
			}
			close OUT;
		}
	}
}

sub quantize_data {
		
	foreach my $n (keys %files) {
		
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
			
			if (($qscore{1}>=$vars{'qscore1min'}) && ($qscore{1}<=$vars{'qscore1max'}) &&  ($qscore{2}<=$vars{'qscore2max'})) {
				push @ratios, $score{2}/$score{1};
			}
			
		}
		my $total_measurements = @ratios;
		my $avg = stdev(@ratios);
		my $norm = 1/$avg;
		printout("Norm factor = $norm based off $total_measurements frags (total $total_frags)\n\n");
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

	foreach my $k (keys %gatc_frag_score) {
		my $score = $gatc_frag_score{$k}{$prot};
		next unless $score >0;
		#~ $score = 0 if $score eq "NA";
		push (@frags, $score);
	}
	
	printout("Sorting frags ...\n");
	my @sorted_frags = sort {$a<=> $b} @frags;
	
	my @quants;
	for (my $q=0.1;$q<=1.05;$q+=0.1) {
		if ($q >= (0-0.01)) {
			push @quants, $q;
		}
	}
	
	printout("Finding quants ...\n");
	my $cut_off = shift(@quants);
	my $count;
	foreach  (@sorted_frags) {
		my $score = $_;
		$count++;
		
		my $prop = $count/@sorted_frags;

		if ($prop>$cut_off) {
			printout("   *** found quant $cut_off at score $score\n");
			$seg{$cut_off}{$prot} = $score;
			$cut_off = shift(@quants);
		}
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

sub help {
	print STDERR "Default variables:\n\n";
	foreach (sort (keys %vars)) {
		print STDERR "$_ = $vars{$_}\n\t$vars_details{$_}\n\n";
	}
	print STDERR "\n";
	exit 1;
}

sub load_gatc_frags {
	printout("\nReading GATC file ...\n");
	open (GATC, "<","$vars{'gatc_frag_file'}") || die "Unable to read GATC file: $!\n";
	@gatc = (<GATC>);
	chomp (@gatc);
	close GATC;
	$gatc_fragments=@gatc;

	foreach my $l (@gatc) {
		my ($chr, $source, $type, $start, $end, $score, $b, $c, $name) = split('\t', $l);
		my $mid = ($start+$end)/2;
		push (@gatc_simple, [$chr, $mid]);
		$gatc_reverse_hash{$mid}=$#gatc_simple;
		push (@{$gatc_chr{$chr}}, $mid);
	}
}