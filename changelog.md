## v1.6.2
*  Fixed: handling of paired-end vs. single-end reads (hotfix for errors introduced in v1.6.1)

## v1.6.1
*  Changed: improved handling of paired-end vs. single-end reads
*  Fixed: better error handling of issues creating the config folder in the user's home directory
*  Fixed: removed the need to specify bowtie2 indices when providing BAM files

## v1.6
*  New feature: automatic grouping and processing of experiments and replicates within the set of input files (if you don't want to use this, run with `--nogroups` to restore old v1.5.3 functionality).  This removes the need for any higher-level scripting of sample processing, and means that a large set of disparate samples can be processed with just a single command. Experimental and replicate group detection relies on three parameters:
    * `--exp_prefix`: the common characters immediately preceding the experiment name (default `_`)
    * `--exp_suffix`: the common characters immediately following the experiment name (if unset, takes the value of `--rep_prefix`)
    * `--rep_prefix`: the common characters that prefix the replicate number (default `_n`)
*  New feature: dry-run with `--n` to check experiment and replicate groupings prior to execution
*  New feature: automatic detection of PE or SE sequencing inputs.  No need to explicitly use CLI flags anymore to specify; `damidseq_pipeline` will detect these automatically.  This should work out of the box for most sequencing facility outputs, but if you need to adjust the detection parameters use `--paired_match` to change the paired-end prefix search string.  Including mixed PE and SE inputs in the same processing run is now completely fine.
*  New feature: use `--datadir` to specify the directory of FASTQs or BAMs to process.
*  New feature: optionally, do not process and overwrite previously generated BAM files when reprocessing a large number of samples.  Use `--ncbam` to enable this (ensure you manually delete any partly-aligned BAM files first if using this; the pipeline applies no sanity checks to detect part-aligned files).
*  New feature: `--catada` flag is an alias for `--just_coverage` (as apparently the built-in CATaDa file generation functionality wasn't very transparent in the old release)
*  Modified: Inline::C normalisation code compilation now saves the compiled C code to the user's .config directory for reuse
*  Fixed: explicitly prevent any form of read extension when using PE sequencing data
*  Fixed: GATC fragment fence-post error corrected (fragments in output bedgraphs no longer have a 1bp overlap)
*  Some minor code refactoring to allow experimental and replicate groupings


## v1.5.3
*  New feature: `--coverage` flag outputs coverage bedGraph track of each sample (useful with Dam-only tracks for CATaDa analysis).  `--just_coverage` exits after coverage tracks are written.
*  New feature: Kernel density normalisation routines sped up ~10x via Inline::C.  This is enabled by default, but will only work if the Inline::C module is installed; damidseq\_pipeline will fall back on the Perl coded routine if this is not present.  If you want to control explicitly, do so via the `--kde_inline` flag.
*  New feature: ChIP-seq-specific processing routines added: `--remdups` removes PCR duplicates, `--unique` only maps unique reads.  Enable all of these with `--chipseq` (also turns of GATC-level binning, GATC-based read extension and will only output coverage tracks not ratios)
*  New feature: New normalisation method, "rawbins".  This writes the normalisation data to TSV format and runs a user-provided script for custom normalisation.  Your script should read the data file and output the normalisation factor.  Advanced users only, use with caution.
*  New feature: consider all occupancy as RPM (reads per million mapped reads), set as default.  Control with `--scores_as_rpm` flag.
*  Changed: bowtie2 command handling refactored more elegantly
*  Added: ability to add custom flags to bowtie2 command with `--bowtie2_add_flags`
*  Fixed: small changes to read mapping code for more accurate mapping
*  Changed: filename filters are now turned off by default (use `--no_file_filters=0` to turn back on)


## v1.4.6
*  New feature: paired-end reads are now handled natively at the fastq level.  Use the --paired flag to enable, and your paired reads should be automagically detected and aligned
*  Changed: BAM file format is now used natively without SAM intermediates at all levels (also fixes the recent samtools version handling bug)
*  Changed: paired-end reads are selected based on the 0x2 bit on the SAM file format bitwise FLAG (reads aligned in proper pair)
*  Fixed: paired-end read counts are enumerated correctly
*  Removed: --keep_sam and --only_sam flags (as SAM format is no longer used)
*  Added: --keep_original_bams flag to retain the non-extended reads BAM file when using single-end sequencing
*  Changed: warning on missing chromosomal identities in GATC alignment made more friendly and less confusing

## v1.4.5
*  Added --ps_debug option to explicitly display pseudocounts calculation

## v1.4.4
*  Fixed issue with some externally-generated BAM files

## v1.4.3
*  Better handling of chromosome names

## v1.4.2
*  Small improvements to SAM file handling

## v1.4.1
*  Fixed issues with samtools v1.3 (this version of samtools introduced backwards incompatibilities when using the 'sort' function.  damidseq_pipeline now checks for the version number and should support all versions of samtools.)

## v1.4
*  New read-extension method: by default, reads are now only extended as far as the next GATC fragment.  Use --extension_method=full to disable this feature and extend every read by the value of --len.
*  Output format is now bedgraph by default.  Use --output_format=gff to restore the previous default.  Changing the default to bedgraph allows users to create TDF tracks directly within the graphical IGV tools, making it easier for end users.
*  Minor code cleanups

## v1.3
*  Major bugfix: reads from the minus strand were not being extended correctly during processing.  The overall impact is minor (correlation between old and new read extension methods is >0.95) but this new method is technically more accurate.
*  added --keep_sam (do not delete the temporary SAM file) option
*  added --only_sam (do not generate BAM files) option (both options are intended for debug purposes only)

## v1.2.7
*  New opition --no_file_filters to prevent any filename trimming/filtering (by default input filenames are trimmed to the first underscore)
*  Small filename issue fixes

## v1.2.6
*  Now uses File::Basename to handle filenames
*  Fixed/cleaned up a number of rare problems with filename handling

## v1.2.5
*  Added --dam and --out_name options.  Use these to set the Dam-only control sample and/or a custom output name
*  Added more sanity checks
*  Minor bugfixes

## v1.2.4
*  Added explicit checks for bowtie2 and samtools executables, and for bowtie2 output

## v1.2.3
*  Fixed a serious error in RPM normalisation calculations (values were inverted) -- please re-run on your samples if you have used this method on them
*  Minor code cleanups

## v1.2.1
*  Added ability to process BAM files generated from paired-end sequencing
*  Cleaner reporting of missing assembly fragments in GATC files
*  Some small bugs fixed
*  General code clean-ups

## v1.2
*  Completely re-written normalisation routine based on kernel density estimation
*  Genomic coverage is now calculated internally rather than using bedtools (uses much less memory, is slightly faster, and drops the requirement for an external binned windows file)
*  Binned window files are no longer required (bins are calculated automatically using the sequence information provided in the BAM headers, and the bin size specified by the --bins command-line option)
*  Better handling of GATC fragment files (should prevent hangs/pauses when creating GATC fragment arrays)
*  Memory optimisation for large files (greatly reduces usage for processing mouse/human data)
*  Added ability to process BAM files directly
*  Much better file-handling all round (now takes sample names directly from filenames by default; the option to use an index.txt file remains but is essentially deprecated)
*  New option: --norm_method=## kde/rpm  "kde" is the default method using kernel density estimation; "rpm" normalises solely on readcounts/million reads only (the "rpm" method is not recommended except for the very rare cases in which a Dam-fusion protein fails to methylate accessible genomic regions, making kde normalisation is inappropriate)
*  Re-written --help output rountines (better formatted and more informative)
*  Ability to read gzipped GATC files
*  Ability to save sets of defaults to enable quick switching between different genomes  (use --save_defaults=## name; use --load_defaults=## name to load; use load_defaults=list to list current available options)
*  New location for config files (in ~/.config/damid_pipeline/).  Existing config file will be migrated automatically
*  Various small bugfixes and code clean-ups

** NB: a number of default parameters have changed with this release.  It is strongly advised to reset all parameters to the default value with --reset_defaults.

## v1.0
*  Initial release
