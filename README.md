# Introduction

[damidseq_pipeline](https://github.com/owenjm/damidseq_pipeline/releases) is a single script that automatically handles sequence alignment, read extension, binned counts, normalisation, pseudocount addition and final ratio file generation. The script uses FASTQ or BAM files as input, and outputs the final log2 ratio files in bedGraph format.

## Features
* Fully automated processing of NGS DamID-seq datasets, from FASTQ input to bedGraph output
* Handles both single- and paired-end datasets
* Can be used with either FASTQ or pre-aligned BAM input files
* Multiple methods of normalisation provided
* As of v1.5.3 and greater, can also handle and process ChIP-seq NGS data

## Citation

If you find this software useful, please cite:

Marshall OJ and Brand AH. (2015) damidseq_pipeline: an automated pipeline for processing DamID sequencing datasets. *Bioinformatics.* 31(20): 3371--3.
([pubmed](http://www.ncbi.nlm.nih.gov/pubmed/26112292); [full text, open access](https://academic.oup.com/bioinformatics/article/31/20/3371/196153))

# Download and installation

[Download the latest version](https://github.com/owenjm/damidseq_pipeline/releases) of the pipeline script and associated files.

Prebuilt GATC fragment files used by the script are available for the following genomes:
* [*Drosophila melanogaster* r5.57](https://github.com/owenjm/damidseq_pipeline/raw/gh-pages/pipeline_gatc_files/Dmel_r5.57.GATC.gff.gz)
* [*D. melanogaster* r6](https://github.com/owenjm/damidseq_pipeline/raw/gh-pages/pipeline_gatc_files/Dmel_BDGP6.GATC.gff.gz)
* [*Mus musculus* GRCm38](https://github.com/owenjm/damidseq_pipeline/raw/gh-pages/pipeline_gatc_files/MmGRCm38.GATC.gff.gz) or
* [Human GRCh38](https://github.com/owenjm/damidseq_pipeline/raw/gh-pages/pipeline_gatc_files/HsGRCh38.GATC.gff.gz).

## Requirements

* [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) v2.1 or above (and appropriate genome indices -- see below) (not required if using pre-aligned BAM files)
* [SAMtools](http://samtools.sourceforge.net) v0.1.9 or above
* a GFF file containing all GATC sites in the genome (a file for the Drosophila genome is provided in the script .zip archive)
* a *nix operating system (e.g. linux, Mac OSX) with Perl v5.10 or greater installed. We recommend using Ubuntu Linux in a virtual machine if using Windows.
* Sequencing data in FASTQ or BAM format

## Installation

1. Extract the pipeline script archive, make the damid_pipeline file executable and place it in your path
1. Install [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
1. Obtain Bowtie 2 indices provided by [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) or [Illumina's iGenome](http://support.illumina.com/sequencing/sequencing_software/igenome.html)

    Alternatively, build the Bowtie 2 index files manually:
    1. Download the latest FASTA genome primary_assembly (or toplevel) file from [Ensembl](http://ftp.ensembl.org/pub/current_fasta/)
        e.g. [the current release for *Mus musculus*](http://ftp.ensembl.org/pub/current_fasta/mus_musculus/dna/)
        
        (alternatively, for *Drosophila*, download from [the Flybase FTP site](ftp://ftp.flybase.net/releases/current/)
    1. Extract the .gz file
    1. Run bowtie2-build in the directory containing the extracted .fasta file. For the examples above:

            bowtie2-build Mus_musculus.GRCm38.dna.primary_assembly.fa GRCm38
            bowtie2-build dmel-all-chromosome-r5.57.fasta dmel_r5.57
1. Install [SAMtools](http://samtools.sourceforge.net)
1. Download a pre-built GATC fragment file for
    * [*D. melanogaster* r5.57](https://github.com/owenjm/damidseq_pipeline/raw/gh-pages/pipeline_gatc_files/Dmel_r5.57.GATC.gff.gz)
    * [*D. melanogaster* r6](https://github.com/owenjm/damidseq_pipeline/raw/gh-pages/pipeline_gatc_files/Dmel_BDGP6.GATC.gff.gz)
    * [*Mus musculus* GRCm38](https://github.com/owenjm/damidseq_pipeline/raw/gh-pages/pipeline_gatc_files/MmGRCm38.GATC.gff.gz) or
    * [Human GRCh38](https://github.com/owenjm/damidseq_pipeline/raw/gh-pages/pipeline_gatc_files/HsGRCh38.GATC.gff.gz).
    
    Alternatively build your own:

    1. Download the FASTA genome sequence, as in step 3 above (no need to extract the gzipped files)
    1. Run the provided [gatc.track.maker.pl](http://github.com/owenjm/damid_pipeline/blob/master/gatc.track.maker.pl?raw=true) script on the fasta sequence, e.g.:

            perl gatc.track.maker.pl --name=dmel_r5.57 dmel-all-chromosome-r5.57.fasta

## Initial setup

In order to run correctly, damidseq_pipeline needs to know the locations of two paths, specified using the following command-line options (it will prompt you for these if you fail to specify them):

1. The directory and basename of the bowtie2 index files (obtained or built in step 3 above)
    (specified with the `--bowtie2_genome_dir` option)
        e.g. in the example above, use

        --bowtie2_genome_dir=[path_to_.bt2_files]/dmel_r5.57
1. The GATC fragment .gff file (downloaded from the pre-built files listed in step 5, or built following the instructions above)
    (specified with the `--gatc_frag_file` option)

In order to setup the pipeline to process the *D. melanogaster* genome, for example, the first-run command would be:

    damidseq_pipeline --gatc_frag_file=path/to/Dmel_r5.57.GATC.gff.gz --bowtie2_genome_dir=path/to/dmel_r5.57/dmel_r.5.57

Once run once with these options and correct values, the paths will be saved for all future runs unless overridden on the command-line. **You do not need to specify these options each time**
{: .notice--info}

(To clear all user-saved values, including these values, run with the `--reset_defaults` option.)

# Using damidseq_pipeline

Run damidseq_pipeline in a directory containing sequencing files in FASTQ or BAM format.  The default behaviour is to process all files in FASTQ format, and if none are found, all files in BAM format.  

Alternatively, individual files may be specified on the command line if the user does not wish to process all available files present in the directory (for example, if the sequencing lane contained multiple replicates).

## Sample names
Sample names are assigned from filenames.  If a single filename being processed begins with "Dam", this will be assigned as the Dam-only control.  

If no sample filename, or multiple filenames, begin with "Dam", use `--dam=[filename]` to specify the Dam-only control sample manually.  If a Dam-only control cannot be automatically determined, damidseq_pipeline will exit and prompt you to specify one.

## Paired-end sequencing files
To process paired-end FASTQ files, use the `--paired` option and the pipeline will search for, and match, paired reads.

BAM files generated from paired-end data are automatically detected and processed, without requiring this option.

## Processing ChIP-seq data
As of v1.5.3, damidseq_pipeline can also handle ChIP-seq data via the `--chipseq` flag.  This option will remove PCR duplicate reads, only process uniquely mapping reads, and output binned coverage tracks in RPM (reads per million mapped reads).

Warning -- do not use this option with DamID-seq data.  DamID-seq is all about the PCR duplicates!
{: .notice--warning}

## Other options
To see all available options, run the script with `--help` command-line option:

    damidseq_pipeline --help

This will give you a list of adjustable parameters and their default and current values if applicable. We recommend keeping these at the default value in most cases; however, these can be modified on the command-line with `--option=value` (no spaces).

## Saving default option values
To save modified values for all future runs, run the script with the parameter you wish to change together with the `--save_defaults` command-line option:

    damidseq_pipeline --save_defaults

If bowtie2 and samtools are not in your path, you can specify these on the command-line also.

## Working with multiple genomes

If the user expects to process data from multiple genomes, separate genome specifications can be saved by using the `--save_defaults=[name]` along with the `--bowtie2_genome_dir` and `--gatc_frag_file` options (and any other custom options that the user wishes to set as default for this genome, e.g. the bin width).  For e.g.:

    damidseq_pipeline --save_defaults=fly --gatc_frag_file=path/to/Dmel_r5.57.GATC.gff.gz --bowtie2_genome_dir=path/to/dmel_r5.57/dmel_r.5.57
    damidseq_pipeline --save_defaults=mouse --bins=500 --gatc_frag_file=path/to/MmGRCm38.GATC.gff.gz --bowtie2_genome_dir=path/to/Mm_GRCm38/GRCm38

Once set up, different genome definitions can be quickly loaded using the --load_defaults=[name] option, e.g.:

    damidseq_pipeline --load_defaults=fly

All currently saved genome definitions can be listed using --load_defaults=list.

## Example dataset

An example set of two small (3000 reads each) fastq.gz files and an index.txt file are provided in the zip archive (as "example.zip"), or you can [download these separately](http://github.com/owenjm/damid_pipeline/blob/master/example.zip?raw=true). Running the pipeline script on these files should successfully produce a polII-vs-Dam.gatc.bedgraph ratio file as output.

# Output and downstream data processing

The final output will be a single ratio file: Sample-vs-Dam.gatc.bedgraph. The .gatc.bedgraph file represents the data at GATC fragment resolution (based on the reference genome) and should be used for all subsequent analysis.

The [bedGraph format](http://genome.ucsc.edu/goldenpath/help/bedgraph.html) is used by default.  The pipeline script can output the final ratio files in [GFF format](http://www.ensembl.org/info/website/upload/gff.html) instead if the `--output_format=gff` command-line switch is used.

### Visualising the DNA binding profiles

The bedgraph output files can be can viewed directly in genome browsers such as [IGV](http://www.broadinstitute.org/software/igv/).  For publication-quality figures, we recommend [pyGenomeTracks](https://pygenometracks.readthedocs.io/).

### Calling significant peaks from the data

The [find_peaks](http://github.com/owenjm/find_peaks) software will process the output .gatc.bedgraph ratio file and call significant peaks present in the dataset.  Please see the find_peaks page for more details.

### Calling transcribed genes from RNA pol II datasets

The [polii.gene.call](http://github.com/owenjm/polii.gene.call) Rscript will call transcribed genes (i.e. gene bodies with significantly enriched pol II occupancy) from the output .gatc.bedgraph file.  Please see the polii.gene.call page for more details.

### Other useful scripts and utilities

A collection of useful R and Perl scripts for comparing and analysing DamID-seq data [is maintained here](http://github.com/owenjm/damid_misc).
