# Introduction

[damidseq_pipeline](https://github.com/owenjm/damidseq_pipeline/releases) is a single script that automatically handles sequence alignment, read extension, binned counts, normalisation, pseudocount addition and final ratio file generation. The script uses FASTQ or BAM files as input, and outputs the final log2 ratio files in bedGraph format.

## Features
* Fully automated processing of all NGS DamID-seq datasets, from FASTQ input to bedGraph output
* Automatic grouping and processing of multiple experimental and replicate batches
* Automatically detects and processs both single- and paired-end datasets
* Can be used with either FASTQ or pre-aligned BAM input files
* Multiple methods of normalisation provided
* Optional generation of outputs for CATaDa processing
* Can also handle and process ChIP-seq or CUT&RUN NGS data
* Free and open-source software maintained by the [Marshall lab](https://marshall-lab.org)

## Citation

If you find this software useful, please cite:

Marshall OJ and Brand AH. (2015) damidseq_pipeline: an automated pipeline for processing DamID sequencing datasets. *Bioinformatics.* 31(20): 3371--3.
([pubmed](http://www.ncbi.nlm.nih.gov/pubmed/26112292); [full text, open access](https://academic.oup.com/bioinformatics/article/31/20/3371/196153))

Please note that damidseq_pipeline has now evolved well beyond the functionality described in that article.

# Download and installation

[Download the latest version](https://github.com/owenjm/damidseq_pipeline/releases) of the pipeline script and associated files.

Prebuilt GATC fragment files used by the script are available for the following genomes:
* [*D. melanogaster* r6](https://github.com/owenjm/damidseq_pipeline/raw/gh-pages/pipeline_gatc_files/Dmel_BDGP6.GATC.gff.gz)
* [*Drosophila melanogaster* r5.57](https://github.com/owenjm/damidseq_pipeline/raw/gh-pages/pipeline_gatc_files/Dmel_r5.57.GATC.gff.gz)
* [*Mus musculus* GRCm38](https://github.com/owenjm/damidseq_pipeline/raw/gh-pages/pipeline_gatc_files/MmGRCm38.GATC.gff.gz) or
* [Human GRCh38](https://github.com/owenjm/damidseq_pipeline/raw/gh-pages/pipeline_gatc_files/HsGRCh38.GATC.gff.gz).

## Requirements

* [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) v2.1 or above (and appropriate genome indices -- see below) (not required if using pre-aligned BAM files)
* [SAMtools](http://samtools.sourceforge.net) v0.1.9 or above
* a GFF file containing all GATC sites in the genome (a file for the Drosophila genome is provided in the script .zip archive)
* a *nix operating system (e.g. linux, Mac OSX) with Perl v5.10 or greater installed. We recommend using Ubuntu Linux in a virtual machine if using Windows.
* Sequencing data in FASTQ or BAM format

## Installation

1. Extract the pipeline script archive, make the `damid_pipeline` file executable and place it in your path
    ```bash
    # Very simple way to do this in a *nix environment,
    # Change to the directory with the extracted files and:
    chmod a+x damidseq_pipeline
    sudo cp damidseq_pipeline /usr/local/bin/
    ```
1. Install [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
1. Obtain Bowtie 2 indices provided by [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) or [Illumina's iGenome](http://support.illumina.com/sequencing/sequencing_software/igenome.html)

    Alternatively, build the Bowtie 2 index files manually:
    1. Download the latest FASTA genome primary_assembly (or toplevel) file from [Ensembl](http://ftp.ensembl.org/pub/current_fasta/)
        e.g. [the current release for *Mus musculus*](http://ftp.ensembl.org/pub/current_fasta/mus_musculus/dna/)
        
        (alternatively, for *Drosophila*, download from [the Flybase FTP site](ftp://ftp.flybase.net/releases/current/))
    1. Extract the .gz file
    1. Run bowtie2-build in the directory containing the extracted .fasta file. For example:

            bowtie2-build Mus_musculus.GRCm38.dna.primary_assembly.fa GRCm38
            bowtie2-build dmel-all-chromosome-r5.57.fasta dmel_r5.57
1. Install [SAMtools](http://samtools.sourceforge.net)
1. Download a pre-built GATC fragment file for
    * [*D. melanogaster* r6](https://github.com/owenjm/damidseq_pipeline/raw/gh-pages/pipeline_gatc_files/Dmel_BDGP6.GATC.gff.gz)
    * [*D. melanogaster* r5.57](https://github.com/owenjm/damidseq_pipeline/raw/gh-pages/pipeline_gatc_files/Dmel_r5.57.GATC.gff.gz) _(... surely nobody's still using release 5, though, right?)_
    * [*Mus musculus* GRCm38](https://github.com/owenjm/damidseq_pipeline/raw/gh-pages/pipeline_gatc_files/MmGRCm38.GATC.gff.gz) or
    * [Human GRCh38](https://github.com/owenjm/damidseq_pipeline/raw/gh-pages/pipeline_gatc_files/HsGRCh38.GATC.gff.gz).
    
    Alternatively, build your own:

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

By default, the pipeline will process all files found in the current working directory.  To process files in a different directory, specify the path with the `--datadir` command-line switch (from v1.6).

Alternatively, individual files may be specified on the command line if the user does not wish to process all available files present in the directory (but there is little reason to do this, as from v1.6 the pipeline can correctly group and process multiple replicates and experiments from the one batch of files).

## Sample names
Sample names are assigned from filenames.  Ideally, start each sample with the name of the protein being profiled: if you do this, everything else will be easy.

If a single filename being processed begins with "Dam", this will be assigned as the Dam-only control. 

If no sample filename, or multiple filenames, begin with "Dam", use `--dam=[filename]` to specify the Dam-only control sample manually.  If a Dam-only control cannot be automatically determined, damidseq_pipeline will exit and prompt you to specify one.  (But, please, save yourself the trouble and just start your filenames with the protein name.  You'll thank me later.)

## Processing multiple experiments and replicates

As of v1.6, damidseq_pipeline will group files into different experiments and replicates, and handle this automatically.  (There's no more need for complex snakemake workarounds or shell-scripted loops if you just want to align and process a set of samples.)

Experimental and replicate group detection relies on at least two parameters:
* `--exp_prefix`: the common characters immediately preceding the experiment name (default `_`)
* `--rep_prefix`: the common characters that prefix the replicate number (default `_n`)

However, if you have additional information between the experiment name and the replicate designation, you can also set:
* `--exp_suffix`: the common characters immediately following the experiment name (takes the value of `--rep_prefix` when unset)

In the Marshall lab, we typically use the following naming format:

  `[protein]_[celltype]-[experiment number]_n[replicate]`

Thus, the following filenames will all work with the default values above:

```
Dpn_NSCs-OM1_n1
Dpn_NSCs-OM1_n2
Dam_NSCs-OM1_n1
Dam_NSCs-OM1_n2
HP1a_Neurons-OM2_n1
Dam_Neurons-OM2_n1
```

But you can designate your own character strings to fit your own naming conventions if you need to, using the command-line switches above.

We think this makes processing much easier (especially these days when you can multiplex 50+ samples in a lane of sequencing).  But, if you prefer to keep things simple and want the old v1.5.3 and earlier functionality back, you can always run with the command-line option `--nogroups` and everything will be like it was before.
{: .notice--info}

## Paired-end sequencing files
As of v1.6, damidseq_pipeline will automatically detect paired-end or single-end reads and process these appropriately.  You can happily mix read types within a single processing run if you want or need to.

What form of sequencing do we suggest?  Paired-end sequencing is always better if you can afford it (and really there should be very little price difference these days).  Paired-ends provides more accurate alignments and allows significantly better alignment within repetitive regions.
{: .notice--info}

If you're using an earlier version of the pipeline and don't want to update, use the `--paired` command-line option and the pipeline will search for, and match, paired reads.  But we recommend updating to the latest version (or if there's a bug preventing you from updating, let us know!).

## Generating coverage bedGraph files for CATaDa

As of v1.6, the `--catada` command-line option will output binned coverage tracks (in RPM, reads per million mapped reads, by default) for individual BAM or FASTQ files and exit (i.e. no ratio files will be generated).  (Earlier versions provided the same functionality with the `--just_coverage` flag, which is still present and works identically to `--catada`.)

To generate ratios _and_ output coverage files, use the `--coverage` flag instead, and get two wishes in one.

Coverage bedGraphs from Dam-only files generated in this way can be used for [Chromatin Accessibilty TaDa (CATaDa)](https://elifesciences.org/articles/32341) processing.

## Processing ChIP-seq or CUT&RUN data
As of v1.5.3 and later, damidseq_pipeline can also handle ChIP-seq or CUT&RUN data via the `--chipseq` flag.  This option will remove PCR duplicate reads, only process uniquely mapping reads, and output binned coverage tracks in RPM (reads per million mapped reads).

Warning -- do not use this option, or attempt to remove PCR duplicates, with DamID-seq data.  **DamID-seq is _all_ about the PCR duplicates!!**  We can't emphasise this enough.
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

## Visualising DNA binding profiles

The bedgraph output files can be can viewed directly in genome browsers such as [IGV](http://www.broadinstitute.org/software/igv/).  For publication-quality figures, we recommend [pyGenomeTracks](https://pygenometracks.readthedocs.io/).

## Calling significant peaks from the data

The [find_peaks](http://github.com/owenjm/find_peaks) software will process the output .gatc.bedgraph ratio file and call significant peaks present in the dataset.  Please see the find_peaks page for more details.

## Calling transcribed genes from RNA pol II datasets

The [polii.gene.call](http://github.com/owenjm/polii.gene.call) Rscript will call transcribed genes (i.e. gene bodies with significantly enriched pol II occupancy) from the output .gatc.bedgraph file.  Please see the polii.gene.call page for more details.

## Comparative protein binding, gene expression and chromatin accessibility analysis

We should have a comprehensive downstream damidBind R package for downstream processing of protein binding, Pol II occupancy and CATaDa chromatin accessibility released shortly.  Watch this space.

# DamID, TaDa, FlyORF-TaDa, NanoDam and CATaDa

No matter what flavour of DamID you're using, if there's a DNA Adenine Methylase enzyme involved, `damidseq_pipeline` is the tool to use.

For every technique other than CATaDa, process your samples as per normal.

For CATaDa, use `--catada` to generate chromatin accessibility coverage bedGraphs.

# Reporting issues, feature requests, and bugs

Please log these via the [damidseq_pipeline GitHub site](https://github.com/owenjm/damidseq_pipeline/).

# DamID protocols and reagents

For our latest lab protocol for Targeted DamID, advice and base plasmid files and sequences, please see our [dedicated page on TaDa](https://marshall-lab.org/tada) for more details.

