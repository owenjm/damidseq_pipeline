### Introduction

Processing DamID-seq data involves extending single-end reads, aligning the reads to the genome and determining the coverage, similar to processing regular ChIP-seq datasets. However, as DamID data is represented as a log2 ratio of (Dam-fusion/Dam), normalisation of the sample and Dam-only control is necessary, and adding pseudocounts to mitigate the effect of background counts when represented as a ratio is highly recommended.

We use a single pipeline script to handle alignment, read extension, binned counts, normalisation, pseudocount addition and final ratio file generation. The script uses .fastq files as input, and outputs the final log2 ratio files in GFF or BEDGRAPH format. These files can easily be converted to TDF for viewing in [IGV](http://www.broadinstitute.org/software/igv/) with the provided [gff2tdf.pl](http://github.com/owenjm/damid_pipeline/blob/master/gff2tdf.pl?raw=true) script (see below).

### Download

Download the latest version of the pipeline script and associated files:
* [As a zipfile](https://github.com/owenjm/damidseq_pipeline/zipball/master)
* [As a tarball](https://github.com/owenjm/damidseq_pipeline/tarball/master)

### Requirements

* [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) v2.1 or above (and appropriate genome indices -- see below) (not required if using pre-aligned BAM files)
* [SAMtools](http://samtools.sourceforge.net) v0.1.9 or above
* a GFF file containing all GATC sites in the genome (a file for the Drosophila genome is provided in the script .zip archive)
* a *nix operating system (e.g. linux, Mac OSX) with Perl v5.10 or greater installed. We recommend using Ubuntu Linux in a virtual machine if using Windows.
* Sequencing data in FASTQ or BAM format

### Installation

1. Unzip the pipeline script zip file, make the damid_pipeline.pl file executable and place it in your path
1. Install [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
1. Obtain Bowtie 2 indices provided by [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) or [Illumina's iGenome](http://support.illumina.com/sequencing/sequencing_software/igenome.html)

    Alternatively, build the Bowtie 2 index files manually:
    1. Download the latest FASTA genome primary_assembly (or toplevel) file from [Ensembl](ftp.ensembl.org/pub/current_fasta/)
        e.g. [the current release for *Mus musculus*](http://ftp.ensembl.org/pub/current_fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz)
        
        (alternatively, for *Drosophila*, download from the [Flybase FTP site](http://ftp.flybase.net/releases/current/)
         e.g. [*D. melanogaster* release 5.57](http://ftp.flybase.net/releases/FB2014_03/dmel_r5.57/fasta/dmel-all-chromosome-r5.57.fasta.gz))
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

### Running the script the first time

In order to run correctly, the script needs to know the locations of two paths, specified using the following command-line options:

1. The directory and basename of the bowtie2 index files (obtained or built in step 3 above)
    (specified with the --bowtie2_genome_dir option)
        e.g. in the example above, use

        --bowtie2_genome_dir=[path_to_.bt2_files]/dmel_r5.57
1. The GATC fragment .gff file (provided in the zip file, or built in step 5 above)
    (specified with the --gatc_frag_file option)

In order to setup the pipeline to process the *D. melanogaster* genome, for example, the first-run command would be:

    damidseq_pipeline.pl --gatc_frag_file=path/to/Dmel_r5.57.GATC.gff.gz --bowtie2_genome_dir=path/to/dmel_r5.57/dmel_r.5.57

If these paths do not already exist and the script is run with these options and correct values, the paths will be saved for all future runs unless overridden on the command-line.

(To clear all user-saved values, run the script with the --reset_defaults option.)

### Running the script

Run the script in a directory with sequencing files in FASTQ or BAM format.  The default behaviour is to process all files in FASTQ format, and if none are found, all files in BAM format.  Alternatively, individual files may be specified on the command line if the user does not wish to process all available files present in the directory (for example, if the sequencing lane contained multiple replicates).

The script will by default determine sample names from the file names, and expects a single filename to start with "Dam" in order to automatically assign the Dam-only control.

To see all available options, run the script with --help command-line option:

    damidseq_pipeline.pl --help

This will give you a list of adjustable parameters and their default and current values if applicable. We recommend keeping these at the default value in most cases; however, these can be modified on the command-line with --option=value (no spaces).

To save modified values for all future runs, run the script with the parameter you wish to change together with the --save_defaults command-line option:

    damidseq_pipeline.pl --save_defaults

If bowtie2 and samtools are not in your path, you can specify these on the command-line also.

### Output

The final output will be two ratio files: Sample-vs-DAM.gff and Sample-vs-DAM.gatc.gff. The .gatc.gff file represents the data at GATC fragment resolution (based on the reference genome) and should be used for all subsequent analysis. The other ratio file contains the coverageBed bins (i.e. 75nt bins by default) and may be useful for data representation.

The [GFF format](http://www.ensembl.org/info/website/upload/gff.html) is used by default.  The pipeline script can output the final ratio files in [BedGraph format](http://genome.ucsc.edu/goldenpath/help/bedgraph.html) instead if the --output_format=bedgraph command-line switch is used.

Either file can be converted to .tdf format for viewing in [IGV](http://www.broadinstitute.org/software/igv/) via the provided [gff2tdf.pl](http://github.com/owenjm/damid_pipeline/blob/master/gff2tdf.pl?raw=true) script.

### Working with multiple genomes

If the user expects to process data from multiple genomes, separate genome specifications can be saved by using the --save_defaults=[name] along with the --bowtie2_genome_dir and --gatc_frag_file options (and any other custom options that the user wishes to set as default for this genome, e.g. the bin width).  For e.g.:

    damidseq_pipeline.pl --save_defaults=fly --gatc_frag_file=path/to/Dmel_r5.57.GATC.gff.gz --bowtie2_genome_dir=path/to/dmel_r5.57/dmel_r.5.57
    damidseq_pipeline.pl --save_defaults=mouse --bins=500 --gatc_frag_file=path/to/MmGRCm38.GATC.gff.gz --bowtie2_genome_dir=path/to/Mm_GRCm38/GRCm38

Once set up, different genome definitions can be quickly loaded using the --load_defaults=[name] option, e.g.:

    damidseq_pipeline.pl --load_defaults=fly

All currently saved genome definitions can be listed using --load_defaults=list.

### Processing FASTQ files with adaptor codes

The damidseq_pipeline will attempt to automatically determine the sample name from the filenames provided.  

Some sequencing facilities, however, may provide only adaptor index IDs rather than sample names.  If sequencing filenames are provided with index IDs rather than names, a file "index.txt" in the working directory can be provided with the format of [index] [sample_name], eg:

    A6 Dam
    A12 polII 

where A6 is the sequencing adaptor index. The sample name cannot contain spaces and there has to be one (and only one) sample called "Dam"; the adaptor index must to be referenced in the fastq filename (e.g. for "A6", either "Index6" or "A006" are expected in the filename). Please see the [provided example](http://github.com/owenjm/damid_pipeline/blob/master/example.zip?raw=true) for an illustration of the index.txt file format and matching file names.

### Example dataset

An example set of two small (3000 reads each) fastq.gz files and an index.txt file are provided in the zip archive (as "example.zip"), or you can [download these separately](http://github.com/owenjm/damid_pipeline/blob/master/example.zip?raw=true). Running the pipeline script on these files should successfully produce a polII-vs-Dam.gff ratio file as output.

