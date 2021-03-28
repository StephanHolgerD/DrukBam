# DrukBam
`DrukBam` is a straight forward program for plotting alignment files without leaving the commandline. It can be used with or without a reference fasta file and allows fast plotting multiple variants or regions of interest.

## reference free
<center><img src="exampleOutput/example_refFreeDirection_19_281367_281468.png" width="100%"/></center>

## including a reference
<center><img src="exampleOutput/example_refFreeDirection_19_281367_281468R.png" width="100%"/></center>

## split by read direction
<center><img src="exampleOutput/example_refFreeDirection_19_281367_281568D.png" width="100%"/></center>


## bigger span
<center><img src="exampleOutput/example_refFreeDirection_19_281067_282568DR.png" width="100%"/></center>




# Usage
<details>
  <summary>DrukBam vcf</summary>

  ```
  usage: DrukBam vcf [-h] -b BAM -v VCF [-p PADDING] [--highlight]
                     [--threads THREADS] [--maxcoverage MAXCOVERAGE]
                     [--direction] [--schematic] [--style STYLE] [--fasta FASTA]
                     [--outputdir OUTPUTDIR] [-i ID] [--chunksize CHUNKSIZE]
                     [--outfmt OUTFMT] [--outlineoff]

  optional arguments:
    -h, --help            show this help message and exit

  required arguments:
    -b BAM, --bam BAM     Pos. sorted and indexed bam file
    -v VCF, --vcf VCF     vcf file with variants of interest
    -p PADDING, --padding PADDING
                          number of nt around the variant
    --highlight           highlight the position of interest

  optional arguments:
    --threads THREADS     number of cpu's to run in paralell, ROI <1000 will
                          always use 1 core
    --maxcoverage MAXCOVERAGE
                          max cov to plot
    --direction           split reads by forward and reverse
    --schematic           plot no nucleotide, recommended for ROI>1000
    --style STYLE         different style options for the plot, provide .ini
                          file
    --fasta FASTA         fasta file for reference related plotting
    --outputdir OUTPUTDIR
                          directory for output
    -i ID, --id ID        output filename
    --chunksize CHUNKSIZE
                          max size of visualized area, can be increases but will
                          sow down calculation
    --outfmt OUTFMT       format of plot, choose between pdf,svg,png
    --outlineoff          plotting of read outline

```
</details>


<details>
  <summary>DrukBam region</summary>

  ```
  usage: DrukBam region [-h] -b BAM -c CHROMOSOME -s START -e END
                      [--threads THREADS] [--maxcoverage MAXCOVERAGE]
                      [--direction] [--schematic] [--style STYLE]
                      [--fasta FASTA] [--outputdir OUTPUTDIR] [-i ID]
                      [--chunksize CHUNKSIZE] [--outfmt OUTFMT] [--outlineoff]

optional arguments:
  -h, --help            show this help message and exit

required arguments:
  -b BAM, --bam BAM     Pos. sorted and indexed bam file
  -c CHROMOSOME, --chromosome CHROMOSOME
                        name of chromosome/contig
  -s START, --start START
                        start of region of interest
  -e END, --end END     end of the region of interest

optional arguments:
  --threads THREADS     number of cpu's to run in paralell, ROI <1000 will
                        always use 1 core
  --maxcoverage MAXCOVERAGE
                        max cov to plot
  --direction           split reads by forward and reverse
  --schematic           plot no nucleotide, recommended for ROI>1000
  --style STYLE         different style options for the plot, provide .ini
                        file
  --fasta FASTA         fasta file for reference related plotting
  --outputdir OUTPUTDIR
                        directory for output
  -i ID, --id ID        output filename
  --chunksize CHUNKSIZE
                        max size of visualized area, can be increases but will
                        sow down calculation
  --outfmt OUTFMT       format of plot, choose between pdf,svg,png
  --outlineoff          plotting of read outline


```
</details>


## Installing
`DrukBam` is available via pypi:
```

```

## Examples:

Samplot requires either BAM files or CRAM files as primary input. If you use
CRAM, you'll also need a reference genome. You can easily acquire a reference genome file with [GGD](https://github.com/gogetdata/ggd-cli), which is also available from conda.

### Basic use case
Using data from NA12878, NA12889, and NA12890 in the
[1000 Genomes Project](http://www.internationalgenome.org/about) (available in the test/data directory of samplot), we will
inspect a possible deletion in NA12878 at 4:115928726-115931880 with respect
to that same region in two unrelated samples NA12889 and NA12890.

The following command will create an image of that region:
```
time samplot plot \
    -n NA12878 NA12889 NA12890 \
    -b samplot/test/data/NA12878_restricted.bam \
      samplot/test/data/NA12889_restricted.bam \
      samplot/test/data/NA12890_restricted.bam \
    -o 4_115928726_115931880.png \
    -c chr4 \
    -s 115928726 \
    -e 115931880 \
    -t DEL

real	0m3.882s
user	0m3.831s
sys	0m0.328s

```

The arguments used above are:

`-n` The names to be shown for each sample in the plot

`-b` The BAM/CRAM files of the samples (space-delimited)

`-o` The name of the output file containing the plot

`-c` The chromosome of the region of interest

`-s` The start location of the region of interest

`-e` The end location of the region of interest

`-t` The type of the variant of interest

This will create an image file named `4_115928726_115931880.png`, shown below:

<img src="/doc/imgs/4_115928726_115931880.png">

### Gene and other genomic feature annotations

Gene annotations (tabixed, gff3 file) and genome features (tabixed, bgzipped, bed file) can be
included in the plots.

Get the gene annotations:
```
wget ftp://ftp.ensembl.org/pub/grch37/release-84/gff3/homo_sapiens/Homo_sapiens.GRCh37.82.gff3.gz
bedtools sort -i Homo_sapiens.GRCh37.82.gff3.gz \
| bgzip -c > Homo_sapiens.GRCh37.82.sort.gff3.gz
tabix Homo_sapiens.GRCh37.82.sort.gff3.gz
```

Get genome annotations, in this case Repeat Masker tracks and a mappability track:
```
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDukeMapabilityUniqueness35bp.bigWig
bigWigToBedGraph wgEncodeDukeMapabilityUniqueness35bp.bigWig wgEncodeDukeMapabilityUniqueness35bp.bed
bgzip wgEncodeDukeMapabilityUniqueness35bp.bed
tabix wgEncodeDukeMapabilityUniqueness35bp.bed.gz

curl http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz \
| bgzip -d -c \
| cut -f 6,7,8,13 \
| bedtools sort -i stdin \
| bgzip -c > rmsk.bed.gz
tabix rmsk.bed.gz
```

Plot:
```
samplot plot \
    -n NA12878 NA12889 NA12890 \
    -b samplot/test/data/NA12878_restricted.bam \
      samplot/test/data/NA12889_restricted.bam \
      samplot/test/data/NA12890_restricted.bam \
    -o 4_115928726_115931880.d100.genes_reps_map.png \
    -c chr4 \
    -s 115928726 \
    -e 115931880 \
    -t DEL \
    -d 100 \
    -T Homo_sapiens.GRCh37.82.sort.gff3.gz \
    -A rmsk.bed.gz wgEncodeDukeMapabilityUniqueness35bp.bed.gz
```

<img src="/doc/imgs/4_115928726_115931880.d100.genes_reps_map.png">

## Generating images from a VCF file
To plot images from structural variant calls in a VCF file, use samplot's
`vcf` subcommand. This accepts a VCF file and the BAM files of samples
you wish to plot, outputting images and an `index.html` page for review.

### Usage
<details>
  <summary> samplot vcf </summary>

  ```
usage: samplot vcf [-h] [--vcf VCF] [-d OUT_DIR] [--ped PED] [--dn_only]
                   [--min_call_rate MIN_CALL_RATE] [--filter FILTER]
                   [-O {png,pdf,eps,jpg}] [--max_hets MAX_HETS]
                   [--min_entries MIN_ENTRIES] [--max_entries MAX_ENTRIES]
                   [--max_mb MAX_MB] [--min_bp MIN_BP]
                   [--important_regions IMPORTANT_REGIONS] -b BAMS [BAMS ...]
                   [--sample_ids SAMPLE_IDS [SAMPLE_IDS ...]]
                   [--command_file COMMAND_FILE] [--format FORMAT] [--gff GFF]
                   [--downsample DOWNSAMPLE] [--manual_run]

optional arguments:
  -h, --help            show this help message and exit
  --vcf VCF, -v VCF     VCF file containing structural variants
  -d OUT_DIR, --out-dir OUT_DIR
                        path to write output PNGs
  --ped PED             path ped (or .fam) file
  --dn_only             plots only putative de novo variants (PED file
                        required)
  --min_call_rate MIN_CALL_RATE
                        only plot variants with at least this call-rate
  --filter FILTER       simple filter that samples must meet. Join multiple
                        filters with '&' and specify --filter multiple times
                        for 'or' e.g. DHFFC < 0.7 & SVTYPE = 'DEL'
  -O {png,pdf,eps,jpg}, --output_type {png,pdf,eps,jpg}
                        type of output figure
  --max_hets MAX_HETS   only plot variants with at most this many
                        heterozygotes
  --min_entries MIN_ENTRIES
                        try to include homref samples as controls to get this
                        many samples in plot
  --max_entries MAX_ENTRIES
                        only plot at most this many heterozygotes
  --max_mb MAX_MB       skip variants longer than this many megabases
  --min_bp MIN_BP       skip variants shorter than this many bases
  --important_regions IMPORTANT_REGIONS
                        only report variants that overlap regions in this bed
                        file
  -b BAMS [BAMS ...], --bams BAMS [BAMS ...]
                        Space-delimited list of BAM/CRAM file names
  --sample_ids SAMPLE_IDS [SAMPLE_IDS ...]
                        Space-delimited list of sample IDs, must have same
                        order as BAM/CRAM file names. BAM RG tag required if
                        this is ommitted.
  --command_file COMMAND_FILE
                        store commands in this file.
  --format FORMAT       comma separated list of FORMAT fields to include in
                        sample plot title
  --gff GFF             genomic regions (.gff with .tbi in same directory)
                        used when building HTML table and table filters
  --downsample DOWNSAMPLE
                        Number of normal reads/pairs to plot
  --manual_run          don't auto-run the samplot plot commands (command_file
                        will be deleted)
  ```
</details>

`samplot vcf` can be used to quickly apply some basic filters to variants. Filters are applied via the `--filter` argument, which may be repeated as many times as desired. Each expression specified with the `--filter` option is applied separately in an OR fashion, which `&` characters may be used within a statement for AND operations.

### Example:
```
samplot vcf \
    --filter "SVTYPE == 'DEL' & SU >= 8" \
    --filter "SVTYPE == 'INV' & SU >= 5" \
    --vcf example.vcf\
    -d test/\
    -O png\
    --important_regions regions.bed\
    -b example.bam > samplot_commands.sh
```
This example will create a directory named test (in the current working directory). A file named `index.html` will be created inside that directory to explore the images created.

**Filters:** The above filters will remove all samples/variants from output except:
* `DUP` variants with at least `SU` of 8
* `INV` variants with `SU` of at least 5

The specific `FORMAT` fields available in your VCF file may be different. I recommend SV VCF annotation with [duphold](https://github.com/brentp/duphold) by [brentp](https://github.com/brentp).

For more complex expression-based VCF filtering, try brentp's [slivar](https://github.com/brentp/slivar), which provides similar but more broad options for filter expressions.

**Region restriction.** Variants can also be filtered by overlap with a set of region (for example, gene coordinates for genes correlated with a disease). The `important_regions` argument provides a BED file of such regions for this example.

**Filtering for de novo SVs**
Using a [PED](https://gatkforums.broadinstitute.org/gatk/discussion/7696/pedigree-ped-files) file with `samplot vcf` allows filtering for variants that may be spontaneous/de novo variants. This filter is a simple Mendelian violation test. If a sample 1) has valid parent IDs in the PED file, 2) has a non-homref genotype (1/0, 0/1, or 1/1 in VCF), 3) passes filters, and 4) both parents have homref genotypes (0/0 in VCF), the sample may have a de novo variant. Filter parameters are not applied to the parents. The sample is plotted along with both parents, which are labeled as father and mother in the image.

Example call with the addition of a PED file:

<pre>
samplot vcf \
    --filter "SVTYPE == 'DEL' & SU >= 8" \
    --filter "SVTYPE == 'INV' & SU >= 5" \
    --vcf example.vcf\
    -d test/\
    -O png\
    <b>--ped family.ped\</b>
    --important_regions regions.bed\
    -b example.bam > samplot_commands.sh
</pre>

**Additional notes.**
* Variants where fewer than 95% of samples have a call (whether reference or alternate) will be excluded by default. This can be altered via the command-line argument `min_call_rate`.
* If you're primarily interested in rare variants, you can use the `max_hets` filter to remove variants that appear in more than `max_hets` samples.
* Large variants can now be plotted easily by samplot through use of `samplot plot`'s `zoom` argument. However, you can still choose to only plot variants larger than a given size using the `max_mb` argument. The `zoom` argument takes an integer parameter and shows only the intervals within +/- that parameter on either side of the breakpoints. A dotted line connects the ends of the variant call bar at the top of the window, showing that the region between breakpoint intervals is not shown.
* By default, if fewer than 6 samples have a variant and additional homref samples are given, control samples will be added from the homref group to reach a total of 6 samples in the plot. This number may be altered using the `min_entries` argument.
* Arguments that are optional in `samplot plot` can by given as arguments to `samplot vcf`. They will be applied to each image generated.


#### CRAM inputs
Samplot also support CRAM input, which requires a reference fasta file for
reading as noted above. Notice that the reference file is not included in this
repository due to size. This time we'll plot an interesting duplication at
X:101055330-101067156.

```
samplot plot \
    -n NA12878 NA12889 NA12890 \
    -b samplot/test/data/NA12878_restricted.cram \
      samplot/test/data/NA12889_restricted.cram \
      samplot/test/data/NA12890_restricted.cram \
    -o cramX_101055330_101067156.png
    -c chrX \
    -s 101055330 \
    -e 101067156 \
    -t DUP \
    -r hg19.fa
```


The arguments used above are the same as those used for the basic use case, with the addition of the following:

`-r` The reference file used for reading CRAM files

#### Plotting without the SV
Samplot can also plot genomic regions that are unrelated to an SV. If you do
not pass the SV type option (`-t`) then the top SV bar will go away and only
the region that is given by `-c` `-s` and `-e` will be displayed.

#### Long read (Oxford nanopore and PacBio) and linked read support
Any alignment that is longer than 1000 bp is treated as a long read, and
the plot design will focus on aligned regions and gaps. Aligned regions are in orange, and gaps follow the same DEL/DUP/INV color code used for short reads. The height of the alignment is based on the size of its largest gap.

<img src="/doc/imgs/longread_del.png">

If the bam file has an MI tag, then the reads will be treated as linked reads.
The plots will be similar to short read plots, but all alignments with the same MI is plotted at the same height according to alignment with the largest gap in the group. A green line connects all alignments in a group.

<img src="/doc/imgs/linkedread_del.png">
