# rescan

## Quick start
### Clone repository
```
git clone https://github.com/rlmcl/rescan
```

### Compile
```
cd rescan
gcc -w -o rescan rescan.c
```

### Install
```
sudo cp rescan /usr/local/bin
```

### Run on SAM-format stream
```
samtools view in.bam [ region ] | rescan [ options ] 
```

### Options
| Argument 				| Type	| Description										|
| --------------------:	| ----: | ------------------------------------------------- |
| `--regions` (`-r`)	| FILE  | file name for bed-format, position-sorted regions |
| `--id` (`-i`)			| STRING| sample ID 										|
| `--chr` (`-c`)		| STRING| chromosome for reporting rescan statistics 		|
| `--start` (`-s`)		| INT   | start position for reporting rescan statistics	|
| `--end` (`-e`)		| INT   | end position for reporting rescan statistics		|
| `--jump` (`-j`)		| INT   | number of bases to jump by in printing output 	|	
| `--distance` (`-d`)	| INT 	| up/downstream distance for searching 				|
| `--maxfrag` (`-m`)	| INT 	| maximum fragment length allowed 					|
| `--minq` (`-q`)		| INT 	| minimum mapping quality for good reads 			|
| `--help` (`-h`)		| -		| get help			 								|

### Regions
As well as controlling data flow into REscan using [SAMtools](http://www.htslib.org/doc/samtools.html) (see *Feed it only the data it needs* under *Caveats and features*, below), you have three options for specifying regions for REscan output:

1. **Don't specify regions**

   REscan will output statistics for every locus encountered, skipping over bases if `-j` (`--jump`) is specified. Output data will correspond to every chromosome encountered in RNAME in the input stream, starting at the lowest position encountered for each chromosome (potentially minus distance \[`-d`/`--distance`\] parameter if read is reverse orientation) and finishing at the highest position for that chromosome.
   
2. **Specify a single region**

   `-c`, `-s` and `-e` (`--chr`, `--start` and `--end`) can be used to instruct REscan only to output data for a single region.
   
3. **Specify multiple regions**

   Using `-r` (`--regions`), a [BED-format](https://en.wikipedia.org/wiki/BED_(file_format)) regions file can be specified. Expected format is tab-delimited, with chromosome, start and end positions and an optional name for each region (eg gene/transcript name). Regions should be position-sorted, with chromosomes in the same order as they appear in the SAM data. Example:  
   `chr6	16299112	16761490	ATXN1`  
   `chr9	27546545	27573866	C9orf72`  
   `chr12	111452268	111599676	ATXN2`  

### Examples

Report REscan statistics for all of chromosome 9 from `input.bam` and store in `output.vcf`:
```
samtools view input.bam chr9 | rescan > output.vcf
```
Stream chromosome 9 data from `input.bam` and report REscan statistics just for the region chr9:27570000-27577000. Specify larger fragment length than default (`-d`):
```
samtools view input.bam chr9 | rescan -c chr9 -s 27570000 -e 27577000 -d 600 > output.vcf
```
Stream from a genome file called `jim_genome.bam` and report REscan statistics for all transcripts with genomic coordinates specified in `transcripts.bed`. Specify sample ID as `Jim`. Only report every 100th base (`-j`). Compress the output with `bgzip`:

```
samtools view jim_genome.bam | rescan -r transcripts.bed -i Jim -j 100 | bgzip > jim_rescan.vcf.gz
```

## About REscan
REscan is a simple tool for counting the number of poorly-paired reads spanning a region, possibly reflecting the presence of a repeat expansion. Results are reported as the fraction of reads with unmapped or distantly-mapped mates. Output is VCF-format with a field `RS`, representing the REscan statistic r<sub>x</sub>/r<sub>t</sub>, where r<sub>x</sub> is the number of poorly-mapped reads pointing towards the locus (represented as `BM` or "badmapped" in the VCF output) and r<sub>t</sub> is the total number of (nearby) reads pointing into the locus (sum of `BM` and `GM` or "goodmapped" in the VCF output).

### Caveats and features

- **REscan reads from stdin**

   This is for simplicity and flexibility, so it can be built into pipelines (eg using `samtools view`), as detailed above.

- **REscan writes to stdout**

   Again for flexibility; it is recommended you pipe into [`bgzip`](http://www.htslib.org/doc/bgzip.html) (block gzip tool that comes with [`tabix`](http://www.htslib.org/doc/tabix.html) -- `sudo apt install tabix`) so that you can index and fast-access regions later (`tabix output.vcf.gz`).

- **REscan expects paired-end, _position-sorted_ data**

   The fundamentals of the REscan statistic rely on mapped reads with poorly-mapped mates, so you should use paired-end sequence data. It won't do anything useful with single-end data.

- **REscan expects only one sample per run**

   You'll get funky results if you stream a merged BAM file. You can set sample ID using the `--id` (`-i`) argument (defaults to `NA`).
   
- **Feed it only the data it needs**

   If you are only generating statistics for certain regions (eg specified using `-r`/`--regions`) then consider also only inputting data relevant to these regions from SAMtools. Bear in mind that for the edges of your regions you will also need SAM data +/- some distance (eg, say, 500 bases) so that you will accurately count reads orientated into your loci of interest. One sensible solution is `samtools view in.bam -L regions1.bed | rescan -r regions2.bed [...]`, where `regions2.bed` contains your regions of interest and `regions1.bed` contains those regions extended by +/- 500 bp.

- **REscan doesn't do variant calling _per se_**

   The REscan statistic is very basic; REscan will not measure the repeat length, for example. For a more comprehensive tool, try ExpansionHunter or one of the many _punSTR_ programmes.

- **It won't just identify repeat expansions**

   There are many ways to end up with poorly-mapped mates as defined by REscan (eg inversions, translocations, CNVs). This programme is called REscan because the original motivation to write it was to identify novel potential repeat expansions. A high statistic at a locus warrants further investigation and should not be considered as definitive evidence for a repeat expansion. You should follow up with other tools and wet-lab validation, if possible.

- **REscan assumes chromosomes are linear**

   It doesn't (currently) do anything fancy with circular chromosomes like the mitochondrial genome so statistics at SAM positions 0 or length(chr) won't be quite right.

- **If your data aren't human, you might need to tweak**

   The longest chromosome length the algorithm can handle is defined in `functions.c` in the `BASELEN` macro (268,435,456, ie 2^28, at time of writing). This is long enough to fit human chromosome 1; if you are working with an organism with longer chromosomes you will need to change this number and recompile.

- **Know your input data**

   It is useful to know a bit about the provenance of your data and how it looks before you use this tool. For example, if you can find it out, it is useful to know the mean fragment length of DNA molecules used in library preparation, and to set `--distance` (`-d`) accordingly. You could also, of course, look at the fragment size distribution in your SAM data in advance to inform this decision.

- **Know your output data**

   Generating these statistics is just the beginning of an exploratory journey; you should look at your data and have a think about how you might clean them up (eg standardising within individuals) and, for interesting loci, cross-reference with structural variant callers. If you still suspect a repeat expansion, have a look at the reference genome for potential motifs.

### License
GNU General Public license v3

### Citation
Pending
