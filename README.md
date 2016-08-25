# rescan

## Quick start
### Clone repository
> git clone https://github.com/rlmcl/rescan

### Compile
> cd rescan

> gcc -w -o rescan rescan.c

### Install
> sudo cp rescan /usr/local/bin

### Run on SAM-format stream
> samtools view in.bam [ region ] | rescan \\ 

>   [ -s start region for reporting ] \\

>   [ -e end region for reporting ] \\

>   [ -d distance parameter ] > out.rescan

## Overview
Rescan is a simple tool for counting the number of poorly-paired reads spanning a region, possibly reflecting the presence of a repeat expansion. Results are reported as the fraction of reads with unmapped or distantly-mapped mates.
