# rescan

## Quick start
### Clone repository
> git clone git clone https://github.com/rlmcl/rescan

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

