void usage()
{
	fprintf(stderr,	"REscan alpha version 0.00\n"
			"Russell McLaughlin\n"
			"------------------\n"
			"Usage:\tsamtools view in.bam [ region ] | rescan [ options ]\n\n"
			"Options:\n"
			"       --start (-s)  INT  : start position for reporting rescan statistics\n"
			"         --end (-e)  INT  : end position for reporting rescan statistics\n"
			"    --distance (-d)  INT  : up/downstream distance for searching\n"
			" --levenshtein (-l) FLOAT : edit distance above which to consider a read bad\n"
			"        --minq (-q)  INT  : minimum mapping quality for good reads\n"
			"        --help (-h)       : print this help message\n"
);
return;
}


