void usage(	char *regionsfile,
			char *id,
			char *userchr,
			int start,
			int end,
			int jump,
			int dist,
			int maxfrag,
			int minq)
{
	char endhelp[20];
	if( regionsfile[0] == NULL )
	{
		strcpy( regionsfile, "currently unspecified" );
	}
	if( end == -1 )
	{
		strcpy(endhelp, "max position in bam" );
	}
	else
	{
		sprintf( endhelp,"%d",end );
	}
	fprintf(stderr,	
			"------------------------------------------\n"
			"REscan version 1.0.0\n"
			"Russell McLaughlin, Trinity College Dublin\n"
			"GNU General Public License v3\n"
			"------------------------------------------\n"
			"Usage:\tsamtools view in.bam [ region ] | rescan [ options ]\n\n"
			"Options:\n"
			"   --regions (-r)   FILE  : file name for bed-format, position-sorted regions (%s)\n"
			"        --id (-i) STRING  : sample ID (%s)\n"
			"       --chr (-c) STRING  : chromosome for reporting rescan statistics (%s)\n"
			"     --start (-s)    INT  : start position for reporting rescan statistics (%d)\n"
			"       --end (-e)    INT  : end position for reporting rescan statistics (%s)\n"
			"      --jump (-j)    INT  : number of bases to jump by in printing output (%d)\n"
			"  --distance (-d)    INT  : up/downstream distance for searching (%d)\n"
			"   --maxfrag (-m)    INT  : maximum fragment length allowed (%d)\n"
			"      --minq (-q)    INT  : minimum mapping quality for good reads (%d)\n"
			"      --help (-h)         : print this help message\n"
			,regionsfile,id,userchr,start,endhelp,jump,dist,maxfrag,minq);
return;
}


