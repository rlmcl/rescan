/*
This file is part of REscan. (https://github.com/rlmcl/rescan)

Copyright 2018-2020 Russell McLaughlin

REscan is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

REscan is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with REscan.  If not, see <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include "errors.c"
#include "functions.c"
#include "help.c"

int main( int argc, char **argv )
{
	char	line[BUFFSIZE],		// lines of SAM-format input
			rdid[RDIDLEN],		// read ID
			lbm[10000][RDIDLEN],// read IDs for tracking left-hand badmapped (lbm) reads
			rname[RNAMELEN],	// RNAME field (ie chromosome name)
			prevrname[RNAMELEN],// RNAME of previous read
			userchr[RNAMELEN],	// User-defined chromosome for outputting
			regionsfile[FNLEN],	// Regions file name
			id[SIDLEN] = "NA",	// Sample ID
			*rnext,				// RNEXT string (reference name for next segment)
			*linep,				// pointer for looping through line elements
			*token;				// for string spliting

	const char *delim ="\t";	// tab delimiter for line inputs
			
	int	start=-1,					// start position for reporting (specify 0 for whole-chromosome)
		end=-1,						// end position for reporting (specify chr length for whole-chromosome)
		userregionspecified = 0,	// to check that all three of -s, -e and -c are provided
		dist = 200,					// upstream/downstream parameter for searching each locus for bad read pairs
		minq = 20,					// minimum mapping quality to consider a good read
		maxfrag = 50000,			// maximum fragment length permitted (distance between read pairs)
		jump = 1,					// jump parameter for skipping bases in output
		maxpos = 0,					// maximum position reached in bam
		minpos = BASELEN,			// minimum position reached in bam
		field,						// field counter to identify right columns from bam
		flag,						// bitwise flag
		pos,						// read position
		mapq,						// mapping quality of read
		pnext,						// pos of next segment
		tlen,						// template length
		startpos,					// start position of read's influence
		endpos,						// end position of read's influence
		direction,					// direction of read's influence
		seqlen,						// length of sequence
		lowestlbm = 0,				// lowest populated slot in left-hand badmapped (lbm) array
		sizeoflbm = 0,				// maximum populated slot in lbm array
		i,							// counter
		option_index=0,				// option counter for getopt
		c;							// character initializer for getopt

	unsigned short *goodmapped,	// goodmapped array
				   *badmapped;	// badmapped array
	goodmapped = (unsigned short *)malloc( BASELEN * sizeof(unsigned short) );
	badmapped  = (unsigned short *)malloc( BASELEN * sizeof(unsigned short) );

	//Get user-defined parameters
	static struct option long_options[] =
	{
		{"add", 1, 0, 0},
		{"start",required_argument,0,'s'},		// start position for reporting rescan statistics
		{"end",required_argument,0,'e'},		// end position for reporting
		{"jump",required_argument,0,'j'},		// jump parameter for skipping bases of output
		{"distance",required_argument,0,'d'},	// up/downstream distance for searching
		{"maxfrag",required_argument,0,'m'},	// maximum fragment length
		{"minq",required_argument,0,'q'},		// minimum mapping quality for good reads
		{"regions",required_argument,0,'r'},	// regions file (bed format)
		{"id",required_argument,0,'i'},			// sample ID to report in VCF
		{"help",no_argument,0,'h'},				// ask for help
		{NULL, 0, NULL, 0}
	};
	
	while ((c = getopt_long(argc, argv, "hs:e:j:d:l:q:r:i:c:", long_options, &option_index)) != -1)
	{
		int this_option_optind = optind ? optind : 1;
		switch (c)
		{
			case 's': start = atoi(optarg); userregionspecified +=1; break;
			case 'e': end = atoi(optarg); userregionspecified +=1; break;
			case 'j': jump = atoi(optarg); break;
			case 'd': dist = atoi(optarg); break;
			case 'm': maxfrag = atoi(optarg); break;
			case 'q': minq = atoi(optarg); break;
			case 'c': strcpy(userchr,optarg); userregionspecified +=1; break;
			case 'r': strcpy(regionsfile,optarg); regionsarespecified = true ; break;
			case 'i': strcpy(id,optarg); break;
			case 'h': usage(regionsfile,id,userchr,start,end,jump,dist,maxfrag,minq); return 0;
		}
	}

	// Has user specified regions file or used -s, -e, -c? If so, parse into linked lists
	checkuserregions(userregionspecified,regionsfile,userchr,start,end);

	// Main part: read stdin (SAM stream) and get to work on the data
	while ( fgets ( line, sizeof line, stdin ) != NULL )
	{
		// initialize useful values
		field = 0;			// for extracting SAM fields
		linep = line;		// pointer for looping through line elements

		if( line[0] == '@' ) { continue; }	// skip header lines, if present

		// extract SAM fields
		while(token = strsep(&linep, delim))
		{
			switch(field)
			{
				case 0 : strcpy( rdid, token ); break;	// rdid now contains read ID
				case 1 : flag = atoi(token); break;		// flag now contains bitwise SAM flag
				case 2 : strcpy( rname, token ); break; // rname now contains reference seq name (chromosome)
				case 3 : pos = atoi(token); break;		// pos now contains mapping position of read
				case 4 : mapq = atoi(token); break;		// mapq now contains mapping quality of read
				case 6 : rnext = token; break;			// rnext now contains chr of next read in pair
				case 7 : pnext = atoi(token); break;	// pnext now contains pos of next read in pair
				case 8 : tlen = atoi(token); break;		// tlen now contains template length
				case 9 : seqlen = strlen(token); break;	// seqlen now contains length of read
			}
			field++;
		}

		// Has chromosome changed? If so, print results for this chromosome
		if( (strcmp(rname, prevrname)) )
		{
			report(prevrname, goodmapped, badmapped, jump, id, argc, argv);
			if( !regionsarespecified )		// No user regions specified; default to min and max bam positions
			{
				resetchr(rname);
			}
		}
		strcpy( prevrname, rname );

		// Determine read orientation
		switch( ( flag & 16 ) )
		{
			case 16 : // segment is reverse complemented
				startpos = pos - dist;
				endpos = pos + seqlen;
				direction = -1;
				break;
			case 0 : // segment is NOT reverse complemented
				startpos = pos;
				endpos = pos + seqlen + dist;
				direction = 1;
				break;
		}
        startpos = ( startpos < 0 ) ? 0 : startpos; // reassign negative startpos (happens at the start of chromosomes)

		if( !regionsarespecified )		// No user regions specified; update min and max bam positions
		{
			uehead->val = ( endpos   > uehead->val ) ?   endpos : uehead->val;	// reassign max pos if current bam position is greater
			ushead->val = ( startpos < ushead->val ) ? startpos : ushead->val;	// reassign min pos if current bam position is lower
		}

		// skip if fragment doesn't overlap region of interest
		if( regionsarespecified && !readisinregions( rname, startpos, endpos ) ) { continue; }

		// Main substance of programme: check if mate is well-aligned
		if( !( flag & 4 ) ) // segment is mapped adequately
		{
			if( mapq >= minq ) // segment is mapped with sufficient quality
			{
				if( !( flag & 8 ) &&							// next segment is mapped adequately
				    ( rnext[0] == '=' )	&&						// next segment is on same chromosome
				    ( ( pnext - pos ) * direction < maxfrag ) )	// next segment is within maxfrag distance
				{
					if( direction == 1 ) // First in pair; increment goodmapped
					{
						increment( goodmapped, startpos, endpos );
					}
					else // 2nd in pair; check mapping qualities
					{
						if( mapq < minq ) // 2nd read is actually bad (didn't know earlier)
						{
							increment( badmapped, startpos, endpos );
							decrement( goodmapped, startpos, endpos );	// molecule was erroneously counted earlier
						}
						else // 2nd read is ok; check lbm for 1st read
						{
							if( checklbm(rdid) )
							{
								increment( badmapped, startpos, endpos );
							}
							else
							{
								// do nothing: this is a fragment with both ends mapped well so was
								// already handled when 1st read was traversed (goodmapped was incremented)
							}
						}
					}
				}
				else // next segment is not mapped adequately, on a different chr or too far away on same chr
				{
					increment( badmapped, startpos, endpos );
				}
			}
			else // segment is not mapped with sufficient quality
			{
				if( direction == 1 ) // First in pair; add read id to lbm (if 2nd: don't need to do anything)
				{
					addtolbm( rdid );
				}
				else // 2nd in pair; just need to check if read id is also in lbm and delete it
				{
					checklbm( rdid );
				}
			}
		}
		else
		{
		}
	}
	// Report final statistics (after last chromosome)
	report(rname,goodmapped,badmapped,jump,id,argc,argv);
	return 0;
}
