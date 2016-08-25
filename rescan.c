#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
//#include "samtools-1.3.1/bam.h"

#define BUFFSIZE 4096

typedef struct
{
	int size;
	void **keys;
	void **values;
} hash;
 
hash *hash_new (int size)
{
	hash *h = calloc(1, sizeof (hash));
	h->keys = calloc(size, sizeof (void *));
	h->values = calloc(size, sizeof (void *));
	h->size = size;
	return h;
}
 
int hash_index (hash *h, void *key)
{
	int i = (int) key % h->size;
	while (h->keys[i] && h->keys[i] != key)
		i = (i + 1) % h->size;
	return i;
}
 
void hash_insert (hash *h, void *key, void *value)
{
	int i = hash_index(h, key);
	h->keys[i] = key;
	h->values[i] = value;
}
 
void *hash_lookup (hash *h, void *key)
{
	int i = hash_index(h, key);
	return h->values[i];
}

void usage()
{
	fprintf(stderr,"REscan development version 0.00\nNo more help than that, I'm afraid\n"
);
return;
}

int main( int argc, char **argv )
{
	int	start=0,	// start position for reporting (specify 0 for whole-chromosome)
		end=250000000,	// end position for reporting (specify chr length for whole-chromosome)
		dist = 200,	// upstream/downstream parameter for searching each locus for bad read pairs
		maxfrag = 2000,	// maximum fragment length permitted (distance between read pairs)
		field,		// field counter to identify right columns from bam
		flag,		// bitwise flag
		pos,		// read position
		pnext,		// pos of next segment
		startpos,	// start position of read's influence
		endpos,		// end position of read's influence
		direction,	// direction of read's influence
		seqlen,		// length of sequence
		hashval,	// holder for current value of good/badmapped hash when looking up or incrementing
		hashval2,	// same as hashval
		i,		// counter
		option_index=0,	// option counter for getopt
		c;		// character initializer for getopt

	char	line[BUFFSIZE],	// lines of SAM-format input
		*rnext,		// RNEXT string (reference name for next segment)
		*token;		// for string spliting
		
	float	result;		// placeholder for results while printing
	
	const char *delim = "\t";// tab delimiter for line inputs

	//Get user-defined parameters
	static struct option long_options[] =
	{
		{"add", 1, 0, 0},
		{"start",required_argument,0,'s'},
		{"end",required_argument,0,'e'},
		{"distance",required_argument,0,'d'},
		{"help",no_argument,0,'h'},
		{NULL, 0, NULL, 0}
	};
	while ((c = getopt_long(argc, argv, "hs:e:d:", long_options, &option_index)) != -1)
	{
		int this_option_optind = optind ? optind : 1;
		switch (c)
		{
			case 'h': usage(); return 0;
			case 's': start = atoi(optarg); break;
			case 'e': end = atoi(optarg); break;
			case 'd': dist = atoi(optarg); break;
		}
	}
	hash 	*goodmapped = hash_new( 250000000 ),//end-start ),
		*badmapped  = hash_new( 250000000 );//end-start );
	

	while ( fgets ( line, sizeof line, stdin ) != NULL )
	{
		// initialize useful values
		field = 0;		// for extracting SAM fieils
		char *linep = line;	// pointer for looping through line elements
		
		while(token = strsep(&linep, delim))
		{
			switch(field)
			{
				case 1 : flag = atoi(token); break;
				case 3 : pos = atoi(token); break;
				case 6 : rnext = token; break;
				case 7 : pnext = atoi(token); break;
				case 9 : seqlen = strlen(token); break;

			}
			field++;
		}
		if( ( flag & 16 ) == 16 ) // segment is reverse complemented
		{
			startpos = pos - dist;
			endpos = pos + seqlen;
			direction = -1;
		}
		else // segment is NOT reverse complemented
		{
			startpos = pos;
			endpos = pos + seqlen + dist;
			direction = 1;
		}
		if( ( flag & 4 ) == 0 ) // segment is mapped adequately
		{
			if( ( flag & 8 ) == 0 ) // next segment is mapped adequately
			{
				if( !strcmp(rnext,"=") ) // next segment is on same chromosome
				{
					if( ( pnext - pos ) * direction < maxfrag ) // next segment is within 2kb
					{
						for( i=startpos; i<=endpos; i++ )
						{
							hashval = hash_lookup( goodmapped, i );
							if( hashval == NULL ) { hash_insert( goodmapped, i, 1 ); }
							else { hash_insert( goodmapped, i, (hashval+1) ); }
						}
					}
					else // next segment is on same chr but too far away (probably rare)
					{
						for( i = startpos; i <= endpos; i++ )
						{
							hashval = hash_lookup( badmapped, i );
							if( hashval == NULL ) { hash_insert( badmapped, i, 1 ); }
							else { hash_insert( badmapped, i, (hashval+1) ); }
						}
					}
				}
				else // next segment is on a different chromosome
				{
					for( i = startpos; i <= endpos; i++ )
					{
						hashval = hash_lookup( badmapped, i );
						if( hashval == NULL ) { hash_insert( badmapped, i, 1 ); }
						else { hash_insert( badmapped, i, (hashval+1) ); }
					}
				}
			}
			else // next segment is not mapped adequately
			{
				for( i = startpos; i <= endpos; i++ )
				{
					hashval = hash_lookup( badmapped, i );
					if( hashval == NULL ) { hash_insert( badmapped, i, 1 ); }
					else { hash_insert( badmapped, i, (hashval+1) ); }
				}
			}
		}
	}
	for( i = start; i < end; i++ )
	{
		hashval = hash_lookup( badmapped, i );
		hashval2= hash_lookup( goodmapped,i );
		if( hashval == NULL ) { hashval = 0; }
		if( hashval2== NULL ) { hashval2= 0; }
		if( hashval > 0 )
		{
			result = (float)hashval/((float)hashval+(float)hashval2);
			fprintf( stdout, "%f\t", result );
		}
		else
		{
			fprintf( stdout, "0\t" );
		}
	}
	fprintf( stdout, "\n" );
	return 0;
}
