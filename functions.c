#define BUFFSIZE	4096		// Buffer for reading lines of input
#define BASELEN		268435456	// Max chromosome size (basepairs; 2^28 = 268,435,456)
#define	READTRK		10000		// Max number of leftmost poorly-mapping reads that can be tracked at once
#define RDIDLEN		255			// Max length of read IDs
#define SIDLEN		255			// Max length of sample IDs
#define RNAMELEN	64			// Max length of reference names (chromsomes)
#define REGNAMELEN	64			// Max length of user-defined region names
#define FNLEN		256			// Max length for regions file name
#define MAXLBM		512			// Max number of reads to hold in memory for leftbadmapped list
#define DEBUG(a) do { printf("%s\n", a); fflush(stdout); } while (0)


// Global bools to keep track of things
bool	thisis1stwrite		= true,	// if first time writing vcf data (print header)
		regionsarespecified	= false;// shorthand for control structures if user provides regions file

int     lbmcounter = 0;

// Linked list structures
struct userstarts 	{ int val;				struct userstarts *next; };
struct userends 	{ int val;				struct userends   *next; };
struct userchrs 	{ char val[RNAMELEN];	struct userchrs   *next; };
struct usernames 	{ char val[REGNAMELEN];	struct usernames  *next; };
struct lbm			{ char val[RDIDLEN];	struct lbm		  *next; };

// Global linked list pointers
//// user-defined regions
struct userstarts *ushead = NULL;	// userstarts head
struct userstarts *uscurr = NULL;	// userstarts current
struct userstarts *usrchk = NULL;	// userstarts region check (when scanning)
struct userstarts *ussrpt = NULL;	// userstarts statistics pointer (when reporting)
struct userends   *uehead = NULL;	// userends head
struct userends   *uecurr = NULL;	// userends current
struct userends   *uerchk = NULL;	// userends region check (when scanning)
struct userends   *uesrpt = NULL;	// userends statistics pointer (when reporting)
struct userchrs   *uchead = NULL;	// userchrs head
struct userchrs   *uccurr = NULL;	// userchrs current
struct userchrs   *ucrchk = NULL;	// userchrs region check (when scanning)
struct userchrs   *ucsrpt =	NULL;	// userchrs statistics pointer (when reporting)
struct usernames  *unhead = NULL;	// usernames head
struct usernames  *uncurr = NULL;	// usernames current
struct usernames  *unrchk = NULL;	// usernames region check (when scanning)
struct usernames  *unsrpt = NULL;	// usernames statistics pointer (when reporting)
//// lbm linked list (lbm = left badmapped)
struct lbm 		  *lbmhead=	NULL;	// left badmapped head
struct lbm 		  *lbmcurr=	NULL;	// left badmapped current


// Check if current read/fragment is in a user-specified region
bool readisinregions( char* rname, int startpos, int endpos)
{
	char prevchr[RNAMELEN];
	if( usrchk == NULL ) // first time checking
	{
		usrchk = ushead;
		uerchk = uehead;
		ucrchk = uchead;
		unrchk = unhead;
	}
	strcpy( prevchr,ucrchk->val );
	while(usrchk != NULL)
	{
		if( !strcmp(rname,ucrchk->val) )
		{
			// is entire fragment 5' of user region? Skip the fragment
			if( endpos<usrchk->val )
			{
				return 0;
			}
			// is fragment within user region?
			else if( startpos<=uerchk->val && endpos>=usrchk->val )
			{
				return 1;
			}
		}
		// Either current user chromosome is different to fragment or fragment is 3' of user region -- skip the region
		usrchk = usrchk->next;
		uerchk = uerchk->next;
		ucrchk = ucrchk->next;
		unrchk = unrchk->next;
	}
	return 0;
}

// write vcf header
void write_vcf_header(char *id, int argc, char **argv)
{
	time_t rawtime;
	struct tm * timeinfo;
	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	fprintf(stdout,	"##fileformat=VCFv4.1\n"
					"##FORMAT=<ID=RS,Number=1,Type=Float,Description=\"REscan statistic\">\n"
					"##FORMAT=<ID=GM,Number=1,Type=Integer,Description=\"Goodmapped reads (REscan): number of reads pointing into locus with well-mapped mate\">\n"
					"##FORMAT=<ID=BM,Number=1,Type=Integer,Description=\"Badmapped reads (REscan): number of reads pointing into locus with poorly-mapped mate\">\n"
					"##INFO=<ID=UserRegion,Number=1,Type=String,Description=\"User-supplied region name for current REscan statistic\">\n"
					"##starttime=%s"
					"##command=",asctime(timeinfo));
	int a;
	for(a = 0; a < argc; a++)
	{
		fprintf(stdout,"%s ",argv[a]);
	}
	fprintf(stdout,"\n");
	fprintf(stdout,"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n",id);
	fflush(stdout);
}

int report(char *rname, unsigned short *goodmapped, unsigned short *badmapped, int jump, char *id, int argc, char **argv)
{
	float result;
	if( thisis1stwrite )
	{
		write_vcf_header(id,argc,argv);
		ussrpt = ushead;
		uesrpt = uehead;
		ucsrpt = uchead;
		unsrpt = unhead;
		thisis1stwrite=false;
		return 0;
	}

	// loop over regions and report
	while(ussrpt != NULL)
	{
		if( !strcmp(rname,ucsrpt->val) )
		{
			// report stats from ussrpt->val to uesrpt->val
			int i;
			for( i = ussrpt->val; i < uesrpt->val ; i+=jump )
			{
				if(goodmapped[i]+badmapped[i]==0) { result = (float)0; }
				else { result = (float)badmapped[i]/((float)goodmapped[i]+(float)badmapped[i]); }
				fprintf( stdout,	"%s\t%d\t"					// Chromsome, position
									".\t.\t.\t.\t.\t"			// variant ID,ref,alt,qual,filter
									"UserRegion=%s\tRS:GM:BM\t"	// info,format
									"%f:%d:%d\n",				// REscan statistic,goodmapped,badmapped	
								ucsrpt->val,				// Chromosome
									i,							// Position
									unsrpt->val,				// User-specified region name
									result,						// REscan statistic
									goodmapped[i],				// goodmapped number
									badmapped[i]				// badmapped number
						);
			}

			ussrpt = ussrpt->next;
			uesrpt = uesrpt->next;
			ucsrpt = ucsrpt->next;
			unsrpt = unsrpt->next;
		}
		else
		{
			return 0;
		}
	}

	// reset goodmapped and badmapped
	memset(goodmapped,0,sizeof(goodmapped));
	memset(badmapped, 0,sizeof( badmapped));

	/////////////////////////////////////////////////////////
	// TODO: should also clear lbm, ie free all of its memory locations
	/////////////////////////////////////////////////////////
}


// add to goodmapped or badmapped array between start and end positions
void increment( unsigned short *a, int *startpos, int *endpos )
{
	int i;
	for( i=startpos; i<=endpos; i++ ) { a[i]++; }
}

// subtract from goodmapped or badmapped array between start and end positions
void decrement( unsigned short *a, int *startpos, int *endpos )
{
	int i;
	for( i=startpos; i<=endpos; i++ ) { a[i]--; }
}

struct lbm* addtolbm(char *val)
{
    lbmcounter++;

	bool create = false;	// toggle to create linked list if this is the first
	if(lbmhead == NULL) { create = true; }

	// if lbm list is bigger than 1,000, start culling the head
	struct lbm *del = (struct lbm*)malloc(sizeof(struct lbm));
	if( lbmcounter >= MAXLBM )
	{
		lbmcounter--;
		lbmhead = lbmhead->next;
		free(del);
//		del = NULL;
	}

	struct lbm *lbmptr = (struct lbm*)malloc(sizeof(struct lbm));
	if(lbmptr == NULL) { diewitherror(LBMERR); }
	strcpy( lbmptr->val, val );
	lbmptr->next = NULL;
	if( create )
	{
		lbmhead = lbmcurr = lbmptr;
		return lbmptr;
	}
	else
	{
		lbmcurr->next = lbmptr;	// add to end of linked list
		lbmcurr = lbmptr;		// not start
		return lbmptr;
	}
}


int checklbm( char *rdid )
{
	// find if read id exists in lbm
	// if so, delete it and return 1
	// else return 0
	struct lbm *current  = lbmhead;
	struct lbm *previous = NULL;

	while( current != NULL )
	{
		if( !strcmp( current->val, rdid ) )
		{
			// delete value from linked list
			if( previous != NULL )			{ previous->next = current->next; }
			else 							{ lbmhead = current->next; }
			if( current == lbmcurr )		{ lbmcurr = previous; }
			else if( current == lbmhead )	{ lbmhead = current->next; }
			free( current );
			current = NULL;

            lbmcounter--;

			// return success at finding rdid
			return 1;
		}
		else
		{
//			fprintf(stdout,"Scenario 2\n");
			previous = current;
			current  = current->next;
		}
	}
	return 0;
}


/// LINKED LISTS
// with help from:
// https://www.thegeekstuff.com/2012/08/c-linked-list-example/

struct userstarts* append_userstarts(int val)
{
	bool create = false;	// toggle to create linked list if this is the first
	if(ushead == NULL) { create = true; }

	struct userstarts *usptr = (struct userstarts*)malloc(sizeof(struct userstarts));
	if(usptr == NULL) { diewitherror(REGSPECERR); }
	usptr->val = val;
	usptr->next = NULL;
	if( create )
	{
		ushead = uscurr = usptr;
		return usptr;
	}
	else
	{
		uscurr->next = usptr;
		uscurr = usptr;
		return usptr;
	}
}

struct userends* append_userends(int val)
{
	bool create = false;	// toggle to create linked list if this is the first
	if(uehead == NULL) { create = true; }

	struct userends *ueptr = (struct userends*)malloc(sizeof(struct userends));
	if(ueptr == NULL) { diewitherror(REGSPECERR); }
	ueptr->val = val;
	ueptr->next = NULL;
	if( create )
	{
		uehead = uecurr = ueptr;
		return ueptr;
	}
	else
	{
		uecurr->next = ueptr;
		uecurr = ueptr;
		return ueptr;
	}
}

struct userchrs* append_userchrs(char val[RNAMELEN])
{
	bool create = false;	// toggle to create linked list if this is the first
	if(uchead == NULL) { create = true; }

	struct userchrs *ucptr = (struct userchrs*)malloc(sizeof(struct userchrs));
	if(ucptr == NULL) { diewitherror(REGSPECERR); }
	strcpy( ucptr->val, val );
	ucptr->next = NULL;
	if( create )
	{
		uchead = uccurr = ucptr;
		return ucptr;
	}
	else
	{
		uccurr->next = ucptr;
		uccurr = ucptr;
		return ucptr;
	}
}

struct usernames* append_usernames(char val[REGNAMELEN])
{
	bool create = false;	// toggle to create linked list if this is the first
	if(unhead == NULL) { create = true; }

	struct usernames *unptr = (struct usernames*)malloc(sizeof(struct usernames));
	if(unptr == NULL) { diewitherror(REGSPECERR); }
	if((val[0] == '\0') || (val[0] == '\n')) { strcpy(val,"region_unnamed\0"); } // Check for empty user region name
	else if(val[strlen(val)-1]=='\n') {val[strlen(val)-1] = '\0';} // If not empty, strip newline from user region name
	strcpy( unptr->val, val );
	unptr->next = NULL;
	if( create )
	{
		unhead = uncurr = unptr;
		return unptr;
	}
	else
	{
		uncurr->next = unptr;
		uncurr = unptr;
		return unptr;
	}
}

// Regions list builder
void buildregionslist( char *regionsfile )
{
	FILE		*regions;			// filehandle for regions file
	char		line[BUFFSIZE],		// lines of regions file
				*linep,				// pointer for looping through line elements
				*token;				// for string spliting
	const char	*delim ="\t";	// tab delimiter for line inputs
	int			field;				// field counter to identify right columns from bam

	regions = fopen( regionsfile, "r" );
	while( fgets( line, sizeof line, regions ) != NULL )
	{
		field = 0;		// for extracting tab-delimited entries (chr start end ID)
		linep = line;	// pointer for looping through line elements
		while(token = strsep(&linep, delim))
		{
			switch( field )
			{
				case 0 : append_userchrs(token); break;				// adds chromosome to userchrs linked list
				case 1 : append_userstarts(atoi(token)); break;		// adds start pos to userstarts linked list
				case 2 : append_userends(atoi(token)); break;		// adds end pos to userends linked list
				case 3 : append_usernames(token); break;			// adds user region name to usernames linked list
			}
			field++;
		}
	}
	fclose( regions );
}

void buildregion( char *userchr, int start, int end )
{
	char temp[REGNAMELEN] = "Userregion\0";
	append_userchrs(userchr);
	append_userstarts(start);
	append_userends(end);
	append_usernames(temp);
}

void checkuserregions(int userregionspecified,char *regionsfile,char *userchr,int start,int end)
{
	// check for well-specified user input
	if( !regionsarespecified && userregionspecified>0 && userregionspecified<3 ) { diewitherror(REGIONSERR); }
	// regions file provided by user
	if( regionsarespecified ) 			{ buildregionslist( regionsfile ); } 
	// -s, -e and -c specified by user
	else if ( userregionspecified==3 ) 	{ buildregion( userchr,start,end ); regionsarespecified = true; } 
	// no regions specified -- start with absurd assignment that gets updated on the fly
	else 								{ buildregion( "",BASELEN,0 ); } 
}

void resetchr( char *rname )
{
	ushead->val = BASELEN;      // reset userstarts | note this is a hack-ey, probably unintuitive way of resetting
	uehead->val = 0;            // reset userends   | these values. If user regions aren't specified we just use the
	strcpy(uchead->val,rname);  // reset userchrs   | user regions linked lists anyway to store the chr and bam positions
	ussrpt = ushead;
	uesrpt = uehead;
	ucsrpt = uchead;
	unsrpt = unhead;
	usrchk = ushead;
	uerchk = uehead;
	ucrchk = uchead;
	unrchk = unhead;
}
