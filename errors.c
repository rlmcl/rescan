#define	E0			"REscan error: "
#define BASELENERR	"out of memory. This is probably because there is a chromosome in your data larger than the current maximum. You can change this by editing the BASELEN definition at the top of functions.c and recompiling."
#define REGIONSERR	"if you use -s, -e or -c, you must specify all three"
#define REGSPECERR	"please check that your regions file is properly formatted"
#define LBMERR		"unspecified problem with reading read IDs (LBM error)"

void diewitherror( char *error )
{
	fprintf( stderr, E0 "%s\n",error );
	exit(1);
}
