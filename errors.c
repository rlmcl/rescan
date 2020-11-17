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
