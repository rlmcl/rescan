#define BUFFSIZE 4096

// define associative array structure ("hash")
typedef struct
{
	int size;
	void **keys;
	void **values;
} hash;

// 
typedef struct
{
	hash	*chrlen,
		*offset,
		*linelength,
		*binlinelength;
} fai;

// function to declare new hash
hash *hash_new (int size)
{
	hash *h = calloc(1, sizeof (hash));
	h->keys = calloc(size, sizeof (void *));
	h->values = calloc(size, sizeof (void *));
	h->size = size;
	return h;
}

// function to get hash index (used in insert and lookup functions)
int hash_index (hash *h, void *key)
{
	int i = (int) key % h->size;
	while (h->keys[i] && h->keys[i] != key)
		i = (i + 1) % h->size;
	return i;
}

// function to insert new key-value pair into hash
void hash_insert (hash *h, void *key, void *value)
{
	int i = hash_index(h, key);
	h->keys[i] = key;
	h->values[i] = value;
}

// function to get element from hash given key
void *hash_lookup (hash *h, void *key)
{
	int i = hash_index(h, key);
	return h->values[i];
}

// function to increment goodmapped or badmapped hashes
void increment (hash *h, int startpos, int endpos)
{
	int i, hashval;
	for( i=startpos; i<=endpos; i++ )
	{
		hashval = hash_lookup( h, i );
		if( hashval == NULL ) { hash_insert( h, i, 1 ); }
		else { hash_insert( h, i, (hashval+1) ); }
	}
}

// function to load fai reference
void load_fai (char *refgenome)
{
	char 	line[128],
		*linep = line,
		*token,
		*sn;
	const char delim = "\t";
	int 	field = 0,
		chrlen,
		offset,
		linelength,
		binlinelength;
	FILE 	*faifile = fopen( refgenome, "r" );
//	fai 	f;
	
//	f.chrlen = hash_new(10);
//	f.offset = hash_new(10);
//	f.linelength = hash_new(10);
//	f.binlinelength = hash_new(10);
	fprintf( stdout, "Hello from the function\n" );	
	while( fgets( line, sizeof line, faifile ) != NULL )
	{
		fprintf( stdout,"%s\n",line );
		while(token = strsep(&linep, delim))
		{
			switch(field)
			{
				case 0 : sn = token; break;
			//	case 1 : hash_insert( f.chrlen, sn, atoi(token) ); break;
			//	case 2 : hash_insert( f.offset, sn, atoi(token) ); break;
			//	case 3 : hash_insert( f.linelength, sn, atoi(token) ); break;
			//	case 4 : hash_insert( f.binlinelength, sn, atoi(token) ); break;
			}
			field++;
		}
	}
//	return f;
}



