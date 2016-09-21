#define BUFFSIZE 4096

// define associative array structure ("hash")
typedef struct
{
	int size;
	void **keys;
	void **values;
} hash;

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


