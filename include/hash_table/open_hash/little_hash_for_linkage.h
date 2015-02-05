/*
 * 
 * CORTEX project contacts:  
 * 		M. Caccamo (mario.caccamo@bbsrc.ac.uk) and 
 * 		Z. Iqbal (zam@well.ox.ac.uk)
 *
 * **********************************************************************
 *
 * This file is part of CORTEX.
 *
 * CORTEX is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CORTEX is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CORTEX.  If not, see <http://www.gnu.org/licenses/>.
 *
 * **********************************************************************
 */
/*
  little_hash_for_linkage.h
*/

#ifndef LINK_HASH_H_
#define LINK_HASH_H_

#include "global.h"
#include "linkage_element.h"
#include "element.h"

typedef struct
{
  short kmer_size;
  long long number_buckets;
  int bucket_size;
  LinkageElement * table; 
  short * next_element; //keeps index of the next free element in bucket 
  long long * collisions;
  long long unique_kmers;
  int max_rehash_tries;
} LinkHashTable;


LinkHashTable * link_hash_table_new(int number_bits, int bucket_size, int max_rehash_tries, short kmer_size);

void link_hash_table_free(LinkHashTable * * little_hash_table);

//if the key is present applies f otherwise adds a new element for kmer
boolean link_hash_table_apply_or_insert(Key key, void (*f)(LinkageElement*), LinkHashTable *);

//applies f to every element of the table
void link_hash_table_traverse(void (*f)(LinkageElement *),LinkHashTable *);

//if the element is not in table create an element with key and adds it
LinkageElement * link_hash_table_find_or_insert(Key key, boolean * found, LinkHashTable * little_hash_table);
LinkageElement * link_hash_table_insert(Key key, LinkHashTable * little_hash_table);

void link_hash_table_print_stats(LinkHashTable *);

long long link_hash_table_get_unique_kmers(LinkHashTable *);

//return entry for kmer
LinkageElement * link_hash_table_find(Key key, LinkHashTable * little_hash_table);

#endif /* LINK_HASH_H_ */
