/*
 * Copyright 2009-2011 Zamin Iqbal and Mario Caccamo
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
  linkage_element.h -  used when genotyping. one colour per allele, plus one
  for reference plus two spare for working
*/

#ifndef LINKAGE_ELEMENT_H_
#define LINKAGE_ELEMENT_H_

#include <stdio.h>
#include <inttypes.h>

#include "global.h"
#include "element.h" // need this for Edges, etc
#include "binary_kmer.h"

//Bubble and PD callers only call biallelic for now
#define MAX_ALLELES_SUPPORTED_FOR_STANDARD_GENOTYPING 2


//type definitions
typedef enum
  {
    one_read   = 0,
    two_read   = 1,
    three_read = 2,
    more_read  = 3
  } LinkStatus;

typedef enum
  {
    just_initiated = 0,
    pvalue_passed  = 1
  } LinkageNodeStatus;

struct linked_node
{
  dBNode* ptr;
  LinkStatus read;
  struct linked_node* next;
};

typedef struct linked_node LinkNode;

typedef struct
{
  unsigned short link_num;
  BinaryKmer kmer;
  LinkageNodeStatus status;
  LinkNode* link;
} LinkageElement;


LinkageElement* new_linkage_element();
void free_linkage_element(LinkageElement** linkage_element);
void linkage_element_initialise_from_normal_element(LinkageElement* e1, 
						    Element* e2);
void linkage_element_assign(LinkageElement* e1, LinkageElement* e2);
boolean db_linkage_node_check_for_flag_ALL_OFF(LinkageElement * node);
BinaryKmer* linkage_element_get_kmer(LinkageElement *);
boolean linkage_element_is_key(Key,LinkageElement);
Key linkage_element_get_key(BinaryKmer*,short kmer_size, Key preallocated_key);
void linkage_element_initialise(LinkageElement *,Key, short kmer_size);
void linkage_element_initialise_kmer_to_zero(LinkageElement * e);
void linkage_element_set_kmer(LinkageElement * e, Key kmer, short kmer_size);
Orientation db_linkage_node_get_orientation(BinaryKmer*, LinkageElement *, short kmer_size);

#endif /* LINKAGE_ELEMENT_H_ */
