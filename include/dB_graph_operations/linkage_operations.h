/*
   linkage_operations.h
*/

#ifndef LINKAGE_OPERATIONS_H_
#define LINKAGE_OPERATIONS_H_

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <limits.h>
#include <string.h>
#include <ctype.h>

// link_operations headers
#include "element.h"
#include "linkage_element.h"
#include "dB_graph.h"
#include "dB_graph_population.h"
#include "event_encoding.h"
#include "open_hash/little_hash_for_linkage.h"


struct passed_node
{
  unsigned short dist;
  dBNode* ptr;
  struct passed_node* next;
};

typedef struct passed_node PrevBranch;

void linkage_element_link_node_coverage_update (LinkNode* l);
void linkage_element_add_a_link(LinkageElement* e1, LinkageElement* e2, dBGraph* db_graph);
boolean if_two_branches_are_connected(dBNode* candidate, dBNode** branch_list, unsigned short* distance, unsigned short* distal_length, 
				      int read_length, int insert_size, int min_sup, dBGraph* db_graph, LinkHashTable* link_list);
void assemble_contigs_from_a_single_source(dBNode* source, Orientation initial_orientation, int read_length, int insert_length, int min_contig_length, int min_sup,
					   dBGraph* db_graph, LinkHashTable* link_list, unsigned int* count, FILE* output, int index);
void assemble_contigs_from_a_single_source_2(dBNode* source, Orientation initial_orientation, int read_length, int insert_length,int min_contig_length, int min_sup,
					     dBGraph* db_graph, LinkHashTable* link_list, unsigned int* count, FILE* output, int index);


#endif /* LINKAGE_OPERATIONS_H_ */
