/*
   print_branches.h
*/

#ifndef PRINT_BRANCHES_H_
#define PRINT_BRANCHES_H_

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <limits.h>
#include <string.h>

// path headers

#include "element.h"
#include "linkage_element.h"
#include "dB_graph.h"
#include "dB_graph_population.h"
#include "event_encoding.h"
#include "open_hash/little_hash_for_linkage.h"
#include "graph_info.h"

int db_graph_get_sequence_for_a_single_branch_for_specific_person_or_pop(dBNode* node, Orientation orientation, dBGraph *db_graph,
									 char* seq, boolean* flag, int index);
void db_graph_print_all_branches_for_specific_person_or_pop(char* out_file_name, dBGraph* db_graph, int index);
boolean db_graph_insert_branching_node_for_specific_person_or_pop(dBNode* node,dBGraph *db_graph,LinkHashTable *link_list,int index);
void db_graph_create_little_hash_table_for_branching_nodes_for_specific_person_or_pop(dBGraph* db_graph, LinkHashTable *link_list, int index);


#endif /* PRINT_BRANCHES_H_ */
