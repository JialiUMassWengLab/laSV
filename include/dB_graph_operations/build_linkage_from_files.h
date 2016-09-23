/*
   build_linkage_from_files.h
*/

#ifndef BUILD_LINKAGE_FROM_FILES_H_
#define BUILD_LINKAGE_FROM_FILES_H_

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <limits.h>
#include <string.h>
#include <ctype.h>

//third party headers
#include <seq_file.h>
#include <string_buffer.h>

// link_operations headers
#include "element.h"
#include "linkage_element.h"
#include "dB_graph.h"
#include "dB_graph_population.h"
#include "event_encoding.h"
#include "file_reader.h"
#include "linkage_operations.h"
#include "open_hash/little_hash_for_linkage.h"


void linkage_build_link_list_from_file(char* input_file_name, dBGraph* db_graph, LinkHashTable* link_list);
void get_nearby_branch_for_a_read(dBNode* node, LinkageElement** candidate_list, dBGraph* db_graph, LinkHashTable* link_list, short* count, int colour);



#endif /* BUILD_LINKAGE_FROM_FILES_H_ */
