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

void load_se_seq_data_into_link_list(
				     const char *file_path,
				     char quality_cutoff, int homopolymer_cutoff, boolean remove_dups_se,
				     char ascii_fq_offset, int colour_index, dBGraph *db_graph, LinkHashTable *link_list,
				     unsigned long long *bad_reads, unsigned long long *dup_reads,
				     unsigned long long *bases_read, unsigned long long *bases_loaded,
				     unsigned long *readlen_count_array, unsigned long readlen_count_array_size);
void load_pe_seq_data_into_link_list(
				     const char *file_path1, const char *file_path2,
				     char quality_cutoff, int homopolymer_cutoff, boolean remove_dups_pe,
				     char ascii_fq_offset, int colour_index, dBGraph *db_graph, LinkHashTable *link_list,
				     unsigned long long *bad_reads, unsigned long long *dup_reads,
				     unsigned long long *bases_read, unsigned long long *bases_loaded,
				     unsigned long *readlen_count_array, unsigned long readlen_count_array_size);
void load_se_filelist_into_link_list(
				     char* se_filelist_path,
				     int qual_thresh, int homopol_limit, boolean remove_dups_se,
				     char ascii_fq_offset, int colour, dBGraph* db_graph, LinkHashTable *link_list,
				     unsigned int *total_files_loaded,
				     unsigned long long *total_bad_reads, unsigned long long *total_dup_reads,
				     unsigned long long *total_bases_read, unsigned long long *total_bases_loaded,
				     unsigned long *readlen_count_array, unsigned long readlen_count_array_size);
void load_pe_filelists_into_link_list(
				      char* pe_filelist_path1, char* pe_filelist_path2,
				      int qual_thresh, int homopol_limit, boolean remove_dups_pe,
				      char ascii_fq_offset, int colour, dBGraph* db_graph, LinkHashTable *link_list,
				      unsigned int *total_file_pairs_loaded,
				      unsigned long long *total_bad_reads, unsigned long long *total_dup_reads,
				      unsigned long long *total_bases_read, unsigned long long *total_bases_loaded,
				      unsigned long *readlen_count_array, unsigned long readlen_count_array_size);





#endif /* BUILD_LINKAGE_FROM_FILES_H_ */
