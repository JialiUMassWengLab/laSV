/*
  build_linkage_from_files.c
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <inttypes.h>

// cortex_var headers
#include "build_linkage_from_files.h"
#define is_base_char(x) ((x) == 'a' || (x) == 'A' || \
                         (x) == 'c' || (x) == 'C' || \
                         (x) == 'g' || (x) == 'G' || \
                         (x) == 't' || (x) == 'T')


void invalid_base_warning1(SeqFile *sf, char b)
{
  warn("Invalid sequence [%c] [path: %s; line: %lu]\n",
       b, seq_get_path(sf), seq_curr_line_number(sf));
}

short _kmer_errors1(SeqFile *sf, char *kmer_str, char* qual_str,
		    short kmer_size, char read_qual,
		    char quality_cutoff, int homopolymer_cutoff,
		    char *skip_homorun_base) // set to != 0 if need to skip a base
{
  short num_bp_to_skip = 0;
  int i;

  for(i = kmer_size-1; i >= 0; i--)
    {
      if(!is_base_char(kmer_str[i]))
	{
	  if(kmer_str[i] != 'n' && kmer_str[i] != 'N')
	    {
	      // Invalid base
	      invalid_base_warning1(sf, kmer_str[i]);
	    }

	  // skip
	  num_bp_to_skip = i+1;
	  break;
	}
    }

  if(read_qual)
    {
      // Find latest poor-quality base position
      for(i = kmer_size-1; i >= 0; i--)
	{
	  // Dev: we don't do any validation on the quality scores
	  if(qual_str[i] <= quality_cutoff)
	    {
	      num_bp_to_skip = i+1;
	      break;
	    }
	}
    }

  if(homopolymer_cutoff != 0)
    {
      int run_length = 1;

      for(i = 1; i < kmer_size; i++)
	{
	  if(kmer_str[i] == kmer_str[i-1])
	    run_length++;
	  else
	    run_length = 1;

	  if(run_length >= homopolymer_cutoff)
	    num_bp_to_skip = MAX(i+1, num_bp_to_skip);
	}

      *skip_homorun_base = (run_length >= homopolymer_cutoff ? kmer_str[kmer_size-1]
			    : 0);
    }

  return num_bp_to_skip;
}

inline char _read_base1(SeqFile *sf, char *b, char *q, char read_qual)
{
  if(!seq_read_base(sf, b))
    {
      return 0;
    }

  if(!is_base_char(*b) && *b != 'n' && *b != 'N')
    {
      // Invalid base
      invalid_base_warning1(sf, *b);

      // Default to 'N'
      *b = 'N';
    }

  if(read_qual && !seq_read_qual(sf, q))
    {
      warn("%s:%d: Couldn't read quality scores [read: %s; path: %s; line: %lu]",
	   __FILE__, __LINE__, seq_get_read_name(sf), seq_get_path(sf),
	   seq_curr_line_number(sf));
    }

  // Convert to upper case
  *b = toupper(*b);

  return 1;
}

inline char _read_k_bases1(SeqFile *sf, char *bases, char *quals, int k,
			   char read_qual)
{
  if(!seq_read_k_bases(sf, bases, k))
    {
      return 0;
    }

  if(read_qual && !seq_read_k_quals(sf, quals, k))
    {
      warn("%s:%d: Couldn't read quality scores [read: %s; path: %s; line: %lu]",
	   __FILE__, __LINE__, seq_get_read_name(sf), seq_get_path(sf),
	   seq_curr_line_number(sf));
    }

  // Convert to upper case
  int i;
  for(i = 0; i < k; i++)
    bases[i] = toupper(bases[i]);

  return 1;
}

char _read_first_kmer1(SeqFile *sf, char *kmer_str, char* qual_str,
		       short kmer_size, char read_qual,
		       char quality_cutoff, int homopolymer_cutoff,
		       char first_base, char first_qual)
{
  short bp_to_skip;

  if(first_base != 0)
    {
      // Already got one base
      kmer_str[0] = first_base;
      qual_str[0] = first_qual;

      if(!_read_k_bases1(sf, kmer_str+1, qual_str+1, kmer_size-1, read_qual))
	return 0;
    }
  else if(!_read_k_bases1(sf, kmer_str, qual_str, kmer_size, read_qual))
    {
      return 0;
    }

  char skip_homorun_base = 0;

  while((bp_to_skip = _kmer_errors1(sf, kmer_str, qual_str, kmer_size, read_qual,
				    quality_cutoff, homopolymer_cutoff,
				    &skip_homorun_base)) > 0)
    {
      //printf("proposed kmer: '%s' bases to skip: %i\n", kmer_str, bp_to_skip);

      while(bp_to_skip == kmer_size)
	{
	  // Re read whole kmer
	  if(!_read_k_bases1(sf, kmer_str, qual_str, kmer_size, read_qual))
	    {
	      return 0;
	    }

	  bp_to_skip = 0;
	  while(bp_to_skip < kmer_size && kmer_str[bp_to_skip] == skip_homorun_base)
	    {
	      bp_to_skip++;
	    }
	}

      if(bp_to_skip > 0)
	{
	  // Skip bp_to_skip bases
	  int bp_to_keep = kmer_size - bp_to_skip;

	  memmove(kmer_str, kmer_str+bp_to_skip, bp_to_keep);
	  memmove(qual_str, qual_str+bp_to_skip, bp_to_keep);

	  if(!_read_k_bases1(sf, kmer_str+bp_to_keep, qual_str+bp_to_keep,
			     bp_to_skip, read_qual))
	    {
	      return 0;
	    }
	}
    }

  return 1;
}

void get_nearby_branch_for_a_read(dBNode* node, LinkageElement** candidate_list, dBGraph* db_graph, LinkHashTable* link_list, short* count, int colour)
{
  dBNode* next = NULL;
  dBNode* tmp_node = NULL;
  Orientation next_orientation, tmp_orientation;
  Nucleotide reverse_edge;

  short i;
  for (i=0; i < 2; i++)
    {
      if (db_node_get_degree(node, colour, i) == 1)
	{
	  short j;
	  for (j=0; j < 4; j++)
	    {
	      if (db_node_edge_exist(node, j, i, colour))
		{
		  tmp_node = db_graph_get_next_node_for_specific_person_or_pop(node,i,&tmp_orientation,j,&reverse_edge,db_graph,colour);
		  break;
		}
	    }
	     
	  BinaryKmer* local_copy_of_kmer=element_get_kmer(tmp_node);
	  BinaryKmer tmp_kmer;
	  LinkageElement* mirror=link_hash_table_find(element_get_key(local_copy_of_kmer,db_graph->kmer_size, &tmp_kmer),link_list);
	  short counter=0;

	  while (db_node_get_degree(tmp_node,colour,tmp_orientation) == 1 && counter < 200 && mirror == NULL)
	    {
	      short k;
	      for (k=0; k < 4; k++)
		{
		  if (db_node_edge_exist(tmp_node, k, tmp_orientation, colour))
		    {
		      next = db_graph_get_next_node_for_specific_person_or_pop(tmp_node,tmp_orientation,&next_orientation,k,&reverse_edge,db_graph,colour);
		      break;
		    }
		}

	      tmp_node = next;
	      tmp_orientation = next_orientation;		  
	      local_copy_of_kmer=element_get_kmer(tmp_node);
	      mirror=link_hash_table_find(element_get_key(local_copy_of_kmer,db_graph->kmer_size, &tmp_kmer),link_list);
	      
	      counter++;
	    }

	  if (mirror != NULL)
	    {
	      candidate_list[*count]=mirror;
	      (*count)++;
	    }
	}
    }

  return;
}


inline void _process_read1(SeqFile *sf, char* kmer_str, char* qual_str,
			   char quality_cutoff, int homopolymer_cutoff,
			   dBGraph *db_graph, LinkHashTable *link_list, 
			   int colour_index, short* count,
			   BinaryKmer curr_kmer, Element* curr_node, LinkageElement *curr_link_node,
			   LinkageElement** candidate_list,
			   unsigned long long *bases_loaded,
			   unsigned long *readlen_count_array,
			   unsigned long readlen_count_array_size)
{
  // Hash table stuff
  BinaryKmer tmp_key;

  short kmer_size = link_list->kmer_size;
  char read_qual = seq_has_quality_scores(sf);

  char base, qual, prev_base;
  unsigned long contig_length;

  int homopol_length = 1;

  char keep_reading = 1;
  char is_first_kmer = 1;

  #ifdef DEBUG_CONTIGS
  printf(">%s\n", seq_get_read_name(sf));
  #endif

  while(keep_reading)
    {
      if(!is_first_kmer)
	{
	  seq_to_binary_kmer(kmer_str, kmer_size, (BinaryKmer*)curr_kmer);

	  linkage_element_get_key((BinaryKmer*)curr_kmer, kmer_size, &tmp_key);
	  curr_link_node = link_hash_table_find(&tmp_key, link_list);

	  if (curr_link_node != NULL)
	    {
	      candidate_list[*count]=curr_link_node;
	      (*count)++;
	    }
	}

    #ifdef DEBUG_CONTIGS
      _print_kmer(curr_kmer, kmer_size);
    #endif

      is_first_kmer = 0;

      contig_length = kmer_size;
      prev_base = kmer_str[kmer_size-1];

      // Set homopol_length for the first kmer in this contig
      if(homopolymer_cutoff != 0)
	{
	  homopol_length = 1;

	  int i;
	  for(i = kmer_size-2; i >= 0 && kmer_str[i] == prev_base; i--)
	    homopol_length++;
	}

      while((keep_reading = _read_base1(sf, &base, &qual, read_qual)))
	{
	  // Check for Ns and low quality scores
	  if(base == 'N' || (read_qual && qual <= quality_cutoff))
	    {
	      keep_reading = _read_first_kmer1(sf, kmer_str, qual_str,
					       kmer_size, read_qual,
					       quality_cutoff, homopolymer_cutoff,
					       0, 0);
	      break;
	    }

	  // Check homopolymer run length
	  if(homopolymer_cutoff != 0)
	    {
	      if(prev_base == base)
		{
		  homopol_length++;

		  if(homopol_length >= homopolymer_cutoff)
		    {
		      // Skip the rest of these bases
		      while((keep_reading = _read_base1(sf, &base, &qual, read_qual)) &&
			    base == prev_base);

		      // Pass on first base that is != prev_base
		      keep_reading = _read_first_kmer1(sf, kmer_str, qual_str,
						       kmer_size, read_qual,
						       quality_cutoff, homopolymer_cutoff,
						       base, qual);
		      break;
		    }
		}
	      else
		{
		  // Reset homopolymer length
		  homopol_length = 1;
		}
	    }

	  // base accepted
	  contig_length++;

	  // Construct new kmer
	  //binary_kmer_assignment_operator(curr_kmer, prev_kmer);
	  Nucleotide nuc = char_to_binary_nucleotide(base);
	  binary_kmer_left_shift_one_base_and_insert_new_base_at_right_end((BinaryKmer*)curr_kmer,
									   nuc,
									   kmer_size);

      #ifdef DEBUG_CONTIGS
	  printf("%c", base);
      #endif

	  // Lookup in db
	  linkage_element_get_key((BinaryKmer*)curr_kmer, kmer_size, &tmp_key);
	  curr_link_node = link_hash_table_find(&tmp_key, link_list);

	  if (curr_link_node != NULL)
	    {
              candidate_list[*count]=curr_link_node;
              (*count)++;
	    }	      
	}

      // Store bases that made it into the graph
      (*bases_loaded) += contig_length;

    #ifdef DEBUG_CONTIGS
      printf(" [%lu]\n", contig_length);
    #endif

      // Store contig length
      if(readlen_count_array != NULL)
	{
	  contig_length = MIN(contig_length, readlen_count_array_size-1);
	  readlen_count_array[contig_length]++;
	}
    }

  if (*count == 0)
    {
      get_nearby_branch_for_a_read(curr_node,candidate_list,db_graph,link_list,count,colour_index);
    }
}

void load_se_seq_data_into_link_list(
				     const char *file_path,
				     char quality_cutoff, int homopolymer_cutoff, boolean remove_dups_se,
				     char ascii_fq_offset, int colour_index, dBGraph *db_graph, LinkHashTable *link_list,
				     unsigned long long *bad_reads, // number of reads that have no good kmers
				     unsigned long long *dup_reads, // number of reads that pcr duplicates
				     unsigned long long *bases_read, // total bases in file
				     unsigned long long *bases_loaded, // bases that make it into the graph
				     unsigned long *readlen_count_array, // histogram of contigs lengths
				     unsigned long readlen_count_array_size) // contigs bigger go in final bin
{
  short kmer_size = db_graph->kmer_size;

  // DEV:
  // First check if this is a cortex binary
  // -> Load binary

  // Open file
  SeqFile *sf = seq_file_open(file_path);

  if(sf == NULL)
    {
      die("Couldn't open single-end sequence file '%s'\n", file_path);
    }

  seq_set_fastq_ascii_offset(sf, ascii_fq_offset);

  // Are we using quality scores
  char read_qual = seq_has_quality_scores(sf);

  // Vars for reading in
  char kmer_str[kmer_size+1];
  char qual_str[kmer_size+1];

  //
  // Hash table variables
  BinaryKmer curr_kmer;

  // Element and dBNode are the same thing
  LinkageElement *curr_link_node = NULL;
  Element *curr_node = NULL;
  Orientation curr_orient = forward;
  BinaryKmer tmp_key;
  boolean curr_found;

  while(seq_next_read(sf))
    {
      //printf("Started seq read: %s\n", seq_get_read_name(sf));

      if(_read_first_kmer1(sf, kmer_str, qual_str, kmer_size, read_qual,
			   quality_cutoff, homopolymer_cutoff, 0, 0))
	{
	  // Check if we want this read
	  char is_dupe = 0;

	  // Look up in table
	  curr_found = false;

	  seq_to_binary_kmer(kmer_str, kmer_size, &curr_kmer);
	  element_get_key(&curr_kmer, kmer_size, &tmp_key);
          curr_node = hash_table_find_or_insert(&tmp_key, &curr_found, db_graph);
          curr_orient = db_node_get_orientation(&curr_kmer, curr_node, kmer_size);

	  linkage_element_get_key(&curr_kmer, kmer_size, &tmp_key);
	  curr_link_node = link_hash_table_find(&tmp_key, link_list);

	  if(remove_dups_se == true)
	    {
	      if(curr_found == true &&
		 db_node_check_read_start(curr_node, curr_orient) == true)
		{
		  (*dup_reads)++;
		  is_dupe = 1;
		}
	      else
		{
		  db_node_set_read_start_status(curr_node, curr_orient);
		}
	    }

	  //debug
	  //char tmpstr[32];
	  //binary_kmer_to_seq(curr_kmer, kmer_size, tmpstr);
	  //printf("First kmer: %s\n", tmpstr);

	  if(!is_dupe)
	    {
	      LinkageElement** candidate_list = (LinkageElement**) calloc(75, sizeof(LinkageElement*));
	      short count=0;

	      if (curr_link_node != NULL) 
		{
		  candidate_list[0] = curr_link_node;
		}

	      _process_read1(sf, kmer_str, qual_str,
			     quality_cutoff, homopolymer_cutoff,
			     db_graph, link_list, colour_index, &count,
			     curr_kmer, curr_node, curr_link_node, 
			     candidate_list, bases_loaded,
			     readlen_count_array, readlen_count_array_size);

	      short j;
	      for (j=0; j < count; j++)
		{
		  short k;
		  for (k=j+1; k < count; k++)
		    {
		      linkage_element_add_a_link(candidate_list[j], candidate_list[k], db_graph);
		      linkage_element_add_a_link(candidate_list[k], candidate_list[j], db_graph);
		    }
		}

	      free(candidate_list);	      
	    }
	}
      else
	{
	  // Couldn't get a single kmer from read
	  (*bad_reads)++;
	}
    }

  // Update with bases read in
  (*bases_read) += seq_total_bases_passed(sf) + seq_total_bases_skipped(sf);

  seq_file_close(sf);
}

void load_pe_seq_data_into_link_list(
				     const char *file_path1, const char *file_path2,
				     char quality_cutoff, int homopolymer_cutoff, boolean remove_dups_pe,
				     char ascii_fq_offset, int colour_index, dBGraph *db_graph, LinkHashTable* link_list,
				     unsigned long long *bad_reads, // number of reads that have no good kmers
				     unsigned long long *dup_reads, // number of reads that are pcr duplicates
				     unsigned long long *bases_read, // total bases in file
				     unsigned long long *bases_loaded, // bases that make it into the graph
				     unsigned long *readlen_count_array, // length of contigs loaded
				     unsigned long readlen_count_array_size) // contigs bigger go in final bin
{
  short kmer_size = db_graph->kmer_size;

  // Open files
  SeqFile *sf1 = seq_file_open(file_path1);
  SeqFile *sf2 = seq_file_open(file_path2);

  if(sf1 == NULL)
    {
      die("Couldn't open paired-end sequence file '%s'\n", file_path1);
    }
  else if(sf2 == NULL)
    {
      die("Couldn't open paired-end sequence file '%s'\n", file_path2);
    }

  seq_set_fastq_ascii_offset(sf1, ascii_fq_offset);
  seq_set_fastq_ascii_offset(sf2, ascii_fq_offset);

  // Are we using quality scores
  char read_qual1 = seq_has_quality_scores(sf1);
  char read_qual2 = seq_has_quality_scores(sf2);

  // Vars for reading in
  char kmer_str1[kmer_size+1], kmer_str2[kmer_size+1];
  char qual_str1[kmer_size+1], qual_str2[kmer_size+1];

  //
  // Hash table variables
  BinaryKmer curr_kmer1, curr_kmer2;

  // Element and dBNode are the same thing
  LinkageElement *curr_link_node1 = NULL, *curr_link_node2 = NULL;
  Element *curr_node1 = NULL, *curr_node2 = NULL;
  Orientation curr_orient1 = forward, curr_orient2 = forward;
  BinaryKmer tmp_key;
  boolean curr_found1, curr_found2;

  while(1)
    {
      char read1 = seq_next_read(sf1);
      char read2 = seq_next_read(sf2);

      if(read1 != read2)
	{
	  die("Paired-end files don't have the same number of reads\n"
          "file1: %s\n"
	      "file2: %s\n",
	      file_path1, file_path2);
	}
      else if(!read1)
	{
	  // Neither have read
	  break;
	}

      read1 = _read_first_kmer1(sf1, kmer_str1, qual_str1, kmer_size, read_qual1,
				quality_cutoff, homopolymer_cutoff, 0, 0);

      read2 = _read_first_kmer1(sf2, kmer_str2, qual_str2, kmer_size, read_qual2,
				quality_cutoff, homopolymer_cutoff, 0, 0);

      if(read1)
	{
	  curr_found1 = false;

	  // Get binary kmer
	  seq_to_binary_kmer(kmer_str1, kmer_size, &curr_kmer1);

	  // Look up first kmer
	  element_get_key(&curr_kmer1, kmer_size, &tmp_key);
	  curr_node1 = hash_table_find_or_insert(&tmp_key, &curr_found1, db_graph);
	  curr_orient1 = db_node_get_orientation(&curr_kmer1, curr_node1, kmer_size);

	  linkage_element_get_key(&curr_kmer1, kmer_size, &tmp_key);
	  curr_link_node1 = link_hash_table_find(&tmp_key,link_list);
	}
      else
	{
	  // Couldn't get a single kmer from read
	  (*bad_reads)++;
	}

      if(read2)
	{
	  curr_found2 = false;

	  // Get binary kmer
	  seq_to_binary_kmer(kmer_str2, kmer_size, &curr_kmer2);

	  // Look up second kmer
	  element_get_key(&curr_kmer2, kmer_size, &tmp_key);
	  curr_node2 = hash_table_find_or_insert(&tmp_key, &curr_found2, db_graph);
	  curr_orient2 = db_node_get_orientation(&curr_kmer2, curr_node2, kmer_size);

	  linkage_element_get_key(&curr_kmer2, kmer_size, &tmp_key);
	  curr_link_node2 = link_hash_table_find(&tmp_key,link_list);
	}
      else
	{
	  // Couldn't get a single kmer from read
	  (*bad_reads)++;
	}

      char is_dupe = 0;

      if(remove_dups_pe == true && (read1 || read2))
	{
	  // Check if all reads that were read in are marked as dupe => ie.
	  //   A) Both reads read in and both are marked as read starts OR
	  //   B) Just one read in and it is marked as read start

	  if((!read1 || (curr_found1 == true &&
			 db_node_check_read_start(curr_node1, curr_orient1) == true)) &&
	     (!read2 || (curr_found2 == true &&
			 db_node_check_read_start(curr_node2, curr_orient2) == true)))
	    {
	      // Both reads are dupes or neither are
	      is_dupe = 1;
	      (*dup_reads) += 2;
	    }
	  else
	    {
	      if(read1)
		db_node_set_read_start_status(curr_node1, curr_orient1);

	      if(read2)
		db_node_set_read_start_status(curr_node2, curr_orient2);
	    }
	}

      /*
    // For debugging
    char tmpstr1[kmer_size+1], tmpstr2[kmer_size+1];
    binary_kmer_to_seq((BinaryKmer*)curr_kmer1, kmer_size, tmpstr1);
    binary_kmer_to_seq((BinaryKmer*)curr_kmer2, kmer_size, tmpstr2);
    printf("Paired first kmers: %s %s\n", read1 ? tmpstr1 : "", read2 ? tmpstr2 : "");
      */

      if(!is_dupe)
	{
	  LinkageElement** candidate_list = (LinkageElement**) calloc(150, sizeof(LinkageElement*));
	  short count=0;

	  if(read1)
	    {
	      if (curr_link_node1 != NULL)
		{
		  candidate_list[count] = curr_link_node1;
		}

	      _process_read1(sf1, kmer_str1, qual_str1,
			     quality_cutoff, homopolymer_cutoff,
			     db_graph, link_list, colour_index, &count,
			     curr_kmer1, curr_node1, curr_link_node1, candidate_list,
			     bases_loaded,
			     readlen_count_array, readlen_count_array_size);
	    }

	  if(read2)
	    {
              if (curr_link_node1 != NULL)
                {
                  candidate_list[count] = curr_link_node1;
                }

	      _process_read1(sf2, kmer_str2, qual_str2,
			     quality_cutoff, homopolymer_cutoff,
			     db_graph, link_list, colour_index, &count,
			     curr_kmer2, curr_node2, curr_link_node2, candidate_list,
			     bases_loaded,
			     readlen_count_array, readlen_count_array_size);
	    }

	  short j;
	  for (j=0; j < count; j++)
	    {
	      short k;
	      for (k=j+1; k < count; k++)
		{
		  linkage_element_add_a_link(candidate_list[j], candidate_list[k], db_graph);
		  linkage_element_add_a_link(candidate_list[k], candidate_list[j], db_graph);
		}
	    }

	  free(candidate_list);
	}
    }

  // Update with bases read in
  (*bases_read) += seq_total_bases_passed(sf1) + seq_total_bases_skipped(sf1) +
    seq_total_bases_passed(sf2) + seq_total_bases_skipped(sf2);

  seq_file_close(sf1);
  seq_file_close(sf2);
}

void load_se_filelist_into_link_list(
				     char* se_filelist_path,
				     int qual_thresh, int homopol_limit, boolean remove_dups_se,
				     char ascii_fq_offset, int colour, dBGraph* db_graph, LinkHashTable* link_list,
				     unsigned int *total_files_loaded,
				     unsigned long long *total_bad_reads, unsigned long long *total_dup_reads,
				     unsigned long long *total_bases_read, unsigned long long *total_bases_loaded,
				     unsigned long *readlen_count_array, unsigned long readlen_count_array_size)
{
qual_thresh += ascii_fq_offset;

// Get absolute path
char absolute_path[PATH_MAX+1];
char* se_filelist_abs_path = realpath(se_filelist_path, absolute_path);

if(se_filelist_abs_path == NULL)
  {
    die("Cannot get absolute path to filelist of SE files: %s\n", se_filelist_path);
  }

/* COMMENT_OUT_DURING_TESTS */
printf("  path: %s\n", se_filelist_abs_path);
/*  */
FILE* se_list_file = fopen(se_filelist_abs_path, "r");

if(se_list_file == NULL)
  {
    die("Cannot open filelist of SE files: %s\n", se_filelist_abs_path);
  }

// Get directory path
StrBuf *dir = file_reader_get_strbuf_of_dir_path(se_filelist_abs_path);

// Stats
unsigned int se_files_loaded = 0;

unsigned long long se_bad_reads = 0;
unsigned long long se_dup_reads = 0;

unsigned long long se_bases_read = 0;
unsigned long long se_bases_loaded = 0;

StrBuf *line = strbuf_new();

while(strbuf_reset_readline(line, se_list_file))
  {
    strbuf_chomp(line);

    if(strbuf_len(line) > 0)
      {
	// Get paths relative to filelist dir
	if(strbuf_get_char(line, 0) != '/')
	  strbuf_insert(line, 0, dir, 0, strbuf_len(dir));

	// Get absolute paths
	char* path_ptr = realpath(line->buff, absolute_path);

	if(path_ptr == NULL)
	  {
	    die("Cannot find sequence file: %s\n",line->buff);
	  }

	load_se_seq_data_into_link_list(path_ptr,
					qual_thresh, homopol_limit, remove_dups_se,
					ascii_fq_offset, colour, db_graph, link_list,
					&se_bad_reads, &se_dup_reads,
					&se_bases_read, &se_bases_loaded,
					readlen_count_array, readlen_count_array_size);

	se_files_loaded++;
      }
  }

 strbuf_free(line);
 strbuf_free(dir);

 // Finished reading single-ended file list
 fclose(se_list_file);

 // Update cumulative stats
 *total_files_loaded += se_files_loaded;
 *total_bad_reads += se_bad_reads;
 *total_dup_reads += se_dup_reads;
 *total_bases_read += se_bases_read;
 *total_bases_loaded += se_bases_loaded;


 // Print SE stats for this set of files
 /* COMMENT_OUT_DURING_TESTS */
 printf("\nNum SE files loaded:%u\n", se_files_loaded);
 printf("\tKmers:%llu\n", hash_table_get_unique_kmers(db_graph));
 printf("\tNumber of bad reads:%llu\n", se_bad_reads);
 printf("\tNumber of dupe reads:%llu\n", se_dup_reads);
 printf("\tSE sequence parsed:%llu\n", se_bases_read);
 printf("\tTotal SE sequence that passed filters:%llu\n", se_bases_loaded);
 /* */
}

void load_pe_filelists_into_link_list(
				      char* pe_filelist_path1, char* pe_filelist_path2,
				      int qual_thresh, int homopol_limit, boolean remove_dups_pe,
				      char ascii_fq_offset, int colour, dBGraph* db_graph, LinkHashTable* link_list,
				      unsigned int *total_file_pairs_loaded,
				      unsigned long long *total_bad_reads, unsigned long long *total_dup_reads,
				      unsigned long long *total_bases_read, unsigned long long *total_bases_loaded,
				      unsigned long *readlen_count_array, unsigned long readlen_count_array_size)
{
  qual_thresh += ascii_fq_offset;

  // Get absolute paths
  char absolute_path1[PATH_MAX+1], absolute_path2[PATH_MAX+1];

  char* pe_filelist_abs_path1 = realpath(pe_filelist_path1, absolute_path1);
  char* pe_filelist_abs_path2 = realpath(pe_filelist_path2, absolute_path2);

  if(pe_filelist_abs_path1 == NULL)
    {
      die("Cannot get absolute path to filelist of PE files: %s\n",
	  pe_filelist_path1);
    }
  else if(pe_filelist_abs_path2 == NULL)
    {
      die("Cannot get absolute path to filelist of PE files: %s\n",
	  pe_filelist_path2);
    }

  // Open files
  FILE* pe_list_file1 = fopen(pe_filelist_abs_path1, "r");
  FILE* pe_list_file2 = fopen(pe_filelist_abs_path2, "r");

  if(pe_list_file1 == NULL)
    {
      die("Cannot open filelist of PE files: %s\n", pe_filelist_abs_path1);
    }
  else if(pe_list_file2 == NULL)
    {
      die("Cannot open filelist of PE files: %s\n", pe_filelist_abs_path2);
    }

  // Get directory paths for filelist files
  StrBuf *dir1 = file_reader_get_strbuf_of_dir_path(pe_filelist_abs_path1);
  StrBuf *dir2 = file_reader_get_strbuf_of_dir_path(pe_filelist_abs_path2);

  // Stats
  unsigned int pe_file_pairs_loaded = 0;

  unsigned long long pe_bad_reads = 0;
  unsigned long long pe_dup_reads = 0;

  unsigned long long pe_bases_read = 0;
  unsigned long long pe_bases_loaded = 0;

  StrBuf *line1 = strbuf_new();
  StrBuf *line2 = strbuf_new();

  while(1)
    {
      char read1 = strbuf_reset_readline(line1, pe_list_file1) > 0;
      char read2 = strbuf_reset_readline(line2, pe_list_file2) > 0;

      if(!read1 && !read2)
	break;

      strbuf_chomp(line1);
      strbuf_chomp(line2);

      if((strbuf_len(line1) == 0) != (strbuf_len(line2) == 0))
	{
	  die("Paired-end files don't have the same number of lines:\n"
          "File 1: %s\n"
	      "File 2: %s", pe_filelist_path1, pe_filelist_path2);
	}
      else if(strbuf_len(line1) > 0)
	{
	  // Get paths relative to filelist dir
	  if(strbuf_get_char(line1, 0) != '/')
	    strbuf_insert(line1, 0, dir1, 0, strbuf_len(dir1));

	  if(strbuf_get_char(line2, 0) != '/')
	    strbuf_insert(line2, 0, dir2, 0, strbuf_len(dir2));

	  // Get absolute paths
	  char* path_ptr1 = realpath(line1->buff, absolute_path1);
	  char* path_ptr2 = realpath(line2->buff, absolute_path2);

	  if(path_ptr1 == NULL)
	    {
	      die("Cannot find sequence file: %s\n", line1->buff);
	    }
	  else if(path_ptr2 == NULL)
	    {
	      die("Cannot find sequence file: %s\n", line2->buff);
	    }

	  /* COMMENT_OUT_DURING_TESTS */
	  // Print file paths
	  printf("Paired-end seq files:\n");

	  printf("  File 1: %s\n", path_ptr1);
	  printf("  File 2: %s\n", path_ptr2);
	  /*   */

	  // Read
	  load_pe_seq_data_into_link_list(path_ptr1, path_ptr2,
					  qual_thresh, homopol_limit, remove_dups_pe,
					  ascii_fq_offset, colour, db_graph, link_list,
					  &pe_bad_reads, &pe_dup_reads,
					  &pe_bases_read, &pe_bases_loaded,
					  readlen_count_array, readlen_count_array_size);

	  pe_file_pairs_loaded++;
	}
    }

  /// Free line buffers
  strbuf_free(line1);
  strbuf_free(line2);

  // Free dir paths
  strbuf_free(dir1);
  strbuf_free(dir2);

  // Finished reading paired-ended file list
  fclose(pe_list_file1);
  fclose(pe_list_file2);

  // Update cumulative stats
  *total_file_pairs_loaded += pe_file_pairs_loaded;
  *total_bad_reads += pe_bad_reads;
  *total_dup_reads += pe_dup_reads;
  *total_bases_read += pe_bases_read;
  *total_bases_loaded += pe_bases_loaded;

  // Print SE stats for this set of files
  /* COMMENT_OUT_DURING_TESTS */
  printf("\nNum PE file pairs loaded:%u\n", pe_file_pairs_loaded);
  printf("\tKmers:%llu\n", hash_table_get_unique_kmers(db_graph));
  printf("\tNumber of bad reads:%llu\n", pe_bad_reads);
  printf("\tNumber of dupe reads:%llu\n", pe_dup_reads);
  printf("\tPE sequence parsed:%llu\n", pe_bases_read);
  printf("\tTotal PE sequence that passed filters:%llu\n", pe_bases_loaded);
}




void linkage_build_link_list_from_file(char* input_file_name, dBGraph* db_graph, LinkHashTable* link_list)
{
  short kmer_size = db_graph->kmer_size;
  FILE* input = fopen(input_file_name, "r");
  if (input == NULL)
    {
      die("linkage_build_link_list_from_file couldn't open file %s!\n\n", input_file_name);
    }

  char line[LINE_MAX]; //LINE_MAX is defined in limits.h
  char curr_seq[kmer_size+1];

  //  FILE* test = fopen("test_file", "w");
  while (fgets(line, LINE_MAX, input) != NULL)
    {
      int length=strlen(line);
      //  fprintf(test, "%d\n", length);

      LinkageElement** curr_list = (LinkageElement**) calloc(150, sizeof(LinkageElement*));
      if (curr_list == NULL)
	{
	  die("Could not allocate for curr_list!\n\n");
	}

      short i=0;
      int start=0;
      while (start+kmer_size < length)
	{
	  strncpy(curr_seq, line+start, kmer_size);
	  curr_seq[kmer_size]='\0';
	  //  fprintf(test, "%s\n", curr_seq);

	  boolean found;
	  BinaryKmer local_copy1, local_copy2;
	  BinaryKmer* curr_kmer = seq_to_binary_kmer(curr_seq,kmer_size,&local_copy1);
	  LinkageElement* new_element = link_hash_table_find_or_insert(element_get_key(curr_kmer,link_list->kmer_size,&local_copy2),&found,link_list);

	  if (new_element != NULL)
	    {
	      curr_list[i] = new_element;
	    }

	  i++;
	  start += kmer_size;
	}

      short j;
      for (j=0; j < i; j++)
	{
	  short k;
	  for (k=j+1; k < i; k++)
	    {
	      linkage_element_add_a_link(curr_list[j], curr_list[k], db_graph);
	      linkage_element_add_a_link(curr_list[k], curr_list[j], db_graph);
	    }
	}

      free(curr_list);
    }

  fclose(input);
  //  fclose(test);
  /*
  long long p;
  for (p=0; p < link_list->number_buckets * link_list->bucket_size; p++)
    {
      if (link_list->table[p].status > 0)
	{
	  short count=0;
	  LinkageElement node = link_list->table[p];
	  if (node.link == NULL) 
	    {
	      printf("%d\n", count);
	    }
	  else 
	    {
	      LinkNode* end=node.link;
	      while (end != NULL)
		{
		  end=end->next;
		  count++;
		}
	      count++;
	      printf("%d\n", count);
	    }
	}
    }
  */
}
