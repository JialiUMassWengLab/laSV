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
