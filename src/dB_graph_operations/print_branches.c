/*
  print_branches.c - print branch sequences
*/

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <limits.h>
#include <string.h>

// cortex_var headers
#include "print_branches.h"
#include "open_hash/little_hash_for_linkage.h"

int db_graph_get_sequence_for_a_single_branch_for_specific_person_or_pop(dBNode* node, Orientation orientation, dBGraph *db_graph, 
									 char* seq, boolean* flag,int index)
{
  if (node == NULL)
    {
      die("Invalid starting node!\n");
    }

  char tmp_seq[db_graph->kmer_size+1];
  char initia[]="";
  strcpy(tmp_seq, initia);

  int length=1;

  BinaryKmer local_copy_of_kmer, tmp_kmer;
  binary_kmer_assignment_operator(local_copy_of_kmer, node->kmer);
  BinaryKmer* rev_kmer = binary_kmer_reverse_complement(&local_copy_of_kmer,db_graph->kmer_size, &tmp_kmer);

  if (orientation == forward) {
    binary_kmer_to_seq(&local_copy_of_kmer,db_graph->kmer_size,tmp_seq);
  }
  else {
    binary_kmer_to_seq(rev_kmer,db_graph->kmer_size,tmp_seq);
  }
  strcpy(seq, tmp_seq);
  db_node_set_status(node, visited);


  Nucleotide reverse_edge=0;
  Orientation next_orientation=0;
  Orientation curr_orientation=0;

  dBNode *node_ptr=NULL, *next_node=NULL;

  if (db_node_get_degree(node, index, orientation) > 1)
    {
      *flag = false;
      return length;
    }

  else if (db_node_get_degree(node, index, orientation) == 0)
    {
      return length;
    }

  else {
    short i;
    for (i=0; i<4; i++)
      {
        if (db_node_edge_exist(node, i, orientation, index)) {
          next_orientation=0;
          reverse_edge=0;
          node_ptr=db_graph_get_next_node_for_specific_person_or_pop(node,orientation,&next_orientation,i,&reverse_edge,db_graph, index);
          char tmp_char = binary_nucleotide_to_char(i);
          strncat(seq, &tmp_char, 1);
	}
      }
  }

  //  if (node_ptr==NULL) {printf("%s\n","Oh!");}
  //  printf("%d\n", length);

  curr_orientation=next_orientation;
  int degree1 = db_node_get_degree(node_ptr, 0, 0);
  int degree2 = db_node_get_degree(node_ptr, 0, 1);

  while ((degree1 == 1) && (degree2 == 1) && (length <= 300))
    {
      //If in a circle
      if (db_node_check_status(node_ptr, visited)) {*flag = false; break;}

      length++;
      db_node_set_status(node_ptr, visited);

      next_node = NULL;
      short j;
      for (j=0; j<4; j++)
        {
          if (db_node_edge_exist(node_ptr, j, curr_orientation, index))
            {
              next_orientation=0;
              reverse_edge=0;
              next_node = db_graph_get_next_node_for_specific_person_or_pop(node_ptr,curr_orientation,&next_orientation,j,&reverse_edge,db_graph, index);
              char tmp_char = binary_nucleotide_to_char(j);
              strncat(seq, &tmp_char, 1);
            }
        }

      if (next_node == NULL) {die("Unexpected disruption when calculating perfect path coverage in colour %d\n", index);}

      node_ptr = next_node;
      curr_orientation=next_orientation;
      degree1 = db_node_get_degree(node_ptr, index, forward);
      degree2 = db_node_get_degree(node_ptr, index, reverse);
    }//while

  if ((curr_orientation == forward) && (degree1 != 1) && (degree2 == 1))
    {
      db_node_set_status(node_ptr, visited);
      length++;
      *flag = false;
    }
  else if ((curr_orientation == reverse) && (degree2 != 1) && (degree1 == 1))
    {
      db_node_set_status(node_ptr, visited);
      length++;
      *flag = false;
    }
  else {seq[strlen(seq)-1]='\0';}

  //  printf("%d\n%s\n", strlen(seq), seq);
  return length;
}


void db_graph_print_all_branches_for_specific_person_or_pop(char* out_file_name, dBGraph* db_graph, int index)
{
  hash_table_traverse(&db_node_action_unset_status_visited_or_visited_and_exists_in_reference, db_graph);

  dBNode *curr_node=NULL, *next_node=NULL;
  int length=0, count=0;
  Orientation next_orientation=0;
  Nucleotide reverse_edge=0;
  char seq[500];
  char initia[]="";
  FILE* output = fopen(out_file_name, "w");

  long long i;
  for(i=0;i<db_graph->number_buckets * db_graph->bucket_size;i++)
    {
      curr_node = &db_graph->table[i];
      if ( (curr_node->coverage[index] > 0) && (!db_node_check_status(curr_node,pruned)) &&
           ((db_node_get_degree(curr_node, index, 0) > 1) || ((db_node_get_degree(curr_node, index, 1) > 1))) )

        {
          short k;
          for(k=0; k < 2; k++)
            {
              int degree = db_node_get_degree(curr_node, index, k);

              if (degree > 1)
                {
                  short j;
                  for(j=0; j < 4; j++)
                    {
                      if (db_node_edge_exist(curr_node,j,k,index))
                        {
                          reverse_edge=0;
                          next_orientation=0;
                          next_node = db_graph_get_next_node_for_specific_person_or_pop(curr_node,k,&next_orientation,j,&reverse_edge,db_graph,index);

                          if (!db_node_check_status(next_node,visited))
                            {
                              boolean flag=true;
                              strcpy(seq, initia);

                              length = db_graph_get_sequence_for_a_single_branch_for_specific_person_or_pop(next_node,next_orientation,db_graph,
													    seq,&flag,index);

                              if (length > 0)
                                {
                                  fprintf(output, "%s", ">b_");
				  short p;
				  for(p=0; p < db_graph->kmer_size; p++)
				    {
				      fprintf(output, "%c", seq[p]);
				    }
				  fprintf(output, "\n");
                                  fprintf(output, "%s\n", seq);

                                  count++;
                                }
                            }
                        }//if edge exist
                    }//for
                }//else if

            }
        }
    }

  fclose(output);
  return;
}


boolean db_graph_insert_branching_node_for_specific_person_or_pop(dBNode* node,dBGraph *db_graph,LinkHashTable *link_list,int index)
{
  if (node == NULL)
    {
      die("Invalid starting node!\n");
    }

  boolean found;
  BinaryKmer* curr_kmer = element_get_kmer(node);
  BinaryKmer local_copy;
  LinkageElement* new_element = link_hash_table_find_or_insert(element_get_key(curr_kmer,link_list->kmer_size,&local_copy),&found,link_list);
  db_node_set_status(node, visited);

  if (found || (new_element == NULL)) {return false;}
  else {return true;}
}


void db_graph_create_little_hash_table_for_branching_nodes_for_specific_person_or_pop(dBGraph* db_graph, LinkHashTable *link_list, int index)
{
  hash_table_traverse(&db_node_action_unset_status_visited_or_visited_and_exists_in_reference, db_graph);

  dBNode *curr_node=NULL, *next_node=NULL;
  Orientation next_orientation;
  Nucleotide reverse_edge;

  long long i;
  for(i=0;i<db_graph->number_buckets * db_graph->bucket_size;i++)
    {
      curr_node = &db_graph->table[i];
      if ( (curr_node->coverage[index] > 0) && (!db_node_check_status(curr_node,pruned)) &&
           ((db_node_get_degree(curr_node, index, 0) > 1) || ((db_node_get_degree(curr_node, index, 1) > 1))) )

        {
          short k;
          for(k=0; k < 2; k++)
            {
              int degree = db_node_get_degree(curr_node, index, k);
              if (degree > 1)
                {
                  short j;
                  for(j=0; j < 4; j++)
                    {
                      if (db_node_edge_exist(curr_node,j,k,index))
			{
                          reverse_edge=0;
                          next_orientation=0;
                          next_node = db_graph_get_next_node_for_specific_person_or_pop(curr_node,k,&next_orientation,j,&reverse_edge,db_graph,index);

                          if (!db_node_check_status(next_node,visited))
                            {
                              boolean if_insert;
                              if_insert = db_graph_insert_branching_node_for_specific_person_or_pop(next_node,db_graph,link_list,index);
			      /*
                              if (if_insert)
                                {
                                  char seq[db_graph->kmer_size];
                                  BinaryKmer local_copy_of_kmer;
                                  binary_kmer_to_seq(&local_copy_of_kmer,db_graph->kmer_size,seq);
                                  printf("%s\n", seq);
                                }
			      */
                            }
                        }//if edge exist
                    }//for
		}//else if

            }
        }
    }

  return;
}
