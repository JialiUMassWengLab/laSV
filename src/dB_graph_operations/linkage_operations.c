/*
  linkage_operations.c
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <inttypes.h>
#include "linkage_operations.h"

void linkage_element_link_node_coverage_update (LinkNode* l)
{
  if (l->read == 0) {l->read=1;}
  else if (l->read == 1) {l->read=2;}
  else if (l->read == 2) {l->read=3;}
}


void linkage_element_add_a_link(LinkageElement* e1, LinkageElement* e2, dBGraph* db_graph)
{
  if ((e1 == NULL) || (e2 == NULL)) 
    {
      die("Trying to add link with null pointer...\n\n");
    }

  if (e1 == e2) {return;}
  
  if (e1->link_num >= 300) {return;}

  BinaryKmer* local_copy_of_kmer=linkage_element_get_kmer(e2);
  BinaryKmer tmp_kmer;
  dBNode* mirror=hash_table_find(element_get_key(local_copy_of_kmer,db_graph->kmer_size, &tmp_kmer),db_graph);

  LinkNode* end=e1->link;
  if (end == NULL)
    {
      LinkNode* new_node = calloc(1,sizeof(LinkNode));
      if (new_node == NULL)
	{
	  die("Couldn't allocate for new LinkNode!\n\n");
	}

      new_node->ptr=mirror;
      new_node->read=0;
      new_node->next=NULL;

      e1->link=new_node;
      e1->link_num++;
    }
  else 
    {
      while ((end->next != NULL)&&(end->ptr != mirror))
	{
	  end=end->next;
	}
      
      if ((end->next == NULL) && (end->ptr != mirror))
	{
	  LinkNode* new_node = calloc(1,sizeof(LinkNode));
	  if (new_node == NULL)
	    {
	      die("Couldn't allocate for new LinkNode!\n\n");
	    }

	  new_node->ptr=mirror;
	  new_node->read=0;
	  new_node->next=NULL;
	  
	  end->next=new_node;
	  e1->link_num++;
	}
      else 
	{
	  linkage_element_link_node_coverage_update(end);
	}
    }

  return;
}


boolean if_two_branches_are_connected(dBNode* candidate, dBNode** branch_list, unsigned short* distance, unsigned short* distal_length, 
				      int read_length, int insert_size, int min_sup, dBGraph* db_graph, LinkHashTable* link_list)
{
  boolean found=false;
  unsigned short current_length=*distal_length;

  if (branch_list[0] == NULL)
    {
      found=true;
      return found;
    }

  int bound = (int) (0.4 * read_length);
  int limit = (int) (1.4 * insert_size);

  unsigned short i=0;
  unsigned short proximal=0;
  while (branch_list[i] != NULL)
    {
      i++;
      if (current_length-distance[i] <= bound) {proximal++;}
    }

  BinaryKmer* local_copy_of_kmer=element_get_kmer(candidate);
  BinaryKmer tmp_kmer;
  LinkageElement* mirror=link_hash_table_find(element_get_key(local_copy_of_kmer,db_graph->kmer_size, &tmp_kmer),link_list);

  //  char tmp_seq[db_graph->kmer_size+1];
  //  printf("%s\n", binary_kmer_to_seq(&tmp_kmer,db_graph->kmer_size,tmp_seq));

  if (mirror == NULL)
    {
      //      die("Couldn't find the branch in link hash table!\nSomething must be wrong!\n\n");
      return false;
    }
    
  unsigned short farthest_link_dist=0;
  short num_proximal_linked_branch_nodes=0;
  LinkNode* curr_link_node=mirror->link;
  while (curr_link_node != NULL)
    {
      short j=i-1;
      while ((j >= 0) && (current_length-distance[j] <= limit))
	{
	  if ((curr_link_node->ptr == branch_list[j]) && (curr_link_node->read >= min_sup))
	    {
	      if (current_length-distance[j] <= bound) {num_proximal_linked_branch_nodes++;}

	      if (j == 0)
		{
		  (*distal_length) = db_graph->kmer_size;
		}
	      else if (distance[j-1] < (*distal_length))
		{
		  (*distal_length) = distance[j-1];
		  farthest_link_dist = distance[j];
		}

	      found=true;
	      break;
	    }
	  else 
	    j--;
	}
      
      curr_link_node = curr_link_node->next;
    }
  /*
  if (j > 8 && ((int) j/4) >= num_distal_linked_branch_nodes)
    {
      found=false;
    }
  */

  if ((proximal > num_proximal_linked_branch_nodes) && (current_length-farthest_link_dist < bound*2))
  //  if (proximal > num_proximal_linked_branch_nodes)
    {
      found=false;
    }

  return found;
}
  


void assemble_contigs_from_a_single_source(dBNode* source, Orientation initial_orientation, int read_length, int insert_size, int min_contig_length, int min_sup,
					   dBGraph* db_graph, LinkHashTable* link_list, unsigned int* count, FILE* output, int index)
{
  typedef struct
  {
    dBNode* ptr;
    Orientation orientation;
    Orientation orien[100];
    unsigned short extend;
    unsigned short revisit;
    dBNode* passed[100];
    dBNode* visited[100];
    unsigned short dist[100];
    char contig[50000];
  } queue_node;

  queue_node* dijkstra_queue;

  void print_one_contig(char* seq, unsigned short start, FILE* filename, unsigned int* num)
  {
    char* str_to_print = seq;

    str_to_print += start - db_graph->kmer_size;
    int length=strlen(str_to_print);
    //    printf("%d\n", length);
    //    printf("%s\n%s\n", seq, str_to_print);

    if (length >= min_contig_length)
      {
	if (length < 24999)
	  {
	    (*num)++;
	    fprintf(filename, "%s%d\t%d\n", ">Contig",*num,length);
	    fprintf(filename, "%s\n", str_to_print);
	  }
	else
	  {
	    char seq1[25000];
	    strncpy(seq1,str_to_print,24998);
	    seq1[24998]='\0';
	    (*num)++;
            fprintf(filename, "%s%d\t%d\n", ">Contig",*num,24998);
            fprintf(filename, "%s\n", seq1);

	    str_to_print += length-24998;
	    strncpy(seq1,str_to_print,24998);
	    seq1[24998]='\0';
            (*num)++;
            fprintf(filename, "%s%d\t%d\n", ">Contig",*num,24998);
            fprintf(filename, "%s\n", seq1);
	  }
      }

    return;
  }

  void assign_queue_node(queue_node* node1, queue_node node2)
  {
    node1->ptr = node2.ptr;
    node1->extend  = node2.extend;
    node1->revisit = node2.revisit;
    node1->orientation = node2.orientation;
    strcpy(node1->contig, node2.contig);

    short k;
    for (k=0; k < 100; k++) 
      {
	node1->visited[k] = node2.visited[k];
	node1->orien[k]   = node2.orien[k];
	node1->passed[k]  = node2.passed[k];
	node1->dist[k]    = node2.dist[k];
      }
  }

  void add_visited_node(queue_node* node, dBNode* previous, boolean* first_visit)
  {
    short visited_ptr=0;
    while ((visited_ptr < 100) && (node->visited[visited_ptr] != NULL)) {visited_ptr++;}

    if (visited_ptr < 100) 
      {
	node->visited[visited_ptr] = node->ptr;
	node->orien[visited_ptr]   = node->orientation;
	node->passed[visited_ptr]  = previous;
	node->dist[visited_ptr]    = strlen(node->contig);
      }
    else 
      {
	//printf("Too many visited nodes\nSet longer array\n\n");
	(*first_visit) = false;
      }

    return;
  }

  void run_along_a_supernode(queue_node* node, dBNode* previous_dB_node, boolean* first_visit, boolean if_set_visited)
  {
    while ((db_node_get_degree(node->ptr,index,node->orientation) == 1) && (*first_visit))
      {
	if (db_node_get_degree(node->ptr,index,opposite_orientation(node->orientation)) > 1)
	  {
	    add_visited_node(node, previous_dB_node, first_visit);
	  }

	dBNode* next=NULL;
	Orientation next_orientation;
	Nucleotide reverse_edge;

	char tmp_char[2];
	short l;
	for (l=0; l < 4; l++)
	  {
	    if (db_node_edge_exist(node->ptr, l, node->orientation, index))
	      {
		next=db_graph_get_next_node_for_specific_person_or_pop(node->ptr,node->orientation,&next_orientation,l,&reverse_edge,db_graph,index);
		tmp_char[0] = binary_nucleotide_to_char(l);
		tmp_char[1] = '\0';
		break;
	      }
	  }

	short k;
	for (k=0; k < 100; k++)
	  {
	    if ((next == node->visited[k])&&(next_orientation == node->orien[k])) {*first_visit=false; break;}
	  }

	if (*first_visit)
	  {
	    previous_dB_node = node->ptr;
	    node->ptr = next;
	    node->orientation = next_orientation;
	    if (db_node_check_status(node->ptr,visited)) 
	      node->revisit++;
	    else 
	      {
		node->revisit=0;
		if (if_set_visited) {db_node_set_status(node->ptr,visited);}
	      }
	    strncat(node->contig, tmp_char, 1);

	    if (strlen(node->contig) >= 49999) 
	      {
		(*first_visit) = false;
		//printf("string length limit reached\n\n");
	      }
	  }//if first_visit
      }//while

    if (db_node_get_degree(node->ptr,index,node->orientation) == 0) 
      {
	(*first_visit) = false;
      }

    return;
  }


  dijkstra_queue = (queue_node*) malloc(60000*sizeof(queue_node));

  if (dijkstra_queue == NULL) {die("Could not allocate a dijkstra queue with 60000 elements\n\n");}

  db_node_set_status(source, visited);

  dijkstra_queue[0].ptr=source;
  dijkstra_queue[0].extend=0;
  dijkstra_queue[0].revisit=0;
  dijkstra_queue[0].orientation=initial_orientation;

  char tmp_seq[db_graph->kmer_size+1];

  BinaryKmer local_copy_of_kmer, tmp_kmer;
  binary_kmer_assignment_operator(local_copy_of_kmer, source->kmer);
  BinaryKmer* rev_kmer = binary_kmer_reverse_complement(&local_copy_of_kmer,db_graph->kmer_size, &tmp_kmer);
  if (initial_orientation == forward) {
    binary_kmer_to_seq(&local_copy_of_kmer,db_graph->kmer_size,tmp_seq);
  }
  else {
    binary_kmer_to_seq(rev_kmer,db_graph->kmer_size,tmp_seq);
  }
  strcpy(dijkstra_queue[0].contig, tmp_seq);

  short k;
  for (k=0; k < 100; k++) 
    {
      dijkstra_queue[0].visited[k]=NULL; 
      dijkstra_queue[0].orien[k]=0;
      dijkstra_queue[0].passed[k]=NULL;
      dijkstra_queue[0].dist[k]=0;
    }

  if (output == NULL)
    {
      die("Couldn't find output file to write contigs\n\n");
    }

  boolean first_visit0=true;
  dBNode* previous_dB_node=NULL;
  run_along_a_supernode(&dijkstra_queue[0], previous_dB_node, &first_visit0, true);

  if (!first_visit0)
    {
      print_one_contig(dijkstra_queue[0].contig, db_graph->kmer_size, output, count);
      free(dijkstra_queue);
      return;
    }


  int begin_ptr=0, end_ptr=0;
  while (begin_ptr != (end_ptr+1) % 60000)
    {
      //add new elements to the queue
      queue_node curr_node;
      assign_queue_node(&curr_node, dijkstra_queue[begin_ptr]);
      /*
      if ((end_ptr >= 59996) && (begin_ptr > 0))
	{
	  int i;
	  for (i=begin_ptr; i <= end_ptr; i++)
	    {
	      assign_queue_node(&dijkstra_queue[i-begin_ptr], dijkstra_queue[i]);
	    }

	  end_ptr -= begin_ptr;
	  begin_ptr = 0;
	}
      */
      short j;
      for (j=0; j<4; j++)
	{
	  if (db_node_edge_exist(curr_node.ptr, j, curr_node.orientation, index))
	    {
	      Orientation tmp_orientation=0;
	      Nucleotide reverse_edge=0;
	      
	      previous_dB_node=curr_node.ptr;
	      dBNode* tmp_dB_node=db_graph_get_next_node_for_specific_person_or_pop(curr_node.ptr,curr_node.orientation,&tmp_orientation,j,&reverse_edge,db_graph,index);
	      char curr_char[2];
	      curr_char[0] = binary_nucleotide_to_char(j);
	      curr_char[1] = '\0';
	      
	      boolean first_visit=true;
	      for (k=0; k < 100; k++)
		{
		  if ((tmp_dB_node == curr_node.visited[k])&&(tmp_orientation == curr_node.orien[k])) {first_visit=false; break;}
		}
	      
	      unsigned short distal_length=strlen(curr_node.contig);
	      if (first_visit && if_two_branches_are_connected(tmp_dB_node,curr_node.passed,curr_node.dist,&distal_length,read_length,insert_size,min_sup,db_graph,link_list))
		{
		  //printf("%d\n", begin_ptr);
		  queue_node tmp_queue_node;
		  assign_queue_node(&tmp_queue_node, curr_node);
		  tmp_queue_node.ptr=tmp_dB_node;
		  tmp_queue_node.orientation=tmp_orientation;
		  strncat(tmp_queue_node.contig,curr_char,1);

		  if (! db_node_check_status(tmp_dB_node, visited))
		    {
		      db_node_set_status(tmp_dB_node, visited);
		      tmp_queue_node.revisit=0;
		    }
		  else tmp_queue_node.revisit++;

		  if (strlen(tmp_queue_node.contig) >= 49999)
		    {
		      first_visit = false;
		      //printf("string length limit reached\n\n");
		    }
		  
		  run_along_a_supernode(&tmp_queue_node, previous_dB_node, &first_visit, true);
		  print_one_contig(tmp_queue_node.contig, distal_length, output, count);
		  
		  if (first_visit && db_node_get_degree(tmp_queue_node.ptr,index,tmp_queue_node.orientation) > 1 && \
		      tmp_queue_node.revisit <= 2000)
		    {
		      if (db_node_get_degree(tmp_queue_node.ptr,index,opposite_orientation(tmp_queue_node.orientation)) > 1)
			{
			  add_visited_node(&tmp_queue_node, previous_dB_node, &first_visit);
			}
		      
		      end_ptr=(end_ptr+1) % 60000;
		      if (end_ptr == begin_ptr) 
			{
			  //int len=strlen(tmp_queue_node.contig);
			  //printf("Too many sequences starting from kmer %s\nSkipping the rest and continue with other kmers...\nCurrent length:%d\n\n", tmp_seq, len);
			  free(dijkstra_queue);
			  return;
			}

		      if (first_visit)
			assign_queue_node(&dijkstra_queue[end_ptr], tmp_queue_node);
		    }
		  
		}//if connected

	      else if (first_visit) 
		{
		  short num_branch_candidate=0;
		  short k=0;
		  while ((curr_node.visited[k] != NULL) && (k < 100)) 
		    {
		      int dist_between_branches=strlen(curr_node.contig)-curr_node.dist[k];
		      if (dist_between_branches < read_length || (dist_between_branches < insert_size && \
								  dist_between_branches > insert_size-2*read_length)) 
			num_branch_candidate++;

		      k++;
		    }
		  k--;

		  if (k < 100)
		    {
		      distal_length=curr_node.dist[k];
		      queue_node tmp_queue_node;
		      assign_queue_node(&tmp_queue_node, curr_node);
		      tmp_queue_node.ptr=tmp_dB_node;
		      tmp_queue_node.orientation=tmp_orientation;
		      strncat(tmp_queue_node.contig,curr_char,1);

		      if (strlen(tmp_queue_node.contig) >= 49999)
			{
			  first_visit = false;
			  //printf("string length limit reached\n\n");
			}
		      
		      //run_along_a_supernode(&tmp_queue_node, previous_dB_node, &first_visit, false);
		      //print_one_contig(tmp_queue_node.contig, distal_length, output, count);

		      //if (tmp_queue_node.extend > 15 || num_branch_candidate > 2 || strlen(curr_node.contig)-distal_length > 1.25 * insert_size)
		      if (tmp_queue_node.extend > 6 || num_branch_candidate > 1)
			{
			  run_along_a_supernode(&tmp_queue_node, previous_dB_node, &first_visit, false);
			  print_one_contig(tmp_queue_node.contig, distal_length, output, count);
			}
		      else
			{
			  if (! db_node_check_status(tmp_dB_node, visited))
			    {
			      //db_node_set_status(tmp_dB_node, visited);
			      tmp_queue_node.revisit=0;
			    }
			  else tmp_queue_node.revisit++;

                          run_along_a_supernode(&tmp_queue_node, previous_dB_node, &first_visit, true);
                          print_one_contig(tmp_queue_node.contig, distal_length, output, count);

			  if (first_visit && db_node_get_degree(tmp_queue_node.ptr,index,tmp_queue_node.orientation) > 1)
			    {
			      if (db_node_get_degree(tmp_queue_node.ptr,index,opposite_orientation(tmp_queue_node.orientation)) > 1)
				{
				  add_visited_node(&tmp_queue_node, previous_dB_node, &first_visit);
				}
			      
			      end_ptr=(end_ptr+1) % 60000;
			      if (end_ptr == begin_ptr)
				{
				  //int len=strlen(tmp_queue_node.contig);
				  //printf("Too many sequences starting from kmer %s\nSkipping the rest and continue with other kmers...\nCurrent length:%d\n\n", tmp_seq, len);
				  free(dijkstra_queue);
				  return;
				}

			      if (first_visit) 
				{
				  tmp_queue_node.extend++;
				  assign_queue_node(&dijkstra_queue[end_ptr], tmp_queue_node);
				}
			    }
			}//else 

		    }
		}//not connected

	    }//if edge_exist
	}//for j

      begin_ptr = (begin_ptr+1) % 60000;
    }//while

  free(dijkstra_queue);
}



//No longer in use
/*
void assemble_contigs_from_a_single_source_2(dBNode* source, Orientation initial_orientation, int read_length, int insert_size, int min_contig_length, int min_sup,
					     dBGraph* db_graph, LinkHashTable* link_list, unsigned int* count, FILE* output, int index)
{
  typedef struct
  {
    dBNode* ptr;
    Orientation orientation;
    Orientation orien[800];
    unsigned short revisit;
    dBNode* passed[800];
    dBNode* visited[800];
    unsigned short dist[800];
    char contig[30000];
  } queue_node;

  void print_one_contig(char* seq, unsigned short start, FILE* filename, unsigned int* num)
  {
    char* str_to_print = seq;

    str_to_print += start - db_graph->kmer_size;
    int length=strlen(str_to_print);
    //    printf("%d\n", length);
    //    printf("%s\n%s\n", seq, str_to_print);

    if (length >= min_contig_length) 
      {
	(*num)++;
	fprintf(filename, "%s%d\t%d\n", ">Contig",*num,length);
	fprintf(filename, "%s\n", str_to_print);
      }

    return;
  }

  void assign_queue_node(queue_node* node1, queue_node node2)
  {
    node1->ptr = node2.ptr;
    node1->revisit = node2.revisit;
    node1->orientation = node2.orientation;
    strcpy(node1->contig, node2.contig);

    short k;
    for (k=0; k < 800; k++) {
      node1->visited[k] = node2.visited[k];
      node1->orien[k]   = node2.orien[k];
      node1->passed[k]  = node2.passed[k];
      node1->dist[k]    = node2.dist[k];
    }
  }

  boolean add_visited_node(queue_node* node, dBNode* previous)
  {
    short visited_ptr=0;
    while ((visited_ptr < 800) && (node->passed[visited_ptr] != NULL)) {visited_ptr++;}
    if (visited_ptr < 800) {
      node->visited[visited_ptr] = node->ptr;
      node->orien[visited_ptr]   = node->orientation;
      node->passed[visited_ptr]  = previous;
      node->dist[visited_ptr]    = strlen(node->contig);
      return true;
    }
    //else {die("Too many visited nodes\nSet longer array\n\n");}
    else {return false;}
  }

  void run_along_a_supernode(queue_node* node, dBNode* previous_dB_node, boolean* first_visit, boolean if_set_visited)
  {
    while ((db_node_get_degree(node->ptr,index,node->orientation) == 1) && (*first_visit))
      {
        if (db_node_get_degree(node->ptr,index,opposite_orientation(node->orientation)) > 1)
          {
            *first_visit=add_visited_node(node, previous_dB_node);
          }

        dBNode* next=NULL;
        Orientation next_orientation;
        Nucleotide reverse_edge;

	char tmp_char[2];
        short l;
        for (l=0; l < 4; l++)
          {
            if (db_node_edge_exist(node->ptr, l, node->orientation, index))
              {
                next=db_graph_get_next_node_for_specific_person_or_pop(node->ptr,node->orientation,&next_orientation,l,&reverse_edge,db_graph,index);
		tmp_char[0] = binary_nucleotide_to_char(l);
                tmp_char[1] = '\0';
                break;
              }
          }

        short k;
        for (k=0; k < 800; k++)
          {
            if ((next == node->visited[k])&&(next_orientation == node->orien[k])) {*first_visit=false; break;}
          }

        if (*first_visit)
          {
            previous_dB_node = node->ptr;
            node->ptr = next;
            node->orientation = next_orientation;
            if (db_node_check_status(node->ptr, visited)) {node->revisit++;}
	    else 
	      {
		node->revisit=0;
		if (if_set_visited) {db_node_set_status(node->ptr, visited);}
	      }
	    strncat(node->contig, tmp_char, 1);
	    if (strlen(node->contig) >= 30000) {(*first_visit) = false;}
          }//if first_visit
      }//while

    if (db_node_get_degree(node->ptr,index,node->orientation) == 0)
      {
	(*first_visit) = false;
      }

    return;
  }


  queue_node* dijkstra_queue = calloc(50000,sizeof(queue_node));

  if (dijkstra_queue == NULL) {die("Could not allocate a dijkstra queue with 50000 elements\n\n");}

  db_node_set_status(source, visited);

  dijkstra_queue[0].ptr=source;
  dijkstra_queue[0].orientation=initial_orientation;
  dijkstra_queue[0].revisit=0;

  char tmp_seq[db_graph->kmer_size+1];

  BinaryKmer local_copy_of_kmer, tmp_kmer;
  binary_kmer_assignment_operator(local_copy_of_kmer, source->kmer);
  BinaryKmer* rev_kmer = binary_kmer_reverse_complement(&local_copy_of_kmer,db_graph->kmer_size, &tmp_kmer);
  if (initial_orientation == forward) {
    binary_kmer_to_seq(&local_copy_of_kmer,db_graph->kmer_size,tmp_seq);
  }
  else {
    binary_kmer_to_seq(rev_kmer,db_graph->kmer_size,tmp_seq);
  }
  strcpy(dijkstra_queue[0].contig, tmp_seq);

  short k;
  for (k=0; k < 800; k++)
    {
      dijkstra_queue[0].visited[k]=NULL;
      dijkstra_queue[0].orien[k]=0;
      dijkstra_queue[0].passed[k]=NULL;
      dijkstra_queue[0].dist[k]=0;
    }

  if (output == NULL)
    {
      die("Couldn't find output file to write contigs\n\n");
    }

  boolean first_visit0=true;
  dBNode* previous_dB_node=NULL;
  run_along_a_supernode(&dijkstra_queue[0], previous_dB_node, &first_visit0, true);

  if (!first_visit0)
    {
      print_one_contig(dijkstra_queue[0].contig, db_graph->kmer_size, output, count);
      free(dijkstra_queue);
      return;
    }


  int begin_ptr=0, end_ptr=0;
  while (begin_ptr <= end_ptr)
    {
      //add new elements to the queue                                                                                                                                       
      queue_node curr_node;
      assign_queue_node(&curr_node, dijkstra_queue[begin_ptr]);
      begin_ptr++;

      if ((end_ptr >= 59996) && (begin_ptr > 0))
        {
          int i;
          for (i=begin_ptr; i <= end_ptr; i++)
            {
              assign_queue_node(&dijkstra_queue[i-begin_ptr], dijkstra_queue[i]);
            }

          end_ptr -= begin_ptr;
          begin_ptr = 0;
        }

      short j;
      for (j=0; j<4; j++)
        {
          if (db_node_edge_exist(curr_node.ptr, j, curr_node.orientation, index))
            {
              Orientation tmp_orientation=0;
              Nucleotide reverse_edge=0;

	      previous_dB_node=curr_node.ptr;
              dBNode* tmp_dB_node=db_graph_get_next_node_for_specific_person_or_pop(curr_node.ptr,curr_node.orientation,&tmp_orientation,j,&reverse_edge,db_graph,index);
              char curr_char[2];
              curr_char[0] = binary_nucleotide_to_char(j);
              curr_char[1] = '\0';

              boolean first_visit=true;
              for (k=0; k < 800; k++)
                {
                  if ((tmp_dB_node == curr_node.visited[k])&&(tmp_orientation == curr_node.orien[k])) {first_visit=false; break;}
                }

	      //if (db_node_check_status(tmp_dB_node,visited))
	      //{
	      //  first_visit = false;
	      //}

	      unsigned short distal_length = strlen(curr_node.contig);
              if (first_visit && if_two_branches_are_connected(tmp_dB_node,curr_node.passed,curr_node.dist,&distal_length,read_length,insert_size,min_sup,db_graph,link_list))
                {
                  queue_node tmp_queue_node;
                  assign_queue_node(&tmp_queue_node, curr_node);
                  tmp_queue_node.ptr=tmp_dB_node;
                  tmp_queue_node.orientation=tmp_orientation;
                  strcat(tmp_queue_node.contig,curr_char);

                  if (strlen(tmp_queue_node.contig) >= 30000)
                    {
                      first_visit = false;
                    }

                  run_along_a_supernode(&tmp_queue_node, previous_dB_node, &first_visit, true);
		  print_one_contig(tmp_queue_node.contig, distal_length, output, count);

                  if (db_node_check_status(tmp_dB_node, visited))
                    {
                      tmp_queue_node.revisit++;
                      if (tmp_queue_node.revisit > (int) (2*insert_size)) {continue;}
                    }
                  else {db_node_set_status(tmp_dB_node, visited); tmp_queue_node.revisit=0;} 

                  if (first_visit && db_node_get_degree(tmp_queue_node.ptr,index,tmp_queue_node.orientation) > 1)
                    {
                      if (db_node_get_degree(tmp_queue_node.ptr,index,opposite_orientation(tmp_queue_node.orientation)) > 1)
                        {
                          add_visited_node(&tmp_queue_node, previous_dB_node);
			}

                      end_ptr++;
                      if (end_ptr > 59999)
                        {
                          //  die("Overflow of dijkstra_queue\nSet longer queue\n\n");
                          char curr_seq[db_graph->kmer_size+1];
                          BinaryKmer* tmp_kmer=element_get_kmer(source);
                          binary_kmer_to_seq(tmp_kmer,db_graph->kmer_size,curr_seq);
                          printf("Too many sequences starting from kmer %s\nSkipping the rest and continue with other kmers...\n\n", curr_seq);
                          free(dijkstra_queue);
                          return;
                        }

                      assign_queue_node(&dijkstra_queue[end_ptr], tmp_queue_node);
                    }

                }//if connected

	      else if (first_visit)
                {
                  short k=0;
                  while ((curr_node.passed[k] != NULL) && (k < 800)) {k++;}
                  k--;
                  if (k < 800)
                    {
                      distal_length=curr_node.dist[k];
                      queue_node tmp_queue_node;
                      assign_queue_node(&tmp_queue_node, curr_node);
                      tmp_queue_node.ptr=tmp_dB_node;
                      tmp_queue_node.orientation=tmp_orientation;
                      strcat(tmp_queue_node.contig,curr_char);

		      //boolean progress=(!db_node_check_status(tmp_dB_node, visited));

		      if (strlen(tmp_queue_node.contig) >= 30000)
                        {
                          first_visit = false;
                        }

                      run_along_a_supernode(&tmp_queue_node, previous_dB_node, &first_visit, false);
		      print_one_contig(tmp_queue_node.contig, distal_length, output, count);

                      if (first_visit && progress && db_node_get_degree(tmp_queue_node.ptr,index,tmp_queue_node.orientation) > 1)
                        {
                          if (db_node_get_degree(tmp_queue_node.ptr,index,opposite_orientation(tmp_queue_node.orientation)) > 1)
                            {
                              add_visited_node(&tmp_queue_node, previous_dB_node);
                            }

                          end_ptr++;
			  if (end_ptr > 59999)
                            {
			      char curr_seq[db_graph->kmer_size+1];
                              BinaryKmer* tmp_kmer=element_get_kmer(source);
                              binary_kmer_to_seq(tmp_kmer,db_graph->kmer_size,curr_seq);
                              printf("Too many sequences starting from kmer %s\nSkipping the rest and continue with other kmers...\n\n", curr_seq);
                              free(dijkstra_queue);
                              return;
                            }

                          assign_queue_node(&dijkstra_queue[end_ptr], tmp_queue_node);
                        }

                    }
		}//else if

            }//if edge_exist
        }//for j
    }//while

  free(dijkstra_queue);
}
*/
