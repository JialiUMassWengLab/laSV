/*
  linkage_element.c - copy of element.c with modifications
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <inttypes.h>

// cortex_var headers
#include "linkage_element.h"

//currently noone calls this in normal use
// In normal use, the priority queue allocates space to put the eloement directly within,
// and calls genotyping_element_initialise
LinkageElement* new_linkage_element()
{
  LinkageElement* e = malloc(sizeof(LinkageElement));

  if (e==NULL)
  {
    die("Unable to allocate a new linkage_element");
  }
  
  binary_kmer_initialise_to_zero(&(e->kmer));
  e->link_num=0;
  e->status=1;
  e->link=NULL;
  
  return e;
}


void free_linkage_element(LinkageElement** element)
{
  LinkNode* current=(*element)->link;
  while (current != NULL)
    {
      LinkNode* next=current->next;
      free(current);
      current=next;
    }

  free(*element);
  *element=NULL;
}

void linkage_element_initialise_from_normal_element(LinkageElement* e1, 
						    Element* e2)
{
  if (e2==NULL)
    {
      e1=NULL;
      return;
    }

  binary_kmer_assignment_operator( (*e1).kmer, (*e2).kmer);
  e1->link_num=0;
  e1->status=1;
  e1->link=NULL;
}


void linkage_element_assign(LinkageElement* e1, LinkageElement* e2)
{
  if ((e1==NULL)||(e2==NULL))
    {
      return;
    }

  binary_kmer_assignment_operator( (*e1).kmer, (*e2).kmer);
  e1->link_num = e2->link_num;
  e1->status   = e2->status;

  LinkNode* current1 = NULL;
  LinkNode* current2 = e2->link;
  while (current2 != NULL)
    {
      LinkNode* new_link = (LinkNode*) malloc(sizeof(LinkNode));
      new_link->ptr  = current2->ptr;
      new_link->read = current2->read;
      new_link->next = NULL;

      if (current1 == NULL) e1->link=new_link;
      else {current1->next=new_link; current1 = new_link;}

      current2 = current2->next;
    }
}


boolean db_linkage_node_check_for_flag_ALL_OFF(LinkageElement * node) {
  if (node->status == 0) {
    return true;
  }
  else {
    return false;
  }
}


//WARNING - this gives you a pointer to a the binary kmer in the node. You could modify contents of the hash table
BinaryKmer* linkage_element_get_kmer(LinkageElement * e){
  return &(e->kmer);
}

boolean linkage_element_is_key(Key key, LinkageElement e){
  //  return key == e.kmer;
  return binary_kmer_comparison_operator(*key, e.kmer);
}

Key linkage_element_get_key(BinaryKmer* kmer, short kmer_size, Key preallocated_key){
  
  BinaryKmer local_rev_kmer;
  binary_kmer_initialise_to_zero(&local_rev_kmer);

  binary_kmer_reverse_complement(kmer,kmer_size, &local_rev_kmer);
  
  if (binary_kmer_less_than(local_rev_kmer,*kmer, kmer_size))
    {
      binary_kmer_assignment_operator(*((BinaryKmer*)preallocated_key),local_rev_kmer);
    }
  else
    {
      binary_kmer_assignment_operator(*((BinaryKmer*)preallocated_key),*kmer);
    }

  return preallocated_key;

}


void linkage_element_initialise(LinkageElement * e, Key kmer, short kmer_size){

  if (e==NULL)
    {
      die("Called elemtn_initialise on NULL ptr");
    }

  BinaryKmer tmp_kmer;
  binary_kmer_initialise_to_zero(&tmp_kmer);
  binary_kmer_assignment_operator( e->kmer, *(linkage_element_get_key(kmer, kmer_size, &tmp_kmer)));

  e->link_num=0;
  e->status=1;
  e->link=NULL;

}


void linkage_element_initialise_kmer_to_zero(LinkageElement * e)
{
  if (e==NULL)
  {
    die("Called genotyping_element_initialise_covgs_and_edges_to_zero on NULL ptr");
  }

  binary_kmer_initialise_to_zero(&(e->kmer));

  e->link_num=0;
  e->status=1;
  e->link=NULL;
}




void linkage_element_set_kmer(LinkageElement * e, Key kmer, short kmer_size)
{
  if (e==NULL)
  {
    die("Called element_set_kmer on NULL ptr");
  }

  BinaryKmer tmp_kmer;
  binary_kmer_initialise_to_zero(&tmp_kmer);
  binary_kmer_assignment_operator( e->kmer, *(linkage_element_get_key(kmer, kmer_size, &tmp_kmer)));
}



Orientation db_linkage_node_get_orientation(BinaryKmer* k, LinkageElement * e, short kmer_size){

  if (binary_kmer_comparison_operator(e->kmer,*k)==true)
    {
      return forward;
    }
  
  BinaryKmer tmp_kmer;

  if (binary_kmer_comparison_operator(e->kmer, *(binary_kmer_reverse_complement(k,kmer_size, &tmp_kmer)))==true)
    {
      return reverse;
    }
  
  char tmpseq1[kmer_size+1];
  char tmpseq2[kmer_size+1];

  die("programming error - you have called db_genotyping_node_get_orientation \n"
      "with a kmer that is neither equal to the kmer in this node, nor its rev comp\n"
      "Arg 1 Kmer is %s and Arg 2 node kmer is %s\n",
      binary_kmer_to_seq(k, kmer_size, tmpseq1),
      binary_kmer_to_seq(&(e->kmer), kmer_size, tmpseq2));  
}



