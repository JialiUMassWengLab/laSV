/* 
  base_encoding/event_encoding.c
*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "event_encoding.h"

//returns Undefined if given non AGCT character
Nucleotide char_to_binary_nucleotide(char c)
{
	switch (c)
	{
	case 'A':
	  return Adenine;
	case 'C':
	  return Cytosine;
	case 'G':
	  return Guanine;
	case 'T':
	  return Thymine;
	case 'a':
	  return Adenine;
	case 'c':
	  return Cytosine;
	case 'g':
	  return Guanine;
	case 't':
	  return Thymine;
	default:
	  return Undefined;
	}
}





Nucleotide reverse_binary_nucleotide(Nucleotide n)
{
  switch (n)
    {
    case Adenine:
      return Thymine;
    case Cytosine:
      return Guanine;
    case Guanine:
      return Cytosine;
    case Thymine:
      return Adenine;
    default:
      die("Calling reverse_binary_nucleotide on non-existent nucleotide %i", n);
    }
}

char binary_nucleotide_to_char(Nucleotide n)
{
  switch (n)
  {
    case Adenine:
      return 'A';
    case Cytosine:
      return 'C';
    case Guanine:
      return 'G';
    case Thymine:
      return 'T';
    default:
      die("Non existent binary nucleotide %d\n",n);
  }
}




boolean good_symbol(char c){
  boolean ret;
  if (c  != 'A' && c != 'a' && 
      c != 'C' && c != 'c' && 
      c != 'G' && c != 'g' && 
      c != 'T' && c != 't' && 
      c != 'N' && c != 'n' 
      ){
    ret = false;
  }	
  else{
    ret =  true;
  }
  
  return ret;
}
