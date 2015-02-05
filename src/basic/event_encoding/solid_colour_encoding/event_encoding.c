/* 
  solid_colour_encoding/event_encoding.c
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
	case '0':
	  return Zero;
	case '1':
	  return One;
	case '2':
	  return Two;
	case '3':
	  return Three;
	default:
	  return Undefined;
	}
}




//this method does basically nothing
Nucleotide reverse_binary_nucleotide(Nucleotide n)
{
  switch (n)
    {
    case Zero:
      return Zero;
    case One:
      return One;
    case Two:
      return Two;
    case Three:
      return Three;
    default:
      die("Calling reverse_binary_nucleotide on non-existent nucleotide %i", n);
    }
}

char binary_nucleotide_to_char(Nucleotide n)
{
	switch (n) {
	case Zero:
	  return '0';
	case One:
	  return '1';
	case Two:
	  return '2';
	case Three:
	  return '3';
	default:
	  die("Non existent binary nucleotide %d\n",n);
	}
}


boolean good_symbol(char c){
  boolean ret;
  if (c  != '0' && c != '1' && 
      c != '2' && c != '3' && 
      c != 'N' && c != 'n' 
      ){
    ret = false;
  }	
  else{
    ret =  true;
  }
  
  return ret;
}
