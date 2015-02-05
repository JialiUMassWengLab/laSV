/*
  global.c
*/

#include <stdlib.h>
#include <stdarg.h> // needed for va_list
#include <stdio.h>
#include <string.h>

#include "global.h"

boolean test_file_existence(char* file)
{
  FILE* fp = fopen(file, "r");
  if(fp == NULL)
  {
    return false;
  }
  else
  {
    fclose(fp);
    return true;
  }
}

// integer comparison: returns negative if a < b
//                                    0 if a == b
//                             positive if a > b
int int_cmp(const void *a, const void *b)
{
  // casting pointer types
  const int *ia = (const int *)a;
  const int *ib = (const int *)b;

  return (*ia  - *ib);
}

void set_string_to_null(char* str, int len)
{
  memset(str, 0, sizeof(char)*len);
}


void die(const char* fmt, ...)
{
  fflush(stdout);

  // Print error
  fprintf(stderr, "Error: ");

  va_list argptr;
  va_start(argptr, fmt);
  vfprintf(stderr, fmt, argptr);
  va_end(argptr);

  // Check if we need to print a newline
  if(*(fmt+strlen(fmt)-1) != '\n')
  {
    fprintf(stderr, "\n");
  }

  exit(EXIT_FAILURE);
}

void warn(const char* fmt, ...)
{
  fflush(stdout);

  // Print warning
  fprintf(stderr, "Warning: ");

  va_list argptr;
  va_start(argptr, fmt);
  vfprintf(stderr, fmt, argptr);
  va_end(argptr);

  // Check if we need to print a newline
  if(*(fmt+strlen(fmt)-1) != '\n')
  {
    fprintf(stderr, "\n");
  }
}

// Placeholder for message() -- a function for standard output
void message(const char* fmt, ...)
{
  va_list argptr;
  va_start(argptr, fmt);
  vfprintf(stderr, fmt, argptr);
  va_end(argptr);
}
