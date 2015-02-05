/*
 seq_file.c
 project: seq_file
 url: https://github.com/noporpoise/seq_file
 author: Isaac Turner <turner.isaac@gmail.com>
 Copyright (C) 20-June-2012
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h> // tolower

#include "sam.h"
#include "string_buffer.h"

#include "seq_file.h"
#include "seq_common.h"

#include "seq_fasta.h"
#include "seq_fastq.h"
#include "seq_plain.h"
#include "seq_sam.h"

const char* seq_file_types[6]
  = {"Unknown", "FASTA", "FASTQ", "Plain", "SAM", "BAM"};

const char* seq_file_types_zipped[6]
  = {"Unknown(zipped)", "FASTA(zipped)", "FASTQ(zipped)", "Plain(zipped)",
     "SAM", "BAM"};

//
// Sequence File reader
//

void _str_to_lower(char *str)
{
  for(; *str != '\0'; str++)
    *str = (char)tolower(*str);
}

char _contains_extension(const char *path, const char *ext, const size_t len)
{
  const char *tmp;

  for(tmp = path; (tmp = strstr(tmp, ext)) != NULL; tmp++)
  {
    // Check extension is followed by end-of-string or separator . or _
    if(*(tmp+len) == '\0' || *(tmp+len) == '.' || *(tmp+len) == '_')
    {
      return 1;
    }
  }

  return 0;
}

void seq_guess_filetype_from_path(const char *path, SeqFileType *file_type,
                                  char *zipped)
{
  #define ARRLEN 20

  const char* exts[ARRLEN] = {".fa",".fasta",
                              ".fq",".fastq",
                              ".faz",".fagz",".fa.gz",".fa.gzip",".fasta.gzip",
                              ".fqz",".fqgz",".fq.gz",".fq.gzip",".fastq.gzip",
                              ".txt",".txtgz",".txt.gz",".txt.gzip",
                              ".sam", ".bam"};

  const SeqFileType types[ARRLEN]
    = {SEQ_FASTA, SEQ_FASTA,
       SEQ_FASTQ, SEQ_FASTQ,
       SEQ_FASTA, SEQ_FASTA, SEQ_FASTA, SEQ_FASTA, SEQ_FASTA,
       SEQ_FASTQ, SEQ_FASTQ, SEQ_FASTQ, SEQ_FASTQ, SEQ_FASTQ,
       SEQ_PLAIN, SEQ_PLAIN, SEQ_PLAIN, SEQ_PLAIN,
       SEQ_SAM, SEQ_BAM};

  const char zips[ARRLEN] = {0, 0,
                             0, 0,
                             1, 1, 1, 1, 1,
                             1, 1, 1, 1, 1,
                             0, 1, 1, 1,
                             0, 0};

  int i;
  size_t strlens[ARRLEN];

  for(i = 0; i < ARRLEN; i++)
    strlens[i] = strlen(exts[i]);

  *file_type = SEQ_UNKNOWN;

  size_t path_len = strlen(path);

  char* strcpy = strdup(path);
  _str_to_lower(strcpy);

  // Check for an extension at the end
  for(i = 0; i < ARRLEN; i++)
  {
    if(path_len >= strlens[i] &&
       strcasecmp(path + path_len - strlens[i], exts[i]) == 0)
    {
      *file_type = types[i];
      *zipped = zips[i];
      free(strcpy);
      return;
    }
  }

  // Check for an extension anywhere in the filename
  for(i = 0; i < ARRLEN; i++)
  {
    if(_contains_extension(path, exts[i], strlens[i]))
    {
      *file_type = types[i];
      *zipped = zips[i];
      free(strcpy);
      return;
    }
  }

  free(strcpy);
}

// Determines file type and opens necessary streams + mallocs memory
void _set_seq_filetype(SeqFile *sf)
{
  // Guess filetype from path
  SeqFileType file_type;
  char zipped;

  seq_guess_filetype_from_path(sf->path, &file_type, &zipped);

  if(file_type == SEQ_SAM)
  {
    // SAM
    sf->sam_file = samopen(sf->path, "r", 0);
    sf->file_type = SEQ_SAM;
    sf->bam = bam_init1();
    return;
  }
  else if(file_type == SEQ_BAM)
  {
    // BAM
    sf->sam_file = samopen(sf->path, "rb", 0);
    sf->file_type = SEQ_BAM;
    sf->bam = bam_init1();
    return;
  }

  // If not SAM or BAM, we can open it and determine its contents -
  // more reliable

  // Open file for the first time
  if(strcmp(sf->path, "-") == 0)
  {
    sf->gz_file = gzdopen(fileno(stdin), "r");
  }
  else
  {
    sf->gz_file = gzopen(sf->path, "r");
  }

  if(sf->gz_file == NULL)
  {
    fprintf(stderr, "Error: Couldn't open gz_file [%s:%i]\n", __FILE__, __LINE__);
    return;
  }

  #ifdef ZLIB_VERNUM
    #if (ZLIB_VERNUM >= 0x1240)
      // Set buffer size to 1Mb
      gzbuffer(sf->gz_file, (unsigned int)1024*1024);
    #endif
  #endif

  int first_char;

  // Move sf->line_number from 0 to 1 on first character
  // Then for each newline, line_number++

  do
  {
    first_char = gzgetc(sf->gz_file);
    sf->line_number++;
  } while (first_char != -1 && (first_char == '\n' || first_char == '\r'));

  if(first_char == -1)
  {
    fprintf(stderr, "seq_file.c warning: empty sequence file\n");
    return;
  }
  else if(first_char == '>')
  {
    // Reading FASTA
    sf->file_type = SEQ_FASTA;
    sf->read_line_start = 1;
  }
  else if(first_char == '@')
  {
    // Reading FASTQ
    sf->file_type = SEQ_FASTQ;
    sf->read_line_start = 1;
    sf->bases_buff = strbuf_new();
  }
  else if(is_base_char(first_char))
  {
    // Plain file
    sf->file_type = SEQ_PLAIN;
    sf->read_line_start = 0;

    if(gzungetc(first_char, sf->gz_file) == -1)
    {
      fprintf(stderr, "seq_file.c error: gzungetc failed\n");
      exit(EXIT_FAILURE);
    }
  }
  else
  {
    fprintf(stderr, "seq_file.c error: unknown filetype starting '%c'\n",
            first_char);
  }
}

SeqFile* _create_default_seq_file(const char* file_path)
{
  SeqFile* sf = (SeqFile*) malloc(sizeof(SeqFile));

  sf->path = file_path;

  sf->gz_file = NULL;
  sf->sam_file = NULL;

  sf->bam = NULL;

  sf->fastq_ascii_offset = 33;

  sf->file_type = SEQ_UNKNOWN;
  sf->read_line_start = 0;

  sf->entry_name = strbuf_new();
  sf->entry_index = 0;

  sf->entry_offset = 0;
  sf->entry_offset_qual = 0;

  sf->entry_read = 0;
  sf->entry_read_qual = 0;

  sf->bases_buff = NULL;

  sf->total_bases_passed = 0;
  sf->total_bases_skipped = 0;

  sf->line_number = 0;

  // For writing
  sf->plain_file = NULL;
  sf->line_wrap = 0;
  sf->curr_line_length = 0;
  sf->write_state = WS_READ_ONLY;

  return sf;
}

SeqFile* seq_file_open(const char* file_path)
{
  SeqFile* sf = _create_default_seq_file(file_path);
  _set_seq_filetype(sf);

  if(sf->file_type == SEQ_UNKNOWN)
  {
    seq_file_close(sf);
    return NULL;
  }
  else
  {
    return sf;
  }
}

SeqFile* seq_file_open_filetype(const char* file_path,
                                SeqFileType file_type)
{
  SeqFile* sf = _create_default_seq_file(file_path);
  sf->file_type = file_type;

  if(file_type == SEQ_FASTQ)
  {
    sf->bases_buff = strbuf_new();
  }

  switch(file_type)
  {
    case SEQ_SAM:
      sf->sam_file = samopen(sf->path, "r", 0);
      sf->bam = bam_init1();
      break;
    case SEQ_BAM:
      sf->sam_file = samopen(sf->path, "rb", 0);
      sf->bam = bam_init1();
      break;
    case SEQ_FASTA:
    case SEQ_FASTQ:
    case SEQ_PLAIN:
      if(strcmp(sf->path, "-") == 0)
      {
        sf->gz_file = gzdopen(fileno(stdin), "r");
      }
      else
      {
        sf->gz_file = gzopen(sf->path, "r");
      }
      break;
    default:
      fprintf(stderr, "seq_file.c warning: invalid SeqFileType in "
                      "function seq_file_open_filetype()\n");
      free(sf);
      return NULL;
  }

  return sf;
}

// Open to write
// file_type must be FASTA, FASTQ
// set gzip to != 0 to turn on gzipping output
// if line_wrap is != 0, sequence lines are wrapped
SeqFile* seq_file_open_write(const char* file_path, SeqFileType file_type,
                             char gzip, unsigned long line_wrap)
{
  if(file_type == SEQ_SAM || file_type == SEQ_BAM)
  {
    fprintf(stderr, "seq_file.c error: cannot write to a SAM or BAM file\n");
    return NULL;
  }
  else if(file_type == SEQ_UNKNOWN)
  {
    fprintf(stderr, "seq_file.c error: cannot open file type SEQ_UNKNOWN\n");
    return NULL;
  }
  else if(file_type == SEQ_PLAIN && line_wrap != 0)
  {
    fprintf(stderr, "seq_file.c warning: cannot set line wrap with 'plain' "
                    "sequence format\n");
    line_wrap = 0;
  }

  gzFile* gz_file = NULL;
  FILE* plain_file = NULL;

  if(gzip)
  {
    if((gz_file = gzopen(file_path, "w")) == NULL)
    {
      return NULL;
    }
  }
  else
  {
    if((plain_file = fopen(file_path, "w")) == NULL)
    {
      return NULL;
    }
  }

  SeqFile* sf = _create_default_seq_file(file_path);

  sf->gz_file = gz_file;
  sf->plain_file = plain_file;
  sf->file_type = file_type;
  sf->line_wrap = line_wrap;

  // Set write state (default is read-only)
  sf->write_state = WS_BEGIN;

  return sf;
}


// Close
// Adds an extra byte (newline) to the end of output (writing) files
size_t seq_file_close(SeqFile* sf)
{
  size_t num_bytes_printed = 0;

  // Add a new line to the end of output files
  if(sf->write_state != WS_READ_ONLY)
  {
    num_bytes_printed = seq_puts(sf, "\n");
    sf->line_number++;
  }

  if(sf->gz_file != NULL)
  {
    gzclose(sf->gz_file);
  }

  if(sf->file_type == SEQ_SAM || sf->file_type == SEQ_BAM)
  {
    bam_destroy1(sf->bam);
    samclose(sf->sam_file);
  }

  if(sf->bases_buff != NULL)
  {
    strbuf_free(sf->bases_buff);
  }

  if(sf->plain_file != NULL)
  {
    fclose(sf->plain_file);
  }

  strbuf_free(sf->entry_name);

  free(sf);

  return num_bytes_printed;
}

SeqFileType seq_file_get_type(const SeqFile* sf)
{
  return sf->file_type;
}

const char* seq_file_get_type_str(const SeqFile* sf)
{
  return seq_file_types[sf->file_type];
}

const char* seq_file_type_str(SeqFileType file_type,
                              char zipped)
{
  return zipped ? seq_file_types_zipped[file_type] : seq_file_types[file_type];
}

// Get a pointer to the file path
const char* seq_get_path(const SeqFile* sf)
{
  return sf->path;
}

// Set FASTQ ASCII offset (also applies to SAM/BAM)
void seq_set_fastq_ascii_offset(SeqFile *sf, char fastq_ascii_offset)
{
  sf->fastq_ascii_offset = fastq_ascii_offset;
}

// Get FASTQ ASCII offset (also applies to SAM/BAM)
char seq_get_fastq_ascii_offset(const SeqFile *sf)
{
  return sf->fastq_ascii_offset;
}

// Get the number of bases read/written so far
unsigned long seq_total_bases_passed(const SeqFile *sf)
{
  return sf->total_bases_passed;
}

// Get the total bases skipped (not read through API) in file so far
unsigned long seq_total_bases_skipped(const SeqFile *sf)
{
  return sf->total_bases_skipped;
}

unsigned long seq_curr_line_number(const SeqFile *sf)
{
  return sf->line_number;
}

char seq_has_quality_scores(const SeqFile *sf)
{
  switch (sf->file_type)
  {
    case SEQ_SAM:
    case SEQ_BAM:
    case SEQ_FASTQ:
      return 1;
    case SEQ_FASTA:
    case SEQ_PLAIN:
      return 0;
    default:
      fprintf(stderr, "seq_file.c: unknown file type [path: %s]\n", sf->path);
      return 0;
  }
}

// Returns 1 if open for writing, 0 otherwise
char seq_is_open_for_write(const SeqFile *sf)
{
  return (sf->write_state != WS_READ_ONLY);
}

// Get the name of the next read
const char* seq_get_read_name(SeqFile *sf)
{
  return sf->entry_name->buff;
}

// Get this read index -- starts from 0
unsigned long seq_get_read_index(SeqFile *sf)
{
  return sf->entry_index;
}

unsigned long seq_get_bases_read(SeqFile *sf)
{
  return sf->entry_offset;
}

unsigned long seq_get_quals_read(SeqFile *sf)
{
  return sf->entry_offset_qual;
}

unsigned long seq_get_length(SeqFile *sf)
{
  if(sf->entry_read)
  {
    fprintf(stderr, "seq_file.c: haven't finished reading sequence - "
                    "seq_get_length() cannot return a length\n");
    return 0;
  }

  return sf->entry_offset;
}



// Returns 1 on success 0 if no more to read
char seq_next_read(SeqFile *sf)
{
  strbuf_reset(sf->entry_name);

  char success;

  switch (sf->file_type)
  {
    case SEQ_FASTA:
      success = seq_next_read_fasta(sf);
      break;
    case SEQ_FASTQ:
      success = seq_next_read_fastq(sf);
      break;
    case SEQ_SAM:
    case SEQ_BAM:
      success = seq_next_read_sam(sf);
      break;
    case SEQ_PLAIN:
      success = seq_next_read_plain(sf);
      break;
    default:
      fprintf(stderr, "seq_file.c: Cannot read from unknown file type "
                      "[file: %s]\n", sf->path);
      return 0;
  }

  sf->entry_offset = 0;
  sf->entry_offset_qual = 0;

  if(success)
  {
    sf->entry_index++;
    sf->entry_read = 1;
    sf->entry_read_qual = seq_has_quality_scores(sf);
  }
  else
  {
    sf->entry_read = 0;
    sf->entry_read_qual = 0;
  }

  return success;
}

/*
 Read a single base from a read
*/

// Read a single base from the current read
// Returns 1 on success, 0 if no more quality scores or run out of bases
char seq_read_base(SeqFile *sf, char *c)
{
  // Check if we have read anything in
  if(!sf->entry_read)
    return 0;

  char success;

  switch (sf->file_type)
  {
    case SEQ_FASTA:
      success = seq_read_base_fasta(sf, c);
      break;
    case SEQ_FASTQ:
      success = seq_read_base_fastq(sf, c);
      break;
    case SEQ_SAM:
    case SEQ_BAM:
      success = seq_read_base_sam(sf, c);
      break;
    case SEQ_PLAIN:
      success = seq_read_base_plain(sf, c);
      break;
    default:
      fprintf(stderr, "seq_file.c: Cannot read from unknown file type "
                      "[file: %s]\n", sf->path);
      return 0;
  }

  if(success)
  {
    sf->entry_offset++;
    sf->total_bases_passed++;
  }

  sf->entry_read = success;

  return success;
}


/*
 Read a single quality score from a read
*/

// Read a single quality score from the current read
// Returns 1 on success, 0 if no more quality scores or run out of bases
char seq_read_qual(SeqFile *sf, char *c)
{
  // Check if we have read anything in
  if(!sf->entry_read_qual)
    return 0;

  char success;

  switch (sf->file_type)
  {
    case SEQ_FASTQ:
      success = seq_read_qual_fastq(sf, c);
      break;
    case SEQ_SAM:
    case SEQ_BAM:
      success = seq_read_qual_sam(sf, c);
      break;
    case SEQ_FASTA:
    case SEQ_PLAIN:
      fprintf(stderr, "seq_file.c: Cannot read from file type without quality "
                      "scores [file: %s; type: %s]\n", sf->path,
                      seq_file_type_str(sf->file_type, 0));
      return 0;
    default:
      fprintf(stderr, "seq_file.c: Cannot read from unknown file type "
                      "[file: %s]\n", sf->path);
      return 0;
  }

  if(success)
    sf->entry_offset_qual++;

  sf->entry_read_qual = success;

  return success;
}

/*
 Read k bases / quality scores from a read
*/

// str must be at least k+1 bytes long.  Null-terminates str at position k+1.
// returns 1 on success, 0 otherwise
char seq_read_k_bases(SeqFile *sf, char* str, int k)
{
  // Check if we have read anything in
  if(!sf->entry_read)
    return 0;

  int i;

  for(i = 0; i < k; i++)
  {
    if(!seq_read_base(sf, str+i))
    {
      sf->entry_read = 0;
      str[0] = '\0';
      return 0;
    }
  }

  str[k] = '\0';
  return 1;
}

char seq_read_k_quals(SeqFile *sf, char* str, int k)
{
  // Check if we have read anything in
  if(!sf->entry_read_qual)
    return 0;

  int i;

  for(i = 0; i < k; i++)
  {
    if(!seq_read_qual(sf, str+i))
    {
      sf->entry_read_qual = 0;
      str[0] = '\0';
      return 0;
    }
  }

  str[k] = '\0';
  return 1;
}


/*
 Read the rest of a read
*/

// returns 1 on success, 0 otherwise
char seq_read_all_bases(SeqFile *sf, StrBuf *sbuf)
{
  // Check if we have read anything in
  if(!sf->entry_read)
    return 0;

  strbuf_reset(sbuf);

  char success;

  switch(sf->file_type)
  {
    case SEQ_SAM:
    case SEQ_BAM:
      success = seq_read_all_bases_sam(sf, sbuf);
      break;
    case SEQ_FASTA:
      success = seq_read_all_bases_fasta(sf, sbuf);
      break;
    case SEQ_FASTQ:
      success = seq_read_all_bases_fastq(sf, sbuf);
      break;
    case SEQ_PLAIN:
      success = seq_read_all_bases_plain(sf, sbuf);
      break;
    default:
      fprintf(stderr, "seq_file.c: tried to read from unknown filetype "
                      "[path: %s]\n", sf->path);
      return 0;
  }

  sf->entry_offset += strbuf_len(sbuf);
  sf->total_bases_passed += strbuf_len(sbuf);

  // read has been exhausted
  sf->entry_read = 0;

  return success;
}


/*
 Read the rest of a read's quality scores
*/

// returns 1 on success, 0 otherwise
char seq_read_all_quals(SeqFile *sf, StrBuf *sbuf)
{
  // Check if we have read anything in
  if(!sf->entry_read_qual)
    return 0;

  strbuf_reset(sbuf);

  char success;

  switch(sf->file_type)
  {
    case SEQ_SAM:
    case SEQ_BAM:
      success = seq_read_all_quals_sam(sf, sbuf);
      break;
    case SEQ_FASTQ:
      success = seq_read_all_quals_fastq(sf, sbuf);
      break;
  default:
      fprintf(stderr, "seq_file.c: tried to read from unknown filetype "
                      "[path: %s]\n", sf->path);
      return 0;
  }

  // Exhausted read quality scores
  sf->entry_offset_qual += strbuf_len(sbuf);
  sf->entry_read_qual = 0;

  return success;
}

/*
 Write to a file. 
 Each function returns the number of bytes written or 0 on failure
*/

size_t seq_file_write_name(SeqFile *sf, const char *name)
{
  size_t num_bytes_printed = 0;

  switch(sf->file_type)
  {
    case SEQ_FASTA:
      num_bytes_printed = seq_file_write_name_fasta(sf, name);
      break;
    case SEQ_FASTQ:
      num_bytes_printed = seq_file_write_name_fastq(sf, name);
      break;
    default:
      fprintf(stderr, "seq_file.c: called seq_file_write_name() with invalid "
                      "file type (%s) [path: %s]\n", seq_file_get_type_str(sf),
                      sf->path);
      return 0;
  }

  sf->write_state = WS_NAME;
  return num_bytes_printed;
}


/*
 Write sequence
*/

size_t seq_file_write_seq(SeqFile *sf, const char *seq)
{
  size_t str_len = strlen(seq);

  unsigned long num_bytes_printed = 0;

  switch(sf->file_type)
  {
    case SEQ_FASTA:
      num_bytes_printed = seq_file_write_seq_fasta(sf, seq, str_len);
      break;
    case SEQ_FASTQ:
      num_bytes_printed = seq_file_write_seq_fastq(sf, seq, str_len);
      break;
    case SEQ_PLAIN:
      num_bytes_printed = seq_file_write_seq_plain(sf, seq);
      break;
    default:
      fprintf(stderr, "seq_file.c: called seq_file_write_seq() with invalid "
                      "file type (%s) [path: %s]\n", seq_file_get_type_str(sf),
                      sf->path);
      return 0;
  }

  sf->total_bases_passed += str_len;

  sf->write_state = WS_SEQ;

  return num_bytes_printed;
}

// Print quality
// Only FASTQ file types are allowed to call this function
size_t seq_file_write_qual(SeqFile *sf, const char *qual)
{
  if(sf->file_type != SEQ_FASTQ)
  {
    fprintf(stderr, "seq_file.c: called seq_file_write_qual() with invalid "
                    "file type (%s) [path: %s]\n", seq_file_get_type_str(sf),
                    sf->path);
    return 0;
  }

  return seq_file_write_qual_fastq(sf, qual);
}
