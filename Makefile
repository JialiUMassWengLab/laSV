#ifndef CC
  CC = gcc
#endif

BIN = bin
TEMP_TEST_DIR = data/tempfiles_can_be_deleted

ifeq ($(MAXK),31)
   BITFIELDS = 1
endif

ifeq ($(MAXK),63)
   BITFIELDS = 2
endif

ifeq ($(MAXK),95)
   BITFIELDS = 3
endif

ifeq ($(MAXK),127)
   BITFIELDS = 4
endif

ifeq ($(MAXK),160)
   BITFIELDS = 5
endif

ifeq ($(MAXK),192)
   BITFIELDS = 6
endif

ifeq ($(MAXK),223)
   BITFIELDS = 7
endif

ifeq ($(MAXK),255)
   BITFIELDS = 8
endif

ifndef BITFIELDS
   BITFIELDS = 1
   MAXK = 31
endif


ifndef NUM_COLS
  NUM_COLS = 1
endif

##ifeq ($(BITFIELDS),0)
#$(error Invalid value for MAXK - either omit or use 32*x-1)
#endif

# Test if running on a mac
UNAME=$(shell uname)
ifeq ($(UNAME),Darwin)
	MAC = 1
endif

# Library paths
IDIR_GSL = libs/gsl-1.15
IDIR_GSL_ALSO = libs/gsl-1.15/gsl
IDIR_STRS = libs/string_buffer
IDIR_SEQ = libs/seq_file
IDIR_BAM = libs/samtools-0.1.18

# Main program includes
IDIR_BASIC = include/basic
IDIR_BASE_ENCODING = ${IDIR_BASIC}/event_encoding/base_encoding
IDIR_COLOUR_ENCODING = ${IDIR_BASIC}/event_encoding/solid_colour_encoding
IDIR_HASH = include/hash_table
IDIR_OPERATION = include/dB_graph_operations
IDIR_OPERATION_CORE = include/dB_graph_operations/core
IDIR_OPERATION_CMD_LINE = include/dB_graph_operations/many_colours

# Correct for zam's paths for CUNIT
#NOT_ZAM=1

ifdef NOT_ZAM
	IDIR_CUNIT = /opt/local/include/CUnit
	LDIR_CUNIT = /opt/local/lib
else
	IDIR_CUNIT = /home/zam/dev/hg/CUnit/CUnit-2.1-0/CUnit/Headers

	ifdef MAC
		IDIR_CUNIT = /Users/zam/dev/hg/laptop_local/repos/CUnit/CUnit-2.1-0/CUnit/Headers
	else
		LDIR_CUNIT = /home/zam/bin/lib
	endif
endif

ifdef MAC
	MACFLAG = -fnested-functions
endif

ARCH = -m64

ifdef 32_BITS
	ARCH =
endif

# Comment out this line to turn off adding the commit version
# (it already checks if hg is installed)
VERSION_STR=$(shell if [ `command -v hg` ]; then echo ' (commit' `hg id --num --id`')'; else echo; fi)

# DEV: Add -DNDEBUG=1 to turn off assert() calls
OPT := $(ARCH) -Wall $(MACFLAG) -DVERSION_STR='"$(VERSION_STR)"' \
       -DNUMBER_OF_BITFIELDS_IN_BINARY_KMER=$(BITFIELDS) \
       -DNUMBER_OF_COLOURS=$(NUM_COLS)

ifdef DEBUG
	OPT := -O0 -g $(OPT)
else
	OPT := -O3 $(OPT)
endif

LIBLIST = -lseqfile -lbam -lstrbuf -lz -lm
TEST_LIBLIST = -lcunit -lncurses $(LIBLIST)

# Add -L/usr/local/lib/ to satisfy some systems that struggle to link libz
LIBINCS = -L/usr/local/lib -I$(IDIR_BAM) \
          -I$(IDIR_SEQ) -I$(IDIR_STRS) \
          -L$(IDIR_BAM) -L$(IDIR_SEQ) -L$(IDIR_STRS)

TEST_LIBINCS = -I$(IDIR_CUNIT) -L$(LDIR_CUNIT) $(LIBINCS)

CFLAGS_BASIC      = -I$(IDIR_BASIC) -I$(IDIR_BASE_ENCODING) $(LIBINCS)
CFLAGS_OPERATION     = -I$(IDIR_OPERATION_CORE) -I$(IDIR_BASIC) -I$(IDIR_HASH) -I$(IDIR_OPERATION) -I$(IDIR_BASE_ENCODING) $(LIBINCS)

OPERATION_OBJ = src/obj/dB_graph_operations/print_branches.o src/obj/dB_graph_operations/linkage_element.o src/obj/dB_graph_operations/little_hash_for_linkage.o src/obj/dB_graph_operations/model_info.o src/obj/dB_graph_operations/global.o src/obj/dB_graph_operations/graph_operations.o src/obj/dB_graph_operations/binary_kmer.o src/obj/dB_graph_operations/element.o src/obj/dB_graph_operations/seq.o src/obj/dB_graph_operations/hash_value.o src/obj/dB_graph_operations/hash_table.o src/obj/dB_graph_operations/dB_graph.o src/obj/dB_graph_operations/dB_graph_population.o src/obj/dB_graph_operations/cmd_line.o src/obj/dB_graph_operations/event_encoding.o src/obj/dB_graph_operations/graph_info.o src/obj/dB_graph_operations/maths.o src/obj/dB_graph_operations/file_reader.o src/obj/dB_graph_operations/linkage_operations.o src/obj/dB_graph_operations/build_linkage_from_files.o


MAXK_AND_TEXT = $(join "", $(MAXK))
NUMCOLS_AND_TEST = $(join "_c", $(NUM_COLS))

dB_graph_operations : $(OPERATION_OBJ)
	mkdir -p $(BIN); $(CC) $(CFLAGS_OPERATION) $(OPT) $(OPT_COLS) -o $(BIN)/dB_graph $(OPERATION_OBJ) $(LIBLIST)


.PHONY : clean
clean :
	rm -rf $(BIN)/*
	rm -rf src/obj

remove_objects:
	rm -rf src/obj/*







#pattern rules


src/obj/dB_graph_operations/%.o : src/dB_graph_operations/%.c include/dB_graph_operations/%.h
	mkdir -p src/obj/dB_graph_operations; $(CC) $(CFLAGS_OPERATION) $(OPT) -c $< -o $@

src/obj/dB_graph_operations/%.o : src/basic/event_encoding/base_encoding/%.c include/basic/event_encoding/base_encoding/%.h
	mkdir -p src/obj/dB_graph_operations; $(CC) $(CFLAGS_OPERATION) $(OPT) -c $< -o $@

src/obj/dB_graph_operations/%.o : src/basic/%.c include/basic/%.h
	mkdir -p src/obj/dB_graph_operations; $(CC) $(CFLAGS_OPERATION) $(OPT) -c $< -o $@

src/obj/dB_graph_operations/%.o : src/hash_table/hash_key/bob_jenkins/%.c include/hash_table/%.h
	mkdir -p src/obj/dB_graph_operations; $(CC) $(CFLAGS_OPERATION) $(OPT) -c $< -o $@

src/obj/dB_graph_operations/%.o : src/hash_table/open_hash/%.c include/hash_table/open_hash/%.h
	mkdir -p src/obj/dB_graph_operations; $(CC) $(CFLAGS_OPERATION) $(OPT) -c $< -o $@
