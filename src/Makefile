.DELETE_ON_ERROR:

SHELL         = /bin/bash

ROOTCONFIG   := root-config
ROOTCFLAGS   := $(shell $(ROOTCONFIG) --cflags)
ROOTLDFLAGS  := $(shell $(ROOTCONFIG) --ldflags)
ROOTLIBS     := $(shell $(ROOTCONFIG) --libs)
ROOTDICTTYPE := $(shell $(ROOTCONFIG) --dicttype)
ROOTCINT     := rootcint
ROOTLIBDIR  := $(shell $(ROOTCONFIG) --libdir)

CXX           = g++
CXXFLAGS      = -O2 -Wall -fPIC $(ROOTCFLAGS)
LD            = g++
LDFLAGS       = -O2 $(ROOTLDFLAGS)
SOFLAGS       = -shared
LIBS          = $(ROOTLIBS)

GENERATOR_OBJ     = TNtuple_generator.o
GENERATOR_SRC     = TNtuple_generator.cxx
GENERATOR_SO      = libTNtuple_generator.so
GENERATOR         = tntuple_generator
GENERATOR_LIBS    = -L$(ROOTLIBDIR) -lCore -lCint -lTree -lHist 

MAIN_OBJ = generatorMain.o
MAIN_SRC = generatorMain.cxx


.PHONY: all lib clean distclean

all: $(GENERATOR)

lib: $(GENERATOR_SO)

$(GENERATOR): $(GENERATOR_SO) $(MAIN_OBJ)
	        $(LD) $(LDFLAGS) $(MAIN_OBJ) $(GENERATOR_OBJ) $(GENERATOR_LIBS) -o $@
		        @echo "$@ done"

$(GENERATOR_SO): $(GENERATOR_OBJ)
	        $(LD) $(SOFLAGS) $(LDFLAGS) $^ -o $@ 


clean:
	        @rm -f $(GENERATOR_OBJ) $(MAIN_OBJ) 

distclean: clean
	        @rm -f $(GENERATOR) $(GENERATOR_SO) tags


TNtuple_Generator.o : TNtuple_generator.h


%.o: %.cxx
	        $(CXX)  $(CXXFLAGS) -c $< -o $@


-include Makefile_me


