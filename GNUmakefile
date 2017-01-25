#############################################################################
# Project : RooRarFit
# Package : RooRarFit
# File    : $Id: GNUmakefile.new,v 1.4 2012/07/04 09:01:35 fwilson Exp $
#
# Description:
#   Standalone makefile for the RooRarFit package.
#
# Instructions:
#   o Make sure the ROOTSYS environment variable is set
#
#   o Run 'make <target>'
#     - Default target makes shared library (libRooRarFit.so) for use in ROOT
#     - Add ROOTBUILD=debug to command line for debug version
#
# Build targets:
#   shlib    - make libRooRarFit.so (default target)
#   include  - make dependency files (*.d)
#   lib      - make libRooRarFit.a
#   bin      - make shlib and rarFit program
#   html     - make html documentation in html subdirectory
#   doc      - make UserGuide documentation in doc subdirectory
#   cleanbin - delete test binaries
#   clean    - delete all intermediate and final build objects
#
# Authors :
#
############################################################################# 
#
ifndef ROOTSYS
  out := $(shell echo "Please define ROOTSYS and retry" >&2)
endif

# Try to setup using default ROOT settings (sometimes it's in test and not etc)
exist := $(wildcard $(ROOTSYS)/test/Makefile.arch)
ifneq ($(strip $(exist)),)
  include $(ROOTSYS)/test/Makefile.arch
else
  include $(ROOTSYS)/etc/Makefile.arch
endif

MFLAGS   = -MM -Wno-deprecated 
WORKDIR  = $(PWD)/tmp
LIBDIR   = $(WORKDIR)
HTMLDOC  = html

ifeq ($(PLATFORM),macosx)
  MYFLAGS = -undefined dynamic_lookup
endif

# -------------------------------------------------------------

PACKAGE = RooRarFit
OBJDIR = $(WORKDIR)/objects
INCDIR = $(WORKDIR)/$(PACKAGE)

VPATH  = $(INCDIR) $(OBJDIR)

INCLUDES += -I.. -I$(INCDIR)/

# ensure pick up right libraries
ROOTCONFIG := $(ROOTSYS)/bin/root-config
ROOTLIBS   := $(shell $(ROOTCONFIG) --libs)

ROOTLIBS += -lMinuit -lHtml -lRooFit -lRooFitCore -lFoam
SKIPRDLIST = $(PACKAGE)_LinkDef.rdl
RDLLIST   = $(filter-out $(SKIPRDLIST),$(wildcard *.rdl))
SKIPLIST  = rarFit.cc runActions.cc
CINTFILE  = $(PACKAGE)Cint.cc
CINTOBJ   = $(OBJDIR)/$(PACKAGE)Cint.o
LIBFILE   = $(LIBDIR)/lib$(PACKAGE).a
SHLIBFILE = $(LIBDIR)/lib$(PACKAGE).$(DllSuf)
BINCCFILE = $(wildcard *.cc)
BINFILE   = $(LIBDIR)/rarFit
BINOFILE  = $(addprefix $(OBJDIR)/,$(patsubst %.cc,%.o,$(BINCCFILE)))

CXXFLAGS := $(shell $(ROOTCONFIG) --cflags)
CXXFLAGS += -Wno-parentheses -fPIC

ifeq ($(ROOTCINT),)
  ROOTCINT = $(ROOTSYS)/bin/rootcint
endif

default: shlib

# List of all includes files to copy from rdl
HHLIST = $(patsubst %.rdl,%.hh,$(RDLLIST))

# List of all source files to build (excluding main.cc)
CCLIST = $(filter-out $(SKIPLIST),$(wildcard *.cc))

# List of all object files to build (excluding main.o)
OLIST = $(addprefix $(OBJDIR)/,$(patsubst %.cc,%.o,$(CCLIST)))

# List of all dependency file to make
DLIST = $(patsubst %.rdl,%.d,$(RDLLIST))

# Implicit rule making all dependency Makefiles included at the end of this makefile
$(INCDIR)/%.d: %.cc $(HHLIST)
	@echo "Making $@"
	@set -e; $(CXX) $(MFLAGS) $(CXXFLAGS) $(INCLUDES) $< \
	          | sed 's/\($(notdir $*)\)\.o[ :]*/\1.o $(notdir $@) : /g' > $@; \
	        [ -s $@ ] || rm -f $@

# Implicit rule copying all RDL to INCDIR/HH
%.hh: %.rdl 
	@mkdir -p $(INCDIR)
	@cp $< $(INCDIR)/$@

$(OBJDIR)/%.o : %.cc 
	@echo "Compiling $<"
	@mkdir -p $(OBJDIR)
	@$(CXX) $(CXXFLAGS) -c $< -o $(OBJDIR)/$(notdir $@) $(INCLUDES)

# Rule to make ROOTCINT output file
$(CINTOBJ): $(RDLLIST) 
	@mkdir -p $(INCDIR)
	@mkdir -p $(OBJDIR)
	@echo "Generating dictionary from $(SKIPRDLIST)"
	@cd $(INCDIR) ; $(ROOTCINT) -f $(CINTFILE) -c $(INCLUDES) $(HHLIST)
	@echo "Compiling $(CINTFILE)"
	@$(CXX) $(CXXFLAGS) -c $(INCDIR)/$(CINTFILE) -o $(CINTOBJ) $(INCLUDES)

# Rule to combine objects into a library
$(LIBFILE): $(OLIST) $(CINTOBJ)
	@echo "Making $(LIBFILE)"
	@rm -f $(LIBFILE)
	@ar q $(LIBFILE) $(OLIST) $(CINTOBJ)
	@ranlib $(LIBFILE)

# Rule to combine objects into a shared library
$(SHLIBFILE): $(OLIST) $(CINTOBJ)
	@echo "Making shlib $(SHLIBFILE)"
	@rm -f $(SHLIBFILE)
	@$(CXX) $(CXXFLAGS) $(MYFLAGS) $(OLIST) $(CINTOBJ) $(SOFLAGS) -o $(SHLIBFILE)

# Rule to build binary
$(BINFILE): $(SHLIBFILE) $(BINOFILE)
	@echo "Making binary $(BINFILE)"
	@rm -f $(BINFILE)
	@$(CXX) $(CXXFLAGS) $(MYFLAGS) $(BINOFILE) $(CINTOBJ) $(ROOTLIBS) -o $(BINFILE)

# Useful build targets
lib     : $(LIBFILE) 
shlib   : $(SHLIBFILE)
bin     : shlib $(BINFILE)

cleanbin :
	@rm -f $(BINFILE)

clean : cleanbin
	@rm -f $(OBJDIR)/* $(INCDIR)/* $(LIBFILE) $(SHLIBFILE)

$(HTMLDOC)/index.html : $(SHLIBFILE)
	@echo "Making HTML documentation in $(HTMLDOC)"
	@( echo 'gSystem->Load("libRooUnfold");'; \
	   echo 'THtml h; h.SetOutputDir("$(HTMLDOC)");'; \
	   echo 'h.MakeAll();';\
	   echo '.q') \
	 | $(ROOTSYS)/bin/root -l -b

html : $(HTMLDOC)/index.html

.PHONY : shlib lib bin default clean cleanbin html doc

-include $(addprefix $(INCDIR)/,$(DLIST))
