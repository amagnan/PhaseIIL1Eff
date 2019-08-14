#Necessary to use shell built-in commands
SHELL=bash

FASTJET_INSTALL=/home/amagnan/SOFTWARE/usr/fastjet-install

USERINCLUDES += -I$(ROOTSYS)/include/
USERINCLUDES += -I$(FASTJET_INSTALL)/include/
#USERINCLUDES += -isystem $(BOOSTSYS)/include

# Define libraries to link
USERLIBS += $(shell root-config --glibs) #-lRooFit -lGenVector #-lTreePlayer -lTMVA
USERLIBS += -Wl,-rpath,$(FASTJET_INSTALL)/lib -lm  -L$(FASTJET_INSTALL)/lib -lfastjettools -lfastjet -lfastjetplugins -lsiscone_spherical -lsiscone
#USERLIBS += -L$(BOOSTSYS)/lib -lboost_regex -lboost_program_options -lboost_filesystem

CXXFLAGS = -Wall -W -std=c++11 
#CXXFLAGS = -Wall -W -O2 -std=c++0x 
LDFLAGS = -shared -Wall -W

CXX=g++

CXXFLAGS += $(USERINCLUDES)
LIBS += $(USERLIBS)

# A list of directories
BASEDIR = $(shell pwd)
EXEDIR = $(BASEDIR)
OBJDIR = $(BASEDIR)
TESTDIR = $(BASEDIR)
OBJ_EXT=o
TEST_EXT=C

# Build a list of srcs and bins to build
EXES=$(wildcard $(BASEDIR)/plotEff.C)
BINS=$(BASEDIR)/plotEff


all: $(BINS)

$(EXEDIR)/plotEff:  $(TESTDIR)/plotEff.C
	$(CXX) -o $@ $(CXXFLAGS) $< $(LIBS)


