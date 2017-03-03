# $Id: Makefile.users.in,v 1.6 2011/03/24 16:20:52 chulwoo Exp $
# Makefile.users.  Generated from Makefile.users.in by configure.
#   This is the makefile for all the commands
#   that are common for all testing applications.
#----------------------------------------------------------------------


SRCDIR = /home/jt2798/cps-jtu-build/public-1.2/build-qmp/../cps_pp-qmp
BUILDDIR = /home/jt2798/cps-jtu-build/public-1.2/build-qmp
QOS = 
QOSLIB = ${QOS}/quser/gcc-lib-user///
CC = mpicc
CXX = mpic++
AS  = as
LDFLAGS =  -fopenmp -L/home/jt2798/cps-jtu-build/public-1.2//local/lib -lqmp -L/home/jt2798/cps-jtu-build/public-1.2//local/lib -lqmp -lqio -llime 

me = $(notdir $(PWD))
BIN = NOARCH.x

#VPATH :=$(SRCDIR)/tests/$(me)



#
# include list for the Columbia code
#
INCLIST = -I${BUILDDIR} -I${SRCDIR}/include  -I/home/jt2798/cps-jtu-build/public-1.2//local/include   

CFLAGS= -g -O2 -fopenmp -O2
CXXFLAGS=  -fopenmp -O2
ASFLAGS= 
DFLAGS +=  -DUSE_OMP -DUSE_QMP -DUSE_QIO

#
# Libraries for the Columbia code
#
# (These are for the scalar version of the code)
#
#

.PHONY: cps clean


LIBLIST =\
  $(BUILDDIR)/cps.a \

#
#  targets
#


all:$(BIN)

.SUFFIXES:
.SUFFIXES:  .o .C .S .c

CSRC :=$(wildcard *.c)
CCSRC :=$(wildcard *.C)
SSRC :=$(wildcard *.S)

COBJ=$(CSRC:.c=.o)
CCOBJ=$(CCSRC:.C=.o)
SOBJ=$(SSRC:.S=.o)

OBJS_SRC = $(SOBJ) $(CCOBJ) $(COBJ)
OBJS := $(notdir $(OBJS_SRC))

$(BIN):  $(OBJS) $(LIBLIST)
	@echo OBJS = $(OBJS)
	$(CXX) $(OBJS) $(LIBLIST) $(LDFLAGS) -o $(BIN)

.c.o:
	$(CC) -o $@ $(CFLAGS) $(DFLAGS) -c $(INCLIST) $<
.C.o:
	$(CXX) -o $@ $(CXXFLAGS) $(DFLAGS) -c $(INCLIST) $<
.S.o:
	$(AS) -o $@ $(ASFLAGS) -c $(INCLIST) $<

cps:
	$(MAKE) -C $(BUILDDIR)

clean:
	rm -f *.dat *.o  $(BIN)
	rm -f ../regressions/*$(me).dat
	rm -f ../regressions/*$(me).checklog
