
# Makefile

SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic
CC = g++

CCOPT = -O9 -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG -DIL_STD -Wall -std=c++11 -Wno-write-strings -g
CPLEXDIR      = /opt/ibm/ILOG/CPLEX_Studio127/cplex
CONCERTDIR    = /opt/ibm/ILOG/CPLEX_Studio127/concert
CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include
CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) -I.
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
LIBS = -L$(CPLEXLIBDIR) -L$(CONCERTLIBDIR)
CCLNFLAGS = -lm -lconcert -lilocplex -lcplex -pthread

all: vc vcroot st stroot bp bproot_dum bproot_psc bp_ccn genclassicinst genrandominst genmuinst gengraph checker


# Tools made by Daniel

st: listcola.cpp
	$(CC) -o $@ $^ -DSTABLEMODEL -DTUNEDPARAMS $(CCFLAGS) $(LIBS) $(CCLNFLAGS)

stroot: listcola.cpp
	$(CC) -o $@ $^ -DSTABLEMODEL -DTUNEDPARAMS -DONLYRELAXATION $(CCFLAGS) $(LIBS) $(CCLNFLAGS)

vc: listcola.cpp
	$(CC) -o $@ $^ -DSYMMETRYRESTR1 -DSYMMETRYRESTR3 $(CCFLAGS) $(LIBS) $(CCLNFLAGS)

vcroot: listcola.cpp
	$(CC) -o $@ $^ -DSYMMETRYRESTR1 -DSYMMETRYRESTR3 -DONLYRELAXATION $(CCFLAGS) $(LIBS) $(CCLNFLAGS)

genclassicinst: genclassicinst.cpp
	$(CC) -o $@ $^ $(CCFLAGS) $(LIBS) $(CCLNFLAGS)

genrandominst: genrandominst.cpp
	$(CC) -o $@ $^ $(CCFLAGS)

genmuinst: genmuinst.cpp
	$(CC) -o $@ $^ $(CCFLAGS) $(LIBS) $(CCLNFLAGS)

gengraph: gengraph.cpp
	$(CC) -o $@ $^ $(CCFLAGS) -lm

checker: checker.cpp
	$(CC) -o $@ $^ $(CCFLAGS) -lm


# Code made by Mauro

bp: main.o bp.o lp.o stable.o graph.o io.o mwis_sewell/wstable.o
	$(CC) -o $@ $^ $(CCFLAGS) $(LIBS) $(CCLNFLAGS)

bp_ccn: main.o bp.o lp_ccn.o stable.o graph.o io.o mwis_sewell/wstable.o
	$(CC) -o $@ $^ $(CCFLAGS) $(LIBS) $(CCLNFLAGS)

bproot_dum: main.o bproot.o lproot.o stable.o graph.o io.o mwis_sewell/wstable.o
	$(CC) -o $@ $^ $(CCFLAGS) $(LIBS) $(CCLNFLAGS)

bproot_psc: main.o bproot.o lprootheur.o stable.o graph.o io.o mwis_sewell/wstable.o
	$(CC) -o $@ $^ $(CCFLAGS) $(LIBS) $(CCLNFLAGS)

main.o: main.cpp bp.h lp.h stable.h graph.h io.h
	$(CC) -c -o $@ $< $(CCFLAGS)

bp.o: bp.cpp bp.h lp.h
	$(CC) -c -o $@ $< $(CCFLAGS)

bproot.o: bp.cpp bp.h lp.h
	$(CC) -c -o $@ $< $(CCFLAGS) -DONLY_RELAXATION

lp.o: lp.cpp lp.h
	$(CC) -c -o $@ $< $(CCFLAGS)

lpheur.o: lp.cpp lp.h
	$(CC) -c -o $@ $< $(CCFLAGS) -DINITIAL_COLUMNS_HEURISTIC

lp_ccn.o: lp.cpp lp.h
	$(CC) -c -o $@ $< $(CCFLAGS) -DINITIAL_COLUMNS_FATHER

lproot.o: lp.cpp lp.h
	$(CC) -c -o $@ $< $(CCFLAGS) -DONLY_RELAXATION

lprootheur.o: lp.cpp lp.h
	$(CC) -c -o $@ $< $(CCFLAGS) -DINITIAL_COLUMNS_HEURISTIC -DONLY_RELAXATION

stable.o: stable.cpp stable.h
	$(CC) -c -o $@ $< $(CCFLAGS)

graph.o: graph.cpp graph.h
	$(CC) -c -o $@ $< $(CCFLAGS)

io.o: io.cpp io.h
	$(CC) -c -o $@ $< $(CCFLAGS)

sewell:
	$(MAKE) -C mwis_sewell


.PHONY: clean

clean:
	$(MAKE) -C mwis_sewell clean
	rm -f *.o
	rm -f vc vcroot st stroot bp bproot_dum bproot_psc bp_ccn genclassicinst genrandominst genmuinst gengraph checker
