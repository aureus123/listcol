
# Makefile

SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic
CC = g++

CCOPT = -O9 -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG -DIL_STD -Wall -std=c++11 -Wno-write-strings
CPLEXDIR      = /opt/ibm/ILOG/CPLEX_Studio127/cplex
CONCERTDIR    = /opt/ibm/ILOG/CPLEX_Studio127/concert
CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include
CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) -I.
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
LIBS = -L$(CPLEXLIBDIR) -L$(CONCERTLIBDIR)
CCLNFLAGS = -lm -lconcert -lilocplex -lcplex -pthread

#all: listcola1 listcola1t listcola2_edge listcola2 listcola2t listcola2t1 listcola2t12 listcola2t13 listcola2t13root genclassicinst genrandominst genmuinst gengraph checker listcolabproot listcolabprootheur listcolabp listcolabpcopy

all: listcola1t listcola2t13 listcola2t13root genclassicinst genrandominst genmuinst gengraph checker listcolabproot listcolabprootheur listcolabp listcolabpcopy

listcola1: listcola.cpp
	$(CC) -o $@ $^ -DSTABLEMODEL $(CCFLAGS) $(LIBS) $(CCLNFLAGS)

listcola1t: listcola.cpp
	$(CC) -o $@ $^ -DSTABLEMODEL -DTUNEDPARAMS $(CCFLAGS) $(LIBS) $(CCLNFLAGS)

listcola2_edge: listcola.cpp
	$(CC) -o $@ $^ -DEDGEINEQ $(CCFLAGS) $(LIBS) $(CCLNFLAGS)

listcola2: listcola.cpp
	$(CC) -o $@ $^ $(CCFLAGS) $(LIBS) $(CCLNFLAGS)

listcola2t: listcola.cpp
	$(CC) -o $@ $^ -DTUNEDPARAMS $(CCFLAGS) $(LIBS) $(CCLNFLAGS)

listcola2t1: listcola.cpp
	$(CC) -o $@ $^ -DTUNEDPARAMS -DSYMMETRYRESTR1 $(CCFLAGS) $(LIBS) $(CCLNFLAGS)

listcola2t12: listcola.cpp
	$(CC) -o $@ $^ -DTUNEDPARAMS -DSYMMETRYRESTR1 -DSYMMETRYRESTR2 $(CCFLAGS) $(LIBS) $(CCLNFLAGS)

listcola2t13: listcola.cpp
	$(CC) -o $@ $^ -DTUNEDPARAMS -DSYMMETRYRESTR1 -DSYMMETRYRESTR3 $(CCFLAGS) $(LIBS) $(CCLNFLAGS)

listcola2t13root: listcola.cpp
	$(CC) -o $@ $^ -DTUNEDPARAMS -DSYMMETRYRESTR1 -DSYMMETRYRESTR3 -DONLYRELAXATION $(CCFLAGS) $(LIBS) $(CCLNFLAGS)


genclassicinst: genclassicinst.cpp
	$(CC) -o $@ $^ $(CCFLAGS) $(LIBS) $(CCLNFLAGS)

genrandominst: genrandominst.cpp
	$(CC) -o $@ $^ $(CCFLAGS) -lm

genmuinst: genmuinst.cpp
	$(CC) -o $@ $^ $(CCFLAGS) -lm

gengraph: gengraph.cpp
	$(CC) -o $@ $^ $(CCFLAGS) -lm

checker: checker.cpp
	$(CC) -o $@ $^ $(CCFLAGS) -lm



listcolabp: main.o bp.o lp.o stable.o graph.o io.o mwis_sewell/wstable.o
	$(CC) -o $@ $^ $(CCFLAGS) $(LIBS) $(CCLNFLAGS)

listcolabpcopy: main.o bp.o lpcopy.o stable.o graph.o io.o mwis_sewell/wstable.o
	$(CC) -o $@ $^ $(CCFLAGS) $(LIBS) $(CCLNFLAGS)

listcolabproot: main.o bproot.o lp.o stable.o graph.o io.o mwis_sewell/wstable.o
	$(CC) -o $@ $^ $(CCFLAGS) $(LIBS) $(CCLNFLAGS)

listcolabprootheur: main.o bproot.o lpheur.o stable.o graph.o io.o mwis_sewell/wstable.o
	$(CC) -o $@ $^ $(CCFLAGS) $(LIBS) $(CCLNFLAGS)

main.o: main.cpp bp.h lp.h stable.h graph.h io.h
	$(CC) -c -o $@ $< $(CCFLAGS)

bp.o: bp.cpp bp.h
	$(CC) -c -o $@ $< $(CCFLAGS)

bproot.o: bp.cpp bp.h
	$(CC) -c -o $@ $< $(CCFLAGS) -DONLY_RELAXATION

lp.o: lp.cpp lp.h
	$(CC) -c -o $@ $< $(CCFLAGS)

lpheur.o: lp.cpp lp.h
	$(CC) -c -o $@ $< $(CCFLAGS) -DINITIAL_COLUMNS_HEURISTIC

lpcopy.o: lp.cpp lp.h
	$(CC) -c -o $@ $< $(CCFLAGS) -DINITIAL_COLUMNS_FATHER

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
	rm -f listcola1 listcola1t listcola2_edge listcola2 listcola2t listcola2t1 listcola2t12 listcola2t13 listcola2t13root genclassicinst genrandominst genmuinst gengraph checker listcolabproot listcolabprootheur listcolabp listcolabpcopy

