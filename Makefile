
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

all: vc vcroot st stroot bproot_dum bproot_psc  \
     genclassicinst gensetcoverinst gensetcoverrailinst  \
     genrandominst genmuinst gengraph checker \
     bp_edg0 bp_edgalt0 bp_edg1 bp_edgalt1 bp_edg2 bp_edgalt2  \
     bp_clr0 bp_clralt0 bp_clr1 bp_clralt1 bp_clr2 bp_clralt2  \
     bp_clr2N3 bp_clralt2N3  \
     bp_clr2N4 bp_clralt2N4  \
     bp_clr2N5 bp_clralt2N5


# Example

mybp0: main.o bp.o mylp.o mygraph.o io.o mwis_sewell/wstable.o
	$(CC) -o $@ $^ $(CCFLAGS) $(LIBS) $(CCLNFLAGS)

mylp.o: lp.cpp lp.h graph.h
	$(CC) -c -o $@ $< $(CCFLAGS) -DINITIAL_COLUMN_STRATEGY=0 -DBRANCHING_STRATEGY=1 -DPREPROCESSING=0 -DVARIABLE_SELECTION=1 -DN_BRANCHS=2

mygraph.o: graph.cpp graph.h
	$(CC) -c -o $@ $< $(CCFLAGS) -DBRANCHING_STRATEGY=1 -DPREPROCESSING=0


# Code made by Mauro

bp_edg0: main.o bp.o lp_edg0.o graph_edg0.o io.o mwis_sewell/wstable.o
	$(CC) -o $@ $^ $(CCFLAGS) $(LIBS) $(CCLNFLAGS)

lp_edg0.o: lp.cpp lp.h graph.h
	$(CC) -c -o $@ $< $(CCFLAGS) -DINITIAL_COLUMN_STRATEGY=0 -DBRANCHING_STRATEGY=0 -DPREPROCESSING=0 -DVARIABLE_SELECTION=0

graph_edg0.o: graph.cpp graph.h
	$(CC) -c -o $@ $< $(CCFLAGS) -DBRANCHING_STRATEGY=0 -DPREPROCESSING=0

bp_edgalt0: main.o bp.o lp_edgalt0.o graph_edg0.o io.o mwis_sewell/wstable.o
	$(CC) -o $@ $^ $(CCFLAGS) $(LIBS) $(CCLNFLAGS)

lp_edgalt0.o: lp.cpp lp.h graph.h
	$(CC) -c -o $@ $< $(CCFLAGS) -DINITIAL_COLUMN_STRATEGY=0 -DBRANCHING_STRATEGY=0 -DPREPROCESSING=0 -DVARIABLE_SELECTION=1


bp_edg1: main.o bp.o lp_edg1.o graph_edg1.o io.o mwis_sewell/wstable.o
	$(CC) -o $@ $^ $(CCFLAGS) $(LIBS) $(CCLNFLAGS)

lp_edg1.o: lp.cpp lp.h graph.h
	$(CC) -c -o $@ $< $(CCFLAGS) -DINITIAL_COLUMN_STRATEGY=0 -DBRANCHING_STRATEGY=0 -DPREPROCESSING=1 -DVARIABLE_SELECTION=0

graph_edg1.o: graph.cpp graph.h
	$(CC) -c -o $@ $< $(CCFLAGS) -DBRANCHING_STRATEGY=0 -DPREPROCESSING=1

bp_edgalt1: main.o bp.o lp_edgalt1.o graph_edg1.o io.o mwis_sewell/wstable.o
	$(CC) -o $@ $^ $(CCFLAGS) $(LIBS) $(CCLNFLAGS)

lp_edgalt1.o: lp.cpp lp.h graph.h
	$(CC) -c -o $@ $< $(CCFLAGS) -DINITIAL_COLUMN_STRATEGY=0 -DBRANCHING_STRATEGY=0 -DPREPROCESSING=1 -DVARIABLE_SELECTION=1


bp_edg2: main.o bp.o lp_edg2.o graph_edg2.o io.o mwis_sewell/wstable.o
	$(CC) -o $@ $^ $(CCFLAGS) $(LIBS) $(CCLNFLAGS)

lp_edg2.o: lp.cpp lp.h graph.h
	$(CC) -c -o $@ $< $(CCFLAGS) -DINITIAL_COLUMN_STRATEGY=0 -DBRANCHING_STRATEGY=0 -DPREPROCESSING=2 -DVARIABLE_SELECTION=0

graph_edg2.o: graph.cpp graph.h
	$(CC) -c -o $@ $< $(CCFLAGS) -DBRANCHING_STRATEGY=0 -DPREPROCESSING=2

bp_edgalt2: main.o bp.o lp_edgalt2.o graph_edg2.o io.o mwis_sewell/wstable.o
	$(CC) -o $@ $^ $(CCFLAGS) $(LIBS) $(CCLNFLAGS)

lp_edgalt2.o: lp.cpp lp.h graph.h
	$(CC) -c -o $@ $< $(CCFLAGS) -DINITIAL_COLUMN_STRATEGY=0 -DBRANCHING_STRATEGY=0 -DPREPROCESSING=2 -DVARIABLE_SELECTION=1


bp_clr0: main.o bp.o lp_clr0.o graph_clr0.o io.o mwis_sewell/wstable.o
	$(CC) -o $@ $^ $(CCFLAGS) $(LIBS) $(CCLNFLAGS)

lp_clr0.o: lp.cpp lp.h graph.h
	$(CC) -c -o $@ $< $(CCFLAGS) -DINITIAL_COLUMN_STRATEGY=0 -DBRANCHING_STRATEGY=1 -DPREPROCESSING=0 -DVARIABLE_SELECTION=0 -DN_BRANCHS=2

graph_clr0.o: graph.cpp graph.h
	$(CC) -c -o $@ $< $(CCFLAGS) -DBRANCHING_STRATEGY=1 -DPREPROCESSING=0

bp_clralt0: main.o bp.o lp_clralt0.o graph_clr0.o io.o mwis_sewell/wstable.o
	$(CC) -o $@ $^ $(CCFLAGS) $(LIBS) $(CCLNFLAGS)

lp_clralt0.o: lp.cpp lp.h graph.h
	$(CC) -c -o $@ $< $(CCFLAGS) -DINITIAL_COLUMN_STRATEGY=0 -DBRANCHING_STRATEGY=1 -DPREPROCESSING=0 -DVARIABLE_SELECTION=1 -DN_BRANCHS=2


bp_clr1: main.o bp.o lp_clr1.o graph_clr1.o io.o mwis_sewell/wstable.o
	$(CC) -o $@ $^ $(CCFLAGS) $(LIBS) $(CCLNFLAGS)

lp_clr1.o: lp.cpp lp.h graph.h
	$(CC) -c -o $@ $< $(CCFLAGS) -DINITIAL_COLUMN_STRATEGY=0 -DBRANCHING_STRATEGY=1 -DPREPROCESSING=1 -DVARIABLE_SELECTION=0 -DN_BRANCHS=2

graph_clr1.o: graph.cpp graph.h
	$(CC) -c -o $@ $< $(CCFLAGS) -DBRANCHING_STRATEGY=1 -DPREPROCESSING=1

bp_clralt1: main.o bp.o lp_clralt1.o graph_clr1.o io.o mwis_sewell/wstable.o
	$(CC) -o $@ $^ $(CCFLAGS) $(LIBS) $(CCLNFLAGS)

lp_clralt1.o: lp.cpp lp.h graph.h
	$(CC) -c -o $@ $< $(CCFLAGS) -DINITIAL_COLUMN_STRATEGY=0 -DBRANCHING_STRATEGY=1 -DPREPROCESSING=1 -DVARIABLE_SELECTION=1 -DN_BRANCHS=2


bp_clr2: main.o bp.o lp_clr2.o graph_clr2.o io.o mwis_sewell/wstable.o
	$(CC) -o $@ $^ $(CCFLAGS) $(LIBS) $(CCLNFLAGS)

lp_clr2.o: lp.cpp lp.h graph.h
	$(CC) -c -o $@ $< $(CCFLAGS) -DINITIAL_COLUMN_STRATEGY=0 -DBRANCHING_STRATEGY=1 -DPREPROCESSING=2 -DVARIABLE_SELECTION=0 -DN_BRANCHS=2

graph_clr2.o: graph.cpp graph.h
	$(CC) -c -o $@ $< $(CCFLAGS) -DBRANCHING_STRATEGY=1 -DPREPROCESSING=2

bp_clralt2: main.o bp.o lp_clralt2.o graph_clr2.o io.o mwis_sewell/wstable.o
	$(CC) -o $@ $^ $(CCFLAGS) $(LIBS) $(CCLNFLAGS)

lp_clralt2.o: lp.cpp lp.h graph.h
	$(CC) -c -o $@ $< $(CCFLAGS) -DINITIAL_COLUMN_STRATEGY=0 -DBRANCHING_STRATEGY=1 -DPREPROCESSING=2 -DVARIABLE_SELECTION=1 -DN_BRANCHS=2


bp_clr2N3: main.o bp.o lp_clr2N3.o graph_clr2.o io.o mwis_sewell/wstable.o
	$(CC) -o $@ $^ $(CCFLAGS) $(LIBS) $(CCLNFLAGS)

lp_clr2N3.o: lp.cpp lp.h graph.h
	$(CC) -c -o $@ $< $(CCFLAGS) -DINITIAL_COLUMN_STRATEGY=0 -DBRANCHING_STRATEGY=1 -DPREPROCESSING=2 -DVARIABLE_SELECTION=0 -DN_BRANCHS=3

bp_clralt2N3: main.o bp.o lp_clralt2N3.o graph_clr2.o io.o mwis_sewell/wstable.o
	$(CC) -o $@ $^ $(CCFLAGS) $(LIBS) $(CCLNFLAGS)

lp_clralt2N3.o: lp.cpp lp.h graph.h
	$(CC) -c -o $@ $< $(CCFLAGS) -DINITIAL_COLUMN_STRATEGY=0 -DBRANCHING_STRATEGY=1 -DPREPROCESSING=2 -DVARIABLE_SELECTION=1 -DN_BRANCHS=3


bp_clr2N4: main.o bp.o lp_clr2N4.o graph_clr2.o io.o mwis_sewell/wstable.o
	$(CC) -o $@ $^ $(CCFLAGS) $(LIBS) $(CCLNFLAGS)

lp_clr2N4.o: lp.cpp lp.h graph.h
	$(CC) -c -o $@ $< $(CCFLAGS) -DINITIAL_COLUMN_STRATEGY=0 -DBRANCHING_STRATEGY=1 -DPREPROCESSING=2 -DVARIABLE_SELECTION=0 -DN_BRANCHS=4

bp_clralt2N4: main.o bp.o lp_clralt2N4.o graph_clr2.o io.o mwis_sewell/wstable.o
	$(CC) -o $@ $^ $(CCFLAGS) $(LIBS) $(CCLNFLAGS)

lp_clralt2N4.o: lp.cpp lp.h graph.h
	$(CC) -c -o $@ $< $(CCFLAGS) -DINITIAL_COLUMN_STRATEGY=0 -DBRANCHING_STRATEGY=1 -DPREPROCESSING=2 -DVARIABLE_SELECTION=1 -DN_BRANCHS=4


bp_clr2N5: main.o bp.o lp_clr2N5.o graph_clr2.o io.o mwis_sewell/wstable.o
	$(CC) -o $@ $^ $(CCFLAGS) $(LIBS) $(CCLNFLAGS)

lp_clr2N5.o: lp.cpp lp.h graph.h
	$(CC) -c -o $@ $< $(CCFLAGS) -DINITIAL_COLUMN_STRATEGY=0 -DBRANCHING_STRATEGY=1 -DPREPROCESSING=2 -DVARIABLE_SELECTION=0 -DN_BRANCHS=5

bp_clralt2N5: main.o bp.o lp_clralt2N5.o graph_clr2.o io.o mwis_sewell/wstable.o
	$(CC) -o $@ $^ $(CCFLAGS) $(LIBS) $(CCLNFLAGS)

lp_clralt2N5.o: lp.cpp lp.h graph.h
	$(CC) -c -o $@ $< $(CCFLAGS) -DINITIAL_COLUMN_STRATEGY=0 -DBRANCHING_STRATEGY=1 -DPREPROCESSING=2 -DVARIABLE_SELECTION=1 -DN_BRANCHS=5


bp: main.o bp.o lp.o graph.o io.o mwis_sewell/wstable.o
	$(CC) -o $@ $^ $(CCFLAGS) $(LIBS) $(CCLNFLAGS)

bp_dum: main.o bp.o lp_dum.o graph.o io.o mwis_sewell/wstable.o
	$(CC) -o $@ $^ $(CCFLAGS) $(LIBS) $(CCLNFLAGS)

bp_psc: main.o bp.o lp_psc.o graph.o io.o mwis_sewell/wstable.o
	$(CC) -o $@ $^ $(CCFLAGS) $(LIBS) $(CCLNFLAGS)

bp_ccn: main.o bp.o lp_ccn.o graph.o io.o mwis_sewell/wstable.o
	$(CC) -o $@ $^ $(CCFLAGS) $(LIBS) $(CCLNFLAGS)

bp_poo: main.o bp.o lp_poo.o graph_poo.o io.o mwis_sewell/wstable.o
	$(CC) -o $@ $^ $(CCFLAGS) $(LIBS) $(CCLNFLAGS)

bp_clr: main.o bp.o lp_clr.o graph_clr.o io.o mwis_sewell/wstable.o
	$(CC) -o $@ $^ $(CCFLAGS) $(LIBS) $(CCLNFLAGS)

bp_ind: main.o bp.o lp_ind.o graph_ind.o io.o mwis_sewell/wstable.o
	$(CC) -o $@ $^ $(CCFLAGS) $(LIBS) $(CCLNFLAGS)

bproot_dum: main.o bproot.o lproot.o graph.o io.o mwis_sewell/wstable.o
	$(CC) -o $@ $^ $(CCFLAGS) $(LIBS) $(CCLNFLAGS)

bproot_psc: main.o bproot.o lprootheur.o graph.o io.o mwis_sewell/wstable.o
	$(CC) -o $@ $^ $(CCFLAGS) $(LIBS) $(CCLNFLAGS)

main.o: main.cpp bp.h lp.h graph.h io.h
	$(CC) -c -o $@ $< $(CCFLAGS)

bp.o: bp.cpp bp.h lp.h
	$(CC) -c -o $@ $< $(CCFLAGS)

bproot.o: bp.cpp bp.h lp.h
	$(CC) -c -o $@ $< $(CCFLAGS) -DONLY_RELAXATION

lp.o: lp.cpp lp.h graph.h
	$(CC) -c -o $@ $< $(CCFLAGS) -DINITIAL_COLUMN_STRATEGY=0 -DBRANCHING_STRATEGY=0

lp_dum.o: lp.cpp lp.h graph.h
	$(CC) -c -o $@ $< $(CCFLAGS) -DINITIAL_COLUMN_STRATEGY=0 -DBRANCHING_STRATEGY=0

lp_psc.o: lp.cpp lp.h graph.h
	$(CC) -c -o $@ $< $(CCFLAGS) -DINITIAL_COLUMN_STRATEGY=2 -DBRANCHING_STRATEGY=0

lp_ccn.o: lp.cpp lp.h graph.h
	$(CC) -c -o $@ $< $(CCFLAGS) -DINITIAL_COLUMN_STRATEGY=1 -DBRANCHING_STRATEGY=0

lp_poo.o: lp.cpp lp.h graph.h
	$(CC) -c -o $@ $< $(CCFLAGS) -DINITIAL_COLUMN_STRATEGY=1 -DBRANCHING_STRATEGY=0 -DSTABLE_POOL

lp_clr.o: lp.cpp lp.h graph.h
	$(CC) -c -o $@ $< $(CCFLAGS) -DINITIAL_COLUMN_STRATEGY=0 -DBRANCHING_STRATEGY=1

lp_ind.o: lp.cpp lp.h graph.h
	$(CC) -c -o $@ $< $(CCFLAGS) -DINITIAL_COLUMN_STRATEGY=0 -DBRANCHING_STRATEGY=2

lproot.o: lp.cpp lp.h graph.h
	$(CC) -c -o $@ $< $(CCFLAGS) -DONLY_RELAXATION -DINITIAL_COLUMN_STRATEGY=0

lprootheur.o: lp.cpp lp.h graph.h
	$(CC) -c -o $@ $< $(CCFLAGS) -DONLY_RELAXATION -DINITIAL_COLUMN_STRATEGY=2

graph.o: graph.cpp graph.h
	$(CC) -c -o $@ $< $(CCFLAGS) -DBRANCHING_STRATEGY=0

graph_poo.o: graph.cpp graph.h
	$(CC) -c -o $@ $< $(CCFLAGS) -DBRANCHING_STRATEGY=0 -DSTABLE_POOL

graph_clr.o: graph.cpp graph.h
	$(CC) -c -o $@ $< $(CCFLAGS) -DBRANCHING_STRATEGY=1

graph_ind.o: graph.cpp graph.h
	$(CC) -c -o $@ $< $(CCFLAGS) -DBRANCHING_STRATEGY=2

io.o: io.cpp io.h
	$(CC) -c -o $@ $< $(CCFLAGS)

sewell:
	$(MAKE) -C mwis_sewell


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

gensetcoverinst: gensetcoverinst.cpp
	$(CC) -o $@ $^ $(CCFLAGS)

gensetcoverrailinst: gensetcoverinst.cpp
	$(CC) -o $@ $^ -DRAIL_INSTANCE $(CCFLAGS)

genrandominst: genrandominst.cpp
	$(CC) -o $@ $^ $(CCFLAGS)

genmuinst: genmuinst.cpp
	$(CC) -o $@ $^ $(CCFLAGS) $(LIBS) $(CCLNFLAGS)

gengraph: gengraph.cpp
	$(CC) -o $@ $^ $(CCFLAGS) -lm

checker: checker.cpp
	$(CC) -o $@ $^ $(CCFLAGS) -lm


.PHONY: clean

clean:
	$(MAKE) -C mwis_sewell clean
	rm -f *.o
	rm -f vc vcroot st stroot bproot_dum bproot_psc bp_edg0 bp_edgalt0 bp_edg1 bp_edgalt1 bp_edg2 bp_edgalt2 bp_clr0 bp_clralt0 bp_clr1 bp_clralt1 bp_clr2 bp_clralt2 bp_clr0N3 bp_clralt0N3 bp_clr1N3 bp_clralt1N3 bp_clr2N3 bp_clralt2N3 bp_clr0N4 bp_clralt0N4 bp_clr1N4 bp_clralt1N4 bp_clr2N4 bp_clralt2N4 bp_clr0N5 bp_clralt0N5 bp_clr1N5 bp_clralt1N5 bp_clr2N5 bp_clralt2N5 genclassicinst gensetcoverinst gensetcoverrailinst genrandominst genmuinst gengraph checker

