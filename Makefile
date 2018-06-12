
# Makefile

SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic
CC = g++

CCOPT = -O9 -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG -DIL_STD -Wall
CPLEXDIR      = /opt/ibm/ILOG/CPLEX_Studio127/cplex
CONCERTDIR    = /opt/ibm/ILOG/CPLEX_Studio127/concert
CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include
CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) -I.
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
LIBS = -L$(CPLEXLIBDIR) -L$(CONCERTLIBDIR)
CCLNFLAGS = -lm -lconcert -lilocplex -lcplex -pthread -std=c++11

main: bp.cpp io.o graph.o lp.o io.o stable.o mwis_sewell/wstable.o
	$(CC) -o $@ $^ $(CCFLAGS) $(LIBS) $(CCLNFLAGS)

lp.o: lp.cpp lp.h graph.o stable.o
	$(CC) -c -o $@ $< $(CCLNFLAGS) $(CCFLAGS)

stable.o: stable.cpp stable.h mwis_sewell/wstable.o
	$(CC) -c -o $@ $< $(CCLNFLAGS)

sewell:
	$(MAKE) -C mwis_sewell

io.o: io.cpp io.h
	$(CC) -c -o $@ $< $(CCLNFLAGS)

graph.o: graph.cpp graph.h
	$(CC) -c -o $@ $< $(CCLNFLAGS)

.PHONY: clean

clean:
	$(MAKE) -C mwis_sewell clean
	rm -f *.o
	rm -f main
