/**
    This file is part of exactcolors.

    exactcolors is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    exactcolors is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with exactcolors.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "mwss.h"
#include "mwss_ext.h"

static const int MAX_NODES=10000;

static int check_ncount(int ncount)
{
   int rval = 0;
   if (ncount > MAX_NODES) {
      fprintf(stderr,"SEWELL_init_graph ncount exceeds MAX_NODES.");
      rval = 1;
   }
   return rval;
}


static int SEWELL_init_graph(MWSSgraphpnt graph,
                             int ncount, int ecount,
                             const int elist[])
{
   int rval = 0;
   int      i, row, col;

   // Initialize the node names and degrees.
   graph->n_nodes = ncount;

   rval = allocate_graph(graph, ncount);
   MWIScheck_rval(rval,"Failed in allocate_graph");

   //MALLOC(node_list, n_nodes+1, tnode);
   for(i = 0; i <= graph->n_nodes; i++) {
      graph->node_list[i].name = i;
      graph->node_list[i].degree = 0;
   }

   // Initialize the adjacency matrix.
   for(row = 1; row <= graph->n_nodes; row++) {
      graph->adj[row][row] = 0;
      for(col = row + 1; col <= graph->n_nodes; col++) {
         graph->adj[row][col] = 0;
         graph->adj[col][row] = 0;
      }
   }

   for (i = 0; i < ecount; ++i) {
      row = elist[2*i]     + 1;
      col = elist[2*i + 1] + 1;
      graph->adj[row][col] = 1;
      graph->adj[col][row] = 1;
   }

   rval = build_graph(graph);
   MWIScheck_rval(rval,"Failed in build_graph");


 CLEANUP:
   return rval;
}

extern
int SEWELL_optimize(int ** newset,
                    int   *nnewset,
                    int ncount, int ecount, const int elist[],
                    NWT nweights[],
                    NWT lower_bound,
                    NWT goal)
{
   double cpu_limit = -1.0;
   return SEWELL_heur (newset,nnewset,ncount,ecount,elist,
                       nweights,lower_bound,goal,cpu_limit);
}
extern
int SEWELL_heur(int ** newset,
                                  int   *nnewset,
                                  int ncount, int ecount, const int elist[],
                                  NWT nweights[],
                                  NWT lower_bound,
                                  NWT goal,
                                  double cpu_limit)
{
   int rval = 0;
   wstable_info   info;
   int i;
   MWSSgraph      graph;
   MWSSdata       data;
   wstable_parameters parms;
   double density =  ((double) ecount) / ((double) (ncount * ( ncount - 1))) * 2.0;
   reset_pointers(&graph, &data, &info);

   if(check_ncount(ncount)) {goto CLEANUP;}

   default_parameters(&parms);
   parms.cpu_limit = cpu_limit;

   if (density < 0.2) {
      parms.clique_cover = 2;
   }

   SEWELL_init_graph(&graph,ncount,ecount,elist);
   for (i = 0; i < ncount; ++i) {
      graph.weight[i+1] = nweights[i];
   }


   rval = initialize_max_wstable(&graph,&info);
   MWIScheck_rval(rval, "Failed in initialize_max_wstable");

   rval = call_max_wstable(&graph,&data,&parms,&info, goal,lower_bound);
   MWIScheck_rval(rval, "Failed in call_max_wstable");

   if (*newset) {free(*newset);}

   *nnewset = data.n_best;
   *newset  = (int*) malloc (data.n_best * sizeof(int));
   if (!*newset) {
      fprintf(stderr,"Failed to allocate *newset");
      rval = 1; goto CLEANUP;
   }

   for (i = 1; i <= data.n_best; ++i) {
      (*newset)[i - 1] = data.best_sol[i]->name - 1;
   }

   if (parms.cpu_limit >= 0 && info.cpu > parms.cpu_limit) {
      printf("cpu_limit of %f seconds exceeded: %f seconds. Solution may not optimum.\n", parms.cpu_limit,info.cpu);
      rval = SEWELL_TIMEOUT;
      goto CLEANUP;
   }



 CLEANUP:
   free_max_wstable(&graph,&data, &info);
   return rval;
}

extern
int SEWELL_node_limit()
{
   return MAX_NODES;
}
