#include "graph.h"
#include "mwis_sewell/mwss.h"
#include <vector>

template <class SpecificMWSS>
class MWSS {

    public:
    MWSS(Graph& G) data(G) {};
    void solve(int k, vector<int>& pi, double goal, vector<int>& stable_set, double& weight) {
        data.solve(k,pi,goal,stable_set,weight);
        return;
    };

    private:
    SpecificMWSS data;

}

class Sewell {

    public:
    Sewell(Graph& G);
    void solve(int k, double goal, vector<int>& stable_set, vector<int>& pi, double& weight);

    private:
    vector<MWSSgraph> Mgraph;   // subgraphs of G (one per color) for the MWSS algorithm
    wstable_parameters Mparms;  // parameters for the MWSS algorithm

};
