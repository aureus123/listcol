#include <vector>
#include <set>
#include <iostream>
#include <climits>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <stack>

using namespace std;

// bye - finish executing and show a message
void bye(string str) {
	cout << str << endl;
	exit(1);
}

class Graph {

    public:

    int vertices;
    int edges;
    int colors;

    vector<vector<int> > adj; // adjacency list
    vector<set<int> > L; // color list

    Graph (char *graph_filename, char *list_filename) {

	    // Open file 
	    FILE *stream = fopen(graph_filename, "rt");
	    if (!stream) bye("Graph file cannot be opened");
	    fscanf(stream, "%d:%d\n", &vertices, &edges);

	    // Do not accept graph of less than 4 vertices or stable sets
	    if (vertices < 4) bye("Number of vertices out range!");
	    if (edges < 1 || edges > vertices*(vertices - 1) / 2) bye("Number of edges out of range!");

        adj.resize(vertices);

	    // Read edges
	    for (int e = 0; e < edges; e++) {
		    int u, v;
		    fscanf(stream, "%d,%d\n", &u, &v);
		    if (u < 0 || u >= v || v >= vertices) {
			    cout << "Error reading edge " << e + 1 << "!" << endl;
			    bye("Bye!");
		    }
            adj[u].push_back(v);
            adj[v].push_back(u);
	    }
	    fclose(stream);

	    // Open file
	    stream = fopen(list_filename, "rt");
	    if (!stream) bye("Cost file cannot be opened");
	    fscanf(stream, "%d:%d\n", &vertices, &colors);


        L.resize(vertices, set<int> ());

	    // Read lists
	    for (int v = 0; v < vertices; v++) {
		    int list_size;
		    fscanf(stream, "%d", &list_size);
		    if (list_size < 1 || list_size > colors) bye("Error reading lists!");

		    // Read a list
		    int last_read = -1;
		    for (int s = 0; s < list_size; s++) {
			    int element;
			    fscanf(stream, "%d", &element);
			    if (element <= last_read || element >= colors) bye("Error reading lists!");
			    last_read = element;
                L[v].insert(element);
		    }
	    }
	    fclose(stream);

    };

};

class Node {

    public:

    vector<int> colored;
    int colored_vertices;
    vector<bool> used_colors;
    int local_bound;

    Node (int n, int c) : colored(n,-1), colored_vertices(0), used_colors(c,false), local_bound(0) {};

};

int main (int argc, char **argv) {

	char* filename = argv[1];
    char arg1[300], arg2[300];
    strncpy(arg1,filename,300);
    strcat(arg1,".graph");
    strncpy(arg2,filename,300);
    strcat(arg2,".list");
    Graph G (arg1, arg2);

    Node* node = new Node (G.vertices, G.colors);
    stack<Node*> st;
    st.push(node);
    int LB = 0;
    int UB = INT_MAX;
    int nodes = 0;

    while (!st.empty()) {

        ++nodes;
        node = st.top();
        st.pop();

        if (node->local_bound >= UB) { // Branch
            delete node;
            continue;
        }

        if (node->colored_vertices == G.vertices) // Branch and maybe update UB
            if (node->local_bound < UB) {
                UB = node->local_bound;
                delete node;
                continue;
            }

        // Select vertex v for branching
        int best_i = -1;
        int best_n1 = 0;
        int best_n2 = 0;
        for (int v = 0; v < G.vertices; ++v) {

            if (node->colored[v] > -1) continue;

            int n1 = 0; //Number of different colors which are unavailable
            for (int k = 0; k < G.colors; ++k) {
                if (G.L[v].find(k) == G.L[v].end())
                    ++n1;
                else
                    for (int u: G.adj[v])
                        if (node->colored[u] == k)
                            ++n1;
            }

            int n2 = 0; // number of uncolored vertices to which vertex v is adjacent
            for (int u: G.adj[v])
                if (node->colored[u] == -1)
                    ++n2;

            if ((n1 > best_n1) || ((n1 == best_n1) && (n2 > best_n2)) || ((n1 == best_n1) && (n2 == best_n2) && (v < best_i)) ) {
                best_i = v;
                best_n1 = n1;
                best_n2 = n2;
            }
        }

        if (best_i == -1) bye ("Branching error");

        // Branch
        int v = best_i;
        for (int k: G.L[v]) {
            bool available = true;
            for (int u: G.adj[v])
                if (node->colored[u] == k) {
                    available = false;
                    break;
                }
            if (available) {
                Node* n = new Node(*node);
                n->colored[v] = k;
                ++n->colored_vertices;
                if (!n->used_colors[k]) {
                    n->used_colors[k] = true;
                    ++n->local_bound;
                }
                st.push(n);
            }
        }


        delete node;

    }

    cout << "Opt = " << UB << endl;
    cout << "Nodes = " << nodes << endl;
    
}
