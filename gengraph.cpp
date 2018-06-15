/*
 * GENGRAPH - Generate a random graph
 * Made in 2018 by Daniel Severin and Mauro Lucci
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>

using namespace std;

/* FUNCTIONS */

/*
 * bye - finish executing and show a message
 */
void bye(char *string)
{
	cout << string << endl;
	exit(1);
}

/*
* urnd - generate numbers with uniform random distribution
*  if flag=false, the interval is [a, b]
*  if flag=true, the interval is [a, b)
*/
float urnd(float a, float b, bool flag)
{
	return a + rand() * (b - a) / (float)(RAND_MAX + (flag ? 1 : 0));
}

/*
 * main - Main program
 */
int main(int argc, char **argv)
{
	cout << "GENGRAPH - Generate a random graph." << endl;

	if (argc < 4) {
		cout << "Usage: gengraph file.graph vertices density" << endl;
		cout << "It generates a random graph with the given vertices and density (in percentage)." << endl;
		cout << "For instances similar to transport ones, use 38.8 as density" << endl;
		bye("Bye!");
	}

	char *filename = argv[1];
	int vertices = atoi(argv[2]);
	if (vertices < 4 || vertices > 10000) bye("Number of vertices out range!");
	float density = atof(argv[3]) / 100.0;
	if (density < 0.0 || density > 1.0) bye("Density out range!");

	int clique_size = vertices * (vertices - 1) / 2;
	int *edge_u = new int[clique_size];
	int *edge_v = new int[clique_size];

	/* generate the random graph and the number of colors */
	int edges = 0;
	for (int u = 0; u < vertices - 1; u++) {
		for (int v = u + 1; v < vertices; v++) {
			if (urnd(0, 1, false) <= density) {
				edge_u[edges] = u;
				edge_v[edges] = v;
				edges++;
			}
		}
	}
	cout << "Statistics:" << endl;
	float realdensity = 100.0 * (float)edges / (float)clique_size;
	cout << "  |V| = " << vertices << ", |E| = " << edges << " (density = " << realdensity << "%)." << endl;

	/* write graph file */
	FILE *stream = fopen(filename, "wt");
	if (!stream) bye("Graph file cannot be created");
	fprintf(stream, "%d:%d\n", vertices, edges);
	for (int e = 0; e < edges; e++) fprintf(stream, "%d,%d\n", edge_u[e], edge_v[e]);
	fclose(stream);

	/* free memory */
	delete[] edge_v;
	delete[] edge_u;
	return 0;
}
