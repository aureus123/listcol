/*
 * GENMUINST - Generate a random instance of the mu-coloring Problem from a graph
 * Made in 2018-2020 by Daniel Severin
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>

#ifdef VISUALC
#include <ctime>
#endif

using namespace std;

/* CONSTANTS */

#define INFDIST 9999999

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
	return a + rand() * (b - a) / ((float)RAND_MAX + (flag ? 1.0 : 0.0));
}

/*
 * main - Main program
 */
int main(int argc, char **argv)
{
	cout << "GENMUINST - Generate a random instance of mu-coloring Problem." << endl;

	if (argc < 2) {
		cout << "Usage: genmuinst file [colors]" << endl;
		cout << "It takes file.graph and generates file.cost and file.list." << endl;
		cout << "If \"colors\" is given, mu(v) = uniform distribution between 2 and colors." << endl;
//		cout << "Otherwise, mu(v) = uniform distribution between 0.5 * |Delta| and |Delta|." << endl;
//		cout << "Otherwise, mu(v) = uniform distribution between 0.1 * |V| and |V|." << endl;
		cout << "Otherwise, mu(v) = uniform distribution between 1 and |V|." << endl;
		bye("Bye!");
	}

	char *filename = argv[1];

	/* open graph file */
	char filename_extension[210];
	strncpy(filename_extension, filename, 200);
	strcat(filename_extension, ".graph");
	FILE *stream = fopen(filename_extension, "rt");
	if (!stream) bye("Graph file cannot be opened");
	int vertices, edges;
	fscanf(stream, "%d:%d\n", &vertices, &edges);

	/* do not accept graph of less than 4 vertices or stable sets */
	if (vertices < 4) bye("Number of vertices out range!");
	if (edges < 1 || edges > vertices*(vertices - 1) / 2) bye("Number of edges out of range!");

	/* read edges for obtaining degree of vertices */
	int *degrees = new int[vertices];
	for (int v = 0; v < vertices; v++) degrees[v] = 0;
	for (int e = 0; e < edges; e++) {
		int u, v;
		fscanf(stream, "%d,%d\n", &u, &v);
		if (u < 0 || u >= v || v >= vertices) {
			cout << "Error reading edge " << e + 1 << "!" << endl;
			bye("Bye!");
		}
		degrees[u]++;
		degrees[v]++;
	}
	fclose(stream);

	/* set number of colors */
	int colors = -1;
	if (argc > 2) {
		colors = atoi(argv[2]);
		if (colors < 3 || colors > 50000) bye("Number of colors out range!");
	}

	cout << "Statistics:" << endl;
	int clique_size = vertices * (vertices - 1) / 2;
	float density = 100.0 * (float)edges / (float)clique_size;
	cout << "  |V| = " << vertices << ", |E| = " << edges << " (density = " << density << "%)." << endl << endl;

	/* generate mu's */
	srand(time(0));
	int *mu = new int[vertices];
	int max_colors_generated = -1;
	for (int v = 0; v < vertices; v++) {
		int m;
		if (colors == -1) {
/*			m = (int)(urnd(0.5 * (float)degrees[v], 1.0 * (float)degrees[v], false) + 0.1);*/
/*			m = (int)(urnd(0.1 * (float)vertices, 1.0 * (float)vertices, false) + 0.1);
			if (m < 2) m = 2; */
			m = (int)(urnd(1.0, 1.0 * (float)vertices, false) + 0.1);
		}
		else m = (int)(urnd(2.0, (float)colors, false) + 0.1);
		if (m > max_colors_generated) max_colors_generated = m;
		mu[v] = m;
	}
	cout << "Colors generated = " << max_colors_generated << endl;

	/* write cost file */
	strncpy(filename_extension, filename, 200);
	strcat(filename_extension, ".cost");
	stream = fopen(filename_extension, "wt");
	if (!stream) bye("Cost file cannot be created");
	fprintf(stream, "%d\n", max_colors_generated);
	for (int k = 0; k < max_colors_generated; k++) fprintf(stream, "1 ");
	fprintf(stream, "\n");
	fclose(stream);

	/* write list file */
	strncpy(filename_extension, filename, 200);
	strcat(filename_extension, ".list");
	stream = fopen(filename_extension, "wt");
	if (!stream) bye("List file cannot be created");
	fprintf(stream, "%d:%d\n", vertices, max_colors_generated);
	for (int v = 0; v < vertices; v++) {
		fprintf(stream, "%d  ", mu[v]);
		for (int s = 0; s < mu[v]; s++) fprintf(stream, "%d ", s);
		fprintf(stream, "\n");
	}
	fclose(stream);

	delete[] degrees;
	return 0;
}
