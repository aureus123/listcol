/*
 * GENSETCOVERINST - Generate an instance for the Set Covering Problem
 * Made in 2021 by Daniel Severin
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>

using namespace std;

#define MAXVERTPERCOLOR 15
#define MAXCOLORPERVERT 200
//#define RAIL_INSTANCE

/* FUNCTIONS */

/*
 * bye - finish executing and show a message
 */
void bye(char *string)
{
	cout << string << endl;
	exit(-1);
}

/*
 * main - Main program
 */
int main(int argc, char **argv)
{
	cout << "GENSETCOVERINST - Generate an instance for the Set Covering Problem." << endl;

	if (argc < 2) {
		cout << "Usage: gensetcoverinst file" << endl;
		cout << "It takes file.txt and generates file.graph, file.cost and file.list." << endl;
		bye("Bye!");
	}

	char *filename = argv[1];

	/* open graph file */
	char filename_extension[210];
	strncpy(filename_extension, filename, 200);
	strcat(filename_extension, ".txt");
	FILE *stream = fopen(filename_extension, "rt");
	if (!stream) bye("SC file cannot be opened");
	int vertices, colors;
	fscanf(stream, "%d %d\n", &vertices, &colors);

	/* do not accept graph of less than 4 vertices */
	if (vertices < 4) bye("Number of vertices out range!");

	/* ask for memory */
	int *weight = new int[colors];
	int *L_size = new int[vertices];
	for (int v = 0; v < vertices; v++) L_size[v] = 0;


#ifndef RAIL_INSTANCE
	/* read weights */
	for (int j = 0; j < colors; j++) {
		int w;
		fscanf(stream, "%d", &w);
		if (w < 1 || w > 100) {
			cout << "Error reading color " << j + 1 << ", weight out of bounds!" << endl;
			bye("Bye!");
		}
		weight[j] = w;
	}

	/* read list of colors */
	bool *col = new bool[colors];
	int **L_set = new int*[vertices];
	for (int v = 0; v < vertices; v++) {
		int size;
		fscanf(stream, "%d", &size);
		if (size < 1 || size > MAXCOLORPERVERT ) {
			cout << "Error reading vertex " << v + 1 << ", size out of bounds!" << endl;
			bye("Bye!");
		}
		for (int j = 0; j < colors; j++) col[j] = false;
		for (int i = 0; i < size; i++) {
			int j;
			fscanf(stream, "%d", &j);
			if (j < 1 || j > colors) {
				cout << "Error reading vertex " << v + 1 << ", color out of bounds!" << endl;
				bye("Bye!");
			}
			col[j - 1] = true;
		}
		L_set[v] = new int[size];
		L_size[v] = size;
		int i = 0;
		for (int j = 0; j < colors; j++) if (col[j]) L_set[v][i++] = j;
		if (i != size) bye("Internal error reading list of colors!");
	}
#else
	/* read Vk's */
	int *Vk_size = new int[colors];
	int **Vk_set = new int*[colors];
	for (int j = 0; j < colors; j++) {
		int w, size;
		fscanf(stream, "%d %d", &w, &size);
		if (w < 1 || w > 10) {
			cout << "Error reading color " << j + 1 << ", weight out of bounds!" << endl;
			bye("Bye!");
		}
		if (size < 1 || size > MAXVERTPERCOLOR ) {
			cout << "Error reading color " << j + 1 << ", size out of bounds!" << endl;
			bye("Bye!");
		}
		weight[j] = w;

		Vk_set[j] = new int[size];
		Vk_size[j] = size;
		for (int i = 0; i < size; i++) {
			int v;
			fscanf(stream, "%d", &v);
			if (v < 1 || v > vertices) {
				cout << "Error reading color " << j + 1 << ", vert out of bounds!" << endl;
				bye("Bye!");
			}
			Vk_set[j][i] = v - 1;
			L_size[v - 1]++;
		}
	}
#endif
	fclose(stream);

	cout << "Statistics:" << endl;
	cout << "  |V| = " << vertices << ", |E| = 0, |C| = " << colors << "." << endl;

	/* write costs to a file */
	strncpy(filename_extension, filename, 200);
	strcat(filename_extension, ".cost");
	stream = fopen(filename_extension, "wt");
	if (!stream) bye("Cost file cannot be created");
	fprintf(stream, "%d\n", colors);
	for (int j = 0; j < colors; j++) fprintf(stream, "%d ", weight[j]);
	fprintf(stream, "\n");
	fclose(stream);

	/* write edge-less graph to a file */
	strncpy(filename_extension, filename, 200);
	strcat(filename_extension, ".graph");
	stream = fopen(filename_extension, "wt");
	if (!stream) bye("Graph file cannot be created");
	fprintf(stream, "%d:0\n", vertices);
	fclose(stream);

	/* write list of vertices to a file */
	strncpy(filename_extension, filename, 200);
	strcat(filename_extension, ".list");
	stream = fopen(filename_extension, "wt");
	if (!stream) bye("List file cannot be created");
	fprintf(stream, "%d:%d\n", vertices, colors);
	
	for (int v = 0; v < vertices; v++) {
		fprintf(stream, "%d  ", L_size[v]);
#ifndef RAIL_INSTANCE
		for (int i = 0; i < L_size[v]; i++) fprintf(stream, "%d ", L_set[v][i]);
#else
		for (int j = 0; j < colors; j++) {
			for (int i = 0; i < Vk_size[j]; i++) {
				if (v == Vk_set[j][i]) fprintf(stream, "%d ", j);
			}
		}
#endif
		fprintf(stream, "\n");
	}
	fclose(stream);

	/* free memory */
#ifndef RAIL_INSTANCE
	for (int v = 0; v < vertices; v++) delete[] L_set[v];
	delete[] col;
#else
	for (int j = 0; j < colors; j++) delete[] Vk_set[j];
	delete[] Vk_size;
#endif
	delete[] L_size;
	delete[] weight;
	return 0;
}
