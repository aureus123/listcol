#include "io.h"
#include <stdio.h>
#include <vector>
#include <utility>
#include <iostream>

/* for linux users: do not define VISUALC */
#ifndef VISUALC
#include <unistd.h>
#include <sys/times.h>
#else
#include <windows.h>
#endif

// set_color - change color of text
void set_color(int color)
{
#ifdef VERBOSE
#ifdef VISUALC
	SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), color);
#endif
#endif
    return;
}

// bye - finish executing and show a message
void bye(string str)
{
#ifdef VERBOSE
	set_color(12);
	cout << str << endl;
	set_color(7);
#endif
	exit(1);
}

// read_graph - read a graph in the following format:
//   in the first line, the number of vertices and edges separated by a colon ":"
//   then, for each line, the endpoints of an edge (u,v) where u < v
//   note that vertices starts from 0, i.e. 0 < v < |V|-1
//   example for a diamond graph:
//     4:5
//     0,1
//     0,2
//     0,3
//     1,2
//     1,3
int read_graph(char *filename, vector<vector<int> >& adj_list)
{
	// Open file 
	FILE *stream = fopen(filename, "rt");
	if (!stream) bye("Graph file cannot be opened");
    int vertices, edges;
	fscanf(stream, "%d:%d\n", &vertices, &edges);

	// Do not accept graph of less than 4 vertices or stable sets
	if (vertices < 4) bye("Number of vertices out range!");
	if (edges < 1 || edges > vertices*(vertices - 1) / 2) bye("Number of edges out of range!");

    adj_list.resize(vertices);

	// Read edges
	for (int e = 0; e < edges; e++) {
		int u, v;
		fscanf(stream, "%d,%d\n", &u, &v);
		if (u < 0 || u >= v || v >= vertices) {
			cout << "Error reading edge " << e + 1 << "!" << endl;
			bye("Bye!");
		}
        adj_list[u].push_back(v);
        adj_list[v].push_back(u);
	}
	fclose(stream);

    return edges;
}

// read_cost - read costs of colors in the following format:
//   in the first line, the number of colors
//   then, the costs of color 0, 1, etc, in succesive order
//   example for C = {0, 1, 2} with costs c_0 = 5, c_1 = 2, c_2 = 8:
//     3
//     5 2 8
void read_cost(char *filename, vector<int>& costs_list)
{
	/* open file */
	FILE *stream = fopen(filename, "rt");
	if (!stream) bye("Cost file cannot be opened");
    int colors;
	fscanf(stream, "%d\n", &colors);

	/* do not accept less than 2 colors */
	if (colors < 2) bye("Number of colors out range!");

    costs_list.resize(colors);

	/* read costs */
	for (int k = 0; k < colors; k++) {
		int ck;
		fscanf(stream, "%d", &ck);
		if (ck < 0) bye("Color cost must be non negative!");
		costs_list[k] = ck;
	}
	fclose(stream);

    return;
}

// read_list - read list of colors per vertex in the following format:
//   in the first line, the number of vertices and colors separated by a colon ":"
//   then, for each line, the cardinal of L(v) followed by the elements of L(v) in increasing order
//   example for |V| = 3, |C| = 5 and L(0) = {1, 2}, L(1) = {0, 2, 3}, L(2) = {0, 1, 4}:
//     3:5
//     2  1 2
//     3  0 2 3
//     3  0 1 4
void read_list(char *filename, int vertices, int colors, vector<vector<int> >& colors_list)
{
	/* open file */
	FILE *stream = fopen(filename, "rt");
	if (!stream) bye("Cost file cannot be opened");

	/* check vv and cc */
	int vv, cc;
	fscanf(stream, "%d:%d\n", &vv, &cc);
	if (vv != vertices) bye("Number of vertices mismatch!");
	if (cc != colors) bye("Number of colors mismatch!");

    colors_list.resize(vertices);

	/* read lists */
	for (int v = 0; v < vertices; v++) {
		int list_size;
		fscanf(stream, "%d", &list_size);
		if (list_size < 1 || list_size > colors) bye("Error reading lists!");
        colors_list[v].resize(list_size);

		/* read a list */
		int last_read = -1;
		for (int s = 0; s < list_size; s++) {
			int element;
			fscanf(stream, "%d", &element);
			if (element <= last_read || element >= colors) bye("Error reading lists!");
			last_read = element;
            colors_list[v][s] = element;
		}
	}
	fclose(stream);

    return;
}
