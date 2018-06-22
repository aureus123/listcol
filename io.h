#ifndef _IO_H_
#define _IO_H_

#include <vector>
#include <string>

using namespace std;

int read_graph(char *filename, vector<vector<int> >& adj_list);
void read_cost(char *filename, vector<int>& costs_list);
void read_list(char *filename, int vertices, int colors, vector<vector<int> >& colors_list);
void set_color(int color);
void bye(string str);

#endif
