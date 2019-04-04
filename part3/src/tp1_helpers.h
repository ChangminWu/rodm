#include <cmath>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <iterator>
#include <vector>
#include <algorithm>

using namespace std;

typedef struct{
	unsigned long v, e;
} graphSize;

typedef struct{
	unsigned long u, v;
} edge;

unsigned long max3(unsigned long a, unsigned long b, unsigned long c);
bool sortvec(const edge& a, const edge& b);
graphSize size_of_graph(char* file_name);
void digraph2graph(char* infile_name, char* outfile_name);


inline unsigned long max3(unsigned long a, unsigned long b, unsigned long c){
    a = (a>b) ? a : b;
    return (a>c) ? a : c;
}

bool sortvec(const edge& a, const edge& b){ 
    if (a.u<b.u) {
		if (a.v < b.v) {
			return 1;
		}
		else {
			return 0;
		}
	}
	else {
		return 0;
	}
	//return a.u < b.u ? a.v < b.v ? 1 : 0 : 0;
}

graphSize size_of_graph(char* file_name){
	FILE *file = fopen(file_name, "r");

	unsigned long s;
	unsigned long t;
	unsigned long nbr_node = 0;
	unsigned long nbr_edge = 0;

	while (fscanf(file, "%lu %lu", &s, &t)==2) {
		nbr_node = max3(nbr_node, s, t);
		nbr_edge += 1;
	}

	nbr_node += 1;

	graphSize gs;
	gs.v = nbr_node;
	gs.e = nbr_edge;
	return gs;
}

void digraph2graph(char* infile_name, char* outfile_name){
	FILE *file_in = fopen(infile_name, "r");
	FILE *file_out = fopen(outfile_name, "w");

	unsigned long s;
	unsigned long t;
	edge e;
	//unordered_set<<unsigned int, unsigned int>> set_edge;
	vector<edge> vec_edge;


	while (fscanf(file_in, "%lu %lu", &s, &t)==2) {
		if (s < t) {
			e.u = s;
			e.v = t;
		}else if (t < s) {
			e.u = t;
			e.v = s;
		}
		//set_edge.insert(e);
		//vec_edge.push_back(e);	
	}
	//vec_edge.assign( set_edge.begin(), set_edge.end() );
	sort(vec_edge.begin(), vec_edge.end(), sortvec); 
	auto last = unique(vec_edge.begin(), vec_edge.end());
	vec_edge.erase(last, vec_edge.end()); 

	for (auto e : vec_edge) {
		fprintf(file_out, "%lu %lu", e.u, e.v);
	}
}

