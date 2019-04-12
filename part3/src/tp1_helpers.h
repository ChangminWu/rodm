#include <cmath>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <iterator>
#include <vector>
#include <set>
#include <algorithm>
#include <list>
#include <map>

using namespace std;
#define NLINKS 100000;

struct GraphSize{
	unsigned long v, e;
};

typedef pair<unsigned long, unsigned long> Edge;
typedef vector<Edge> EdgeList;
typedef vector<bool> AdjVec; 
typedef vector<AdjVec> AdjMatrix;
typedef vector<unsigned long> AdjList;

struct AdjArray{
	AdjList ind, l;
};

unsigned long max3(unsigned long a, unsigned long b, unsigned long c);
GraphSize size_of_graph(char* file_name);
void digraph2graph(char* infile_name, char* outfile_name);


inline unsigned long max3(unsigned long a, unsigned long b, unsigned long c){
    a = (a>b) ? a : b;
    return (a>c) ? a : c;
}

GraphSize size_of_graph(char* file_name){
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

	GraphSize gs;
	gs.v = nbr_node;
	gs.e = nbr_edge;
	fclose(file);
	return gs;
}

void digraph2graph(char* infile_name, char* outfile_name){
	FILE *file_in = fopen(infile_name, "r");
	FILE *file_out = fopen(outfile_name, "w");

	unsigned long s;
	unsigned long t;
	EdgeList vec_edge;
	map<unsigned long, unsigned long> node_reindex;
	unsigned long max_idx = 0;
	
	while (fscanf(file_in, "%lu %lu", &s, &t)==2) {
		if (s != t) {
			if (!node_reindex.count(s)){
				node_reindex[s] = max_idx;
				max_idx += 1;
			}

			if (!node_reindex.count(t)){
				node_reindex[t] = max_idx;
				max_idx += 1;
			}
		}

		if (node_reindex[s] < node_reindex[t]) {
			vec_edge.push_back(make_pair(node_reindex[s], node_reindex[t]));
		}else if (node_reindex[t] < node_reindex[s]) {
			vec_edge.push_back(make_pair(node_reindex[t], node_reindex[s]));
		}	
	}

	sort(vec_edge.begin(), vec_edge.end(), [](const pair<unsigned long, unsigned long> &a, const pair<unsigned long, unsigned long> &b)
											{return a.first==b.first?a.second<b.second:a.first<b.first;}); 
	
	auto last = unique(vec_edge.begin(), vec_edge.end(), [](const pair<unsigned long, unsigned long> &a, const pair<unsigned long, unsigned long> &b)
											{return a.first==b.first && a.second==b.second;});
	
	vec_edge.erase(last, vec_edge.end()); 

	for (auto& e : vec_edge) {
		//cout << e.first << "; " << e.second << endl;
		fprintf(file_out, "%lu %lu\n", e.first, e.second);
	}
	fclose(file_in);
	fclose(file_out);
}

map<unsigned long, unsigned long> degree_list(char* infile_name, char* outfile_name){
	FILE *file_in = fopen(infile_name, "r");
	FILE *file_out = fopen(outfile_name, "w");

	unsigned long s;
	unsigned long t;
	map<unsigned long, unsigned long> degrees;
	
	while (fscanf(file_in, "%lu %lu", &s, &t)==2) {
		if (degrees.count(s)){
			degrees[s] += 1;
		} else {degrees[s] = 1; }
		if (degrees.count(t)){
			degrees[t] += 1;
		} else {degrees[t] = 1; }
	}

	for (auto& node: degrees) {
		fprintf(file_out, "%lu %lu\n", node.first, node.second);
	}
	fclose(file_in);
	fclose(file_out);
	return degrees;	
}

unsigned long degree_sum(char* file_name, map<unsigned long, unsigned long> degrees){
	FILE *file = fopen(file_name, "r");

	unsigned long s;
	unsigned long t;
	unsigned long Q = 0;

	while (fscanf(file, "%lu %lu", &s, &t)==2) {
		Q += degrees[s] * degrees[t];
	}
	fclose(file);
	return Q;
}

void degree_dist(char* file_name, map<unsigned long, unsigned long> degrees){
	FILE *file = fopen(file_name, "w");

	unsigned long deg;
	map<unsigned long, unsigned long> degree_map;

	for(auto& node: degrees){
		deg = node.second;
		if (degree_map.count(deg)){
			degree_map[deg] += 1;
		}
		else {
			degree_map[deg] = 1;
		}
	}

	for (auto& pair: degree_map){
		fprintf(file, "%lu %lu\n", pair.first, pair.second);
	}
	fclose(file);
}

EdgeList load_edge_list(char* file_name){
	FILE *file = fopen(file_name, "r");
	unsigned long el = NLINKS;

	EdgeList g(el);
	unsigned long s;
	unsigned long t;
	unsigned long nbr_edge = 0;

	while (fscanf(file, "%lu %lu", &s, &t)==2) {
		g[nbr_edge] = make_pair(s,t);
		nbr_edge += 1;
		if (nbr_edge >= el) {
			el += NLINKS;
			g.resize(el);
		}
	}
	fclose(file);
	return g;
}

AdjMatrix load_adj_matrix(char* file_name, unsigned long nbr_node){
	FILE *file = fopen(file_name, "r");
	
	unsigned long s,t;
	AdjMatrix adj_matrix(nbr_node, AdjVec(nbr_node, 0));

	while (fscanf(file, "%lu %lu", &s, &t)==2) {
		adj_matrix[s][t] = 1;
		adj_matrix[t][s] = 1;
	}
	fclose(file);
	return adj_matrix;
}

AdjArray load_adj_array(char* file_name, map<unsigned long, unsigned long> degrees) {
	FILE *file = fopen(file_name, "r");
	unsigned long s,t;
	vector<unsigned long> cumdeg_index(degrees.size()+1);

	cumdeg_index[0] = 0;
	for (auto& node: degrees) {
		cumdeg_index[node.first+1] = cumdeg_index[node.first] + node.second;
	}
	
	AdjArray adj_array;
	AdjList curr_ind(cumdeg_index.size(),0);
	AdjList adj_list(cumdeg_index.back()); 
	
	while (fscanf(file, "%lu %lu", &s, &t)==2) {
		adj_list[cumdeg_index[s] + curr_ind[s]] = t;
		adj_list[cumdeg_index[t] + curr_ind[t]] = s;
		curr_ind[s] += 1;
		curr_ind[t] += 1;
	}
	fclose(file);
	adj_array.ind = cumdeg_index;
	adj_array.l   = adj_list;
	return adj_array;
}

vector<unsigned long> connected_component(AdjArray adj_array) {
	unsigned long nbr_node = adj_array.ind.size()-1;
	set<unsigned long> unvisited;
	vector<unsigned long> unvisited_vect;
	unsigned long size_connected_component;
	vector<unsigned long> size_vect;
	for (unsigned long i=0; i<nbr_node; i++) {
		unvisited.insert(i);
	}
	unsigned long s = rand() % nbr_node;     
	list<unsigned long> queue;
	unvisited.erase(s);
	//unvisited.erase(remove(unvisited.begin(), unvisited.end(), s), unvisited.end());
	queue.push_back(s); 
	while (!unvisited.empty()) {
		size_connected_component = 0;
		while (!queue.empty()) {
			s = queue.front(); 
			queue.pop_front();
			size_connected_component += 1;
			for (unsigned long j = adj_array.ind[s]; j < adj_array.ind[s+1]; j++) {
				if (unvisited.find(adj_array.l[j])!=unvisited.end()) {
					unvisited.erase(adj_array.l[j]);
					queue.push_back(adj_array.l[j]); 
				}
			}
		}
		size_vect.push_back(size_connected_component);
		if (unvisited.size()>0) {
			unvisited_vect.assign(unvisited.begin(), unvisited.end());
			unsigned long new_s = rand() %  unvisited_vect.size();
			unvisited.erase(unvisited_vect[new_s]);
			queue.push_back(unvisited_vect[new_s]);
			unvisited_vect.clear();
		}
	} 
	return size_vect;
}

double average(unsigned short array[], int array_size){
	unsigned long sum=0;
	double avg;
	for (int i=0; i<array_size; i++) {
		sum += array[i];
	}
	avg = sum / static_cast<double>(array_size);
	return avg;
}

double diameter_estimation_1(AdjArray adj_array, int init_node) {
	unsigned long nbr_node = adj_array.ind.size()-1;
	unsigned long init_s, s;
	set<unsigned long> unvisited;
	unsigned short diameter[nbr_node];
	list<unsigned long> queue;
	unsigned short diameters[init_node];


	for (int i = 0; i < init_node; i++) {
		for (unsigned long j=0; j<nbr_node; j++) {
			unvisited.insert(j);
		}
		init_s = rand() % nbr_node;
		diameter[init_s] = 0;
		unvisited.erase(init_s);
		queue.push_back(init_s);
		while (!queue.empty()) {
			s = queue.front(); 
			queue.pop_front();
			for (unsigned long k = adj_array.ind[s]; k < adj_array.ind[s+1]; k++) {
				if (unvisited.find(adj_array.l[k])!=unvisited.end()) {
					diameter[adj_array.l[k]] = diameter[s] + 1; 
					unvisited.erase(adj_array.l[k]);
					queue.push_back(adj_array.l[k]); 
				}
			}
		}
		diameters[i] = *max_element(diameter, diameter+nbr_node); 
		unvisited.clear();
	}
	return 2.0*average(diameters, init_node);
}

unsigned long diameter_estimation_2(AdjArray adj_array) {
	unsigned long nbr_node = adj_array.ind.size()-1;
	unsigned long init_s, s;
	
	set<unsigned long> unvisited;
	unsigned short diameter[nbr_node];
	
	list<unsigned long> queue;
	unsigned long diat_estim;

	init_s = rand() % nbr_node;
	for (int i=0; i<2; i++) {
		for (unsigned long j=0; j<nbr_node; j++) {
			unvisited.insert(j);
		}
		diameter[init_s] = 0;
		unvisited.erase(init_s);
		queue.push_back(init_s);
		while (!queue.empty()) {
			s = queue.front(); 
			queue.pop_front();
			for (unsigned long k = adj_array.ind[s]; k < adj_array.ind[s+1]; k++) {
				if (unvisited.find(adj_array.l[k])!=unvisited.end()) {
					diameter[adj_array.l[k]] = diameter[s] + 1; 
					unvisited.erase(adj_array.l[k]);
					queue.push_back(adj_array.l[k]); 
				}
			}
		}
		init_s = s;	
	}
	diat_estim = *max_element(diameter, diameter+nbr_node); // diameter[s];
	return diat_estim;
}

unsigned long tiangles(EdgeList edge_list, AdjArray adj_array, map<unsigned long, unsigned long> degrees){
	vector<pair<unsigned long, unsigned long>> sorted_degree;
	vector<unsigned long> reindex(degrees.size());
	unsigned long u,v,w,count;
	vector<unsigned long> uneighbor, vneighbor, wintsct;

	for (unsigned long i=0; i<degrees.size(); i++){
		sorted_degree.push_back(make_pair(i, degrees[i]));
	}
	
	sort(sorted_degree.begin(), sorted_degree.end(), [](const pair<unsigned long, unsigned long> &a, const pair<unsigned long, unsigned long> &b){return a.second>b.second;}); 
	
	for (unsigned long i=0; i<degrees.size(); i++){
		reindex[sorted_degree[i].first] = i;
	}
	
	count = 0;
	for (auto& edge: edge_list) {
		u = edge.first;
		v = edge.second;
		for (unsigned long i = adj_array.ind[u]; i < adj_array.ind[u+1]; i++) {
			w = adj_array.l[i];
			if (reindex[w]<reindex[u] && reindex[w]>reindex[v]){
				uneighbor.push_back(w);
			}else if (reindex[w]>reindex[u] && reindex[v]>reindex[w]) {
				uneighbor.push_back(w);
			}
		}
		for (unsigned long i = adj_array.ind[v]; i < adj_array.ind[v+1]; i++) {
			w = adj_array.l[i];
			if (reindex[w]<reindex[v] && reindex[u]<reindex[w]){
				vneighbor.push_back(w);
			}else if (reindex[w]>reindex[v] && reindex[u]>reindex[w]) {
				vneighbor.push_back(w);
			}
		}
		sort(uneighbor.begin(), uneighbor.end());
		sort(vneighbor.begin(), vneighbor.end());
		set_intersection(uneighbor.begin(),uneighbor.end(),
                          vneighbor.begin(), vneighbor.end(),
                          back_inserter(wintsct));
		//cout << uneighbor.size() << "; " << vneighbor.size() << "; " << wintsct.size() << endl;
		count += wintsct.size();
		uneighbor.clear();
		vneighbor.clear();
		wintsct.clear();
	}
	return count;
}
