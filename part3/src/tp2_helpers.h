#include <iostream>
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
#include <numeric>

using namespace std;

struct GraEdgeList{
    vector< pair <unsigned long, unsigned long> > edge_list;
    map< unsigned long, int > degrees;
}; 

struct CoreNumber{
    map< unsigned long, int > cores;
    int graph_core;
};

GraEdgeList load_graph_kcore(char* file_graph) {
    FILE *file = fopen(file_graph, "r");
    unsigned long s,t;
    vector< pair <unsigned long, unsigned long> > edges;
    map<unsigned long, int> degree_list;

    while (fscanf(file, "%lu %lu", &s, &t)==2) {
        if (s != t) {
            edges.push_back(make_pair(s, t));
        }
    }
    fclose(file);

    sort(edges.begin(), edges.end(), [](const pair<unsigned long, unsigned long> &a, const pair<unsigned long, unsigned long> &b)
                                        {return a.first==b.first?a.second<b.second:a.first<b.first;});

    for (auto&edge : edges) {
        if (!degree_list.count(edge.first)) {
            degree_list[edge.first] = 1;
        } else {degree_list[edge.first] += 1;}
        
        if (!degree_list.count(edge.second)) {
            degree_list[edge.second] = 1;
        } else {degree_list[edge.second] += 1;}
    }

    GraEdgeList graph;
    graph.degrees = degree_list; 
    graph.edge_list = edges;    
    degree_list.clear();
    edges.clear();

    return graph;
}

CoreNumber kcore(GraEdgeList graph){
    unsigned long nbr_node = graph.degrees.size();
    vector < pair<unsigned long , int> > degree_vec;
    int max_core = 0;
    map<unsigned long, int> core_list; 
    unsigned long s;

    for (auto& x: graph.degrees) {
        core_list[x.first] = 0;
    }

    // for (auto itr = graph.degrees.begin(); itr != graph.degrees.end(); ++itr) {
    //     degree_vec.push_back(*itr);
    // }
    copy(graph.degrees.begin(), graph.degrees.end(), back_inserter<vector < pair <unsigned long , int > > > (degree_vec));

    sort(degree_vec.begin(), degree_vec.end(), [](const pair<unsigned long, int> &a, const pair<unsigned long, int> &b)
                                        {return a.second<b.second;});

    unsigned long iter = nbr_node;

    while (!graph.degrees.empty()) {
        s = degree_vec[0].first;
        
        max_core = max(max_core, degree_vec[0].second);
        core_list[s] = max_core;

        for (auto &edge: graph.edge_list) {
            if (edge.first == s) {
                graph.edge_list.erase(remove(graph.edge_list.begin(), graph.edge_list.end(), edge), graph.edge_list.end());
            } else if (edge.second == s) {
                graph.edge_list.erase(remove(graph.edge_list.begin(), graph.edge_list.end(), edge), graph.edge_list.end());
            }
        }
        graph.degrees.clear();
        degree_vec.clear();

        for (auto& edge : graph.edge_list) {
            if (!graph.degrees.count(edge.first)) {
                graph.degrees[edge.first] = 1;
            } else { graph.degrees[edge.first] += 1;}
        
            if (!graph.degrees.count(edge.second)) {
                graph.degrees[edge.second] = 1;
            } else { graph.degrees[edge.second] += 1;}
        }

        for (auto& x: core_list){
            if (!graph.degrees.count(x.first) && x.second == 0) {
                core_list[x.first] = max_core;
            }
        }

        copy(graph.degrees.begin(), graph.degrees.end(), back_inserter<vector < pair <unsigned long , int > > > (degree_vec));
        sort(degree_vec.begin(), degree_vec.end(), [](const pair<unsigned long, int> &a, const pair<unsigned long, int> &b)
                                        {return a.second<b.second;});
    }
    
    CoreNumber core_numbers;
    core_numbers.cores = core_list;
    core_numbers.graph_core = max_core;
    return core_numbers;
}

void post_process(char* file_name, map< unsigned long, int > degree_list, map<unsigned long, int> core_list) {
    FILE *file = fopen(file_name, "w");
    vector<unsigned long> nodes;
    for (auto& x: degree_list) {
        nodes.push_back(x.first);
    }

    for (auto& node: nodes) {
        fprintf(file, "%lu %d %d", node, degree_list[node], core_list[node]);
    }
    fclose(file);
}

// void load_graph(char* file_graph, char* file_name, char* new_file_name) {
//     FILE *gfile = fopen(file_graph, "r");
//     FILE *nfile = fopen(file_name, "r");
//     FILE *nnfile = fopen(new_file_name, "w");
//     unsigned long s, t;
//     char str[100];
//     vector< pair <unsigned long, unsigned long> > edge_list;
//     map<unsigned long, unsigned long> reindex;
//     unsigned long max_idx = 0;

//     // while (fscanf(gfile, "%lu %lu", &s, &t)==2) {
//     //     if (!reindex.count(s)) {
//     //         reindex[s] = max_idx;
//     //         max_idx += 1;
//     //     }

//     //     if (!reindex.count(t)) {
//     //         reindex[t] = max_idx;
//     //         max_idx += 1;
//     //     }

//     //     edge_list.push_back(make_pair(reindex[s], reindex[t]));
    
//     // }
//     // fclose(gfile);
//     // cout << max_idx << endl;
//     // AdjMatrix adj_matrix (max_idx, vector<bool> (max_idx, 0) );
//     // for (auto& edge: edge_list) {
//     //     adj_matrix[edge.first][edge.second] = 1;
//     // }

//     while (fscanf(nfile, "%lu %s", &t, &str)==2) {
//         cout << str << endl;
//         fprintf(nnfile, "%lu %s\n", t, str);
//     }

//     fclose(nfile);
//     fclose(nnfile);
// }

// TransMatrix init_transition(AdjMatrix adj_matrix) {
//     unsigned long nbr_node = adj_matrix.size();
//     unsigned long sum_of_row;

//     TransMatrix transition(nbr_node, vector<double> (nbr_node, 0.0) );
//     for (unsigned long i=0; i<nbr_node; i++) {
//         sum_of_row = accumulate(adj_matrix[i].begin(), adj_matrix[i].end(), 0);
//         if (sum_of_row) {
//             for (unsigned long j=0; j<nbr_node; j++) {
//                 transition[i][j] = static_cast<double>(adj_matrix[i][j]) / static_cast<double>(sum_of_row);
//             }
//         } else {
//             for (unsigned long j=0; j<nbr_node; j++) {
//                 transition[i][j] = 1.0 / static_cast<double>(nbr_node);
//             }
//         }
//         //sum_of_row = 0;
//         //for_each(adj_matrix[i].begin(), adj_matrix[i].end(), [&] (unsigned short n) {sum_of_row += n;});
//     }
//     return transition;
// }

// vector<double> matvec_prod(TransMatrix transition, vector<double> parray) {
//     unsigned long nbr_node = parray.size();
//     vector<double> sum(nbr_node, 0.0);
//     for (unsigned long j=0; j<nbr_node; j++) {
//         double column_sum = 0.0;
//         for (unsigned long i=0; i<nbr_node; i++) {
//             column_sum += transition[i][j] * parray[i];
//         }
//         sum[j] = column_sum;
//     }
//     return sum;
// }

// vector<double> normalize(vector<double> parray) {
//     unsigned long nbr_node = parray.size();
//     double norm = 0.0;
//     for (unsigned long i=0; i < nbr_node; i++) {
//         norm += abs(parray[i]);
//     }
//     for (unsigned long i=0; i < nbr_node; i++) {
//         parray[i] += (1.0-norm) / static_cast<double>(nbr_node);
//     }
//     return parray;
// } 

// vector<double> vecscalar_prod(vector<double> parray, double alpha) {
//     unsigned long nbr_node = parray.size();
//     for (unsigned long i=0; i<nbr_node; i++) {
//         parray[i] *= alpha;
//     }
//     return parray;
// }

// double difference(vector<double> parray1, vector<double> parray2) {
//     unsigned long nbr_node = parray1.size();
//     double sum = 0.0;
//     for (unsigned long i=0; i<nbr_node; i++) {
//         sum += abs(parray1[i]-parray2[i]);
//     }
//     return sum;
// }

// vector<double> pagerank(AdjMatrix adj_matrix, double alpha, double tol, unsigned int max_iter, char* file_name) {
//     unsigned long nbr_node = adj_matrix.size();
//     double sum;
//     vector<double> pvalue(nbr_node, 1.0 / nbr_node);
//     vector<double> ident(nbr_node, 1.0 / nbr_node);
//     vector<double> new_pvalue;
//     FILE *file = fopen(file_name, "w");

//     TransMatrix transition = init_transition(adj_matrix); 
//     for (unsigned int iter=0;  iter<max_iter; iter++) {
//         sum = 0.0;
//         new_pvalue = matvec_prod(transition, pvalue);
//         for (unsigned long i=0; i<nbr_node; i++) {
//             new_pvalue[i] = new_pvalue[i] * (1.0-alpha) + alpha * (1.0 /nbr_node);
//         }
//         new_pvalue = normalize(new_pvalue);
//         for (unsigned long i=0; i<nbr_node; i++) {
//             sum += new_pvalue[i] * new_pvalue[i];
//         }
//         fprintf(file, "%u %f\n", iter, sum);

//         if (difference(new_pvalue, pvalue)<tol) {
//             break;
//         }
//         pvalue = new_pvalue;
//         new_pvalue.clear();
//     }
//     fclose(file);
//     return pvalue;
// }

// void post_process(char* infile_name, char* outfile_name, vector<double> pvalue){
//     FILE *ifile = fopen(infile_name, "r");
//     FILE *ofile = fopen(outfile_name, "w");
//     unsigned long s;
//     unsigned long nbr_node = pvalue.size();
//     char str[100];
//     vector<string> pagename;
//     vector< pair <unsigned long, double> >  pvalue_id(nbr_node);

//     while (fscanf(ifile, "%lu %s", &s, &str)==2) {
//         pagename.push_back(str);
//     }
//     fclose(ifile);

//     for (unsigned long i=0; i<nbr_node; i++) {
//         pvalue_id[i] = make_pair(i, pvalue[i]);
//     }

//     sort(pvalue_id.begin(), pvalue_id.end(), [](const pair<unsigned long, double> &a, const pair<unsigned long, double> &b)
// 											{return a.second<b.second;});

//     fprintf(ofile, "Five Pages with lowest pagerank values:\n");
//     for (int i=0; i<5; i++) {
//         fprintf(ofile, "%s %f\n", pagename[pvalue_id[i].first], pvalue_id[i].second);
//     } 

//     fprintf(ofile, "Five Pages with highest pagerank values:\n");
//     for (int i=nbr_node-1; i>nbr_node-6; i--) {
//         fprintf(ofile, "%s %f\n", pagename[pvalue_id[i].first], pvalue_id[i].second);
//     } 

//     fclose(ofile);
// }



