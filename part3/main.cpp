#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <iterator>
#include <vector>
//#include "src/tp1_helpers.h"
#include "src/tp2_helpers.h"

using namespace std;

int main(int argc, char** argv) {   
    /* TP1 */
    
    /* declare variables */
    // GraphSize gs;
    // time_t t1,t2;
    // map<unsigned long, unsigned long> degrees;
    // unsigned long Q;
    // EdgeList g1;
    // AdjMatrix g2;
    // AdjArray g3;
    // vector<unsigned long> size_comp;
    // unsigned long diam, num_tri;

    /* exercise 2 */
    //t1 = time(NULL);
    //gs = size_of_graph(argv[1]);
    //cout << "Number of nodes: " << gs.v << endl;
    //cout << "Number of edges: " << gs.e << endl;
    //t2 = time(NULL);
    //printf("- Overall time for checking graph size = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));

    /* exercise 3 */
    //digraph2graph(argv[1], argv[2]);

    /* check new size of graph (undirected with sorted edges, self-loop removed and re-index of nodes) */
    // gs = size_of_graph(argv[2]);
    // cout << "Number of new nodes: " << gs.v << endl;
    // cout << "Number of new edges: " << gs.e << endl;

    /* exercise 4 */
    // degrees = degree_list(argv[2], argv[3]);

    /* exercise 5 */
    // t1 = time(NULL);
    // Q = degree_sum(argv[2], degrees);
    // t2 = time(NULL);
    // cout << "Final Sum: " << Q << endl; 
    // printf("- Overall time for Q value = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));

    /* exercise 6 */
    // degree_dist(argv[4], degrees);

    /* exercise 7 */
    // t1 = time(NULL);
    // g1 = load_edge_list(argv[2]);
    // g1.clear();
    // t2 = time(NULL);
    // printf("- Overall time for EdgeList = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
    
    // t1 = time(NULL);
    // g2 = load_adj_matrix(argv[2], gs.v);
    // g2.clear();
    // t2 = time(NULL);
    // printf("- Overall time for AdjMatrix = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));

    // t1 = time(NULL);
    // g3 = load_adj_array(argv[2], degrees);
    // g3.ind.clear();
    // g3.l.clear();
    // t2 = time(NULL);
    // printf("- Overall time for AdjList = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
    
    /* exercise 8 */
    // size_comp = connected_component(g3);
    // for (auto&size: size_comp) {
    //      cout << size << "; ";
    // }
    // cout << gs.v << endl;
    
    /* exercise 8 */
    // diam = diameter_estimation_2(g3);
    // cout << diam << endl;
    
    /* exercise 9 */
    // num_tri = tiangles(g1, g3, degrees);
    // cout << num_tri << endl;
    
    /* ---------------------------------------------------------------------------------- */
    /* TP2 */

    /* declare variables */
    double alpha = 0.15;
    double tol = 0.0000001;
    unsigned int max_iter = 100;
    vector<double> pvalue;
    CoreNumber core_number;
    GraEdgeList graph;
    time_t t1,t2;

    /*exercise 3: Kcore*/
    t1 = time(NULL);

    graph = load_graph_kcore(argv[1]);
    cout << graph.degrees.size() << endl;

    core_number = kcore(graph);
    t2 = time(NULL);
    cout << core_number.graph_core << endl;
    cout << core_number.cores.size() << endl;
    printf("- Overall time for KCore = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
    
    /*exercise 1: pagerank*/
    // load_graph(argv[1], argv[2], argv[3]);
    // pvalue = pagerank(adj_matrix, alpha, tol, max_iter, argv[4]);
    // post_process(argv[3], argv[5], pvalue);
    return 0;
}

