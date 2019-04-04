#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <iterator>
#include <vector>
#include "src/tp1_helpers.h"

using namespace std;

int main(int argc, char** argv) {   
    graphSize gs;
    time_t t1,t2,t3;

    t1 = time(NULL);

    gs = size_of_graph(argv[1]);
    cout << "Number of nodes: " << gs.v << endl;
    cout << "Number of edges: " << gs.e << endl;

    t2 = time(NULL);
    printf("- Overall time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));

    digraph2graph(argv[1], argv[2]);

    gs = size_of_graph(argv[2]);
    cout << "Number of nodes: " << gs.v << endl;
    cout << "Number of edges: " << gs.e << endl;

    t3 = time(NULL);
	printf("- Overall time = %ldh%ldm%lds\n",(t3-t2)/3600,((t3-t2)%3600)/60,((t3-t2)%60));

    return 0;
}