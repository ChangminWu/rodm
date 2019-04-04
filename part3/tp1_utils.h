#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <chrono>

using namespace std;

typedef struct{
    unsigned long v;
    unsigned long e;
} gStats;

unsigned long max3(unsigned long a, unsigned long b, unsigned long c);
gStats number_of_nodes(ifstream& file);
void clean_data(ifstream& file_in, ofstream& file_out);

inline unsigned long max3(unsigned long a, unsigned long b, unsigned long c){
    a = (a>b) ? a : b;
    return (a>c) ? a : c;
}

gStats number_of_nodes(ifstream& file) {
    char ch;
    unsigned long nbr_node = 0;
    unsigned long nbr_edge = 0;
    unsigned long s;
    unsigned long t;
    string line;

    auto start = chrono::high_resolution_clock::now();
    while (getline(file, line)) {
        if (line.length() != 0 && line[0] != '#') {
            istringstream is(line);
            is >> s >> t;
            nbr_node = max3(nbr_node, s, t);
            nbr_edge += 1;
        }
    }

    file.close();

    cout << "Number of nodes: " << nbr_node << "\n";
    cout << "Number of edges: " << nbr_edge << "\n";
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end - start;
    cout << "Taking time: " << elapsed.count() << " s" << endl;

    gStats gs;
    gs.v = nbr_node;
    gs.e = nbr_edge; 

    return gs;
}

void clean_data(ifstream& file_in, ofstream& file_out) {
    char ch;
    unsigned long s;
    unsigned long t;
    string line;

    while (getline(file_in, line)) {
        if (line.length() != 0 && line[0] != '#') {
            istringstream is(line);
            is >> s >> t;
            if (s < t) {
                file_out << s << t << "\n";
            }
        } else {
            file_out << line;
        }
    }
    file_in.close();
    file_out.close();
}








