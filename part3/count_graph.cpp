#include <stdio.h>
#include <stdlib.h>
#include <time.h>

using namespace std;

#define NLINKS 100000000;

typedef struct{
    unsigned long s;
    unsigned long t;
} edge;

typedef struct{
    unsigned long v;
    unsigned long e;
    edge *edges;
} edgelist;


inline unsigned long max3(unsigned long a, unsigned long b, unsigned long c) {
    a = (a>b) ? a : b;
    return (a>c) ? a : c;
}


edgelist* read_edgelist(char* input) {
    unsigned long el = NLINKS;
    FILE *file = fopen(input, "r");

    edgelist *g = (edgelist*) malloc(sizeof(edgelist));
    g->v = 0;
    g->e = 0;
    g->edges = (edge*) malloc(sizeof(edge)*el);

    while (fscanf(file, "%lu %lu", &(g->edges[g->e].s), &(g->edges[g->e].t))==2){
        g->v = max3(g->v, g->edges[g->e].s, g->edges[g->e].t);
        if (g->e++==el) {
            el += NLINKS;
            g->edges = (edge*) realloc(g->edges, el*sizeof(edge));
        }
    }
    fclose(file);

    g->v++;
    g->edges = (edge*) realloc(g->edges, g->e*sizeof(edge));
    
    return g;
}

void free_edgelist(edgelist *g) {
    free(g->edges);
    free(g);
}

int main(int argc, char** argv) {
    edgelist* g;
    time_t t1, t2;

    t1 = time(NULL);

    printf("Reading edgelist from file %s\n", argv[1]);
    g = read_edgelist(argv[1]);

    printf("Number of nodes: %lu\n", g->v);
    printf("Number of edges: %lu\n", g->e);

    free_edgelist(g);

    t2 = time(NULL);

    printf("- Overall time = %ldh%ldm%lds\n", (t2-t1)/3600, ((t2-t1)%3600)/60, ((t2-t1)%60));
    return 0;   
}



