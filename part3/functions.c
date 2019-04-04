#include <stdlib.h>
#include <stdio.h>
#include <time.h>//to estimate the runing time
#include "functions.h"
#include <math.h>  // for abs value fct

#define NLINKS 10 //maximum number of edges for memory allocation, will increase if needed
#define NPOINTS 200 //maximum number of edges for memory allocation, will increase if needed

//compute the maximum of three unsigned long
inline unsigned long max3(unsigned long a, unsigned long b, unsigned long c){
	a=(a>b) ? a : b;
	return (a>c) ? a : c;
}

//Exercise 2
unsigned long ReadNbNodes(FILE *file){
    
    // reset pointer to file start
    rewind(file);

    // loop over file, saving biggest node number
    unsigned long nbNodes = 0;  
    unsigned long leftNode; unsigned long rightNode;  
    int count=0;
    while (fscanf(file,"%lu %lu", &leftNode, &rightNode)==2) {
        nbNodes = max3(nbNodes, leftNode, rightNode);
        count=count+1;
    }
    printf("There are %lu nodes\n", nbNodes + 1);
    printf("There are %d edges \n", count);
    // reset pointer to file start (for the next person to come along)
    rewind(file);

    return nbNodes + 1; // because counting starts at 0 
}

//Exercise 3, DÉSOLÉ, C'EST UN ESSAIS, ÇA MARCHE PAS
void CleaningData(FILE *file){
     FILE *NewFile = fopen("New_email-Eu-core.txt","w") ;
//     FILE *NewFile2 = fopen("NewPrime_email-Eu-core.txt","w") ;

     // loop over file, saving biggest node number
     unsigned long leftNode; unsigned long rightNode;
     unsigned long leftPrime; unsigned long rightPrime;
     int count=0;
     int limit=1;

     while(limit<=25571){
         count=0;
         rewind(file);
         while (count<=limit) {
             fscanf(file,"%lu %lu", &leftNode, &rightNode);
             count=count+1;
         }

         while (fscanf(file,"%lu %lu", &leftPrime, &rightPrime)==2){
             if((leftNode==rightPrime&&rightNode==leftPrime)||(rightNode==rightPrime&&leftNode==leftPrime)){
//                fprintf(NewFile, "%lu %lu \n", leftNode, rightNode );
             }
             else{
                fprintf(NewFile, "%lu %lu\n", leftNode, rightNode);
//                fprintf(NewFile, "%lu %lu\n", leftPrime, rightPrime);
                break;
             }
         }

         limit=limit+1;
       }
}


//Exercise 4
unsigned long* GenerateDegreeAssign(FILE *file, unsigned long nbNodes){
 
    // reset pointer to file start
    rewind(file);

    // read through the file to create degree assignment
    unsigned long *degrees = malloc(sizeof(unsigned long) * nbNodes);

    for (unsigned long nodeidx=0; nodeidx<nbNodes; nodeidx++) {degrees[nodeidx]=0;};  // initialize

    unsigned long leftNode; unsigned long rightNode;  
    while (fscanf(file,"%lu %lu", &leftNode, &rightNode)==2) {  // update degree count
        ++degrees[leftNode]; ++degrees[rightNode];}

    // create and fill file to store degree assignment
    FILE *degree_assignment = fopen("degree_assignment.txt", "w");

    if (degree_assignment == NULL)  // if problem opening
    {
        printf("Error opening degree assignment file!\n");
        exit(1);
    }
    
    unsigned long Node;
    for (Node=0; Node<nbNodes; Node++) {
        fprintf(degree_assignment, "%lu %lu\n", Node, degrees[Node]);
        };

    // free memory, close / reset pointer to file
    fclose(degree_assignment);
    rewind(file);

    return(degrees);

}

//Exercise 5
void PrintSumDegreePairs(FILE* file, unsigned long nbNodes){

    clock_t begin = clock();

    // open degree assignment file and commit content to a table
    FILE *degree_assignment = fopen("degree_assignment.txt","r") ; // reading file

    // if (degree_assignment == NULL) { // in case misread
    //     return -1;
    // }

    unsigned long *degrees = malloc(sizeof(unsigned long) * nbNodes);
    unsigned long Node; unsigned long Degree;
    while (fscanf(degree_assignment,"%lu %lu", &Node, &Degree)==2) {  // update degree count
        degrees[Node] = Degree;}

    // go through edgelist file and evaluate product of edge-degrees
    unsigned long LeftNode; unsigned long RightNode;
    unsigned long sum_prod_degree = 0;
    while (fscanf(file,"%lu %lu", &LeftNode, &RightNode)==2) {  // update degree count
        unsigned long a = degrees[LeftNode];
        unsigned long b = degrees[RightNode];
        sum_prod_degree += a*b;}

    // stop clock
    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

    printf("The sum-product of edge-degrees is %lu calculated in %f seconds\n", sum_prod_degree, time_spent);

    // free memory, close / reset pointer to file
    fclose(degree_assignment);
    rewind(file);
    free(degrees);
}

//Exercise 6
void GenerateDegreeHist(unsigned long* degrees, unsigned long nbNodes){

    // write degree histogram in file
    FILE *degree_hist = fopen("degree_hist.txt", "w");
    if (degree_hist == NULL)  // if problem opening
    {
        printf("Error opening degree histogram file!\n");
        exit(1);
    }

    // find max degree
    unsigned long maxdegree = 0;
    for (unsigned idx=0; idx < nbNodes; idx++) {
        if (degrees[idx] > maxdegree) {maxdegree = degrees[idx];}
    }

    // degree histogram
    for (unsigned long degree=0; degree <= maxdegree; degree++) {

        unsigned long nbnodesvalued = 0;

        for (unsigned long idx=0; idx < nbNodes; idx++) {
            
            if (degrees[idx]==degree) {++nbnodesvalued;}

        }

    fprintf(degree_hist, "%lu %lu\n", degree, nbnodesvalued);

    }

    // close file
    fclose(degree_hist);
    
}

//Exercise 7.1
unsigned long* ReadAsTable(FILE* file){

    rewind(file);

    // Determine maximum index
    unsigned long Elem1; unsigned long Elem2; unsigned MaxIdx = 0;
    while (fscanf(file,"%lu %lu", &Elem1, &Elem2)==2) {
       ++MaxIdx;}

    rewind(file);

    // Read as Table
    unsigned long *Table = malloc(sizeof(unsigned long) * MaxIdx);
    while (fscanf(file,"%lu %lu", &Elem1, &Elem2)==2) {
        Table[Elem1] = Elem2;}

    rewind(file);

    return (Table);
}

//Exercise 7.2
unsigned long** CreateAdjacencyMatrix(FILE* file, unsigned long nbNodes){

    rewind(file);

    // Create Adjacency Matrix
    unsigned long **AdjacencyMatrix = malloc(sizeof(unsigned long *) * nbNodes);
    for (unsigned long line=0; line<nbNodes; line++){
        AdjacencyMatrix[line] = malloc(sizeof(unsigned long) * nbNodes); 
    }

    // Initialize values to 0
    for (unsigned long i=0; i<nbNodes; i++) {
        for (unsigned long j=0; j<nbNodes; j++){
            AdjacencyMatrix[i][j]=0;
        };
    }

    // Go through edges and activate them with 1 
    unsigned long LeftNode; unsigned long RightNode;
    while (fscanf(file,"%lu %lu", &LeftNode, &RightNode)==2) {
        AdjacencyMatrix[LeftNode][RightNode] = 1;
    }

    rewind(file);
    return(AdjacencyMatrix);
}


void FreeAdjacencyMatrix(unsigned long** AdjacencyMatrix, unsigned long nbNodes){

    for (unsigned long line=0; line<nbNodes; line++){
    free(AdjacencyMatrix[line]);
    }

	free(AdjacencyMatrix);
}

//Exercise 7.3
adjacencyarray* CreateAdjacencyArray(FILE* file, unsigned long directed){


    rewind(file);

    printf("Computing Adjacency Array...\n");
    clock_t begin = clock();


    // Create and Initialize Adjacency Array
    adjacencyarray *AdjacencyArray = malloc(sizeof(adjacencyarray));

    AdjacencyArray->nbnodes = 0;
    AdjacencyArray->maxcapacity = NPOINTS;
    AdjacencyArray->nodes = malloc(sizeof(node) * NPOINTS);
    for (unsigned long i=0; i<NPOINTS; i++) {

        if (directed==0) {  // initialize for undirected graph
            AdjacencyArray->nodes[i].neighbors = malloc(sizeof(unsigned long));
            AdjacencyArray->nodes[i].degree = 0;
            AdjacencyArray->nodes[i].active = 1;
            AdjacencyArray->nodes[i].maxcapacity = 1;
        }

        if (directed==1) {  // initialize for directed graph
            AdjacencyArray->nodes[i].neighborsin = malloc(sizeof(unsigned long));
            AdjacencyArray->nodes[i].neighborsout = malloc(sizeof(unsigned long));
            AdjacencyArray->nodes[i].degreein = 0;
            AdjacencyArray->nodes[i].degreeout = 0;
            AdjacencyArray->nodes[i].degree = 0;
            AdjacencyArray->nodes[i].maxcapacityin = 1;
            AdjacencyArray->nodes[i].maxcapacityout = 1;

        }

    }

    // Go through edges, Fill the Adjacency Array
    unsigned long leftNode; unsigned long rightNode;  
    unsigned long old_node_maxcapacity = AdjacencyArray->maxcapacity;
    unsigned long counter = 0;



    while (fscanf(file,"%lu %lu", &leftNode, &rightNode)==2) {


        // ++counter; 
        // printf("ok for step %lu \n", counter);


        // Nodes: update nb nodes 
        AdjacencyArray->nbnodes = max3(AdjacencyArray->nbnodes, leftNode, rightNode); 
        
        old_node_maxcapacity = AdjacencyArray->maxcapacity;

        // Nodes: reallocate memory for new nodes
        if (AdjacencyArray->nbnodes + 1 >= AdjacencyArray->maxcapacity) {  // because indexing starts at 0
            AdjacencyArray->maxcapacity = AdjacencyArray->nbnodes + 1;
            AdjacencyArray->nodes = realloc(AdjacencyArray->nodes, sizeof(node) * AdjacencyArray->maxcapacity);
        }

        
        // Nodes: initialize new nodes
        for (unsigned long i=old_node_maxcapacity; i<AdjacencyArray->maxcapacity; i++) {
        
            if (directed==0) {  // initialize for undirected graph
                AdjacencyArray->nodes[i].neighbors = malloc(sizeof(unsigned long));
                AdjacencyArray->nodes[i].degree = 0;
                AdjacencyArray->nodes[i].active = 1;
                AdjacencyArray->nodes[i].maxcapacity = 1;
            }

            if (directed==1) {  // initialize for directed graph
                AdjacencyArray->nodes[i].neighborsin = malloc(sizeof(unsigned long));
                AdjacencyArray->nodes[i].neighborsout = malloc(sizeof(unsigned long));
                AdjacencyArray->nodes[i].degreein = 0;
                AdjacencyArray->nodes[i].degreeout = 0;
                AdjacencyArray->nodes[i].degree = 0;
                AdjacencyArray->nodes[i].maxcapacityin = 1;
                AdjacencyArray->nodes[i].maxcapacityout = 1;

            }
        }
        

        old_node_maxcapacity = AdjacencyArray->maxcapacity ; 

        
        // Edges: update degree
        ++AdjacencyArray->nodes[leftNode].degree;
        ++AdjacencyArray->nodes[rightNode].degree;

        if (directed==1) {
            ++AdjacencyArray->nodes[rightNode].degreein;
            ++AdjacencyArray->nodes[leftNode].degreeout;
        }

        // printf("Nb nodes %lu vs Max capacity %lu \n", AdjacencyArray->nbnodes, AdjacencyArray->maxcapacity);
        // printf("Left node %lu degree %lu vs Left node capacity %lu \n", leftNode, AdjacencyArray->nodes[leftNode].degree, AdjacencyArray->nodes[leftNode].maxcapacity);
        // printf("Right node %lu degree %lu vs Right node capacity %lu \n", rightNode, AdjacencyArray->nodes[rightNode].degree, AdjacencyArray->nodes[leftNode].maxcapacity);


        // Edges: reallocate memory for new neighbors

        if (directed==1) {

            if (AdjacencyArray->nodes[leftNode].degreeout + 1 >= AdjacencyArray->nodes[leftNode].maxcapacityout){
                AdjacencyArray->nodes[leftNode].maxcapacityout *= 2;
                AdjacencyArray->nodes[leftNode].neighborsout = realloc(AdjacencyArray->nodes[leftNode].neighborsout, sizeof(unsigned long) * AdjacencyArray->nodes[leftNode].maxcapacityout);
            }
            AdjacencyArray->nodes[leftNode].neighborsout[AdjacencyArray->nodes[leftNode].degreeout - 1] = rightNode; // for indexing

            if (AdjacencyArray->nodes[rightNode].degreein + 1  >= AdjacencyArray->nodes[rightNode].maxcapacityin){
                AdjacencyArray->nodes[rightNode].maxcapacityin *= 2;
                AdjacencyArray->nodes[rightNode].neighborsin = realloc(AdjacencyArray->nodes[rightNode].neighborsin, sizeof(unsigned long) * AdjacencyArray->nodes[rightNode].maxcapacityin);
            }
            AdjacencyArray->nodes[rightNode].neighborsin[AdjacencyArray->nodes[rightNode].degreein - 1] = leftNode; // for indexing

        }


        if (directed==0) {

            if (AdjacencyArray->nodes[leftNode].degree + 1 >= AdjacencyArray->nodes[leftNode].maxcapacity){
                AdjacencyArray->nodes[leftNode].maxcapacity *= 2;
                AdjacencyArray->nodes[leftNode].neighbors = realloc(AdjacencyArray->nodes[leftNode].neighbors, sizeof(unsigned long) * AdjacencyArray->nodes[leftNode].maxcapacity);
            }
            AdjacencyArray->nodes[leftNode].neighbors[AdjacencyArray->nodes[leftNode].degree - 1] = rightNode; // for indexing

            if (AdjacencyArray->nodes[rightNode].degree + 1  >= AdjacencyArray->nodes[rightNode].maxcapacity){
                AdjacencyArray->nodes[rightNode].maxcapacity *= 2;
                AdjacencyArray->nodes[rightNode].neighbors = realloc(AdjacencyArray->nodes[rightNode].neighbors, sizeof(unsigned long) * AdjacencyArray->nodes[rightNode].maxcapacity);
            }
            AdjacencyArray->nodes[rightNode].neighbors[AdjacencyArray->nodes[rightNode].degree - 1] = leftNode; // for indexing

        }


    }


    ++AdjacencyArray->nbnodes;

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Adjacency Array computation time: %f seconds\n", time_spent);
    printf("This graph has %lu nodes\n", AdjacencyArray->nbnodes);

    rewind(file);

    return AdjacencyArray;}


//Exercise 8
fifo* FIFOCreate(unsigned long maxcapacity){
    fifo* FIFO = malloc(sizeof(fifo));
    FIFO->cursor = -1;
    FIFO->maxcapacity = maxcapacity;
    FIFO->elements = malloc(sizeof(unsigned long) * maxcapacity);
    return (FIFO);
}

// Auxiliary Functions for FIFO
void FIFOAdd(fifo* FIFO, unsigned long element){

    // reallocate memory if FIFO at maximum capacity
    if (FIFO->maxcapacity == FIFO->cursor - 1) {
        FIFO->elements = realloc(FIFO->elements, sizeof(unsigned long) * FIFO->maxcapacity * 2);
        FIFO->maxcapacity *= 2;
    }

    // add element, update cursor
    FIFO->elements[FIFO->cursor + 1] = element;
    ++FIFO->cursor;
}

unsigned long FIFOPop(fifo* FIFO){

    // stop if FIFO is empty
    if (FIFO->cursor == -1) {
        printf("FIFO is empty"); return(-1);
    }

    // remove and return last element, update cursor
    --FIFO->cursor;
    return FIFO->elements[FIFO->cursor + 1];
}



//TP 2

//Exercise 1
// power iteration method, requires directed Adjacency Array in input
double* CreatePageRank(adjacencyarray* AdjacencyArray, double alpha, unsigned long t) {  

    printf("Computing Page Rank...");
    clock_t begin = clock();


    // initialization
    double* pagerank = malloc(sizeof(double) * AdjacencyArray->nbnodes);

    double uniformstateproba = (double) 1/AdjacencyArray->nbnodes;

    for (unsigned pageidx = 0; pageidx< AdjacencyArray->nbnodes; pageidx++) {
        pagerank[pageidx] = uniformstateproba;
    }


    for (unsigned long timecount=0; timecount < t; timecount++) {

        // Matrix-Vector product : Transition Array x Current State
        for (unsigned long linenode=0; linenode < AdjacencyArray->nbnodes; linenode++) {

            double pagerankupdate = 0;

            for (unsigned long columnidx=0; columnidx < AdjacencyArray->nodes[linenode].degreein; columnidx++){

                unsigned long columnnode = AdjacencyArray->nodes[linenode].neighborsin[columnidx];

                double T_ij = (double) 1/AdjacencyArray->nodes[columnnode].degreeout;

                pagerankupdate += T_ij * pagerank[columnnode];

            }

            pagerank[linenode] = pagerankupdate;
        
        }

        // Current State reweighted with Uniform State configuration
        for (unsigned long linenode=0; linenode < AdjacencyArray->nbnodes; linenode++) {

            pagerank[linenode] = (1 - alpha) * pagerank[linenode] + alpha * uniformstateproba;

        }

        // Current State gets Normalized
        double norm1 = 0;

        for (unsigned long linenode=0; linenode < AdjacencyArray->nbnodes; linenode++) {
            norm1 += fabs(pagerank[linenode]);
        }

        double corrective = (double) (1 - norm1) / AdjacencyArray->nbnodes;

        for (unsigned long linenode=0; linenode < AdjacencyArray->nbnodes; linenode++) {
            pagerank[linenode] += corrective;
        }

        // printf("page 3 rank is %f \n", pagerank[3]);

    }

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Page Rank computation time: %f seconds\n", time_spent);

    return pagerank;

}


//Excercise 3
void CreatekCore(adjacencyarray* AdjacencyArray, unsigned long write) {

    // k core decomposition
    printf("Computing k-Core... \n");
    clock_t begin = clock();

    // initialization
    unsigned long i = AdjacencyArray->nbnodes;
    unsigned long c = 0;

    // write decomposition in file
    FILE *kcore = fopen("kcore_decomposition.txt", "w") ; // reading file
    if (kcore == NULL)  // if problem opening
    {
        printf("Error opening kcore decomposition file!\n");
        exit(1);
    }

    
    while (i > 0) {

        if (i % 1000 == 0) {
            printf("%lu steps left \n", i);
        }

        // get minimum degree
        unsigned long mindegree = AdjacencyArray->nbnodes;
        for (unsigned idx = 0; idx < AdjacencyArray->nbnodes; idx++){
            if (AdjacencyArray->nodes[idx].degree < mindegree && AdjacencyArray->nodes[idx].active==1) {
                mindegree = AdjacencyArray->nodes[idx].degree;};
        }
        
        // get a minimum degree node
        unsigned long mindegreenode = 0;
        for (unsigned idx = 0; idx < AdjacencyArray->nbnodes; idx++){
            if (AdjacencyArray->nodes[idx].degree == mindegree && AdjacencyArray->nodes[idx].active==1) {
                mindegreenode = idx; break;};
        }  

        // get core value for that node
        c = max3(0, c, mindegree);

        // write in file
        fprintf(kcore, "%lu %lu %lu %lu \n", mindegreenode, i, mindegree, c);
        // prints node idx, its new rank, its degree, its core value

        // update graph: the current node's neighbors
        for (unsigned neighbor = 0; neighbor < AdjacencyArray->nodes[mindegreenode].degree; neighbor++){
            --AdjacencyArray->nodes[neighbor].degree;
        }

        // update graph: the current node
        free(AdjacencyArray->nodes[mindegreenode].neighbors);
        AdjacencyArray->nodes[mindegreenode].active = 0;

        // update i
        --i;
    }

    fclose(kcore);

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("k-Core decomposition computation time: %f seconds\n", time_spent);
    printf("Core value of the graph: %lu\n", c);
}



//TP 3
void GenerateRandomizedGraph(double p, double q){

    // cluster assignment
    unsigned long* clusters = malloc(sizeof(unsigned long) * 400);
    for (unsigned long cluster=0; cluster<4; cluster++){
        for (unsigned long idx= cluster*100; idx< (cluster+1)*100; idx++){
            clusters[idx] = cluster;
        }
    }

    // create file to write in
    FILE* edgelist = fopen("randomizedgraph.txt", "w") ;  
    if (edgelist == NULL) {  // in case misread
        printf("Error opening degree histogram file!\n");
        exit(1);    
    }

    // random number placeholder
    double r;


    // generates randomized edgelist
    for (unsigned long leftNode=0; leftNode<399; leftNode++){
        for (unsigned long rightNode=leftNode+1; rightNode<400; rightNode++){  // be careful not to overlap! otherwise SegFault...

            // case: same cluster
            if (clusters[leftNode]==clusters[rightNode]){
                // generates edge with proba p
                r = (float)rand()/(float)(RAND_MAX/1);
                if (r < p){
                    fprintf(edgelist, "%lu   %lu\n", leftNode%400, rightNode%400);
                }
            }

            // case: different cluster
            else {
                // generates edge with proba q
                r = (float)rand()/(float)(RAND_MAX/1);
                if (r < q){
                    fprintf(edgelist, "%lu   %lu\n", leftNode%400, rightNode%400);
                }
            }

        }
    }


    fclose(edgelist);

}

void shuffle400(unsigned long* array){
    
    // based on https://en.wikipedia.org/wiki/Fisher–Yates_shuffle
    // for integer array
    
    unsigned long j;
    unsigned long r;  // random nb

    for (unsigned i=0; i<=398; i++){
        
        // find a 'good' random pair-node to permute with
        r = rand()%(399-i);
        j = i + r;

        // permute
        unsigned long temp1 = array[i]; 
        unsigned long temp2 = array[j];

        array[i] = temp2;
        array[j] = temp1;

    }

}
