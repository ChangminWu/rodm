#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>  // for abs value fct
#include <time.h>//to estimate the runing time
#include "functions.h"


#define NLINKS 10 //maximum number of edges for memory allocation, will increase if needed
#define NPOINTS 200 //maximum number of edges for memory allocation, will increase if needed

/*
To compile:
gcc p1.c functions.c -O9 -o out

To execute:
./out wiki_edgelist.txt

Prior structure:
email, amazon, livejournal, orkut, friendster txt files with no heading
p1 repository
*/

// Implement free?


int main(int argc,char** argv){

    // p1ex5
    time_t t1,t2;
    t1=time(NULL);
    // email (substitutable for another undirected graph)
    printf("For lr... \n");
    FILE* file = fopen("com-lj.ungraph.txt","r") ;  // reading file, substitute filename here to work with other graphs
    if (file == NULL) {  // in case misread
        return -1;
    }  

    CleaningData(file);

    t2=time(NULL);

    printf("- Overall time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));

    //unsigned long nbNodes = ReadNbNodes(file);

    //unsigned long* degrees = GenerateDegreeAssign(file, nbNodes);

    //GenerateDegreeHist(degrees, nbNodes);

    //PrintSumDegreePairs(file, nbNodes);

    //adjacencyarray* AdjacencyArray = CreateAdjacencyArray(file, 0);



    // // Implementing BFSColor

    // unsigned long *Coloring = malloc(sizeof(unsigned long) * nbNodes); 
    // for (unsigned long idx = 0; idx < nbNodes; idx++) {Coloring[idx]=0;} // Color 0 = not marked
    // fifo *FIFO = FIFOCreate(nbNodes);

    // unsigned long StartNode = 0; unsigned long CurrentColor = 0;

    // for (unsigned long idx = 0; idx < nbNodes; idx++) {
        
    //     if (Coloring[idx]==0) {

    //         StartNode = idx;
    //         ++CurrentColor; 
    //         FIFOAdd(FIFO, StartNode);
    //         Coloring[StartNode] = CurrentColor;

    //         while (FIFO->cursor != -1){  // to indicate FIFO not empty 

    //             unsigned long CurrentNode = FIFOPop(FIFO);

    //             // printf("CurrentNode: %lu \n", CurrentNode);
    //             /eadNbNod/ printf("FIFO Cursor: %lu ... FIFO MaxCap: %lu \n", FIFO->cursor, FIFO->maxcapacity);

    //             for (unsigned long idx = 0; idx < AdjacencyArray->nodes[CurrentNode].degree; idx++){
                    
    //                 unsigned long neighbor = AdjacencyArray->nodes[CurrentNode].neighbors[idx];

    //                 // if neighbor not marked, add to FIFO and mark 
    //                 if (Coloring[neighbor]==0) {
    //                     FIFOAdd(FIFO, neighbor);
    //                     Coloring[neighbor] = CurrentColor;
    //                 }

    //             }
    //         }

    //     }
    
    // }

    // FILE *ConnectedComponents = fopen("connected_components.txt", "w") ; // reading file
    // if (ConnectedComponents == NULL)  // if problem opening
    // {
    //     printf("Error opening connected components file!\n");
    //     exit(1);
    // }
    
    // for (unsigned long idx=0; idx < nbNodes; idx++) {
    //     fprintf(ConnectedComponents, "%lu %lu\n", idx, Coloring[idx]);
    //     };

    // fclose(ConnectedComponents);






    return 0;

}
