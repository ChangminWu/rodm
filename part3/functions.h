
// structure for FIFO

typedef struct {
    unsigned long cursor ;
    unsigned long maxcapacity;
    unsigned long *elements ;
} fifo;



// structures for Adjacency Array

typedef struct {
    unsigned long* neighbors;
    unsigned long* neighborsin;
    unsigned long* neighborsout;
    unsigned long degree;
    unsigned long degreein;
    unsigned long degreeout;
    unsigned long active;
    unsigned long maxcapacity;
    unsigned long maxcapacityin;
    unsigned long maxcapacityout;
} node;

typedef struct {
    unsigned long nbnodes;
    unsigned long maxcapacity;
    node *nodes;
} adjacencyarray;


inline unsigned long max3(unsigned long a, unsigned long b, unsigned long c);

unsigned long ReadNbNodes(FILE *file);

void CleaningData(FILE *file);

unsigned long* GenerateDegreeAssign(FILE *file, unsigned long nbNodes);

void GenerateDegreeHist(unsigned long* degrees, unsigned long nbNodes);

unsigned long* ReadAsTable(FILE* file);

void PrintSumDegreePairs(FILE* file, unsigned long nbNodes);

unsigned long** CreateAdjacencyMatrix(FILE* file, unsigned long nbNodes);

adjacencyarray* CreateAdjacencyArray(FILE* file, unsigned long directed);

double* CreatePageRank(adjacencyarray* AdjacencyArray, double alpha, unsigned long t);

fifo* FIFOCreate(unsigned long maxcapacity);

void FIFOAdd(fifo* FIFO, unsigned long element);

unsigned long FIFOPop(fifo* FIFO);

void CreatekCore(adjacencyarray* AdjacencyArray, unsigned long write);

void GenerateRandomizedGraph(double p, double q);

void shuffle400(unsigned long* array);


