#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <iostream>
using namespace std;

template <typename Datatype>
class shared_data
{
public:
    Datatype *table; // * localtable;
    MPI_Win wintable;
    int _noderank;
    shared_data(MPI_Comm nodecomm, int noderank, int tablesize = 1) : _noderank(noderank)
    {
        int localtablesize = 0;
        MPI_Aint winsize;
        int windisp;
        if (noderank == 0)
            localtablesize = tablesize;
        MPI_Win_allocate_shared(localtablesize * sizeof(Datatype), sizeof(Datatype),
                                MPI_INFO_NULL, nodecomm, &table, &wintable);
        if (noderank != 0)
        {
            MPI_Win_shared_query(wintable, 0, &winsize, &windisp, &table);
        }
    }
    Datatype &operator()() { return &table; }
    Datatype &operator[](int n) { return table[n]; }
};

int main(void)
{

    int nodesize, noderank;
    int size, rank;

    MPI_Comm allcomm, nodecomm;

    allcomm = MPI_COMM_WORLD;

    MPI_Init(NULL, NULL);

    MPI_Comm_size(allcomm, &size);
    MPI_Comm_rank(allcomm, &rank);

    // Create node-local communicator

    MPI_Comm_split_type(allcomm, MPI_COMM_TYPE_SHARED, rank,
                        MPI_INFO_NULL, &nodecomm);

    MPI_Comm_size(nodecomm, &nodesize);
    MPI_Comm_rank(nodecomm, &noderank);

    shared_data<int> table(nodecomm, noderank, nodesize);
    MPI_Barrier(nodecomm);
    if (noderank == 0)
    {
        for (int i = 0; i < nodesize; i++)
        {
            table[i] = rank * nodesize + i;
        }
    }
    MPI_Barrier(nodecomm);
    printf("rank %d, noderank %d, table[%d] = %d\n", rank, noderank, noderank, table[noderank]);

    MPI_Finalize();
}
