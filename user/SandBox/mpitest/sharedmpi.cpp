#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <iostream>
#include <string>
using namespace std;

class MPI_shared_com
{
public:
    MPI_Comm nodecomm;
    int nodesize, noderank;
    MPI_shared_com(MPI_Comm allcomm, int rank) {
        MPI_Comm_split_type(allcomm, MPI_COMM_TYPE_SHARED, rank, MPI_INFO_NULL, &nodecomm);
        MPI_Comm_size(nodecomm, &nodesize);
        MPI_Comm_rank(nodecomm, &noderank);
    }
};
template<typename Datatype>
class shared_data
{
public:
    Datatype* table;// * localtable;
    MPI_Win wintable;
    int _noderank;
    shared_data(MPI_shared_com shd_c, int tablesize = 1) :_noderank(shd_c.noderank)
    {
        int localtablesize = 0;
        MPI_Aint winsize;
        int windisp;
        if (shd_c.noderank == 0) localtablesize = tablesize;
        MPI_Win_allocate_shared(localtablesize * sizeof(Datatype), sizeof(Datatype),
            MPI_INFO_NULL, shd_c.nodecomm, &table, &wintable);
        if (shd_c.noderank != 0)
        {
            MPI_Win_shared_query(wintable, 0, &winsize, &windisp, &table);
        }
    }
    void lock_all()
    {
        MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, MPI_MODE_NOCHECK, wintable);
        //MPI_Win_lock_all(MPI_MODE_NOCHECK, wintable);
    }
    void unlock_all()
    {
        MPI_Win_unlock(0, wintable);
        //MPI_Win_unlock_all(wintable);
    }
    Datatype& operator()() { return table[0]; }
    Datatype& operator[](int n) { return table[n]; }
};

template<typename Data>
class shared_linklist
{
public:
    struct linknode
    {
        Data data;
        int next;
    };

    struct supportdata
    {
        int list_head, list_tail;
        int mempool_head, tablesize;
    };
    shared_data<supportdata> sdata;
    shared_data<linknode> mempool;
    shared_linklist(MPI_shared_com shd_c, int tablesize) :mempool(shd_c, tablesize), sdata(shd_c)
    {
        sdata().list_head = -1;
        sdata().list_tail = -1;
        sdata().mempool_head = 0;
        sdata().tablesize = tablesize;
    }
    void lock_all()
    {
        sdata.lock_all();
        mempool.lock_all();
    }
    void unlock_all()
    {
        mempool.unlock_all();
        sdata.unlock_all();
    }
    bool push(const Data& input)
    {
        //lock_all();
        if (sdata().mempool_head >= sdata().tablesize)
        {
            return false;
        }
        if (sdata().list_head >= 0)
        {
            linknode temp;
            temp.data = input;
            temp.next = -1;
            mempool.lock_all();

            //MPI_Put( /* data on origin: */   &temp, sizeof(linknode), MPI_BYTE,
                /* data on target: */  // &mempool[sdata().mempool_head], 1, sizeof(linknode), MPI_BYTE,mempool.wintable);
            MPI_Put( /* data on origin: */   &temp, sizeof(linknode), MPI_BYTE, 0,
                /* data on target: */   sdata().mempool_head * sizeof(linknode), 1, sizeof(linknode), MPI_BYTE,
                mempool.wintable);
            temp.data = mempool[sdata().list_tail].data;
            temp.next = sdata().mempool_head;
            MPI_Put( /* data on origin: */   &temp, sizeof(linknode), MPI_BYTE, 0,
                /* data on target: */   sdata().list_tail * sizeof(linknode), 1, sizeof(linknode), MPI_BYTE,
                mempool.wintable);
            //mempool[sdata().list_tail].next = sdata().mempool_head;
            mempool.unlock_all();
            //mempool[sdata().mempool_head].data = input;
            //mempool[sdata().mempool_head].next = -1;
            sdata().list_tail = sdata().mempool_head;
            sdata().mempool_head++;
        }
        else
        {
            sdata().list_head = 0;
            sdata().list_tail = 0;
            sdata().mempool_head = 1;
            mempool[0].data = input;
            mempool[0].next = -1;
        }
        //unlock_all();
    }
    //void exe() { if (_noderank == 0)table[0].data = 100; }
    linknode& operator[](int n) { return mempool[n]; }


};

int main(void)
{
    int i, flag;

    int nodesize, noderank;
    int size, rank, irank;
    int tablesize, localtablesize;
    int* table, * localtable;
    int* model;

    MPI_Comm allcomm, nodecomm;

    char verstring[MPI_MAX_LIBRARY_VERSION_STRING];
    char nodename[MPI_MAX_PROCESSOR_NAME];

    MPI_Aint winsize;
    int windisp;
    int* winptr;

    int version, subversion, verstringlen, nodestringlen;

    allcomm = MPI_COMM_WORLD;

    MPI_Win wintable;



    MPI_Init(NULL, NULL);

    MPI_Comm_size(allcomm, &size);
    MPI_Comm_rank(allcomm, &rank);

    MPI_Get_processor_name(nodename, &nodestringlen);

    MPI_Get_version(&version, &subversion);
    MPI_Get_library_version(verstring, &verstringlen);

    if (rank == 0)
    {
        printf("Version %d, subversion %d\n", version, subversion);
        printf("Library <%s>\n", verstring);
    }

    // Create node-local communicator

    MPI_Comm_split_type(allcomm, MPI_COMM_TYPE_SHARED, rank,
        MPI_INFO_NULL, &nodecomm);

    MPI_Comm_size(nodecomm, &nodesize);
    MPI_Comm_rank(nodecomm, &noderank);

    // Only rank 0 on a node actually allocates memory
    tablesize = nodesize;
    localtablesize = 0;

    if (noderank == 0) localtablesize = tablesize;

    // debug info

    printf("Rank %d of %d, rank %d of %d in node <%s>, localtablesize %d\n",
        rank, size, noderank, nodesize, nodename, localtablesize);


    MPI_Win_allocate_shared(localtablesize * sizeof(int), sizeof(int),
        MPI_INFO_NULL, nodecomm, &localtable, &wintable);

    MPI_Win_get_attr(wintable, MPI_WIN_MODEL, &model, &flag);

    if (1 != flag)
    {
        printf("Attribute MPI_WIN_MODEL not defined\n");
    }
    else
    {
        if (MPI_WIN_UNIFIED == *model)
        {
            if (rank == 0) printf("Memory model is MPI_WIN_UNIFIED\n");
        }
        else
        {
            if (rank == 0) printf("Memory model is *not* MPI_WIN_UNIFIED\n");

            MPI_Finalize();
            return 1;
        }
    }

    // need to get local pointer valid for table on rank 0

    table = localtable;

    if (noderank != 0)
    {
        MPI_Win_shared_query(wintable, 0, &winsize, &windisp, &table);
    }

    // All table pointers should now point to copy on noderank 0

    // Initialise table on rank 0 with appropriate synchronisation

    MPI_Win_fence(0, wintable);

    //table[noderank] =  noderank;

    if (noderank == 0)
    {
        for (i = 0; i < tablesize; i++)
        {
            table[i] = rank * tablesize + i;
        }
    }


    MPI_Win_fence(0, wintable);

    // Check we did it right

    printf("rank %d, noderank %d, table[%d] = %d\n",
        rank, noderank, 0, table[0]);
    /*
for (i = 0; i < tablesize; i++)
{
    printf("rank %d, noderank %d, table[%d] = %d\n",
        rank, noderank, i, table[i]);
}*/
    MPI_shared_com shd_comm(MPI_COMM_WORLD, rank);
    shared_linklist<int> lltest(shd_comm, 100000);
    //shared_data<int> lltest(shd_comm, 100000);
    //cout << "haha" << endl;
    //MPI_Win_fence(0, wintable);
/*    if (noderank == 0)
    {
        lltest.push(100);
        lltest.push(10);
        lltest.push(10);
        lltest.push(105);
        lltest.push(1);
    }*/
    //MPI_Win_fence(0, wintable);
    for (int i = 0;i < 3;i++)
    {
        lltest.push(rank * 10 + i);
    }
    //MPI_Barrier(MPI_COMM_WORLD);
    //MPI_Win_fence(0, wintable);
    //printf("rank %d, noderank %d, table[%d] = %d\n", rank, noderank, 0, lltest[0].data);

    /*if (noderank == 1)
    {
        printf("rank %d, noderank %d\n", rank, noderank);
        for (int i = 0;i != -1;i = lltest[i].next)
        {
            //printf("i=%d,data= %d,next=%d\n", i, lltest[i].data, lltest[i].next);
            printf(" %d", lltest[i].data);
        }
        printf("\n");
    }*/
    MPI_Win_fence(0, wintable);

    /*printf("rank %d, noderank %d\n", rank, noderank);
    for (int i = 0;i != -1;i = lltest[i].next)
    {
        //printf("i=%d,data= %d,next=%d\n", i, lltest[i].data, lltest[i].next);
        printf(" %d", lltest[i].data);
    }*/
    string str = string("rank") + to_string(rank) + string("\n");
    for (int i = 0;i < 18;i++)
    {
        //printf("i=%d,data= %d,next=%d\n", i, lltest[i].data, lltest[i].next);
        str += string(" ") + to_string(lltest[i].data) + string(",") + to_string(lltest[i].next);

    }
    cout << str << endl;






    MPI_Finalize();
}
