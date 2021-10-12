#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <iostream>
#include <atomic>
using namespace std;
class mpi_manager
{
public:
    MPI_Comm comm;
    int size, rank;
};
class mpi_node_manager : public mpi_manager
{
public:
    mpi_node_manager(MPI_Comm allcomm) {
        int rank_t = 0;;
        MPI_Comm_rank(allcomm, &rank);
        MPI_Comm_split_type(allcomm, MPI_COMM_TYPE_SHARED, rank_t, MPI_INFO_NULL, &comm);
        MPI_Comm_size(comm, &size);
        MPI_Comm_rank(comm, &rank);
    }
};
class mpi_sync
{
public:
    int sync_, sync_r_;
    //MPI_Comm comm;
    mpi_manager& manager;
    mpi_sync(mpi_manager& manager_in) : sync_(0), sync_r_(0), manager(manager_in) {}
    void sync() {
        MPI_Reduce(&sync_, &sync_r_, 1, MPI_INT, MPI_MAX, 0, manager.comm);
        MPI_Barrier(MPI_COMM_WORLD);
    }
};
class my_mutex {
    atomic_flag& flag;
public:
    my_mutex(atomic_flag* ptr, mpi_sync& sync_in) :flag(*ptr) { if (sync_in.manager.rank == 0) flag.clear(); sync_in.sync(); }
    void lock() {
        while (flag.test_and_set(std::memory_order_acquire));
    }
    void unlock() {
        flag.clear(memory_order_release);
    }
};
template<typename Datatype>
class shared_data
{
public:
    Datatype* table;// * localtable;
    MPI_Win wintable;
    mpi_manager& manger;
    int size;
    shared_data(mpi_manager& manger_in, int size_in = 1) : manger(manger_in), size(size_in)
    {
        int localsize = 0;
        //int* model;
        //int flag;
        MPI_Aint winsize;
        int windisp;
        if (manger.rank == 0) localsize = size;
        MPI_Win_allocate_shared(localsize * sizeof(Datatype), sizeof(Datatype),
            MPI_INFO_NULL, manger.comm, &table, &wintable);
        //table = localtable;
        //MPI_Win_get_attr(wintable, MPI_WIN_MODEL, &model, &flag);
        if (manger.rank != 0)
        {
            MPI_Win_shared_query(wintable, 0, &winsize, &windisp, &table);
        }
    }
    Datatype& operator()() { return *table; }
    Datatype& operator[](int n) { return table[n]; }
    Datatype* ptr() { return table; }
};
class mpi_mutex {
    shared_data<atomic_flag> flag_shared;

    my_mutex mutex_;
public:
    mpi_mutex(mpi_manager& manger_in, mpi_sync& sync_in) :flag_shared(manger_in), mutex_(flag_shared.ptr(), sync_in) {}
    //my_mutex(void* ptr, mpi_sync& sync_in) :flag(*(atomic_flag*)ptr) { if (sync_in.manager.rank == 0) flag.clear(); sync_in.sync(); }
    void lock() {
        mutex_.lock();
    }
    void unlock() {
        mutex_.unlock();
    }
};
int main(void)
{
    int i, flag;

    //int nodesize, noderank;
    int size, rank;//, irank;
    int tablesize, localtablesize;
    int* table, * localtable;
    std::atomic<int>* table_a;
    int* model;

    MPI_Comm allcomm, nodecomm;

    //char verstring[MPI_MAX_LIBRARY_VERSION_STRING];
    //char nodename[MPI_MAX_PROCESSOR_NAME];

    MPI_Aint winsize;
    int windisp;
    //int* winptr;

    //int version, subversion, verstringlen, nodestringlen;

    allcomm = MPI_COMM_WORLD;

    MPI_Win wintable;



    MPI_Init(NULL, NULL);

    MPI_Comm_size(allcomm, &size);
    MPI_Comm_rank(allcomm, &rank);


    // Create node-local communicator

    MPI_Comm_split_type(allcomm, MPI_COMM_TYPE_SHARED, rank,
        MPI_INFO_NULL, &nodecomm);
    mpi_node_manager node_manager(allcomm);

    //MPI_Comm_size(nodecomm, &nodesize);
    //MPI_Comm_rank(nodecomm, &noderank);

    // Only rank 0 on a node actually allocates memory
    tablesize = 2;
    localtablesize = 0;

    if (node_manager.rank == 0) localtablesize = tablesize;

    // debug info

    //printf("Rank %d of %d, rank %d of %d in node <%s>, localtablesize %d\n", rank, size, noderank, nodesize, nodename, localtablesize);


    MPI_Win_allocate_shared(localtablesize * sizeof(int), sizeof(int),
        MPI_INFO_NULL, node_manager.comm, &localtable, &wintable);

    MPI_Win_get_attr(wintable, MPI_WIN_MODEL, &model, &flag);
    table = localtable;
    // need to get local pointer valid for table on rank 0

    if (node_manager.rank != 0)
    {
        MPI_Win_shared_query(wintable, 0, &winsize, &windisp, &table);
    }

    // All table pointers should now point to copy on noderank 0

    // Initialise table on rank 0 with appropriate synchronisation
    table_a = (std::atomic<int>*)table;
    mpi_sync sync0(node_manager);
    //my_mutex mutexmy((atomic_flag*)&table[1], sync0);
    mpi_mutex mutexmy(node_manager, sync0);
    MPI_Win_fence(0, wintable);

    //table[noderank] =  noderank;
    //mutexmy.lock();
    table_a[0] = 0;
    sync0.sync();
    for (int i = 0; i < 1000; i++)
    {
        mutexmy.lock();
        table[0]= table[0]+2;
        mutexmy.unlock();
        //table_a[0]++;
    }

    sync0.sync();
    cout << *table_a << endl;
    
    
    //mutexmy.unlock();

    
    //MPI_Reduce(&sync, &sync_r, 1, MPI_INT, MPI_MAX, 0, nodecomm);
    //MPI_Barrier(MPI_COMM_WORLD);
    //cout << sizeof(atomic_flag) << endl;
    //if (rank == 0)
    //    cout << "-----------------" << endl;
    //sync0.sync();
    //MPI_Reduce(&sync, &sync_r, 1, MPI_INT, MPI_MAX, 0, nodecomm);
    //MPI_Barrier(MPI_COMM_WORLD);
    //cout << "before" << endl;
    //MPI_Reduce(&sync, &sync_r, 1, MPI_INT, MPI_MAX, 0, nodecomm);
    //MPI_Barrier(MPI_COMM_WORLD);
    //sync0.sync();
    //if (rank == 0)
    //    cout << "-----------------" << endl;
    //sync0.sync();
    //MPI_Reduce(&sync, &sync_r, 1, MPI_INT, MPI_MAX, 0, nodecomm);
    //MPI_Barrier(MPI_COMM_WORLD);
    //cout << "after" << endl;



    MPI_Win_fence(0, wintable);

    // Check we did it right
/*
    if (noderank == 0)
        printf("rank %d, noderank %d, table[%d] = %d\n",
            rank, noderank, 0, table[0]);
            */

            /*
        for (i = 0; i < tablesize; i++)
        {
            printf("rank %d, noderank %d, table[%d] = %d\n",
                rank, noderank, i, table[i]);
        }*/



    MPI_Finalize();
}
