#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <atomic>
#include <chrono>
#include <iostream>
#include <random>
using namespace std;
class mpi_manager {
   public:
    MPI_Comm comm;
    int size, rank;
};
class mpi_node_manager : public mpi_manager {
   public:
    mpi_node_manager(MPI_Comm allcomm) {
        int rank_t = 0;
        ;
        MPI_Comm_rank(allcomm, &rank);
        MPI_Comm_split_type(allcomm, MPI_COMM_TYPE_SHARED, rank_t, MPI_INFO_NULL, &comm);
        MPI_Comm_size(comm, &size);
        MPI_Comm_rank(comm, &rank);
    }
};
class mpi_sync {
    int sync_, sync_r_;

   public:
    mpi_manager& manager;
    mpi_sync(mpi_manager& manager_in)
        : sync_(0), sync_r_(0), manager(manager_in) {}
    void sync() {
        MPI_Reduce(&sync_, &sync_r_, 1, MPI_INT, MPI_MAX, 0, manager.comm);
        MPI_Barrier(MPI_COMM_WORLD);
    }
};
class my_mutex {
    atomic_flag& flag;

   public:
    my_mutex(atomic_flag* ptr, mpi_sync& sync_in)
        : flag(*ptr) {
        if (sync_in.manager.rank == 0)
            flag.clear();
        sync_in.sync();
    }
    void lock() {
        while (flag.test_and_set(std::memory_order_acquire))
            ;
    }
    void unlock() { flag.clear(memory_order_release); }
};

class my_room {
    atomic_flag& flag;
    atomic<int>& p_in_room;

   public:
    my_room(atomic_flag* ptr, atomic<int>* ptr2, mpi_sync& sync_in)
        : flag(*ptr), p_in_room(*ptr2) {
        if (sync_in.manager.rank == 0) {
            flag.clear();
            p_in_room = 0;
        }
        sync_in.sync();
    }
    void lock() {
        while (flag.test_and_set(std::memory_order_acquire))
            ;
    }
    void unlock() { flag.clear(memory_order_release); }
    void enter() {
        while (flag.test_and_set(std::memory_order_acquire))
            ;
        p_in_room++;
        flag.clear(memory_order_release);
    }
    void exit() { p_in_room--; }
    void give_back() {
        flag.clear(memory_order_release);
        p_in_room--;
    }
    void require_empty_room() {
        while (flag.test_and_set(std::memory_order_acquire))
            ;
        while (p_in_room != 0)
            ;
        p_in_room++;
    }
};
template <typename Datatype>
class shared_data {
    Datatype* table;
    MPI_Win wintable;
    mpi_manager& manager;
    int size;

   public:
    shared_data(mpi_manager& manager_in, int size_in = 1)
        : manager(manager_in), size(size_in) {
        int localsize = 0, windisp;
        MPI_Aint winsize;
        if (manager.rank == 0)
            localsize = size;
        MPI_Win_allocate_shared(localsize * sizeof(Datatype), sizeof(Datatype), MPI_INFO_NULL, manager.comm, &table, &wintable);
        if (manager.rank != 0) {
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
    mpi_mutex(mpi_manager& manager_in, mpi_sync& sync_in)
        : flag_shared(manager_in), mutex_(flag_shared.ptr(), sync_in) {}
    void lock() { mutex_.lock(); }
    void unlock() { mutex_.unlock(); }
};

class mpi_room {
    shared_data<atomic_flag> flag_shared;
    shared_data<atomic<int>> a_int_shared;

   public:
    my_room room_;
    mpi_room(mpi_manager& manager_in, mpi_sync& sync_in)
        : flag_shared(manager_in), a_int_shared(manager_in), room_(flag_shared.ptr(), a_int_shared.ptr(), sync_in) {}
    void lock() { room_.lock(); }
    void unlock() { room_.unlock(); }
    void enter() { room_.enter(); }
    void exit() { room_.exit(); }
    void give_back() { room_.give_back(); }
    void require_empty_room() { room_.require_empty_room(); }
};
template <typename Data>
class shared_linklist {
   public:
    struct linknode {
        Data data;
        atomic<int> next;
    };

    struct supportdata {
        atomic<int> list_head, list_tail;
        atomic<int> size, capacity;
        atomic<int> spare_head;
        // atomic<int> mempool_head;
    };

    shared_data<linknode> mempool;
    shared_data<supportdata> sdata;
    shared_data<int> memo_spare;
    mpi_mutex mutex_;
    shared_linklist(mpi_manager& manager_in, int capacity_in, mpi_sync& sync_in)
        : mempool(manager_in, capacity_in), sdata(manager_in), memo_spare(manager_in, capacity_in), mutex_(manager_in, sync_in) {
        if (manager_in.rank == 0) {
            sdata().list_head = -1;
            sdata().list_tail = -1;
            sdata().capacity = capacity_in;
            sdata().spare_head = 0;
            sdata().size = 0;
            for (int i = 0; i < capacity_in; i++) {
                memo_spare[i] = i;
                mempool[i].next = -2;
            }
        }
        sync_in.sync();
    }
    bool insert(const Data& input, int position = -1) {
        mutex_.lock();
        if (sdata().size >= sdata().capacity ||
            (position >= 0 && mempool[position].next == -2)) {
            mutex_.unlock();
            return false;
        }

        atomic<int> loc;
        loc = memo_spare[sdata().spare_head];
        sdata().spare_head++;
        mempool[loc].data = input;
        if (position >= 0) {
            mempool[loc].next.store(mempool[position].next);
        } else if (position == -1) {
            mempool[loc].next.store(sdata().list_head);
        }

        if (sdata().size > 0) {
            if (position == -1) {
                sdata().list_head.store(loc);
            } else {
                mempool[position].next.store(loc);
            }
            if (mempool[loc].next == -1)
                sdata().list_tail.store(loc);
            sdata().size++;
            // sdata().mempool_head++;
        } else {
            mempool[loc].next.store(-1);
            sdata().list_head.store(loc);
            sdata().list_tail.store(loc);
            sdata().size++;
            // sdata().mempool_head = 1;
        }
        mutex_.unlock();
        return true;
    }

    bool remove(int position = -1) {
        mutex_.lock();
        if (sdata().size <= 0 ||
            (position >= 0 &&
             (mempool[position].next == -1 || mempool[position].next == -2))) {
            mutex_.unlock();
            return false;
        }
        atomic<int> loc;
        if (position == -1) {
            loc.store(sdata().list_head);
        } else {
            loc.store(mempool[position].next);
        }

        sdata().size--;
        if (position == -1) {
            sdata().list_head.store(mempool[sdata().list_head].next);
        } else {
            mempool[position].next.store(mempool[mempool[position].next].next);
        }

        sdata().spare_head--;
        memo_spare[sdata().spare_head] = loc;
        mempool[loc].next = -2;

        mutex_.unlock();
        return true;
    }
    bool push(const Data& input) {
        mutex_.lock();
        if (sdata().size >= sdata().capacity) {
            mutex_.unlock();
            return false;
        }

        atomic<int> loc;
        loc = memo_spare[sdata().spare_head];
        sdata().spare_head++;
        mempool[loc].data = input;
        mempool[loc].next = -1;

        if (sdata().size > 0) {
            // sdata().mempool_head++;
            mempool[sdata().list_tail].next.store(loc);
            sdata().list_tail.store(loc);
            sdata().size++;
            // sdata().mempool_head++;
        } else {
            sdata().list_head.store(loc);
            sdata().list_tail.store(loc);
            sdata().size++;
            // sdata().mempool_head = 1;
        }
        mutex_.unlock();
        return true;
    }
    // void exe() { if (_noderank == 0)table[0].data = 100; }
    void lock() { mutex_.lock(); }
    void unlock() { mutex_.unlock(); }
    int head() { return sdata().list_head; }
    int size() { return sdata().size; }
    linknode& operator[](int n) { return mempool[n]; }
};
int main(void) {
    std::atomic<int>* table_a;

    MPI_Comm allcomm = MPI_COMM_WORLD;

    MPI_Init(NULL, NULL);

    mpi_node_manager node_manager(allcomm);
    shared_data<atomic<int>> int_share(node_manager);
    table_a = (std::atomic<int>*)(int_share.ptr());
    mpi_sync sync0(node_manager);
    mpi_mutex mutexmy(node_manager, sync0);
    mpi_room roommy(node_manager, sync0);

    // int_share() = 0;
    /*for (int i = 0; i < 1000000; i++) {
        mutexmy.lock();
        int_share()++;
        //(*table_a)++;
        mutexmy.unlock();
    }*/
    //sync0.sync();
    //cout << int_share() << endl;

    shared_linklist<int> lltest(node_manager, 10000, sync0);
    //bool ordered = false;
    sync0.sync();
    for (int i = 0; i < (node_manager.rank+1) * 10 + 10; i++) {
        lltest.push(node_manager.rank * 100 + i);
    }
    sync0.sync();
    if (node_manager.rank == 0) {
        //cout << "aaa!" << endl;
        // mutexmy.lock();

        string str = string(":");
        int j = 0;
        for (int i = lltest.head(); i != -1; i = lltest[i].next) {
            // printf("i=%d,data= %d,next=%d\n", i, lltest[i].data, lltest[i].next);
            //str += string(" (") + to_string(lltest[i].data) + string(",") + to_string(lltest[i].next) + string(")");
            str += string(" ") + to_string(lltest[i].data);

            j++;
            //cout << lltest[i].data << endl;
        }
        cout << "rank" << node_manager.rank << ":" << j << "," << lltest.size() << str
             << endl;
    }
    int_share() = 1;
    sync0.sync();
    /*
    for (int i = 0; i < 1000; i++) {
        for (int j = 0; j < 1000; j++) {
            if ((j + node_manager.rank) % 100 == 23) {
                roommy.require_empty_room();
                for (int k = 1; k < 30; k++)
                    int_share() = k * int_share();
                int_share() = 1;
                roommy.give_back();
            } else {
                roommy.enter();
                int x = int_share();
                if (x != 1) {
                    cout << "error:int=" << x << endl;
                }
                roommy.exit();
            }
        }
    }
*/

    int save, found_i, insert_pos;
    bool found_disordered = false, success = false, ordered = false, insert = false, found_pos = false;
    while (!ordered) {
        if (found_disordered) {
            roommy.require_empty_room();
            int inext = lltest[found_i].next;
            if (inext >= 0 && lltest[found_i].data > lltest[inext].data) {
                save = lltest[inext].data;
                found_pos = lltest.remove(found_i);
            }
            found_disordered = false;
            //cout << "here:" << success << endl;
            roommy.give_back();
        } else if (found_pos) {
            roommy.enter();
            int i = lltest.head();
            if (save <= lltest[i].data) {
                insert_pos = -1;
                found_pos = false;
            } else {
                for (; i != -1; i = lltest[i].next) {
                    int inext = lltest[i].next;
                    if (inext >=0 && save <= lltest[inext].data) {
                        insert_pos = i;
                        found_pos = false;
                        break;
                    }
                }
                if (found_pos == true) {
                    insert_pos = -2;
                    found_pos = false;
                }
            }
            insert = true;
            roommy.exit();
        } else if (insert) {
            roommy.require_empty_room();
            if (insert_pos == -2) {
                lltest.push(save);
                insert = false;
            } else if (insert_pos == -1) {
                lltest.insert(save);
                insert = false;
            } else {
                int inext = lltest[insert_pos].next;
                if (inext >= 0 && save <= lltest[inext].data) {
                    found_pos = !lltest.insert(save, insert_pos);
                    insert = false;
                } else {
                    found_pos = true;
                    insert = false;
                }
            }

            roommy.give_back();
        } else {
            roommy.enter();
            int i;
            for (i = lltest.head(); i != -1; i = lltest[i].next) {
                if (lltest[i].next >= 0) {
                    int inext = lltest[i].next;
                    if (lltest[i].data > lltest[inext].data) {
                        found_i = i;
                        found_disordered = true;
                        break;
                    }
                }
            }
            if (i == -1)
                ordered = true;
            roommy.exit();
        }
    }

    sync0.sync();
    if (node_manager.rank == 0) {
        //cout << "aaa!" << endl;
        // mutexmy.lock();

        string str = string(":");
        int j = 0;
        for (int i = lltest.head(); i != -1; i = lltest[i].next) {
            // printf("i=%d,data= %d,next=%d\n", i, lltest[i].data, lltest[i].next);
            //str += string(" (") + to_string(lltest[i].data) + string(",") + to_string(lltest[i].next) + string(")");
            str += string(" ") + to_string(lltest[i].data);

            j++;
            //cout << lltest[i].data << endl;
        }
        cout << "rank" << node_manager.rank << ":" << j << "," << lltest.size() << str
             << endl;
    }
    sync0.sync();
    //cout << "rank" << node_manager.rank << ":" << save << endl;

    // cout << str << endl;
    // std::cout.flush();
    // mutexmy.unlock();
    MPI_Finalize();
}
