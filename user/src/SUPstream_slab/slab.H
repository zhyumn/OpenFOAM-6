/* ref: https://github.com/bbu/userland-slab-allocator */

#ifndef __GNUC__
#error Can be compiled only with GCC.
#endif

#pragma once

#include <stdint.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <stddef.h>
#include <sys/mman.h>
#include <math.h>
#include <assert.h>
#include "create_memory.H"
#include <time.h>

#define HY_DEBUG

namespace Foam
{
    struct slab_header;
    typedef SharedPointer<slab_header, 0> SharedPointer_SLAB;

    struct slab_header
    {
        SharedPointer_SLAB prev, next;
        uint64_t slots;
        uintptr_t refcount;
        SharedPointer_SLAB page;
        uint8_t data[];
    };

    class slab_chain
    {
        size_t itemsize, itemcount;
        size_t slabsize, pages_per_alloc;
        uint64_t initial_slotmask, empty_slotmask;
        uintptr_t alignment_mask;
        SharedPointer_SLAB partial, empty, full;

    public:
        //MemPool* pmempool;
        void init(size_t);
        SHAREDPOINTER alloc(MemPool &);
        void free(SHAREDPOINTER, MemPool &);
        void traverse(void (*)(const void *));
        void destroy(MemPool &);

#ifndef NDEBUG
        int is_valid();
#endif
        void report();
        void report(int &a, int &b, int &c);
        void stats(FILE *const out);
        void dump(FILE *const out);
        void props(FILE *const out);
    };

    class Slab
    {
    public:
        MemPool mempool;
        //slab_chain* s;
        Foam::SUPstream::mpi_manager &manager_;
        size_t nslab_chain;
        Foam::SUPstream::shared_data<slab_chain> s;
        Slab(Foam::SUPstream::mpi_manager &manager_in, size_t maximum_size, size_t memsize) : mempool(manager_in, memsize),
                                                                                              manager_(manager_in),
                                                                                              nslab_chain(ceil(log2(maximum_size))),
                                                                                              s(manager_in, nslab_chain + 1)
        {
            MemStart<0>::start = (size_t)mempool.memory;
            //nslab_chain = ceil(log2(maximum_size));
            //s = new slab_chain[nslab_chain + 1];
            if (manager_in.rank == 0)
            {
                size_t isize = 1;
                for (size_t i = 0; i <= nslab_chain; i++)
                {
                    //s[i].pmempool = &mempool;
                    //std::cout << "isize= !!!!!!!!!!!!!!!!!!!!!!!!!!!"<<isize<<std::endl;
                    s[i].init(isize);
                    isize *= 2;
                }
            }
        }
        void d_report()
        {
            mempool.report();
            int to = 0, us = 0, fr = 0;
            for (size_t i = 0; i <= nslab_chain; i++)
            {
                s[i].report();
                s[i].report(to, us, fr);
            }
            std::cout << "Slab" << std::endl;
            std::cout << "Total:" << to << " Used:" << us << " Free:" << fr << std::endl;
        }
        void report()
        {
            mempool.report();
            int to = 0, us = 0, fr = 0;
            for (size_t i = 0; i <= nslab_chain; i++)
            {
                s[i].report(to, us, fr);
            }
            std::cout << "Slab" << std::endl;
            std::cout << "Total:" << to << " Used:" << us << " Free:" << fr << std::endl;
        }
        SHAREDPOINTER alloc(size_t size)
        {
            //std::cout << "alloc size= !!!!!!!!!!!!!!!!!!!!!!!!!!!"<<size<<std::endl;
            size_t index = ceil(log2(size));
            //std::cout << "alloc index= !!!!!!!!!!!!!!!!!!!!!!!!!!!"<<index<<std::endl;
            //std::cout << "alloc nslab_chain= !!!!!!!!!!!!!!!!!!!!!!!!!!!"<<nslab_chain<<std::endl;
            assert(index <= nslab_chain);
            SHAREDPOINTER ret = s[index].alloc(mempool);
            assert(ret != sptr_NULL);
            return ret;
        }
        void free(size_t size, SHAREDPOINTER addr)
        {
            size_t index = ceil(log2(size));
            //std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!free index= " << index << "," << size << std::endl;
            assert(index <= nslab_chain);
            s[index].free(addr, mempool);
        }
        ~Slab()
        {
            if (manager_.rank == 0)
                for (size_t i = 0; i <= nslab_chain; i++)
                {
                    s[i].destroy(mempool);
                }
        }
    };
}