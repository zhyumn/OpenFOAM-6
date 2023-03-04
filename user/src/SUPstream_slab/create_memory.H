#ifndef CREATE_MEMORY_HPP_
#define CREATE_MEMORY_HPP_

#include <stddef.h>
#include <unistd.h>
#include <math.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <cstring>
#include <stdio.h>
#include "SUPstream.H"
#define SHAREDPOINTER size_t

namespace Foam
{
    extern SHAREDPOINTER sptr_NULL;
    extern size_t pagesize;
    template <size_t I>
    struct MemStart
    {
        static size_t start;
    };

    template <class T, size_t I = 0>
    struct SharedPointer
    {
        SHAREDPOINTER offset;
        SharedPointer() : offset(sptr_NULL){};
        SharedPointer(T *ptr) : offset((size_t)ptr - MemStart<I>::start){};
        SharedPointer(SHAREDPOINTER offset_in) : offset(offset_in){};
        SharedPointer(SharedPointer<T, I> &sptr) : offset(sptr.offset){};
        T &operator*() const
        {
            return *((T *)(MemStart<I>::start + offset));
        }
        T *operator->()
        {
            return (T *)(MemStart<I>::start + offset);
        }
        T *operator()()
        {
            return (T *)(MemStart<I>::start + offset);
        }

        const T *operator()() const
        {
            return (const T *)(MemStart<I>::start + offset);
        }
        operator SHAREDPOINTER()
        {
            return offset;
        }
        SHAREDPOINTER operator=(SHAREDPOINTER &a)
        {
            return offset = a;
        }
        SHAREDPOINTER operator=(const SHAREDPOINTER &a)
        {
            return offset = a;
        }
        SHAREDPOINTER operator=(SharedPointer<T, I> &a)
        {
            return offset = a.offset;
        }
        SHAREDPOINTER operator+(const SHAREDPOINTER &a)
        {
            return offset + a;
        }
        bool isNULL()
        {
            return offset == -1;
        }
    };

    class MemPool
    {
    public:
        size_t npage;
        size_t memsize;
        Foam::SUPstream::shared_data<uint8_t> memory_MPI;
        Foam::SUPstream::shared_data<size_t> page_info;
        void *memory;

        //MemPool() :memory(NULL), page_info(NULL), npage(0), memsize(0) {}
        MemPool(Foam::SUPstream::mpi_manager &manager_in, size_t size);

        void create_memory(size_t size);

        SHAREDPOINTER get_page(size_t size);

        int free_page(SHAREDPOINTER, size_t);

        void free_memory();
        void report();
        ~MemPool() { free_memory(); }
    };
}

#endif //!CREATE_MEMORT_HPP_
