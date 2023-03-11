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
        //SharedPointer(T *ptr) : offset((size_t)ptr - MemStart<I>::start){};
        SharedPointer(T *ptr) : offset(reinterpret_cast<size_t>(ptr) - MemStart<I>::start){};
        SharedPointer(SHAREDPOINTER offset_in) : offset(offset_in){};
        SharedPointer(SharedPointer<T, I> &sptr) : offset(sptr.offset){};
        inline T &operator*() const
        {
            return *reinterpret_cast<T *>(MemStart<I>::start + offset);
            //return *((T *)(MemStart<I>::start + offset));
        }
        inline T *operator->()
        {
            return reinterpret_cast<T *>(MemStart<I>::start + offset);
            //return (T *)(MemStart<I>::start + offset);
        }
        inline T *operator()()
        {
            return reinterpret_cast<T *>(MemStart<I>::start + offset);
            //return (T *)(MemStart<I>::start + offset);
        }

        inline const T *operator()() const
        {
            return reinterpret_cast<const T *>(MemStart<I>::start + offset);
            //return (const T *)(MemStart<I>::start + offset);
        }
        inline operator SHAREDPOINTER()
        {
            return offset;
        }
        inline SHAREDPOINTER operator=(SHAREDPOINTER &a)
        {
            return offset = a;
        }
        inline SHAREDPOINTER operator=(const SHAREDPOINTER &a)
        {
            return offset = a;
        }
        inline SHAREDPOINTER operator=(const SharedPointer<T, I> &a)
        {
            return offset = a.offset;
        }
        inline SHAREDPOINTER operator+(const SHAREDPOINTER &a)
        {
            return offset + a;
        }
        inline bool operator==(const SharedPointer<T, I> &a)
        {
            return this->offset == a.offset;
        }
        inline bool isNULL()
        {
            return offset == sptr_NULL;
        }
        inline bool notNULL()
        {
            return offset != sptr_NULL;
        }
        inline void setNULL()
        {
            offset = sptr_NULL;
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
