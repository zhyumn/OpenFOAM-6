#include "create_memory.H"
#define UNLIKELY(exp) __builtin_expect(exp, 0)

namespace Foam
{
    size_t pagesize = (size_t)sysconf(_SC_PAGESIZE);
    //SharedPointer<void*> sptr_NULL((size_t)-1);
    SHAREDPOINTER sptr_NULL = -1;
    template <>
    size_t MemStart<0>::start = 0;
    //size_t SharedPointer<void*, 0>::start = 0;

    //MemPool(size_t size) :memory(NULL), page_info(NULL), npage(0), memsize(0) { create_memory(size); }
    MemPool::MemPool(Foam::SUPstream::mpi_manager &manager_in, size_t size) : npage(ceil(size / (float)pagesize)),
                                                                              memsize(npage * pagesize),
                                                                              memory_MPI(manager_in, memsize),
                                                                              page_info(manager_in, npage),
                                                                              memory(memory_MPI.ptr())
    {
    }

    SHAREDPOINTER MemPool::get_page(size_t size)
    {
        size_t n_page = ceil(size / (float)pagesize);
        size_t empty_page = -1;
        bool flag = false;

        size_t j = 0;
        for (size_t i = 0; i + n_page <= npage; i += n_page) //n_page)
        {
            for (j = 0; j < n_page && page_info[i + j] == 0; j++)
                ;

            if (j == n_page)
            {
                empty_page = i;
                flag = true;
                break;
            }
        }
        if (flag == false)
        {
            return sptr_NULL;
        }
        for (size_t i = 0; i < n_page; i++)
        {
            page_info[empty_page + i] = 1;
        }
        return (SHAREDPOINTER)(pagesize * empty_page);
    }

    int MemPool::free_page(SHAREDPOINTER addr, size_t size)
    {
        size_t n_page = ceil(size / (float)pagesize);
        memset((char *)memory + addr, 0, n_page * pagesize);
        size_t offset = addr / pagesize;
        for (size_t i = 0; i < n_page; i++)
            page_info[offset + i] = 0;
        return 0;
    }

    void MemPool::free_memory()
    {
        if (memory == NULL)
            return;
        //if (UNLIKELY(munmap(memory, memsize) == -1))
        //    perror("free_memory");
    }

    void MemPool::report()
    {
        int sum = 0;
        for (size_t i = 0; i < npage; i++)
        {
            if (page_info[i] == 0)
                sum++;
        }
        sum *= pagesize;
        std::cout << "Mempool" << std::endl;
        std::cout << "Total:" << memsize << " Avai:" << sum << " used:" << memsize - sum << std::endl;
    }
}