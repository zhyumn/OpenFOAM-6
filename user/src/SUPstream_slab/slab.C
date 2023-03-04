/* ref: https://github.com/bbu/userland-slab-allocator */

#include "slab.H"




#define SLAB_DUMP_COLOURED

#ifdef SLAB_DUMP_COLOURED
# define GRAY(s)   "\033[1;30m" s "\033[0m"
# define RED(s)    "\033[0;31m" s "\033[0m"
# define GREEN(s)  "\033[0;32m" s "\033[0m"
# define YELLOW(s) "\033[1;33m" s "\033[0m"
#else
# define GRAY(s)   s
# define RED(s)    s
# define GREEN(s)  s
# define YELLOW(s) s
#endif

#define SLOTS_ALL_ZERO ((uint64_t) 0)
#define SLOTS_FIRST ((uint64_t) 1)
#define FIRST_FREE_SLOT(s) ((size_t) __builtin_ctzll(s))
#define FREE_SLOTS(s) ((size_t) __builtin_popcountll(s))
#define ONE_USED_SLOT(slots, empty_slotmask) \
    ( \
        ( \
            (~(slots) & (empty_slotmask))       & \
            ((~(slots) & (empty_slotmask)) - 1)   \
        ) == SLOTS_ALL_ZERO \
    )

#define POWEROF2(x) ((x) != 0 && ((x) & ((x) - 1)) == 0)

#define LIKELY(exp) __builtin_expect(exp, 1)
#define UNLIKELY(exp) __builtin_expect(exp, 0)

namespace Foam
{

#ifndef NDEBUG
    int slab_chain::is_valid()
    {
        assert(POWEROF2(pagesize));
        //std::cout<<"this->slabsize="<<this->slabsize<<std::endl;
        assert(POWEROF2(this->slabsize));
        assert(POWEROF2(this->pages_per_alloc));

        assert(this->itemcount >= 2 && this->itemcount <= 64);
        assert(this->itemsize >= 1 && this->itemsize <= SIZE_MAX);
        assert(this->pages_per_alloc >= pagesize);
        assert(this->pages_per_alloc >= this->slabsize);

        assert(offsetof(struct slab_header, data) +
            this->itemsize * this->itemcount <= this->slabsize);

        assert(this->empty_slotmask == ~SLOTS_ALL_ZERO >> (64 - this->itemcount));
        assert(this->initial_slotmask == (this->empty_slotmask ^ SLOTS_FIRST));
        assert(this->alignment_mask == ~(this->slabsize - 1));

        SharedPointer_SLAB  heads[] =
        { this->full, this->empty, this->partial };

        for (size_t head = 0; head < 3; ++head) {
            SharedPointer_SLAB prev, slab;

            for (slab = heads[head]; slab != sptr_NULL; slab = slab->next) {
                if (prev == sptr_NULL)
                    assert(slab->prev == sptr_NULL);
                else
                    assert(slab->prev == prev);

                switch (head) {
                case 0:
                    assert(slab->slots == SLOTS_ALL_ZERO);
                    break;

                case 1:
                    assert(slab->slots == this->empty_slotmask);
                    break;

                case 2:
                    assert((slab->slots & ~this->empty_slotmask) == SLOTS_ALL_ZERO);
                    assert(FREE_SLOTS(slab->slots) >= 1);
                    assert(FREE_SLOTS(slab->slots) < this->itemcount);
                    break;
                }

                if (slab->refcount == 0) {
                    assert((uintptr_t)slab % this->slabsize == 0);

                    if (this->slabsize >= pagesize)
                        assert((uintptr_t)slab->page % this->slabsize == 0);
                    else
                        assert((uintptr_t)slab->page % pagesize == 0);
                }
                else {
                    if (this->slabsize >= pagesize)
                        assert((uintptr_t)slab % this->slabsize == 0);
                    else
                        assert((uintptr_t)slab % pagesize == 0);
                }

                prev = slab;
            }
        }

        return 1;
    }
#endif

    void slab_chain::init(const size_t itemsize)
    {
        assert(itemsize >= 1 && itemsize <= SIZE_MAX);
        assert(POWEROF2(pagesize));

        this->itemsize = itemsize;

        const size_t data_offset = offsetof(struct slab_header, data);
        const size_t least_slabsize = data_offset + 64 * this->itemsize;
        this->slabsize = (size_t)1 << (size_t)ceil(log2(least_slabsize));
        this->itemcount = 64;

        if (this->slabsize - least_slabsize != 0) {
            const size_t shrinked_slabsize = this->slabsize >> 1;

            if (data_offset < shrinked_slabsize&&
                shrinked_slabsize - data_offset >= 2 * this->itemsize) {

                this->slabsize = shrinked_slabsize;
                this->itemcount = (shrinked_slabsize - data_offset) / this->itemsize;
            }
        }

        this->pages_per_alloc = this->slabsize > pagesize ?
            this->slabsize : pagesize;

        this->empty_slotmask = ~SLOTS_ALL_ZERO >> (64 - this->itemcount); //warning 64?
        this->initial_slotmask = this->empty_slotmask ^ SLOTS_FIRST;
        this->alignment_mask = ~(this->slabsize - 1);
        this->partial = sptr_NULL;
        this->empty = sptr_NULL;
        this->full = sptr_NULL;
        //std::cout<<"here!!!!!!!!!!!!!init"<<std::endl;
        assert(is_valid());
    }

    SHAREDPOINTER slab_chain::alloc(MemPool& mempool)
    {
        //std::cout<<"here!!!!!!!!!!!!!alloc start"<<std::endl;
        assert(is_valid());

        if (LIKELY(this->partial != sptr_NULL)) {
            /* found a partial slab, locate the first free slot */
            register const size_t slot = FIRST_FREE_SLOT(this->partial->slots);
            this->partial->slots ^= SLOTS_FIRST << slot;

            if (UNLIKELY(this->partial->slots == SLOTS_ALL_ZERO)) {
                /* slab has become full, change state from partial to full */
                SharedPointer_SLAB tmp = this->partial;

                /* skip first slab from partial list */
                if (LIKELY((this->partial = this->partial->next) != sptr_NULL))
                    this->partial->prev = sptr_NULL;

                if (LIKELY((tmp->next = this->full) != sptr_NULL))
                    this->full->prev = tmp;

                this->full = tmp;
                //std::cout<<"here!!!!!!!!!!!!!alloc end"<<std::endl;
                return (SHAREDPOINTER)(this->full->data + slot * this->itemsize - MemStart<0>::start);
            }
            else {
                //std::cout<<"here!!!!!!!!!!!!!alloc end"<<std::endl;
                return (SHAREDPOINTER)(this->partial->data + slot * this->itemsize - MemStart<0>::start);
            }
        }
        else if (LIKELY((this->partial = this->empty) != sptr_NULL)) {
            /* found an empty slab, change state from empty to partial */
            if (LIKELY((this->empty = this->empty->next) != sptr_NULL))
                this->empty->prev = sptr_NULL;

            this->partial->next = sptr_NULL;

            /* slab is located either at the beginning of page, or beyond */
            UNLIKELY(this->partial->refcount != 0) ?
                this->partial->refcount++ : this->partial->page->refcount++;

            this->partial->slots = this->initial_slotmask;
            //std::cout<<"here!!!!!!!!!!!!!alloc end"<<std::endl;
            return (SHAREDPOINTER)(this->partial->data - MemStart<0>::start);
        }
        else {
            /* no empty or partial slabs available, create a new one */


            this->partial = mempool.get_page(this->pages_per_alloc);

            if (UNLIKELY(this->partial == sptr_NULL))
               {
                //std::cout<<"here!!!!!!!!!!!!!alloc end"<<std::endl;
                 return perror("get_page failed"), this->partial = sptr_NULL;}

            SharedPointer_SLAB prev;

            const char* const page_end =
                (char*)(this->partial + this->pages_per_alloc);

            union {
                const char* c;
                SharedPointer_SLAB  s;
            } curr = {
                .c = (const char*)(this->partial + this->slabsize)
            };

            __builtin_prefetch(this->partial(), 1);

            this->partial->prev = sptr_NULL;
            this->partial->next = sptr_NULL;
            this->partial->refcount = 1;
            this->partial->slots = this->initial_slotmask;

            if (LIKELY(curr.c != page_end)) {
                curr.s->prev = sptr_NULL;
                curr.s->refcount = 0;
                curr.s->page = this->partial;
                curr.s->slots = this->empty_slotmask;
                this->empty = curr.s;
                prev = curr.s;

                while (LIKELY((curr.c += this->slabsize) != page_end)) {
                    prev->next = curr.s;
                    curr.s->prev = prev;
                    curr.s->refcount = 0;
                    curr.s->page = this->partial;
                    curr.s->slots = this->empty_slotmask;
                    prev = curr.s;
                }

                prev->next = sptr_NULL;
            }
            //std::cout<<"here!!!!!!!!!!!!!alloc end"<<std::endl;
            return (SHAREDPOINTER)(this->partial->data - MemStart<0>::start);
        }

        /* unreachable */
    }

    void slab_chain::free(SHAREDPOINTER addr,MemPool& mempool)
    {
        //std::cout<<"here!!!!!!!!!!!!!free"<<std::endl;
        assert(is_valid());
        assert(addr != sptr_NULL);

        SharedPointer_SLAB slab((uintptr_t)addr & this->alignment_mask);

        register const int slot = ((char*)addr - (char*)slab.offset -
            offsetof(struct slab_header, data)) / this->itemsize;

        if (UNLIKELY(slab->slots == SLOTS_ALL_ZERO)) {
            /* target slab is full, change state to partial */
            slab->slots = SLOTS_FIRST << slot;

            if (LIKELY(slab != this->full)) {
                if (LIKELY((slab->prev->next = slab->next) != sptr_NULL))
                    slab->next->prev = slab->prev;

                slab->prev = sptr_NULL;
            }
            else if (LIKELY((this->full = this->full->next) != sptr_NULL)) {
                this->full->prev = sptr_NULL;
            }

            slab->next = this->partial;

            if (LIKELY(this->partial != sptr_NULL))
                this->partial->prev = slab;

            this->partial = slab;
        }
        else if (UNLIKELY(ONE_USED_SLOT(slab->slots, this->empty_slotmask))) {
            /* target slab is partial and has only one filled slot */
            if (UNLIKELY(slab->refcount == 1 || (slab->refcount == 0 &&
                slab->page->refcount == 1))) {

                /* unmap the whole page if this slab is the only partial one */
                if (LIKELY(slab != this->partial)) {
                    if (LIKELY((slab->prev->next = slab->next) != sptr_NULL))
                        slab->next->prev = slab->prev;
                }
                else if (LIKELY((this->partial = this->partial->next) != sptr_NULL)) {
                    this->partial->prev = sptr_NULL;
                }

                SharedPointer_SLAB page = UNLIKELY(slab->refcount != 0) ? slab : slab->page;
                const char* const page_end = (char*)(page + this->pages_per_alloc);
                char found_head = 0;

                union {
                    char* c;
                    SharedPointer_SLAB s;
                } s = {
                .c = (char*)sptr_NULL
                };

                for (s.s = page; s.c != page_end; s.c += this->slabsize) {
                    if (UNLIKELY(s.s == this->empty))
                        found_head = 1;
                    else if (UNLIKELY(s.s == slab))
                        continue;
                    else if (LIKELY((s.s->prev->next = s.s->next) != sptr_NULL))
                        s.s->next->prev = s.s->prev;
                }

                if (UNLIKELY(found_head && (this->empty = this->empty->next) != sptr_NULL))
                    this->empty->prev = sptr_NULL;
                if (UNLIKELY(mempool.free_page(page, this->pages_per_alloc) == -1))
                    perror("free_page");

            }
            else {
                slab->slots = this->empty_slotmask;

                if (LIKELY(slab != this->partial)) {
                    if (LIKELY((slab->prev->next = slab->next) != sptr_NULL))
                        slab->next->prev = slab->prev;

                    slab->prev = sptr_NULL;
                }
                else if (LIKELY((this->partial = this->partial->next) != sptr_NULL)) {
                    this->partial->prev = sptr_NULL;
                }

                slab->next = this->empty;

                if (LIKELY(this->empty != sptr_NULL))
                    this->empty->prev = slab;

                this->empty = slab;

                UNLIKELY(slab->refcount != 0) ?
                    slab->refcount-- : slab->page->refcount--;
            }
        }
        else {
            /* target slab is partial, no need to change state */
            slab->slots |= SLOTS_FIRST << slot;
        }
    }

    void slab_chain::traverse(void (*fn)(const void*))
    {
        //std::cout<<"here!!!!!!!!!!!!!traverse"<<std::endl;
        assert(fn != NULL);
        assert(is_valid());

        SharedPointer_SLAB slab;
        SharedPointer_SLAB item, end;
        const size_t data_offset = offsetof(struct slab_header, data);

        for (slab = this->partial; slab != sptr_NULL; slab = slab->next) {
            item = slab + data_offset;
            end = item + this->itemcount * this->itemsize;
            uint64_t mask = SLOTS_FIRST;

            do {
                if (!(slab->slots & mask))
                    fn(item());

                mask <<= 1;
            } while ((item.offset += this->itemsize) != end.offset);
        }

        for (slab = this->full; slab != sptr_NULL; slab = slab->next) {
            item = slab + data_offset;
            end = item + this->itemcount * this->itemsize;

            do fn(item());
            while ((item.offset += this->itemsize) != end.offset);
        }
    }

    void slab_chain::destroy(MemPool& mempool)
    {
        //std::cout<<"here!!!!!!!!!!!!!destroy"<<std::endl;
        assert(is_valid());

        SharedPointer_SLAB heads[] = { this->partial, this->empty, this->full };
        SharedPointer_SLAB pages_head, pages_tail;

        for (size_t i = 0; i < 3; ++i) {
            SharedPointer_SLAB slab = heads[i];

            while (slab != sptr_NULL) {
                if (slab->refcount != 0) {
                    SharedPointer_SLAB page = slab;
                    slab = slab->next;

                    if (UNLIKELY(pages_head == sptr_NULL))
                        pages_head = page;
                    else
                        pages_tail->next = page;

                    pages_tail = page;
                }
                else {
                    slab = slab->next;
                }
            }
        }

        if (LIKELY(pages_head != sptr_NULL)) {
            pages_tail->next = sptr_NULL;
            SharedPointer_SLAB page = pages_head;


            do {
                SharedPointer_SLAB target = page;
                page = page->next;

                if (UNLIKELY(mempool.free_page(target, this->pages_per_alloc) == -1))
                    perror("free_page");
            } while (page != sptr_NULL);
        }
    }

    void slab_chain::dump(FILE* const out)
    {
                //std::cout<<"here!!!!!!!!!!!!!dump"<<std::endl;
        assert(out != NULL);
        assert(is_valid());

        SharedPointer_SLAB heads[] =
        { this->partial, this->empty, this->full };

        const char* labels[] = { "part", "empt", "full" };

        for (size_t i = 0; i < 3; ++i) {
            SharedPointer_SLAB slab = heads[i];

            fprintf(out,
                YELLOW("%6s ") GRAY("|%2d%13s|%2d%13s|%2d%13s|%2d%13s") "\n",
                labels[i], 64, "", 48, "", 32, "", 16, "");

            unsigned long long total = 0, row;

            for (row = 1; slab != sptr_NULL; slab = slab->next, ++row) {
                const unsigned used = this->itemcount - FREE_SLOTS(slab->slots);
                fprintf(out, GRAY("%6llu "), row);

                for (int k = 63; k >= 0; --k) {
                    fprintf(out, slab->slots & (SLOTS_FIRST << k) ? GREEN("1") :
                        ((size_t)k >= this->itemcount ? GRAY("0") : RED("0")));
                }

                fprintf(out, RED(" %8u") "\n", used);
                total += used;
            }

            fprintf(out,
                GREEN("%6s ") GRAY("^%15s^%15s^%15s^%15s") YELLOW(" %8llu") "\n\n",
                "", "", "", "", "", total);
        }
    }

    void slab_chain::stats(FILE* const out)
    {
                        //std::cout<<"here!!!!!!!!!!!!!stats"<<std::endl;
        assert(out != NULL);
        assert(is_valid());

        long long unsigned
            total_nr_slabs = 0,
            total_used_slots = 0,
            total_free_slots = 0;

        float occupancy;

        SharedPointer_SLAB heads[] =
        { this->partial, this->empty, this->full };

        const char* labels[] = { "Partial", "Empty", "Full" };

        fprintf(out, "%8s %17s %17s %17s %17s\n", "",
            "Slabs", "Used", "Free", "Occupancy");

        for (size_t i = 0; i < 3; ++i) {
            long long unsigned nr_slabs = 0, used_slots = 0, free_slots = 0;
            SharedPointer_SLAB slab;

            for (slab = heads[i]; slab != sptr_NULL; slab = slab->next) {
                nr_slabs++;
                used_slots += this->itemcount - FREE_SLOTS(slab->slots);
                free_slots += FREE_SLOTS(slab->slots);
            }

            occupancy = used_slots + free_slots ?
                100 * (float)used_slots / (used_slots + free_slots) : 0.0;

            fprintf(out, "%8s %17llu %17llu %17llu %16.2f%%\n",
                labels[i], nr_slabs, used_slots, free_slots, occupancy);

            total_nr_slabs += nr_slabs;
            total_used_slots += used_slots;
            total_free_slots += free_slots;
        }

        occupancy = total_used_slots + total_free_slots ?
            100 * (float)total_used_slots / (total_used_slots + total_free_slots) :
            0.0;

        fprintf(out, "%8s %17llu %17llu %17llu %16.2f%%\n", "Total",
            total_nr_slabs, total_used_slots, total_free_slots, occupancy);
    }

    void slab_chain::report()
    {
                        //std::cout<<"here!!!!!!!!!!!!!stats"<<std::endl;
        assert(is_valid());

        long long unsigned
            total_nr_slabs = 0,
            total_used_slots = 0,
            total_free_slots = 0;

        float occupancy;

        SharedPointer_SLAB heads[] =
        { this->partial, this->empty, this->full };

        const char* labels[] = { "Partial", "Empty", "Full" };
        std::cout<< "itemsize:" << this->itemsize<<std::endl;
        std::cout<< "Slabs Used Free Occupancy" <<std::endl;
        //fprintf(out, "%8s %17s %17s %17s %17s\n", "",
        //    "Slabs", "Used", "Free", "Occupancy");

        for (size_t i = 0; i < 3; ++i) {
            long long unsigned nr_slabs = 0, used_slots = 0, free_slots = 0;
            SharedPointer_SLAB slab;

            for (slab = heads[i]; slab != sptr_NULL; slab = slab->next) {
                nr_slabs++;
                used_slots += this->itemcount - FREE_SLOTS(slab->slots);
                free_slots += FREE_SLOTS(slab->slots);
            }

            occupancy = used_slots + free_slots ?
                100 * (float)used_slots / (used_slots + free_slots) : 0.0;

            std::cout<< labels[i]<< " " << nr_slabs<< " " << used_slots<< " " << free_slots<< " " << occupancy <<std::endl;
            //fprintf(out, "%8s %17llu %17llu %17llu %16.2f%%\n",
            //    labels[i], nr_slabs, used_slots, free_slots, occupancy);

            total_nr_slabs += nr_slabs;
            total_used_slots += used_slots;
            total_free_slots += free_slots;
        }

        occupancy = total_used_slots + total_free_slots ?
            100 * (float)total_used_slots / (total_used_slots + total_free_slots) :
            0.0;
        std::cout<<  "Total"<< " " << total_nr_slabs<< " " <<  total_used_slots<< " " << total_free_slots<< " " << occupancy <<std::endl;
        std::cout<< std::endl;
        //fprintf(out, "%8s %17llu %17llu %17llu %16.2f%%\n", "Total",
        //    total_nr_slabs, total_used_slots, total_free_slots, occupancy);
    }

    void slab_chain::report(int &total, int &used,int &free_)
    {
                        //std::cout<<"here!!!!!!!!!!!!!stats"<<std::endl;
        assert(is_valid());

        long long unsigned
            total_nr_slabs = 0,
            total_used_slots = 0,
            total_free_slots = 0;

        float occupancy;

        SharedPointer_SLAB heads[] =
        { this->partial, this->empty, this->full };

        const char* labels[] = { "Partial", "Empty", "Full" };
        //std::cout<< "itemsize:" << this->itemsize<<std::endl;
        //std::cout<< "Slabs Used Free Occupancy" <<std::endl;
        //fprintf(out, "%8s %17s %17s %17s %17s\n", "",
        //    "Slabs", "Used", "Free", "Occupancy");

        for (size_t i = 0; i < 3; ++i) {
            long long unsigned nr_slabs = 0, used_slots = 0, free_slots = 0;
            SharedPointer_SLAB slab;

            for (slab = heads[i]; slab != sptr_NULL; slab = slab->next) {
                nr_slabs++;
                used_slots += this->itemcount - FREE_SLOTS(slab->slots);
                free_slots += FREE_SLOTS(slab->slots);
            }

            occupancy = used_slots + free_slots ?
                100 * (float)used_slots / (used_slots + free_slots) : 0.0;

            //std::cout<< labels[i]<< " " << nr_slabs<< " " << used_slots<< " " << free_slots<< " " << occupancy <<std::endl;
            //fprintf(out, "%8s %17llu %17llu %17llu %16.2f%%\n",
            //    labels[i], nr_slabs, used_slots, free_slots, occupancy);

            total_nr_slabs += nr_slabs;
            total_used_slots += used_slots;
            total_free_slots += free_slots;
        }
        total += total_nr_slabs* this->pages_per_alloc;
        used += total_used_slots * this->itemsize;
        free_ += total_free_slots *  this->itemsize;
        occupancy = total_used_slots + total_free_slots ?
            100 * (float)total_used_slots / (total_used_slots + total_free_slots) :
            0.0;
        //std::cout<<  "Total"<< " " << total_nr_slabs<< " " <<  total_used_slots<< " " << total_free_slots<< " " << occupancy <<std::endl;
        //std::cout<< std::endl;
        //fprintf(out, "%8s %17llu %17llu %17llu %16.2f%%\n", "Total",
        //    total_nr_slabs, total_used_slots, total_free_slots, occupancy);
    }

    void slab_chain::props(FILE* const out)
    {
                        //std::cout<<"here!!!!!!!!!!!!!props"<<std::endl;
        assert(out != NULL);
        assert(is_valid());

        fprintf(out,
            "%18s: %8zu\n"
            "%18s: %8zu = %.2f * (%zu pagesize)\n"
            "%18s: %8zu = (%zu offset) + (%zu itemcount) * (%zu itemsize)\n"
            "%18s: %8zu = (%zu slabsize) - (%zu total)\n"
            "%18s: %8zu = %zu * (%zu pagesize)\n"
            "%18s: %8zu = (%zu alloc) / (%zu slabsize)\n",

            "pagesize",
            pagesize,

            "slabsize",
            this->slabsize, (float)this->slabsize / pagesize, pagesize,

            "total",
            offsetof(struct slab_header, data) + this->itemcount * this->itemsize,
            offsetof(struct slab_header, data), this->itemcount, this->itemsize,

            "waste per slab",
            this->slabsize - offsetof(struct slab_header, data) -
            this->itemcount * this->itemsize, this->slabsize,
            offsetof(struct slab_header, data) + this->itemcount * this->itemsize,

            "pages per alloc",
            this->pages_per_alloc, this->pages_per_alloc / pagesize,
            pagesize,

            "slabs per alloc",
            this->pages_per_alloc / this->slabsize, this->pages_per_alloc,
            this->slabsize
        );
    }
}


