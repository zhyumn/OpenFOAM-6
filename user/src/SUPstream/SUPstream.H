/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::SUPstream

Description

SourceFiles
    SUPstreamI.H
    SUPstream.C
    SUPstreamIO.C

\*---------------------------------------------------------------------------*/

#ifndef SUPstream_H
#define SUPstream_H

#include "UPstream.H"

#include <mpi.h>
#include <atomic>
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Forward declaration of classes
    class Istream;
    class Ostream;

    // Forward declaration of friend functions and operators
    class SUPstream;
    Istream& operator>>(Istream&, SUPstream&);
    Ostream& operator<<(Ostream&, const SUPstream&);

    /*---------------------------------------------------------------------------*\
                         Class SUPstream Declaration
\*---------------------------------------------------------------------------*/

    class SUPstream
        : public UPstream
    {
        // Private data

        //- Description of data_
        scalar data_;

        // Private Member Functions

        //- Disallow default bitwise copy construct
        SUPstream(const SUPstream&);

        //- Disallow default bitwise assignment
        void operator=(const SUPstream&);

    public:
        // Public classes

        class mpi_manager
        {
        public:
            MPI_Comm comm;
            int size, rank;
        };
        class mpi_node_manager : public mpi_manager
        {
        public:
            mpi_node_manager(MPI_Comm allcomm)
            {
                int rank_t = 0;
                MPI_Comm_rank(allcomm, &rank);
                MPI_Comm_split_type(allcomm, MPI_COMM_TYPE_SHARED, rank_t, MPI_INFO_NULL, &comm);
                MPI_Comm_size(comm, &size);
                MPI_Comm_rank(comm, &rank);
            }
            mpi_node_manager() {};
            void init(MPI_Comm allcomm)
            {
                int rank_t = 0;
                MPI_Comm_rank(allcomm, &rank);
                MPI_Comm_split_type(allcomm, MPI_COMM_TYPE_SHARED, rank_t, MPI_INFO_NULL, &comm);
                MPI_Comm_size(comm, &size);
                MPI_Comm_rank(comm, &rank);
            }
        };

        class mpi_sync
        {
            int sync_, sync_r_;

        public:
            mpi_manager& manager;
            mpi_sync(mpi_manager& manager_in)
                : sync_(0), sync_r_(0), manager(manager_in) {}
            void sync()
            {
                MPI_Reduce(&sync_, &sync_r_, 1, MPI_INT, MPI_MAX, 0, manager.comm);
                MPI_Barrier(MPI_COMM_WORLD);
            }
            void operator()()
            {
                sync();
            }
        };
        class my_mutex
        {
            std::atomic_flag& flag;

        public:
            my_mutex(std::atomic_flag* ptr, mpi_sync& sync_in)
                : flag(*ptr)
            {
                if (sync_in.manager.rank == 0)
                    flag.clear();
                sync_in.sync();
            }
            void lock()
            {
                while (flag.test_and_set(std::memory_order_acquire))
                    ;
            }
            void unlock() { flag.clear(std::memory_order_release); }
        };


        class my_room
        {
            std::atomic_flag& flag;
            std::atomic<int>& p_in_room;

        public:
            my_room(std::atomic_flag* ptr, std::atomic<int>* ptr2, mpi_sync& sync_in)
                : flag(*ptr), p_in_room(*ptr2)
            {
                if (sync_in.manager.rank == 0)
                {
                    flag.clear();
                    p_in_room = 0;
                }
                sync_in.sync();
            }
            void lock()
            {
                while (flag.test_and_set(std::memory_order_acquire))
                    ;
            }
            void unlock() { flag.clear(std::memory_order_release); }
            void enter()
            {
                while (flag.test_and_set(std::memory_order_acquire))
                    ;
                p_in_room++;
                flag.clear(std::memory_order_release);
            }
            void exit() { p_in_room--; }
            void give_back()
            {
                flag.clear(std::memory_order_release);
                p_in_room--;
            }
            void require_empty_room()
            {
                while (flag.test_and_set(std::memory_order_acquire))
                    ;
                while (p_in_room != 0)
                    ;
                p_in_room++;
            }
        };
        template <typename Datatype>
        class shared_data
        {
            Datatype* table;
            MPI_Win wintable;
            mpi_manager& manager;
            int size;

        public:
            shared_data(mpi_manager& manager_in, int size_in = 1)
                : manager(manager_in), size(size_in)
            {
                int localsize = 0, windisp;
                MPI_Aint winsize;
                if (manager.rank == 0)
                    localsize = size;
                MPI_Win_allocate_shared(localsize * sizeof(Datatype), sizeof(Datatype), MPI_INFO_NULL, manager.comm, &table, &wintable);
                if (manager.rank != 0)
                {
                    MPI_Win_shared_query(wintable, 0, &winsize, &windisp, &table);
                }
            }
            Datatype& operator()() { return *table; }
            Datatype& operator[](int n) { return table[n]; }
            Datatype* ptr() { return table; }
        };
        class mpi_mutex
        {
            shared_data<std::atomic_flag> flag_shared;
            my_mutex mutex_;

        public:
            mpi_mutex(mpi_manager& manager_in, mpi_sync& sync_in)
                : flag_shared(manager_in), mutex_(flag_shared.ptr(), sync_in) {}
            void lock() { mutex_.lock(); }
            void unlock() { mutex_.unlock(); }
        };

        class mpi_room
        {
            shared_data<std::atomic_flag> flag_shared;
            shared_data<std::atomic<int>> a_int_shared;

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

    public:
        // Static data members

        //- Static data staticData
        static scalar staticData;

        static mpi_node_manager node_manager;

        static mpi_sync Sync;

        // Static member functions
        static void node_init();

        // Constructors

        // Selectors

        //- Select null constructed

        //- Destructor

        // Member Functions

        // Access

        // Check

        // Edit

        // Write

        // Member Operators

        // Friend Functions

        // Friend Operators

        // IOstream Operators

        friend Istream& operator>>(Istream&, SUPstream&);
        friend Ostream& operator<<(Ostream&, const SUPstream&);
    };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}  // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#include "SUPstreamI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
