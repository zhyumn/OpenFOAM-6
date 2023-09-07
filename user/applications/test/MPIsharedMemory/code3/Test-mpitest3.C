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

Application
    

Description

\*---------------------------------------------------------------------------*/

#include "SUPstream.H"
#include "slab.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#include "setRootCase.H"
#include "createTime.H"
    if (SUPstream::parRun())
    {
        Info << "\nMPIrun !\n"
             << endl;

        SUPstream::node_init();
        double t1, t2;

        SUPstream::shared_data<int> a(SUPstream::node_manager);

        SUPstream::Sync();

        if (SUPstream::node_manager.rank == 0)
        {
            a() = 100;
        }
        SUPstream::Sync();

        Pout << "a=" << a() << endl;

        SUPstream::Sync();
        a() = 0;
        SUPstream::Sync();
        t1 = MPI_Wtime();
        for (int i = 0; i < 10000; i++)
        {
            a()++;
        }
        SUPstream::Sync();
        t2 = MPI_Wtime();
        Pout << "test1, a=" << a() << ", correct value should be " << SUPstream::node_manager.size * 10000 << ", run time " << t2 - t1 << endl;

        SUPstream::Sync();
        a() = 0;
        SUPstream::Sync();
        std::atomic<int> *atomic_a = reinterpret_cast<std::atomic<int> *>(&a());
        t1 = MPI_Wtime();
        for (int i = 0; i < 10000; i++)
        {
            atomic_a->operator++();
        }
        SUPstream::Sync();
        t2 = MPI_Wtime();
        Pout << "test2, a=" << a() << ", correct value should be " << SUPstream::node_manager.size * 10000 << ", run time " << t2 - t1 << endl;

        SUPstream::Sync();
        a() = 0;
        SUPstream::Sync();
        SUPstream::mpi_mutex mutex(SUPstream::node_manager, SUPstream::Sync);
        t1 = MPI_Wtime();
        for (int i = 0; i < 10000; i++)
        {
            mutex.lock();
            a()++;
            mutex.unlock();
        }
        SUPstream::Sync();
        t2 = MPI_Wtime();
        Pout << "test3, a=" << a() << ", correct value should be " << SUPstream::node_manager.size * 10000 << ", run time " << t2 - t1 << endl;
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
         << "  ClockTime = " << runTime.elapsedClockTime() << " s"
         << nl << endl;

    Info << "End\n"
         << endl;

    return 0;
}

// ************************************************************************* //
