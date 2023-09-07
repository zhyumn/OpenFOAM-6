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

struct Node
{
    int value;
    SharedPointer<Node> next;
};
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
        pslab = new Foam::Slab(SUPstream::node_manager, sizeof(Node) * 100, sizeof(Node));
        SUPstream::shared_data<SharedPointer<Node>> head(SUPstream::node_manager);
        head() = sptr_NULL;
        SUPstream::Sync();
        if (SUPstream::node_manager.rank == 0)
        {
            SharedPointer<Node> tail(sptr_NULL);
            head() = pslab->alloc(sizeof(Node));
            tail = head();

            tail->value = 0;

            for (int i = 1; i < 10; i++)
            {
                tail->next = pslab->alloc(sizeof(Node));
                tail = tail->next;
                tail->value = i;
                tail->next = sptr_NULL;
            }
        }
        SUPstream::Sync();
        if (SUPstream::node_manager.rank == 1)
        {
            SharedPointer<Node> ptr(head());
            while (ptr.notNULL())
            {
                Pout << ptr->value << ",";
                ptr = ptr->next;
            }
            Pout << endl;
        }
        SUPstream::Sync();
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
