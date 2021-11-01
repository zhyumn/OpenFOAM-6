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
    mpiFoam

Description

\*---------------------------------------------------------------------------*/
#include "ISAT.H"
//#include "SUPstream.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#include "setRootCase.H"
    if (SUPstream::parRun())
    {
        Info << "\nMPIrun !\n"
             << endl;
        SUPstream::node_init();
        SUPstream::shared_data<unsigned int> int_share(SUPstream::node_manager);
        std::atomic<unsigned int>* table_a=(std::atomic<unsigned int>*)(int_share.ptr());
        int_share()=0;
        SUPstream::Sync();
        for (int i = 0; i < 10000000; i++)
        {
            //mutexmy.lock();
            //int_share()++;
            (*table_a)++;
            //mutexmy.unlock();
        }
        SUPstream::Sync();
        Pout<<" a = "<<int_share()<<endl;

        ISAT test(SUPstream::node_manager,10);

        Info<<"Test! "<<test.tleaf.v<<endl;
        Info<<"Test! "<<test.tleaf.node_->v<<endl;
        Info<<"Test! "<<test.tleaf.node_->leaf_->v<<endl;
    }
#include "createTime.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    

    Info << nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
         << "  ClockTime = " << runTime.elapsedClockTime() << " s"
         << nl << endl;

    Info << "End\n"
         << endl;

    return 0;
}

// ************************************************************************* //
