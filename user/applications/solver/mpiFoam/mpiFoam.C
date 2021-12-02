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
#include "functionTabulation.H"
//#include "SUPstream.H"
#include "fvCFD.H"
#include "psiReactionThermo.H"
#include "thermoPhysicsTypes.H"
#include <cstdlib>
#include <ctime>
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char* argv[])
{
#include "setRootCase.H"
#include "createTime.H"
    if (SUPstream::parRun())
    {
        Info << "\nMPIrun !\n"
            << endl;
        SUPstream::node_init();
        SUPstream::shared_data<unsigned int> int_share(SUPstream::node_manager);
        std::atomic<unsigned int>* table_a = reinterpret_cast<std::atomic<unsigned int>*>(int_share.ptr());
        int_share() = 0;
        SUPstream::Sync();
        for (int i = 0; i < 1000; i++)//10000000
        {
            //mutexmy.lock();
            //int_share()++;
            (*table_a)++;
            //mutexmy.unlock();
        }
        SUPstream::Sync();
        Pout << " a = " << int_share() << endl;

        functionTabulation test;
        //parallelISAT_chem<psiReactionThermo, constGasHThermoPhysics> test(SUPstream::node_manager, 1000000, SUPstream::Sync);

        //Info<<"Test! "<<test.tleaf.v<<endl;
        //Info<<"Test! "<<test.tleaf.node_->v<<endl;
        //Info<<"Test! "<<test.tleaf.node_->leaf_->v<<endl;
        SUPstream::Sync();
        //Pout << SUPstream::node_manager.rank << endl;
        //if (SUPstream::node_manager.rank == 0)
        //label index;
        scalar input;
        scalar out_approx;
        scalar out_exact;
        scalar maxerror = 0;
        if (1)
        {
            //std::srand(std::time(nullptr)*SUPstream::node_manager.rank);
            std::srand(SUPstream::node_manager.rank);
            //Pout << SUPstream::node_manager.rank << endl;

            for (int i = 0;i < 1000;i++)
            {
                int random_variable = std::rand();

                input = (random_variable % 1000) / 100.0;
                //index = test.binaryTreeSearch(input);
                out_approx = test.calculate(input);
                out_exact = test.func(input);
                Pout << input << "," << out_exact << "," << out_approx << ",error = " << out_exact - out_approx << endl;
                if (fabs(out_exact - out_approx) > maxerror)
                    maxerror = fabs(out_exact - out_approx);

            }
            Pout << "maxerror=" << maxerror << endl;

            //test.insert(2);

        }
        SUPstream::Sync();
        if (SUPstream::node_manager.rank == 0)
        {   
            test.writeDot("tree.gv");
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
