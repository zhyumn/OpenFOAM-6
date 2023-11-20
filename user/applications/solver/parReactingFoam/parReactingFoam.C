/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    reactingFoam

Description
    Solver for combustion with chemical reactions.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "turbulentFluidThermoModel.H"
#include "psiReactionThermo.H"
#include "CombustionModel.H"
#include "multivariateScheme.H"
#include "pimpleControl.H"
#include "pressureControl.H"
#include "fvOptions.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "slab.H"
#include "SUPstream.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#include "postProcess.H"

#include "setRootCaseLists.H"
#include "createTime.H"

    if (SUPstream::parRun())
    {
        size_t memorySize_G = runTime.controlDict().subDict("memorySize").lookupOrDefault<label>("G", 0);
        size_t memorySize_M = runTime.controlDict().subDict("memorySize").lookupOrDefault<label>("M", 300);
        size_t memorySize_K = runTime.controlDict().subDict("memorySize").lookupOrDefault<label>("K", 0);
        size_t memorySize_B = runTime.controlDict().subDict("memorySize").lookupOrDefault<label>("B", 0);
/*         size_t memorySize = memorySize_G;
        memorySize *= 1024;
        memorySize += memorySize_M;
        memorySize *= 1024;
        memorySize += memorySize_K;
        memorySize *= 1024;
        memorySize += memorySize_B; */
        size_t memorySize = memorySize_G * (1<<30) + memorySize_M * (1<<20) +  memorySize_K * (1<<10) + memorySize_B;
        //size_t memorySize = runTime.controlDict().lookupOrDefault<label>("memorySize", 60000*5000);
        size_t maxMemBlock = runTime.controlDict().lookupOrDefault<label>("maxMemBlock", 60000);

        //FatalErrorInFunction << "memorySize = " << memorySize << ",maxMemBlock = " << maxMemBlock << exit(FatalError);

        SUPstream::node_init();
        //Foam::Slab slab(SUPstream::node_manager, 1000, 100000);
        //size_t memsize =60000;
        //memsize *=5000;
        //pslab = new Foam::Slab(SUPstream::node_manager, 60000, memsize);
        pslab = new Foam::Slab(SUPstream::node_manager, maxMemBlock, memorySize);
    }
    //FatalErrorInFunction  << exit(FatalError);

    /*
    {
        if (SUPstream::node_manager.rank == 0)
        {
            int a0_0 = pslab->alloc(536);
            int a0_1 = pslab->alloc(52464);
            pslab->free(52464, a0_1);
            pslab->free(536, a0_0);
        }
        SUPstream::Sync();
        int a1_0;
        int a1_1;
        if (SUPstream::node_manager.rank == 1)
        {
            a1_0 = pslab->alloc(536);
            a1_1 = pslab->alloc(52464);
        }
        SUPstream::Sync();
        if (SUPstream::node_manager.rank == 0)
        {
            pslab->free(52464, a1_1);
            pslab->free(536, a1_0);
        }
        SUPstream::Sync();
        int a2_0, a2_1;
        if (SUPstream::node_manager.rank == 1)
        {
            a2_0 = pslab->alloc(536);
            a2_1 = pslab->alloc(52464);
        }
        if (SUPstream::node_manager.rank == 0)
        {
            pslab->free(52464, a2_1);
            pslab->free(536, a2_0);
        }
        SUPstream::Sync();
        if (SUPstream::node_manager.rank == 0)
        {
            int a3_0 = pslab->alloc(536);
            int a3_1 = pslab->alloc(52464);
        }
        SUPstream::Sync();
        if (SUPstream::node_manager.rank == 0)
        {
            int a3_2 = pslab->alloc(536);
            int a3_3 = pslab->alloc(52464);
            pslab->free(52464, a3_2);
            int a3_4 = pslab->alloc(536);
            int a3_5 = pslab->alloc(52464);
        }
    }*/
#include "createMesh.H"
#include "createControl.H"
#include "createTimeControls.H"
#include "initContinuityErrs.H"
#include "createFields.H"
#include "createFieldRefs.H"

    turbulence->validate();

    if (!LTS)
    {
#include "compressibleCourantNo.H"
#include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << "\nStarting time loop\n"
         << endl;

    while (runTime.run())
    {
#include "readTimeControls.H"

        if (LTS)
        {
#include "setRDeltaT.H"
        }
        else
        {
#include "compressibleCourantNo.H"
#include "setDeltaT.H"
        }

        runTime++;

        Info << "Time = " << runTime.timeName() << nl << endl;

#include "rhoEqn.H"

        while (pimple.loop())
        {
#include "UEqn.H"
#include "YEqn.H"
#include "EEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                if (pimple.consistent())
                {
#include "pcEqn.H"
                }
                else
                {
#include "pEqn.H"
                }
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        rho = thermo.rho();

        runTime.write();

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << endl;
    }

    Info << "End\n"
         << endl;

    return 0;
}

// ************************************************************************* //
