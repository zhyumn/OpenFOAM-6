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

\*---------------------------------------------------------------------------*/

#include "makeChemistrySolverTypesPT.H"

#include "thermoPhysicsTypes.H"
#include "psiReactionThermo.H"
#include "rhoReactionThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Chemistry solvers based on sensibleEnthalpy
    makeChemistrySolverTypesPT(psiReactionThermo, constGasHThermoPhysics);
    makeChemistrySolverTypesPT(psiReactionThermo, gasHThermoPhysics);
    makeChemistrySolverTypesPT
    (
        psiReactionThermo,
        constIncompressibleGasHThermoPhysics
    );
    makeChemistrySolverTypesPT
    (
        psiReactionThermo,
        incompressibleGasHThermoPhysics
    );
    makeChemistrySolverTypesPT(psiReactionThermo, icoPoly8HThermoPhysics);
    makeChemistrySolverTypesPT(psiReactionThermo, constFluidHThermoPhysics);
    makeChemistrySolverTypesPT
    (
        psiReactionThermo,
        constAdiabaticFluidHThermoPhysics
    );
    makeChemistrySolverTypesPT(psiReactionThermo, constHThermoPhysics);

    makeChemistrySolverTypesPT(rhoReactionThermo, constGasHThermoPhysics);
    makeChemistrySolverTypesPT(rhoReactionThermo, gasHThermoPhysics);
    makeChemistrySolverTypesPT
    (
        rhoReactionThermo,
        constIncompressibleGasHThermoPhysics
    );
    makeChemistrySolverTypesPT
    (
        rhoReactionThermo,
        incompressibleGasHThermoPhysics
    );
    makeChemistrySolverTypesPT(rhoReactionThermo, icoPoly8HThermoPhysics);
    makeChemistrySolverTypesPT(rhoReactionThermo, constFluidHThermoPhysics);
    makeChemistrySolverTypesPT
    (
        rhoReactionThermo,
        constAdiabaticFluidHThermoPhysics
    );
    makeChemistrySolverTypesPT(rhoReactionThermo, constHThermoPhysics);

    // Chemistry solvers based on sensibleInternalEnergy
    makeChemistrySolverTypesPT(psiReactionThermo, constGasEThermoPhysics);
    makeChemistrySolverTypesPT(psiReactionThermo, gasEThermoPhysics);
    makeChemistrySolverTypesPT
    (
        psiReactionThermo,
        constIncompressibleGasEThermoPhysics
    );
    makeChemistrySolverTypesPT
    (
        psiReactionThermo,
        incompressibleGasEThermoPhysics
    );
    makeChemistrySolverTypesPT(psiReactionThermo, icoPoly8EThermoPhysics);
    makeChemistrySolverTypesPT(psiReactionThermo, constFluidEThermoPhysics);
    makeChemistrySolverTypesPT
    (
        psiReactionThermo,
        constAdiabaticFluidEThermoPhysics
    );
    makeChemistrySolverTypesPT(psiReactionThermo, constEThermoPhysics);

    makeChemistrySolverTypesPT(rhoReactionThermo, constGasEThermoPhysics);
    makeChemistrySolverTypesPT(rhoReactionThermo, gasEThermoPhysics);
    makeChemistrySolverTypesPT
    (
        rhoReactionThermo,
        constIncompressibleGasEThermoPhysics
    );
    makeChemistrySolverTypesPT
    (
        rhoReactionThermo,
        incompressibleGasEThermoPhysics
    );
    makeChemistrySolverTypesPT(rhoReactionThermo, icoPoly8EThermoPhysics);
    makeChemistrySolverTypesPT(rhoReactionThermo, constFluidEThermoPhysics);
    makeChemistrySolverTypesPT
    (
        rhoReactionThermo,
        constAdiabaticFluidEThermoPhysics
    );
    makeChemistrySolverTypesPT(rhoReactionThermo, constEThermoPhysics);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
