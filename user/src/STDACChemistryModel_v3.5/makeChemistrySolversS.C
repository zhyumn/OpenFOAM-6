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

#include "makeChemistrySolverTypesS.H"

#include "thermoPhysicsTypes.H"
#include "psiReactionThermo.H"
#include "rhoReactionThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Chemistry solvers based on sensibleEnthalpy
    makeChemistrySolverTypesS(psiReactionThermo, constGasHThermoPhysics);
    makeChemistrySolverTypesS(psiReactionThermo, gasHThermoPhysics);
    makeChemistrySolverTypesS
    (
        psiReactionThermo,
        constIncompressibleGasHThermoPhysics
    );
    makeChemistrySolverTypesS
    (
        psiReactionThermo,
        incompressibleGasHThermoPhysics
    );
    makeChemistrySolverTypesS(psiReactionThermo, icoPoly8HThermoPhysics);
    makeChemistrySolverTypesS(psiReactionThermo, constFluidHThermoPhysics);
    makeChemistrySolverTypesS
    (
        psiReactionThermo,
        constAdiabaticFluidHThermoPhysics
    );
    makeChemistrySolverTypesS(psiReactionThermo, constHThermoPhysics);

    makeChemistrySolverTypesS(rhoReactionThermo, constGasHThermoPhysics);
    makeChemistrySolverTypesS(rhoReactionThermo, gasHThermoPhysics);
    makeChemistrySolverTypesS
    (
        rhoReactionThermo,
        constIncompressibleGasHThermoPhysics
    );
    makeChemistrySolverTypesS
    (
        rhoReactionThermo,
        incompressibleGasHThermoPhysics
    );
    makeChemistrySolverTypesS(rhoReactionThermo, icoPoly8HThermoPhysics);
    makeChemistrySolverTypesS(rhoReactionThermo, constFluidHThermoPhysics);
    makeChemistrySolverTypesS
    (
        rhoReactionThermo,
        constAdiabaticFluidHThermoPhysics
    );
    makeChemistrySolverTypesS(rhoReactionThermo, constHThermoPhysics);

    // Chemistry solvers based on sensibleInternalEnergy
    makeChemistrySolverTypesS(psiReactionThermo, constGasEThermoPhysics);
    makeChemistrySolverTypesS(psiReactionThermo, gasEThermoPhysics);
    makeChemistrySolverTypesS
    (
        psiReactionThermo,
        constIncompressibleGasEThermoPhysics
    );
    makeChemistrySolverTypesS
    (
        psiReactionThermo,
        incompressibleGasEThermoPhysics
    );
    makeChemistrySolverTypesS(psiReactionThermo, icoPoly8EThermoPhysics);
    makeChemistrySolverTypesS(psiReactionThermo, constFluidEThermoPhysics);
    makeChemistrySolverTypesS
    (
        psiReactionThermo,
        constAdiabaticFluidEThermoPhysics
    );
    makeChemistrySolverTypesS(psiReactionThermo, constEThermoPhysics);

    makeChemistrySolverTypesS(rhoReactionThermo, constGasEThermoPhysics);
    makeChemistrySolverTypesS(rhoReactionThermo, gasEThermoPhysics);
    makeChemistrySolverTypesS
    (
        rhoReactionThermo,
        constIncompressibleGasEThermoPhysics
    );
    makeChemistrySolverTypesS
    (
        rhoReactionThermo,
        incompressibleGasEThermoPhysics
    );
    makeChemistrySolverTypesS(rhoReactionThermo, icoPoly8EThermoPhysics);
    makeChemistrySolverTypesS(rhoReactionThermo, constFluidEThermoPhysics);
    makeChemistrySolverTypesS
    (
        rhoReactionThermo,
        constAdiabaticFluidEThermoPhysics
    );
    makeChemistrySolverTypesS(rhoReactionThermo, constEThermoPhysics);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
