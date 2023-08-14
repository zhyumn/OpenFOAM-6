/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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

#include "makeChemistryReductionMethodsPT.H"

#include "thermoPhysicsTypes.H"

#include "psiReactionThermo.H"
#include "rhoReactionThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Chemistry solvers based on sensibleEnthalpy
    makeChemistryReductionMethodsPT(psiReactionThermo, constGasHThermoPhysics);
    makeChemistryReductionMethodsPT(psiReactionThermo, gasHThermoPhysics);
    makeChemistryReductionMethodsPT
    (
        psiReactionThermo,
        constIncompressibleGasHThermoPhysics
    );
    makeChemistryReductionMethodsPT
    (
        psiReactionThermo,
        incompressibleGasHThermoPhysics
    );
    makeChemistryReductionMethodsPT(psiReactionThermo, icoPoly8HThermoPhysics);
    makeChemistryReductionMethodsPT(psiReactionThermo, constFluidHThermoPhysics);
    makeChemistryReductionMethodsPT
    (
        psiReactionThermo,
        constAdiabaticFluidHThermoPhysics
    );
    makeChemistryReductionMethodsPT(psiReactionThermo, constHThermoPhysics);

    makeChemistryReductionMethodsPT(rhoReactionThermo, constGasHThermoPhysics);
    makeChemistryReductionMethodsPT(rhoReactionThermo, gasHThermoPhysics);
    makeChemistryReductionMethodsPT
    (
        rhoReactionThermo,
        constIncompressibleGasHThermoPhysics
    );
    makeChemistryReductionMethodsPT
    (
        rhoReactionThermo,
        incompressibleGasHThermoPhysics
    );
    makeChemistryReductionMethodsPT(rhoReactionThermo, icoPoly8HThermoPhysics);
    makeChemistryReductionMethodsPT(rhoReactionThermo, constFluidHThermoPhysics);
    makeChemistryReductionMethodsPT
    (
        rhoReactionThermo,
        constAdiabaticFluidHThermoPhysics
    );
    makeChemistryReductionMethodsPT(rhoReactionThermo, constHThermoPhysics);


    // Chemistry solvers based on sensibleInternalEnergy
    makeChemistryReductionMethodsPT(psiReactionThermo, constGasEThermoPhysics);
    makeChemistryReductionMethodsPT(psiReactionThermo, gasEThermoPhysics);
    makeChemistryReductionMethodsPT
    (
        psiReactionThermo,
        constIncompressibleGasEThermoPhysics
    );
    makeChemistryReductionMethodsPT
    (
        psiReactionThermo,
        incompressibleGasEThermoPhysics
    );
    makeChemistryReductionMethodsPT(psiReactionThermo, icoPoly8EThermoPhysics);
    makeChemistryReductionMethodsPT(psiReactionThermo, constFluidEThermoPhysics);
    makeChemistryReductionMethodsPT
    (
        psiReactionThermo,
        constAdiabaticFluidEThermoPhysics
    );
    makeChemistryReductionMethodsPT(psiReactionThermo, constEThermoPhysics);

    makeChemistryReductionMethodsPT(rhoReactionThermo, constGasEThermoPhysics);
    makeChemistryReductionMethodsPT(rhoReactionThermo, gasEThermoPhysics);
    makeChemistryReductionMethodsPT
    (
        rhoReactionThermo,
        constIncompressibleGasEThermoPhysics
    );
    makeChemistryReductionMethodsPT
    (
        rhoReactionThermo,
        incompressibleGasEThermoPhysics
    );
    makeChemistryReductionMethodsPT(rhoReactionThermo, icoPoly8EThermoPhysics);
    makeChemistryReductionMethodsPT(rhoReactionThermo, constFluidEThermoPhysics);
    makeChemistryReductionMethodsPT
    (
        rhoReactionThermo,
        constAdiabaticFluidEThermoPhysics
    );
    makeChemistryReductionMethodsPT(rhoReactionThermo, constEThermoPhysics);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
