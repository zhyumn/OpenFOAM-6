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

#include "makeChemistryTabulationMethodsPT.H"

#include "thermoPhysicsTypes.H"

#include "psiReactionThermo.H"
#include "rhoReactionThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Chemistry solvers based on sensibleEnthalpy
    makeChemistryTabulationMethodsPT(psiReactionThermo, constGasHThermoPhysics);
    makeChemistryTabulationMethodsPT(psiReactionThermo, gasHThermoPhysics);
    makeChemistryTabulationMethodsPT
    (
        psiReactionThermo,
        constIncompressibleGasHThermoPhysics
    );
    makeChemistryTabulationMethodsPT
    (
        psiReactionThermo,
        incompressibleGasHThermoPhysics
    );
    makeChemistryTabulationMethodsPT(psiReactionThermo, icoPoly8HThermoPhysics);
    makeChemistryTabulationMethodsPT(psiReactionThermo, constFluidHThermoPhysics);
    makeChemistryTabulationMethodsPT
    (
        psiReactionThermo,
        constAdiabaticFluidHThermoPhysics
    );
    makeChemistryTabulationMethodsPT(psiReactionThermo, constHThermoPhysics);

    makeChemistryTabulationMethodsPT(rhoReactionThermo, constGasHThermoPhysics);
    makeChemistryTabulationMethodsPT(rhoReactionThermo, gasHThermoPhysics);
    makeChemistryTabulationMethodsPT
    (
        rhoReactionThermo,
        constIncompressibleGasHThermoPhysics
    );
    makeChemistryTabulationMethodsPT
    (
        rhoReactionThermo,
        incompressibleGasHThermoPhysics
    );
    makeChemistryTabulationMethodsPT(rhoReactionThermo, icoPoly8HThermoPhysics);
    makeChemistryTabulationMethodsPT(rhoReactionThermo, constFluidHThermoPhysics);
    makeChemistryTabulationMethodsPT
    (
        rhoReactionThermo,
        constAdiabaticFluidHThermoPhysics
    );
    makeChemistryTabulationMethodsPT(rhoReactionThermo, constHThermoPhysics);

    // Chemistry solvers based on sensibleInternalEnergy

    makeChemistryTabulationMethodsPT(psiReactionThermo, constGasEThermoPhysics);
    makeChemistryTabulationMethodsPT(psiReactionThermo, gasEThermoPhysics);
    makeChemistryTabulationMethodsPT
    (
        psiReactionThermo,
        constIncompressibleGasEThermoPhysics
    );
    makeChemistryTabulationMethodsPT
    (
        psiReactionThermo,
        incompressibleGasEThermoPhysics
    );
    makeChemistryTabulationMethodsPT(psiReactionThermo, icoPoly8EThermoPhysics);
    makeChemistryTabulationMethodsPT(psiReactionThermo, constFluidEThermoPhysics);
    makeChemistryTabulationMethodsPT
    (
        psiReactionThermo,
        constAdiabaticFluidEThermoPhysics
    );
    makeChemistryTabulationMethodsPT(psiReactionThermo, constEThermoPhysics);

    makeChemistryTabulationMethodsPT(rhoReactionThermo, constGasEThermoPhysics);
    makeChemistryTabulationMethodsPT(rhoReactionThermo, gasEThermoPhysics);
    makeChemistryTabulationMethodsPT
    (
        rhoReactionThermo,
        constIncompressibleGasEThermoPhysics
    );
    makeChemistryTabulationMethodsPT
    (
        rhoReactionThermo,
        incompressibleGasEThermoPhysics
    );
    makeChemistryTabulationMethodsPT(rhoReactionThermo, icoPoly8EThermoPhysics);
    makeChemistryTabulationMethodsPT(rhoReactionThermo, constFluidEThermoPhysics);
    makeChemistryTabulationMethodsPT
    (
        rhoReactionThermo,
        constAdiabaticFluidEThermoPhysics
    );
    makeChemistryTabulationMethodsPT(rhoReactionThermo, constEThermoPhysics);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
