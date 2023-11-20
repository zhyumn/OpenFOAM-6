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

#include "makeChemistryReductionMethodSs.H"

#include "thermoPhysicsTypes.H"

#include "psiReactionThermo.H"
#include "rhoReactionThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Chemistry solvers based on sensibleEnthalpy
    makeChemistryReductionMethodSs(psiReactionThermo, constGasHThermoPhysics);
    makeChemistryReductionMethodSs(psiReactionThermo, gasHThermoPhysics);
    makeChemistryReductionMethodSs
    (
        psiReactionThermo,
        constIncompressibleGasHThermoPhysics
    );
    makeChemistryReductionMethodSs
    (
        psiReactionThermo,
        incompressibleGasHThermoPhysics
    );
    makeChemistryReductionMethodSs(psiReactionThermo, icoPoly8HThermoPhysics);
    makeChemistryReductionMethodSs(psiReactionThermo, constFluidHThermoPhysics);
    makeChemistryReductionMethodSs
    (
        psiReactionThermo,
        constAdiabaticFluidHThermoPhysics
    );
    makeChemistryReductionMethodSs(psiReactionThermo, constHThermoPhysics);

    makeChemistryReductionMethodSs(rhoReactionThermo, constGasHThermoPhysics);
    makeChemistryReductionMethodSs(rhoReactionThermo, gasHThermoPhysics);
    makeChemistryReductionMethodSs
    (
        rhoReactionThermo,
        constIncompressibleGasHThermoPhysics
    );
    makeChemistryReductionMethodSs
    (
        rhoReactionThermo,
        incompressibleGasHThermoPhysics
    );
    makeChemistryReductionMethodSs(rhoReactionThermo, icoPoly8HThermoPhysics);
    makeChemistryReductionMethodSs(rhoReactionThermo, constFluidHThermoPhysics);
    makeChemistryReductionMethodSs
    (
        rhoReactionThermo,
        constAdiabaticFluidHThermoPhysics
    );
    makeChemistryReductionMethodSs(rhoReactionThermo, constHThermoPhysics);


    // Chemistry solvers based on sensibleInternalEnergy
    makeChemistryReductionMethodSs(psiReactionThermo, constGasEThermoPhysics);
    makeChemistryReductionMethodSs(psiReactionThermo, gasEThermoPhysics);
    makeChemistryReductionMethodSs
    (
        psiReactionThermo,
        constIncompressibleGasEThermoPhysics
    );
    makeChemistryReductionMethodSs
    (
        psiReactionThermo,
        incompressibleGasEThermoPhysics
    );
    makeChemistryReductionMethodSs(psiReactionThermo, icoPoly8EThermoPhysics);
    makeChemistryReductionMethodSs(psiReactionThermo, constFluidEThermoPhysics);
    makeChemistryReductionMethodSs
    (
        psiReactionThermo,
        constAdiabaticFluidEThermoPhysics
    );
    makeChemistryReductionMethodSs(psiReactionThermo, constEThermoPhysics);

    makeChemistryReductionMethodSs(rhoReactionThermo, constGasEThermoPhysics);
    makeChemistryReductionMethodSs(rhoReactionThermo, gasEThermoPhysics);
    makeChemistryReductionMethodSs
    (
        rhoReactionThermo,
        constIncompressibleGasEThermoPhysics
    );
    makeChemistryReductionMethodSs
    (
        rhoReactionThermo,
        incompressibleGasEThermoPhysics
    );
    makeChemistryReductionMethodSs(rhoReactionThermo, icoPoly8EThermoPhysics);
    makeChemistryReductionMethodSs(rhoReactionThermo, constFluidEThermoPhysics);
    makeChemistryReductionMethodSs
    (
        rhoReactionThermo,
        constAdiabaticFluidEThermoPhysics
    );
    makeChemistryReductionMethodSs(rhoReactionThermo, constEThermoPhysics);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
