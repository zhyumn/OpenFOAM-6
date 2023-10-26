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

#include "makeChemistryTabulationMethodSs.H"

#include "thermoPhysicsTypes.H"

#include "psiReactionThermo.H"
#include "rhoReactionThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Chemistry solvers based on sensibleEnthalpy
    makeChemistryTabulationMethodSs(psiReactionThermo, constGasHThermoPhysics);
    makeChemistryTabulationMethodSs(psiReactionThermo, gasHThermoPhysics);
    makeChemistryTabulationMethodSs
    (
        psiReactionThermo,
        constIncompressibleGasHThermoPhysics
    );
    makeChemistryTabulationMethodSs
    (
        psiReactionThermo,
        incompressibleGasHThermoPhysics
    );
    makeChemistryTabulationMethodSs(psiReactionThermo, icoPoly8HThermoPhysics);
    makeChemistryTabulationMethodSs(psiReactionThermo, constFluidHThermoPhysics);
    makeChemistryTabulationMethodSs
    (
        psiReactionThermo,
        constAdiabaticFluidHThermoPhysics
    );
    makeChemistryTabulationMethodSs(psiReactionThermo, constHThermoPhysics);

    makeChemistryTabulationMethodSs(rhoReactionThermo, constGasHThermoPhysics);
    makeChemistryTabulationMethodSs(rhoReactionThermo, gasHThermoPhysics);
    makeChemistryTabulationMethodSs
    (
        rhoReactionThermo,
        constIncompressibleGasHThermoPhysics
    );
    makeChemistryTabulationMethodSs
    (
        rhoReactionThermo,
        incompressibleGasHThermoPhysics
    );
    makeChemistryTabulationMethodSs(rhoReactionThermo, icoPoly8HThermoPhysics);
    makeChemistryTabulationMethodSs(rhoReactionThermo, constFluidHThermoPhysics);
    makeChemistryTabulationMethodSs
    (
        rhoReactionThermo,
        constAdiabaticFluidHThermoPhysics
    );
    makeChemistryTabulationMethodSs(rhoReactionThermo, constHThermoPhysics);

    // Chemistry solvers based on sensibleInternalEnergy

    makeChemistryTabulationMethodSs(psiReactionThermo, constGasEThermoPhysics);
    makeChemistryTabulationMethodSs(psiReactionThermo, gasEThermoPhysics);
    makeChemistryTabulationMethodSs
    (
        psiReactionThermo,
        constIncompressibleGasEThermoPhysics
    );
    makeChemistryTabulationMethodSs
    (
        psiReactionThermo,
        incompressibleGasEThermoPhysics
    );
    makeChemistryTabulationMethodSs(psiReactionThermo, icoPoly8EThermoPhysics);
    makeChemistryTabulationMethodSs(psiReactionThermo, constFluidEThermoPhysics);
    makeChemistryTabulationMethodSs
    (
        psiReactionThermo,
        constAdiabaticFluidEThermoPhysics
    );
    makeChemistryTabulationMethodSs(psiReactionThermo, constEThermoPhysics);

    makeChemistryTabulationMethodSs(rhoReactionThermo, constGasEThermoPhysics);
    makeChemistryTabulationMethodSs(rhoReactionThermo, gasEThermoPhysics);
    makeChemistryTabulationMethodSs
    (
        rhoReactionThermo,
        constIncompressibleGasEThermoPhysics
    );
    makeChemistryTabulationMethodSs
    (
        rhoReactionThermo,
        incompressibleGasEThermoPhysics
    );
    makeChemistryTabulationMethodSs(rhoReactionThermo, icoPoly8EThermoPhysics);
    makeChemistryTabulationMethodSs(rhoReactionThermo, constFluidEThermoPhysics);
    makeChemistryTabulationMethodSs
    (
        rhoReactionThermo,
        constAdiabaticFluidEThermoPhysics
    );
    makeChemistryTabulationMethodSs(rhoReactionThermo, constEThermoPhysics);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
