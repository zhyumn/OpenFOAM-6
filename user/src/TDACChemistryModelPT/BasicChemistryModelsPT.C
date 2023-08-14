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

InClass
    Foam::psiChemistryModel

Description
    Creates chemistry model instances templated on the type of thermodynamics

\*---------------------------------------------------------------------------*/

#include "makeChemistryModel.H"

#include "psiReactionThermo.H"
#include "rhoReactionThermo.H"

//#include "StandardChemistryModel.H"
#include "TDACChemistryModelPT.H"
#include "thermoPhysicsTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Make base types
    //makeChemistryModel(psiReactionThermo);
    //makeChemistryModel(rhoReactionThermo);

    // Chemistry moldels based on sensibleEnthalpy

    makeChemistryModelType(
        TDACChemistryModelPT,
        psiReactionThermo,
        constGasHThermoPhysics);

    makeChemistryModelType(
        TDACChemistryModelPT,
        psiReactionThermo,
        gasHThermoPhysics);

    makeChemistryModelType(
        TDACChemistryModelPT,
        psiReactionThermo,
        constIncompressibleGasHThermoPhysics);

    makeChemistryModelType(
        TDACChemistryModelPT,
        psiReactionThermo,
        incompressibleGasHThermoPhysics);

    makeChemistryModelType(
        TDACChemistryModelPT,
        psiReactionThermo,
        icoPoly8HThermoPhysics);

    makeChemistryModelType(
        TDACChemistryModelPT,
        psiReactionThermo,
        constFluidHThermoPhysics);

    makeChemistryModelType(
        TDACChemistryModelPT,
        psiReactionThermo,
        constAdiabaticFluidHThermoPhysics);

    makeChemistryModelType(
        TDACChemistryModelPT,
        psiReactionThermo,
        constHThermoPhysics);

    makeChemistryModelType(
        TDACChemistryModelPT,
        rhoReactionThermo,
        constGasHThermoPhysics);

    makeChemistryModelType(
        TDACChemistryModelPT,
        rhoReactionThermo,
        gasHThermoPhysics);

    makeChemistryModelType(
        TDACChemistryModelPT,
        rhoReactionThermo,
        constIncompressibleGasHThermoPhysics);

    makeChemistryModelType(
        TDACChemistryModelPT,
        rhoReactionThermo,
        incompressibleGasHThermoPhysics);

    makeChemistryModelType(
        TDACChemistryModelPT,
        rhoReactionThermo,
        icoPoly8HThermoPhysics);

    makeChemistryModelType(
        TDACChemistryModelPT,
        rhoReactionThermo,
        constFluidHThermoPhysics);

    makeChemistryModelType(
        TDACChemistryModelPT,
        rhoReactionThermo,
        constAdiabaticFluidHThermoPhysics);

    makeChemistryModelType(
        TDACChemistryModelPT,
        rhoReactionThermo,
        constHThermoPhysics);

    // Chemistry moldels based on sensibleInternalEnergy

    makeChemistryModelType(
        TDACChemistryModelPT,
        psiReactionThermo,
        constGasEThermoPhysics);

    makeChemistryModelType(
        TDACChemistryModelPT,
        psiReactionThermo,
        gasEThermoPhysics);

    makeChemistryModelType(
        TDACChemistryModelPT,
        psiReactionThermo,
        constIncompressibleGasEThermoPhysics);

    makeChemistryModelType(
        TDACChemistryModelPT,
        psiReactionThermo,
        incompressibleGasEThermoPhysics);

    makeChemistryModelType(
        TDACChemistryModelPT,
        psiReactionThermo,
        icoPoly8EThermoPhysics);

    makeChemistryModelType(
        TDACChemistryModelPT,
        psiReactionThermo,
        constFluidEThermoPhysics);

    makeChemistryModelType(
        TDACChemistryModelPT,
        psiReactionThermo,
        constAdiabaticFluidEThermoPhysics);

    makeChemistryModelType(
        TDACChemistryModelPT,
        psiReactionThermo,
        constEThermoPhysics);

    makeChemistryModelType(
        TDACChemistryModelPT,
        rhoReactionThermo,
        constGasEThermoPhysics);

    makeChemistryModelType(
        TDACChemistryModelPT,
        rhoReactionThermo,
        gasEThermoPhysics);

    makeChemistryModelType(
        TDACChemistryModelPT,
        rhoReactionThermo,
        constIncompressibleGasEThermoPhysics);

    makeChemistryModelType(
        TDACChemistryModelPT,
        rhoReactionThermo,
        incompressibleGasEThermoPhysics);

    makeChemistryModelType(
        TDACChemistryModelPT,
        rhoReactionThermo,
        icoPoly8EThermoPhysics);

    makeChemistryModelType(
        TDACChemistryModelPT,
        rhoReactionThermo,
        constFluidEThermoPhysics);

    makeChemistryModelType(
        TDACChemistryModelPT,
        rhoReactionThermo,
        constAdiabaticFluidEThermoPhysics);

    makeChemistryModelType(
        TDACChemistryModelPT,
        rhoReactionThermo,
        constEThermoPhysics);
}

// ************************************************************************* //
