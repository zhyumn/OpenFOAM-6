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
#include "STDACChemistryModel.H"
#include "thermoPhysicsTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Make base types
    //makeChemistryModel(psiReactionThermo);
    //makeChemistryModel(rhoReactionThermo);

    // Chemistry moldels based on sensibleEnthalpy

    makeChemistryModelType(
        STDACChemistryModel,
        psiReactionThermo,
        constGasHThermoPhysics);

    makeChemistryModelType(
        STDACChemistryModel,
        psiReactionThermo,
        gasHThermoPhysics);

    makeChemistryModelType(
        STDACChemistryModel,
        psiReactionThermo,
        constIncompressibleGasHThermoPhysics);

    makeChemistryModelType(
        STDACChemistryModel,
        psiReactionThermo,
        incompressibleGasHThermoPhysics);

    makeChemistryModelType(
        STDACChemistryModel,
        psiReactionThermo,
        icoPoly8HThermoPhysics);

    makeChemistryModelType(
        STDACChemistryModel,
        psiReactionThermo,
        constFluidHThermoPhysics);

    makeChemistryModelType(
        STDACChemistryModel,
        psiReactionThermo,
        constAdiabaticFluidHThermoPhysics);

    makeChemistryModelType(
        STDACChemistryModel,
        psiReactionThermo,
        constHThermoPhysics);

    makeChemistryModelType(
        STDACChemistryModel,
        rhoReactionThermo,
        constGasHThermoPhysics);

    makeChemistryModelType(
        STDACChemistryModel,
        rhoReactionThermo,
        gasHThermoPhysics);

    makeChemistryModelType(
        STDACChemistryModel,
        rhoReactionThermo,
        constIncompressibleGasHThermoPhysics);

    makeChemistryModelType(
        STDACChemistryModel,
        rhoReactionThermo,
        incompressibleGasHThermoPhysics);

    makeChemistryModelType(
        STDACChemistryModel,
        rhoReactionThermo,
        icoPoly8HThermoPhysics);

    makeChemistryModelType(
        STDACChemistryModel,
        rhoReactionThermo,
        constFluidHThermoPhysics);

    makeChemistryModelType(
        STDACChemistryModel,
        rhoReactionThermo,
        constAdiabaticFluidHThermoPhysics);

    makeChemistryModelType(
        STDACChemistryModel,
        rhoReactionThermo,
        constHThermoPhysics);

    // Chemistry moldels based on sensibleInternalEnergy

    makeChemistryModelType(
        STDACChemistryModel,
        psiReactionThermo,
        constGasEThermoPhysics);

    makeChemistryModelType(
        STDACChemistryModel,
        psiReactionThermo,
        gasEThermoPhysics);

    makeChemistryModelType(
        STDACChemistryModel,
        psiReactionThermo,
        constIncompressibleGasEThermoPhysics);

    makeChemistryModelType(
        STDACChemistryModel,
        psiReactionThermo,
        incompressibleGasEThermoPhysics);

    makeChemistryModelType(
        STDACChemistryModel,
        psiReactionThermo,
        icoPoly8EThermoPhysics);

    makeChemistryModelType(
        STDACChemistryModel,
        psiReactionThermo,
        constFluidEThermoPhysics);

    makeChemistryModelType(
        STDACChemistryModel,
        psiReactionThermo,
        constAdiabaticFluidEThermoPhysics);

    makeChemistryModelType(
        STDACChemistryModel,
        psiReactionThermo,
        constEThermoPhysics);

    makeChemistryModelType(
        STDACChemistryModel,
        rhoReactionThermo,
        constGasEThermoPhysics);

    makeChemistryModelType(
        STDACChemistryModel,
        rhoReactionThermo,
        gasEThermoPhysics);

    makeChemistryModelType(
        STDACChemistryModel,
        rhoReactionThermo,
        constIncompressibleGasEThermoPhysics);

    makeChemistryModelType(
        STDACChemistryModel,
        rhoReactionThermo,
        incompressibleGasEThermoPhysics);

    makeChemistryModelType(
        STDACChemistryModel,
        rhoReactionThermo,
        icoPoly8EThermoPhysics);

    makeChemistryModelType(
        STDACChemistryModel,
        rhoReactionThermo,
        constFluidEThermoPhysics);

    makeChemistryModelType(
        STDACChemistryModel,
        rhoReactionThermo,
        constAdiabaticFluidEThermoPhysics);

    makeChemistryModelType(
        STDACChemistryModel,
        rhoReactionThermo,
        constEThermoPhysics);
}

// ************************************************************************* //
