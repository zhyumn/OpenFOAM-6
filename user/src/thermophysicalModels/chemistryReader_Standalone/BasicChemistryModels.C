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

//#include "makeChemistryModel.H"

#include "psiReactionThermo.H"
#include "rhoReactionThermo.H"

#include "StandardChemistryModelS.H"
#include "TDACChemistryModel.H"
#include "thermoPhysicsTypes.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#define makeChemistryModel(Comp)                                               \
                                                                               \
    typedef BasicChemistryModel<Comp> BasicChemistryModel##Comp;               \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        BasicChemistryModel##Comp,                                             \
        "BasicChemistryModel<"#Comp">",                                        \
        0                                                                      \
    );                                                                         \
                                                                               \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        BasicChemistryModel##Comp,                                             \
        thermo                                                                 \
    );


#define makeChemistryModelType(SS, Comp, Thermo)                               \
                                                                               \
    typedef SS<Comp, Thermo> SS##Comp##Thermo;                                 \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        SS##Comp##Thermo,                                                      \
        (#SS"<"#Comp"," + Thermo::typeName() + ">").c_str(),                   \
        0                                                                      \
    );
namespace Foam
{
    // Make base types
    //makeChemistryModel(psiReactionThermo);
    //makeChemistryModel(rhoReactionThermo);

    // Chemistry moldels based on sensibleEnthalpy


    makeChemistryModelType
    (
        StandardChemistryModelS,
        psiReactionThermo,
        gasHThermoPhysics
    );

 
}

// ************************************************************************* //
