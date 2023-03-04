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

#ifndef makeParallelIAST_H
#define makeParallelIAST_H

#include "thermoPhysicsTypes.H"
#include "psiReactionThermo.H"
#include "chemistryTabulationMethod.H"
#include "parallelISAT_chem.H"
#include "makeChemistryTabulationMethods.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#define makeParallelISATMethod(Comp, Thermo)                        \
                                                                               \
    typedef chemistryTabulationMethods::parallelISAT_chem<Comp, Thermo> parallelISAT_chem##Comp##Thermo;     \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        parallelISAT_chem##Comp##Thermo,                                       \
        ("parallelISAT<" + word(Comp::typeName_())                                      \
      + "," + Thermo::typeName() + ">").c_str(),                               \
        0                                                                      \
    );                                                                         \
                                                                               \
    chemistryTabulationMethod<Comp, Thermo>::                                  \
        adddictionaryConstructorToTable<parallelISAT_chem##Comp##Thermo>                      \
        add##chemistryTabulationMethods##parallelISAT_chem##Comp##Thermo##ConstructorToTable_;
namespace Foam
{

    makeParallelISATMethod
    (
        psiReactionThermo,
        gasHThermoPhysics
    );
/*
    makeChemistryTabulationMethod
    (
        parallelISAT_chem,
        psiReactionThermo,
        gasHThermoPhysics
    );
    */
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
