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
#include "chemistryTabulationMethodS.H"
#include "ISAT_manager.H"
#include "makeChemistryTabulationMethodSs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#define makeParallelISATMethod(Comp, Thermo)                        \
                                                                               \
    typedef chemistryTabulationMethodSs::ISAT_manager<Comp, Thermo> ISAT_manager##Comp##Thermo;     \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        ISAT_manager##Comp##Thermo,                                       \
        ("parallelISAT<" + word(Comp::typeName_())                                      \
      + "," + Thermo::typeName() + ">").c_str(),                               \
        0                                                                      \
    );                                                                         \
                                                                               \
    chemistryTabulationMethodS<Comp, Thermo>::                                  \
        adddictionaryConstructorToTable<ISAT_manager##Comp##Thermo>                      \
        add##chemistryTabulationMethodSs##ISAT_manager##Comp##Thermo##ConstructorToTable_;
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
        ISAT_manager,
        psiReactionThermo,
        gasHThermoPhysics
    );
    */
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
