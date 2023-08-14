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
#include "chemistryTabulationMethodPT.H"
//#include "parallelISAT_chem.H"
#include "ISATPT.H"
#include "makeChemistryTabulationMethodsPT.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#define makeISATMethodPT(Comp, Thermo)                        \
                                                                               \
    typedef chemistryTabulationMethodsPT::ISAT<Comp, Thermo> ISAT##Comp##Thermo;     \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        ISAT##Comp##Thermo,                                       \
        ("parallelISAT<" + word(Comp::typeName_())                                      \
      + "," + Thermo::typeName() + ">").c_str(),                               \
        0                                                                      \
    );                                                                         \
                                                                               \
    chemistryTabulationMethodPT<Comp, Thermo>::                                  \
        adddictionaryConstructorToTable<ISAT##Comp##Thermo>                      \
        add##chemistryTabulationMethodsPT##ISAT##Comp##Thermo##ConstructorToTable_;
namespace Foam
{
    makeISATMethodPT
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
