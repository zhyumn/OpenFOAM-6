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

#ifndef makeChemistryReductionMethodsPT_H
#define makeChemistryReductionMethodsPT_H

#include "chemistryReductionMethodPT.H"

#include "noChemistryReductionPT.H"
#include "DACS.H"
#include "DRGS.H"
#include "DRGEPS.H"
#include "EFAS.H"
#include "PFAS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeChemistryReductionMethodPT(SS, Comp, Thermo)                         \
                                                                               \
    typedef chemistryReductionMethodsPT::SS<Comp, Thermo> SS##Comp##Thermo;      \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        SS##Comp##Thermo,                                                      \
        (#SS"<" + word(Comp::typeName_())                                      \
      + "," + Thermo::typeName() + ">").c_str(),                               \
        0                                                                      \
    );                                                                         \
                                                                               \
    chemistryReductionMethodPT<Comp, Thermo>::                             \
        adddictionaryConstructorToTable<SS##Comp##Thermo>                      \
        add##chemistryReductionMethodsPT##SS##Comp##Thermo##ConstructorToTable_;


#define makeChemistryReductionMethodsPT(CompChemModel, Thermo)                   \
                                                                               \
    typedef chemistryReductionMethodPT<CompChemModel, Thermo>                    \
        chemistryReductionMethodPT##CompChemModel##Thermo;                       \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        chemistryReductionMethodPT##CompChemModel##Thermo,                       \
        (word(chemistryReductionMethodPT##CompChemModel##Thermo::typeName_()) +  \
        '<' + word(CompChemModel::typeName_()) + "," + Thermo::typeName() + '>'\
        ).c_str(),                                                             \
        0                                                                      \
    );                                                                         \
                                                                               \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        chemistryReductionMethodPT##CompChemModel##Thermo,                       \
        dictionary                                                             \
    );                                                                         \
                                                                               \
    makeChemistryReductionMethodPT                                               \
    (                                                                          \
        none,                                                                  \
        CompChemModel,                                                         \
        Thermo                                                                 \
    );                                                                         \
                                                                               \
    makeChemistryReductionMethodPT                                               \
    (                                                                          \
        DAC,                                                                   \
        CompChemModel,                                                         \
        Thermo                                                                 \
    );                                                                         \
                                                                               \
    makeChemistryReductionMethodPT                                               \
    (                                                                          \
        DRG,                                                                   \
        CompChemModel,                                                         \
        Thermo                                                                 \
    );                                                                         \
                                                                               \
    makeChemistryReductionMethodPT                                               \
    (                                                                          \
        DRGEP,                                                                 \
        CompChemModel,                                                         \
        Thermo                                                                 \
    );                                                                         \
                                                                               \
    makeChemistryReductionMethodPT                                               \
    (                                                                          \
        EFA,                                                                   \
        CompChemModel,                                                         \
        Thermo                                                                 \
    );                                                                         \
                                                                               \
    makeChemistryReductionMethodPT                                               \
    (                                                                          \
        PFA,                                                                   \
        CompChemModel,                                                         \
        Thermo                                                                 \
    );                                                                         \


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
