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

#ifndef makeChemistryTabulationMethodsPT_H
#define makeChemistryTabulationMethodsPT_H

#include "chemistryTabulationMethodPT.H"
#include "noChemistryTabulationPT.H"
#include "ISATPT.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeChemistryTabulationMethodPT(SS, Comp, Thermo)                        \
                                                                               \
    typedef chemistryTabulationMethodsPT::SS<Comp, Thermo> SS##Comp##Thermo;     \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        SS##Comp##Thermo,                                                      \
        (#SS"<" + word(Comp::typeName_())                                      \
      + "," + Thermo::typeName() + ">").c_str(),                               \
        0                                                                      \
    );                                                                         \
                                                                               \
    chemistryTabulationMethodPT<Comp, Thermo>::                                  \
        adddictionaryConstructorToTable<SS##Comp##Thermo>                      \
        add##chemistryTabulationMethodsPT##SS##Comp##Thermo##ConstructorToTable_;


#define makeChemistryTabulationMethodsPT(CompChemModel, Thermo)                  \
                                                                               \
    typedef chemistryTabulationMethodPT<CompChemModel, Thermo>                   \
        chemistryTabulationMethodPT##CompChemModel##Thermo;                      \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        chemistryTabulationMethodPT##CompChemModel##Thermo,                      \
        (word(chemistryTabulationMethodPT##CompChemModel##Thermo::typeName_()) + \
        '<' + word(CompChemModel::typeName_()) + "," + Thermo::typeName() + '>'\
        ).c_str(),                                                             \
        0                                                                      \
    );                                                                         \
                                                                               \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        chemistryTabulationMethodPT##CompChemModel##Thermo,                      \
        dictionary                                                             \
    );                                                                         \
                                                                               \
    makeChemistryTabulationMethodPT                                              \
    (                                                                          \
        none,                                                                  \
        CompChemModel,                                                         \
        Thermo                                                                 \
    );                                                                         \
                                                                               \
    makeChemistryTabulationMethodPT                                              \
    (                                                                          \
        ISAT,                                                                  \
        CompChemModel,                                                         \
        Thermo                                                                 \
    );


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
