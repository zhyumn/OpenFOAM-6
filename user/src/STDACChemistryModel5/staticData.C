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

//#include "parISATleaf.H"
//#include "SVD.H"
//#include "IOstream.H"
#include "staticData.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    namespace chemistryTabulationMethodSs
    {
        bool static_data::mechRedActive = true;
        label static_data::nAdditionalEqns_ = 1;
        bool static_data::variableTimeStep = true;
        label static_data::nEqns = 0;
        label static_data::completeSpaceSize_ = 0;
        scalar static_data::tolerance_ = 0.0;
        label static_data::idT_ = 0;
        label static_data::idp_ = 0;
        label static_data::iddeltaT_ = 0;
    }
}
// ************************************************************************* //
