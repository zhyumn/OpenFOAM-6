/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
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

Description
    Template for use with codeStream.

\*---------------------------------------------------------------------------*/

#include "dictionary.H"
#include "Ostream.H"
#include "Pstream.H"
#include "unitConversion.H"

//{{{ begin codeInclude
#line 25 "/panfs/jay/groups/25/suo-yang/srini237/OpenFOAM/OpenFOAM-6/user/tutorials/Nike/TML_reacting_Turb/0/CO2.#codeStream"
#include "fvCFD.H"
//}}} end codeInclude

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    void codeStream_7e566b90866c2428f8ed2dce995f08ee2da1d685
    (
        Ostream& os,
        const dictionary& dict
    )
    {
//{{{ begin code
        #line 39 "/panfs/jay/groups/25/suo-yang/srini237/OpenFOAM/OpenFOAM-6/user/tutorials/Nike/TML_reacting_Turb/0/CO2.#codeStream"
const IOdictionary & d = static_cast<const IOdictionary&>(dict);
            const fvMesh& mesh = refCast<const fvMesh>(d.db());
            const scalar PI = 3.1415926535897932384626;

            double thickness = 5.95e-7;
            double Htop=0.0623;//0.05;
            double Hbottom=0.0001;//0.0001;

            scalarField CO2(mesh.nCells());
            forAll(CO2, i)
            {
                const scalar z = mesh.C()[i][2];
                
                    CO2[i] =0.5*(Htop+Hbottom)+0.5*(Htop-Hbottom)*erf(sqrt(PI)*z/thickness);
                
            }
            CO2.writeEntry( "", os);
//}}} end code
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

