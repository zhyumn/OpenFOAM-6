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
#line 25 "/panfs/jay/groups/25/suo-yang/srini237/OpenFOAM/OpenFOAM-6/user/tutorials/Nike/TemporalMixingLayer_2d_12-1/0/N2.#codeStream"
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
    void codeStream_80df80d53b21306ff4879a258001646cb3b8705b
    (
        Ostream& os,
        const dictionary& dict
    )
    {
//{{{ begin code
        #line 39 "/panfs/jay/groups/25/suo-yang/srini237/OpenFOAM/OpenFOAM-6/user/tutorials/Nike/TemporalMixingLayer_2d_12-1/0/N2.#codeStream"
const IOdictionary & d = static_cast<const IOdictionary&>(dict);
            const fvMesh& mesh = refCast<const fvMesh>(d.db());
            const scalar PI = 3.1415926535897932384626;
            
            double dx=2.159e-5;
            double thickness =dx/29.16;//5.901;
            double Ntop=0.001;//0.9499;
            double Nbottom=0.999;//0.0001;

            double WMC=86.1754;
            double WMN=28.0134;

            scalarField N2(mesh.nCells());
            forAll(N2, i)
            {
                const scalar z = mesh.C()[i][2];
                double a = 0.5*(1+erf(sqrt(PI)*z/thickness));
                N2[i] =0.5*(Ntop+Nbottom)+0.5*(Ntop-Nbottom)*erf(sqrt(PI)*z/thickness);
                double X=(Ntop*a+Nbottom*(1-a));
                //N2[i] = WMN*X/(WMN*X+WMC*(1-X));
                
            }
            N2.writeEntry( "", os);
//}}} end code
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

