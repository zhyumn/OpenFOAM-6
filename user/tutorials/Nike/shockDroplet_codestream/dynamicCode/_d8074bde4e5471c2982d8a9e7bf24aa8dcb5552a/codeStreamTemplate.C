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
#line 25 "/scratch.nike/srini237/OpenFOAM/OpenFOAM-6/user/tutorials/Nike/shockDroplet_codestream/0/N2.#codeStream"
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
    void codeStream_d8074bde4e5471c2982d8a9e7bf24aa8dcb5552a
    (
        Ostream& os,
        const dictionary& dict
    )
    {
//{{{ begin code
        #line 39 "/scratch.nike/srini237/OpenFOAM/OpenFOAM-6/user/tutorials/Nike/shockDroplet_codestream/0/N2.#codeStream"
const IOdictionary & d = static_cast<const IOdictionary&>(dict);
            const fvMesh& mesh = refCast<const fvMesh>(d.db());

            double r = 0.125;
            double x_center = 0.5;
            double y_center = 0.5;
            double thickness = 1.0/15;
            double omega=0.01;
            double Nlow=0.0001;
            double Nhigh=0.9499;

            scalarField N2(mesh.nCells(), Nhigh);
            forAll(N2, i)
            {
                const scalar x = mesh.C()[i][0];
                const scalar y = mesh.C()[i][1];
                //const scalar z = mesh.C()[i][2];
                scalar r_tmp=sqrt(sqr(y - y_center)+sqr(x-x_center));
                if (r_tmp < r-thickness/2)
                {
                    N2[i] = Nlow;
                }
                else if (r_tmp < r+thickness/2)
                {
                    N2[i] =0.5*(Nlow+Nhigh)+0.5*(-Nlow+Nhigh)*tanh((r_tmp-r)/omega);
                }
            }
            N2.writeEntry( "", os);
//}}} end code
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

