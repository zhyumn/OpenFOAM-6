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
#line 25 "/panfs/jay/groups/25/suo-yang/srini237/OpenFOAM/OpenFOAM-6/user/tutorials/Nike/TemporalMixingLayer_2d/0/U.#codeStream"
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
    void codeStream_a04e1528bb38d395ab167a15ec3e27e58069520e
    (
        Ostream& os,
        const dictionary& dict
    )
    {
//{{{ begin code
        #line 39 "/panfs/jay/groups/25/suo-yang/srini237/OpenFOAM/OpenFOAM-6/user/tutorials/Nike/TemporalMixingLayer_2d/0/U.#codeStream"
const IOdictionary &d = static_cast<const IOdictionary &>(dict);
            const fvMesh &mesh = refCast<const fvMesh>(d.db());
            const scalar PI = 3.1415926535897932384626;

            //double delta = 0.2*5.901;
            double lambda1 = 4.3375e-6;//delta *7.29;
            double delta = 5.95e-7;//lambda1;///7.29;
            double lambda3 = 0.6*lambda1;

            double Utop = 108.126;
            double Ubottom = -87.86;
            double F2D = 0.05;//0.05;
            double dU = mag(Utop - Ubottom);
            double Gamma3 = lambda1*dU;
            
            
            double A[4] = {1.0, 0.5, 0.35, 0.0};
            double B[4] = {1.0, 0.0, 0.0, 0.025};
            
            double alpha = 2 * PI / lambda1;
            double alphaS = alpha * delta;

            double a1, a2;

            vectorField U(mesh.nCells());
            forAll(U, i)
            {
                //const scalar x = mesh.C()[i][0];
                const scalar z = mesh.C()[i][2];

                U[i][0] = 0.5 * (Utop + Ubottom) + 0.5 * (Utop - Ubottom) * erf(sqrt(PI) * z / delta);
                U[i][1] = 0;
                U[i][2] = 0;

            }

            scalar zS;
            scalar sqrtPI = sqrt(PI);
            for (int i = 0; i < 4; i++)
            {
                forAll(U, j)
                {
                    const scalar x = mesh.C()[j][0];
                    const scalar z = mesh.C()[j][2];
                    zS = z / delta;
                    a1 =  0.25 * exp(sqr(alphaS) * 0.25 / PI) * (erf(zS * sqrt(PI) + alphaS * 0.5 / sqrtPI) - 1);
                    a2 = -0.25 * exp(sqr(alphaS) * 0.25 / PI) * (erf(zS * sqrt(PI) - alphaS * 0.5 / sqrtPI) + 1);
                    U[j][0] += F2D * lambda1 * dU / Gamma3 * A[i] * delta * sin(alpha * x) * (a2 * exp(-z * alpha) - a1 * exp(z * alpha));//* delta
                    U[j][2] += F2D * lambda1 * dU / Gamma3 * A[i]  * cos(alpha * x) * (a2 * exp(-z * alpha) + a1 * exp(z * alpha));//* delta
                }
                alpha *= 0.5;
                alphaS *= 0.5;
            }

            /*
            double Gamma1 = lambda3*delta;
            double alpha2 = 2 * PI / lambda3;
            double alpha2S = alpha2 * delta;
            for (int i = 0; i < 4; i++)
            {
                forAll(U, j)
                {
                    const scalar x = mesh.C()[j][0];
                    const scalar y = mesh.C()[j][1];
                    const scalar z = mesh.C()[j][2];
                    zS = z / delta;
                    a1 =  0.25 * exp(sqr(alpha2S) * 0.25 / PI) * (erf(zS * sqrt(PI) + alpha2S * 0.5 / sqrtPI) - 1);
                    a2 = -0.25 * exp(sqr(alpha2S) * 0.25 / PI) * (erf(zS * sqrt(PI) - alpha2S * 0.5 / sqrtPI) + 1);
                    U[j][0] += F3D * lambda3 * dU / Gamma1 * B[i] * delta * sin(alpha2 * y) * (a2 * exp(-z * alpha) - a1 * exp(z * alpha));
                    U[j][1] += F3D * lambda3 * dU / Gamma1 * B[i] * delta * cos(alpha2 * y) * (a2 * exp(-z * alpha) + a1 * exp(z * alpha));
                }
                alpha *= 0.5;
                alphaS *= 0.5;
            }
            */

            U.writeEntry("", os);
//}}} end code
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

