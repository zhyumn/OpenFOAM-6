/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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

#include "parallelISAT_1D.H"
namespace Foam
{
    parallelISAT_1D* ISAT_1D::pISAT = nullptr;
    std::ostream& operator<<(std::ostream& out, ISAT_1D::leafData& A)
    {
        out << A.v << ", " << A.Rv;
        return out;
    }
    std::ostream& operator<<(std::ostream& out, ISAT_1D::nodeData& A)
    {
        out << A.v;
        return out;
    }
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

//const dataType Foam::parallelISAT_1D::staticData();


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //




Foam::parallelISAT_1D::parallelISAT_1D()
    :
    parallelISAT(SUPstream::node_manager, 100000, SUPstream::Sync),
    tolerance_(1e-1)
{
    Foam::ISAT_1D::pISAT = this;
}






// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //




// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::parallelISAT_1D::~parallelISAT_1D()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::ISAT_1D::outputType Foam::ISAT_1D::leafData::func(const inputType& x)
{
    return sin(x);
}
Foam::ISAT_1D::gradientType Foam::ISAT_1D::leafData::gradFunc(const inputType& x)
{
    return cos(x);
}

void Foam::ISAT_1D::leafData::set(const inputType& x)
{
    v = x; Rv = func(x); computeA(v, Rv, A); EOA = Foam::ISAT_1D::pISAT->tolerance() / max(fabs(A), 1);
}


void Foam::ISAT_1D::leafData::set(const inputType& x, const outputType& y)
{
    v = x; Rv = y; computeA(x, y, A); EOA = Foam::ISAT_1D::pISAT->tolerance() / max(fabs(A), 1);
}
template<typename ...Args>
void Foam::ISAT_1D::leafData::computeA(const inputType& x, const outputType& y, gradientType& A, Args ...args)
{
    A = gradFunc(x);
}
bool Foam::ISAT_1D::leafData::checkSolution(const inputType& x, const outputType& y)
{
    outputType y_approx;
    retrieve(x, y_approx);
    return fabs(y_approx - y) < pISAT->tolerance();
}

// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


// ************************************************************************* //
