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

#include "ISATleaf.H"
#include "SVD.H"
#include "IOstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// Defined as static to be able to dynamically change it during simulations
// (all chemPoints refer to the same object)
/*
void Foam::ISATleaf::print(Foam::Ostream& OFout, int a)
{
    for (int i = 0;i < a;i++)
        OFout << " ";
    OFout << "*" << value_ << endl;
}
*/

Foam::ISATleaf::ISATleaf(int n_in, int n_out,
    const scalarList& v,
    ISATNode* node) :
    node_(node), value_(n_in),
    data_(n_out), A_(n_in, n_out),
    EOA_(n_in, n_in), numRetrieve_(0), pTimeTagList_(nullptr)
{
    forAll(value_, i)
    {
        value_[i] = v[i];
    }
}
Foam::ISATleaf::ISATleaf(int n_in, int n_out, const scalarList& v, ISATNode* node, const scalarList& data_in
) : node_(node), value_(n_in), data_(n_out), A_(n_in, n_out), EOA_(n_in, n_in), numRetrieve_(0), pTimeTagList_(nullptr)
{
    forAll(value_, i)
    {
        value_[i] = v[i];
    }
    forAll(data_, i)
    {
        data_[i] = data_in[i];
    }

}
bool Foam::ISATleaf::inEOA(const scalarList& point, const scalarRectangularMatrix& scaleIn)
{
    scalarRectangularMatrix dx(value_.size(), 1);
    forAll(value_, i)
    {
        dx[i][0] = (point[i] - value_[i]) / scaleIn[i][i];
    }
    
    return ((dx.T()) * EOA_ * dx)[0][0] <= 1.0;

}


void Foam::ISATleaf::eval(const scalarList& value, scalarList& ret)
{
    scalarRectangularMatrix dx(1, value_.size()), retm;
    ret = data_;
    for (int i = 0;i < value_.size();i++)
        dx[0][i] = value[i] - value_[i];
    retm = dx * A_;
    for (int i = 0;i < data_.size();i++)
        ret[i] += retm[0][i];
}

void Foam::ISATleaf::grow(const scalarList& point, const scalarRectangularMatrix& scaleIn)
{
    scalarRectangularMatrix dx(value_.size(), 1);
    double pbp;
    for (int i = 0;i < value_.size();i++)
        dx[i][0] = (point[i] - value_[i]) / scaleIn[i][i];
    pbp = (dx.T() * EOA_ * dx)[0][0];
    EOA_ = EOA_ + EOA_ * dx * (dx.T()) * EOA_ * (1 - pbp) / sqr(pbp);
}


void Foam::ISATleaf::increaseNumRetrieve()
{
    this->numRetrieve_++;
}



void Foam::ISATleaf::resetNumRetrieve()
{
    this->numRetrieve_ = 0;
}
// ************************************************************************* //
