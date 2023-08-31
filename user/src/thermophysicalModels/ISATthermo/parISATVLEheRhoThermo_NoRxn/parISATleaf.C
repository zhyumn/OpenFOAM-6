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

#include "parISATleaf.H"
#include "SVD.H"
#include "IOstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
Foam::parISATleaf::parISATleaf(int n_in, int n_out,
                               const scalarList &v) : node_(), numRetrieve_(0) //, pTimeTagList_(nullptr)
{
    value_.init(n_in);
    data_.init(n_out);
    A_.init(n_in, n_out);
    EOA_.init(n_in, n_in);
    forAll(value_, i)
    {
        value_[i] = v[i];
    }
}
Foam::parISATleaf::parISATleaf(int n_in, int n_out,
                               const scalarList &v,
                               SharedPointer<parISATNode> &node) : node_(node), numRetrieve_(0) //, pTimeTagList_(nullptr)
{
    value_.init(n_in);
    data_.init(n_out);
    A_.init(n_in, n_out);
    EOA_.init(n_in, n_in);
    forAll(value_, i)
    {
        value_[i] = v[i];
    }
}

void Foam::parISATleaf::init(int n_in, int n_out)
{
    lastUsed = 0;
    node_.offset = sptr_NULL;
    numRetrieve_ = 0;
    value_.init(n_in);
    data_.init(n_out);
    A_.init(n_in, n_out);
    EOA_.init(n_in, n_in);
    forAll(value_, i)
    {
        value_[i] = 0;
    }
}

void Foam::parISATleaf::reuse()
{
    lastUsed = 0;
    node_.offset = sptr_NULL;
    numRetrieve_ = 0;
    forAll(value_, i)
    {
        value_[i] = 0;
    }
}
Foam::parISATleaf::parISATleaf(int n_in, int n_out,
                               const scalarList &v,
                               SharedPointer<parISATNode> &node, const scalarList &data_in) : node_(node), numRetrieve_(0) //, pTimeTagList_(nullptr)
{
    value_.init(n_in);
    data_.init(n_out);
    A_.init(n_in, n_out);
    EOA_.init(n_in, n_in);
    forAll(value_, i)
    {
        value_[i] = v[i];
    }
    forAll(data_, i)
    {
        data_[i] = data_in[i];
    }
}
bool Foam::parISATleaf::inEOA(const scalarList &point, const scalarRectangularMatrix &scaleIn)
{
    static scalarList dx(value_.size());
    if (value_.size() > 4)
    {
        FatalErrorInFunction << " value_.size()=" << value_.size()
                             //<< ",offset=" << xx
                             << exit(FatalError);
    }
    forAll(value_, i)
    {
        dx[i] = (point[i]);
    }
    forAll(value_, i)
    {
        dx[i] = 1 / scaleIn[i][i];
    }
    forAll(value_, i)
    {
        dx[i] = (point[i] - value_[i]) / scaleIn[i][i];
    }
    scalar ret = 0;

    forAll(value_, i)
    {
        forAll(value_, j)
        {

            ret += dx[i] * dx[j] * EOA_(i, j);
        }
    }

    return ret <= 1.0;
}

void Foam::parISATleaf::eval(const scalarList &value, scalarList &ret)
{
    static scalarList dx(value_.size()), retm;
    for (int i = 0; i < value_.size(); i++)
        ret[i] = data_[i];
    for (int i = 0; i < value_.size(); i++)
        dx[i] = value[i] - value_[i];
    //retm = dx * A_;
    for (int i = 0; i < data_.size(); i++)
    {
        for (int j = 0; j < value_.size(); j++)
            ret[i] += dx[j] * A_(j, i);
    }
}

void Foam::parISATleaf::grow(const scalarList &point, const scalarRectangularMatrix &scaleIn)
{
    static scalarRectangularMatrix dx(value_.size(), 1);
    static scalarRectangularMatrix EOA_tmp(EOA_.n_in_, EOA_.n_out_);
    for (int i = 0; i < EOA_.n_in_; i++)
        for (int j = 0; j < EOA_.n_out_; j++)
        {
            EOA_tmp[i][j] = EOA_(i, j);
        }
    double pbp;
    for (int i = 0; i < value_.size(); i++)
        dx[i][0] = (point[i] - value_[i]) / scaleIn[i][i];
    pbp = (dx.T() * EOA_tmp * dx)[0][0];
    EOA_tmp = EOA_tmp + EOA_tmp * dx * (dx.T()) * EOA_tmp * (1 - pbp) / sqr(pbp);
    for (int i = 0; i < EOA_.n_in_; i++)
        for (int j = 0; j < EOA_.n_out_; j++)
        {
            EOA_(i, j) = EOA_tmp[i][j];
        }
}

void Foam::parISATleaf::increaseNumRetrieve()
{
    this->numRetrieve_++;
}

void Foam::parISATleaf::resetNumRetrieve()
{
    this->numRetrieve_ = 0;
}
// ************************************************************************* //
