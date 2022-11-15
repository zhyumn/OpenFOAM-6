/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2018 OpenFOAM Foundation
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

#include "PengRobinson.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class Specie>
Foam::PengRobinson<Specie>::PengRobinson(
    const dictionary &dict)
    : Specie(dict),
      Tc_(readScalar(dict.subDict("equationOfState").lookup("Tc"))),
      Vc_(1/readScalar(dict.subDict("equationOfState").lookup("rhoc"))),
      Zc_(1.0),
      Pc_(readScalar(dict.subDict("equationOfState").lookup("Pc"))),
      omega_(readScalar(dict.subDict("equationOfState").lookup("omega")))
      //Hig_phase_(dict.subDict("equationOfState").lookup("Hig_phase")),
      //Hig2_phase_(dict.subDict("equationOfState").lookup("Hig2_phase"))
{
    Zc_ = Pc_ * Vc_ / (RR * Tc_);
    sqrt_rPc_ = sqrt(1 / Pc_);
    sqrt_rTc_ = sqrt(1 / Tc_);
    coef_ = sqrt(0.45724) * sqrt_rPc_ * Tc_;

    if (omega_ > 0.49)
    {
        kappa_ = omega_ * (omega_ * (omega_ * (0.016666) - 0.164423) + 1.48503) + 0.379642;
    }
    else
    {
        kappa_ = omega_ * (omega_ * (-0.26992) + 1.54226) + 0.37464;
    }
    Aa_ = coef_ * (1.0 + kappa_);
    Ab_ = -coef_ * kappa_ * sqrt_rTc_;

    dAdTa_ = -0.45724 * Tc_ * sqr(kappa_) / Pc_;
    dAdTb_ = 3 * 0.45724 * kappa_ * (1 + kappa_) / Pc_ * Tc_ * sqrt(Tc_);
    dAdTc_ = -2 * 0.45724 * sqr(Tc_) * sqr(1 + kappa_) / Pc_;

    dAdPa_ = 0.45724 * Tc_ * sqr(kappa_) / Pc_;
    dAdPb_ = -2 * 0.45724 * kappa_ * (1 + kappa_) / Pc_ * Tc_ * sqrt(Tc_);
    dAdPc_ = 0.45724 * sqr(Tc_) * sqr(1 + kappa_) / Pc_;

    aa_ = RR * 1.0e-03 * Aa_;
    ab_ = RR * 1.0e-03 * Ab_;

    dadta_ = 0.45724 * sqr(RR * 1.0e-03) * Tc_ * sqr(kappa_) / Pc_;
    dadtb_ = -0.45724 * sqr(RR * 1.0e-03) * kappa_ * (1 + kappa_) / Pc_ * Tc_ * sqrt(Tc_);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class Specie>
void Foam::PengRobinson<Specie>::write(Ostream &os) const
{
    Specie::write(os);
    dictionary dict("PengRobinson");
    dict.add("Tc", Tc_);
    dict.add("Vc", Vc_);
    dict.add("Zc", Zc_);
    dict.add("Pc", Pc_);
    dict.add("omega", omega_);
    //dict.add("Hig_coef", Hig_phase_);
    //dict.add("Hig2_coef", Hig2_phase_);

    os << indent << dict.dictName() << dict;
}

// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template <class Specie>
Foam::Ostream &Foam::operator<<(
    Ostream &os,
    const PengRobinson<Specie> &pg)
{
    pg.write(os);
    return os;
}

// ************************************************************************* //
