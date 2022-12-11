/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "chungTransportMixture.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class ThermoMixture>
Foam::chungTransportMixture<ThermoMixture>::chungTransportMixture(const dictionary &dict, PtrList<SingleThermoType> &speciesData)
    : ThermoMixture(dict, speciesData)
{
}

template <class ThermoMixture>
inline Foam::chungTransportMixture<ThermoMixture>::chungTransportMixture(
    const word &name,
    PtrList<SingleThermoType> &speciesData,
    const speciesTable &specieNames,
    const dictionary &thermoDict)
    : ThermoMixture(name, speciesData, specieNames, thermoDict),
      w_global(this->N_),
      sigma_global(this->N_),
      e_k0_global(this->N_),
      kij_global(this->N_),
      mwij_global(this->N_),
      pow3_sigma_global(this->N_),
      sqr_sigma_global(this->N_),
      sqrt_mwij_global(this->N_),
      mu_M_global(this->N_),
      Mmd_global(this->N_),
      sigmd_global(this->N_),
      sqrt_Mmd_global(this->N_),
      sqr_sigmd_global(this->N_),
      mw_binary_M_global(this->N_),
      e_k_mix_M_global(this->N_),
      w_mix_M_global(this->N_),
      Dimix_opt_list1(this->N_),
      Dimix_opt_list2(this->N_),
      Dimix_opt_list3(this->N_),
      kappa_mu_opt_list1(this->N_)
{
    forAll(this->X_, i)
    {
        w_global[i].resize(this->N_);
        sigma_global[i].resize(this->N_);
        e_k0_global[i].resize(this->N_);
        kij_global[i].resize(this->N_);
        mwij_global[i].resize(this->N_);
        pow3_sigma_global[i].resize(this->N_);
        sqr_sigma_global[i].resize(this->N_);
        sqrt_mwij_global[i].resize(this->N_);
        mu_M_global[i].resize(this->N_);
        Mmd_global[i].resize(this->N_);
        sigmd_global[i].resize(this->N_);
        sqrt_Mmd_global[i].resize(this->N_);
        sqr_sigmd_global[i].resize(this->N_);
        Dimix_opt_list1[i].resize(this->N_);
        Dimix_opt_list2[i].resize(this->N_);
        kappa_mu_opt_list1[i].resize(this->N_);
        mw_binary_M_global[i].resize(this->N_);
        e_k_mix_M_global[i].resize(this->N_);
        w_mix_M_global[i].resize(this->N_);
    }
}

template <class ThermoMixture>
inline void Foam::chungTransportMixture<ThermoMixture>::chung_init()
{

    forAll(this->X_, spid)
    {
        forAll(this->X_, spjd)
        {
            w_global[spid][spjd] = 0.50 * ((*this)[spid].omega_ + (*this)[spjd].omega_);
            sigma_global[spid][spjd] = sqrt((0.8090 * cbrt((*this)[spid].Vc_ * 1.0e+03)) * (0.8090 * cbrt((*this)[spjd].Vc_ * 1.0e+03)));
            e_k0_global[spid][spjd] = sqrt((*this)[spid].Tc_ * (*this)[spjd].Tc_) / 1.2593;

            kij_global[spid][spjd] = sqrt((*this)[spid].kappa_ * (*this)[spid].kappa_);
            mwij_global[spid][spjd] = 2.0 * (*this)[spid].W() * (*this)[spjd].W() / ((*this)[spid].W() + (*this)[spjd].W());
            pow3_sigma_global[spid][spjd] = pow3(sigma_global[spid][spjd]);
            sqr_sigma_global[spid][spjd] = sqr(sigma_global[spid][spjd]);
            sqrt_mwij_global[spid][spjd] = sqrt(mwij_global[spid][spjd]);
            mu_M_global[spid][spjd] = sqr((*this)[spid].mu_) * sqr((*this)[spjd].mu_) / (pow3_sigma_global[spid][spjd] * e_k0_global[spid][spjd]);
            Mmd_global[spid][spjd] = 1 / (*this)[spid].W() + 1 / (*this)[spjd].W();
            sigmd_global[spid][spjd] = cbrt((*this)[spid].sigmvi_) + cbrt((*this)[spjd].sigmvi_);
            sqrt_Mmd_global[spid][spjd] = sqrt(Mmd_global[spid][spjd]);
            sqr_sigmd_global[spid][spjd] = sqr(sigmd_global[spid][spjd]);
            mw_binary_M_global[spid][spjd] = e_k0_global[spid][spjd] * sqr_sigma_global[spid][spjd] * sqrt_mwij_global[spid][spjd];
            e_k_mix_M_global[spid][spjd] = e_k0_global[spid][spjd] * pow3_sigma_global[spid][spjd];
            w_mix_M_global[spid][spjd] = w_global[spid][spjd] * pow3_sigma_global[spid][spjd];
        }
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class ThermoMixture>
void Foam::chungTransportMixture<ThermoMixture>::chungTransportMixture::write(Ostream &os) const
{
    os << this->name() << endl;
    os << token::BEGIN_BLOCK << incrIndent << nl;

    ThermoMixture::write(os);

    dictionary dict("transport");
    os << indent << dict.dictName() << dict;

    os << decrIndent << token::END_BLOCK << nl;
}

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template <class ThermoMixture>
Foam::Ostream &Foam::operator<<(Ostream &os, const chungTransportMixture<ThermoMixture> &ct)
{
    ct.write(os);
    return os;
}

// ************************************************************************* //
