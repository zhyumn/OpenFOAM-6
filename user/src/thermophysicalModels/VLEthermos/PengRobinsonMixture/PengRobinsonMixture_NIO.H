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

/*---------------------------------------------------------------------------*/
/* PengRobinson EOS                                                          */
/* Thermal solver                                                            */
/* Developed by   Ping Yi,  03/01/2019,   UMN                                */
/* mix means each phase mixture; ave means two-phase mixture                 */
/*---------------------------------------------------------------------------*/

#include "IFstream.H"
#include "IOmanip.H"
#include "OFstream.H"
#include "PengRobinson.H"
#include "PengRobinsonMixture.H"
#include "mathematicalConstants.H"
#include "specie.H"
#include "string.H"
#include "tmp.H"

using namespace Foam;
// * * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * //
/*
template <class ThermoMixture>
inline bool Foam::PengRobinsonMixture<ThermoMixture>::Ln_fugacityCoefficient_opt_NIO(
    scalar p, //input
    scalar T,
    const scalarList &Xinp,
    scalarList &fugcoef,
    scalarList &B_,
    scalarList &Tsqrt_AbyP_,
    scalar &Zmix,
    scalar &dg,
    const label pflag) const
{

    bool retflag = true;
    scalar Bmixz = 0.0;
    forAll(Xinp, i)
    {
        Bmixz += Xinp[i] * B_[i];
    }
    scalar Amixz = 0.0;
    static scalarList Amix2(this->N_);

    scalar rsqrT_P = p / (T * T);
    scalar temp, temp2;

    forAll(Xinp, i)
    {
        temp = 0.0;
        Amix2[i] = 0.0;
        forAll(Xinp, j)
        {
            temp2 = Xinp[j] * Tsqrt_AbyP_[j] * oneMbico[i][j];

            Amix2[i] += temp2;
            temp += temp2;
        }
        temp2 = Tsqrt_AbyP_[i];
        Amix2[i] *= rsqrT_P * temp2;
        Amixz += temp * Xinp[i] * temp2;
    }
    Amixz *= rsqrT_P;

    //A_A2(p, T, Xinp, Amixz, Amix2);

    const scalar small8 = 1.e-08; //PY: need to be noticed

    //scalar Z0, Zmix0;
    Zmix = Z_gibbs_mix(p, T, Amixz, Bmixz, Xinp);
    //std::tie(Z0, Zmix0) = Z_opt2(p, T, Amixz, Xinp);
    //Zmix = (pflag == 0) ? Z0 : Zmix0;

    const scalar cc0 = 2.0 * sqrt(2.0);
    const scalar cc1 = 1.0 + sqrt(2.0);
    const scalar cc2 = 1.0 - sqrt(2.0);

    scalar zfactor00 = Zmix;
    scalar ZMB = zfactor00 - Bmixz;

    if (ZMB < 0.0)
    {
        //FatalErrorInFunction << "ZMB < 0.0" << exit(FatalError);
        zfactor00 = Bmixz + small8;
        ZMB = small8;
        retflag = false;
    }
    scalar logZMB = log(ZMB);
    scalar logX = log((zfactor00 + cc1 * Bmixz) / (zfactor00 + cc2 * Bmixz));

    scalar coe_tt = Amixz / (cc0 * Bmixz) * logX;
    scalar coe_b = (zfactor00 - 1.0) / Bmixz + coe_tt / Bmixz;
    //scalar coe_0 = -logZMB;
    scalar coe_Amix2 = -coe_tt * 2.0 / Amixz;
    dg = (Zmix - 1 - logZMB) * RR * 1e-3 * T + RR * 1e-3 * T * (-Amixz / (cc0 * Bmixz)) * logX;

    forAll(Xinp, spid)
    {
        fugcoef[spid] = B_[spid] * coe_b + coe_Amix2 * Amix2[spid] - logZMB;
    }
    return retflag;
}*/

// ************************************************************************* //
