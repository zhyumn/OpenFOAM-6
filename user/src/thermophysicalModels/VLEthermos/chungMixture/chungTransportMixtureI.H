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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class ThermoMixture>
inline Foam::chungTransportMixture<ThermoMixture>::chungTransportMixture(
    const ThermoMixture& t)
    : ThermoMixture(t)
{
}

template <class ThermoMixture>
inline Foam::chungTransportMixture<ThermoMixture>::chungTransportMixture(
    const word& name,
    const chungTransportMixture& ct)
    : ThermoMixture(name, ct)
{
}

template <class ThermoMixture>
inline Foam::autoPtr<Foam::chungTransportMixture<ThermoMixture>>
Foam::chungTransportMixture<ThermoMixture>::clone() const
{
    return autoPtr<chungTransportMixture<ThermoMixture>>(
        new chungTransportMixture<ThermoMixture>(*this));
}

template <class ThermoMixture>
inline Foam::autoPtr<Foam::chungTransportMixture<ThermoMixture>>
Foam::chungTransportMixture<ThermoMixture>::New(
    const dictionary& dict)
{
    return autoPtr<chungTransportMixture<ThermoMixture>>(
        new chungTransportMixture<ThermoMixture>(dict));
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class ThermoMixture>
inline Foam::scalar Foam::chungTransportMixture<ThermoMixture>::mu(
    scalar p,
    scalar T,
    const label pf,
    const scalarList* Xinp) const
{
    const scalarList* Xp = &this->X_;
    if (Xinp != NULL)
    {
        Xp = Xinp;
    }
    const scalarList& X = *Xp;

    scalarListList kij(this->N_);
    scalarListList mwij(this->N_);
    scalarListList w(this->N_);
    scalarListList sigma(this->N_);
    scalarListList e_k0(this->N_);
    forAll(X, i)
    {
        kij[i].resize(this->N_);
        mwij[i].resize(this->N_);
        w[i].resize(this->N_);
        sigma[i].resize(this->N_);
        e_k0[i].resize(this->N_);
    }

    scalarList k(this->N_, Zero); //need to be sure
    scalarList e_k(this->N_, Zero);
    scalarList mw_sp3(this->N_, Zero);
    scalarList Vc_sp3(this->N_, Zero);

    forAll(X, spid)
    {
        e_k[spid] = (*this)[spid].Tc_ / 1.2593;
        mw_sp3[spid] = (*this)[spid].W();           //g/mol
        Vc_sp3[spid] = (*this)[spid].Vc_ * 1.0e+03; //Vc is m3/kmol, Vc_sp3 is cm3/mol
    }
    forAll(X, spid)
    {
        forAll(X, spjd)
        {
            w[spid][spjd] = 0.50 * ((*this)[spid].omega_ + (*this)[spjd].omega_);
            sigma[spid][spjd] = sqrt((0.8090 * pow(Vc_sp3[spid], 1.0 / 3.0)) * (0.8090 * pow(Vc_sp3[spjd], 1.0 / 3.0)));
            //sigma(i,j)= bip_kij(i,j)*sqrt((0.809d0*vc(i)**(1.d0/3.d0))*(0.809d0*vc(j)**(1.d0/3.d0)))
            e_k0[spid][spjd] = sqrt(e_k[spid] * e_k[spjd]);
            //e_k0(i,j)= bip_kij(i,j)*sqrt(e_k(i)*e_k(j))
            //kij[spid][spjd] = sqrt(k[spid] * k[spjd]);
            kij[spid][spjd] = sqrt((*this)[spid].kappa_ * (*this)[spid].kappa_);
            mwij[spid][spjd] = 2.0 * mw_sp3[spid] * mw_sp3[spjd] / (mw_sp3[spid] + mw_sp3[spjd]);
        }
    }

    scalar sigma_mix3 = 0.0;
    scalar sigma_mix2 = 0.0;

    //scalar sigma_mix_23 = sigma_mix2/sigma_mix3;
    scalar e_k_mix = 0.0;
    scalar w_mix = 0.0;
    scalar k_mix = 0.0;
    scalar mu4 = 0.0;

    scalar mw_binary = 0.0;
    forAll(X, spid)
    {
        forAll(X, spjd)
        {
            sigma_mix3 += X[spid] * X[spjd] * pow3(sigma[spid][spjd]);
            sigma_mix2 += X[spid] * X[spjd] * sqr(sigma[spid][spjd]);
            e_k_mix += X[spid] * X[spjd] * e_k0[spid][spjd] * pow3(sigma[spid][spjd]);
        }
    }

    forAll(X, spid)
    {
        forAll(X, spjd)
        {

            w_mix += X[spid] * X[spjd] * w[spid][spjd] * pow3(sigma[spid][spjd]);
            mu4 += X[spid] * X[spjd] * sqr((*this)[spid].mu_) * sqr((*this)[spjd].mu_) / (pow3(sigma[spid][spjd]) * e_k0[spid][spjd]);
            k_mix += X[spid] * X[spjd] * kij[spid][spjd];
            mw_binary += (X[spid] * X[spjd] * e_k0[spid][spjd] * sqr(sigma[spid][spjd]) * sqrt(mwij[spid][spjd])); // (e_k_mix * sigma_mix2);
        }
    }
    e_k_mix /= sigma_mix3;
    w_mix /= sigma_mix3;
    mu4 *= sigma_mix3 * e_k_mix;
    scalar sigma_mix = pow(sigma_mix3, 1.0 / 3);
    scalar sigma_mix2_dev = sigma_mix * sigma_mix;

    mw_binary /= (e_k_mix * sigma_mix2_dev);
    mw_binary = sqr(mw_binary);

    scalar volu_mix = sigma_mix3 / pow3(0.8090); //Eq.(9-5.43);Eq.(9-4.8) bug??
    //scalar volu_mix= pow3(sigma_mix3/0.8090);
    scalar tempc_mix = 1.25930 * e_k_mix;
    scalar mur_mix4 = pow4(131.30) * mu4 / sqr(tempc_mix * volu_mix); //volu_mix in cm3/mol, mu4 in debyes
    scalar t_mix = T / e_k_mix;
    scalar Fcm = 1.0 - 0.2756 * w_mix + 0.0590350 * mur_mix4 + k_mix;
    const scalar A = 1.16145, B = 0.14874, C = 0.52487, D = 0.77320, E = 2.16178;
    const scalar F = 2.43787, G = -6.435e-4, H = 7.27371, S = 18.0323, W = -0.76830;

    scalar omegav = A * pow(t_mix, -B) + C * (exp(-D * t_mix)) + E * exp(-F * t_mix)
        + G * pow(t_mix, B) * sin(S * pow(t_mix, W) - H);
    //scalar omegav = 1.161450 * (pow(t_mix, -0.148740) + 0.524870 * exp(-0.773200 * t_mix) + 2.161780 * exp(-2.437870 * t_mix)) + (-6.435e-4) * pow(t_mix, 0.148740) * sin(18.0323 * pow(t_mix, -0.7683) - 7.27371); //???GT...

    scalar VcM = pow(sigma_mix / 0.809, 3); //✓


    scalar eta0 = (4.0785e-5) * sqrt(mw_binary * T) * Fcm / (omegav * pow(VcM, 2.0 / 3.0)); //unit P ✓


    const scalar  a0[10] = { 6.32402, 0.0012102, 5.28346, 6.62263, 19.74540,
                        -1.89992, 24.27450, 0.79716, -0.23816, 0.068629 };
    const scalar a1[10] = { 50.41190, -0.0011536, 254.209, 38.0957, 7.63034,
                       -12.53670, 3.44945, 1.11764, 0.067695,  0.34793 };
    const scalar a2[10] = { -51.68, -0.0062571, -168.481, -8.46414, -14.35440,
                        4.98529, -11.29130, 0.012348, -0.8163, 0.59256 };
    const scalar a3[10] = { 1189.02, 0.037283, 3898.270, 31.4178, 31.52670,
                             -18.1507, 69.3464, -4.1161, 4.02528, -0.72663 };//-4.11661,69.3466

    scalar MA[10] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    for (int i = 0; i < 10; i++)
    {
        // Notice that MA[0] = A1, A[1] = A2, ..., A[9] = A10 in eq (11) 
        MA[i] = a0[i] + a1[i] * w_mix + a2[i] * mu4 + a3[i] * k_mix;
    }

    scalar rho = this->rho(p, T, pf, Xinp) / this->W(Xinp);
    //scalar aa=1
    scalar Y = 1e-6 * rho * VcM / 6.0;
    scalar G1 = (1.0 - 0.5 * Y) / pow(1.0 - Y, 3.0);
    scalar G2 = (MA[0] * (1.0 - exp(-MA[3] * Y)) / Y + MA[1] * G1 * exp(MA[4] * Y) + MA[2] * G1) /
        (MA[0] * MA[3] + MA[1] + MA[2]);
    scalar etak = eta0 * (1.0 / G2 + MA[5] * Y);
    scalar etap = ((40.785e-6 / sqrt(t_mix)) * sqrt(mw_binary * T) / pow(VcM, 2.0 / 3.0))
        * MA[6] * Y * Y * G2 * exp(MA[7] + MA[8] / t_mix + MA[9] / sqr(t_mix));
    //1. mixture phases viscosity 40.785 ? 26.69 there is error in 9-5.24
    scalar finalMu = 0.1 * (etak + etap);
    return finalMu;
}


template <class ThermoMixture>
inline Foam::scalar Foam::chungTransportMixture<ThermoMixture>::mu_rho(
    scalar p,
    scalar T,
    scalar rho_in,
    const scalarList& Xinp) const
{

    scalarListList kij(this->N_);
    scalarListList mwij(this->N_);
    scalarListList w(this->N_);
    scalarListList sigma(this->N_);
    scalarListList e_k0(this->N_);
    forAll(Xinp, i)
    {
        kij[i].resize(this->N_);
        mwij[i].resize(this->N_);
        w[i].resize(this->N_);
        sigma[i].resize(this->N_);
        e_k0[i].resize(this->N_);
    }

    scalarList k(this->N_, Zero); //need to be sure
    scalarList e_k(this->N_, Zero);
    scalarList mw_sp3(this->N_, Zero);
    scalarList Vc_sp3(this->N_, Zero);

    forAll(Xinp, spid)
    {
        e_k[spid] = (*this)[spid].Tc_ / 1.2593;
        mw_sp3[spid] = (*this)[spid].W();           //g/mol
        Vc_sp3[spid] = (*this)[spid].Vc_ * 1.0e+03; //Vc is m3/kmol, Vc_sp3 is cm3/mol
    }
    forAll(Xinp, spid)
    {
        forAll(Xinp, spjd)
        {
            w[spid][spjd] = 0.50 * ((*this)[spid].omega_ + (*this)[spjd].omega_);
            sigma[spid][spjd] = sqrt((0.8090 * pow(Vc_sp3[spid], 1.0 / 3.0)) * (0.8090 * pow(Vc_sp3[spjd], 1.0 / 3.0)));
            //sigma(i,j)= bip_kij(i,j)*sqrt((0.809d0*vc(i)**(1.d0/3.d0))*(0.809d0*vc(j)**(1.d0/3.d0)))
            e_k0[spid][spjd] = sqrt(e_k[spid] * e_k[spjd]);
            //e_k0(i,j)= bip_kij(i,j)*sqrt(e_k(i)*e_k(j))
            //kij[spid][spjd] = sqrt(k[spid] * k[spjd]);
            kij[spid][spjd] = sqrt((*this)[spid].kappa_ * (*this)[spid].kappa_);
            mwij[spid][spjd] = 2.0 * mw_sp3[spid] * mw_sp3[spjd] / (mw_sp3[spid] + mw_sp3[spjd]);
        }
    }

    scalar sigma_mix3 = 0.0;
    scalar sigma_mix2 = 0.0;

    //scalar sigma_mix_23 = sigma_mix2/sigma_mix3;
    scalar e_k_mix = 0.0;
    scalar w_mix = 0.0;
    scalar k_mix = 0.0;
    scalar mu4 = 0.0;

    scalar mw_binary = 0.0;
    forAll(Xinp, spid)
    {
        forAll(Xinp, spjd)
        {
            sigma_mix3 += Xinp[spid] * Xinp[spjd] * pow3(sigma[spid][spjd]);
            sigma_mix2 += Xinp[spid] * Xinp[spjd] * sqr(sigma[spid][spjd]);
            e_k_mix += Xinp[spid] * Xinp[spjd] * e_k0[spid][spjd] * pow3(sigma[spid][spjd]);
        }
    }

    forAll(Xinp, spid)
    {
        forAll(Xinp, spjd)
        {

            w_mix += Xinp[spid] * Xinp[spjd] * w[spid][spjd] * pow3(sigma[spid][spjd]);
            mu4 += Xinp[spid] * Xinp[spjd] * sqr((*this)[spid].mu_) * sqr((*this)[spjd].mu_) / (pow3(sigma[spid][spjd]) * e_k0[spid][spjd]);
            k_mix += Xinp[spid] * Xinp[spjd] * kij[spid][spjd];
            mw_binary += (Xinp[spid] * Xinp[spjd] * e_k0[spid][spjd] * sqr(sigma[spid][spjd]) * sqrt(mwij[spid][spjd])); // (e_k_mix * sigma_mix2);
        }
    }
    e_k_mix /= sigma_mix3;
    w_mix /= sigma_mix3;
    mu4 *= sigma_mix3 * e_k_mix;
    scalar sigma_mix = pow(sigma_mix3, 1.0 / 3);
    scalar sigma_mix2_dev = sigma_mix * sigma_mix;

    mw_binary /= (e_k_mix * sigma_mix2_dev);
    mw_binary = sqr(mw_binary);

    scalar volu_mix = sigma_mix3 / pow3(0.8090); //Eq.(9-5.43);Eq.(9-4.8) bug??
    //scalar volu_mix= pow3(sigma_mix3/0.8090);
    scalar tempc_mix = 1.25930 * e_k_mix;
    scalar mur_mix4 = pow4(131.30) * mu4 / sqr(tempc_mix * volu_mix); //volu_mix in cm3/mol, mu4 in debyes
    scalar t_mix = T / e_k_mix;
    scalar Fcm = 1.0 - 0.2756 * w_mix + 0.0590350 * mur_mix4 + k_mix;
    const scalar A = 1.16145, B = 0.14874, C = 0.52487, D = 0.77320, E = 2.16178;
    const scalar F = 2.43787, G = -6.435e-4, H = 7.27371, S = 18.0323, W = -0.76830;

    scalar omegav = A * pow(t_mix, -B) + C * (exp(-D * t_mix)) + E * exp(-F * t_mix)
        + G * pow(t_mix, B) * sin(S * pow(t_mix, W) - H);
    //scalar omegav = 1.161450 * (pow(t_mix, -0.148740) + 0.524870 * exp(-0.773200 * t_mix) + 2.161780 * exp(-2.437870 * t_mix)) + (-6.435e-4) * pow(t_mix, 0.148740) * sin(18.0323 * pow(t_mix, -0.7683) - 7.27371); //???GT...

    scalar VcM = pow(sigma_mix / 0.809, 3); //✓


    scalar eta0 = (4.0785e-5) * sqrt(mw_binary * T) * Fcm / (omegav * pow(VcM, 2.0 / 3.0)); //unit P ✓


    const scalar  a0[10] = { 6.32402, 0.0012102, 5.28346, 6.62263, 19.74540,
                        -1.89992, 24.27450, 0.79716, -0.23816, 0.068629 };
    const scalar a1[10] = { 50.41190, -0.0011536, 254.209, 38.0957, 7.63034,
                       -12.53670, 3.44945, 1.11764, 0.067695,  0.34793 };
    const scalar a2[10] = { -51.68, -0.0062571, -168.481, -8.46414, -14.35440,
                        4.98529, -11.29130, 0.012348, -0.8163, 0.59256 };
    const scalar a3[10] = { 1189.02, 0.037283, 3898.270, 31.4178, 31.52670,
                             -18.1507, 69.3464, -4.1161, 4.02528, -0.72663 };//-4.11661,69.3466

    scalar MA[10] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    for (int i = 0; i < 10; i++)
    {
        // Notice that MA[0] = A1, A[1] = A2, ..., A[9] = A10 in eq (11) 
        MA[i] = a0[i] + a1[i] * w_mix + a2[i] * mu4 + a3[i] * k_mix;
    }

    scalar rho = rho_in / this->W(Xinp);
    //scalar aa=1
    scalar Y = 1e-6 * rho * VcM / 6.0;
    scalar G1 = (1.0 - 0.5 * Y) / pow(1.0 - Y, 3.0);
    scalar G2 = (MA[0] * (1.0 - exp(-MA[3] * Y)) / Y + MA[1] * G1 * exp(MA[4] * Y) + MA[2] * G1) /
        (MA[0] * MA[3] + MA[1] + MA[2]);
    scalar etak = eta0 * (1.0 / G2 + MA[5] * Y);
    scalar etap = ((40.785e-6 / sqrt(t_mix)) * sqrt(mw_binary * T) / pow(VcM, 2.0 / 3.0))
        * MA[6] * Y * Y * G2 * exp(MA[7] + MA[8] / t_mix + MA[9] / sqr(t_mix));
    //1. mixture phases viscosity 40.785 ? 26.69 there is error in 9-5.24
    scalar finalMu = 0.1 * (etak + etap);
    return finalMu;
}

//- Return mass Diffusivity of specie ith [m^2/s] based on Fick's law using 
//  binary Diffusion coefficients based on Fuller and Takahashi models
template <class Specie>
inline Foam::scalar Foam::chungTransportMixture<Specie>::Dimix
(
    const scalar p,
    const scalar T,
    label speciei
) const
{

    //const scalarList* Xp = &this->X_;
    //if (Xinp != NULL)
    //{
    //    Xp = Xinp;
    //}
    const scalarList& X = this->X_;
    scalar sum = 0.0;
    scalar Dimix = 0.0;
    scalarList Ymd(this->N_, Zero);
    scalar sumY = 0.0;
    forAll(X, spid)
    {
        Ymd[spid] = X[spid] * (*this)[spid].W();
        sumY += X[spid] * (*this)[spid].W();
    }
    forAll(X, spid)
    {
        Ymd[spid] /= sumY;
    }
    scalarListList Pcmd(this->N_);
    scalarListList Tcmd(this->N_);
    scalarListList Mmd(this->N_);
    scalarListList sigmd(this->N_);

    forAll(X, i)
    {
        Pcmd[i].resize(this->N_);
        Tcmd[i].resize(this->N_);
        Mmd[i].resize(this->N_);
        sigmd[i].resize(this->N_);
    }

    forAll(X, spid)
    {
        forAll(X, spjd)
        {
            if (X[spid] + X[spjd] != 0)
            {
                Pcmd[spid][spjd] = (X[spid] * (*this)[spid].Pc_ + X[spjd] * (*this)[spjd].Pc_) / (X[spid] + X[spjd]);
                Tcmd[spid][spjd] = (X[spid] * (*this)[spid].Tc_ + X[spjd] * (*this)[spjd].Tc_) / (X[spid] + X[spjd]);
            }
            else
            {
                Tcmd[spid][spjd] = 1e-40;
                Pcmd[spid][spjd] = 1e-40;
            }
            Mmd[spid][spjd] = 1 / (*this)[spid].W() + 1 / (*this)[spjd].W();
            sigmd[spid][spjd] = pow((*this)[spid].sigmvi_, 1.0 / 3) + pow((*this)[spjd].sigmvi_, 1.0 / 3);
        }
    }
    forAll(X, j)
    {
        if (speciei != j)
        {
            sum = sum + X[j] /
                (
                    phi(p / (Pcmd[j][speciei]), T / (Tcmd[j][speciei]))
                    * 0.001 * pow(T, 1.75) * pow(Mmd[j][speciei], 0.5)
                    / ((p / 101325) * pow(sigmd[j][speciei], 2))
                    );
        }
    }
    //Dimix = (1-Ymd_[speciei])/sum;
    if (sum == 0.0) { Dimix = 0; }
    else { Dimix = (1 - Ymd[speciei]) / (sum); }

    return Dimix / 1e4;
}

//- Return coefficient phi based on Takahashi model
template <class Specie>
inline Foam::scalar Foam::chungTransportMixture<Specie>::phi
(
    scalar PRR,
    scalar TRR
) const
{
    const scalar A01 = 0.038042, B01 = 1.52267;
    const scalar A02 = 0.067433, B02 = 2.16794;
    const scalar A03 = 0.098317, B03 = 2.42910;
    const scalar A04 = 0.137610, B04 = 2.77605;
    const scalar A05 = 0.175081, B05 = 2.98256;
    const scalar A06 = 0.216376, B06 = 3.11384;
    const scalar A08 = 0.314051, B08 = 3.50264;
    //    const scalar C08 = 0.000000, E08 = 0.00000; //unused
    const scalar A10 = 0.385736, B10 = 3.07773;
    const scalar C10 = 0.141211, E10 = 13.45454;
    const scalar A12 = 0.514553, B12 = 3.54744;
    const scalar C12 = 0.278407, E12 = 14.00000;
    const scalar A14 = 0.599184, B14 = 3.61216;
    const scalar C14 = 0.372683, E14 = 10.00900;
    const scalar A16 = 0.557725, B16 = 3.41882;
    const scalar C16 = 0.504894, E16 = 8.57519;
    const scalar A18 = 0.593007, B18 = 3.18415;
    const scalar C18 = 0.678469, E18 = 10.37483;
    const scalar A20 = 0.696001, B20 = 3.37660;
    const scalar C20 = 0.665702, E20 = 11.21674;
    const scalar A25 = 0.790770, B25 = 3.27984;
    const scalar C25 = 0.000000, E25 = 0.00000;
    const scalar A30 = 0.502100, B30 = 2.39031;
    const scalar C30 = 0.602907, E30 = 6.19043;
    const scalar A40 = 0.837452, B40 = 3.23513;
    const scalar C40 = 0.000000, E40 = 0.00000;
    const scalar A50 = 0.890390, B50 = 3.13001;

    scalar Y = 0, YA = 0, YB = 0;
    TRR = (TRR < 1.0) ? 1.0 : TRR;
    TRR = (TRR > 5.0) ? 5.0 : TRR;
    PRR = (PRR > 5.0) ? 5.0 : PRR;

    if (PRR < 0.1)
    {
        Y = 1.01 * (1.0 - A01 * pow(TRR, -B01));
        Y = 1.0 + (Y - 1.0) / (0.1 - 0.0) * (PRR - 0.0);
    }
    else if ((PRR >= 0.1) && (PRR < 0.2))
    {
        YA = 1.01 * (1.0 - A01 * pow(TRR, -B01));
        YB = 1.01 * (1.0 - A02 * pow(TRR, -B02));
        Y = YA + (YB - YA) / (0.2 - 0.1) * (PRR - 0.1);
    }
    else if ((PRR >= 0.2) && (PRR < 0.3))
    {
        YA = 1.01 * (1.0 - A02 * pow(TRR, -B02));
        YB = 1.01 * (1.0 - A03 * pow(TRR, -B03));
        Y = YA + (YB - YA) / (0.3 - 0.2) * (PRR - 0.2);
    }
    else if ((PRR >= 0.3) && (PRR < 0.4))
    {
        YA = 1.01 * (1.0 - A03 * pow(TRR, -B03));
        YB = 1.01 * (1.0 - A04 * pow(TRR, -B04));
        Y = YA + (YB - YA) / (0.4 - 0.3) * (PRR - 0.3);
    }
    else if ((PRR >= 0.4) && (PRR < 0.5))
    {
        YA = 1.01 * (1.0 - A04 * pow(TRR, -B04));
        YB = 1.01 * (1.0 - A05 * pow(TRR, -B05));
        Y = YA + (YB - YA) / (0.5 - 0.4) * (PRR - 0.4);
    }
    else if ((PRR >= 0.5) && (PRR < 0.6))
    {
        YA = 1.01 * (1.0 - A05 * pow(TRR, -B05));
        YB = 1.01 * (1.0 - A06 * pow(TRR, -B06));
        Y = YA + (YB - YA) / (0.6 - 0.5) * (PRR - 0.5);
    }
    else if ((PRR >= 0.6) && (PRR < 0.8))
    {
        YA = 1.01 * (1.0 - A06 * pow(TRR, -B06));
        YB = 1.01 * (1.0 - A08 * pow(TRR, -B08));
        Y = YA + (YB - YA) / (0.8 - 0.6) * (PRR - 0.6);
    }
    else if ((PRR >= 0.8) && (PRR < 1.0))
    {
        YA = 1.01 * (1.0 - A08 * pow(TRR, -B08));
        YB = 1.02 * (1.0 - A10 * pow(TRR, -B10)) * (1.0 - C10 * pow(TRR, -E10));
        Y = YA + (YB - YA) / (1.0 - 0.8) * (PRR - 0.8);
    }
    else if ((PRR >= 1.0) && (PRR < 1.2))
    {
        YA = 1.02 * (1.0 - A10 * pow(TRR, -B10)) * (1.0 - C10 * pow(TRR, -E10));
        YB = 1.02 * (1.0 - A12 * pow(TRR, -B12)) * (1.0 - C12 * pow(TRR, -E12));
        Y = YA + (YB - YA) / (1.2 - 1.0) * (PRR - 1.0);
    }
    else if ((PRR >= 1.2) && (PRR < 1.4))
    {
        YA = 1.02 * (1.0 - A12 * pow(TRR, -B12)) * (1.0 - C12 * pow(TRR, -E12));
        YB = 1.02 * (1.0 - A14 * pow(TRR, -B14)) * (1.0 - C14 * pow(TRR, -E14));
        Y = YA + (YB - YA) / (1.4 - 1.2) * (PRR - 1.2);
    }
    else if ((PRR >= 1.4) && (PRR < 1.6))
    {
        YA = 1.02 * (1.0 - A14 * pow(TRR, -B14)) * (1.0 - C14 * pow(TRR, -E14));
        YB = 1.02 * (1.0 - A16 * pow(TRR, -B16)) * (1.0 - C16 * pow(TRR, -E16));
        Y = YA + (YB - YA) / (1.6 - 1.4) * (PRR - 1.4);
    }
    else if ((PRR >= 1.6) && (PRR < 1.8))
    {
        YA = 1.02 * (1.0 - A16 * pow(TRR, -B16)) * (1.0 - C16 * pow(TRR, -E16));
        YB = 1.03 * (1.0 - A18 * pow(TRR, -B18)) * (1.0 - C18 * pow(TRR, -E18));
        Y = YA + (YB - YA) / (1.8 - 1.6) * (PRR - 1.6);
    }
    else if ((PRR >= 1.8) && (PRR < 2.0))
    {
        YA = 1.03 * (1.0 - A18 * pow(TRR, -B18)) * (1.0 - C18 * pow(TRR, -E18));
        YB = 1.03 * (1.0 - A20 * pow(TRR, -B20)) * (1.0 - C20 * pow(TRR, -E20));
        Y = YA + (YB - YA) / (2.0 - 1.8) * (PRR - 1.8);
    }
    else if ((PRR >= 2.0) && (PRR < 2.5))
    {
        YA = 1.03 * (1.0 - A20 * pow(TRR, -B20)) * (1.0 - C20 * pow(TRR, -E20));
        YB = 1.04 * (1.0 - A25 * pow(TRR, -B25)) * (1.0 - C25 * pow(TRR, -E25));
        Y = YA + (YB - YA) / (2.5 - 2.0) * (PRR - 2.0);
    }
    else if ((PRR >= 2.5) && (PRR < 3.0))
    {
        YA = 1.04 * (1.0 - A25 * pow(TRR, -B25)) * (1.0 - C25 * pow(TRR, -E25));
        YB = 1.05 * (1.0 - A30 * pow(TRR, -B30)) * (1.0 - C30 * pow(TRR, -E30));
        Y = YA + (YB - YA) / (3.0 - 2.5) * (PRR - 2.5);
    }
    else if ((PRR >= 3.0) && (PRR < 4.0))
    {
        YA = 1.05 * (1.0 - A30 * pow(TRR, -B30)) * (1.0 - C30 * pow(TRR, -E30));
        YB = 1.06 * (1.0 - A40 * pow(TRR, -B40)) * (1.0 - C40 * pow(TRR, -E40));
        Y = YA + (YB - YA) / (4.0 - 3.0) * (PRR - 3.0);
    }
    else if ((PRR >= 4.0) && (PRR < 5.0))
    {
        YA = 1.06 * (1.0 - A40 * pow(TRR, -B40)) * (1.0 - C40 * pow(TRR, -E40));
        YB = 1.07 * (1.0 - A50 * pow(TRR, -B50));
        Y = YA + (YB - YA) / (5.0 - 4.0) * (PRR - 4.0);
    }
    else if (PRR >= 5.0)
    {
        Y = 1.07 * (1.0 - A50 * pow(TRR, -B50));
    }

    if (Y == 0) { Y = 1e-16; }
    return Y;
}

template <class Specie>
inline Foam::scalar Foam::chungTransportMixture<Specie>::kappa(
    scalar p,
    scalar T,
    const label pf,
    const scalarList* Xinp) const
{

    const scalarList* Xp = &this->X_;
    if (Xinp != NULL)
    {
        Xp = Xinp;
    }
    const scalarList& X = *Xp;

    scalarListList kij(this->N_);
    scalarListList mwij(this->N_);
    scalarListList w(this->N_);
    scalarListList sigma(this->N_);
    scalarListList e_k0(this->N_);
    forAll(X, i)
    {
        kij[i].resize(this->N_);
        mwij[i].resize(this->N_);
        w[i].resize(this->N_);
        sigma[i].resize(this->N_);
        e_k0[i].resize(this->N_);
    }

    scalarList k(this->N_, Zero); //need to be sure
    scalarList e_k(this->N_, Zero);
    scalarList mw_sp3(this->N_, Zero);
    scalarList Vc_sp3(this->N_, Zero);

    forAll(X, spid)
    {
        e_k[spid] = (*this)[spid].Tc_ / 1.2593;
        mw_sp3[spid] = (*this)[spid].W();           //g/mol
        Vc_sp3[spid] = (*this)[spid].Vc_ * 1.0e+03; //Vc is m3/kmol, Vc_sp3 is cm3/mol
    }
    forAll(X, spid)
    {
        forAll(X, spjd)
        {
            w[spid][spjd] = 0.50 * ((*this)[spid].omega_ + (*this)[spjd].omega_);
            sigma[spid][spjd] = sqrt((0.8090 * pow(Vc_sp3[spid], 1.0 / 3.0)) * (0.8090 * pow(Vc_sp3[spjd], 1.0 / 3.0)));
            //sigma(i,j)= bip_kij(i,j)*sqrt((0.809d0*vc(i)**(1.d0/3.d0))*(0.809d0*vc(j)**(1.d0/3.d0)))
            e_k0[spid][spjd] = sqrt(e_k[spid] * e_k[spjd]);
            //e_k0(i,j)= bip_kij(i,j)*sqrt(e_k(i)*e_k(j))
            //kij[spid][spjd] = sqrt(k[spid] * k[spjd]);
            kij[spid][spjd] = sqrt((*this)[spid].kappa_ * (*this)[spid].kappa_);
            mwij[spid][spjd] = 2.0 * mw_sp3[spid] * mw_sp3[spjd] / (mw_sp3[spid] + mw_sp3[spjd]);
        }
    }

    scalar sigma_mix3 = 0.0;
    scalar sigma_mix2 = 0.0;

    //scalar sigma_mix_23 = sigma_mix2/sigma_mix3;
    scalar e_k_mix = 0.0;
    scalar w_mix = 0.0;
    scalar k_mix = 0.0;
    scalar mu4 = 0.0;

    scalar mw_binary = 0.0;
    forAll(X, spid)
    {
        forAll(X, spjd)
        {
            sigma_mix3 += X[spid] * X[spjd] * pow3(sigma[spid][spjd]);
            sigma_mix2 += X[spid] * X[spjd] * sqr(sigma[spid][spjd]);
            e_k_mix += X[spid] * X[spjd] * e_k0[spid][spjd] * pow3(sigma[spid][spjd]);
        }
    }
    scalar sigma_mix = pow(sigma_mix3, 1.0 / 3);
    scalar sigma_mix2_dev = sigma_mix * sigma_mix;

    forAll(X, spid)
    {
        forAll(X, spjd)
        {

            w_mix += X[spid] * X[spjd] * w[spid][spjd] * pow3(sigma[spid][spjd]);
            mu4 += X[spid] * X[spjd] * sqr((*this)[spid].mu_) * sqr((*this)[spjd].mu_) / (pow3(sigma[spid][spjd]) * e_k0[spid][spjd]);
            k_mix += X[spid] * X[spjd] * kij[spid][spjd];
            mw_binary += (X[spid] * X[spjd] * e_k0[spid][spjd] * sqr(sigma[spid][spjd]) * sqrt(mwij[spid][spjd])); // (e_k_mix * sigma_mix2);
        }
    }
    e_k_mix /= sigma_mix3;
    w_mix /= sigma_mix3;
    mu4 *= sigma_mix3 * e_k_mix;
    mw_binary /= (e_k_mix * sigma_mix2_dev);
    mw_binary = sqr(mw_binary);

    scalar volu_mix = sigma_mix3 / pow3(0.8090); //Eq.(9-5.43);Eq.(9-4.8) bug??
    //scalar volu_mix= pow3(sigma_mix3/0.8090);
    scalar tempc_mix = 1.25930 * e_k_mix;
    //scalar mur_mix4 = pow4(131.30) * mu4 / sqr(tempc_mix * volu_mix); //volu_mix in cm3/mol, mu4 in debyes
    //scalar t_mix = T / e_k_mix;
    //scalar omegav = 1.161450 * (pow(t_mix, -0.148740) + 0.524870 * exp(-0.773200 * t_mix) + 2.161780 * exp(-2.437870 * t_mix)) + (-6.435e-4) * pow(t_mix, 0.148740) * sin(18.0323 * pow(t_mix, -0.7683) - 7.27371); //???GT...

    //1. mixture phases viscosity 40.785 ? 26.69 there is error in 9-5.24
    //viscos = 1e-07 * 26.69 * Fcm * sqrt(mw_binary * T) / (sigma_mix2 * omegav);

    scalar mw_phase = this->W(Xinp);
    scalar Cv_phase = this->Cp(p, T, pf, Xinp) - this->CpMCv(p, T, pf, Xinp); //J/kgK //warning 
    //scalar Cv_phase = this->Cp(p, T);
    scalar alpha = Cv_phase * mw_phase / (RR * 1.0e-03) - 1.50;
    scalar tempr = T / tempc_mix;

    //compute specific volume, m3/mol
    //scalar amix_phase = this->a(p,T,Xinp), bmix_phase = this->b(p,T,Xinp), dadt_phase = this->dadt(p,T,Xinp), dadt2 = this->d2adt2(p,T,Xinp);
    scalar sv_phase = this->V(p, T, pf, Xinp); //m3/mol

    //2, heat conductivity,W/m.K
    scalar beta = 0.78620 - 0.71090 * w_mix + 1.31680 * sqr(w_mix);
    scalar mw_binary0 = mw_binary / 1.0e+3;
    scalar zeta = 2.00 + 10.50 * sqr(tempr);
    scalar psi = 1.0 + alpha * ((0.2150 + 0.282880 * alpha - 1.0610 * beta + 0.266650 * zeta) / (0.63660 + beta * zeta + 1.0610 * alpha * beta));

    scalar qq = 3.5860 * 1e-03 * sqrt(tempc_mix / mw_binary0) / pow(volu_mix, 2.0 / 3.0);
    scalar yy = volu_mix / (6.0e+06 * sv_phase);
    scalar G1 = (1.0 - 0.50 * yy) / pow3(1.0 - yy);

    // Table 10-3
    //scalar a[7] = { 2.41660, -5.0924e-01, 6.61070, 1.4543e+01, 7.9274e-01, -5.86340, 9.1098e+01 };
    //scalar bb[7] = { 7.4824e-01, -1.50940, 5.62070, -8.91390, 8.2019e-01, 1.2801e+01, 1.2811e+02 };
    //scalar c[7] = { -9.1858e-01, -4.9991e+01, 6.4760e+01, -5.63790, -6.9369e-01, 9.58930, -5.4217e+01 };
    //scalar d[7] = { 1.2712e+02, 6.9983e+01, 2.7039e+01, 7.4344e+01, 6.31730, 6.5529e+01, 5.2381e+02 };
    scalar mur_mix4 = pow4(131.30) * mu4 / sqr(tempc_mix * volu_mix);
    const scalar b0[7] = { 2.41657, -0.50924, 6.61069, 14.54250, 0.79274,
                          -5.8634, 81.171 };
    const scalar b1[7] = { 0.74824, -1.50936, 5.62073, -8.91387, 0.82019,
                          12.8005, 114.1580 };
    const scalar b2[7] = { -0.91858, -49.99120, 64.7599, -5.63794, -0.69369,
                           9.58926, -60.841 };
    const scalar b3[7] = { 121.721, 69.9834, 27.0389, 74.3435, 6.31734,
                          -65.5292, 466.7750 };

    scalar MB[7] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    for (int j = 0; j < 7; j++)
    {
        // Notice that MB[0] = B1, B[1] = B2, ..., B[6] = B7 in eq (13) 
        MB[j] = b0[j] + b1[j] * w_mix + b2[j] * mur_mix4 + b3[j] * k_mix;
    }
    //scalar B_[8] = { 0.0 };
    //for (label ii = 0; ii < 7; ii++)
    //{
    //    B_[ii] = a[ii] + bb[ii] * w_mix + c[ii] * mu4 + d[ii] * k_mix;
    //}

    //scalar G2 = ((B_[1] / yy) * (1.0 - exp(-B_[4] * yy)) + B_[2] * G1 * exp(B_[5] * yy) + B_[3] * G1) / (B_[1] * B_[4] + B_[2] + B_[3]);
    //Info<< "G2* * * *"<<G2<<"* * "<< yy <<"* * "<< volu_mix <<"**"<<sv_phase << endl;
    //Info<< "comp_phase* * * *"<<comp_phase[0] <<"* * "<< comp_phase[1] << endl;

    scalar t_mix = T / e_k_mix;
    scalar Fcm = 1.0 - 0.2756 * w_mix + 0.0590350 * mur_mix4 + k_mix;
    const scalar A = 1.16145, B = 0.14874, C = 0.52487, D = 0.77320, E = 2.16178;
    const scalar F = 2.43787, G = -6.435e-4, H = 7.27371, S = 18.0323, W = -0.76830;

    scalar omegav = A * pow(t_mix, -B) + C * (exp(-D * t_mix)) + E * exp(-F * t_mix)
        + G * pow(t_mix, B) * sin(S * pow(t_mix, W) - H);
    scalar VcM = pow(sigma_mix / 0.809, 3); //✓
    scalar eta0 = (4.0785e-5) * sqrt(mw_binary * T) * Fcm / (omegav * pow(VcM, 2.0 / 3.0)); //unit P ✓
    scalar lamda0 = 7.452 * (eta0 / mw_binary) * psi;

    scalar rho = this->rho(p, T, pf, Xinp) / this->W(Xinp);

    scalar Y = 1e-6 * rho * VcM / 6.0;
    scalar H2 = (MB[0] * (1.0 - exp(-MB[3] * Y)) / Y + MB[1] * G1 * exp(MB[4] * Y) + MB[2] * G1) /
        (MB[0] * MB[3] + MB[1] + MB[2]);
    scalar lamdak = lamda0 * (1.0 / H2 + MB[5] * Y);
    scalar Tr = T / tempc_mix;
    scalar lamdap = (3.039e-4 * pow(tempc_mix / mw_binary, 1.0 / 2.0) / pow(VcM, 2.0 / 3.0)) * MB[6] * pow(Y, 2.0) * H2 * pow(Tr, 1.0 / 2.0);
    //return 31.20 * mu(p, T, Xinp) * psi * (1.0 / G2 + B_[6] * yy) / mw_binary0 + qq * B_[7] * sqr(yy) * sqrt(tempr) * G2;
    scalar kappa = 418.6798 * (lamdak + lamdap);
    return kappa;
    //printf("kappa_phase_thcond_phase %5.8f\n\n",thcond_phase);
}


template <class Specie>
inline Foam::scalar Foam::chungTransportMixture<Specie>::kappa_rho(
    scalar p,
    scalar T,
    scalar rho_in,
    const scalarList& Xinp) const
{

    scalarListList kij(this->N_);
    scalarListList mwij(this->N_);
    scalarListList w(this->N_);
    scalarListList sigma(this->N_);
    scalarListList e_k0(this->N_);
    forAll(Xinp, i)
    {
        kij[i].resize(this->N_);
        mwij[i].resize(this->N_);
        w[i].resize(this->N_);
        sigma[i].resize(this->N_);
        e_k0[i].resize(this->N_);
    }

    scalarList k(this->N_, Zero); //need to be sure
    scalarList e_k(this->N_, Zero);
    scalarList mw_sp3(this->N_, Zero);
    scalarList Vc_sp3(this->N_, Zero);

    forAll(Xinp, spid)
    {
        e_k[spid] = (*this)[spid].Tc_ / 1.2593;
        mw_sp3[spid] = (*this)[spid].W();           //g/mol
        Vc_sp3[spid] = (*this)[spid].Vc_ * 1.0e+03; //Vc is m3/kmol, Vc_sp3 is cm3/mol
    }
    forAll(Xinp, spid)
    {
        forAll(Xinp, spjd)
        {
            w[spid][spjd] = 0.50 * ((*this)[spid].omega_ + (*this)[spjd].omega_);
            sigma[spid][spjd] = sqrt((0.8090 * pow(Vc_sp3[spid], 1.0 / 3.0)) * (0.8090 * pow(Vc_sp3[spjd], 1.0 / 3.0)));
            //sigma(i,j)= bip_kij(i,j)*sqrt((0.809d0*vc(i)**(1.d0/3.d0))*(0.809d0*vc(j)**(1.d0/3.d0)))
            e_k0[spid][spjd] = sqrt(e_k[spid] * e_k[spjd]);
            //e_k0(i,j)= bip_kij(i,j)*sqrt(e_k(i)*e_k(j))
            //kij[spid][spjd] = sqrt(k[spid] * k[spjd]);
            kij[spid][spjd] = sqrt((*this)[spid].kappa_ * (*this)[spid].kappa_);
            mwij[spid][spjd] = 2.0 * mw_sp3[spid] * mw_sp3[spjd] / (mw_sp3[spid] + mw_sp3[spjd]);
        }
    }

    scalar sigma_mix3 = 0.0;
    scalar sigma_mix2 = 0.0;

    //scalar sigma_mix_23 = sigma_mix2/sigma_mix3;
    scalar e_k_mix = 0.0;
    scalar w_mix = 0.0;
    scalar k_mix = 0.0;
    scalar mu4 = 0.0;

    scalar mw_binary = 0.0;
    forAll(Xinp, spid)
    {
        forAll(Xinp, spjd)
        {
            sigma_mix3 += Xinp[spid] * Xinp[spjd] * pow3(sigma[spid][spjd]);
            sigma_mix2 += Xinp[spid] * Xinp[spjd] * sqr(sigma[spid][spjd]);
            e_k_mix += Xinp[spid] * Xinp[spjd] * e_k0[spid][spjd] * pow3(sigma[spid][spjd]);
        }
    }
    scalar sigma_mix = pow(sigma_mix3, 1.0 / 3);
    scalar sigma_mix2_dev = sigma_mix * sigma_mix;

    forAll(Xinp, spid)
    {
        forAll(Xinp, spjd)
        {

            w_mix += Xinp[spid] * Xinp[spjd] * w[spid][spjd] * pow3(sigma[spid][spjd]);
            mu4 += Xinp[spid] * Xinp[spjd] * sqr((*this)[spid].mu_) * sqr((*this)[spjd].mu_) / (pow3(sigma[spid][spjd]) * e_k0[spid][spjd]);
            k_mix += Xinp[spid] * Xinp[spjd] * kij[spid][spjd];
            mw_binary += (Xinp[spid] * Xinp[spjd] * e_k0[spid][spjd] * sqr(sigma[spid][spjd]) * sqrt(mwij[spid][spjd])); // (e_k_mix * sigma_mix2);
        }
    }
    e_k_mix /= sigma_mix3;
    w_mix /= sigma_mix3;
    mu4 *= sigma_mix3 * e_k_mix;
    mw_binary /= (e_k_mix * sigma_mix2_dev);
    mw_binary = sqr(mw_binary);

    scalar volu_mix = sigma_mix3 / pow3(0.8090); //Eq.(9-5.43);Eq.(9-4.8) bug??
    //scalar volu_mix= pow3(sigma_mix3/0.8090);
    scalar tempc_mix = 1.25930 * e_k_mix;
    //scalar mur_mix4 = pow4(131.30) * mu4 / sqr(tempc_mix * volu_mix); //volu_mix in cm3/mol, mu4 in debyes
    //scalar t_mix = T / e_k_mix;
    //scalar omegav = 1.161450 * (pow(t_mix, -0.148740) + 0.524870 * exp(-0.773200 * t_mix) + 2.161780 * exp(-2.437870 * t_mix)) + (-6.435e-4) * pow(t_mix, 0.148740) * sin(18.0323 * pow(t_mix, -0.7683) - 7.27371); //???GT...

    //1. mixture phases viscosity 40.785 ? 26.69 there is error in 9-5.24
    //viscos = 1e-07 * 26.69 * Fcm * sqrt(mw_binary * T) / (sigma_mix2 * omegav);

    scalar mw_phase = this->W(Xinp);
    //scalar Cv_phase = this->Cp(p, T, pf, Xinp) - this->CpMCv(p, T, pf, Xinp); //J/kgK //warning 
    //scalar Cv_phase = Cv_in; //J/kgK //warning 
    scalar Cv_phase = this->dHidealdT(p,T,Xinp)- RR * 1e-3;
    //scalar Cv_phase = this->Cp(p, T);
    scalar alpha = Cv_phase * mw_phase / (RR * 1.0e-03) - 1.50;
    scalar tempr = T / tempc_mix;

    //compute specific volume, m3/mol
    //scalar amix_phase = this->a(p,T,Xinp), bmix_phase = this->b(p,T,Xinp), dadt_phase = this->dadt(p,T,Xinp), dadt2 = this->d2adt2(p,T,Xinp);
    //scalar sv_phase = this->V(p, T, pf, Xinp); //m3/mol
    scalar rho = rho_in / this->W(Xinp);
    //scalar sv_phase = V_in; //m3/mol

    //2, heat conductivity,W/m.K
    scalar beta = 0.78620 - 0.71090 * w_mix + 1.31680 * sqr(w_mix);
    scalar mw_binary0 = mw_binary / 1.0e+3;
    scalar zeta = 2.00 + 10.50 * sqr(tempr);
    scalar psi = 1.0 + alpha * ((0.2150 + 0.282880 * alpha - 1.0610 * beta + 0.266650 * zeta) / (0.63660 + beta * zeta + 1.0610 * alpha * beta));

    scalar qq = 3.5860 * 1e-03 * sqrt(tempc_mix / mw_binary0) / pow(volu_mix, 2.0 / 3.0);
    //scalar yy = volu_mix / (6.0e+06 * sv_phase);
    scalar yy = rho * volu_mix/ 6.0e+06;
    scalar G1 = (1.0 - 0.50 * yy) / pow3(1.0 - yy);

    // Table 10-3
    //scalar a[7] = { 2.41660, -5.0924e-01, 6.61070, 1.4543e+01, 7.9274e-01, -5.86340, 9.1098e+01 };
    //scalar bb[7] = { 7.4824e-01, -1.50940, 5.62070, -8.91390, 8.2019e-01, 1.2801e+01, 1.2811e+02 };
    //scalar c[7] = { -9.1858e-01, -4.9991e+01, 6.4760e+01, -5.63790, -6.9369e-01, 9.58930, -5.4217e+01 };
    //scalar d[7] = { 1.2712e+02, 6.9983e+01, 2.7039e+01, 7.4344e+01, 6.31730, 6.5529e+01, 5.2381e+02 };
    scalar mur_mix4 = pow4(131.30) * mu4 / sqr(tempc_mix * volu_mix);
    const scalar b0[7] = { 2.41657, -0.50924, 6.61069, 14.54250, 0.79274,
                          -5.8634, 81.171 };
    const scalar b1[7] = { 0.74824, -1.50936, 5.62073, -8.91387, 0.82019,
                          12.8005, 114.1580 };
    const scalar b2[7] = { -0.91858, -49.99120, 64.7599, -5.63794, -0.69369,
                           9.58926, -60.841 };
    const scalar b3[7] = { 121.721, 69.9834, 27.0389, 74.3435, 6.31734,
                          -65.5292, 466.7750 };

    scalar MB[7] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    for (int j = 0; j < 7; j++)
    {
        // Notice that MB[0] = B1, B[1] = B2, ..., B[6] = B7 in eq (13) 
        MB[j] = b0[j] + b1[j] * w_mix + b2[j] * mur_mix4 + b3[j] * k_mix;
    }
    //scalar B_[8] = { 0.0 };
    //for (label ii = 0; ii < 7; ii++)
    //{
    //    B_[ii] = a[ii] + bb[ii] * w_mix + c[ii] * mu4 + d[ii] * k_mix;
    //}

    //scalar G2 = ((B_[1] / yy) * (1.0 - exp(-B_[4] * yy)) + B_[2] * G1 * exp(B_[5] * yy) + B_[3] * G1) / (B_[1] * B_[4] + B_[2] + B_[3]);
    //Info<< "G2* * * *"<<G2<<"* * "<< yy <<"* * "<< volu_mix <<"**"<<sv_phase << endl;
    //Info<< "comp_phase* * * *"<<comp_phase[0] <<"* * "<< comp_phase[1] << endl;

    scalar t_mix = T / e_k_mix;
    scalar Fcm = 1.0 - 0.2756 * w_mix + 0.0590350 * mur_mix4 + k_mix;
    const scalar A = 1.16145, B = 0.14874, C = 0.52487, D = 0.77320, E = 2.16178;
    const scalar F = 2.43787, G = -6.435e-4, H = 7.27371, S = 18.0323, W = -0.76830;

    scalar omegav = A * pow(t_mix, -B) + C * (exp(-D * t_mix)) + E * exp(-F * t_mix)
        + G * pow(t_mix, B) * sin(S * pow(t_mix, W) - H);
    scalar VcM = pow(sigma_mix / 0.809, 3); //✓
    scalar eta0 = (4.0785e-5) * sqrt(mw_binary * T) * Fcm / (omegav * pow(VcM, 2.0 / 3.0)); //unit P ✓
    scalar lamda0 = 7.452 * (eta0 / mw_binary) * psi;

    //scalar rho = rho_in / this->W(Xinp);

    scalar Y = 1e-6 * rho * VcM / 6.0;
    scalar H2 = (MB[0] * (1.0 - exp(-MB[3] * Y)) / Y + MB[1] * G1 * exp(MB[4] * Y) + MB[2] * G1) /
        (MB[0] * MB[3] + MB[1] + MB[2]);
    scalar lamdak = lamda0 * (1.0 / H2 + MB[5] * Y);
    scalar Tr = T / tempc_mix;
    scalar lamdap = (3.039e-4 * pow(tempc_mix / mw_binary, 1.0 / 2.0) / pow(VcM, 2.0 / 3.0)) * MB[6] * pow(Y, 2.0) * H2 * pow(Tr, 1.0 / 2.0);
    //return 31.20 * mu(p, T, Xinp) * psi * (1.0 / G2 + B_[6] * yy) / mw_binary0 + qq * B_[7] * sqr(yy) * sqrt(tempr) * G2;
    scalar kappa = 418.6798 * (lamdak + lamdap);
    return kappa;
    //printf("kappa_phase_thcond_phase %5.8f\n\n",thcond_phase);
}

template <class ThermoMixture>
inline Foam::scalar Foam::chungTransportMixture<ThermoMixture>::alphah(
    const scalar p,
    const scalar T) const
{
    return kappa(p, T, 1) / this->Cp(p, T);
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //
template <class Specie>
Foam::chungTransportMixture<Specie> Foam::operator==(
    const chungTransportMixture<Specie>& a,
    const chungTransportMixture<Specie>& b)
{
    return chungTransportMixture<Specie>(a);
}

template <class Specie>
Foam::chungTransportMixture<Specie> Foam::operator*(
    const scalar s,
    const chungTransportMixture<Specie>& b)
{
    return chungTransportMixture<Specie>(b);
}

template <class Specie>
inline void Foam::chungTransportMixture<Specie>::operator+=(const chungTransportMixture<Specie>& pg)
{
    return;
}

// ************************************************************************* //
