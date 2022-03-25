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

#include "IFstream.H"
#include "IOmanip.H"
#include "OFstream.H"
#include "VLE.H"
#include "labelList.H"
#include "mathematicalConstants.H"
#include "scalarMatrices.H"
#include "specie.H"
#include "string.H"
#include <random>

//#define TPN TPn_flash
//#define TPN TPn_flash_New_TPD

// naming, G:gas, L:liq, xxx_ value,
// * * * * * * * * * * * * * * * * Basic properties  * * * * * * * * * * * * * *
// //

// volume fracion of gas phase
template<class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::alpha(scalar p, scalar T, solution& sol) const
{
    scalar sv_G_ = this->V(p, T, 1, &sol.X_gas()); // gas phase specific volume
    scalar sv_L_ = this->V(p, T, 0, &sol.X_liq()); // liquid phase specific volume

    scalar alpha_G_ = sol.vaporfra * sv_G_ /
        (sol.vaporfra * sv_G_ + (1.0 - sol.vaporfra) * sv_L_);
    if (alpha_G_ > 1.0) // 0.9999999)
    {
        FatalErrorInFunction << "alpha_G_=" << alpha_G_ << exit(FatalError);
        alpha_G_ = 0.9999999;
    }
    else if (alpha_G_ < 0.0) // 1.0e-07)
    {
        FatalErrorInFunction << "alpha_G_=" << alpha_G_ << exit(FatalError);
        alpha_G_ = 1.0e-07;
    }

    return alpha_G_;
}

//-mixing rule
// density kg/m^3
template<class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::rho(scalar p, scalar T) const
{
    // autoPtr<solution> sol(TPn_flash(p, T));
    scalar rho_M_ = p * this->W() / (Z(p, T) * RR * 1.0e-03 * T);
    return rho_M_;
}

// density kg/m^3
template<class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::rho(scalar p, scalar T, solution& sol)
const // TODO test rho!!! and drhodX
{
    // scalar alpha = this->alpha(p, T, sol);

    // scalar ZGas = this->ThermoMixture::Z(p, T, 1, &(sol.X_gas()));
    // scalar ZLiq = this->ThermoMixture::Z(p, T, 0, &(sol.X_liq()));
    // scalar ZMixture = ZGas * alpha + ZLiq * (1.0 - alpha);
    // scalar rhoMixture = p * this->W() / (ZMixture * RR * 1.0e-03 * T);
    // return rhoMixture;

    scalar rho_M_ = p * this->W() / (Z(p, T, sol) * RR * 1.0e-03 * T);
    return rho_M_;
}

// compressiblity factor
template<class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::Z(scalar p, scalar T) const
{
    autoPtr<solution> sol(TPN(p, T));
    return Z(p, T, sol());
}
// compressiblity factor
template<class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::Z(scalar p, scalar T, solution& sol) const
{
    scalar vf_ = sol.vaporfra;
    scalar Z_G_ = this->ThermoMixture::Z(p, T, 1, &(sol.X_gas()));
    scalar Z_L_ = this->ThermoMixture::Z(p, T, 0, &(sol.X_liq()));
    // return (vf * ZGas * ZGas + (1 - vf) * ZLiq * ZLiq) / (vf * ZGas + (1 - vf)
    // * ZLiq);
    return vf_ * Z_G_ + (1 - vf_) * Z_L_;
}

template<class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::psi(scalar p, scalar T) const
{
    return rho(p, T) / p;
}

template<class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::psi(scalar p, scalar T, autoPtr<solution>& sol) const
{
    return rho(p, T, sol()) / p;
}
template<class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::vaporfra(scalar p, scalar T) const
{
    // auto aa=TPN(p, T);
    autoPtr<solution> sol(TPN(p, T));
    // FatalErrorInFunction
    //    << "DEBUG:"
    //    << exit(FatalError);
    // return aa().vaporfra;
    return sol().vaporfra;
}

template<class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::vaporfra(scalar p,
    scalar T,
    autoPtr<solution>& sol) const
{
    return sol().vaporfra;
}

// S might be wrong
template<class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::S(scalar p, scalar T) const
{
    autoPtr<solution> sol(TPN(p, T));
    // label pfg = 1, pfl = 0;
    scalar mw_gas = this->W(&(sol().X_gas())); // kg/mol
    scalar mw_mixture = this->W();
    scalar ygas = sol().vaporfra * mw_gas / mw_mixture;
    return this->ThermoMixture::S(p, T, 1, &(sol().X_gas())) * ygas +
        this->ThermoMixture::S(p, T, 0, &(sol().X_liq())) * (1.0 - ygas);
}

// * * * * * * * * * * * * * * * * Partial Derivative  * * * * * * * * * * * * *
// * //

// *** system (T,P,X) *** //
// X is mole number, but sum(X)=1

// vi: In mixture, mole fraction of specie i, gas phase.
template<class ThermoMixture>
inline autoPtr<scalarList>
Foam::VLE<ThermoMixture>::dvidT(scalar p, scalar T, solution& sol) const
{
    scalarSquareMatrix A(this->X_.size(), Zero);
    scalarList& b = *(new scalarList(this->X_.size()));
    scalar vf_ = sol.vaporfra;
    autoPtr<scalarList> ddT_Ln_fugcoef_G_, ddT_Ln_fugcoef_L_;
    List<autoPtr<scalarList>> ddxi_Ln_fugcoef_G_(this->X_.size()),
        ddxi_Ln_fugcoef_L_(this->X_.size());
    ddT_Ln_fugcoef_G_.reset(
        (this->ddT_Ln_fugacityCoefficient(p, T, 1, &sol.X_gas())).ptr());
    for (int i = 0; i < this->X_.size(); i++)
        ddxi_Ln_fugcoef_G_[i].reset(
            (this->ddxi_Ln_fugacityCoefficient(p, T, i, 1, &sol.X_gas())).ptr());
    ddT_Ln_fugcoef_L_.reset(
        (this->ddT_Ln_fugacityCoefficient(p, T, 0, &sol.X_liq())).ptr());
    for (int i = 0; i < this->X_.size(); i++)
        ddxi_Ln_fugcoef_L_[i].reset(
            (this->ddxi_Ln_fugacityCoefficient(p, T, i, 0, &sol.X_liq())).ptr());
    for (int i = 0; i < this->X_.size(); i++) {
        for (int j = 0; j < this->X_.size(); j++) {
            A(i, j) = 1;
            if (i == j)
                A(i, j) -= this->X_[i] / (sol.X_gas()[i] * sol.X_liq()[i]);
            for (int k = 0; k < this->X_.size(); k++) {
                if (j == k)
                    A(i, j) +=
                    (1 - vf_) * ddxi_Ln_fugcoef_G_[k]()[i] * (sol.X_gas()[k] - 1);
                else
                    A(i, j) += (1 - vf_) * ddxi_Ln_fugcoef_G_[k]()[i] * sol.X_gas()[k];
                if (j == k)
                    A(i, j) += vf_ * ddxi_Ln_fugcoef_L_[k]()[i] * (sol.X_liq()[k] - 1);
                else
                    A(i, j) += vf_ * ddxi_Ln_fugcoef_L_[k]()[i] * sol.X_liq()[k];
            }
        }
    }
    for (int i = 0; i < this->X_.size(); i++)
        b[i] = (1 - vf_) * vf_ * (ddT_Ln_fugcoef_G_()[i] - ddT_Ln_fugcoef_L_()[i]);
    LUsolve(A, b);
    return autoPtr<scalarList>(&b);
}
template<class ThermoMixture>
inline autoPtr<scalarList>
Foam::VLE<ThermoMixture>::dvidP(scalar p, scalar T, solution& sol) const
{
    scalarSquareMatrix A(this->X_.size(), Zero);
    scalarList& b = *(new scalarList(this->X_.size()));
    scalar vf_ = sol.vaporfra;
    autoPtr<scalarList> ddP_Ln_fugcoef_G_, ddP_Ln_fugcoef_L_;
    List<autoPtr<scalarList>> ddxi_Ln_fugcoef_G_(this->X_.size()),
        ddxi_Ln_fugcoef_L_(this->X_.size());

    ddP_Ln_fugcoef_G_.reset(
        (this->ddP_Ln_fugacityCoefficient(p, T, 1, &sol.X_gas())).ptr());

    for (int i = 0; i < this->X_.size(); i++)
        ddxi_Ln_fugcoef_G_[i].reset(
            (this->ddxi_Ln_fugacityCoefficient(p, T, i, 1, &sol.X_gas())).ptr());

    ddP_Ln_fugcoef_L_.reset(
        (this->ddP_Ln_fugacityCoefficient(p, T, 0, &sol.X_liq())).ptr());

    for (int i = 0; i < this->X_.size(); i++)
        ddxi_Ln_fugcoef_L_[i].reset(
            (this->ddxi_Ln_fugacityCoefficient(p, T, i, 0, &sol.X_liq())).ptr());
    for (int i = 0; i < this->X_.size(); i++) {
        for (int j = 0; j < this->X_.size(); j++) {
            A(i, j) = 1;
            if (i == j)
                A(i, j) -= this->X_[i] / (sol.X_gas()[i] * sol.X_liq()[i]);
            for (int k = 0; k < this->X_.size(); k++) {
                if (j == k)
                    A(i, j) +=
                    (1 - vf_) * ddxi_Ln_fugcoef_G_[k]()[i] * (sol.X_gas()[k] - 1);
                else
                    A(i, j) += (1 - vf_) * ddxi_Ln_fugcoef_G_[k]()[i] * sol.X_gas()[k];
                if (j == k)
                    A(i, j) += vf_ * ddxi_Ln_fugcoef_L_[k]()[i] * (sol.X_liq()[k] - 1);
                else
                    A(i, j) += vf_ * ddxi_Ln_fugcoef_L_[k]()[i] * sol.X_liq()[k];
            }
        }
    }
    for (int i = 0; i < this->X_.size(); i++)
        b[i] = (1 - vf_) * vf_ * (ddP_Ln_fugcoef_G_()[i] - ddP_Ln_fugcoef_L_()[i]);
    LUsolve(A, b);
    return autoPtr<scalarList>(&b);
}
template<class ThermoMixture>
inline autoPtr<scalarList>
Foam::VLE<ThermoMixture>::dvidXi(scalar p,
    scalar T,
    label di,
    solution& sol) const
{
    scalarSquareMatrix A(this->X_.size(), Zero);
    scalarList& b = *(new scalarList(this->X_.size()));

    scalar vf_ = sol.vaporfra;
    if (vf_ >= 1-1e-6 || vf_ <= 1e-6)
    {
        for (int k = 0; k < this->X_.size(); k++)
        {
            b[k] = 0;
        }
        if (vf_ >= 1-1e-6)
        {
            b[di] = 1.0;
        }
        return autoPtr<scalarList>(&b);
    }
    autoPtr<scalarList> ddT_Ln_fugcoef_G_, ddT_Ln_fugcoef_L_;
    List<autoPtr<scalarList>> ddxi_Ln_fugcoef_G_(this->X_.size()),
        ddxi_Ln_fugcoef_L_(this->X_.size());

    ddT_Ln_fugcoef_G_.reset(
        (this->ddT_Ln_fugacityCoefficient(p, T, 1, &sol.X_gas())).ptr());

    for (int i = 0; i < this->X_.size(); i++)
        ddxi_Ln_fugcoef_G_[i].reset(
            (this->ddxi_Ln_fugacityCoefficient(p, T, i, 1, &sol.X_gas())).ptr());

    ddT_Ln_fugcoef_L_.reset(
        (this->ddT_Ln_fugacityCoefficient(p, T, 0, &sol.X_liq())).ptr());

    for (int i = 0; i < this->X_.size(); i++)
        ddxi_Ln_fugcoef_L_[i].reset(
            (this->ddxi_Ln_fugacityCoefficient(p, T, i, 0, &sol.X_liq())).ptr());
    for (int i = 0; i < this->X_.size(); i++) {
        for (int j = 0; j < this->X_.size(); j++) {
            A(i, j) = -1.0 / (vf_ * (1 - vf_));
            if (i == j)
                A(i, j) +=
                1 / (sol.X_liq()[i] * (1 - vf_)) + 1 / (sol.X_gas()[i] * vf_);
            for (int k = 0; k < this->X_.size(); k++) {
                if (j == k)
                    A(i, j) += ddxi_Ln_fugcoef_G_[k]()[i] * (1 - sol.X_gas()[k]) / vf_;
                else
                    A(i, j) += ddxi_Ln_fugcoef_G_[k]()[i] * (0 - sol.X_gas()[k]) / vf_;
                if (j == k)
                    A(i, j) +=
                    ddxi_Ln_fugcoef_L_[k]()[i] * (1 - sol.X_liq()[k]) / (1 - vf_);
                else
                    A(i, j) +=
                    ddxi_Ln_fugcoef_L_[k]()[i] * (0 - sol.X_liq()[k]) / (1 - vf_);
            }
        }
    }
    for (int k = 0; k < this->X_.size(); k++) {
        b[k] = 0;
        for (int j = 0; j < this->X_.size(); j++) {
            if (j == di) {
                for (int s = 0; s < this->X_.size(); s++) {
                    if (s == j)
                        b[k] +=
                        ddxi_Ln_fugcoef_L_[s]()[k] * (1 - sol.X_liq()[s]) / (1 - vf_);
                    else
                        b[k] +=
                        ddxi_Ln_fugcoef_L_[s]()[k] * (0 - sol.X_liq()[s]) / (1 - vf_);
                }

                if (k == j) {
                    b[k] += 1 / sol.X_liq()[k] * (1 - sol.X_liq()[k]) / (1 - vf_);
                }
                else {
                    b[k] += 1 / sol.X_liq()[k] * (0 - sol.X_liq()[k]) / (1 - vf_);
                }
            }
        }
    }
    LUsolve(A, b);
    return autoPtr<scalarList>(&b);
}

// vf
template<class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::dvfdT(scalar p, scalar T, solution& sol) const
{
    autoPtr<scalarList> dvidT_(this->dvidT(p, T, sol));
    scalar dvfdT_ = 0;
    for (int i = 0; i < this->X_.size(); i++) {
        dvfdT_ += dvidT_()[i];
    }
    return dvfdT_;
}
template<class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::dvfdP(scalar p, scalar T, solution& sol) const
{
    autoPtr<scalarList> dvidP_(this->dvidP(p, T, sol));

    scalar dvfdP_ = 0;
    for (int i = 0; i < this->X_.size(); i++) {
        dvfdP_ += dvidP_()[i];
    }
    return dvfdP_;
}
template<class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::dvfdXi(scalar p,
    scalar T,
    label di,
    solution& sol) const
{
    autoPtr<scalarList> dvidXi_(this->dvidXi(p, T, di, sol));
    scalar dvfdXi_ = -sol.vaporfra;
    for (int i = 0; i < this->X_.size(); i++) {
        dvfdXi_ += dvidXi_()[i];
    }
    return dvfdXi_;
}

// Z mixture
template<class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::dZdT(scalar p, scalar T, solution& sol) const
{
    autoPtr<scalarList> dXdT_G_(this->dXdT_gas(p, T, sol));
    autoPtr<scalarList> dXdT_L_(this->dXdT_liq(p, T, sol));
    scalar dvfdT_ = this->dvfdT(p, T, sol);
    scalar vf_ = sol.vaporfra;
    scalar Z_G_ = this->ThermoMixture::Z(p, T, 1, &(sol.X_gas()));
    scalar Z_L_ = this->ThermoMixture::Z(p, T, 0, &(sol.X_liq()));
    scalar dZdT_G_ = this->ThermoMixture::dZdT(p, T, 1, &(sol.X_gas()));
    scalar dZdT_L_ = this->ThermoMixture::dZdT(p, T, 0, &(sol.X_liq()));
    for (int i = 0; i < this->X_.size(); i++) {
        dZdT_G_ +=
            this->ThermoMixture::dZdxi(p, T, i, 1, &(sol.X_gas())) * dXdT_G_()[i];
        dZdT_L_ +=
            this->ThermoMixture::dZdxi(p, T, i, 0, &(sol.X_liq())) * dXdT_L_()[i];
    }
    // scalar temp = (vf * ZGas + (1 - vf) * ZLiq);
    // return (dvfdT * ZGas * ZGas + vf * 2 * ZGas * dZgasdT + 2 * (1 - vf) * ZLiq
    // * dZliqdT - dvfdT * ZLiq * ZLiq) / temp - (vf * ZGas * ZGas + (1 - vf) *
    // ZLiq * ZLiq) / (temp * temp) * (dvfdT * ZGas + vf * dZgasdT + (1 - vf) *
    // dZliqdT - dvfdT * ZLiq);
    return dvfdT_ * (Z_G_ - Z_L_) + vf_ * dZdT_G_ + (1 - vf_) * dZdT_L_;
}
template<class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::dZdP(scalar p, scalar T, solution& sol) const
{
    autoPtr<scalarList> dXdP_G_(this->dXdP_gas(p, T, sol));
    autoPtr<scalarList> dXdP_L_(this->dXdP_liq(p, T, sol));
    scalar dvfdP_ = this->dvfdP(p, T, sol);
    scalar vf_ = sol.vaporfra;
    scalar Z_G_ = this->ThermoMixture::Z(p, T, 1, &(sol.X_gas()));
    scalar Z_L_ = this->ThermoMixture::Z(p, T, 0, &(sol.X_liq()));
    scalar dZdP_G_ = this->ThermoMixture::dZdP(p, T, 1, &(sol.X_gas()));
    scalar dZdP_L_ = this->ThermoMixture::dZdP(p, T, 0, &(sol.X_liq()));
    for (int i = 0; i < this->X_.size(); i++) {
        dZdP_G_ +=
            this->ThermoMixture::dZdxi(p, T, i, 1, &(sol.X_gas())) * dXdP_G_()[i];
        dZdP_L_ +=
            this->ThermoMixture::dZdxi(p, T, i, 0, &(sol.X_liq())) * dXdP_L_()[i];
    }
    // scalar temp = (vf * ZGas + (1 - vf) * ZLiq);
    // return (dvfdt * ZGas * ZGas + vf * 2 * ZGas * dZgasdT + 2 * (1 - vf) * ZLiq
    // * dZliqdT - dvfdt * ZLiq * ZLiq) / temp - (vf * ZGas * ZGas + (1 - vf) *
    // ZLiq * ZLiq) / (temp * temp) * (dvfdt * ZGas + vf * dZgasdT + (1 - vf) *
    // dZliqdT - dvfdt * ZLiq);
    return dvfdP_ * (Z_G_ - Z_L_) + vf_ * dZdP_G_ + (1 - vf_) * dZdP_L_;
}
template<class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::dZdXi(scalar p,
    scalar T,
    label di,
    solution& sol) const
{
    autoPtr<scalarList> dXdXi_G_(this->dXdXi_gas(p, T, di, sol));
    autoPtr<scalarList> dXdXi_L_(this->dXdXi_liq(p, T, di, sol));
    scalar dvfdXi_ = this->dvfdXi(p, T, di, sol);
    scalar vf_ = sol.vaporfra;
    scalar Z_G_ = this->ThermoMixture::Z(p, T, 1, &(sol.X_gas()));
    scalar Z_L_ = this->ThermoMixture::Z(p, T, 0, &(sol.X_liq()));
    scalar dZdXi_G_ = 0;
    scalar dZdXi_L_ = 0;
    for (int i = 0; i < this->X_.size(); i++) {
        dZdXi_G_ +=
            this->ThermoMixture::dZdxi(p, T, i, 1, &(sol.X_gas())) * dXdXi_G_()[i];
        dZdXi_L_ +=
            this->ThermoMixture::dZdxi(p, T, i, 0, &(sol.X_liq())) * dXdXi_L_()[i];
    }
    // scalar temp = (vf * ZGas + (1 - vf) * ZLiq);
    // return (dvfdt * ZGas * ZGas + vf * 2 * ZGas * dZgasdT + 2 * (1 - vf) * ZLiq
    // * dZliqdT - dvfdt * ZLiq * ZLiq) / temp - (vf * ZGas * ZGas + (1 - vf) *
    // ZLiq * ZLiq) / (temp * temp) * (dvfdt * ZGas + vf * dZgasdT + (1 - vf) *
    // dZliqdT - dvfdt * ZLiq);
    return dvfdXi_ * (Z_G_ - Z_L_) + vf_ * dZdXi_G_ + (1 - vf_) * dZdXi_L_;
}

// rho
template<class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::drhodT(scalar p, scalar T, solution& sol) const
{
    scalar Z_M_ = Z(p, T, sol);
    scalar dZdT_M_ = dZdT(p, T, sol);
    scalar drhodT_M_ = -p * this->W() / (Z_M_ * Z_M_ * RR * 1.0e-03 * T * T) *
        (dZdT_M_ * T + Z_M_);
    return drhodT_M_;
}
template<class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::drhodP(scalar p, scalar T, solution& sol) const
{
    scalar Z_M_ = Z(p, T, sol);
    scalar dZdP_M_ = dZdP(p, T, sol);
    scalar drhodP_M_ = this->W() / (Z_M_ * RR * 1.0e-03 * T) -
        p * this->W() / (Z_M_ * Z_M_ * RR * 1.0e-03 * T) * dZdP_M_;
    return drhodP_M_;
}
template<class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::drhodXi(scalar p,
    scalar T,
    label di,
    solution& sol) const
{
    scalar Z_M_ = Z(p, T, sol);
    scalar dZdXi_M_ = dZdXi(p, T, di, sol);
    scalar drhodXi_M_ =
        p * this->dWdXi(di) / (Z_M_ * RR * 1.0e-03 * T) -
        p * this->W() / (Z_M_ * Z_M_ * RR * 1.0e-03 * T) * dZdXi_M_;
    return drhodXi_M_;
}

// Es
template<class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::dEdT(scalar p, scalar T, solution& sol) const
{
    scalar rho_M_ = rho(p, T, sol);
    return Cp_Hs(p, T, sol) + p / (rho_M_ * rho_M_) * drhodT(p, T, sol);
}

template<class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::dEdP(scalar p, scalar T, solution& sol) const
{
    scalar rho_M_ = rho(p, T, sol);
    return dHsdP(p, T, sol) - 1 / rho_M_ +
        p / (rho_M_ * rho_M_) * drhodP(p, T, sol);
}

template<class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::dEdXi(scalar p, scalar T, label di, solution& sol) const
{
    scalar rho_M_ = rho(p, T, sol);
    return dHsdXi(p, T, di, sol) + p / (rho_M_ * rho_M_) * drhodXi(p, T, di, sol);
}

// *** system (E,rho,X) *** //

template<class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::dTdE_rhoX(scalar p, scalar T, solution& sol) const
{
    scalar dEdT_ = dEdT(p, T, sol);
    scalar dEdP_ = dEdP(p, T, sol);
    scalar drhodT_ = drhodT(p, T, sol);
    scalar drhodP_ = drhodP(p, T, sol);
    scalar det = dEdT_ * drhodP_ - dEdP_ * drhodT_;
    return drhodP_ / det;
}

template<class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::dTdrho_EX(scalar p, scalar T, solution& sol) const
{
    scalar dEdT_ = dEdT(p, T, sol);
    scalar dEdP_ = dEdP(p, T, sol);
    scalar drhodT_ = drhodT(p, T, sol);
    scalar drhodP_ = drhodP(p, T, sol);
    scalar det = dEdT_ * drhodP_ - dEdP_ * drhodT_;
    return -dEdP_ / det;
}

template<class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::dPdE_rhoX(scalar p, scalar T, solution& sol) const
{
    scalar dEdT_ = dEdT(p, T, sol);
    scalar dEdP_ = dEdP(p, T, sol);
    scalar drhodT_ = drhodT(p, T, sol);
    scalar drhodP_ = drhodP(p, T, sol);
    scalar det = dEdT_ * drhodP_ - dEdP_ * drhodT_;
    return -drhodT_ / det;
}

template<class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::dPdrho_EX(scalar p, scalar T, solution& sol) const
{
    scalar dEdT_ = dEdT(p, T, sol);
    scalar dEdP_ = dEdP(p, T, sol);
    scalar drhodT_ = drhodT(p, T, sol);
    scalar drhodP_ = drhodP(p, T, sol);
    scalar det = dEdT_ * drhodP_ - dEdP_ * drhodT_;
    return dEdT_ / det;
}


template<class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::dTdXi_Erho(scalar p, scalar T, label di, solution& sol) const
{
    scalar dEdT_ = dEdT(p, T, sol);
    scalar dEdP_ = dEdP(p, T, sol);
    scalar drhodT_ = drhodT(p, T, sol);
    scalar drhodP_ = drhodP(p, T, sol);
    scalar dEdXi_ = dEdXi(p, T, di, sol);
    scalar drhodXi_ = drhodXi(p, T, di, sol);
    scalar det = dEdT_ * drhodP_ - dEdP_ * drhodT_;
    return -(drhodP_ * dEdXi_ - dEdP_ * drhodXi_) / det;
}
template<class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::dPdXi_Erho(scalar p, scalar T, label di, solution& sol) const
{
    scalar dEdT_ = dEdT(p, T, sol);
    scalar dEdP_ = dEdP(p, T, sol);
    scalar drhodT_ = drhodT(p, T, sol);
    scalar drhodP_ = drhodP(p, T, sol);
    scalar dEdXi_ = dEdXi(p, T, di, sol);
    scalar drhodXi_ = drhodXi(p, T, di, sol);
    scalar det = dEdT_ * drhodP_ - dEdP_ * drhodT_;
    return -(-drhodT_ * dEdXi_ + dEdT_ * drhodXi_) / det;
}
template<class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::dvfdE_rhoX(scalar p, scalar T, solution& sol) const
{
    return dvfdT(p, T, sol) * dTdE_rhoX(p, T, sol) + dvfdP(p, T, sol) * dPdE_rhoX(p, T, sol);
}
template<class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::dvfdrho_EX(scalar p, scalar T, solution& sol) const
{
    return dvfdT(p, T, sol) * dTdrho_EX(p, T, sol) + dvfdP(p, T, sol) * dPdrho_EX(p, T, sol);
}
template<class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::dvfdXi_Erho(scalar p, scalar T, label di, solution& sol) const
{
    return dvfdXi(p, T, di, sol) + dvfdT(p, T, sol) * dTdXi_Erho(p, T, di, sol) + dvfdP(p, T, sol) * dPdXi_Erho(p, T, di, sol);
}
// * * * * * * * * * * * * * * * *  Derivative  properties * * * * * * * * * * *
// * * * //

// isothermal compressibility
template<class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::kappaT(scalar p, scalar T, solution& sol) const
{
    scalar rho_ = rho(p, T, sol);
    scalar drhodP_ = drhodP(p, T, sol);
    return drhodP_ / rho_;
}

// isobaric expansivity
template<class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::alphaP(scalar p, scalar T, solution& sol) const
{
    scalar rho_ = rho(p, T, sol);
    scalar drhodT_ = drhodT(p, T, sol);
    return -drhodT_ / rho_;
}
// isentropic compressibility
template<class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::kappaS(scalar p, scalar T, solution& sol) const
{
    scalar rho_ = rho(p, T, sol);
    scalar kappaT_ = kappaT(p, T, sol);
    scalar Cp_ = Cp(p, T, sol);
    scalar alphaP_ = alphaP(p, T, sol);
    return kappaT_ - T * alphaP_ * alphaP_ / (rho_ * Cp_);
}
// speed of sound
template<class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::c(scalar p, scalar T, solution& sol) const
{
    scalar rho_ = rho(p, T, sol);
    scalar kappaS_ = kappaS(p, T, sol);
    return sqrt(1.0 / (kappaS_ * rho_));
}
template<class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::c(scalar p, scalar T) const
{
    autoPtr<solution> sol(TPN(p, T));
    return c(p, T, sol());
}
// ************************************************************************* //
