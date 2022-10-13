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
template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::alpha(scalar p, scalar T, solution &sol) const
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
template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::rho(scalar p, scalar T) const
{
    // autoPtr<solution> sol(TPn_flash(p, T));

    scalar rho_M_ = p * this->W() / (Z(p, T) * RR * 1.0e-03 * T);
    return rho_M_;
}

template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::rrho(scalar p, scalar T) const
{
    // autoPtr<solution> sol(TPn_flash(p, T));

    scalar rrho_M_ = (Z(p, T) * RR * 1.0e-03 * T) / (p * this->W());
    return rrho_M_;
}
template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::rrho(scalar p, scalar T, solution &sol) const
{
    // autoPtr<solution> sol(TPn_flash(p, T));
    scalar Z_ = Z(p, T, sol);
    scalar rrho_M_ = (Z_ * RR * 1.0e-03 * T) / (p * this->W());
    return rrho_M_;
}

// density kg/m^3
template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::rho(scalar p, scalar T, solution &sol)
    const // TODO test rho!!! and drhodX
{
    // scalar alpha = this->alpha(p, T, sol);

    // scalar ZGas = this->ThermoMixture::Z(p, T, 1, &(sol.X_gas()));
    // scalar ZLiq = this->ThermoMixture::Z(p, T, 0, &(sol.X_liq()));
    // scalar ZMixture = ZGas * alpha + ZLiq * (1.0 - alpha);
    // scalar rhoMixture = p * this->W() / (ZMixture * RR * 1.0e-03 * T);
    // return rhoMixture;
    scalar Z_ = Z(p, T, sol);
    scalar rho_M_ = p * this->W() / (Z_ * RR * 1.0e-03 * T);
    return rho_M_;
}

template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::rho_noVLE(scalar p, scalar T, solution &sol)
    const // TODO test rho!!! and drhodX
{

    scalar Z_ = this->ThermoMixture::Z(p, T, sol.X_gas(), 1);
    scalar rho_M_ = p * this->W() / (Z_ * RR * 1.0e-03 * T);
    return rho_M_;
}

// compressiblity factor
template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::Z(scalar p, scalar T) const
{
    autoPtr<solution> sol(TPN(p, T));
    return Z(p, T, sol());
}
// compressiblity factor
template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::Z(scalar p, scalar T, solution &sol) const
{
    scalar vf_ = sol.vaporfra;
    scalar Z_G_ = this->ThermoMixture::Z(p, T, sol.X_gas(), 1);
    scalar Z_L_ = this->ThermoMixture::Z(p, T, sol.X_liq(), 0);
    // return (vf * ZGas * ZGas + (1 - vf) * ZLiq * ZLiq) / (vf * ZGas + (1 - vf)
    // * ZLiq);
    return vf_ * Z_G_ + (1 - vf_) * Z_L_;
}

template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::psi(scalar p, scalar T) const
{
    return rho(p, T) / p;
}

template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::psi(scalar p, scalar T, autoPtr<solution> &sol) const
{
    return rho(p, T, sol()) / p;
}
template <class ThermoMixture>
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

template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::vaporfra(scalar p,
                                   scalar T,
                                   autoPtr<solution> &sol) const
{
    return sol().vaporfra;
}

// S might be wrong
template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::S(scalar p, scalar T) const
{
    autoPtr<solution> sol(TPN(p, T));
    // label pfg = 1, pfl = 0;
    scalar mw_gas = this->W(sol().X_gas()); // kg/mol
    scalar mw_mixture = this->W();
    scalar ygas = sol().vaporfra * mw_gas / mw_mixture;
    return this->ThermoMixture::S(p, T, sol().X_gas(), 1) * ygas +
           this->ThermoMixture::S(p, T, sol().X_liq(), 0) * (1.0 - ygas);
}

// * * * * * * * * * * * * * * * * Partial Derivative  * * * * * * * * * * * * *
// * //

// *** system (T,P,X) *** //
// X is mole number, but sum(X)=1

template <class ThermoMixture>
inline autoPtr<scalarSquareMatrix>
Foam::VLE<ThermoMixture>::A_dvid_(scalar p, scalar T, solution &sol) const
{
    scalarSquareMatrix &A = *(new scalarSquareMatrix(this->X_.size(), Zero));
    scalar vf_ = sol.vaporfra;

    List<autoPtr<scalarList>> ddxi_Ln_fugcoef_G_(this->X_.size()),
        ddxi_Ln_fugcoef_L_(this->X_.size());

    for (int i = 0; i < this->X_.size(); i++)
        ddxi_Ln_fugcoef_G_[i].reset(
            (this->ddxi_Ln_fugacityCoefficient(p, T, i, sol.X_gas(), 1)).ptr());

    for (int i = 0; i < this->X_.size(); i++)
        ddxi_Ln_fugcoef_L_[i].reset(
            (this->ddxi_Ln_fugacityCoefficient(p, T, i, sol.X_liq(), 0)).ptr());
    for (int i = 0; i < this->X_.size(); i++)
    {
        for (int j = 0; j < this->X_.size(); j++)
        {
            A(i, j) = 1;
            if (i == j)
                A(i, j) -= this->X_[i] / (sol.X_gas()[i] * sol.X_liq()[i]);
            for (int k = 0; k < this->X_.size(); k++)
            {
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

    return autoPtr<scalarSquareMatrix>(&A);
}

// vi: In mixture, mole fraction of specie i, gas phase.
template <class ThermoMixture>
inline autoPtr<scalarList>
Foam::VLE<ThermoMixture>::dvidT(scalar p, scalar T, solution &sol) const
{
    scalarSquareMatrix A(this->X_.size(), Zero);
    scalarList &b = *(new scalarList(this->X_.size()));
    scalar vf_ = sol.vaporfra;
    int n_sp = 0;
    for (int i = 0; i < this->X_.size(); i++)
    {
        if (this->X_[i] > 0)
            n_sp++;
    }
    if (n_sp == 1 || vf_ == 0 || vf_ == 1)
    {
        for (int i = 0; i < this->X_.size(); i++)
            b[i] = 0;
        return autoPtr<scalarList>(&b);
    }
    autoPtr<scalarList> ddT_Ln_fugcoef_G_, ddT_Ln_fugcoef_L_;
    List<autoPtr<scalarList>> ddxi_Ln_fugcoef_G_(this->X_.size()),
        ddxi_Ln_fugcoef_L_(this->X_.size());
    ddT_Ln_fugcoef_G_.reset(
        (this->ddT_Ln_fugacityCoefficient(p, T, sol.X_gas(), 1)).ptr());
    for (int i = 0; i < this->X_.size(); i++)
        ddxi_Ln_fugcoef_G_[i].reset(
            (this->ddxi_Ln_fugacityCoefficient(p, T, i, sol.X_gas(), 1)).ptr());
    ddT_Ln_fugcoef_L_.reset(
        (this->ddT_Ln_fugacityCoefficient(p, T, sol.X_liq(), 0)).ptr());
    for (int i = 0; i < this->X_.size(); i++)
        ddxi_Ln_fugcoef_L_[i].reset(
            (this->ddxi_Ln_fugacityCoefficient(p, T, i, sol.X_liq(), 0)).ptr());
    for (int i = 0; i < this->X_.size(); i++)
    {
        for (int j = 0; j < this->X_.size(); j++)
        {
            A(i, j) = 1;
            if (i == j)
                A(i, j) -= this->X_[i] / (sol.X_gas()[i] * sol.X_liq()[i]);
            for (int k = 0; k < this->X_.size(); k++)
            {
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

template <class ThermoMixture>
inline autoPtr<scalarList>
Foam::VLE<ThermoMixture>::dvidT(scalar p, scalar T, solution &sol, scalarSquareMatrix &A, labelList &pivotIndices) const
{
    scalarList &b = *(new scalarList(this->X_.size()));
    scalar vf_ = sol.vaporfra;
    int n_sp = 0;
    for (int i = 0; i < this->X_.size(); i++)
    {
        if (this->X_[i] > 0)
            n_sp++;
    }
    if (n_sp == 1 || vf_ == 0 || vf_ == 1)
    {
        for (int i = 0; i < this->X_.size(); i++)
            b[i] = 0;
        return autoPtr<scalarList>(&b);
    }
    autoPtr<scalarList> ddT_Ln_fugcoef_G_, ddT_Ln_fugcoef_L_;

    ddT_Ln_fugcoef_G_.reset(
        (this->ddT_Ln_fugacityCoefficient(p, T, sol.X_gas(), 1)).ptr());

    ddT_Ln_fugcoef_L_.reset(
        (this->ddT_Ln_fugacityCoefficient(p, T, sol.X_liq(), 0)).ptr());

    for (int i = 0; i < this->X_.size(); i++)
        b[i] = (1 - vf_) * vf_ * (ddT_Ln_fugcoef_G_()[i] - ddT_Ln_fugcoef_L_()[i]);
    //LUsolve(A, b);
    LUBacksubstitute(A, pivotIndices, b);
    return autoPtr<scalarList>(&b);
}
template <class ThermoMixture>
inline autoPtr<scalarList>
Foam::VLE<ThermoMixture>::dvidP(scalar p, scalar T, solution &sol) const
{
    scalarSquareMatrix A(this->X_.size(), Zero);
    scalarList &b = *(new scalarList(this->X_.size()));
    scalar vf_ = sol.vaporfra;
    int n_sp = 0;
    for (int i = 0; i < this->X_.size(); i++)
    {
        if (this->X_[i] > 0)
            n_sp++;
    }
    if (n_sp == 1)
    {
        for (int i = 0; i < this->X_.size(); i++)
            b[i] = 0;
        return autoPtr<scalarList>(&b);
    }
    autoPtr<scalarList> ddP_Ln_fugcoef_G_, ddP_Ln_fugcoef_L_;
    List<autoPtr<scalarList>> ddxi_Ln_fugcoef_G_(this->X_.size()),
        ddxi_Ln_fugcoef_L_(this->X_.size());

    ddP_Ln_fugcoef_G_.reset(
        (this->ddP_Ln_fugacityCoefficient(p, T, sol.X_gas(), 1)).ptr());

    for (int i = 0; i < this->X_.size(); i++)
        ddxi_Ln_fugcoef_G_[i].reset(
            (this->ddxi_Ln_fugacityCoefficient(p, T, i, sol.X_gas(), 1)).ptr());

    ddP_Ln_fugcoef_L_.reset(
        (this->ddP_Ln_fugacityCoefficient(p, T, sol.X_liq(), 0)).ptr());

    for (int i = 0; i < this->X_.size(); i++)
        ddxi_Ln_fugcoef_L_[i].reset(
            (this->ddxi_Ln_fugacityCoefficient(p, T, i, sol.X_liq(), 0)).ptr());
    for (int i = 0; i < this->X_.size(); i++)
    {
        for (int j = 0; j < this->X_.size(); j++)
        {
            A(i, j) = 1;
            if (i == j)
                A(i, j) -= this->X_[i] / (sol.X_gas()[i] * sol.X_liq()[i]);
            for (int k = 0; k < this->X_.size(); k++)
            {
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

template <class ThermoMixture>
inline autoPtr<scalarList>
Foam::VLE<ThermoMixture>::dvidP(scalar p, scalar T, solution &sol, scalarSquareMatrix &A, labelList &pivotIndices) const
{

    scalarList &b = *(new scalarList(this->X_.size()));
    scalar vf_ = sol.vaporfra;
    int n_sp = 0;
    for (int i = 0; i < this->X_.size(); i++)
    {
        if (this->X_[i] > 0)
            n_sp++;
    }
    if (n_sp == 1)
    {
        for (int i = 0; i < this->X_.size(); i++)
            b[i] = 0;
        return autoPtr<scalarList>(&b);
    }
    autoPtr<scalarList> ddP_Ln_fugcoef_G_, ddP_Ln_fugcoef_L_;

    ddP_Ln_fugcoef_G_.reset(
        (this->ddP_Ln_fugacityCoefficient(p, T, sol.X_gas(), 1)).ptr());

    ddP_Ln_fugcoef_L_.reset(
        (this->ddP_Ln_fugacityCoefficient(p, T, sol.X_liq(), 0)).ptr());

    for (int i = 0; i < this->X_.size(); i++)
        b[i] = (1 - vf_) * vf_ * (ddP_Ln_fugcoef_G_()[i] - ddP_Ln_fugcoef_L_()[i]);
    //LUsolve(A, b);
    LUBacksubstitute(A, pivotIndices, b);
    return autoPtr<scalarList>(&b);
}
template <class ThermoMixture>
inline autoPtr<scalarList>
Foam::VLE<ThermoMixture>::dvidXi(scalar p,
                                 scalar T,
                                 label di,
                                 solution &sol) const
{
    scalarSquareMatrix A(this->X_.size(), Zero);
    scalarList &b = *(new scalarList(this->X_.size()));

    scalar vf_ = sol.vaporfra;
    if (vf_ >= 1 - 1e-6 || vf_ <= 1e-6)
    {
        for (int k = 0; k < this->X_.size(); k++)
        {
            b[k] = 0;
        }
        if (vf_ >= 1 - 1e-6)
        {
            b[di] = 1.0;
        }
        return autoPtr<scalarList>(&b);
    }
    autoPtr<scalarList> ddT_Ln_fugcoef_G_, ddT_Ln_fugcoef_L_;
    List<autoPtr<scalarList>> ddxi_Ln_fugcoef_G_(this->X_.size()),
        ddxi_Ln_fugcoef_L_(this->X_.size());

    ddT_Ln_fugcoef_G_.reset(
        (this->ddT_Ln_fugacityCoefficient(p, T, sol.X_gas(), 1)).ptr());

    for (int i = 0; i < this->X_.size(); i++)
        ddxi_Ln_fugcoef_G_[i].reset(
            (this->ddxi_Ln_fugacityCoefficient(p, T, i, sol.X_gas(), 1)).ptr());

    ddT_Ln_fugcoef_L_.reset(
        (this->ddT_Ln_fugacityCoefficient(p, T, sol.X_liq(), 0)).ptr());

    for (int i = 0; i < this->X_.size(); i++)
        ddxi_Ln_fugcoef_L_[i].reset(
            (this->ddxi_Ln_fugacityCoefficient(p, T, i, sol.X_liq(), 0)).ptr());
    for (int i = 0; i < this->X_.size(); i++)
    {
        for (int j = 0; j < this->X_.size(); j++)
        {
            A(i, j) = -1.0 / (vf_ * (1 - vf_));
            if (i == j)
                A(i, j) +=
                    1 / (sol.X_liq()[i] * (1 - vf_)) + 1 / (sol.X_gas()[i] * vf_);
            for (int k = 0; k < this->X_.size(); k++)
            {
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
    for (int k = 0; k < this->X_.size(); k++)
    {
        b[k] = 0;
        for (int j = 0; j < this->X_.size(); j++)
        {
            if (j == di)
            {
                for (int s = 0; s < this->X_.size(); s++)
                {
                    if (s == j)
                        b[k] +=
                            ddxi_Ln_fugcoef_L_[s]()[k] * (1 - sol.X_liq()[s]) / (1 - vf_);
                    else
                        b[k] +=
                            ddxi_Ln_fugcoef_L_[s]()[k] * (0 - sol.X_liq()[s]) / (1 - vf_);
                }

                if (k == j)
                {
                    b[k] += 1 / sol.X_liq()[k] * (1 - sol.X_liq()[k]) / (1 - vf_);
                }
                else
                {
                    b[k] += 1 / sol.X_liq()[k] * (0 - sol.X_liq()[k]) / (1 - vf_);
                }
            }
        }
    }
    LUsolve(A, b);
    return autoPtr<scalarList>(&b);
}

// vf
template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::dvfdT(scalar p, scalar T, solution &sol) const
{
    autoPtr<scalarList> dvidT_(this->dvidT(p, T, sol));
    scalar dvfdT_ = 0;
    for (int i = 0; i < this->X_.size(); i++)
    {
        dvfdT_ += dvidT_()[i];
    }
    return dvfdT_;
}

template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::dvfdT(scalar p, scalar T, solution &sol, scalarList &dvidT_) const
{
    //autoPtr<scalarList> dvidT_(this->dvidT(p, T, sol));
    scalar dvfdT_ = 0;
    for (int i = 0; i < this->X_.size(); i++)
    {
        dvfdT_ += dvidT_[i];
    }
    return dvfdT_;
}
template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::dvfdP(scalar p, scalar T, solution &sol) const
{
    autoPtr<scalarList> dvidP_(this->dvidP(p, T, sol));

    scalar dvfdP_ = 0;
    for (int i = 0; i < this->X_.size(); i++)
    {
        dvfdP_ += dvidP_()[i];
    }
    return dvfdP_;
}

template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::dvfdP(scalar p, scalar T, solution &sol, scalarList &dvidP_) const
{
    //autoPtr<scalarList> dvidP_(this->dvidP(p, T, sol));

    scalar dvfdP_ = 0;
    for (int i = 0; i < this->X_.size(); i++)
    {
        dvfdP_ += dvidP_[i];
    }
    return dvfdP_;
}
template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::dvfdXi(scalar p,
                                 scalar T,
                                 label di,
                                 solution &sol) const
{
    autoPtr<scalarList> dvidXi_(this->dvidXi(p, T, di, sol));
    scalar dvfdXi_ = -sol.vaporfra;
    for (int i = 0; i < this->X_.size(); i++)
    {
        dvfdXi_ += dvidXi_()[i];
    }
    return dvfdXi_;
}

// Z mixture
template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::dZdT(scalar p, scalar T, solution &sol) const
{
    autoPtr<scalarList> dvidt;
    dvidt.reset((this->dvidT(p, T, sol)).ptr());
    autoPtr<scalarList> dXdT_G_(this->dXdT_gas(p, T, sol, dvidt));
    autoPtr<scalarList> dXdT_L_(this->dXdT_liq(p, T, sol, dvidt));
    scalar dvfdT_ = this->dvfdT(p, T, sol, dvidt());
    scalar vf_ = sol.vaporfra;
    scalar Z_G_ = this->ThermoMixture::Z(p, T, sol.X_gas(), 1);
    scalar Z_L_ = this->ThermoMixture::Z(p, T, sol.X_liq(), 0);
    scalar dZdT_G_ = this->ThermoMixture::dZdT(p, T, sol.X_gas(), 1);
    scalar dZdT_L_ = this->ThermoMixture::dZdT(p, T, sol.X_liq(), 0);
    for (int i = 0; i < this->X_.size(); i++)
    {
        dZdT_G_ +=
            this->ThermoMixture::dZdxi(p, T, i, sol.X_gas(), 1) * dXdT_G_()[i];
        dZdT_L_ +=
            this->ThermoMixture::dZdxi(p, T, i, sol.X_liq(), 0) * dXdT_L_()[i];
    }
    // scalar temp = (vf * ZGas + (1 - vf) * ZLiq);
    // return (dvfdT * ZGas * ZGas + vf * 2 * ZGas * dZgasdT + 2 * (1 - vf) * ZLiq
    // * dZliqdT - dvfdT * ZLiq * ZLiq) / temp - (vf * ZGas * ZGas + (1 - vf) *
    // ZLiq * ZLiq) / (temp * temp) * (dvfdT * ZGas + vf * dZgasdT + (1 - vf) *
    // dZliqdT - dvfdT * ZLiq);
    return dvfdT_ * (Z_G_ - Z_L_) + vf_ * dZdT_G_ + (1 - vf_) * dZdT_L_;
}

template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::dZdT(scalar p, scalar T, solution &sol, scalarSquareMatrix &A, labelList &pivotIndices) const
{
    autoPtr<scalarList> dvidt;
    dvidt.reset((this->dvidT(p, T, sol, A, pivotIndices)).ptr());
    autoPtr<scalarList> dXdT_G_(this->dXdT_gas(p, T, sol, dvidt));
    autoPtr<scalarList> dXdT_L_(this->dXdT_liq(p, T, sol, dvidt));
    scalar dvfdT_ = this->dvfdT(p, T, sol, dvidt());
    scalar vf_ = sol.vaporfra;
    scalar Z_G_ = this->ThermoMixture::Z(p, T, sol.X_gas(), 1);
    scalar Z_L_ = this->ThermoMixture::Z(p, T, sol.X_liq(), 0);
    scalar dZdT_G_ = this->ThermoMixture::dZdT(p, T, sol.X_gas(), 1);
    scalar dZdT_L_ = this->ThermoMixture::dZdT(p, T, sol.X_liq(), 0);
    for (int i = 0; i < this->X_.size(); i++)
    {
        dZdT_G_ +=
            this->ThermoMixture::dZdxi(p, T, i, sol.X_gas(), 1) * dXdT_G_()[i];
        dZdT_L_ +=
            this->ThermoMixture::dZdxi(p, T, i, sol.X_liq(), 0) * dXdT_L_()[i];
    }
    return dvfdT_ * (Z_G_ - Z_L_) + vf_ * dZdT_G_ + (1 - vf_) * dZdT_L_;
}
template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::dZdP(scalar p, scalar T, solution &sol) const
{
    autoPtr<scalarList> dvidp;
    dvidp.reset((this->dvidP(p, T, sol)).ptr());
    autoPtr<scalarList> dXdP_G_(this->dXdP_gas(p, T, sol, dvidp));
    autoPtr<scalarList> dXdP_L_(this->dXdP_liq(p, T, sol, dvidp));
    scalar dvfdP_ = this->dvfdP(p, T, sol, dvidp());
    scalar vf_ = sol.vaporfra;
    scalar Z_G_ = this->ThermoMixture::Z(p, T, sol.X_gas(), 1);
    scalar Z_L_ = this->ThermoMixture::Z(p, T, sol.X_liq(), 0);
    scalar dZdP_G_ = this->ThermoMixture::dZdP(p, T, sol.X_gas(), 1);
    scalar dZdP_L_ = this->ThermoMixture::dZdP(p, T, sol.X_liq(), 0);
    for (int i = 0; i < this->X_.size(); i++)
    {
        dZdP_G_ +=
            this->ThermoMixture::dZdxi(p, T, i, sol.X_gas(), 1) * dXdP_G_()[i];
        dZdP_L_ +=
            this->ThermoMixture::dZdxi(p, T, i, sol.X_liq(), 0) * dXdP_L_()[i];
    }
    // scalar temp = (vf * ZGas + (1 - vf) * ZLiq);
    // return (dvfdt * ZGas * ZGas + vf * 2 * ZGas * dZgasdT + 2 * (1 - vf) * ZLiq
    // * dZliqdT - dvfdt * ZLiq * ZLiq) / temp - (vf * ZGas * ZGas + (1 - vf) *
    // ZLiq * ZLiq) / (temp * temp) * (dvfdt * ZGas + vf * dZgasdT + (1 - vf) *
    // dZliqdT - dvfdt * ZLiq);
    //Info << "dZdP_G_=" << dZdP_G_ << endl;
    //Info << "dZdP_L_=" << dZdP_L_ << endl;
    return dvfdP_ * (Z_G_ - Z_L_) + vf_ * dZdP_G_ + (1 - vf_) * dZdP_L_;
}

template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::dZdP(scalar p, scalar T, solution &sol, scalarSquareMatrix &A, labelList &pivotIndices) const
{
    autoPtr<scalarList> dvidp;
    dvidp.reset((this->dvidP(p, T, sol, A, pivotIndices)).ptr());
    autoPtr<scalarList> dXdP_G_(this->dXdP_gas(p, T, sol, dvidp));
    autoPtr<scalarList> dXdP_L_(this->dXdP_liq(p, T, sol, dvidp));
    scalar dvfdP_ = this->dvfdP(p, T, sol, dvidp());
    scalar vf_ = sol.vaporfra;
    scalar Z_G_ = this->ThermoMixture::Z(p, T, sol.X_gas(), 1);
    scalar Z_L_ = this->ThermoMixture::Z(p, T, sol.X_liq(), 0);
    scalar dZdP_G_ = this->ThermoMixture::dZdP(p, T, sol.X_gas(), 1);
    scalar dZdP_L_ = this->ThermoMixture::dZdP(p, T, sol.X_liq(), 0);
    for (int i = 0; i < this->X_.size(); i++)
    {
        dZdP_G_ +=
            this->ThermoMixture::dZdxi(p, T, i, sol.X_gas(), 1) * dXdP_G_()[i];
        dZdP_L_ +=
            this->ThermoMixture::dZdxi(p, T, i, sol.X_liq(), 0) * dXdP_L_()[i];
    }

    return dvfdP_ * (Z_G_ - Z_L_) + vf_ * dZdP_G_ + (1 - vf_) * dZdP_L_;
}
template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::dZdXi(scalar p,
                                scalar T,
                                label di,
                                solution &sol) const
{
    autoPtr<scalarList> dXdXi_G_(this->dXdXi_gas(p, T, di, sol));
    autoPtr<scalarList> dXdXi_L_(this->dXdXi_liq(p, T, di, sol));
    scalar dvfdXi_ = this->dvfdXi(p, T, di, sol);
    scalar vf_ = sol.vaporfra;
    scalar Z_G_ = this->ThermoMixture::Z(p, T, sol.X_gas(), 1);
    scalar Z_L_ = this->ThermoMixture::Z(p, T, sol.X_liq(), 0);
    scalar dZdXi_G_ = 0;
    scalar dZdXi_L_ = 0;
    for (int i = 0; i < this->X_.size(); i++)
    {
        dZdXi_G_ +=
            this->ThermoMixture::dZdxi(p, T, i, sol.X_gas(), 1) * dXdXi_G_()[i];
        dZdXi_L_ +=
            this->ThermoMixture::dZdxi(p, T, i, sol.X_liq(), 0) * dXdXi_L_()[i];
    }
    // scalar temp = (vf * ZGas + (1 - vf) * ZLiq);
    // return (dvfdt * ZGas * ZGas + vf * 2 * ZGas * dZgasdT + 2 * (1 - vf) * ZLiq
    // * dZliqdT - dvfdt * ZLiq * ZLiq) / temp - (vf * ZGas * ZGas + (1 - vf) *
    // ZLiq * ZLiq) / (temp * temp) * (dvfdt * ZGas + vf * dZgasdT + (1 - vf) *
    // dZliqdT - dvfdt * ZLiq);
    return dvfdXi_ * (Z_G_ - Z_L_) + vf_ * dZdXi_G_ + (1 - vf_) * dZdXi_L_;
}

// rho
template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::drhodT(scalar p, scalar T) const
{
    autoPtr<solution> sol(TPN(p, T));
    scalar Z_M_ = Z(p, T, sol());
    scalar dZdT_M_ = dZdT(p, T, sol());
    scalar drhodT_M_ = -p * this->W() / (Z_M_ * Z_M_ * RR * 1.0e-03 * T * T) *
                       (dZdT_M_ * T + Z_M_);
    return drhodT_M_;
}

template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::drrhodT(scalar p, scalar T) const
{
    autoPtr<solution> sol(TPN(p, T));
    scalar Z_M_ = Z(p, T, sol());
    scalar dZdT_M_ = dZdT(p, T, sol());
    scalar drhodT_M_ = -p * this->W() / (Z_M_ * Z_M_ * RR * 1.0e-03 * T * T) *
                       (dZdT_M_ * T + Z_M_);
    return -1 / sqr(rho(p, T)) * drhodT(p, T);
}

template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::drhodT(scalar p, scalar T, solution &sol) const
{
    scalar Z_M_ = Z(p, T, sol);
    scalar dZdT_M_ = dZdT(p, T, sol);
    scalar drhodT_M_ = -p * this->W() / (Z_M_ * Z_M_ * RR * 1.0e-03 * T * T) *
                       (dZdT_M_ * T + Z_M_);
    return drhodT_M_;
}

template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::drhodT_noVLE(scalar p, scalar T, solution &sol) const
{
    scalar Z_M_ = this->ThermoMixture::Z(p, T, sol.X_gas(), 1);
    scalar dZdT_M_ = this->ThermoMixture::dZdT(p, T, sol.X_gas(), 1);
    scalar drhodT_M_ = -p * this->W() / (Z_M_ * Z_M_ * RR * 1.0e-03 * T * T) *
                       (dZdT_M_ * T + Z_M_);
    return drhodT_M_;
}
template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::drhodP(scalar p, scalar T, solution &sol) const
{
    scalar Z_M_ = Z(p, T, sol);
    scalar dZdP_M_ = dZdP(p, T, sol);

    scalar drhodP_M_ = this->W() / (Z_M_ * RR * 1.0e-03 * T) -
                       p * this->W() / (Z_M_ * Z_M_ * RR * 1.0e-03 * T) * dZdP_M_;
    return drhodP_M_;
}

template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::drhodP_noVLE(scalar p, scalar T, solution &sol) const
{
    scalar Z_M_ = this->ThermoMixture::Z(p, T, sol.X_gas(), 1);
    scalar dZdP_M_ = this->ThermoMixture::dZdP(p, T, sol.X_gas(), 1);

    scalar drhodP_M_ = this->W() / (Z_M_ * RR * 1.0e-03 * T) -
                       p * this->W() / (Z_M_ * Z_M_ * RR * 1.0e-03 * T) * dZdP_M_;
    return drhodP_M_;
}
template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::drhodXi(scalar p,
                                  scalar T,
                                  label di,
                                  solution &sol) const
{
    scalar Z_M_ = Z(p, T, sol);
    scalar dZdXi_M_ = dZdXi(p, T, di, sol);
    scalar drhodXi_M_ =
        p * this->dWdXi(di) / (Z_M_ * RR * 1.0e-03 * T) -
        p * this->W() / (Z_M_ * Z_M_ * RR * 1.0e-03 * T) * dZdXi_M_;
    return drhodXi_M_;
}

// Es
template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::dEdT(scalar p, scalar T, solution &sol) const
{
    scalar rho_M_ = rho(p, T, sol);
    return Cp_Hs(p, T, sol) + p / (rho_M_ * rho_M_) * drhodT(p, T, sol);
}

template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::dEdP(scalar p, scalar T, solution &sol) const
{
    scalar rho_M_ = rho(p, T, sol);
    return dHsdP(p, T, sol) - 1 / rho_M_ +
           p / (rho_M_ * rho_M_) * drhodP(p, T, sol);
}

template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::dEdXi(scalar p, scalar T, label di, solution &sol) const
{
    scalar rho_M_ = rho(p, T, sol);
    return dHsdXi(p, T, di, sol) + p / (rho_M_ * rho_M_) * drhodXi(p, T, di, sol);
}

// *** system (E,rho,X) *** //

template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::dTdE_rhoX(scalar p, scalar T, solution &sol) const
{
    scalar dEdT_ = dEdT(p, T, sol);
    scalar dEdP_ = dEdP(p, T, sol);
    scalar drhodT_ = drhodT(p, T, sol);
    scalar drhodP_ = drhodP(p, T, sol);
    scalar det = dEdT_ * drhodP_ - dEdP_ * drhodT_;
    return drhodP_ / det;
}

template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::dTdrho_EX(scalar p, scalar T, solution &sol) const
{
    scalar dEdT_ = dEdT(p, T, sol);
    scalar dEdP_ = dEdP(p, T, sol);
    scalar drhodT_ = drhodT(p, T, sol);
    scalar drhodP_ = drhodP(p, T, sol);
    scalar det = dEdT_ * drhodP_ - dEdP_ * drhodT_;
    return -dEdP_ / det;
}

template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::dPdE_rhoX(scalar p, scalar T, solution &sol) const
{
    scalar dEdT_ = dEdT(p, T, sol);
    scalar dEdP_ = dEdP(p, T, sol);
    scalar drhodT_ = drhodT(p, T, sol);
    scalar drhodP_ = drhodP(p, T, sol);
    scalar det = dEdT_ * drhodP_ - dEdP_ * drhodT_;
    return -drhodT_ / det;
}

template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::dPdrho_EX(scalar p, scalar T, solution &sol) const
{
    scalar dEdT_ = dEdT(p, T, sol);
    scalar dEdP_ = dEdP(p, T, sol);
    scalar drhodT_ = drhodT(p, T, sol);
    scalar drhodP_ = drhodP(p, T, sol);
    scalar det = dEdT_ * drhodP_ - dEdP_ * drhodT_;
    return dEdT_ / det;
}

template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::dTdXi_Erho(scalar p, scalar T, label di, solution &sol) const
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
template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::dPdXi_Erho(scalar p, scalar T, label di, solution &sol) const
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
template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::dvfdE_rhoX(scalar p, scalar T, solution &sol) const
{
    return dvfdT(p, T, sol) * dTdE_rhoX(p, T, sol) + dvfdP(p, T, sol) * dPdE_rhoX(p, T, sol);
}
template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::dvfdrho_EX(scalar p, scalar T, solution &sol) const
{
    return dvfdT(p, T, sol) * dTdrho_EX(p, T, sol) + dvfdP(p, T, sol) * dPdrho_EX(p, T, sol);
}
template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::dvfdXi_Erho(scalar p, scalar T, label di, solution &sol) const
{
    return dvfdXi(p, T, di, sol) + dvfdT(p, T, sol) * dTdXi_Erho(p, T, di, sol) + dvfdP(p, T, sol) * dPdXi_Erho(p, T, di, sol);
}
// * * * * * * * * * * * * * * * *  Derivative  properties * * * * * * * * * * *
// * * * //

// isothermal compressibility
template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::kappaT(scalar p, scalar T, solution &sol) const
{
    scalar rho_ = rho(p, T, sol);
    scalar drhodP_ = drhodP(p, T, sol);
    return drhodP_ / rho_;
}

template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::kappaT_noVLE(scalar p, scalar T, solution &sol) const
{
    scalar rho_ = rho_noVLE(p, T, sol);
    scalar drhodP_ = drhodP_noVLE(p, T, sol);
    return drhodP_ / rho_;
}

// isobaric expansivity
template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::alphaP(scalar p, scalar T, solution &sol) const
{
    scalar rho_ = rho(p, T, sol);
    scalar drhodT_ = drhodT(p, T, sol);
    return -drhodT_ / rho_;
}
template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::alphaP_noVLE(scalar p, scalar T, solution &sol) const
{
    scalar rho_ = rho_noVLE(p, T, sol);
    scalar drhodT_ = drhodT_noVLE(p, T, sol);
    return -drhodT_ / rho_;
}
// isentropic compressibility
template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::kappaS(scalar p, scalar T, solution &sol) const
{
    scalar rho_ = rho(p, T, sol);
    scalar kappaT_ = kappaT(p, T, sol);
    scalar Cp_ = Cp(p, T, sol);
    scalar alphaP_ = alphaP(p, T, sol);
    return kappaT_ - T * alphaP_ * alphaP_ / (rho_ * Cp_);
}

template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::kappaS_noVLE(scalar p, scalar T, solution &sol) const
{
    scalar rho_ = rho_noVLE(p, T, sol);
    scalar kappaT_ = kappaT_noVLE(p, T, sol);
    scalar Cp_ = this->dHadT(p, T, sol.X_gas(), 1);
    scalar alphaP_ = alphaP_noVLE(p, T, sol);
    return kappaT_ - T * alphaP_ * alphaP_ / (rho_ * Cp_);
}
// speed of sound
template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::c(scalar p, scalar T, solution &sol) const
{
    scalar rho_ = rho(p, T, sol);
    scalar kappaS_ = kappaS(p, T, sol);
    if (kappaS_ < 0)
    {
        //return c_noVLE(p, T, sol);
        FatalErrorInFunction << " kappaS_ =" << kappaS_
                             << ",p=" << p << ",T=" << T << ",X=" << this->X_ << exit(FatalError);
    }
    return sqrt(1.0 / (kappaS_ * rho_));
}

template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::c_opt(scalar p, scalar T, scalar rho_, solution &sol) const
{
    scalar W_ = this->W();
    scalar Z_M_ = Z(p, T, sol);
    autoPtr<scalarSquareMatrix> pA(A_dvid_(p, T, sol));
    labelList pivotIndices(pA().m());
    LUDecompose(pA(), pivotIndices);
    scalar dZdP_M_ = dZdP(p, T, sol, pA(), pivotIndices);
    scalar drhodP_ = W_ / (Z_M_ * RR * 1.0e-03 * T) - p * W_ / (Z_M_ * Z_M_ * RR * 1.0e-03 * T) * dZdP_M_;
    scalar kappaT_ = drhodP_ / rho_;
    scalar Cp_ = Cp(p, T, sol, pA(), pivotIndices);
    scalar dZdT_M_ = dZdT(p, T, sol, pA(), pivotIndices);
    scalar drhodT_ = -p * W_ / (Z_M_ * Z_M_ * RR * 1.0e-03 * T * T) * (dZdT_M_ * T + Z_M_);
    scalar alphaP_ = -drhodT_ / rho_;
    scalar kappaS_ = kappaT_ - T * alphaP_ * alphaP_ / (rho_ * Cp_);
    if (kappaS_ < 0)
    {
        //return c_noVLE(p, T, sol);
        FatalErrorInFunction << " kappaS_ =" << kappaS_
                             << ",p=" << p << ",T=" << T << ",X=" << this->X_ << exit(FatalError);
    }
    return sqrt(1.0 / (kappaS_ * rho_));
}

template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::c_opt_full(scalar p, scalar T, scalar rho_, solution &sol) const
{
    const scalar sqrt2 = sqrt(2.0);
    const scalar cc0 = 2 * sqrt(2.0);
    const scalar cc1 = 1.0 + sqrt(2.0);
    const scalar cc2 = 1.0 - sqrt(2.0);

    const scalar small8 = 1.e-08; //PY: need to be noticed

    scalar W_ = this->W();
    scalar vf_ = sol.vaporfra;

    scalar A_G_ = this->ThermoMixture::A(p, T, sol.X_gas());
    scalar A_L_ = this->ThermoMixture::A(p, T, sol.X_liq());

    scalar B_G_ = this->ThermoMixture::B(p, T, sol.X_gas());
    scalar B_L_ = this->ThermoMixture::B(p, T, sol.X_liq());

    scalar B_G_2 = B_G_ * B_G_;
    scalar B_L_2 = B_L_ * B_L_;

    scalar Z_G_ = this->ThermoMixture::Z(p, T, A_G_, B_G_, sol.X_gas(), 1);
    scalar Z_L_ = this->ThermoMixture::Z(p, T, A_L_, B_L_, sol.X_liq(), 0);
    scalar Z_M_ = vf_ * Z_G_ + (1 - vf_) * Z_L_;

    scalar Z_G_2 = Z_G_ * Z_G_;
    scalar Z_L_2 = Z_L_ * Z_L_;

    scalar a2_G_ = B_G_ - 1.0;
    scalar a2_L_ = B_L_ - 1.0;

    scalar a1_G_ = A_G_ - 2.0 * B_G_ - 3.0 * B_G_ * B_G_;
    scalar a1_L_ = A_L_ - 2.0 * B_L_ - 3.0 * B_L_ * B_L_;

    scalar dAdP_G_ = this->dAdP(p, T, sol.X_gas()), dBdP_G_ = this->dBdP(p, T, sol.X_gas());
    scalar dAdP_L_ = this->dAdP(p, T, sol.X_liq()), dBdP_L_ = this->dBdP(p, T, sol.X_liq());

    scalar da2dP_G_ = dBdP_G_;
    scalar da2dP_L_ = dBdP_L_;

    scalar da1dP_G_ = dAdP_G_ - 2.0 * dBdP_G_ - 6.0 * B_G_ * dBdP_G_;
    scalar da1dP_L_ = dAdP_L_ - 2.0 * dBdP_L_ - 6.0 * B_L_ * dBdP_L_;

    scalar da0dP_G_ = -(dAdP_G_ * B_G_ + A_G_ * dBdP_G_) + 2 * B_G_ * dBdP_G_ + 3 * B_G_2 * dBdP_G_;
    scalar da0dP_L_ = -(dAdP_L_ * B_L_ + A_L_ * dBdP_L_) + 2 * B_L_ * dBdP_L_ + 3 * B_L_2 * dBdP_L_;

    scalar dZdP_G_0 = -(da2dP_G_ * Z_G_2 + da1dP_G_ * Z_G_ + da0dP_G_) / (3 * Z_G_2 + 2 * Z_G_ * a2_G_ + a1_G_);
    scalar dZdP_L_0 = -(da2dP_L_ * Z_L_2 + da1dP_L_ * Z_L_ + da0dP_L_) / (3 * Z_L_2 + 2 * Z_L_ * a2_L_ + a1_L_);

    scalar dAdT_G_ = this->ThermoMixture::dAdT(p, T, sol.X_gas()), dBdT_G_ = this->ThermoMixture::dBdT(p, T, sol.X_gas());
    scalar dAdT_L_ = this->ThermoMixture::dAdT(p, T, sol.X_liq()), dBdT_L_ = this->ThermoMixture::dBdT(p, T, sol.X_liq());

    scalar da2dT_G_ = dBdT_G_;
    scalar da2dT_L_ = dBdT_L_;

    scalar da1dT_G_ = dAdT_G_ - 2.0 * dBdT_G_ - 6.0 * B_G_ * dBdT_G_;
    scalar da1dT_L_ = dAdT_L_ - 2.0 * dBdT_L_ - 6.0 * B_L_ * dBdT_L_;

    scalar da0dT_G_ = -(dAdT_G_ * B_G_ + A_G_ * dBdT_G_) + 2 * B_G_ * dBdT_G_ + 3 * B_G_2 * dBdT_G_;
    scalar da0dT_L_ = -(dAdT_L_ * B_L_ + A_L_ * dBdT_L_) + 2 * B_L_ * dBdT_L_ + 3 * B_L_2 * dBdT_L_;

    scalar dZdT_G_0 = -(da2dT_G_ * Z_G_2 + da1dT_G_ * Z_G_ + da0dT_G_) / (3 * Z_G_2 + 2 * Z_G_ * a2_G_ + a1_G_);
    scalar dZdT_L_0 = -(da2dT_L_ * Z_L_2 + da1dT_L_ * Z_L_ + da0dT_L_) / (3 * Z_L_2 + 2 * Z_L_ * a2_L_ + a1_L_);

    scalarSquareMatrix A(this->X_.size(), Zero);

    List<autoPtr<scalarList>> ddxi_Ln_fugcoef_G_(this->X_.size()),
        ddxi_Ln_fugcoef_L_(this->X_.size());

    scalarList Bi(this->N_);
    scalarList dBidP(this->N_);
    scalarList dBidT(this->N_);
    for (int i = 0; i < this->X_.size(); i++)
    {
        Bi[i] = (*this)[i].B(p, T);
        dBidP[i] = (*this)[i].dBdP(p, T);
        dBidT[i] = (*this)[i].dBdT(p, T);
    }

    scalarList dzdxi_G_(this->N_);
    scalarList dzdxi_L_(this->N_);
    scalarList dAdxi_G_(this->N_);
    scalarList dAdxi_L_(this->N_);
    scalarList dBdxi_G_(this->N_);
    scalarList dBdxi_L_(this->N_);
    for (int i = 0; i < this->X_.size(); i++)
    {
        dAdxi_G_[i] = this->dAdxi(p, T, i, sol.X_gas());
        dBdxi_G_[i] = this->dBdxi(p, T, i, sol.X_gas());
        dAdxi_L_[i] = this->dAdxi(p, T, i, sol.X_liq());
        dBdxi_L_[i] = this->dBdxi(p, T, i, sol.X_liq());

        scalar da2dxi_G_ = dBdxi_G_[i];
        scalar da2dxi_L_ = dBdxi_L_[i];

        scalar da1dxi_G_ = dAdxi_G_[i] - 2.0 * dBdxi_G_[i] - 6.0 * B_G_ * dBdxi_G_[i];
        scalar da1dxi_L_ = dAdxi_L_[i] - 2.0 * dBdxi_L_[i] - 6.0 * B_L_ * dBdxi_L_[i];

        scalar da0dxi_G_ = -(dAdxi_G_[i] * B_G_ + A_G_ * dBdxi_G_[i]) + 2 * B_G_ * dBdxi_G_[i] + 3 * B_G_ * B_G_ * dBdxi_G_[i];
        scalar da0dxi_L_ = -(dAdxi_L_[i] * B_L_ + A_L_ * dBdxi_L_[i]) + 2 * B_L_ * dBdxi_L_[i] + 3 * B_L_ * B_L_ * dBdxi_L_[i];

        dzdxi_G_[i] = -(da2dxi_G_ * Z_G_2 + da1dxi_G_ * Z_G_ + da0dxi_G_) / (3 * Z_G_2 + 2 * Z_G_ * a2_G_ + a1_G_);
        dzdxi_L_[i] = -(da2dxi_L_ * Z_L_2 + da1dxi_L_ * Z_L_ + da0dxi_L_) / (3 * Z_L_2 + 2 * Z_L_ * a2_L_ + a1_L_);
    }

    scalar xxx_G_ = ((Z_G_ + cc1 * B_G_) / (Z_G_ + cc2 * B_G_));
    scalar xxx_L_ = ((Z_L_ + cc1 * B_L_) / (Z_L_ + cc2 * B_L_));
    scalar logxxx_G_ = log(xxx_G_);
    scalar logxxx_L_ = log(xxx_L_);

    autoPtr<scalarList> pA2_G_(this->A2(p, T, sol.X_gas()));
    const scalarList &A2_G_ = pA2_G_();
    autoPtr<scalarList> pA2_L_(this->A2(p, T, sol.X_liq()));
    const scalarList &A2_L_ = pA2_L_();

    for (int i = 0; i < this->X_.size(); i++)
    {
        scalarList &ddxi_Ln_fugcoef = *(new scalarList(this->N_));
        ddxi_Ln_fugcoef_G_[i].reset(&ddxi_Ln_fugcoef);
        scalar dAmixzdxi = dAdxi_G_[i];
        scalar dBmixzdxi = dBdxi_G_[i];

        autoPtr<scalarList> pdAmix2dxi(this->dA2dxi(p, T, i, sol.X_gas()));
        const scalarList &dAmix2dxi = pdAmix2dxi();

        scalar zfactor00 = Z_G_;
        scalar ZMB = zfactor00 - B_G_;

        if (ZMB <= 0.0)
        {
            FatalErrorInFunction << "ZMB <= 0.0" << exit(FatalError);
            zfactor00 = B_G_ + small8;
            ZMB = small8;
        }
        forAll(this->X_, spid)
        {
            ddxi_Ln_fugcoef[spid] = -Bi[spid] / (B_G_2)*dBmixzdxi * (Z_G_ - 1) + Bi[spid] / B_G_ * dzdxi_G_[i];
            ddxi_Ln_fugcoef[spid] -= (dzdxi_G_[i] - dBmixzdxi) / (Z_G_ - B_G_);
            ddxi_Ln_fugcoef[spid] -= (dAmix2dxi[spid] / (sqrt2 * B_G_) - A2_G_[spid] * dBmixzdxi / (sqrt2 * B_G_2) - (dAmixzdxi * Bi[spid]) / (cc0 * B_G_2) + (A_G_ * Bi[spid] * dBmixzdxi) / (sqrt2 * B_G_2 * B_G_)) * logxxx_G_;
            ddxi_Ln_fugcoef[spid] -= (A2_G_[spid] / (sqrt2 * B_G_) - A_G_ * Bi[spid] / (cc0 * B_G_2)) * ((dzdxi_G_[i] + cc1 * dBmixzdxi) / (Z_G_ + cc1 * B_G_) - (dzdxi_G_[i] + cc2 * dBmixzdxi) / (Z_G_ + cc2 * B_G_));
        }
    }

    for (int i = 0; i < this->X_.size(); i++)
    {
        scalarList &ddxi_Ln_fugcoef = *(new scalarList(this->N_));
        ddxi_Ln_fugcoef_L_[i].reset(&ddxi_Ln_fugcoef);
        scalar dAmixzdxi = dAdxi_L_[i];
        scalar dBmixzdxi = dBdxi_L_[i];

        autoPtr<scalarList> pdAmix2dxi(this->dA2dxi(p, T, i, sol.X_liq()));
        const scalarList &dAmix2dxi = pdAmix2dxi();

        scalar zfactor00 = Z_L_;
        scalar ZMB = zfactor00 - B_L_;

        if (ZMB <= 0.0)
        {
            FatalErrorInFunction << "ZMB <= 0.0" << exit(FatalError);
            zfactor00 = B_L_ + small8;
            ZMB = small8;
        }
        forAll(this->X_, spid)
        {
            ddxi_Ln_fugcoef[spid] = -Bi[spid] / (B_L_2)*dBmixzdxi * (Z_L_ - 1) + Bi[spid] / B_L_ * dzdxi_L_[i];
            ddxi_Ln_fugcoef[spid] -= (dzdxi_L_[i] - dBmixzdxi) / (Z_L_ - B_L_);
            ddxi_Ln_fugcoef[spid] -= (dAmix2dxi[spid] / (sqrt2 * B_L_) - A2_L_[spid] * dBmixzdxi / (sqrt2 * B_L_2) - (dAmixzdxi * Bi[spid]) / (cc0 * B_L_2) + (A_L_ * Bi[spid] * dBmixzdxi) / (sqrt2 * B_L_2 * B_L_)) * logxxx_L_;
            ddxi_Ln_fugcoef[spid] -= (A2_L_[spid] / (sqrt2 * B_L_) - A_L_ * Bi[spid] / (cc0 * B_L_2)) * ((dzdxi_L_[i] + cc1 * dBmixzdxi) / (Z_L_ + cc1 * B_L_) - (dzdxi_L_[i] + cc2 * dBmixzdxi) / (Z_L_ + cc2 * B_L_));
        }
    }
    for (int i = 0; i < this->X_.size(); i++)
    {
        for (int j = 0; j < this->X_.size(); j++)
        {
            A(i, j) = 1;
            if (i == j)
                A(i, j) -= this->X_[i] / (sol.X_gas()[i] * sol.X_liq()[i]);
            for (int k = 0; k < this->X_.size(); k++)
            {
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

    labelList pivotIndices(A.m());
    LUDecompose(A, pivotIndices);

    autoPtr<scalarList> dvidp(new scalarList(this->X_.size()));

    int n_sp = 0;
    for (int i = 0; i < this->X_.size(); i++)
    {
        if (this->X_[i] > 0)
            n_sp++;
    }
    if (n_sp == 1)
    {
        for (int i = 0; i < this->X_.size(); i++)
            dvidp()[i] = 0;
    }
    else
    {
        autoPtr<scalarList> ddP_Ln_fugcoef_G_(new scalarList(this->N_)), ddP_Ln_fugcoef_L_(new scalarList(this->N_));
        {

            autoPtr<scalarList> pdAmix2dP(this->dA2dP(p, T, sol.X_gas()));
            const scalarList &dAmix2dP = pdAmix2dP();

            scalar zfactor00 = Z_G_;
            scalar ZMB = zfactor00 - B_G_;

            if (ZMB <= 0.0)
            {
                FatalErrorInFunction << "ZMB <= 0.0" << exit(FatalError);
                zfactor00 = B_G_ + small8;
                ZMB = small8;
            }
            forAll(this->X_, spid)
            {
                ddP_Ln_fugcoef_G_()[spid] = dBidP[spid] / B_G_ * (Z_G_ - 1) - Bi[spid] / (B_G_2)*dBdP_G_ * (Z_G_ - 1) + Bi[spid] / B_G_ * dZdP_G_0;
                ddP_Ln_fugcoef_G_()[spid] -= (dZdP_G_0 - dBdP_G_) / (Z_G_ - B_G_);
                ddP_Ln_fugcoef_G_()[spid] -= (dAmix2dP[spid] / (sqrt2 * B_G_) - A2_G_[spid] * dBdP_G_ / (sqrt2 * B_G_2) - (dAdP_G_ * Bi[spid] + A_G_ * dBidP[spid]) / (cc0 * B_G_2) + (A_G_ * Bi[spid] * dBdP_G_) / (sqrt2 * B_G_2 * B_G_)) * logxxx_G_;
                ddP_Ln_fugcoef_G_()[spid] -= (A2_G_[spid] / (sqrt2 * B_G_) - A_G_ * Bi[spid] / (cc0 * B_G_2)) * ((dZdP_G_0 + cc1 * dBdP_G_) / (Z_G_ + cc1 * B_G_) - (dZdP_G_0 + cc2 * dBdP_G_) / (Z_G_ + cc2 * B_G_));
            }
        }

        {

            autoPtr<scalarList> pdAmix2dP(this->dA2dP(p, T, sol.X_liq()));
            const scalarList &dAmix2dP = pdAmix2dP();

            scalar zfactor00 = Z_L_;
            scalar ZMB = zfactor00 - B_L_;

            if (ZMB <= 0.0)
            {
                FatalErrorInFunction << "ZMB <= 0.0" << exit(FatalError);
                zfactor00 = B_L_ + small8;
                ZMB = small8;
            }
            forAll(this->X_, spid)
            {
                ddP_Ln_fugcoef_L_()[spid] = dBidP[spid] / B_L_ * (Z_L_ - 1) - Bi[spid] / (B_L_2)*dBdP_L_ * (Z_L_ - 1) + Bi[spid] / B_L_ * dZdP_L_0;
                ddP_Ln_fugcoef_L_()[spid] -= (dZdP_L_0 - dBdP_L_) / (Z_L_ - B_L_);
                ddP_Ln_fugcoef_L_()[spid] -= (dAmix2dP[spid] / (sqrt2 * B_L_) - A2_L_[spid] * dBdP_L_ / (sqrt2 * B_L_2) - (dAdP_L_ * Bi[spid] + A_L_ * dBidP[spid]) / (cc0 * B_L_2) + (A_L_ * Bi[spid] * dBdP_L_) / (sqrt2 * B_L_2 * B_L_)) * logxxx_L_;
                ddP_Ln_fugcoef_L_()[spid] -= (A2_L_[spid] / (sqrt2 * B_L_) - A_L_ * Bi[spid] / (cc0 * B_L_2)) * ((dZdP_L_0 + cc1 * dBdP_L_) / (Z_L_ + cc1 * B_L_) - (dZdP_L_0 + cc2 * dBdP_L_) / (Z_L_ + cc2 * B_L_));
            }
        }

        for (int i = 0; i < this->X_.size(); i++)
            dvidp()[i] = (1 - vf_) * vf_ * (ddP_Ln_fugcoef_G_()[i] - ddP_Ln_fugcoef_L_()[i]);

        LUBacksubstitute(A, pivotIndices, dvidp());
    }

    autoPtr<scalarList> dXdP_G_(this->dXdP_gas(p, T, sol, dvidp));
    autoPtr<scalarList> dXdP_L_(this->dXdP_liq(p, T, sol, dvidp));
    scalar dvfdP_ = this->dvfdP(p, T, sol, dvidp());

    scalar dZdP_G_ = dZdP_G_0; //this->ThermoMixture::dZdP(p, T, sol.X_gas(), 1);
    scalar dZdP_L_ = dZdP_L_0; //this->ThermoMixture::dZdP(p, T, sol.X_liq(), 0);

    for (int i = 0; i < this->X_.size(); i++)
    {
        dZdP_G_ += dzdxi_G_[i] * dXdP_G_()[i];
        dZdP_L_ += dzdxi_L_[i] * dXdP_L_()[i];
    }

    scalar dZdP_M_ = dvfdP_ * (Z_G_ - Z_L_) + vf_ * dZdP_G_ + (1 - vf_) * dZdP_L_;

    scalar drhodP_ = W_ / (Z_M_ * RR * 1.0e-03 * T) - p * W_ / (Z_M_ * Z_M_ * RR * 1.0e-03 * T) * dZdP_M_;
    scalar kappaT_ = drhodP_ / rho_;

    autoPtr<scalarList> dvidt(new scalarList(this->X_.size()));

    if (n_sp == 1 || vf_ == 0 || vf_ == 1)
    {
        for (int i = 0; i < this->X_.size(); i++)
            dvidt()[i] = 0;
    }
    else
    {
        autoPtr<scalarList> ddT_Ln_fugcoef_G_(new scalarList(this->N_)), ddT_Ln_fugcoef_L_(new scalarList(this->N_));
        {

            autoPtr<scalarList> pdAmix2dT(this->dA2dT(p, T, sol.X_gas()));
            const scalarList &dAmix2dT = pdAmix2dT();

            scalar zfactor00 = Z_G_;
            scalar ZMB = zfactor00 - B_G_;

            if (ZMB <= 0.0)
            {
                FatalErrorInFunction << "ZMB <= 0.0" << exit(FatalError);
                zfactor00 = B_G_ + small8;
                ZMB = small8;
            }
            forAll(this->X_, spid)
            {
                ddT_Ln_fugcoef_G_()[spid] = dBidT[spid] / B_G_ * (Z_G_ - 1) - Bi[spid] / (B_G_2)*dBdT_G_ * (Z_G_ - 1) + Bi[spid] / B_G_ * dZdT_G_0;
                ddT_Ln_fugcoef_G_()[spid] -= (dZdT_G_0 - dBdT_G_) / (Z_G_ - B_G_);
                ddT_Ln_fugcoef_G_()[spid] -= (dAmix2dT[spid] / (sqrt2 * B_G_) - A2_G_[spid] * dBdT_G_ / (sqrt2 * B_G_2) - (dAdT_G_ * Bi[spid] + A_G_ * dBidT[spid]) / (cc0 * B_G_2) + (A_G_ * Bi[spid] * dBdT_G_) / (sqrt2 * B_G_2 * B_G_)) * logxxx_G_;
                ddT_Ln_fugcoef_G_()[spid] -= (A2_G_[spid] / (sqrt2 * B_G_) - A_G_ * Bi[spid] / (cc0 * B_G_2)) * ((dZdT_G_0 + cc1 * dBdT_G_) / (Z_G_ + cc1 * B_G_) - (dZdT_G_0 + cc2 * dBdT_G_) / (Z_G_ + cc2 * B_G_));
            }
        }

        {

            autoPtr<scalarList> pdAmix2dT(this->dA2dT(p, T, sol.X_liq()));
            const scalarList &dAmix2dT = pdAmix2dT();

            scalar zfactor00 = Z_L_;
            scalar ZMB = zfactor00 - B_L_;

            if (ZMB <= 0.0)
            {
                FatalErrorInFunction << "ZMB <= 0.0" << exit(FatalError);
                zfactor00 = B_L_ + small8;
                ZMB = small8;
            }
            forAll(this->X_, spid)
            {
                ddT_Ln_fugcoef_L_()[spid] = dBidT[spid] / B_L_ * (Z_L_ - 1) - Bi[spid] / (B_L_2)*dBdT_L_ * (Z_L_ - 1) + Bi[spid] / B_L_ * dZdT_L_0;
                ddT_Ln_fugcoef_L_()[spid] -= (dZdT_L_0 - dBdT_L_) / (Z_L_ - B_L_);
                ddT_Ln_fugcoef_L_()[spid] -= (dAmix2dT[spid] / (sqrt2 * B_L_) - A2_L_[spid] * dBdT_L_ / (sqrt2 * B_L_2) - (dAdT_L_ * Bi[spid] + A_L_ * dBidT[spid]) / (cc0 * B_L_2) + (A_L_ * Bi[spid] * dBdT_L_) / (sqrt2 * B_L_2 * B_L_)) * logxxx_L_;
                ddT_Ln_fugcoef_L_()[spid] -= (A2_L_[spid] / (sqrt2 * B_L_) - A_L_ * Bi[spid] / (cc0 * B_L_2)) * ((dZdT_L_0 + cc1 * dBdT_L_) / (Z_L_ + cc1 * B_L_) - (dZdT_L_0 + cc2 * dBdT_L_) / (Z_L_ + cc2 * B_L_));
            }
        }

        for (int i = 0; i < this->X_.size(); i++)
            dvidt()[i] = (1 - vf_) * vf_ * (ddT_Ln_fugcoef_G_()[i] - ddT_Ln_fugcoef_L_()[i]);

        LUBacksubstitute(A, pivotIndices, dvidt());
    }

    autoPtr<scalarList> dXdT_G_(this->dXdT_gas(p, T, sol, dvidt));
    autoPtr<scalarList> dXdT_L_(this->dXdT_liq(p, T, sol, dvidt));
    scalar dvfdT_ = this->dvfdT(p, T, sol, dvidt());

    scalar dZdT_G_ = dZdT_G_0;
    scalar dZdT_L_ = dZdT_L_0;
    for (int i = 0; i < this->X_.size(); i++)
    {
        dZdT_G_ += dzdxi_G_[i] * dXdT_G_()[i];
        dZdT_L_ += dzdxi_L_[i] * dXdT_L_()[i];
    }

    scalar dZdT_M_ = dvfdT_ * (Z_G_ - Z_L_) + vf_ * dZdT_G_ + (1 - vf_) * dZdT_L_;

    scalar W_G_ = this->W(sol.X_gas());
    scalar W_L_ = this->W(sol.X_liq());

    scalar dmw_gasdT = this->dWdT_gas(p, T, sol, dvidt);

    scalar ygas = sol.vaporfra * W_G_ / W_;

    scalarList dHadxi_G_(this->N_);
    scalarList dHadxi_L_(this->N_);

    scalar ddxi_dH_G_, dH_G_;
    scalar ddxi_dH_L_, dH_L_;

    for (int i = 0; i < this->X_.size(); i++)
    {
        scalar dBmixdxi = Bi[i];

        if (xxx_G_ > 0)
        {
            dH_G_ = (Z_G_ - 1) * RR * 1e-3 * T + RR * 1e-3 * T * ((T * dAdT_G_ + A_G_) / (cc0 * B_G_)) * logxxx_G_;

            ddxi_dH_G_ = dzdxi_G_[i] * RR * 1e-3 * T;

            ddxi_dH_G_ += RR * 1e-3 * T * ((T * this->d2AdTdxi(p, T, i, sol.X_gas()) + this->dAdxi(p, T, i, sol.X_gas())) / (cc0 * B_G_)) * logxxx_G_;
            ddxi_dH_G_ -= RR * 1e-3 * T * ((T * dAdT_G_ + A_G_) / (cc0 * B_G_2) * dBmixdxi) * logxxx_G_;
            ddxi_dH_G_ += RR * 1e-3 * T * ((T * dAdT_G_ + A_G_) / (cc0 * B_G_)) * ((dzdxi_G_[i] + cc1 * dBmixdxi) / (Z_G_ + cc1 * B_G_) - (dzdxi_G_[i] + cc2 * dBmixdxi) / (Z_G_ + cc2 * B_G_));
            ddxi_dH_G_ = ddxi_dH_G_ / W_G_ - dH_G_ / (W_G_ * W_G_) * (*this)[i].W() * 1e-3;

            dHadxi_G_[i] = ddxi_dH_G_ + this->dHidealdxi(p, T, i, sol.X_gas());
        }
        else
        {
            ddxi_dH_G_ = dzdxi_G_[i] * RR * 1e-3 * T;
            dH_G_ = (Z_G_ - 1) * RR * 1e-3 * T;
            ddxi_dH_G_ = ddxi_dH_G_ / W_G_ - dH_G_ / (W_G_ * W_G_) * (*this)[i].W() * 1e-3;

            dHadxi_G_[i] = ddxi_dH_G_ + this->dHidealdxi(p, T, i, sol.X_gas());
        }
        if (xxx_L_ > 0)
        {
            dH_L_ = (Z_L_ - 1) * RR * 1e-3 * T + RR * 1e-3 * T * ((T * dAdT_L_ + A_L_) / (cc0 * B_L_)) * logxxx_L_;

            ddxi_dH_L_ = dzdxi_L_[i] * RR * 1e-3 * T;

            ddxi_dH_L_ += RR * 1e-3 * T * ((T * this->d2AdTdxi(p, T, i, sol.X_liq()) + this->dAdxi(p, T, i, sol.X_liq())) / (cc0 * B_L_)) * logxxx_L_;
            ddxi_dH_L_ -= RR * 1e-3 * T * ((T * dAdT_L_ + A_L_) / (cc0 * B_L_2) * dBmixdxi) * logxxx_L_;
            ddxi_dH_L_ += RR * 1e-3 * T * ((T * dAdT_L_ + A_L_) / (cc0 * B_L_)) * ((dzdxi_L_[i] + cc1 * dBmixdxi) / (Z_L_ + cc1 * B_L_) - (dzdxi_L_[i] + cc2 * dBmixdxi) / (Z_L_ + cc2 * B_L_));
            ddxi_dH_L_ = ddxi_dH_L_ / W_L_ - dH_L_ / (W_L_ * W_L_) * (*this)[i].W() * 1e-3;

            dHadxi_L_[i] = ddxi_dH_L_ + this->dHidealdxi(p, T, i, sol.X_liq());
        }
        else
        {
            ddxi_dH_L_ = dzdxi_L_[i] * RR * 1e-3 * T;
            dH_L_ = (Z_L_ - 1) * RR * 1e-3 * T;
            ddxi_dH_L_ = ddxi_dH_L_ / W_L_ - dH_L_ / (W_L_ * W_L_) * (*this)[i].W() * 1e-3;

            dHadxi_L_[i] = ddxi_dH_L_ + this->dHidealdxi(p, T, i, sol.X_liq());
        }
    }

    scalar dygasdT = (dvfdT_ * W_G_ + sol.vaporfra * dmw_gasdT) / W_;

    scalar bmix_G_ = this->b(p, T, sol.X_gas());

    scalar ddT_dH_G_;

    if (xxx_G_ > 0)
    {
        ddT_dH_G_ = dZdT_G_0 * RR * 1e-3 * T + (Z_G_ - 1) * RR * 1e-3;
        ddT_dH_G_ += ((T * this->d2adt2(p, T, sol.X_gas())) / (cc0 * bmix_G_)) * logxxx_G_;
        ddT_dH_G_ += ((T * this->dadt(p, T, sol.X_gas()) - this->a(p, T, sol.X_gas())) / (cc0 * bmix_G_)) * ((dZdT_G_0 - cc1 * B_G_ / T) / (Z_G_ + cc1 * B_G_) - (dZdT_G_0 - cc2 * B_G_ / T) / (Z_G_ + cc2 * B_G_));
    }
    else
    {
        ddT_dH_G_ = dZdT_G_0 * RR * 1e-3 * T + (Z_G_ - 1) * RR * 1e-3;
    }

    ddT_dH_G_ = ddT_dH_G_ / W_G_; // J/kg
                                  // Based on Eq.12, Internal energy is estimated based on ideal(close to absolute 0 T) + departure

    scalar bmix_L_ = this->b(p, T, sol.X_liq());

    scalar ddT_dH_L_;

    if (xxx_L_ > 0)
    {
        ddT_dH_L_ = dZdT_L_0 * RR * 1e-3 * T + (Z_L_ - 1) * RR * 1e-3;
        ddT_dH_L_ += ((T * this->d2adt2(p, T, sol.X_liq())) / (cc0 * bmix_L_)) * logxxx_L_;
        ddT_dH_L_ += ((T * this->dadt(p, T, sol.X_liq()) - this->a(p, T, sol.X_liq())) / (cc0 * bmix_L_)) * ((dZdT_L_0 - cc1 * B_L_ / T) / (Z_L_ + cc1 * B_L_) - (dZdT_L_0 - cc2 * B_L_ / T) / (Z_L_ + cc2 * B_L_));
    }
    else
    {
        ddT_dH_L_ = dZdT_L_0 * RR * 1e-3 * T + (Z_L_ - 1) * RR * 1e-3;
    }

    ddT_dH_L_ = ddT_dH_L_ / W_L_; // J/kg

    scalar dHadT_gas = ddT_dH_G_ + this->dHidealdT(p, T, sol.X_gas());

    scalar dHadT_liq = ddT_dH_L_ + this->dHidealdT(p, T, sol.X_liq());

    for (int i = 0; i < this->X_.size(); i++)
    {
        dHadT_gas += dHadxi_G_[i] * dXdT_G_()[i];
    }

    for (int i = 0; i < this->X_.size(); i++)
    {
        dHadT_liq += dHadxi_L_[i] * dXdT_L_()[i];
    }
    scalar Ha_G_ = dH_G_ / W_G_ + this->Hideal(p, T, sol.X_gas()), Ha_L_ = dH_L_ / W_L_ + this->Hideal(p, T, sol.X_liq());

    scalar Cp_ = Ha_G_ * dygasdT + dHadT_gas * ygas - Ha_L_ * dygasdT + dHadT_liq * (1.0 - ygas);

    //scalar Cp_ = this->ThermoMixture::Ha(p, T, sol.X_gas(), 1) * dygasdT + dHadT_gas * ygas - this->ThermoMixture::Ha(p, T, sol.X_liq(), 0) * dygasdT + dHadT_liq * (1.0 - ygas);

    scalar drhodT_ = -p * W_ / (Z_M_ * Z_M_ * RR * 1.0e-03 * T * T) * (dZdT_M_ * T + Z_M_);
    scalar alphaP_ = -drhodT_ / rho_;
    scalar kappaS_ = kappaT_ - T * alphaP_ * alphaP_ / (rho_ * Cp_);
    if (kappaS_ < 0)
    {
        //return c_noVLE(p, T, sol);
        FatalErrorInFunction << " kappaS_ =" << kappaS_ << ",p=" << p << ",T=" << T << ",X=" << this->X_ << exit(FatalError);
    }
    return sqrt(1.0 / (kappaS_ * rho_));
}

template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::c_noVLE(scalar p, scalar T, solution &sol) const
{
    scalar rho_ = rho_noVLE(p, T, sol);
    scalar kappaS_ = kappaS_noVLE(p, T, sol);
    return sqrt(1.0 / (kappaS_ * rho_));
}
template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::c(scalar p, scalar T) const
{
    autoPtr<solution> sol(TPN(p, T));
    return c(p, T, sol());
}
template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::c_noVLE(scalar p, scalar T) const
{
    autoPtr<solution> sol(TPN(p, T));
    return c_noVLE(p, T, sol());
}
// ************************************************************************* //
