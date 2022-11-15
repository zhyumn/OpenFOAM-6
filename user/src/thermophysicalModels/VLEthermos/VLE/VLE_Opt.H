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

// ************************************************************************* //

template <class ThermoMixture>
inline Foam::autoPtr<scalarRectangularMatrix> Foam::VLE<ThermoMixture>::dTHvf_G_rhoY_dXrhoP(scalar p, scalar T, solution &sol) const
{
    int nsp = this->X_.size();
    //scalar temp;
    scalarRectangularMatrix &grad = *(new scalarRectangularMatrix(nsp + 2, 4 + nsp));

    const VLE<ThermoMixture> &This = *this;
    scalar vf = sol.vaporfra;
    scalarList &X_gas = sol.X_gas();
    scalarList &X_liq = sol.X_liq();
    scalar rho_value = rho(p, T, sol);

    autoPtr<scalarSquareMatrix> pA(A_dvid_(p, T, sol));
    labelList pivotIndices(pA().m());
    LUDecompose(pA(), pivotIndices);

    autoPtr<scalarList> dvidP_value(this->dvidP(p, T, sol, pA(), pivotIndices));
    autoPtr<scalarList> dvidT_value(this->dvidT(p, T, sol, pA(), pivotIndices));
    List<autoPtr<scalarList>> dvidXi_value(nsp);
    for (int i = 0; i < nsp; i++)
    {
        dvidXi_value[i] = this->dvidXi(p, T, i, sol);
    }
    scalar dvfdP_value = 0, dvfdT_value = 0;
    scalarList dvfdXi_value(nsp);
    scalarList dHsdxi_gas_value(nsp), dHsdxi_liq_value(nsp);
    for (int i = 0; i < nsp; i++)
    {
        dvfdT_value += dvidT_value()[i];
        dvfdP_value += dvidP_value()[i];
        dvfdXi_value[i] = -vf;
        for (int j = 0; j < nsp; j++)
        {
            dvfdXi_value[i] += dvidXi_value[i]()[j];
        }
        dHsdxi_gas_value[i] = this->dHsdxi(p, T, i, X_gas, 1);
        dHsdxi_liq_value[i] = this->dHsdxi(p, T, i, X_liq, 0);
    }

    scalar dHsdT_gas = this->ThermoMixture::dHsdT(p, T, X_gas, 1);
    scalar dHsdT_liq = this->ThermoMixture::dHsdT(p, T, X_liq, 0);
    scalar dHsdP_gas = this->ThermoMixture::dHsdP(p, T, X_gas, 1);
    scalar dHsdP_liq = this->ThermoMixture::dHsdP(p, T, X_liq, 0);
    scalarList dHsdXi_gas(nsp);
    scalarList dHsdXi_liq(nsp);

    scalarList dXdT_gas_value(nsp);
    scalarList dXdT_liq_value(nsp);
    scalarList dXdP_gas_value(nsp);
    scalarList dXdP_liq_value(nsp);
    if (vf <= 1e-6)
    {
        for (int i = 0; i < nsp; i++)
        {
            dXdT_gas_value[i] = 0;
            dXdT_liq_value[i] = -(dvidT_value()[i] - X_liq[i] * dvfdT_value) / (1 - vf);
            dXdP_gas_value[i] = 0;
            dXdP_liq_value[i] = -(dvidP_value()[i] - X_liq[i] * dvfdP_value) / (1 - vf);
        }
    }
    else if (vf >= 1 - 1e-6)
    {
        for (int i = 0; i < nsp; i++)
        {
            dXdT_gas_value[i] = (dvidT_value()[i] - X_gas[i] * dvfdT_value) / vf;
            dXdT_liq_value[i] = 0;
            dXdP_gas_value[i] = (dvidP_value()[i] - X_gas[i] * dvfdP_value) / vf;
            dXdP_liq_value[i] = 0;
        }
    }
    else
    {
        for (int i = 0; i < nsp; i++)
        {
            dXdT_gas_value[i] = (dvidT_value()[i] - X_gas[i] * dvfdT_value) / vf;
            dXdT_liq_value[i] = -(dvidT_value()[i] - X_liq[i] * dvfdT_value) / (1 - vf);
            dXdP_gas_value[i] = (dvidP_value()[i] - X_gas[i] * dvfdP_value) / vf;
            dXdP_liq_value[i] = -(dvidP_value()[i] - X_liq[i] * dvfdP_value) / (1 - vf);
        }
    }

    List<scalarList> dXdXi_gas_value(nsp);
    List<scalarList> dXdXi_liq_value(nsp);

    if (vf <= 1e-6)
    {
        for (int i = 0; i < nsp; i++)
        {
            dXdXi_gas_value[i].resize(nsp);
            dXdXi_liq_value[i].resize(nsp);
            for (int j = 0; j < nsp; j++)
            {
                dXdXi_gas_value[i][j] = 0;
                if (i == j)
                {
                    dXdXi_liq_value[i][j] = 1 / (1 - vf) * (1 - dvidXi_value[i]()[j]) - X_liq[j] / (1 - vf) * (1 - dvfdXi_value[i] - vf);
                }
                else
                {
                    dXdXi_liq_value[i][j] = 1 / (1 - vf) * (-dvidXi_value[i]()[j]) - X_liq[j] / (1 - vf) * (1 - dvfdXi_value[i] - vf);
                }
            }
        }
    }
    else if (vf >= 1 - 1e-6)
    {
        for (int i = 0; i < nsp; i++)
        {
            dXdXi_gas_value[i].resize(nsp);
            dXdXi_liq_value[i].resize(nsp);
            for (int j = 0; j < nsp; j++)
            {
                dXdXi_gas_value[i][j] = (dvidXi_value[i]()[j] - X_gas[j] * (dvfdXi_value[i] + vf)) / vf;
                dXdXi_liq_value[i][j] = 0;
            }
        }
    }
    else
    {
        for (int i = 0; i < nsp; i++)
        {
            dXdXi_gas_value[i].resize(nsp);
            dXdXi_liq_value[i].resize(nsp);
            for (int j = 0; j < nsp; j++)
            {
                dXdXi_gas_value[i][j] = (dvidXi_value[i]()[j] - X_gas[j] * (dvfdXi_value[i] + vf)) / vf;
                if (i == j)
                {
                    dXdXi_liq_value[i][j] = 1 / (1 - vf) * (1 - dvidXi_value[i]()[j]) - X_liq[j] / (1 - vf) * (1 - dvfdXi_value[i] - vf);
                }
                else
                {
                    dXdXi_liq_value[i][j] = 1 / (1 - vf) * (-dvidXi_value[i]()[j]) - X_liq[j] / (1 - vf) * (1 - dvfdXi_value[i] - vf);
                }
            }
        }
    }

    for (int i = 0; i < nsp; i++)
    {
        dHsdT_gas += dHsdxi_gas_value[i] * dXdT_gas_value[i];
        dHsdT_liq += dHsdxi_liq_value[i] * dXdT_liq_value[i];
        dHsdP_gas += dHsdxi_gas_value[i] * dXdP_gas_value[i];
        dHsdP_liq += dHsdxi_liq_value[i] * dXdP_liq_value[i];
        dHsdXi_gas[i] = 0;
        dHsdXi_liq[i] = 0;
        for (int j = 0; j < nsp; j++)
        {
            dHsdXi_gas[i] += dHsdxi_gas_value[j] * dXdXi_gas_value[i][j];
            dHsdXi_liq[i] += dHsdxi_liq_value[j] * dXdXi_liq_value[i][j];
        }
    }

    scalar mw_gas = this->W(X_gas); //kg/mol
    scalar mw_mixture = this->W();
    scalar ygas = vf * mw_gas / mw_mixture;
    scalar dmw_gasdT = 0;
    scalar dmw_gasdP = 0;
    for (int i = 0; i < nsp; i++)
    {
        dmw_gasdT += dXdT_gas_value[i] * This[i].W() * 1.0e-03;
        dmw_gasdP += dXdP_gas_value[i] * This[i].W() * 1.0e-03;
    }
    //Info << "dmw_gasdT" << dmw_gasdT << endl;
    scalarList dmw_gasdXi(nsp);
    //Info << "mw_gas" << This[0].W() * 1.0e-03 << endl;
    for (int i = 0; i < nsp; i++)
    {
        dmw_gasdXi[i] = 0;
        for (int j = 0; j < nsp; j++)
        {
            dmw_gasdXi[i] += dXdXi_gas_value[i][j] * This[j].W() * 1.0e-03;
        }
    }
    scalar dygasdT = (dvfdT_value * mw_gas + vf * dmw_gasdT) / mw_mixture;
    scalar dygasdP = (dvfdP_value * mw_gas + vf * dmw_gasdP) / mw_mixture;
    scalarList dygasdXi(nsp);
    for (int i = 0; i < nsp; i++)
    {
        dygasdXi[i] = ((dvfdXi_value[i] + vf) * mw_gas + vf * dmw_gasdXi[i]) / mw_mixture - ygas / mw_mixture * This[i].W() * 1.0e-03;
    }

    scalar Hs_gas = this->ThermoMixture::Hs(p, T, X_gas, 1);
    scalar Hs_liq = this->ThermoMixture::Hs(p, T, X_liq, 0);
    scalar Cp_Hs_value = Hs_gas * dygasdT + dHsdT_gas * ygas - Hs_liq * dygasdT + dHsdT_liq * (1.0 - ygas);
    scalar dHsdP_value = Hs_gas * dygasdP + dHsdP_gas * ygas - Hs_liq * dygasdP + dHsdP_liq * (1.0 - ygas);
    scalarList dHsdXi_value(nsp);
    for (int i = 0; i < nsp; i++)
    {
        dHsdXi_value[i] = Hs_gas * dygasdXi[i] + dHsdXi_gas[i] * ygas - Hs_liq * dygasdXi[i] + dHsdXi_liq[i] * (1.0 - ygas);
    }

    scalar ZMixture = Z(p, T, sol);
    scalar ZGas = this->ThermoMixture::Z(p, T, X_gas, 1);
    scalar ZLiq = this->ThermoMixture::Z(p, T, X_liq, 0);
    scalar dZgasdT = this->ThermoMixture::dZdT(p, T, X_gas, 1);
    scalar dZliqdT = this->ThermoMixture::dZdT(p, T, X_liq, 0);
    scalar dZgasdP = this->ThermoMixture::dZdP(p, T, X_gas, 1);
    scalar dZliqdP = this->ThermoMixture::dZdP(p, T, X_liq, 0);
    scalarList dZgasdXi(nsp), dZliqdXi(nsp);
    scalarList dZdxi_gas_value(nsp), dZdxi_liq_value(nsp);
    for (int i = 0; i < nsp; i++)
    {
        dZdxi_gas_value[i] = this->ThermoMixture::dZdxi(p, T, i, X_gas, 1);
        dZdxi_liq_value[i] = this->ThermoMixture::dZdxi(p, T, i, X_liq, 0);
    }
    for (int i = 0; i < nsp; i++)
    {
        dZgasdT += dZdxi_gas_value[i] * dXdT_gas_value[i];
        dZliqdT += dZdxi_liq_value[i] * dXdT_liq_value[i];
        dZgasdP += dZdxi_gas_value[i] * dXdP_gas_value[i];
        dZliqdP += dZdxi_liq_value[i] * dXdP_liq_value[i];
        dZgasdXi[i] = 0;
        dZliqdXi[i] = 0;
        for (int j = 0; j < nsp; j++)
        {
            dZgasdXi[i] += dZdxi_gas_value[j] * dXdXi_gas_value[i][j];
            dZliqdXi[i] += dZdxi_liq_value[j] * dXdXi_liq_value[i][j];
        }
    }
    //temp = (vf * ZGas + (1 - vf) * ZLiq);
    scalar dZMixturedT = dvfdT_value * (ZGas - ZLiq) + vf * dZgasdT + (1 - vf) * dZliqdT;
    scalar dZMixturedP = dvfdP_value * (ZGas - ZLiq) + vf * dZgasdP + (1 - vf) * dZliqdP;
    scalarList dZMixturedXi(nsp);

    for (int i = 0; i < nsp; i++)
    {
        //dZMixturedXi[i] = (dvfdXi_value[i] * ZGas * ZGas + vf * 2 * ZGas * dZgasdXi[i] + 2 * (1 - vf) * ZLiq * dZliqdXi[i] - dvfdXi_value[i] * ZLiq * ZLiq) / temp - (vf * ZGas * ZGas + (1 - vf) * ZLiq * ZLiq) / (temp * temp) * (dvfdXi_value[i] * ZGas + vf * dZgasdXi[i] + (1 - vf) * dZliqdXi[i] - dvfdXi_value[i] * ZLiq);
        dZMixturedXi[i] = dvfdXi_value[i] * (ZGas - ZLiq) + vf * dZgasdXi[i] + (1 - vf) * dZliqdXi[i];
    }

    scalar drhodT_value = -p * mw_mixture / (ZMixture * ZMixture * RR * 1.0e-03 * T * T) * (dZMixturedT * T + ZMixture);
    scalar drhodP_value = mw_mixture / (ZMixture * RR * 1.0e-03 * T) - p * mw_mixture / (ZMixture * ZMixture * RR * 1.0e-03 * T) * dZMixturedP;
    scalarList drhodXi_value(nsp);
    for (int i = 0; i < nsp; i++)
    {
        drhodXi_value[i] = p * (This[i].W() * 1.0e-03 - mw_mixture) / (ZMixture * RR * 1.0e-03 * T) - p * mw_mixture / (ZMixture * ZMixture * RR * 1.0e-03 * T) * dZMixturedXi[i];
    }

    scalar drhodT_G_value = -p * mw_gas / (ZGas * ZGas * RR * 1.0e-03 * T * T) * (dZgasdT * T + ZGas) + p / (ZGas * RR * 1.0e-03 * T) * dmw_gasdT;
    scalar drhodP_G_value = mw_gas / (ZGas * RR * 1.0e-03 * T) - p * mw_gas / (ZGas * ZGas * RR * 1.0e-03 * T) * dZgasdP + p / (ZGas * RR * 1.0e-03 * T) * dmw_gasdP;

    scalarList drhodXi_G_value(nsp);
    for (int i = 0; i < nsp; i++)
    {
        drhodXi_G_value[i] = -p * mw_gas / (ZGas * ZGas * RR * 1.0e-03 * T) * dZgasdXi[i] + p / (ZGas * RR * 1.0e-03 * T) * dmw_gasdXi[i];
    }

    scalarList dYdT_G_value(nsp), dYdP_G_value(nsp);
    for (int i = 0; i < nsp; i++)
    {
        dYdT_G_value[i] = dXdT_gas_value[i] * This[i].W() * 1.0e-03 / mw_gas - sol.X_gas()[i] * This[i].W() * 1.0e-03 / (mw_gas * mw_gas) * dmw_gasdT;
        dYdP_G_value[i] = dXdP_gas_value[i] * This[i].W() * 1.0e-03 / mw_gas - sol.X_gas()[i] * This[i].W() * 1.0e-03 / (mw_gas * mw_gas) * dmw_gasdP;
    }

    List<scalarList> dYdXi_G_value(nsp);
    for (int i = 0; i < nsp; i++)
    {
        dYdXi_G_value[i].resize(nsp);
        for (int j = 0; j < nsp; j++)
        {
            dYdXi_G_value[i][j] = dXdXi_gas_value[i][j] * This[j].W() * 1.0e-03 / mw_gas - sol.X_gas()[j] * This[j].W() * 1.0e-03 / (mw_gas * mw_gas) * dmw_gasdXi[i];
        }
    }

    for (int i = 0; i < nsp; i++)
    {
        grad[i][0] = -drhodXi_value[i] / drhodT_value;
        grad[i][1] = -Cp_Hs_value * drhodXi_value[i] / drhodT_value + dHsdXi_value[i];
        grad[i][2] = -dvfdT_value * drhodXi_value[i] / drhodT_value + dvfdXi_value[i];
        grad[i][3] = -drhodT_G_value * drhodXi_value[i] / drhodT_value + drhodXi_G_value[i];
        for (int j = 0; j < nsp; j++)
        {
            grad[i][4 + j] = -dYdT_G_value[j] * drhodXi_value[i] / drhodT_value + dYdXi_G_value[i][j];
        }
    }
    grad[nsp][0] = 1 / drhodT_value;
    grad[nsp][1] = Cp_Hs_value / drhodT_value;
    grad[nsp][2] = dvfdT_value / drhodT_value;
    grad[nsp][3] = drhodT_G_value / drhodT_value;
    for (int j = 0; j < nsp; j++)
    {
        grad[nsp][4 + j] = dYdT_G_value[j] / drhodT_value;
    }

    grad[nsp + 1][0] = -drhodP_value / drhodT_value;
    grad[nsp + 1][1] = -Cp_Hs_value * drhodP_value / drhodT_value + dHsdP_value;
    grad[nsp + 1][2] = -dvfdT_value * drhodP_value / drhodT_value + dvfdP_value;
    grad[nsp + 1][3] = -drhodT_G_value * drhodP_value / drhodT_value + drhodP_G_value;
    for (int j = 0; j < nsp; j++)
    {
        grad[nsp + 1][4 + j] = -dYdT_G_value[j] * drhodP_value / drhodT_value + dYdP_G_value[j];
    }

    return autoPtr<scalarRectangularMatrix>(&grad);
}

template <class ThermoMixture>
inline Foam::autoPtr<scalarRectangularMatrix> Foam::VLE<ThermoMixture>::dErhovf_G_rhoY_dXTP(scalar p, scalar T, solution &sol) const
{
    int nsp = this->X_.size();
    //scalar temp;
    scalarRectangularMatrix &grad = *(new scalarRectangularMatrix(nsp + 2, 4 + nsp));

    const VLE<ThermoMixture> &This = *this;
    scalar vf = sol.vaporfra;
    scalarList &X_gas = sol.X_gas();
    scalarList &X_liq = sol.X_liq();
    scalar rho_value = rho(p, T, sol);

    autoPtr<scalarSquareMatrix> pA(A_dvid_(p, T, sol));
    labelList pivotIndices(pA().m());
    LUDecompose(pA(), pivotIndices);

    autoPtr<scalarList> dvidP_value(this->dvidP(p, T, sol, pA(), pivotIndices));
    autoPtr<scalarList> dvidT_value(this->dvidT(p, T, sol, pA(), pivotIndices));
    List<autoPtr<scalarList>> dvidXi_value(nsp);
    for (int i = 0; i < nsp; i++)
    {
        dvidXi_value[i] = this->dvidXi(p, T, i, sol);
    }
    scalar dvfdP_value = 0, dvfdT_value = 0;
    scalarList dvfdXi_value(nsp);
    scalarList dHsdxi_gas_value(nsp), dHsdxi_liq_value(nsp);

    for (int i = 0; i < nsp; i++)
    {
        dvfdT_value += dvidT_value()[i];
        dvfdP_value += dvidP_value()[i];
        dvfdXi_value[i] = -vf;
        for (int j = 0; j < nsp; j++)
        {
            dvfdXi_value[i] += dvidXi_value[i]()[j];
        }
        dHsdxi_gas_value[i] = this->dHsdxi(p, T, i, X_gas, 1);
        dHsdxi_liq_value[i] = this->dHsdxi(p, T, i, X_liq, 0);
    }

    scalar dHsdT_gas = this->ThermoMixture::dHsdT(p, T, X_gas, 1);
    scalar dHsdT_liq = this->ThermoMixture::dHsdT(p, T, X_liq, 0);
    scalar dHsdP_gas = this->ThermoMixture::dHsdP(p, T, X_gas, 1);
    scalar dHsdP_liq = this->ThermoMixture::dHsdP(p, T, X_liq, 0);
    scalarList dHsdXi_gas(nsp);
    scalarList dHsdXi_liq(nsp);

    scalarList dXdT_gas_value(nsp);
    scalarList dXdT_liq_value(nsp);
    scalarList dXdP_gas_value(nsp);
    scalarList dXdP_liq_value(nsp);
    if (vf <= 1e-6)
    {
        for (int i = 0; i < nsp; i++)
        {
            dXdT_gas_value[i] = 0;
            dXdT_liq_value[i] = -(dvidT_value()[i] - X_liq[i] * dvfdT_value) / (1 - vf);
            dXdP_gas_value[i] = 0;
            dXdP_liq_value[i] = -(dvidP_value()[i] - X_liq[i] * dvfdP_value) / (1 - vf);
        }
    }
    else if (vf >= 1 - 1e-6)
    {
        for (int i = 0; i < nsp; i++)
        {
            dXdT_gas_value[i] = (dvidT_value()[i] - X_gas[i] * dvfdT_value) / vf;
            dXdT_liq_value[i] = 0;
            dXdP_gas_value[i] = (dvidP_value()[i] - X_gas[i] * dvfdP_value) / vf;
            dXdP_liq_value[i] = 0;
        }
    }
    else
    {
        for (int i = 0; i < nsp; i++)
        {
            dXdT_gas_value[i] = (dvidT_value()[i] - X_gas[i] * dvfdT_value) / vf;
            dXdT_liq_value[i] = -(dvidT_value()[i] - X_liq[i] * dvfdT_value) / (1 - vf);
            dXdP_gas_value[i] = (dvidP_value()[i] - X_gas[i] * dvfdP_value) / vf;
            dXdP_liq_value[i] = -(dvidP_value()[i] - X_liq[i] * dvfdP_value) / (1 - vf);
        }
    }

    List<scalarList> dXdXi_gas_value(nsp);
    List<scalarList> dXdXi_liq_value(nsp);

    if (vf <= 1e-6)
    {
        for (int i = 0; i < nsp; i++)
        {
            dXdXi_gas_value[i].resize(nsp);
            dXdXi_liq_value[i].resize(nsp);
            for (int j = 0; j < nsp; j++)
            {
                dXdXi_gas_value[i][j] = 0;
                if (i == j)
                {
                    dXdXi_liq_value[i][j] = 1 / (1 - vf) * (1 - dvidXi_value[i]()[j]) - X_liq[j] / (1 - vf) * (1 - dvfdXi_value[i] - vf);
                }
                else
                {
                    dXdXi_liq_value[i][j] = 1 / (1 - vf) * (-dvidXi_value[i]()[j]) - X_liq[j] / (1 - vf) * (1 - dvfdXi_value[i] - vf);
                }
            }
        }
    }
    else if (vf >= 1 - 1e-6)
    {
        for (int i = 0; i < nsp; i++)
        {
            dXdXi_gas_value[i].resize(nsp);
            dXdXi_liq_value[i].resize(nsp);
            for (int j = 0; j < nsp; j++)
            {
                dXdXi_gas_value[i][j] = (dvidXi_value[i]()[j] - X_gas[j] * (dvfdXi_value[i] + vf)) / vf;
                dXdXi_liq_value[i][j] = 0;
            }
        }
    }
    else
    {
        for (int i = 0; i < nsp; i++)
        {
            dXdXi_gas_value[i].resize(nsp);
            dXdXi_liq_value[i].resize(nsp);
            for (int j = 0; j < nsp; j++)
            {
                dXdXi_gas_value[i][j] = (dvidXi_value[i]()[j] - X_gas[j] * (dvfdXi_value[i] + vf)) / vf;
                if (i == j)
                {
                    dXdXi_liq_value[i][j] = 1 / (1 - vf) * (1 - dvidXi_value[i]()[j]) - X_liq[j] / (1 - vf) * (1 - dvfdXi_value[i] - vf);
                }
                else
                {
                    dXdXi_liq_value[i][j] = 1 / (1 - vf) * (-dvidXi_value[i]()[j]) - X_liq[j] / (1 - vf) * (1 - dvfdXi_value[i] - vf);
                }
            }
        }
    }

    for (int i = 0; i < nsp; i++)
    {
        dHsdT_gas += dHsdxi_gas_value[i] * dXdT_gas_value[i];
        dHsdT_liq += dHsdxi_liq_value[i] * dXdT_liq_value[i];
        dHsdP_gas += dHsdxi_gas_value[i] * dXdP_gas_value[i];
        dHsdP_liq += dHsdxi_liq_value[i] * dXdP_liq_value[i];
        dHsdXi_gas[i] = 0;
        dHsdXi_liq[i] = 0;
        for (int j = 0; j < nsp; j++)
        {
            dHsdXi_gas[i] += dHsdxi_gas_value[j] * dXdXi_gas_value[i][j];
            dHsdXi_liq[i] += dHsdxi_liq_value[j] * dXdXi_liq_value[i][j];
        }
    }

    scalar mw_gas = this->W(X_gas); //kg/mol
    scalar mw_mixture = this->W();
    scalar ygas = vf * mw_gas / mw_mixture;
    scalar dmw_gasdT = 0;
    scalar dmw_gasdP = 0;
    for (int i = 0; i < nsp; i++)
    {
        dmw_gasdT += dXdT_gas_value[i] * This[i].W() * 1.0e-03;
        dmw_gasdP += dXdP_gas_value[i] * This[i].W() * 1.0e-03;
    }
    scalarList dmw_gasdXi(nsp);
    for (int i = 0; i < nsp; i++)
    {
        dmw_gasdXi[i] = 0;
        for (int j = 0; j < nsp; j++)
        {
            dmw_gasdXi[i] += dXdXi_gas_value[i][j] * This[j].W() * 1.0e-03;
        }
    }
    scalar dygasdT = (dvfdT_value * mw_gas + vf * dmw_gasdT) / mw_mixture;
    scalar dygasdP = (dvfdP_value * mw_gas + vf * dmw_gasdP) / mw_mixture;
    scalarList dygasdXi(nsp);
    for (int i = 0; i < nsp; i++)
    {
        dygasdXi[i] = ((dvfdXi_value[i] + vf) * mw_gas + vf * dmw_gasdXi[i]) / mw_mixture - ygas / mw_mixture * This[i].W() * 1.0e-03;
    }

    scalar Hs_gas = this->ThermoMixture::Hs(p, T, X_gas, 1);
    scalar Hs_liq = this->ThermoMixture::Hs(p, T, X_liq, 0);
    scalar Cp_Hs_value = Hs_gas * dygasdT + dHsdT_gas * ygas - Hs_liq * dygasdT + dHsdT_liq * (1.0 - ygas);
    scalar dHsdP_value = Hs_gas * dygasdP + dHsdP_gas * ygas - Hs_liq * dygasdP + dHsdP_liq * (1.0 - ygas);
    scalarList dHsdXi_value(nsp);
    for (int i = 0; i < nsp; i++)
    {
        dHsdXi_value[i] = Hs_gas * dygasdXi[i] + dHsdXi_gas[i] * ygas - Hs_liq * dygasdXi[i] + dHsdXi_liq[i] * (1.0 - ygas);
    }

    scalar ZMixture = Z(p, T, sol);
    scalar ZGas = this->ThermoMixture::Z(p, T, X_gas, 1);
    scalar ZLiq = this->ThermoMixture::Z(p, T, X_liq, 0);
    scalar dZgasdT = this->ThermoMixture::dZdT(p, T, X_gas, 1);
    scalar dZliqdT = this->ThermoMixture::dZdT(p, T, X_liq, 0);
    scalar dZgasdP = this->ThermoMixture::dZdP(p, T, X_gas, 1);
    scalar dZliqdP = this->ThermoMixture::dZdP(p, T, X_liq, 0);
    scalarList dZgasdXi(nsp), dZliqdXi(nsp);
    scalarList dZdxi_gas_value(nsp), dZdxi_liq_value(nsp);
    for (int i = 0; i < nsp; i++)
    {
        dZdxi_gas_value[i] = this->ThermoMixture::dZdxi(p, T, i, X_gas, 1);
        dZdxi_liq_value[i] = this->ThermoMixture::dZdxi(p, T, i, X_liq, 0);
    }
    for (int i = 0; i < nsp; i++)
    {
        dZgasdT += dZdxi_gas_value[i] * dXdT_gas_value[i];
        dZliqdT += dZdxi_liq_value[i] * dXdT_liq_value[i];
        dZgasdP += dZdxi_gas_value[i] * dXdP_gas_value[i];
        dZliqdP += dZdxi_liq_value[i] * dXdP_liq_value[i];
        dZgasdXi[i] = 0;
        dZliqdXi[i] = 0;
        for (int j = 0; j < nsp; j++)
        {
            dZgasdXi[i] += dZdxi_gas_value[j] * dXdXi_gas_value[i][j];
            dZliqdXi[i] += dZdxi_liq_value[j] * dXdXi_liq_value[i][j];
        }
    }
    //temp = (vf * ZGas + (1 - vf) * ZLiq);
    scalar dZMixturedT = dvfdT_value * (ZGas - ZLiq) + vf * dZgasdT + (1 - vf) * dZliqdT;
    scalar dZMixturedP = dvfdP_value * (ZGas - ZLiq) + vf * dZgasdP + (1 - vf) * dZliqdP;
    scalarList dZMixturedXi(nsp);

    for (int i = 0; i < nsp; i++)
    {
        dZMixturedXi[i] = dvfdXi_value[i] * (ZGas - ZLiq) + vf * dZgasdXi[i] + (1 - vf) * dZliqdXi[i];
    }

    scalar drhodT_value = -p * mw_mixture / (ZMixture * ZMixture * RR * 1.0e-03 * T * T) * (dZMixturedT * T + ZMixture);
    scalar drhodP_value = mw_mixture / (ZMixture * RR * 1.0e-03 * T) - p * mw_mixture / (ZMixture * ZMixture * RR * 1.0e-03 * T) * dZMixturedP;
    scalarList drhodXi_value(nsp);
    for (int i = 0; i < nsp; i++)
    {
        drhodXi_value[i] = p * (This[i].W() * 1.0e-03 - mw_mixture) / (ZMixture * RR * 1.0e-03 * T) - p * mw_mixture / (ZMixture * ZMixture * RR * 1.0e-03 * T) * dZMixturedXi[i];
    }

    scalar drhodT_G_value = -p * mw_gas / (ZGas * ZGas * RR * 1.0e-03 * T * T) * (dZgasdT * T + ZGas) + p / (ZGas * RR * 1.0e-03 * T) * dmw_gasdT;
    scalar drhodP_G_value = mw_gas / (ZGas * RR * 1.0e-03 * T) - p * mw_gas / (ZGas * ZGas * RR * 1.0e-03 * T) * dZgasdP + p / (ZGas * RR * 1.0e-03 * T) * dmw_gasdP;

    scalarList drhodXi_G_value(nsp);
    for (int i = 0; i < nsp; i++)
    {
        //drhodXi_G_value[i] = p * (This[i].W() * 1.0e-03 - mw_gas) / (ZGas * RR * 1.0e-03 * T) - p * mw_gas / (ZGas * ZGas * RR * 1.0e-03 * T) * dZgasdXi[i] + p / (ZGas * RR * 1.0e-03 * T) * dmw_gasdXi[i];
        drhodXi_G_value[i] = -p * mw_gas / (ZGas * ZGas * RR * 1.0e-03 * T) * dZgasdXi[i] + p / (ZGas * RR * 1.0e-03 * T) * dmw_gasdXi[i];
    }

    for (int i = 0; i < nsp; i++)
    {
        grad[i][0] = dHsdXi_value[i] + p / (rho_value * rho_value) * drhodXi_value[i];
        grad[i][1] = drhodXi_value[i];
        grad[i][2] = dvfdXi_value[i];
        grad[i][3] = drhodXi_G_value[i];

        for (int j = 0; j < nsp; j++)
        {
            grad[i][4 + j] = dXdXi_gas_value[i][j] * This[j].W() * 1.0e-03 / mw_gas - sol.X_gas()[j] * This[j].W() * 1.0e-03 / (mw_gas * mw_gas) * dmw_gasdXi[i];
        }
    }

    grad[nsp][0] = Cp_Hs_value + p / (rho_value * rho_value) * drhodT_value;
    grad[nsp][1] = drhodT_value;
    grad[nsp][2] = dvfdT_value;
    grad[nsp][3] = drhodT_G_value;
    //Info << "dmw_gasdT=" << dmw_gasdT<<",Y0=" << sol.X_gas()[0] * This[0].W() * 1.0e-03 / mw_gas<<endl;
    for (int j = 0; j < nsp; j++)
    {
        grad[nsp][4 + j] = dXdT_gas_value[j] * This[j].W() * 1.0e-03 / mw_gas - sol.X_gas()[j] * This[j].W() * 1.0e-03 / (mw_gas * mw_gas) * dmw_gasdT;
    }
    grad[nsp + 1][0] = dHsdP_value - 1 / rho_value + p / (rho_value * rho_value) * drhodP_value;
    grad[nsp + 1][1] = drhodP_value;
    grad[nsp + 1][2] = dvfdP_value;
    grad[nsp + 1][3] = drhodP_G_value;
    for (int j = 0; j < nsp; j++)
    {
        grad[nsp + 1][4 + j] = dXdP_gas_value[j] * This[j].W() * 1.0e-03 / mw_gas - sol.X_gas()[j] * This[j].W() * 1.0e-03 / (mw_gas * mw_gas) * dmw_gasdP;
    }
    return autoPtr<scalarRectangularMatrix>(&grad);
}

template <class ThermoMixture>
inline Foam::autoPtr<scalarRectangularMatrix> Foam::VLE<ThermoMixture>::dTrhovfc_G_rhoY_dXHsP(scalar p, scalar T, solution &sol) const
{
    int nsp = this->X_.size();
    scalar temp;
    scalarRectangularMatrix &grad = *(new scalarRectangularMatrix(nsp + 2, 5 + nsp));

    const VLE<ThermoMixture> &This = *this;
    scalar vf = sol.vaporfra;
    scalarList &X_gas = sol.X_gas();
    scalarList &X_liq = sol.X_liq();
    scalar rho_value = rho(p, T, sol);

    autoPtr<scalarSquareMatrix> pA(A_dvid_(p, T, sol));
    labelList pivotIndices(pA().m());
    LUDecompose(pA(), pivotIndices);

    autoPtr<scalarList> dvidP_value(this->dvidP(p, T, sol, pA(), pivotIndices));
    autoPtr<scalarList> dvidT_value(this->dvidT(p, T, sol, pA(), pivotIndices));
    List<autoPtr<scalarList>> dvidXi_value(nsp);
    for (int i = 0; i < nsp; i++)
    {
        dvidXi_value[i] = this->dvidXi(p, T, i, sol);
    }
    scalar dvfdP_value = 0, dvfdT_value = 0;
    scalarList dvfdXi_value(nsp);
    scalarList dHsdxi_gas_value(nsp), dHsdxi_liq_value(nsp);

    for (int i = 0; i < nsp; i++)
    {
        dvfdT_value += dvidT_value()[i];
        dvfdP_value += dvidP_value()[i];
        dvfdXi_value[i] = -vf;
        for (int j = 0; j < nsp; j++)
        {
            dvfdXi_value[i] += dvidXi_value[i]()[j];
        }
        dHsdxi_gas_value[i] = this->dHsdxi(p, T, i, X_gas, 1);
        dHsdxi_liq_value[i] = this->dHsdxi(p, T, i, X_liq, 0);
    }

    scalar dHsdT_gas = this->ThermoMixture::dHsdT(p, T, X_gas, 1);
    scalar dHsdT_liq = this->ThermoMixture::dHsdT(p, T, X_liq, 0);
    scalar dHsdP_gas = this->ThermoMixture::dHsdP(p, T, X_gas, 1);
    scalar dHsdP_liq = this->ThermoMixture::dHsdP(p, T, X_liq, 0);
    scalarList dHsdXi_gas(nsp);
    scalarList dHsdXi_liq(nsp);

    scalarList dXdT_gas_value(nsp);
    scalarList dXdT_liq_value(nsp);
    scalarList dXdP_gas_value(nsp);
    scalarList dXdP_liq_value(nsp);
    if (vf <= 1e-6)
    {
        for (int i = 0; i < nsp; i++)
        {
            dXdT_gas_value[i] = 0;
            dXdT_liq_value[i] = -(dvidT_value()[i] - X_liq[i] * dvfdT_value) / (1 - vf);
            dXdP_gas_value[i] = 0;
            dXdP_liq_value[i] = -(dvidP_value()[i] - X_liq[i] * dvfdP_value) / (1 - vf);
        }
    }
    else if (vf >= 1 - 1e-6)
    {
        for (int i = 0; i < nsp; i++)
        {
            dXdT_gas_value[i] = (dvidT_value()[i] - X_gas[i] * dvfdT_value) / vf;
            dXdT_liq_value[i] = 0;
            dXdP_gas_value[i] = (dvidP_value()[i] - X_gas[i] * dvfdP_value) / vf;
            dXdP_liq_value[i] = 0;
        }
    }
    else
    {
        for (int i = 0; i < nsp; i++)
        {
            dXdT_gas_value[i] = (dvidT_value()[i] - X_gas[i] * dvfdT_value) / vf;
            dXdT_liq_value[i] = -(dvidT_value()[i] - X_liq[i] * dvfdT_value) / (1 - vf);
            dXdP_gas_value[i] = (dvidP_value()[i] - X_gas[i] * dvfdP_value) / vf;
            dXdP_liq_value[i] = -(dvidP_value()[i] - X_liq[i] * dvfdP_value) / (1 - vf);
        }
    }

    List<scalarList> dXdXi_gas_value(nsp);
    List<scalarList> dXdXi_liq_value(nsp);

    if (vf <= 1e-6)
    {
        for (int i = 0; i < nsp; i++)
        {
            dXdXi_gas_value[i].resize(nsp);
            dXdXi_liq_value[i].resize(nsp);
            for (int j = 0; j < nsp; j++)
            {
                dXdXi_gas_value[i][j] = 0;
                if (i == j)
                {
                    dXdXi_liq_value[i][j] = 1 / (1 - vf) * (1 - dvidXi_value[i]()[j]) - X_liq[j] / (1 - vf) * (1 - dvfdXi_value[i] - vf);
                }
                else
                {
                    dXdXi_liq_value[i][j] = 1 / (1 - vf) * (-dvidXi_value[i]()[j]) - X_liq[j] / (1 - vf) * (1 - dvfdXi_value[i] - vf);
                }
            }
        }
    }
    else if (vf >= 1 - 1e-6)
    {
        for (int i = 0; i < nsp; i++)
        {
            dXdXi_gas_value[i].resize(nsp);
            dXdXi_liq_value[i].resize(nsp);
            for (int j = 0; j < nsp; j++)
            {
                dXdXi_gas_value[i][j] = (dvidXi_value[i]()[j] - X_gas[j] * (dvfdXi_value[i] + vf)) / vf;
                dXdXi_liq_value[i][j] = 0;
            }
        }
    }
    else
    {
        for (int i = 0; i < nsp; i++)
        {
            dXdXi_gas_value[i].resize(nsp);
            dXdXi_liq_value[i].resize(nsp);
            for (int j = 0; j < nsp; j++)
            {
                dXdXi_gas_value[i][j] = (dvidXi_value[i]()[j] - X_gas[j] * (dvfdXi_value[i] + vf)) / vf;
                if (i == j)
                {
                    dXdXi_liq_value[i][j] = 1 / (1 - vf) * (1 - dvidXi_value[i]()[j]) - X_liq[j] / (1 - vf) * (1 - dvfdXi_value[i] - vf);
                }
                else
                {
                    dXdXi_liq_value[i][j] = 1 / (1 - vf) * (-dvidXi_value[i]()[j]) - X_liq[j] / (1 - vf) * (1 - dvfdXi_value[i] - vf);
                }
            }
        }
    }

    for (int i = 0; i < nsp; i++)
    {
        dHsdT_gas += dHsdxi_gas_value[i] * dXdT_gas_value[i];
        dHsdT_liq += dHsdxi_liq_value[i] * dXdT_liq_value[i];
        dHsdP_gas += dHsdxi_gas_value[i] * dXdP_gas_value[i];
        dHsdP_liq += dHsdxi_liq_value[i] * dXdP_liq_value[i];
        dHsdXi_gas[i] = 0;
        dHsdXi_liq[i] = 0;
        for (int j = 0; j < nsp; j++)
        {
            dHsdXi_gas[i] += dHsdxi_gas_value[j] * dXdXi_gas_value[i][j];
            dHsdXi_liq[i] += dHsdxi_liq_value[j] * dXdXi_liq_value[i][j];
        }
    }

    scalar mw_gas = this->W(X_gas); //kg/mol
    scalar mw_mixture = this->W();
    scalar ygas = vf * mw_gas / mw_mixture;
    scalar dmw_gasdT = 0;
    scalar dmw_gasdP = 0;
    for (int i = 0; i < nsp; i++)
    {
        dmw_gasdT += dXdT_gas_value[i] * This[i].W() * 1.0e-03;
        dmw_gasdP += dXdP_gas_value[i] * This[i].W() * 1.0e-03;
    }
    scalarList dmw_gasdXi(nsp);
    for (int i = 0; i < nsp; i++)
    {
        dmw_gasdXi[i] = 0;
        for (int j = 0; j < nsp; j++)
        {
            dmw_gasdXi[i] += dXdXi_gas_value[i][j] * This[j].W() * 1.0e-03;
        }
    }
    scalar dygasdT = (dvfdT_value * mw_gas + vf * dmw_gasdT) / mw_mixture;
    scalar dygasdP = (dvfdP_value * mw_gas + vf * dmw_gasdP) / mw_mixture;
    scalarList dygasdXi(nsp);
    for (int i = 0; i < nsp; i++)
    {
        dygasdXi[i] = ((dvfdXi_value[i] + vf) * mw_gas + vf * dmw_gasdXi[i]) / mw_mixture - ygas / mw_mixture * This[i].W() * 1.0e-03;
    }

    scalar Hs_gas = this->ThermoMixture::Hs(p, T, X_gas, 1);
    scalar Hs_liq = this->ThermoMixture::Hs(p, T, X_liq, 0);
    scalar Cp_Hs_value = Hs_gas * dygasdT + dHsdT_gas * ygas - Hs_liq * dygasdT + dHsdT_liq * (1.0 - ygas);
    scalar dHsdP_value = Hs_gas * dygasdP + dHsdP_gas * ygas - Hs_liq * dygasdP + dHsdP_liq * (1.0 - ygas);
    scalarList dHsdXi_value(nsp);
    for (int i = 0; i < nsp; i++)
    {
        dHsdXi_value[i] = Hs_gas * dygasdXi[i] + dHsdXi_gas[i] * ygas - Hs_liq * dygasdXi[i] + dHsdXi_liq[i] * (1.0 - ygas);
    }

    scalar ZMixture = Z(p, T, sol);
    scalar ZGas = this->ThermoMixture::Z(p, T, X_gas, 1);
    scalar ZLiq = this->ThermoMixture::Z(p, T, X_liq, 0);
    scalar dZgasdT = this->ThermoMixture::dZdT(p, T, X_gas, 1);
    scalar dZliqdT = this->ThermoMixture::dZdT(p, T, X_liq, 0);
    scalar dZgasdP = this->ThermoMixture::dZdP(p, T, X_gas, 1);
    scalar dZliqdP = this->ThermoMixture::dZdP(p, T, X_liq, 0);
    scalarList dZgasdXi(nsp), dZliqdXi(nsp);
    scalarList dZdxi_gas_value(nsp), dZdxi_liq_value(nsp);
    for (int i = 0; i < nsp; i++)
    {
        dZdxi_gas_value[i] = this->ThermoMixture::dZdxi(p, T, i, X_gas, 1);
        dZdxi_liq_value[i] = this->ThermoMixture::dZdxi(p, T, i, X_liq, 0);
    }
    for (int i = 0; i < nsp; i++)
    {
        dZgasdT += dZdxi_gas_value[i] * dXdT_gas_value[i];
        dZliqdT += dZdxi_liq_value[i] * dXdT_liq_value[i];
        dZgasdP += dZdxi_gas_value[i] * dXdP_gas_value[i];
        dZliqdP += dZdxi_liq_value[i] * dXdP_liq_value[i];
        dZgasdXi[i] = 0;
        dZliqdXi[i] = 0;
        for (int j = 0; j < nsp; j++)
        {
            dZgasdXi[i] += dZdxi_gas_value[j] * dXdXi_gas_value[i][j];
            dZliqdXi[i] += dZdxi_liq_value[j] * dXdXi_liq_value[i][j];
        }
    }
    temp = (vf * ZGas + (1 - vf) * ZLiq);
    scalar dZMixturedT = dvfdT_value * (ZGas - ZLiq) + vf * dZgasdT + (1 - vf) * dZliqdT;
    scalar dZMixturedP = dvfdP_value * (ZGas - ZLiq) + vf * dZgasdP + (1 - vf) * dZliqdP;
    scalarList dZMixturedXi(nsp);

    for (int i = 0; i < nsp; i++)
    {
        dZMixturedXi[i] = dvfdXi_value[i] * (ZGas - ZLiq) + vf * dZgasdXi[i] + (1 - vf) * dZliqdXi[i];
    }

    scalar drhodT_value = -p * mw_mixture / (ZMixture * ZMixture * RR * 1.0e-03 * T * T) * (dZMixturedT * T + ZMixture);
    scalar drhodP_value = mw_mixture / (ZMixture * RR * 1.0e-03 * T) - p * mw_mixture / (ZMixture * ZMixture * RR * 1.0e-03 * T) * dZMixturedP;
    scalarList drhodXi_value(nsp);
    for (int i = 0; i < nsp; i++)
    {
        drhodXi_value[i] = p * (This[i].W() * 1.0e-03 - mw_mixture) / (ZMixture * RR * 1.0e-03 * T) - p * mw_mixture / (ZMixture * ZMixture * RR * 1.0e-03 * T) * dZMixturedXi[i];
    }

    scalar drhodT_G_value = -p * mw_gas / (ZGas * ZGas * RR * 1.0e-03 * T * T) * (dZgasdT * T + ZGas) + p / (ZGas * RR * 1.0e-03 * T) * dmw_gasdT;
    scalar drhodP_G_value = mw_gas / (ZGas * RR * 1.0e-03 * T) - p * mw_gas / (ZGas * ZGas * RR * 1.0e-03 * T) * dZgasdP + p / (ZGas * RR * 1.0e-03 * T) * dmw_gasdP;

    scalarList drhodXi_G_value(nsp);
    for (int i = 0; i < nsp; i++)
    {
        //drhodXi_G_value[i] = p * (This[i].W() * 1.0e-03 - mw_gas) / (ZGas * RR * 1.0e-03 * T) - p * mw_gas / (ZGas * ZGas * RR * 1.0e-03 * T) * dZgasdXi[i] + p / (ZGas * RR * 1.0e-03 * T) * dmw_gasdXi[i];
        drhodXi_G_value[i] = -p * mw_gas / (ZGas * ZGas * RR * 1.0e-03 * T) * dZgasdXi[i] + p / (ZGas * RR * 1.0e-03 * T) * dmw_gasdXi[i];
    }

    scalarList dYdT_G_value(nsp), dYdP_G_value(nsp);
    for (int i = 0; i < nsp; i++)
    {
        dYdT_G_value[i] = dXdT_gas_value[i] * This[i].W() * 1.0e-03 / mw_gas - sol.X_gas()[i] * This[i].W() * 1.0e-03 / (mw_gas * mw_gas) * dmw_gasdT;
        dYdP_G_value[i] = dXdP_gas_value[i] * This[i].W() * 1.0e-03 / mw_gas - sol.X_gas()[i] * This[i].W() * 1.0e-03 / (mw_gas * mw_gas) * dmw_gasdP;
    }

    List<scalarList> dYdXi_G_value(nsp);
    for (int i = 0; i < nsp; i++)
    {
        dYdXi_G_value[i].resize(nsp);
        for (int j = 0; j < nsp; j++)
        {
            dYdXi_G_value[i][j] = dXdXi_gas_value[i][j] * This[j].W() * 1.0e-03 / mw_gas - sol.X_gas()[j] * This[j].W() * 1.0e-03 / (mw_gas * mw_gas) * dmw_gasdXi[i];
        }
    }

    for (int i = 0; i < nsp; i++)
    {
        grad[i][0] = -dHsdXi_value[i] / Cp_Hs_value;
        grad[i][1] = -drhodT_value * dHsdXi_value[i] / Cp_Hs_value + drhodXi_value[i];
        grad[i][2] = -dvfdT_value * dHsdXi_value[i] / Cp_Hs_value + dvfdXi_value[i];
        grad[i][3] = 0;
        grad[i][4] = -drhodT_G_value * dHsdXi_value[i] / Cp_Hs_value + drhodXi_G_value[i];
        for (int j = 0; j < nsp; j++)
        {
            grad[i][5 + j] = (dXdT_gas_value[j] * This[j].W() * 1.0e-03 / mw_gas - sol.X_gas()[j] * This[j].W() * 1.0e-03 / (mw_gas * mw_gas) * dmw_gasdT) * dHsdXi_value[i] / Cp_Hs_value +
                             dXdXi_gas_value[i][j] * This[j].W() * 1.0e-03 / mw_gas - sol.X_gas()[j] * This[j].W() * 1.0e-03 / (mw_gas * mw_gas) * dmw_gasdXi[i];
        }
    }
    grad[nsp][0] = 1 / Cp_Hs_value;
    grad[nsp][1] = drhodT_value / Cp_Hs_value;
    grad[nsp][2] = dvfdT_value / Cp_Hs_value;
    grad[nsp][3] = 0;
    grad[nsp][4] = drhodT_G_value / Cp_Hs_value;
    for (int j = 0; j < nsp; j++)
    {
        grad[nsp][5 + j] = (dXdT_gas_value[j] * This[j].W() * 1.0e-03 / mw_gas - sol.X_gas()[j] * This[j].W() * 1.0e-03 / (mw_gas * mw_gas) * dmw_gasdT) / Cp_Hs_value;
    }

    grad[nsp + 1][0] = -dHsdP_value / Cp_Hs_value;
    grad[nsp + 1][1] = -drhodT_value * dHsdP_value / Cp_Hs_value + drhodP_value;
    grad[nsp + 1][2] = -dvfdT_value * dHsdP_value / Cp_Hs_value + dvfdP_value;
    grad[nsp + 1][3] = 0;
    grad[nsp + 1][4] = -drhodT_G_value * dHsdP_value / Cp_Hs_value + drhodP_G_value;
    for (int j = 0; j < nsp; j++)
    {
        grad[nsp + 1][5 + j] = -(dXdT_gas_value[j] * This[j].W() * 1.0e-03 / mw_gas - sol.X_gas()[j] * This[j].W() * 1.0e-03 / (mw_gas * mw_gas) * dmw_gasdT) * dHsdP_value / Cp_Hs_value +
                               dXdP_gas_value[j] * This[j].W() * 1.0e-03 / mw_gas - sol.X_gas()[j] * This[j].W() * 1.0e-03 / (mw_gas * mw_gas) * dmw_gasdP;
    }
    return autoPtr<scalarRectangularMatrix>(&grad);
}