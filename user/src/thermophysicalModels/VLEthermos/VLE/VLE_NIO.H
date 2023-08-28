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

//#define HY_DEBUG

// * * * * * * * * * * * * * * * * Basic properties  * * * * * * * * * * * * * *

template <class ThermoMixture>
inline void Foam::VLE<ThermoMixture>::TPn_flash_Matheis_NIO(
    scalar p, //input
    scalar T, //input
    const scalarList &Xinp,
    Foam::scalarList &K_init,
    bool &success,
    VLE_data &ret_solution,
    scalar &dg) const
{
    const int itermax_outer = 5000;
    const int itermax_inner = 100;

    const double tol_inner = 1e-9;
    const double tol_outer = 1e-8;
    success = false;

    scalar &Z_l = ret_solution.Liq.Z;
    scalar &Z_g = ret_solution.Gas.Z;

    scalar dgl, dgv;

    scalarList &K = ret_solution.K;
    K = K_init;
    scalar &vf = ret_solution.vaporFraction;

    //autoPtr<scalarList> Ln_fugcoef_liq(new scalarList(this->N_));
    //autoPtr<scalarList> Ln_fugcoef_gas(new scalarList(this->N_));
    scalarList Ln_fugcoef_liq(this->N_);
    scalarList Ln_fugcoef_gas(this->N_);
    scalarList F(this->N_);
    scalarList F_old(this->N_);
    scalarList B_(this->N_);
    scalarList Tsqrt_AbyP_(this->N_);
    scalar sqrtT = sqrt(T);
    scalarList Ktmp(this->N_);
    scalarList R_(this->N_);
    scalarList lnK(this->N_);
    scalar relax = 1.0;
    for (int i = 0; i < this->N_; i++)
    {
        lnK[i] = log(K[i]);
    }
    forAll(Xinp, i)
    {
        B_[i] = (*this)[i].B(p, T);
        Tsqrt_AbyP_[i] = (*this)[i].Tsqrt_AbyP(sqrtT);
    }

    int iter_outer = 0;
    double eps_outer = 0;
    while (iter_outer <= itermax_outer)
    {

        /*Info << '[' << K[0];
        for (int i = 1; i < this->N_; i++)
        {
            Info << ',' << K[i];
        }
        Info << "]," << endl;
        */
        //Info << '[' << K[0] << ',' << K[1] << "]," << endl;
        scalar sum1 = 0, sum2 = 0;
        for (int i = 0; i < this->N_; i++)
        {
            sum1 += Xinp[i] * K[i];
            sum2 += Xinp[i] / K[i];
        }
        /*if (iter_outer == 0 && (sum1 < 1 || sum2 < 1))
        {
            K[0] = 1 / Xinp[0];
            K[1] = Xinp[1];
            continue;
        }*/
        if (sum1 >= 1.0 && sum2 >= 1)
        {
            double eps_inner = 1.0;
            int iter_inner = 1;
            double vf_min = 0, vf_max = 1;

            for (int i = 0; i < this->N_; i++)
            {
                if (K[i] > 1)
                    vf_min = max(0.0, (Xinp[i] * K[i] - 1) / (K[i] - 1));
                else
                    vf_max = min(1.0, (1 - Xinp[i]) / (1 - K[i]));
            }
            vf = 0.5 * (vf_min + vf_max);

            double RR = 0, dRR = 0, dvf = 0, vf_n = 0;
            while (eps_inner >= tol_inner && iter_inner <= itermax_inner)
            {
                RR = 0;
                dRR = 0;
                for (int i = 0; i < this->N_; i++)
                {
                    RR += Xinp[i] * (K[i] - 1) / (1 + vf * (K[i] - 1));
                    dRR -= Xinp[i] * (K[i] - 1) * (K[i] - 1) / ((1 + vf * (K[i] - 1)) * (1 + vf * (K[i] - 1)));
                }
                if (RR > 0)
                    vf_min = vf;
                else
                    vf_max = vf;

                dvf = -RR / dRR;
                vf_n = vf + dvf;

                if (vf_n < vf_min || vf_n > vf_max)
                    vf_n = 0.5 * (vf_min + vf_max);

                eps_inner = fabs((vf_n - vf) / vf);
                vf = vf_n;

                iter_inner++;
            }

            for (int i = 0; i < this->N_; i++)
            {
                ret_solution.Liq.X[i] = Xinp[i] / (1 + vf * (K[i] - 1));
                ret_solution.Gas.X[i] = K[i] * ret_solution.Liq.X[i];
            }
        }
        else if (sum2 < 1)
        {
            vf = 1;
            double sum_ = 0;
            for (int i = 0; i < this->N_; i++)
            {
                //K[i] = sqr(K[i]);
                ret_solution.Liq.X[i] = Xinp[i] / K[i];
                ret_solution.Gas.X[i] = Xinp[i];
                sum_ += ret_solution.Liq.X[i];
            }
            for (int i = 0; i < this->N_; i++)
            {
                ret_solution.Liq.X[i] /= sum_;
            }
            Z_g = this->ThermoMixture::Z_gibbs_mix(p, T, ret_solution.Gas.X);
            //Z_g = this->ThermoMixture::Z(p, T, ret_solution.Gas.X, 1);
            //Z_l = this->ThermoMixture::Z(p, T, ret_solution.Liq.X, 0);
            Z_l = this->ThermoMixture::Z_gibbs_mix(p, T, ret_solution.Liq.X);
            break;
        }
        else if (sum1 <= 1)
        {
            vf = 0;
            double sum_ = 0;
            for (int i = 0; i < this->N_; i++)
            {
                //K[i] = sqr(K[i]);
                ret_solution.Liq.X[i] = Xinp[i];
                ret_solution.Gas.X[i] = Xinp[i] * K[i];
                sum_ += ret_solution.Gas.X[i];
            }
            for (int i = 0; i < this->N_; i++)
            {
                ret_solution.Gas.X[i] /= sum_;
            }
            Z_g = this->ThermoMixture::Z_gibbs_mix(p, T, ret_solution.Gas.X);
            //Z_g = this->ThermoMixture::Z(p, T, ret_solution.Gas.X, 1);
            //Z_l = this->ThermoMixture::Z(p, T, ret_solution.Liq.X, 0);
            Z_l = this->ThermoMixture::Z_gibbs_mix(p, T, ret_solution.Liq.X);
            break;
        }
        else
            FatalErrorInFunction << "TPn_flash_Matheis" << exit(FatalError);

        if (iter_outer > 0 && (vf > 1 - 1e-8 || vf < 1e-8))
        {
            if (vf > 1 - 1e-8)
                vf = 1.0;
            else
                vf = 0.0;
            for (int i = 0; i < this->N_; i++)
            {
                ret_solution.Liq.X[i] = Xinp[i];
                ret_solution.Gas.X[i] = Xinp[i];
            }
            Z_g = this->ThermoMixture::Z_gibbs_mix(p, T, ret_solution.Gas.X);
            Z_l = this->ThermoMixture::Z_gibbs_mix(p, T, ret_solution.Liq.X);
            //Z_g = this->ThermoMixture::Z(p, T, ret_solution.Gas.X, 1);
            //Z_l = this->ThermoMixture::Z(p, T, ret_solution.Liq.X, 0);
            success = true;
            break;
        }

        this->Ln_fugacityCoefficient_opt(p, T, ret_solution.Liq.X, Ln_fugcoef_liq, B_, Tsqrt_AbyP_, Z_l, dgl, 0);

        this->Ln_fugacityCoefficient_opt(p, T, ret_solution.Gas.X, Ln_fugcoef_gas, B_, Tsqrt_AbyP_, Z_g, dgv, 1);
        dg = dgl * (1 - vf) + dgv * vf;

        eps_outer = 0;
        scalarList R_(this->N_);

        for (int i = 0; i < this->N_; i++)
        {
            //scalar logK=log(K[i]);
            //F_old[i] = F[i];
            F[i] = Ln_fugcoef_liq[i] - Ln_fugcoef_gas[i] - lnK[i]; //log(K[i]);
            if (relax < 1)
                R_[i] = exp(F[i] * relax * 2);
            else
                R_[i] = exp(F[i] * relax);

            eps_outer += F[i] * F[i];
        }
        eps_outer = sqrt(eps_outer);
        //if (relax < 1 || (iter_outer > 300 && relax < 4) || (iter_outer > 600 && relax < 4) || (iter_outer > 1000 && relax < 4))
        //    relax *= 4;
        //else
        relax *= 2;
        if (eps_outer < tol_outer)
        {
            success = true;
            break;
        }
        else
        {

            //scalar relax = 2;
            scalar sum1 = 0, sum2 = 0;

            //relax *= 2;

            do
            {
                relax /= 2;
                for (int i = 0; i < this->N_; i++)
                {
                    //Ktmp[i] = exp(Ln_fugcoef_liq()[i] - Ln_fugcoef_gas()[i]);
                    //Ktmp[i] = K[i] * exp(F[i] * relax);
                    Ktmp[i] = K[i] * R_[i];
                }
                sum1 = 0;
                sum2 = 0;
                for (int i = 0; i < this->N_; i++)
                {
                    sum1 += Xinp[i] * Ktmp[i];
                    sum2 += Xinp[i] / Ktmp[i];
                    R_[i] = sqrt(R_[i]);
                }
            } while ((sum1 < 1.0 || sum2 < 1) && relax > 1e-10);
            for (int i = 0; i < this->N_; i++)
            {
                //Ktmp[i] = exp(Ln_fugcoef_liq()[i] - Ln_fugcoef_gas()[i]);
                K[i] = Ktmp[i];
                lnK[i] += relax * F[i];
                //K[i] = exp(Ln_fugcoef_liq()[i] - Ln_fugcoef_gas()[i]);
            }
        }
        iter_outer++;
    }
    //if (iter_outer > itermax_outer||)
    //    success = true;

    return;
}
template <class ThermoMixture>
inline void Foam::VLE<ThermoMixture>::TPn_flash_New_TPD_Tudisco_NIO // TPn flash with new IO
    (
        scalar p,               //input
        scalar T,               //input
        const scalarList &Xinp, //input
        VLE_data &ret_solution  //output
        ) const
{
#ifdef HY_DEBUG
    Info << setprecision(50) << "p=" << p << ",T=" << T << ",X=" << Xinp << endl;
    p = 47510858.779683121;
    T = 567.0820033742425;
    this->X_[0] = 0.15081333856095205;
    this->X_[1] = 0.024826764817530642;
    this->X_[2] = 0.82435989662151732;
#endif
    if (noVLE) //TODO: update all values
    {
        forAll(Xinp, i)
        {
            ret_solution.Gas.X[i] = Xinp[i];
            ret_solution.Liq.X[i] = Xinp[i];
            ret_solution.K[i] = 1;
        }
        ret_solution.vaporFraction = 1;
        NotImplemented;
        return;
    }

    if (Xinp.size() == 1) //TODO: update all values
    {
        autoPtr<scalarList> fugcoef_liq(this->Ln_fugacityCoefficient(p, T, Xinp, 0).ptr());
        autoPtr<scalarList> fugcoef_gas(this->Ln_fugacityCoefficient(p, T, Xinp, 1).ptr());
        ret_solution.K[0] = 1;
        ret_solution.Gas.X[0] = 1;
        ret_solution.Liq.X[0] = 1;
        if (fugcoef_liq()[0] > fugcoef_gas()[0])
        {
            ret_solution.vaporFraction = 1;
        }
        else
        {
            ret_solution.vaporFraction = 0;
        }
        NotImplemented;
        return;
    }

    bool success = false;
    bool stable = false, isVapor = false;
    scalar gmin = DBL_MAX;
    const scalar delta = 1e-13;
    const scalar be_negative = -1.E-9;
    autoPtr<List<scalarList>> p_K_init;
    autoPtr<solution> sol_min;
    autoPtr<scalarList> p_TPD_star;
    scalarList K_init(Xinp.size());
    scalar maxK = -1, minK = 1e10;
    int maxI = 0, minI = 0;
    if (inputK == false)
    {

        //Info << "Wil=";
        for (int i = 0; i < Xinp.size(); i++)
        {

            K_init[i] = (*this)[i].Pc_ / p * exp(5.373 * (1.0 + (*this)[i].omega_) * (1.0 - (*this)[i].Tc_ / T));
            //Info << K_init[i] << " ,";
            if (K_init[i] > maxK)
            {
                maxK = K_init[i];
                maxI = i;
            }
            if (K_init[i] < minK)
            {
                minK = K_init[i];
                minI = i;
            }

            K_init[i] *= max((*this)[i].Pc_ / p, 1);
        }
        //Info << endl;
        scalar rat, mid;
        if (maxK > 1e4)
        {
            rat = sqrt(maxK / minK);
            mid = sqrt(sqrt(maxK * minK));
        }
        else
        {
            rat = sqrt(maxK / minK) * 10;
            mid = sqrt(sqrt(maxK * minK));
        }
        //scalar rat = sqrt(maxK / minK), mid = sqrt(sqrt(maxK * minK));
        //scalar rat = maxK / minK, mid = sqrt(maxK * minK);
        for (int i = 0; i < Xinp.size(); i++)
        {
            if (K_init[i] > mid)
                K_init[i] *= rat;
            else
                K_init[i] /= rat;
            K_init[i] /= mid;
        }
    }
    if (inputK == true)
    {
        for (int i = 0; i < Xinp.size(); i++)
        {
            K_init[i] = Kinit[i];
        }
    }
    scalar ggas1, gliq1, g1 = -1e10; //, Zliq1, Zgas1;
    TPn_flash_Matheis_NIO(p, T, Xinp, K_init, success, ret_solution, g1);
    //autoPtr<solution> sol = TPn_flash_Matheis(p, T, Xinp, stable, isVapor, K_init, success, Zliq1, Zgas1, g1);
    if (success)
    {
        ggas1 = this->Gideal_Mole(p, T, ret_solution.Gas.X);
        gliq1 = this->Gideal_Mole(p, T, ret_solution.Liq.X);
        g1 += ret_solution.vaporFraction * ggas1 + (1 - ret_solution.vaporFraction) * gliq1;
        /*
        ggas1 = this->Gideal_Mole(p, T, sol().X_gas());
        gliq1 = this->Gideal_Mole(p, T, sol().X_liq());
        g1 += sol().vaporfra * ggas1 + (1 - sol().vaporfra) * gliq1;
*/

        if (ret_solution.Liq.Z > 0.87 && ret_solution.Gas.Z > 0.95 && ret_solution.vaporFraction < 0.99999 && ret_solution.vaporFraction > 0.00001)
        {
            g1 += 10000;
        }

        /*
        if (Zliq1 > 0.87 && Zgas1 > 0.95 && sol().vaporfra < 0.99999 && sol().vaporfra > 0.00001)
        {
            g1 += 10000;
        }
        */

        if (ret_solution.Liq.Z > ret_solution.Gas.Z)
        {
            scalarList Xliq_temp = ret_solution.Liq.X;
            scalarList Xgas_temp = ret_solution.Gas.X;
            scalar Liq_Z = ret_solution.Liq.Z;
            scalar Gas_Z = ret_solution.Gas.Z;
            ret_solution.Liq.Z = Gas_Z;
            ret_solution.Gas.Z = Liq_Z;
            scalar alpha = 1.0 - ret_solution.vaporFraction;
            forAll(Xinp, i)
            {
                ret_solution.Liq.X[i] = Xgas_temp[i];
                ret_solution.Gas.X[i] = Xliq_temp[i];
                ret_solution.K[i] = 1.0 / ret_solution.K[i];
            }
            ret_solution.vaporFraction = alpha;
        }
    }

    if (Xinp.size() > 2)
    {
        for (int i = 0; i < Xinp.size(); i++)
        {
            if (i == maxI || i == minI)
                continue;
            K_init[i] = 1 / K_init[i];
        }

        bool success2 = false;
        scalar ggas2, gliq2, g2 = -1e10; //, Zgas2, Zliq2;
        VLE_data solution2(Xinp.size());
        TPn_flash_Matheis_NIO(p, T, Xinp, K_init, success2, solution2, g2);
        //autoPtr<solution> sol2 = TPn_flash_Matheis(p, T, Xinp, stable, isVapor, K_init, success2, Zliq2, Zgas2, g2);
        if (success2)
        {
            ggas2 = this->Gideal_Mole(p, T, solution2.Gas.X);
            gliq2 = this->Gideal_Mole(p, T, solution2.Liq.X);
            g2 += solution2.vaporFraction * ggas2 + (1 - solution2.vaporFraction) * gliq2;

            if (solution2.Liq.Z > 0.87 && solution2.Gas.Z > 0.95 && solution2.vaporFraction < 0.99999 && solution2.vaporFraction > 0.00001)
            {
                g2 += 10000;
            }
            //if (sol2().vaporfra < 0.000001 || sol2().vaporfra > 0.99999)
            //    g2 -= 0.1;
            //Info << "G=" << g2 << endl;

            if (solution2.Liq.Z > solution2.Gas.Z)
            {
                scalarList Xliq_temp = solution2.Liq.X;
                scalarList Xgas_temp = solution2.Gas.X;
                scalar Liq_Z = solution2.Liq.Z;
                scalar Gas_Z = solution2.Gas.Z;
                solution2.Liq.Z = Gas_Z;
                solution2.Gas.Z = Liq_Z;
                scalar alpha = 1.0 - solution2.vaporFraction;
                forAll(Xinp, i)
                {
                    solution2.Liq.X[i] = Xgas_temp[i];
                    solution2.Gas.X[i] = Xliq_temp[i];
                    solution2.K[i] = 1.0 / solution2.K[i];
                }
                solution2.vaporFraction = alpha;
            }
        }
        if (success2 && success)
        {
            if (g1 < g2)
            {
                return;
            }
            else
            {
                forAll(Xinp, i)
                {
                    ret_solution.Gas.X[i] = solution2.Gas.X[i];
                    ret_solution.Liq.X[i] = solution2.Liq.X[i];
                    ret_solution.K[i] = solution2.K[i];
                }
                ret_solution.vaporFraction = solution2.vaporFraction;
                ret_solution.Gas.Z = solution2.Gas.Z;
                ret_solution.Liq.Z = solution2.Liq.Z;

                /*                 forAll(Xinp, i)
                {
                    ret_solution.Gas.X[i] = sol2().X_gas()[i];
                    ret_solution.Liq.X[i] = sol2().X_liq()[i];
                    ret_solution.K[i] = sol2().equalconstant()[i];
                }
                ret_solution.vaporFraction = sol2().vaporfra;
                ret_solution.Gas.Z = this->ThermoMixture::Z(p, T, ret_solution.Gas.X, 1);
                ret_solution.Liq.Z = this->ThermoMixture::Z(p, T, ret_solution.Liq.X, 0); */
                return;
            }
        }
        else if (success)
        {
            /*             forAll(Xinp, i)
            {
                ret_solution.Gas.X[i] = sol().X_gas()[i];
                ret_solution.Liq.X[i] = sol().X_liq()[i];
                ret_solution.K[i] = sol().equalconstant()[i];
            }
            ret_solution.vaporFraction = sol().vaporfra; */
            return;
        }
        else if (success2)
        {

            forAll(Xinp, i)
            {
                ret_solution.Gas.X[i] = solution2.Gas.X[i];
                ret_solution.Liq.X[i] = solution2.Liq.X[i];
                ret_solution.K[i] = solution2.K[i];
            }
            ret_solution.vaporFraction = solution2.vaporFraction;
            ret_solution.Gas.Z = solution2.Gas.Z;
            ret_solution.Liq.Z = solution2.Liq.Z;
            /*             forAll(Xinp, i)
            {
                ret_solution.Gas.X[i] = sol2().X_gas()[i];
                ret_solution.Liq.X[i] = sol2().X_liq()[i];
                ret_solution.K[i] = sol2().equalconstant()[i];
            }
            ret_solution.vaporFraction = sol2().vaporfra;
            ret_solution.Gas.Z = this->ThermoMixture::Z(p, T, ret_solution.Gas.X, 1);
            ret_solution.Liq.Z = this->ThermoMixture::Z(p, T, ret_solution.Liq.X, 0); */
            return;
        }
    }
    else
    {
        if (success)
        {
            /*             forAll(Xinp, i)
            {
                ret_solution.Gas.X[i] = sol().X_gas()[i];
                ret_solution.Liq.X[i] = sol().X_liq()[i];
                ret_solution.K[i] = sol().equalconstant()[i];
            }
            ret_solution.vaporFraction = sol().vaporfra; */
            return;
        }
    }

    /*if (sol().vaporfra < 1e-8||sol().vaporfra > 0.99999)
    {
        scalarList *liq = sol().X_gas.ptr();
        scalarList *gas = sol().X_liq.ptr();
        sol().X_gas.reset(gas);
        sol().X_liq.reset(liq);
        sol().vaporfra = 1;
    }*/
    /*
    if(sol().vaporfra<0.99&&sol().vaporfra>0.01)
    {
        //Info<<sol().vaporfra<<endl;
        //FatalErrorInFunction << "K=!! "<< exit(FatalError);
        FatalErrorInFunction << "K= "
                << sol().equalconstant()  << " \n vf=" << sol().vaporfra << exit(FatalError);
    }
    */

    //return sol;
    if (success)
    {
        /*         forAll(Xinp, i)
        {
            ret_solution.Gas.X[i] = sol().X_gas()[i];
            ret_solution.Liq.X[i] = sol().X_liq()[i];
            ret_solution.K[i] = sol().equalconstant()[i];
        }
        ret_solution.vaporFraction = sol().vaporfra; */
        return;
    }
    //scalar tttt;
    scalar tmp_g, tmp_Zgas, tmp_Zliq;
    tie(stable, p_K_init, isVapor, p_TPD_star) = solveTPD_BFGS_v2(p, T);
    labelList order;
    sortedOrder(p_TPD_star(), order);
    if (p_TPD_star()[order[0]] < be_negative)
    {
        autoPtr<solution> sol2 = TPn_flash_Matheis(p, T, Xinp, stable, isVapor, p_K_init()[order[0]], success, tmp_Zliq, tmp_Zgas, tmp_g);
        forAll(Xinp, i)
        {
            ret_solution.Gas.X[i] = sol2().X_gas()[i];
            ret_solution.Liq.X[i] = sol2().X_liq()[i];
            ret_solution.K[i] = sol2().equalconstant()[i];
        }
        ret_solution.vaporFraction = sol2().vaporfra;
        ret_solution.Gas.Z = this->ThermoMixture::Z(p, T, ret_solution.Gas.X, 1);
        ret_solution.Liq.Z = this->ThermoMixture::Z(p, T, ret_solution.Liq.X, 0);
        return;
        if (success)
        {
            forAll(Xinp, i)
            {
                ret_solution.Gas.X[i] = sol2().X_gas()[i];
                ret_solution.Liq.X[i] = sol2().X_liq()[i];
                ret_solution.K[i] = sol2().equalconstant()[i];
            }
            ret_solution.vaporFraction = sol2().vaporfra;
            ret_solution.Gas.Z = this->ThermoMixture::Z(p, T, ret_solution.Gas.X, 1);
            ret_solution.Liq.Z = this->ThermoMixture::Z(p, T, ret_solution.Liq.X, 0);
            return;
        }
    }
    else
    {
        /*         forAll(Xinp, i)
        {
            ret_solution.Gas.X[i] = sol().X_gas()[i];
            ret_solution.Liq.X[i] = sol().X_liq()[i];
            ret_solution.K[i] = sol().equalconstant()[i];
        }
        ret_solution.vaporFraction = sol().vaporfra; */

        return;
    }

    for (int i = 0; i < p_TPD_star->size(); i++)
    {

        if (i == 0 || (p_TPD_star()[order[i]] < be_negative && fabs(p_TPD_star()[order[i]] - p_TPD_star()[order[i - 1]]) > delta))
        {
            scalar g;
            autoPtr<solution> sol = TPn_flash_Matheis(p, T, Xinp, stable, isVapor, p_K_init()[order[i]], success, tmp_Zliq, tmp_Zgas, g);
            //scalar ggas = this->G_Mole(p, T, sol().X_gas());
            //scalar gliq = this->G_Mole(p, T, sol().X_liq());
            //scalar g = sol().vaporfra * ggas + (1 - sol().vaporfra) * gliq;

            if (g < gmin)
            {
                gmin = g;
                autoPtrReset(sol_min, sol);
                //sol_min.reset(sol);
            }
        }
    }
    forAll(Xinp, i)
    {
        ret_solution.Gas.X[i] = sol_min().X_gas()[i];
        ret_solution.Liq.X[i] = sol_min().X_liq()[i];
        ret_solution.K[i] = sol_min().equalconstant()[i];
    }
    ret_solution.vaporFraction = sol_min().vaporfra;
    ret_solution.Gas.Z = this->ThermoMixture::Z(p, T, ret_solution.Gas.X, 1);
    ret_solution.Liq.Z = this->ThermoMixture::Z(p, T, ret_solution.Liq.X, 0);

    return;
}
template <class ThermoMixture>
inline autoPtr<scalarList>
Foam::VLE<ThermoMixture>::dvidP(scalar p, scalar T, VLE_data &sol) const
{
    scalarSquareMatrix A(this->X_.size(), Zero);
    scalarList &b = *(new scalarList(this->X_.size()));
    scalar vf_ = sol.vaporFraction;
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
        (this->ddP_Ln_fugacityCoefficient(p, T, sol.Gas.X, 1)).ptr());

    for (int i = 0; i < this->X_.size(); i++)
        ddxi_Ln_fugcoef_G_[i].reset(
            (this->ddxi_Ln_fugacityCoefficient(p, T, i, sol.Gas.X, 1)).ptr());

    ddP_Ln_fugcoef_L_.reset(
        (this->ddP_Ln_fugacityCoefficient(p, T, sol.Liq.X, 0)).ptr());

    for (int i = 0; i < this->X_.size(); i++)
        ddxi_Ln_fugcoef_L_[i].reset(
            (this->ddxi_Ln_fugacityCoefficient(p, T, i, sol.Liq.X, 0)).ptr());
    for (int i = 0; i < this->X_.size(); i++)
    {
        for (int j = 0; j < this->X_.size(); j++)
        {
            A(i, j) = 1;
            if (i == j)
                A(i, j) -= this->X_[i] / (sol.Gas.X[i] * sol.Liq.X[i]);
            for (int k = 0; k < this->X_.size(); k++)
            {
                if (j == k)
                    A(i, j) +=
                        (1 - vf_) * ddxi_Ln_fugcoef_G_[k]()[i] * (sol.Gas.X[k] - 1);
                else
                    A(i, j) += (1 - vf_) * ddxi_Ln_fugcoef_G_[k]()[i] * sol.Gas.X[k];
                if (j == k)
                    A(i, j) += vf_ * ddxi_Ln_fugcoef_L_[k]()[i] * (sol.Liq.X[k] - 1);
                else
                    A(i, j) += vf_ * ddxi_Ln_fugcoef_L_[k]()[i] * sol.Liq.X[k];
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
Foam::VLE<ThermoMixture>::dvidT(scalar p, scalar T, VLE_data &sol) const
{
    scalarSquareMatrix A(this->X_.size(), Zero);
    scalarList &b = *(new scalarList(this->X_.size()));
    scalar vf_ = sol.vaporFraction;
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
        (this->ddT_Ln_fugacityCoefficient(p, T, sol.Gas.X, 1)).ptr());
    for (int i = 0; i < this->X_.size(); i++)
        ddxi_Ln_fugcoef_G_[i].reset(
            (this->ddxi_Ln_fugacityCoefficient(p, T, i, sol.Gas.X, 1)).ptr());
    ddT_Ln_fugcoef_L_.reset(
        (this->ddT_Ln_fugacityCoefficient(p, T, sol.Liq.X, 0)).ptr());
    for (int i = 0; i < this->X_.size(); i++)
        ddxi_Ln_fugcoef_L_[i].reset(
            (this->ddxi_Ln_fugacityCoefficient(p, T, i, sol.Liq.X, 0)).ptr());
    for (int i = 0; i < this->X_.size(); i++)
    {
        for (int j = 0; j < this->X_.size(); j++)
        {
            A(i, j) = 1;
            if (i == j)
                A(i, j) -= this->X_[i] / (sol.Gas.X[i] * sol.Liq.X[i]);
            for (int k = 0; k < this->X_.size(); k++)
            {
                if (j == k)
                    A(i, j) +=
                        (1 - vf_) * ddxi_Ln_fugcoef_G_[k]()[i] * (sol.Gas.X[k] - 1);
                else
                    A(i, j) += (1 - vf_) * ddxi_Ln_fugcoef_G_[k]()[i] * sol.Gas.X[k];
                if (j == k)
                    A(i, j) += vf_ * ddxi_Ln_fugcoef_L_[k]()[i] * (sol.Liq.X[k] - 1);
                else
                    A(i, j) += vf_ * ddxi_Ln_fugcoef_L_[k]()[i] * sol.Liq.X[k];
            }
        }
    }
    for (int i = 0; i < this->X_.size(); i++)
        b[i] = (1 - vf_) * vf_ * (ddT_Ln_fugcoef_G_()[i] - ddT_Ln_fugcoef_L_()[i]);
    LUsolve(A, b);
    return autoPtr<scalarList>(&b);
}

template <class ThermoMixture>
inline tuple<scalar, scalar, scalar, scalar> Foam::VLE<ThermoMixture>::dTPdErho_NIO(scalar p, scalar T, VLE_D1_data &ret_solution) const
{
    int nsp = this->X_.size();

    const VLE<ThermoMixture> &This = *this;

    scalar &vf = ret_solution.vaporFraction;
    scalarList &X_gas = ret_solution.Gas.X;
    scalarList &X_liq = ret_solution.Liq.X;

    scalar &ZGas = ret_solution.Gas.Z;
    scalar &ZLiq = ret_solution.Liq.Z;

    //scalar ZMixture = vf * ret_solution.Gas.Z + (1 - vf) * ret_solution.Liq.Z;

    scalar rho_value = rho_NIO(p, T, ret_solution); //p * mw_mixture / (ZMixture * RR * 1.0e-03 * T);

    autoPtr<scalarList> dvidP_value(this->dvidP(p, T, ret_solution));
    autoPtr<scalarList> dvidT_value(this->dvidT(p, T, ret_solution));

    scalar dvfdP_value = 0, dvfdT_value = 0;

    scalarList dHsdxi_gas_value(nsp), dHsdxi_liq_value(nsp);
    for (int i = 0; i < nsp; i++)
    {
        dvfdT_value += dvidT_value()[i];
        dvfdP_value += dvidP_value()[i];
        dHsdxi_gas_value[i] = this->dHsdxi(p, T, i, ZGas, X_gas, 1);
        dHsdxi_liq_value[i] = this->dHsdxi(p, T, i, ZLiq, X_liq, 0);
    }
    scalar dHsdT_gas = this->ThermoMixture::dHsdT(p, T, ZGas, X_gas, 1);
    scalar dHsdT_liq = this->ThermoMixture::dHsdT(p, T, ZLiq, X_liq, 0);
    scalar dHsdP_gas = this->ThermoMixture::dHsdP(p, T, ZGas, X_gas, 1);
    scalar dHsdP_liq = this->ThermoMixture::dHsdP(p, T, ZLiq, X_liq, 0);

    scalarList dXdT_gas_value(nsp);
    scalarList dXdT_liq_value(nsp);
    scalarList dXdP_gas_value(nsp);
    scalarList dXdP_liq_value(nsp);

    if (vf == 0)
    {
        for (int i = 0; i < nsp; i++)
        {
            dXdT_gas_value[i] = 0;
            dXdT_liq_value[i] = -(dvidT_value()[i] - X_liq[i] * dvfdT_value) / (1 - vf);
            dXdP_gas_value[i] = 0;
            dXdP_liq_value[i] = -(dvidP_value()[i] - X_liq[i] * dvfdP_value) / (1 - vf);
        }
    }
    else if (vf == 1)
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

    for (int i = 0; i < nsp; i++)
    {
        dHsdT_gas += dHsdxi_gas_value[i] * dXdT_gas_value[i];
        dHsdT_liq += dHsdxi_liq_value[i] * dXdT_liq_value[i];
        dHsdP_gas += dHsdxi_gas_value[i] * dXdP_gas_value[i];
        dHsdP_liq += dHsdxi_liq_value[i] * dXdP_liq_value[i];
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
    scalar dygasdT = (dvfdT_value * mw_gas + vf * dmw_gasdT) / mw_mixture;
    scalar dygasdP = (dvfdP_value * mw_gas + vf * dmw_gasdP) / mw_mixture;

    scalar Hs_gas = this->ThermoMixture::Hs(p, T, ZGas, X_gas, 1);
    scalar Hs_liq = this->ThermoMixture::Hs(p, T, ZLiq, X_liq, 0);
    scalar Cp_Hs_value = Hs_gas * dygasdT + dHsdT_gas * ygas - Hs_liq * dygasdT + dHsdT_liq * (1.0 - ygas);
    scalar dHsdP_value = Hs_gas * dygasdP + dHsdP_gas * ygas - Hs_liq * dygasdP + dHsdP_liq * (1.0 - ygas);

    scalar ZMixture = Z_NIO(p, T, ret_solution);
    //scalar ZGas = this->ThermoMixture::Z(p, T, X_gas, 1);
    //scalar ZLiq = this->ThermoMixture::Z(p, T, X_liq, 0);

    scalar dZgasdT = this->ThermoMixture::dZdT(p, T, ZGas, X_gas, 1);
    scalar dZliqdT = this->ThermoMixture::dZdT(p, T, ZLiq, X_liq, 0);
    scalar dZgasdP = this->ThermoMixture::dZdP(p, T, ZGas, X_gas, 1);
    scalar dZliqdP = this->ThermoMixture::dZdP(p, T, ZLiq, X_liq, 0);
    scalarList dZgasdXi(nsp), dZliqdXi(nsp);
    scalarList dZdxi_gas_value(nsp), dZdxi_liq_value(nsp);
    for (int i = 0; i < nsp; i++)
    {
        dZdxi_gas_value[i] = this->ThermoMixture::dZdxi(p, T, i, ZGas, X_gas, 1);
        dZdxi_liq_value[i] = this->ThermoMixture::dZdxi(p, T, i, ZLiq, X_liq, 0);
    }
    for (int i = 0; i < nsp; i++)
    {
        dZgasdT += dZdxi_gas_value[i] * dXdT_gas_value[i];
        dZliqdT += dZdxi_liq_value[i] * dXdT_liq_value[i];
        dZgasdP += dZdxi_gas_value[i] * dXdP_gas_value[i];
        dZliqdP += dZdxi_liq_value[i] * dXdP_liq_value[i];
    }

    scalar dZMixturedT = dvfdT_value * (ZGas - ZLiq) + vf * dZgasdT + (1 - vf) * dZliqdT;
    scalar dZMixturedP = dvfdP_value * (ZGas - ZLiq) + vf * dZgasdP + (1 - vf) * dZliqdP;

    scalar drhodT_value = -p * mw_mixture / (ZMixture * ZMixture * RR * 1.0e-03 * T * T) * (dZMixturedT * T + ZMixture);
    scalar drhodP_value = mw_mixture / (ZMixture * RR * 1.0e-03 * T) - p * mw_mixture / (ZMixture * ZMixture * RR * 1.0e-03 * T) * dZMixturedP;

    scalar dEdT_ = Cp_Hs_value + p / (rho_value * rho_value) * drhodT_value;
    scalar dEdP_ = dHsdP_value - 1 / rho_value + p / (rho_value * rho_value) * drhodP_value;

    scalar det = dEdT_ * drhodP_value - dEdP_ * drhodT_value;

    scalar dTdE_rhoX_value = drhodP_value / det;
    scalar dPdE_rhoX_value = -drhodT_value / det;

    scalar dTdrho_EX_value = -dEdP_ / det;
    scalar dPdrho_EX_value = dEdT_ / det;

    ret_solution.dZdT = dZMixturedT;
    ret_solution.dZdP = dZMixturedP;
    ret_solution.dHsdT = Cp_Hs_value;
    ret_solution.dHsdP = dHsdP_value;
    ret_solution.dvfdT = dvfdT_value;
    ret_solution.dvfdP = dvfdP_value;

    return make_tuple(dTdE_rhoX_value, dPdE_rhoX_value, dTdrho_EX_value, dPdrho_EX_value);
}
template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::Z_NIO(scalar p, scalar T, VLE_data &sol) const
{
    scalar &vf_ = sol.vaporFraction;
    //scalar Z_G_0 = this->ThermoMixture::Z(p, T, sol.Gas.X, 1);
    //scalar Z_L_0 = this->ThermoMixture::Z(p, T, sol.Liq.X, 0);
    scalar &Z_G_ = sol.Gas.Z;
    scalar &Z_L_ = sol.Liq.Z;
    return vf_ * Z_G_ + (1 - vf_) * Z_L_;
}
template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::rho_NIO(scalar p, scalar T, VLE_data &sol)
    const
{

    scalar Z_ = Z_NIO(p, T, sol);
    scalar rho_M_ = p * this->W() / (Z_ * RR * 1.0e-03 * T);
    return rho_M_;
}
template <class ThermoMixture>
inline Foam::scalar Foam::VLE<ThermoMixture>::Ha_NIO(scalar p, scalar T, VLE_data &sol) const
{
    scalar mw_gas = this->W(sol.Gas.X); //kg/mol
    scalar mw_mixture = this->W();
    scalar ygas = sol.vaporFraction * mw_gas / mw_mixture;
    return this->ThermoMixture::Ha(p, T, sol.Gas.Z, sol.Gas.X, 1) * ygas + this->ThermoMixture::Ha(p, T, sol.Liq.Z, sol.Liq.X, 0) * (1.0 - ygas);
}
template <class ThermoMixture>
inline Foam::scalar Foam::VLE<ThermoMixture>::Hs_NIO(
    const scalar p,
    const scalar T, VLE_data &sol) const
{
    return Ha_NIO(p, T, sol) - Hc();
}
template <class ThermoMixture>
inline Foam::scalar Foam::VLE<ThermoMixture>::Es_NIO(
    const scalar p,
    const scalar T,
    VLE_data &sol) const
{
    return Hs_NIO(p, T, sol) - p / this->rho_NIO(p, T, sol);
}

template <class ThermoMixture>
inline std::tuple<Foam::scalar, Foam::scalar> Foam::VLE<ThermoMixture>::UVn_flash_NIO(
    scalar e,
    scalar rho,
    scalar p0,
    scalar T0,
    VLE_D1_data &retSol)
{
    if (T0 < 0)
    {
        FatalErrorInFunction
            << "Negative initial temperature T0: " << T0
            << abort(FatalError);
    }
    const int maxIter_ = 100;
    //const scalar tol_ = 1.0e-4;
    scalar T1 = T0;
    scalar p1 = p0;
    //auto sol(this->TPN(p1, T1));
    // scalar e11 = Es_vle(p1, T1, sol());      //, e21, e12;
    // scalar rho11 = this->rho(p1, T1, sol()); //, rho21, rho12;
    //scalar T2 = T1 + 1, p2 = p1 + 10; //(this->*limit)(T0+1);
    //scalar Ttol = T0 * tol_, ptol = p0 * tol_;
    scalar de = 0, drho = 0; //, dedT = 0, dedp = 0, drhodT = 0, drhodp = 0;
    //scalar de_old = 1e10, drho_old = 1e10;
    //scalar dete;
    int iter = 0;
    scalar maxratio = 0.3;
    scalar dTdE, dPdE, dTdrho, dPdrho;
    if (p0 < 0)
    {
        p0 = 1e5;
    }

    scalar t_f;
    this->TPn_flash_New_TPD_Tudisco_NIO(p1, T1, this->X_, retSol);
    scalar e11 = Es_NIO(p1, T1, retSol);
    scalar rho11 = this->rho_NIO(p1, T1, retSol);
    scalar minf = sqr(e11 - e) + sqr(rho11 - rho);
    std::tie(dTdE, dPdE, dTdrho, dPdrho) = this->dTPdErho_NIO(p1, T1, retSol);
    //std::tie(dTdE, dPdE, dTdrho, dPdrho) = this->dTPdErho(p1, T1, sol());

    scalar dT = 1e-2;
    scalar dp = 1;

    do
    {

        de = (e11 - e);
        drho = (rho11 - rho);
        dT = -(dTdE * de + dTdrho * drho);
        dp = -(dPdE * de + dPdrho * drho);
        if (mag(dT) > maxratio * T1)
        {
            dp *= maxratio * T1 / mag(dT);
            dT = sign(dT) * maxratio * T1;
        }
        if (mag(dp) > maxratio * p1)
        {
            dT *= maxratio * p1 / mag(dp);
            dp = sign(dp) * maxratio * p1;
        }

        //T1 = T1 > 1 ? T1 : 1;
        //p1 = p1 > 100 ? p1 : 100;
        T1 = T1 + dT;
        p1 = p1 + dp;
        TPn_flash_New_TPD_Tudisco_NIO(p1, T1, this->X_, retSol);
        //retSol.reset(this->TPN(p1, T1).ptr());
        e11 = Es_NIO(p1, T1, retSol);
        rho11 = this->rho_NIO(p1, T1, retSol);
        t_f = sqr(e11 - e) + sqr(rho11 - rho);

        if (t_f > minf)
        {
            //maxratio *= 0.9;
            if (maxratio > 1e-5)
            {
                maxratio *= 0.9;
            }
        }

        std::tie(dTdE, dPdE, dTdrho, dPdrho) = this->dTPdErho_NIO(p1, T1, retSol);
        if (t_f < minf)
        {
            minf = t_f;
        }

        if (iter++ > 5 * maxIter_)
        {
            /* WarningInFunction << setprecision(50)
                              << "Maximum number of iterations exceeded: " << iter - 1
                              << " error:" << sqrt(minf) //mag(y2[1] - e) + mag(y2[0] - rho)
                              << " e: " << e
                              << " rho: " << rho
                              << " T0: " << T0
                              << "  p0: " << p0
                              << " T1: " << T1
                              << "  p1: " << p1
                              << " X: " << this->X_
                              << setprecision(6)
                              << endl;
                              */
            //<< abort(FatalError);
            //continue;
            break;
        }
    } while ((mag(e11 - e) > 1e-2 || mag(rho11 - rho) > 1e-4));

    return std::make_tuple(p1, T1);
}

template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::c_NIO(scalar p, scalar T, scalar rho_, VLE_D1_data &sol) const
{

    scalar W_ = this->W();
    scalar vf_ = sol.vaporFraction;

    scalar A_G_ = this->ThermoMixture::A(p, T, sol.Gas.X);
    scalar A_L_ = this->ThermoMixture::A(p, T, sol.Liq.X);

    scalar B_G_ = this->ThermoMixture::B(p, T, sol.Gas.X);
    scalar B_L_ = this->ThermoMixture::B(p, T, sol.Liq.X);

    scalar Z_G_ = this->ThermoMixture::Z(p, T, A_G_, B_G_, sol.Gas.X, 1);
    scalar Z_L_ = this->ThermoMixture::Z(p, T, A_L_, B_L_, sol.Liq.X, 0);
    scalar Z_M_ = vf_ * Z_G_ + (1 - vf_) * Z_L_;

    scalar dZdP_M_ = sol.dZdP;

    scalar drhodP_ = W_ / (Z_M_ * RR * 1.0e-03 * T) - p * W_ / (Z_M_ * Z_M_ * RR * 1.0e-03 * T) * dZdP_M_;
    scalar kappaT_ = drhodP_ / rho_;

    scalar dZdT_M_ = sol.dZdT;

    scalar Cp_ = sol.dHsdT;

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
inline autoPtr<scalarSquareMatrix>
Foam::VLE<ThermoMixture>::A_dvid_NIO(scalar p, scalar T, VLE_data &sol) const
{
    scalarSquareMatrix &A = *(new scalarSquareMatrix(this->X_.size(), Zero));
    scalar vf_ = sol.vaporFraction;

    List<autoPtr<scalarList>> ddxi_Ln_fugcoef_G_(this->X_.size()),
        ddxi_Ln_fugcoef_L_(this->X_.size());

    for (int i = 0; i < this->X_.size(); i++)
        ddxi_Ln_fugcoef_G_[i].reset(
            (this->ddxi_Ln_fugacityCoefficient(p, T, i, sol.Gas.X, 1)).ptr());

    for (int i = 0; i < this->X_.size(); i++)
        ddxi_Ln_fugcoef_L_[i].reset(
            (this->ddxi_Ln_fugacityCoefficient(p, T, i, sol.Liq.X, 0)).ptr());
    for (int i = 0; i < this->X_.size(); i++)
    {
        for (int j = 0; j < this->X_.size(); j++)
        {
            A(i, j) = 1;
            if (i == j)
                A(i, j) -= this->X_[i] / (sol.Gas.X[i] * sol.Liq.X[i]);
            for (int k = 0; k < this->X_.size(); k++)
            {
                if (j == k)
                    A(i, j) +=
                        (1 - vf_) * ddxi_Ln_fugcoef_G_[k]()[i] * (sol.Gas.X[k] - 1);
                else
                    A(i, j) += (1 - vf_) * ddxi_Ln_fugcoef_G_[k]()[i] * sol.Gas.X[k];
                if (j == k)
                    A(i, j) += vf_ * ddxi_Ln_fugcoef_L_[k]()[i] * (sol.Liq.X[k] - 1);
                else
                    A(i, j) += vf_ * ddxi_Ln_fugcoef_L_[k]()[i] * sol.Liq.X[k];
            }
        }
    }

    return autoPtr<scalarSquareMatrix>(&A);
}

template <class ThermoMixture>
inline autoPtr<scalarList>
Foam::VLE<ThermoMixture>::dvidP_NIO(scalar p, scalar T, VLE_D1_data &sol, scalarSquareMatrix &A, labelList &pivotIndices) const
{

    scalarList &b = *(new scalarList(this->X_.size()));
    scalar vf_ = sol.vaporFraction;
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
        (this->ddP_Ln_fugacityCoefficient(p, T, sol.Gas.X, 1)).ptr());

    ddP_Ln_fugcoef_L_.reset(
        (this->ddP_Ln_fugacityCoefficient(p, T, sol.Liq.X, 0)).ptr());

    for (int i = 0; i < this->X_.size(); i++)
        b[i] = (1 - vf_) * vf_ * (ddP_Ln_fugcoef_G_()[i] - ddP_Ln_fugcoef_L_()[i]);
    //LUsolve(A, b);
    LUBacksubstitute(A, pivotIndices, b);
    return autoPtr<scalarList>(&b);
}

template <class ThermoMixture>
inline autoPtr<scalarList>
Foam::VLE<ThermoMixture>::dvidT_NIO(scalar p, scalar T, VLE_D1_data &sol, scalarSquareMatrix &A, labelList &pivotIndices) const
{
    scalarList &b = *(new scalarList(this->X_.size()));
    scalar vf_ = sol.vaporFraction;
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
        (this->ddT_Ln_fugacityCoefficient(p, T, sol.Gas.X, 1)).ptr());

    ddT_Ln_fugcoef_L_.reset(
        (this->ddT_Ln_fugacityCoefficient(p, T, sol.Liq.X, 0)).ptr());

    for (int i = 0; i < this->X_.size(); i++)
        b[i] = (1 - vf_) * vf_ * (ddT_Ln_fugcoef_G_()[i] - ddT_Ln_fugcoef_L_()[i]);
    //LUsolve(A, b);
    LUBacksubstitute(A, pivotIndices, b);
    return autoPtr<scalarList>(&b);
}

template <class ThermoMixture>
inline autoPtr<scalarList>
Foam::VLE<ThermoMixture>::dvidXi_NIO(scalar p, scalar T, label di, VLE_D1_data &sol) const
{
    scalarSquareMatrix A(this->X_.size(), Zero);
    scalarList &b = *(new scalarList(this->X_.size()));

    scalar vf_ = sol.vaporFraction;
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
        (this->ddT_Ln_fugacityCoefficient(p, T, sol.Gas.X, 1)).ptr());

    for (int i = 0; i < this->X_.size(); i++)
        ddxi_Ln_fugcoef_G_[i].reset(
            (this->ddxi_Ln_fugacityCoefficient(p, T, i, sol.Gas.X, 1)).ptr());

    ddT_Ln_fugcoef_L_.reset(
        (this->ddT_Ln_fugacityCoefficient(p, T, sol.Liq.X, 0)).ptr());

    for (int i = 0; i < this->X_.size(); i++)
        ddxi_Ln_fugcoef_L_[i].reset(
            (this->ddxi_Ln_fugacityCoefficient(p, T, i, sol.Liq.X, 0)).ptr());
    for (int i = 0; i < this->X_.size(); i++)
    {
        for (int j = 0; j < this->X_.size(); j++)
        {
            A(i, j) = -1.0 / (vf_ * (1 - vf_));
            if (i == j)
                A(i, j) +=
                    1 / (sol.Liq.X[i] * (1 - vf_)) + 1 / (sol.Gas.X[i] * vf_);
            for (int k = 0; k < this->X_.size(); k++)
            {
                if (j == k)
                    A(i, j) += ddxi_Ln_fugcoef_G_[k]()[i] * (1 - sol.Gas.X[k]) / vf_;
                else
                    A(i, j) += ddxi_Ln_fugcoef_G_[k]()[i] * (0 - sol.Gas.X[k]) / vf_;
                if (j == k)
                    A(i, j) +=
                        ddxi_Ln_fugcoef_L_[k]()[i] * (1 - sol.Liq.X[k]) / (1 - vf_);
                else
                    A(i, j) +=
                        ddxi_Ln_fugcoef_L_[k]()[i] * (0 - sol.Liq.X[k]) / (1 - vf_);
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
                            ddxi_Ln_fugcoef_L_[s]()[k] * (1 - sol.Liq.X[s]) / (1 - vf_);
                    else
                        b[k] +=
                            ddxi_Ln_fugcoef_L_[s]()[k] * (0 - sol.Liq.X[s]) / (1 - vf_);
                }

                if (k == j)
                {
                    b[k] += 1 / sol.Liq.X[k] * (1 - sol.Liq.X[k]) / (1 - vf_);
                }
                else
                {
                    b[k] += 1 / sol.Liq.X[k] * (0 - sol.Liq.X[k]) / (1 - vf_);
                }
            }
        }
    }
    LUsolve(A, b);
    return autoPtr<scalarList>(&b);
}

template <class ThermoMixture>
inline void Foam::VLE<ThermoMixture>::dTPvfcdXErho_NIO(scalar p, scalar T, VLE_D1_data &sol, scalarRectangularMatrix &out) const
{
    int nsp = this->X_.size();
    scalarRectangularMatrix &grad = out; //*(new scalarRectangularMatrix(nsp + 2, 4));

    const VLE<ThermoMixture> &This = *this;
    scalar vf = sol.vaporFraction;
    scalarList &X_gas = sol.Gas.X;
    scalarList &X_liq = sol.Liq.X;
    scalar rho_value = rho_NIO(p, T, sol);

    List<autoPtr<scalarList>> dvidXi_value(nsp);
    for (int i = 0; i < nsp; i++)
    {
        dvidXi_value[i] = this->dvidXi_NIO(p, T, i, sol);
    }
    scalar dvfdP_value = sol.dvfdP;
    scalar dvfdT_value = sol.dvfdT;
    scalarList dvfdXi_value(nsp);
    scalarList dHsdxi_gas_value(nsp), dHsdxi_liq_value(nsp);
    for (int i = 0; i < nsp; i++)
    {
        dvfdXi_value[i] = -vf;
        for (int j = 0; j < nsp; j++)
        {
            dvfdXi_value[i] += dvidXi_value[i]()[j];
        }
        dHsdxi_gas_value[i] = this->dHsdxi(p, T, i, X_gas, 1);
        dHsdxi_liq_value[i] = this->dHsdxi(p, T, i, X_liq, 0);
    }

    scalarList dHsdXi_gas(nsp);
    scalarList dHsdXi_liq(nsp);

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

    scalarList dmw_gasdXi(nsp);
    for (int i = 0; i < nsp; i++)
    {
        dmw_gasdXi[i] = 0;
        for (int j = 0; j < nsp; j++)
        {
            dmw_gasdXi[i] += dXdXi_gas_value[i][j] * This[j].W() * 1.0e-03;
        }
    }

    scalarList dygasdXi(nsp);
    for (int i = 0; i < nsp; i++)
    {
        dygasdXi[i] = ((dvfdXi_value[i] + vf) * mw_gas + vf * dmw_gasdXi[i]) / mw_mixture - ygas / mw_mixture * This[i].W() * 1.0e-03;
    }

    scalar Hs_gas = this->ThermoMixture::Hs(p, T, X_gas, 1);
    scalar Hs_liq = this->ThermoMixture::Hs(p, T, X_liq, 0);
    //scalar Cp_Hs_value = sol.dHsdT;
    //scalar dHsdP_value = sol.dHsdP;
    scalarList dHsdXi_value(nsp);
    for (int i = 0; i < nsp; i++)
    {
        dHsdXi_value[i] = Hs_gas * dygasdXi[i] + dHsdXi_gas[i] * ygas - Hs_liq * dygasdXi[i] + dHsdXi_liq[i] * (1.0 - ygas);
    }

    scalar ZMixture = Z_NIO(p, T, sol);
    scalar ZGas = this->ThermoMixture::Z(p, T, X_gas, 1);
    scalar ZLiq = this->ThermoMixture::Z(p, T, X_liq, 0);

    scalarList dZgasdXi(nsp), dZliqdXi(nsp);
    scalarList dZdxi_gas_value(nsp), dZdxi_liq_value(nsp);
    for (int i = 0; i < nsp; i++)
    {
        dZdxi_gas_value[i] = this->ThermoMixture::dZdxi(p, T, i, X_gas, 1);
        dZdxi_liq_value[i] = this->ThermoMixture::dZdxi(p, T, i, X_liq, 0);
    }
    for (int i = 0; i < nsp; i++)
    {
        dZgasdXi[i] = 0;
        dZliqdXi[i] = 0;
        for (int j = 0; j < nsp; j++)
        {
            dZgasdXi[i] += dZdxi_gas_value[j] * dXdXi_gas_value[i][j];
            dZliqdXi[i] += dZdxi_liq_value[j] * dXdXi_liq_value[i][j];
        }
    }

    //scalar dZMixturedT = sol.dZdT;
    //scalar dZMixturedP = sol.dZdP;
    scalarList dZMixturedXi(nsp);

    for (int i = 0; i < nsp; i++)
    {
        //dZMixturedXi[i] = (dvfdXi_value[i] * ZGas * ZGas + vf * 2 * ZGas * dZgasdXi[i] + 2 * (1 - vf) * ZLiq * dZliqdXi[i] - dvfdXi_value[i] * ZLiq * ZLiq) / temp - (vf * ZGas * ZGas + (1 - vf) * ZLiq * ZLiq) / (temp * temp) * (dvfdXi_value[i] * ZGas + vf * dZgasdXi[i] + (1 - vf) * dZliqdXi[i] - dvfdXi_value[i] * ZLiq);
        dZMixturedXi[i] = dvfdXi_value[i] * (ZGas - ZLiq) + vf * dZgasdXi[i] + (1 - vf) * dZliqdXi[i];
    }

    scalar drhodT_value = -p * mw_mixture / (ZMixture * ZMixture * RR * 1.0e-03 * T * T) * (sol.dZdT * T + ZMixture);                        //dZMixturedT
    scalar drhodP_value = mw_mixture / (ZMixture * RR * 1.0e-03 * T) - p * mw_mixture / (ZMixture * ZMixture * RR * 1.0e-03 * T) * sol.dZdP; //dZMixturedP
    scalarList drhodXi_value(nsp);
    for (int i = 0; i < nsp; i++)
    {
        drhodXi_value[i] = p * (This[i].W() * 1.0e-03 - mw_mixture) / (ZMixture * RR * 1.0e-03 * T) - p * mw_mixture / (ZMixture * ZMixture * RR * 1.0e-03 * T) * dZMixturedXi[i];
    }

    List<scalarList> dYdXi_G_value(nsp);
    for (int i = 0; i < nsp; i++)
    {
        dYdXi_G_value[i].resize(nsp);
        for (int j = 0; j < nsp; j++)
        {
            dYdXi_G_value[i][j] = dXdXi_gas_value[i][j] * This[j].W() * 1.0e-03 / mw_gas - X_gas[j] * This[j].W() * 1.0e-03 / (mw_gas * mw_gas) * dmw_gasdXi[i];
        }
    }

    scalar dEdT_ = sol.dHsdT + p / (rho_value * rho_value) * drhodT_value;                 //Cp_Hs_value
    scalar dEdP_ = sol.dHsdP - 1 / rho_value + p / (rho_value * rho_value) * drhodP_value; //dHsdP_value
    scalarList dEdXi_(nsp);
    for (int i = 0; i < nsp; i++)
    {
        dEdXi_[i] = dHsdXi_value[i] + p / (rho_value * rho_value) * drhodXi_value[i];
        //scalar drhodXi_ = drhodXi(p, T, i, sol);
    }
    scalar det = dEdT_ * drhodP_value - dEdP_ * drhodT_value;

    for (int i = 0; i < nsp; i++)
    {
        scalar dTdXi_Erho_value = -(drhodP_value * dEdXi_[i] - dEdP_ * drhodXi_value[i]) / det;
        scalar dPdXi_Erho_value = -(-drhodT_value * dEdXi_[i] + dEdT_ * drhodXi_value[i]) / det;
        //grad[i][0] = -(drhodP_value * dEdXi_[i] - dEdP_ * drhodXi_value[i]) / det;
        //grad[i][1] = -(-drhodT_value * dEdXi_[i] + dEdT_ * drhodXi_value[i]) / det;
        grad[i][0] = dTdXi_Erho_value;
        grad[i][1] = dPdXi_Erho_value;
        grad[i][2] = dvfdXi_value[i] + dvfdT_value * dTdXi_Erho_value + dvfdP_value * dPdXi_Erho_value;
        grad[i][3] = 0;
    }

    scalar dTdE_rhoX_value = drhodP_value / det;
    scalar dPdE_rhoX_value = -drhodT_value / det;
    scalar dvfdE_rhoX_value = dvfdT_value * dTdE_rhoX_value + dvfdP_value * dPdE_rhoX_value;

    scalar dTdrho_EX_value = -dEdP_ / det;
    scalar dPdrho_EX_value = dEdT_ / det;
    scalar dvfdrho_EX_value = dvfdT_value * dTdrho_EX_value + dvfdP_value * dPdrho_EX_value;

    grad[nsp][0] = dTdE_rhoX_value;
    grad[nsp][1] = dPdE_rhoX_value;
    grad[nsp][2] = dvfdE_rhoX_value;
    grad[nsp][3] = 0;

    grad[nsp + 1][0] = dTdrho_EX_value;
    grad[nsp + 1][1] = dPdrho_EX_value;
    grad[nsp + 1][2] = dvfdrho_EX_value;
    grad[nsp + 1][3] = 0;

    //return autoPtr<scalarRectangularMatrix>(&grad);
}

template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::c_full_NIO(scalar p, scalar T, scalar rho_, VLE_D1_data &sol) const
{
    const scalar Tstd = 298.15;
    const scalar sqrt2 = sqrt(2.0);
    const scalar cc0 = 2 * sqrt(2.0);
    const scalar cc1 = 1.0 + sqrt(2.0);
    const scalar cc2 = 1.0 - sqrt(2.0);

    const scalar small8 = 1.e-08; //PY: need to be noticed

    scalarList &X_gas = sol.Gas.X;
    scalarList &X_liq = sol.Liq.X;

    scalar W_ = this->W();
    scalar vf_ = sol.vaporFraction;

    scalar A_G_ = this->ThermoMixture::A(p, T, X_gas);
    scalar A_L_ = this->ThermoMixture::A(p, T, X_liq);

    scalar B_G_ = this->ThermoMixture::B(p, T, X_gas);
    scalar B_L_ = this->ThermoMixture::B(p, T, X_liq);

    scalar B_G_2 = B_G_ * B_G_;
    scalar B_L_2 = B_L_ * B_L_;

    scalar Z_G_ = this->ThermoMixture::Z(p, T, A_G_, B_G_, X_gas, 1);
    scalar Z_L_ = this->ThermoMixture::Z(p, T, A_L_, B_L_, X_liq, 0);
    scalar Z_M_ = vf_ * Z_G_ + (1 - vf_) * Z_L_;

    scalar Z_G_2 = Z_G_ * Z_G_;
    scalar Z_L_2 = Z_L_ * Z_L_;

    scalar a2_G_ = B_G_ - 1.0;
    scalar a2_L_ = B_L_ - 1.0;

    scalar a1_G_ = A_G_ - 2.0 * B_G_ - 3.0 * B_G_ * B_G_;
    scalar a1_L_ = A_L_ - 2.0 * B_L_ - 3.0 * B_L_ * B_L_;

    scalar dAdP_G_ = this->dAdP(p, T, X_gas), dBdP_G_ = this->dBdP(p, T, X_gas);
    scalar dAdP_L_ = this->dAdP(p, T, X_liq), dBdP_L_ = this->dBdP(p, T, X_liq);

    scalar da2dP_G_ = dBdP_G_;
    scalar da2dP_L_ = dBdP_L_;

    scalar da1dP_G_ = dAdP_G_ - 2.0 * dBdP_G_ - 6.0 * B_G_ * dBdP_G_;
    scalar da1dP_L_ = dAdP_L_ - 2.0 * dBdP_L_ - 6.0 * B_L_ * dBdP_L_;

    scalar da0dP_G_ = -(dAdP_G_ * B_G_ + A_G_ * dBdP_G_) + 2 * B_G_ * dBdP_G_ + 3 * B_G_2 * dBdP_G_;
    scalar da0dP_L_ = -(dAdP_L_ * B_L_ + A_L_ * dBdP_L_) + 2 * B_L_ * dBdP_L_ + 3 * B_L_2 * dBdP_L_;

    scalar dZdP_G_0 = -(da2dP_G_ * Z_G_2 + da1dP_G_ * Z_G_ + da0dP_G_) / (3 * Z_G_2 + 2 * Z_G_ * a2_G_ + a1_G_);
    scalar dZdP_L_0 = -(da2dP_L_ * Z_L_2 + da1dP_L_ * Z_L_ + da0dP_L_) / (3 * Z_L_2 + 2 * Z_L_ * a2_L_ + a1_L_);

    scalar dAdT_G_ = this->ThermoMixture::dAdT(p, T, X_gas), dBdT_G_ = this->ThermoMixture::dBdT(p, T, X_gas);
    scalar dAdT_L_ = this->ThermoMixture::dAdT(p, T, X_liq), dBdT_L_ = this->ThermoMixture::dBdT(p, T, X_liq);

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
    /*static scalarList *ddxi_Ln_fugcoef_G_ = NULL;

    static scalarList *ddxi_Ln_fugcoef_L_ = NULL;

    if (ddxi_Ln_fugcoef_G_ == NULL)
    {
        ddxi_Ln_fugcoef_G_ = new scalarList[this->X_.size()];
        ddxi_Ln_fugcoef_L_ = new scalarList[this->X_.size()];

        for(int i =0;i<this->X_.size();i++)
        {
            ddxi_Ln_fugcoef_G_[i].resize(this->X_.size());
            ddxi_Ln_fugcoef_L_[i].resize(this->X_.size());
        }
    }
*/
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
        dAdxi_G_[i] = this->dAdxi(p, T, i, X_gas);
        dBdxi_G_[i] = this->dBdxi(p, T, i, X_gas);
        dAdxi_L_[i] = this->dAdxi(p, T, i, X_liq);
        dBdxi_L_[i] = this->dBdxi(p, T, i, X_liq);

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

    autoPtr<scalarList> pA2_G_(this->A2(p, T, X_gas));
    const scalarList &A2_G_ = pA2_G_();
    autoPtr<scalarList> pA2_L_(this->A2(p, T, X_liq));
    const scalarList &A2_L_ = pA2_L_();

    for (int i = 0; i < this->X_.size(); i++)
    {
        scalarList &ddxi_Ln_fugcoef = *(new scalarList(this->N_));
        ddxi_Ln_fugcoef_G_[i].reset(&ddxi_Ln_fugcoef);
        //scalarList &ddxi_Ln_fugcoef = ddxi_Ln_fugcoef_G_[i];

        scalar dAmixzdxi = dAdxi_G_[i];
        scalar dBmixzdxi = dBdxi_G_[i];

        autoPtr<scalarList> pdAmix2dxi(this->dA2dxi(p, T, i, X_gas));
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
        //scalarList &ddxi_Ln_fugcoef = ddxi_Ln_fugcoef_L_[i];

        scalar dAmixzdxi = dAdxi_L_[i];
        scalar dBmixzdxi = dBdxi_L_[i];

        autoPtr<scalarList> pdAmix2dxi(this->dA2dxi(p, T, i, X_liq));
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
                A(i, j) -= this->X_[i] / (X_gas[i] * X_liq[i]);
            for (int k = 0; k < this->X_.size(); k++)
            {
                if (j == k)
                    A(i, j) +=
                        (1 - vf_) * ddxi_Ln_fugcoef_G_[k]()[i] * (X_gas[k] - 1);
                else
                    A(i, j) += (1 - vf_) * ddxi_Ln_fugcoef_G_[k]()[i] * X_gas[k];
                if (j == k)
                    A(i, j) += vf_ * ddxi_Ln_fugcoef_L_[k]()[i] * (X_liq[k] - 1);
                else
                    A(i, j) += vf_ * ddxi_Ln_fugcoef_L_[k]()[i] * X_liq[k];
            }
        }
    }

    labelList pivotIndices(A.m());
    LUDecompose(A, pivotIndices);

    //autoPtr<scalarList> dvidp(new scalarList(this->X_.size()));
    scalarList dvidp(this->X_.size());

    int n_sp = 0;
    for (int i = 0; i < this->X_.size(); i++)
    {
        if (this->X_[i] > 0)
            n_sp++;
    }
    if (n_sp == 1)
    {
        for (int i = 0; i < this->X_.size(); i++)
            dvidp[i] = 0;
    }
    else
    {
        autoPtr<scalarList> ddP_Ln_fugcoef_G_(new scalarList(this->N_)), ddP_Ln_fugcoef_L_(new scalarList(this->N_));
        {

            autoPtr<scalarList> pdAmix2dP(this->dA2dP(p, T, X_gas));
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

            autoPtr<scalarList> pdAmix2dP(this->dA2dP(p, T, X_liq));
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
            dvidp[i] = (1 - vf_) * vf_ * (ddP_Ln_fugcoef_G_()[i] - ddP_Ln_fugcoef_L_()[i]);

        LUBacksubstitute(A, pivotIndices, dvidp);
    }

    autoPtr<scalarList> dXdP_G_(this->dXdP_gas_NIO(p, T, sol, dvidp));
    autoPtr<scalarList> dXdP_L_(this->dXdP_liq_NIO(p, T, sol, dvidp));
    scalar dvfdP_ = this->dvfdP_NIO(p, T, sol, dvidp);

    scalar dZdP_G_ = dZdP_G_0; //this->ThermoMixture::dZdP(p, T, X_gas, 1);
    scalar dZdP_L_ = dZdP_L_0; //this->ThermoMixture::dZdP(p, T, X_liq, 0);

    for (int i = 0; i < this->X_.size(); i++)
    {
        dZdP_G_ += dzdxi_G_[i] * dXdP_G_()[i];
        dZdP_L_ += dzdxi_L_[i] * dXdP_L_()[i];
    }

    scalar dZdP_M_ = dvfdP_ * (Z_G_ - Z_L_) + vf_ * dZdP_G_ + (1 - vf_) * dZdP_L_;

    scalar drhodP_ = W_ / (Z_M_ * RR * 1.0e-03 * T) - p * W_ / (Z_M_ * Z_M_ * RR * 1.0e-03 * T) * dZdP_M_;
    scalar kappaT_ = drhodP_ / rho_;

    //autoPtr<scalarList> dvidt(new scalarList(this->X_.size()));
    scalarList dvidt(this->X_.size());

    if (n_sp == 1 || vf_ == 0 || vf_ == 1)
    {
        for (int i = 0; i < this->X_.size(); i++)
            dvidt[i] = 0;
    }
    else
    {
        autoPtr<scalarList> ddT_Ln_fugcoef_G_(new scalarList(this->N_)), ddT_Ln_fugcoef_L_(new scalarList(this->N_));
        {

            autoPtr<scalarList> pdAmix2dT(this->dA2dT(p, T, X_gas));
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

            autoPtr<scalarList> pdAmix2dT(this->dA2dT(p, T, X_liq));
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
            dvidt[i] = (1 - vf_) * vf_ * (ddT_Ln_fugcoef_G_()[i] - ddT_Ln_fugcoef_L_()[i]);

        LUBacksubstitute(A, pivotIndices, dvidt);
    }

    autoPtr<scalarList> dXdT_G_(this->dXdT_gas_NIO(p, T, sol, dvidt));
    autoPtr<scalarList> dXdT_L_(this->dXdT_liq_NIO(p, T, sol, dvidt));
    scalar dvfdT_ = this->dvfdT_NIO(p, T, sol, dvidt);

    scalar dZdT_G_ = dZdT_G_0;
    scalar dZdT_L_ = dZdT_L_0;
    for (int i = 0; i < this->X_.size(); i++)
    {
        dZdT_G_ += dzdxi_G_[i] * dXdT_G_()[i];
        dZdT_L_ += dzdxi_L_[i] * dXdT_L_()[i];
    }

    scalar dZdT_M_ = dvfdT_ * (Z_G_ - Z_L_) + vf_ * dZdT_G_ + (1 - vf_) * dZdT_L_;

    scalar W_G_ = this->W(X_gas);
    scalar W_L_ = this->W(X_liq);

    scalar dmw_gasdT = this->dWdT_gas_NIO(p, T, sol, dvidt);

    scalar ygas = vf_ * W_G_ / W_;

    scalarList dHadxi_G_(this->N_);
    scalarList dHadxi_L_(this->N_);

    scalar ddxi_dH_G_ = 0, dH_G_ = 0;
    scalar ddxi_dH_L_ = 0, dH_L_ = 0;

    for (int i = 0; i < this->X_.size(); i++)
    {
        scalar dBmixdxi = Bi[i];

        if (xxx_G_ > 0)
        {
            dH_G_ = (Z_G_ - 1) * RR * 1e-3 * T + RR * 1e-3 * T * ((T * dAdT_G_ + A_G_) / (cc0 * B_G_)) * logxxx_G_;

            ddxi_dH_G_ = dzdxi_G_[i] * RR * 1e-3 * T;

            ddxi_dH_G_ += RR * 1e-3 * T * ((T * this->d2AdTdxi(p, T, i, X_gas) + this->dAdxi(p, T, i, X_gas)) / (cc0 * B_G_)) * logxxx_G_;
            ddxi_dH_G_ -= RR * 1e-3 * T * ((T * dAdT_G_ + A_G_) / (cc0 * B_G_2) * dBmixdxi) * logxxx_G_;
            ddxi_dH_G_ += RR * 1e-3 * T * ((T * dAdT_G_ + A_G_) / (cc0 * B_G_)) * ((dzdxi_G_[i] + cc1 * dBmixdxi) / (Z_G_ + cc1 * B_G_) - (dzdxi_G_[i] + cc2 * dBmixdxi) / (Z_G_ + cc2 * B_G_));
            ddxi_dH_G_ = ddxi_dH_G_ / W_G_ - dH_G_ / (W_G_ * W_G_) * (*this)[i].W() * 1e-3;

            dHadxi_G_[i] = ddxi_dH_G_ + this->dHidealdxi(p, T, i, X_gas) - this->dHidealdxi(1, Tstd, i, X_gas);
        }
        else
        {
            ddxi_dH_G_ = dzdxi_G_[i] * RR * 1e-3 * T;
            dH_G_ = (Z_G_ - 1) * RR * 1e-3 * T;
            ddxi_dH_G_ = ddxi_dH_G_ / W_G_ - dH_G_ / (W_G_ * W_G_) * (*this)[i].W() * 1e-3;

            dHadxi_G_[i] = ddxi_dH_G_ + this->dHidealdxi(p, T, i, X_gas) - this->dHidealdxi(1, Tstd, i, X_gas);
        }
        if (xxx_L_ > 0)
        {
            dH_L_ = (Z_L_ - 1) * RR * 1e-3 * T + RR * 1e-3 * T * ((T * dAdT_L_ + A_L_) / (cc0 * B_L_)) * logxxx_L_;

            ddxi_dH_L_ = dzdxi_L_[i] * RR * 1e-3 * T;

            ddxi_dH_L_ += RR * 1e-3 * T * ((T * this->d2AdTdxi(p, T, i, X_liq) + this->dAdxi(p, T, i, X_liq)) / (cc0 * B_L_)) * logxxx_L_;
            ddxi_dH_L_ -= RR * 1e-3 * T * ((T * dAdT_L_ + A_L_) / (cc0 * B_L_2) * dBmixdxi) * logxxx_L_;
            ddxi_dH_L_ += RR * 1e-3 * T * ((T * dAdT_L_ + A_L_) / (cc0 * B_L_)) * ((dzdxi_L_[i] + cc1 * dBmixdxi) / (Z_L_ + cc1 * B_L_) - (dzdxi_L_[i] + cc2 * dBmixdxi) / (Z_L_ + cc2 * B_L_));
            ddxi_dH_L_ = ddxi_dH_L_ / W_L_ - dH_L_ / (W_L_ * W_L_) * (*this)[i].W() * 1e-3;

            dHadxi_L_[i] = ddxi_dH_L_ + this->dHidealdxi(p, T, i, X_liq) - this->dHidealdxi(1, Tstd, i, X_liq);
        }
        else
        {
            ddxi_dH_L_ = dzdxi_L_[i] * RR * 1e-3 * T;
            dH_L_ = (Z_L_ - 1) * RR * 1e-3 * T;
            ddxi_dH_L_ = ddxi_dH_L_ / W_L_ - dH_L_ / (W_L_ * W_L_) * (*this)[i].W() * 1e-3;

            dHadxi_L_[i] = ddxi_dH_L_ + this->dHidealdxi(p, T, i, X_liq) - this->dHidealdxi(1, Tstd, i, X_liq);
        }
    }

    scalar dygasdT = (dvfdT_ * W_G_ + vf_ * dmw_gasdT) / W_;
    scalar dmw_gasdP = this->dWdP_gas_NIO(p, T, sol, dvidp);
    scalar dygasdP = (dvfdP_ * W_G_ + vf_ * dmw_gasdP) / W_;

    scalar bmix_G_ = this->b(p, T, X_gas);

    scalar ddT_dH_G_, ddP_dH_G_;
    if (xxx_G_ > 0)
    {
        ddT_dH_G_ = dZdT_G_0 * RR * 1e-3 * T + (Z_G_ - 1) * RR * 1e-3;
        ddT_dH_G_ += ((T * this->d2adt2(p, T, X_gas)) / (cc0 * bmix_G_)) * logxxx_G_;
        ddT_dH_G_ += ((T * this->dadt(p, T, X_gas) - this->a(p, T, X_gas)) / (cc0 * bmix_G_)) * ((dZdT_G_0 - cc1 * B_G_ / T) / (Z_G_ + cc1 * B_G_) - (dZdT_G_0 - cc2 * B_G_ / T) / (Z_G_ + cc2 * B_G_));

        ddP_dH_G_ = dZdP_G_0 * RR * 1e-3 * T;
        ddP_dH_G_ += ((T * this->dadt(p, T, X_gas) - this->a(p, T, X_gas)) / (cc0 * bmix_G_)) * ((dZdP_G_0 + cc1 * B_G_ / p) / (Z_G_ + cc1 * B_G_) - (dZdP_G_0 + cc2 * B_G_ / p) / (Z_G_ + cc2 * B_G_));
    }
    else
    {
        ddP_dH_G_ = dZdP_G_0 * RR * 1e-3 * T;
        ddT_dH_G_ = dZdT_G_0 * RR * 1e-3 * T + (Z_G_ - 1) * RR * 1e-3;
    }

    ddT_dH_G_ = ddT_dH_G_ / W_G_; // J/kg
    // Based on Eq.12, Internal energy is estimated based on ideal(close to absolute 0 T) + departure

    ddP_dH_G_ = ddP_dH_G_ / W_G_;

    scalar bmix_L_ = this->b(p, T, X_liq);

    scalar ddT_dH_L_, ddP_dH_L_;

    if (xxx_L_ > 0)
    {
        ddT_dH_L_ = dZdT_L_0 * RR * 1e-3 * T + (Z_L_ - 1) * RR * 1e-3;
        ddT_dH_L_ += ((T * this->d2adt2(p, T, X_liq)) / (cc0 * bmix_L_)) * logxxx_L_;
        ddT_dH_L_ += ((T * this->dadt(p, T, X_liq) - this->a(p, T, X_liq)) / (cc0 * bmix_L_)) * ((dZdT_L_0 - cc1 * B_L_ / T) / (Z_L_ + cc1 * B_L_) - (dZdT_L_0 - cc2 * B_L_ / T) / (Z_L_ + cc2 * B_L_));

        ddP_dH_L_ = dZdP_L_0 * RR * 1e-3 * T;
        ddP_dH_L_ += ((T * this->dadt(p, T, X_liq) - this->a(p, T, X_liq)) / (cc0 * bmix_L_)) * ((dZdP_L_0 + cc1 * B_L_ / p) / (Z_L_ + cc1 * B_L_) - (dZdP_L_0 + cc2 * B_L_ / p) / (Z_L_ + cc2 * B_L_));
    }
    else
    {
        ddP_dH_L_ = dZdP_L_0 * RR * 1e-3 * T;
        ddT_dH_L_ = dZdT_L_0 * RR * 1e-3 * T + (Z_L_ - 1) * RR * 1e-3;
    }

    ddT_dH_L_ = ddT_dH_L_ / W_L_; // J/kg

    ddP_dH_L_ = ddP_dH_L_ / W_L_;

    scalar dHadT_gas = ddT_dH_G_ + this->dHidealdT(p, T, X_gas);
    scalar dHadP_gas = ddP_dH_G_;

    scalar dHadT_liq = ddT_dH_L_ + this->dHidealdT(p, T, X_liq);
    scalar dHadP_liq = ddP_dH_L_;

    for (int i = 0; i < this->X_.size(); i++)
    {
        dHadT_gas += dHadxi_G_[i] * dXdT_G_()[i];
        dHadP_gas += dHadxi_G_[i] * dXdP_G_()[i];
    }

    for (int i = 0; i < this->X_.size(); i++)
    {
        dHadT_liq += dHadxi_L_[i] * dXdT_L_()[i];
        dHadP_liq += dHadxi_L_[i] * dXdP_L_()[i];
    }
    //const scalar Tstd = 298.15;
    scalar Ha_G_ = dH_G_ / W_G_ + this->Hideal(p, T, X_gas) - this->Hideal(1, Tstd, X_gas), Ha_L_ = dH_L_ / W_L_ + this->Hideal(p, T, X_liq) - this->Hideal(1, Tstd, X_liq);

    scalar Cp_ = Ha_G_ * dygasdT + dHadT_gas * ygas - Ha_L_ * dygasdT + dHadT_liq * (1.0 - ygas);
    scalar dHsdP_ = Ha_G_ * dygasdP + dHadP_gas * ygas - Ha_L_ * dygasdP + dHadP_liq * (1.0 - ygas);

    //scalar Cp_ = this->ThermoMixture::Ha(p, T, X_gas, 1) * dygasdT + dHadT_gas * ygas - this->ThermoMixture::Ha(p, T, X_liq, 0) * dygasdT + dHadT_liq * (1.0 - ygas);
    //Info << "2 Ha_G_= " << dHsdP_ << endl;
    scalar drhodT_ = -p * W_ / (Z_M_ * Z_M_ * RR * 1.0e-03 * T * T) * (dZdT_M_ * T + Z_M_);
    scalar alphaP_ = -drhodT_ / rho_;
    scalar kappaS_ = kappaT_ - T * alphaP_ * alphaP_ / (rho_ * Cp_);
    sol.dZdT = dZdT_M_;
    sol.dZdP = dZdP_M_;
    sol.dHsdT = Cp_;
    sol.dHsdP = dHsdP_;
    sol.dvfdT = dvfdT_;
    sol.dvfdP = dvfdP_;
    if (kappaS_ < 0)
    {
        //return c_noVLE(p, T, sol);
        FatalErrorInFunction << " kappaS_ =" << kappaS_ << ",p=" << p << ",T=" << T << ",X=" << this->X_ << exit(FatalError);
    }
    return sqrt(1.0 / (kappaS_ * rho_));
}

template <class ThermoMixture>
inline Foam::autoPtr<Foam::scalarList> Foam::VLE<ThermoMixture>::dXdP_gas_NIO(scalar p, scalar T, VLE_D1_data &sol, scalarList &dvidt) const
{
    scalarList &dXdT = *(new scalarList(this->X_.size()));
    //autoPtr<scalarList> dvidt;
    //dvidt.reset((this->dvidP(p, T, sol)).ptr());
    scalar vf = sol.vaporFraction;
    if (vf <= 0.0)
    {
        for (int i = 0; i < this->X_.size(); i++)
        {
            dXdT[i] = 0;
        }
    }
    else
    {
        scalar rvf = 1 / vf;
        for (int i = 0; i < this->X_.size(); i++)
        {
            dXdT[i] = 0;
            for (int j = 0; j < this->X_.size(); j++)
            {
                if (i == j)
                    dXdT[i] -= rvf * (sol.Gas.X[i] - 1) * dvidt[j];
                else
                    dXdT[i] -= rvf * sol.Gas.X[i] * dvidt[j];
            }
        }
    }
    return autoPtr<scalarList>(&dXdT);
}

template <class ThermoMixture>
inline Foam::autoPtr<Foam::scalarList> Foam::VLE<ThermoMixture>::dXdP_liq_NIO(scalar p, scalar T, VLE_D1_data &sol, scalarList &dvidt) const
{
    scalarList &dXdT = *(new scalarList(this->X_.size()));
    scalar vf = sol.vaporFraction;
    if (vf == 1)
    {
        for (int i = 0; i < this->X_.size(); i++)
        {
            dXdT[i] = 0;
        }
        return autoPtr<scalarList>(&dXdT);
    }
    //autoPtr<scalarList> dvidt;
    //dvidt.reset((this->dvidP(p, T, sol)).ptr());
    //scalar vf = sol.vaporfra;
    for (int i = 0; i < this->X_.size(); i++)
    {
        dXdT[i] = 0;
        for (int j = 0; j < this->X_.size(); j++)
        {
            if (i == j)
                dXdT[i] += 1.0 / (1 - vf) * (sol.Liq.X[i] - 1) * dvidt[j];
            else
                dXdT[i] += 1.0 / (1 - vf) * sol.Liq.X[i] * dvidt[j];
        }
    }
    return autoPtr<scalarList>(&dXdT);
}

template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::dvfdP_NIO(scalar p, scalar T, VLE_D1_data &sol, scalarList &dvidP_) const
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
inline Foam::autoPtr<Foam::scalarList> Foam::VLE<ThermoMixture>::dXdT_gas_NIO(scalar p, scalar T, VLE_D1_data &sol, scalarList &dvidt) const
{
    scalarList &dXdT = *(new scalarList(this->X_.size()));
    //autoPtr<scalarList> dvidt;
    //dvidt.reset((this->dvidT(p, T, sol)).ptr());
    scalar vf = sol.vaporFraction;
    if (vf <= 0.0)
    {
        for (int i = 0; i < this->X_.size(); i++)
        {
            dXdT[i] = 0;
        }
    }
    else
    {
        scalar rvf = 1 / vf;
        for (int i = 0; i < this->X_.size(); i++)
        {
            dXdT[i] = 0;
            for (int j = 0; j < this->X_.size(); j++)
            {
                if (i == j)
                    dXdT[i] -= rvf * (sol.Gas.X[i] - 1) * dvidt[j];
                else
                    dXdT[i] -= rvf * sol.Gas.X[i] * dvidt[j];
            }
        }
    }
    return autoPtr<scalarList>(&dXdT);
}

template <class ThermoMixture>
inline Foam::autoPtr<Foam::scalarList> Foam::VLE<ThermoMixture>::dXdT_liq_NIO(scalar p, scalar T, VLE_D1_data &sol, scalarList &dvidt) const
{
    scalarList &dXdT = *(new scalarList(this->X_.size()));
    //autoPtr<scalarList> dvidt;
    //dvidt.reset((this->dvidT(p, T, sol)).ptr());
    scalar vf = sol.vaporFraction;
    if (vf >= 1.0)
    {
        for (int i = 0; i < this->X_.size(); i++)
        {
            dXdT[i] = 0;
        }
    }
    else
    {
        for (int i = 0; i < this->X_.size(); i++)
        {
            dXdT[i] = 0;

            for (int j = 0; j < this->X_.size(); j++)
            {
                if (i == j)
                    dXdT[i] += 1.0 / (1 - vf) * (sol.Liq.X[i] - 1) * dvidt[j];
                else
                    dXdT[i] += 1.0 / (1 - vf) * sol.Liq.X[i] * dvidt[j];
            }
        }
    }
    return autoPtr<scalarList>(&dXdT);
}

template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::dvfdT_NIO(scalar p, scalar T, VLE_D1_data &sol, scalarList &dvidT_) const
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
inline Foam::scalar Foam::VLE<ThermoMixture>::dWdT_gas_NIO(scalar p, scalar T, VLE_D1_data &sol, scalarList &dvidt) const
{
    scalar dWdT = 0;
    const VLE<ThermoMixture> &This = *this;
    autoPtr<scalarList> dXdT;

    dXdT.reset((this->dXdT_gas_NIO(p, T, sol, dvidt)).ptr());
    for (int i = 0; i < this->X_.size(); i++)
    {
        dWdT += dXdT()[i] * This[i].W() * 1.0e-03;
    }
    return dWdT;
}

template <class ThermoMixture>
inline Foam::scalar Foam::VLE<ThermoMixture>::dWdP_gas_NIO(scalar p, scalar T, VLE_D1_data &sol, scalarList &dvidt) const
{
    scalar dWdT = 0;
    const VLE<ThermoMixture> &This = *this;
    autoPtr<scalarList> dXdT;

    dXdT.reset((this->dXdP_gas_NIO(p, T, sol, dvidt)).ptr());
    for (int i = 0; i < this->X_.size(); i++)
    {
        dWdT += dXdT()[i] * This[i].W() * 1.0e-03;
    }
    return dWdT;
}

template <class ThermoMixture>
inline void Foam::VLE<ThermoMixture>::dErhovfcdXTP_NIO(scalar p, scalar T, VLE_D1_data &sol, scalarRectangularMatrix &out) const
{
    int nsp = this->X_.size();
    //scalar temp;
    scalarRectangularMatrix &grad = out; //*(new scalarRectangularMatrix(nsp + 2, 4));

    const VLE<ThermoMixture> &This = *this;
    scalar vf = sol.vaporFraction;
    scalarList &X_gas = sol.Gas.X;
    scalarList &X_liq = sol.Liq.X;
    scalar rho_value = rho_NIO(p, T, sol);

    List<autoPtr<scalarList>> dvidXi_value(nsp);
    for (int i = 0; i < nsp; i++)
    {
        dvidXi_value[i] = this->dvidXi_NIO(p, T, i, sol);
    }
    scalar dvfdP_value = sol.dvfdP, dvfdT_value = sol.dvfdT;
    scalarList dvfdXi_value(nsp);
    scalarList dHsdxi_gas_value(nsp), dHsdxi_liq_value(nsp);
    for (int i = 0; i < nsp; i++)
    {
        dvfdXi_value[i] = -vf;
        for (int j = 0; j < nsp; j++)
        {
            dvfdXi_value[i] += dvidXi_value[i]()[j];
        }
        dHsdxi_gas_value[i] = this->dHsdxi(p, T, i, X_gas, 1);
        dHsdxi_liq_value[i] = this->dHsdxi(p, T, i, X_liq, 0);
    }

    scalarList dHsdXi_gas(nsp);
    scalarList dHsdXi_liq(nsp);

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

    scalarList dmw_gasdXi(nsp);
    for (int i = 0; i < nsp; i++)
    {
        dmw_gasdXi[i] = 0;
        for (int j = 0; j < nsp; j++)
        {
            dmw_gasdXi[i] += dXdXi_gas_value[i][j] * This[j].W() * 1.0e-03;
        }
    }

    scalarList dygasdXi(nsp);
    for (int i = 0; i < nsp; i++)
    {
        dygasdXi[i] = ((dvfdXi_value[i] + vf) * mw_gas + vf * dmw_gasdXi[i]) / mw_mixture - ygas / mw_mixture * This[i].W() * 1.0e-03;
    }

    scalar Hs_gas = this->ThermoMixture::Hs(p, T, X_gas, 1);
    scalar Hs_liq = this->ThermoMixture::Hs(p, T, X_liq, 0);
    scalar Cp_Hs_value = sol.dHsdT;
    scalar dHsdP_value = sol.dHsdP;
    scalarList dHsdXi_value(nsp);
    for (int i = 0; i < nsp; i++)
    {
        dHsdXi_value[i] = Hs_gas * dygasdXi[i] + dHsdXi_gas[i] * ygas - Hs_liq * dygasdXi[i] + dHsdXi_liq[i] * (1.0 - ygas);
    }

    scalar ZMixture = Z_NIO(p, T, sol);
    scalar ZGas = this->ThermoMixture::Z(p, T, X_gas, 1);
    scalar ZLiq = this->ThermoMixture::Z(p, T, X_liq, 0);
    scalarList dZgasdXi(nsp), dZliqdXi(nsp);
    scalarList dZdxi_gas_value(nsp), dZdxi_liq_value(nsp);
    for (int i = 0; i < nsp; i++)
    {
        dZdxi_gas_value[i] = this->ThermoMixture::dZdxi(p, T, i, X_gas, 1);
        dZdxi_liq_value[i] = this->ThermoMixture::dZdxi(p, T, i, X_liq, 0);
    }
    for (int i = 0; i < nsp; i++)
    {
        dZgasdXi[i] = 0;
        dZliqdXi[i] = 0;
        for (int j = 0; j < nsp; j++)
        {
            dZgasdXi[i] += dZdxi_gas_value[j] * dXdXi_gas_value[i][j];
            dZliqdXi[i] += dZdxi_liq_value[j] * dXdXi_liq_value[i][j];
        }
    }
    //temp = (vf * ZGas + (1 - vf) * ZLiq);
    scalar dZMixturedT = sol.dZdT;
    scalar dZMixturedP = sol.dZdP;
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

    for (int i = 0; i < nsp; i++)
    {
        grad[i][0] = dHsdXi_value[i] + p / (rho_value * rho_value) * drhodXi_value[i];
        grad[i][1] = drhodXi_value[i];
        grad[i][2] = dvfdXi_value[i];
        grad[i][3] = 0;
    }
    grad[nsp][0] = Cp_Hs_value + p / (rho_value * rho_value) * drhodT_value;
    grad[nsp][1] = drhodT_value;
    grad[nsp][2] = dvfdT_value;
    grad[nsp][3] = 0;
    grad[nsp + 1][0] = dHsdP_value - 1 / rho_value + p / (rho_value * rho_value) * drhodP_value;
    grad[nsp + 1][1] = drhodP_value;
    grad[nsp + 1][2] = dvfdP_value;
    grad[nsp + 1][3] = 0;
    //return autoPtr<scalarRectangularMatrix>(&grad);
}

template <class ThermoMixture>
inline Foam::scalar Foam::VLE<ThermoMixture>::PV_nflash_NIO(
    scalar p,
    scalar rho,
    scalar T0,
    VLE_D1_data &sol) const
{

    if (T0 < 0)
    {
        FatalErrorInFunction
            << "Negative initial temperature T0: " << T0
            << abort(FatalError);
    }
    const label maxIter_ = 100;
    scalar T1 = T0;
    this->TPn_flash_New_TPD_Tudisco_NIO(p, T1, this->X_, sol);
    scalar f = 1 / rho;

    //sol.reset(this->TPN(p, T1).ptr());
    scalar F1 = rrho_NIO(p, T1, sol), F2 = F1;
    scalar T2 = T0 + 1;
    scalar T2r = T2;
    const scalar tol_ = 1.0e-4;
    scalar Ttol = T0 * tol_;
    scalar dF = 0, dfdt = 0, dFold = 1e100;
    int iter = 0;
    scalar relax = 1;
    bool inSection = false;
    bool bisec = false;
    double Tm, Fm;
    do
    {
        F1 = F2;
        this->TPn_flash_New_TPD_Tudisco_NIO(p, T2, this->X_, sol);
        T2r = T2;
        //sol.reset(this->TPN(p, T2).ptr());
        F2 = rrho_NIO(p, T2, sol);
        dF = (F2 - f);

        if (fabs(dF) > fabs(dFold))
            relax *= 0.85;

        dFold = dF;
        inSection = ((F2 - f) * (F1 - f) < 0);

        dfdt = (F2 - F1) / (T2 - T1); //(this->*dFdT)(p, Test);

        if (!inSection)
        {
            T1 = T2;
            if (mag(relax * dF / dfdt) > 30)
                T2 = T1 - sign(relax * dF / dfdt) * 30;
            else
                T2 = T1 - relax * dF / dfdt;
        }
        else
        {
            do
            {
                dfdt = (F2 - F1) / (T2 - T1);
                if (bisec)
                {
                    Tm = (T2 + T1) / 2;
                }
                else
                {
                    Tm = T2 - (F2 - f) / dfdt;
                }
                //sol.reset(this->TPN(p, Tm).ptr());
                this->TPn_flash_New_TPD_Tudisco_NIO(p, Tm, this->X_, sol);
                T2r = Tm;
                Fm = rrho_NIO(p, Tm, sol);
                if ((Fm - f) * (F1 - f) < 0)
                {
                    bisec = ((T2 - Tm) / (T2 - T1) < 0.1);
                    T2 = Tm;
                    F2 = Fm;
                }
                else
                {
                    bisec = ((T1 - Tm) / (T1 - T2) < 0.1);
                    T1 = Tm;
                    F1 = Fm;
                }
                if (iter++ > maxIter_)
                {
                    FatalErrorInFunction << setprecision(50)
                                         << "Maximum number of iterations exceeded: " << maxIter_
                                         << " f: " << f
                                         << " p: " << p
                                         << " T0: " << T0
                                         << "  Tg: " << T2
                                         << " X: " << this->X_
                                         << abort(FatalError);
                }
            } while (mag(T2 - T1) > Ttol && mag(rho - 1 / Fm) > 1e-4);
            T2 = Tm;
            F2 = Fm;
        }
        //Tnew =
        //   (this->*limit)(Test - ((this->*F)(p, Test) - f) / (this->*dFdT)(p, Test));

        if (T2 == 0)
        {
            FatalErrorInFunction << setprecision(50)
                                 << "T2 == 0: " << maxIter_
                                 << " f: " << f
                                 << " p: " << p
                                 << " T0: " << T0
                                 << "  Tg: " << T2
                                 << " X: " << this->X_
                                 << abort(FatalError);
        }
        if (iter++ > maxIter_)
        {
            FatalErrorInFunction << setprecision(50)
                                 << "Maximum number of iterations exceeded: " << maxIter_
                                 << " f: " << f
                                 << " p: " << p
                                 << " T0: " << T0
                                 << "  Tg: " << T2
                                 << " X: " << this->X_
                                 << abort(FatalError);
        }

    } while (mag(T2 - T1) > Ttol && mag(1 / F2 - rho) > 1e-4);
    return T2r;
}

template <class ThermoMixture>
inline Foam::scalar Foam::VLE<ThermoMixture>::EP_nflash_NIO(
    scalar E,
    scalar p,
    scalar T0,
    VLE_D1_data &sol) const
{

    if (T0 < 0)
    {
        FatalErrorInFunction
            << "Negative initial temperature T0: " << T0
            << abort(FatalError);
    }
    const label maxIter_ = 100;
    scalar T1 = T0;
    this->TPn_flash_New_TPD_Tudisco_NIO(p, T1, this->X_, sol);
    scalar f = E;

    //sol.reset(this->TPN(p, T1).ptr());
    scalar F1 = Es_NIO(p, T1, sol), F2 = F1;
    scalar T2 = T0 + 1;
    scalar T2r = T2;
    const scalar tol_ = 1.0e-4;
    scalar Ttol = T0 * tol_;
    scalar dF = 0, dfdt = 0, dFold = 1e100;
    int iter = 0;
    scalar relax = 1;
    bool inSection = false;
    bool bisec = false;
    double Tm, Fm;
    do
    {
        F1 = F2;
        this->TPn_flash_New_TPD_Tudisco_NIO(p, T2, this->X_, sol);
        T2r = T2;
        //sol.reset(this->TPN(p, T2).ptr());
        F2 = Es_NIO(p, T2, sol);
        dF = (F2 - f);

        if (fabs(dF) > fabs(dFold))
            relax *= 0.85;

        dFold = dF;
        inSection = ((F2 - f) * (F1 - f) < 0);

        dfdt = (F2 - F1) / (T2 - T1); //(this->*dFdT)(p, Test);

        if (!inSection)
        {
            T1 = T2;
            if (mag(relax * dF / dfdt) > 30)
                T2 = T1 - sign(relax * dF / dfdt) * 30;
            else
                T2 = T1 - relax * dF / dfdt;
        }
        else
        {
            do
            {
                dfdt = (F2 - F1) / (T2 - T1);
                if (bisec)
                {
                    Tm = (T2 + T1) / 2;
                }
                else
                {
                    Tm = T2 - (F2 - f) / dfdt;
                }
                //sol.reset(this->TPN(p, Tm).ptr());
                this->TPn_flash_New_TPD_Tudisco_NIO(p, Tm, this->X_, sol);
                T2r = Tm;
                Fm = Es_NIO(p, Tm, sol);
                if ((Fm - f) * (F1 - f) < 0)
                {
                    bisec = ((T2 - Tm) / (T2 - T1) < 0.1);
                    T2 = Tm;
                    F2 = Fm;
                }
                else
                {
                    bisec = ((T1 - Tm) / (T1 - T2) < 0.1);
                    T1 = Tm;
                    F1 = Fm;
                }
                if (iter++ > maxIter_)
                {
                    FatalErrorInFunction << setprecision(50)
                                         << "Maximum number of iterations exceeded: " << maxIter_
                                         << " f: " << f
                                         << " p: " << p
                                         << " T0: " << T0
                                         << "  Tg: " << T2
                                         << " X: " << this->X_
                                         << abort(FatalError);
                }
            } while (mag(T2 - T1) > Ttol && mag(f - Fm) > 1e-8 * f);
            T2 = Tm;
        }
        //Tnew =
        //   (this->*limit)(Test - ((this->*F)(p, Test) - f) / (this->*dFdT)(p, Test));

        if (T2 == 0)
        {
            FatalErrorInFunction << setprecision(50)
                                 << "T2 == 0: " << maxIter_
                                 << " f: " << f
                                 << " p: " << p
                                 << " T0: " << T0
                                 << "  Tg: " << T2
                                 << " X: " << this->X_
                                 << abort(FatalError);
        }
        if (iter++ > maxIter_)
        {
            FatalErrorInFunction << setprecision(50)
                                 << "Maximum number of iterations exceeded: " << maxIter_
                                 << " f: " << f
                                 << " p: " << p
                                 << " T0: " << T0
                                 << "  Tg: " << T2
                                 << " X: " << this->X_
                                 << abort(FatalError);
        }

    } while (mag(T2 - T1) > Ttol && mag(f - F2) > 1e-6 * f);
    return T2r;
}

template <class ThermoMixture>
inline Foam::scalar
Foam::VLE<ThermoMixture>::rrho_NIO(scalar p, scalar T, VLE_D1_data &sol) const
{
    // autoPtr<solution> sol(TPn_flash(p, T));
    scalar Z_ = Z_NIO(p, T, sol);
    scalar rrho_M_ = (Z_ * RR * 1.0e-03 * T) / (p * this->W());
    return rrho_M_;
}

template <class ThermoMixture>
inline Foam::autoPtr<scalarRectangularMatrix> Foam::VLE<ThermoMixture>::dTHvfcdXrhoP_NIO(scalar p, scalar T, VLE_D1_data &sol) const
{
    int nsp = this->X_.size();
    scalarRectangularMatrix &grad = *(new scalarRectangularMatrix(nsp + 2, 4));

    const VLE<ThermoMixture> &This = *this;
    scalar vf = sol.vaporFraction;
    scalarList &X_gas = sol.Gas.X;
    scalarList &X_liq = sol.Liq.X;
    scalar rho_value = rho_NIO(p, T, sol);

    List<autoPtr<scalarList>> dvidXi_value(nsp);
    for (int i = 0; i < nsp; i++)
    {
        dvidXi_value[i] = this->dvidXi_NIO(p, T, i, sol);
    }
    scalar dvfdP_value = sol.dvfdP, dvfdT_value = sol.dvfdT;
    scalarList dvfdXi_value(nsp);
    scalarList dHsdxi_gas_value(nsp), dHsdxi_liq_value(nsp);
    for (int i = 0; i < nsp; i++)
    {
        dvfdXi_value[i] = -vf;
        for (int j = 0; j < nsp; j++)
        {
            dvfdXi_value[i] += dvidXi_value[i]()[j];
        }
        dHsdxi_gas_value[i] = this->dHsdxi(p, T, i, X_gas, 1);
        dHsdxi_liq_value[i] = this->dHsdxi(p, T, i, X_liq, 0);
    }

    scalarList dHsdXi_gas(nsp);
    scalarList dHsdXi_liq(nsp);

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

    scalarList dmw_gasdXi(nsp);
    for (int i = 0; i < nsp; i++)
    {
        dmw_gasdXi[i] = 0;
        for (int j = 0; j < nsp; j++)
        {
            dmw_gasdXi[i] += dXdXi_gas_value[i][j] * This[j].W() * 1.0e-03;
        }
    }
    scalarList dygasdXi(nsp);
    for (int i = 0; i < nsp; i++)
    {
        dygasdXi[i] = ((dvfdXi_value[i] + vf) * mw_gas + vf * dmw_gasdXi[i]) / mw_mixture - ygas / mw_mixture * This[i].W() * 1.0e-03;
    }

    scalar Hs_gas = this->ThermoMixture::Hs(p, T, X_gas, 1);
    scalar Hs_liq = this->ThermoMixture::Hs(p, T, X_liq, 0);
    scalar Cp_Hs_value = sol.dHsdT;
    scalar dHsdP_value = sol.dHsdP;
    scalarList dHsdXi_value(nsp);
    for (int i = 0; i < nsp; i++)
    {
        dHsdXi_value[i] = Hs_gas * dygasdXi[i] + dHsdXi_gas[i] * ygas - Hs_liq * dygasdXi[i] + dHsdXi_liq[i] * (1.0 - ygas);
    }

    scalar ZMixture = Z_NIO(p, T, sol);
    scalar ZGas = this->ThermoMixture::Z(p, T, X_gas, 1);
    scalar ZLiq = this->ThermoMixture::Z(p, T, X_liq, 0);

    scalarList dZgasdXi(nsp), dZliqdXi(nsp);
    scalarList dZdxi_gas_value(nsp), dZdxi_liq_value(nsp);
    for (int i = 0; i < nsp; i++)
    {
        dZdxi_gas_value[i] = this->ThermoMixture::dZdxi(p, T, i, X_gas, 1);
        dZdxi_liq_value[i] = this->ThermoMixture::dZdxi(p, T, i, X_liq, 0);
    }
    for (int i = 0; i < nsp; i++)
    {
        dZgasdXi[i] = 0;
        dZliqdXi[i] = 0;
        for (int j = 0; j < nsp; j++)
        {
            dZgasdXi[i] += dZdxi_gas_value[j] * dXdXi_gas_value[i][j];
            dZliqdXi[i] += dZdxi_liq_value[j] * dXdXi_liq_value[i][j];
        }
    }

    scalar dZMixturedT = sol.dZdT;
    scalar dZMixturedP = sol.dZdP;
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

    for (int i = 0; i < nsp; i++)
    {
        grad[i][0] = -drhodXi_value[i] / drhodT_value;
        grad[i][1] = -Cp_Hs_value * drhodXi_value[i] / drhodT_value + dHsdXi_value[i];
        grad[i][2] = -dvfdT_value * drhodXi_value[i] / drhodT_value + dvfdXi_value[i];
        grad[i][3] = 0;
    }
    grad[nsp][0] = 1 / drhodT_value;
    grad[nsp][1] = Cp_Hs_value / drhodT_value;
    grad[nsp][2] = dvfdT_value / drhodT_value;
    grad[nsp][3] = 0;

    grad[nsp + 1][0] = -drhodP_value / drhodT_value;
    grad[nsp + 1][1] = -Cp_Hs_value * drhodP_value / drhodT_value + dHsdP_value;
    grad[nsp + 1][2] = -dvfdT_value * drhodP_value / drhodT_value + dvfdP_value;
    grad[nsp + 1][3] = 0;

    return autoPtr<scalarRectangularMatrix>(&grad);
}

template <class ThermoMixture>
inline void Foam::VLE<ThermoMixture>::dTEvfcdXrhoP_NIO(scalar p, scalar T, VLE_D1_data &sol, scalarRectangularMatrix &out) const
{
    int nsp = this->X_.size();
    scalarRectangularMatrix &grad = out;

    const VLE<ThermoMixture> &This = *this;
    scalar vf = sol.vaporFraction;
    scalarList &X_gas = sol.Gas.X;
    scalarList &X_liq = sol.Liq.X;
    scalar rho_value = rho_NIO(p, T, sol);

    List<autoPtr<scalarList>> dvidXi_value(nsp);
    for (int i = 0; i < nsp; i++)
    {
        dvidXi_value[i] = this->dvidXi_NIO(p, T, i, sol);
    }
    scalar dvfdP_value = sol.dvfdP, dvfdT_value = sol.dvfdT;
    scalarList dvfdXi_value(nsp);
    scalarList dHsdxi_gas_value(nsp), dHsdxi_liq_value(nsp);
    for (int i = 0; i < nsp; i++)
    {
        dvfdXi_value[i] = -vf;
        for (int j = 0; j < nsp; j++)
        {
            dvfdXi_value[i] += dvidXi_value[i]()[j];
        }
        dHsdxi_gas_value[i] = this->dHsdxi(p, T, i, X_gas, 1);
        dHsdxi_liq_value[i] = this->dHsdxi(p, T, i, X_liq, 0);
    }

    scalarList dHsdXi_gas(nsp);
    scalarList dHsdXi_liq(nsp);

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

    scalarList dmw_gasdXi(nsp);
    for (int i = 0; i < nsp; i++)
    {
        dmw_gasdXi[i] = 0;
        for (int j = 0; j < nsp; j++)
        {
            dmw_gasdXi[i] += dXdXi_gas_value[i][j] * This[j].W() * 1.0e-03;
        }
    }
    scalarList dygasdXi(nsp);
    for (int i = 0; i < nsp; i++)
    {
        dygasdXi[i] = ((dvfdXi_value[i] + vf) * mw_gas + vf * dmw_gasdXi[i]) / mw_mixture - ygas / mw_mixture * This[i].W() * 1.0e-03;
    }

    scalar Hs_gas = this->ThermoMixture::Hs(p, T, X_gas, 1);
    scalar Hs_liq = this->ThermoMixture::Hs(p, T, X_liq, 0);
    scalar Cp_Hs_value = sol.dHsdT;
    scalar dHsdP_value = sol.dHsdP;
    scalarList dHsdXi_value(nsp);
    for (int i = 0; i < nsp; i++)
    {
        dHsdXi_value[i] = Hs_gas * dygasdXi[i] + dHsdXi_gas[i] * ygas - Hs_liq * dygasdXi[i] + dHsdXi_liq[i] * (1.0 - ygas);
    }

    scalar ZMixture = Z_NIO(p, T, sol);
    scalar ZGas = this->ThermoMixture::Z(p, T, X_gas, 1);
    scalar ZLiq = this->ThermoMixture::Z(p, T, X_liq, 0);

    scalarList dZgasdXi(nsp), dZliqdXi(nsp);
    scalarList dZdxi_gas_value(nsp), dZdxi_liq_value(nsp);
    for (int i = 0; i < nsp; i++)
    {
        dZdxi_gas_value[i] = this->ThermoMixture::dZdxi(p, T, i, X_gas, 1);
        dZdxi_liq_value[i] = this->ThermoMixture::dZdxi(p, T, i, X_liq, 0);
    }
    for (int i = 0; i < nsp; i++)
    {
        dZgasdXi[i] = 0;
        dZliqdXi[i] = 0;
        for (int j = 0; j < nsp; j++)
        {
            dZgasdXi[i] += dZdxi_gas_value[j] * dXdXi_gas_value[i][j];
            dZliqdXi[i] += dZdxi_liq_value[j] * dXdXi_liq_value[i][j];
        }
    }

    scalar dZMixturedT = sol.dZdT;
    scalar dZMixturedP = sol.dZdP;
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

    for (int i = 0; i < nsp; i++)
    {
        grad[i][0] = -drhodXi_value[i] / drhodT_value;
        grad[i][1] = -Cp_Hs_value * drhodXi_value[i] / drhodT_value + dHsdXi_value[i];
        grad[i][2] = -dvfdT_value * drhodXi_value[i] / drhodT_value + dvfdXi_value[i];
        grad[i][3] = 0;
    }
    grad[nsp][0] = 1 / drhodT_value;
    grad[nsp][1] = Cp_Hs_value / drhodT_value + p / sqr(rho_value);
    grad[nsp][2] = dvfdT_value / drhodT_value;
    grad[nsp][3] = 0;

    //Info << "drhodT=" << drhodT_value << endl;
    //Info << "drhodP=" << drhodP_value << endl;

    grad[nsp + 1][0] = -drhodP_value / drhodT_value;
    grad[nsp + 1][1] = -Cp_Hs_value * drhodP_value / drhodT_value + dHsdP_value - 1 / rho_value;
    grad[nsp + 1][2] = -dvfdT_value * drhodP_value / drhodT_value + dvfdP_value;
    grad[nsp + 1][3] = 0;

    //return autoPtr<scalarRectangularMatrix>(&grad);
}

template <class ThermoMixture>
inline Foam::autoPtr<scalarRectangularMatrix> Foam::VLE<ThermoMixture>::dTrhovfcdXEP_NIO(scalar p, scalar T, VLE_D1_data &sol) const
{
    int nsp = this->X_.size();
    scalarRectangularMatrix &grad = *(new scalarRectangularMatrix(nsp + 2, 4));

    const VLE<ThermoMixture> &This = *this;
    scalar vf = sol.vaporFraction;
    scalarList &X_gas = sol.Gas.X;
    scalarList &X_liq = sol.Liq.X;
    scalar rho_value = rho_NIO(p, T, sol);

    List<autoPtr<scalarList>> dvidXi_value(nsp);
    for (int i = 0; i < nsp; i++)
    {
        dvidXi_value[i] = this->dvidXi_NIO(p, T, i, sol);
    }
    scalar dvfdP_value = sol.dvfdP, dvfdT_value = sol.dvfdT;
    scalarList dvfdXi_value(nsp);
    scalarList dHsdxi_gas_value(nsp), dHsdxi_liq_value(nsp);
    for (int i = 0; i < nsp; i++)
    {
        dvfdXi_value[i] = -vf;
        for (int j = 0; j < nsp; j++)
        {
            dvfdXi_value[i] += dvidXi_value[i]()[j];
        }
        dHsdxi_gas_value[i] = this->dHsdxi(p, T, i, X_gas, 1);
        dHsdxi_liq_value[i] = this->dHsdxi(p, T, i, X_liq, 0);
    }

    scalarList dHsdXi_gas(nsp);
    scalarList dHsdXi_liq(nsp);

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

    scalarList dmw_gasdXi(nsp);
    for (int i = 0; i < nsp; i++)
    {
        dmw_gasdXi[i] = 0;
        for (int j = 0; j < nsp; j++)
        {
            dmw_gasdXi[i] += dXdXi_gas_value[i][j] * This[j].W() * 1.0e-03;
        }
    }
    scalarList dygasdXi(nsp);
    for (int i = 0; i < nsp; i++)
    {
        dygasdXi[i] = ((dvfdXi_value[i] + vf) * mw_gas + vf * dmw_gasdXi[i]) / mw_mixture - ygas / mw_mixture * This[i].W() * 1.0e-03;
    }

    scalar Hs_gas = this->ThermoMixture::Hs(p, T, X_gas, 1);
    scalar Hs_liq = this->ThermoMixture::Hs(p, T, X_liq, 0);
    scalar Cp_Hs_value = sol.dHsdT;
    scalar dHsdP_value = sol.dHsdP;
    scalarList dHsdXi_value(nsp);
    for (int i = 0; i < nsp; i++)
    {
        dHsdXi_value[i] = Hs_gas * dygasdXi[i] + dHsdXi_gas[i] * ygas - Hs_liq * dygasdXi[i] + dHsdXi_liq[i] * (1.0 - ygas);
    }

    scalar ZMixture = Z_NIO(p, T, sol);
    scalar ZGas = this->ThermoMixture::Z(p, T, X_gas, 1);
    scalar ZLiq = this->ThermoMixture::Z(p, T, X_liq, 0);

    scalarList dZgasdXi(nsp), dZliqdXi(nsp);
    scalarList dZdxi_gas_value(nsp), dZdxi_liq_value(nsp);
    for (int i = 0; i < nsp; i++)
    {
        dZdxi_gas_value[i] = this->ThermoMixture::dZdxi(p, T, i, X_gas, 1);
        dZdxi_liq_value[i] = this->ThermoMixture::dZdxi(p, T, i, X_liq, 0);
    }
    for (int i = 0; i < nsp; i++)
    {
        dZgasdXi[i] = 0;
        dZliqdXi[i] = 0;
        for (int j = 0; j < nsp; j++)
        {
            dZgasdXi[i] += dZdxi_gas_value[j] * dXdXi_gas_value[i][j];
            dZliqdXi[i] += dZdxi_liq_value[j] * dXdXi_liq_value[i][j];
        }
    }

    scalar dZMixturedT = sol.dZdT;
    scalar dZMixturedP = sol.dZdP;
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
    scalarList dEsdXi_value(nsp);
    for (int i = 0; i < nsp; i++)
    {
        dEsdXi_value[i] = dHsdXi_value[i] + p / sqr(rho_value) * drhodXi_value[i];
    }

    scalar dEsdT = Cp_Hs_value + p / sqr(rho_value) * drhodT_value;
    scalar dEsdP = dHsdP_value - 1 / rho_value + p / sqr(rho_value) * drhodP_value;
    for (int i = 0; i < nsp; i++)
    {
        grad[i][0] = -dEsdXi_value[i] / dEsdT;
        grad[i][1] = -drhodT_value * dEsdXi_value[i] / dEsdT + drhodXi_value[i];
        grad[i][2] = -dvfdT_value * dEsdXi_value[i] / dEsdT + dvfdXi_value[i];
        grad[i][3] = 0;
    }
    grad[nsp][0] = 1 / dEsdT;
    grad[nsp][1] = drhodT_value / dEsdT;
    grad[nsp][2] = dvfdT_value / dEsdT;
    grad[nsp][3] = 0;

    //Info << "drhodT=" << drhodT_value << endl;
    //Info << "drhodP=" << drhodP_value << endl;

    grad[nsp + 1][0] = -dEsdP / dEsdT;
    grad[nsp + 1][1] = -drhodT_value * dEsdP / dEsdT + dHsdP_value;
    grad[nsp + 1][2] = -dvfdT_value * dEsdP / dEsdT + dvfdP_value;
    grad[nsp + 1][3] = 0;

    return autoPtr<scalarRectangularMatrix>(&grad);
}
// ************************************************************************* //
