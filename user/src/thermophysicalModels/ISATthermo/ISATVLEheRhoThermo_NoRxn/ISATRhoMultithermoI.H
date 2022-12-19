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

#include "ISATRhoMultithermo.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template <class ThermoMixture, template <class> class Type>
inline Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::ISATRhoMultithermo(
    const ThermoMixture &sp)
    : ThermoMixture(sp)
{
}
template <class ThermoMixture, template <class> class Type>
inline Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::ISATRhoMultithermo(
    const word &name,
    PtrList<SingleThermoType> &speciesData,
    const speciesTable &specieNames,
    const dictionary &thermoDict,
    ISATRhoMultithermo<ThermoMixture, Type> *pin)
    : ThermoMixture(name, speciesData, specieNames, thermoDict),
      ISATDict_(
          IFstream(
              fileName(thermoDict.lookup("ISATDictFile")).expand())()),
      temp_p(pin),
      Erhovfcfun_XTP(temp_p ? temp_p : this),
      THvfcfun_XrhoP(temp_p ? temp_p : this),
      TPvfcfun_XErho(temp_p ? temp_p : this),
      Erhovfctree_XTP(speciesData.size() + 2, 4, Erhovfcfun_XTP, "Erhovfctree_XTP", ISATDict_),
      THvfctree_XrhoP(speciesData.size() + 2, 4, THvfcfun_XrhoP, "THvfctree_XrhoP", ISATDict_),
      TPvfctree_XErho(speciesData.size() + 2, 4, TPvfcfun_XErho, "TPvfctree_XErho", ISATDict_)
{
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::T(
    scalar f,
    scalar p,
    scalar T0,
    scalar (ISATRhoMultithermo<ThermoMixture, Type>::*F)(const scalar, const scalar) const,
    scalar (ISATRhoMultithermo<ThermoMixture, Type>::*dFdT)(const scalar, const scalar) const,
    scalar (ISATRhoMultithermo<ThermoMixture, Type>::*limit)(const scalar) const) const
{

    if (T0 < 0)
    {
        FatalErrorInFunction
            << "Negative initial temperature T0: " << T0
            << abort(FatalError);
    }
    //static int call_n=0;
    //static int inter_n=0;
    scalar T1 = T0;
    scalar F1 = (this->*F)(p, T1), F2 = F1;
    scalar T2 = (this->*limit)(T0 + 1);
    if (T2 == T1)
        T2 = (this->*limit)(T0 - 1);
    // if(F1>f)
    //     T2=(this->*limit)(T0-30);
    //scalar F1=(this->*F)(p, T1);
    //scalar F2=(this->*F)(p, T2);
    //scalar Tnew = T0;
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
        F2 = (this->*F)(p, T2);
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
                T2 = (this->*limit)(T1 - sign(relax * dF / dfdt) * 30);
            else
                T2 = (this->*limit)(T1 - relax * dF / dfdt);
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
                Fm = (this->*F)(p, Tm);
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
            } while (mag(T2 - T1) > Ttol && mag(f - Fm) > 0.0001);
            T2 = Tm;
        }
        //Tnew =
        //   (this->*limit)(Test - ((this->*F)(p, Test) - f) / (this->*dFdT)(p, Test));

        if (T2 == 0)
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

    } while (mag(T2 - T1) > Ttol && mag(f - F2) > 0.0001);
    return T2;
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::Tsol(
    scalar f,
    scalar p,
    scalar T0,
    scalar (ISATRhoMultithermo<ThermoMixture, Type>::*F)(const scalar, const scalar, typename ISATRhoMultithermo<ThermoMixture, Type>::solution &) const,
    scalar (ISATRhoMultithermo<ThermoMixture, Type>::*dFdT)(const scalar, const scalar)
        const,
    scalar (ISATRhoMultithermo<ThermoMixture, Type>::*limit)(const scalar) const,
    autoPtr<typename ISATRhoMultithermo<ThermoMixture, Type>::solution> &sol) const
{

    if (T0 < 0)
    {
        FatalErrorInFunction
            << "Negative initial temperature T0: " << T0
            << abort(FatalError);
    }
    //static int call_n=0;
    //static int inter_n=0;
    scalar T1 = T0;
    sol.reset(this->TPN(p, T1).ptr());
    scalar F1 = (this->*F)(p, T1, sol()), F2 = F1;
    scalar T2 = (this->*limit)(T0 + 1);
    if (T2 == T1)
        T2 = (this->*limit)(T0 - 1);
    scalar T2N = T2;
    // if(F1>f)
    //     T2=(this->*limit)(T0-30);
    //scalar F1=(this->*F)(p, T1);
    //scalar F2=(this->*F)(p, T2);
    //scalar Tnew = T0;
    scalar Ttol = T0 * tol_;
    scalar dF = 0, dfdt = 0, dFold = 1e100;
    int iter = 0;
    scalar relax = 1;
    bool inSection = false;
    bool bisec = false;
    double Tm, Fm;
    do
    {
        T2 = T2N;
        F1 = F2;
        sol.reset(this->TPN(p, T2).ptr());
        F2 = (this->*F)(p, T2, sol());
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
                T2N = (this->*limit)(T1 - sign(relax * dF / dfdt) * 30);
            else
                T2N = (this->*limit)(T1 - relax * dF / dfdt);
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
                sol.reset(this->TPN(p, Tm).ptr());
                Fm = (this->*F)(p, Tm, sol());
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
                                 << "Maximum number of iterations exceeded: " << maxIter_
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

    } while (mag(T2 - T1) > Ttol && mag(f - F2) > 1e-8 * f);
    return T2;
}

template <class ThermoMixture, template <class> class Type>
inline std::tuple<Foam::scalar, Foam::scalar> Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::TrhoEs(
    scalar e,
    scalar rho,
    scalar p0,
    scalar T0,
    scalar (ISATRhoMultithermo<ThermoMixture, Type>::*Fe)(const scalar, const scalar, typename ISATRhoMultithermo<ThermoMixture, Type>::solution &) const,
    scalar (ISATRhoMultithermo<ThermoMixture, Type>::*Frho)(const scalar, const scalar, typename ISATRhoMultithermo<ThermoMixture, Type>::solution &) const,
    autoPtr<typename ISATRhoMultithermo<ThermoMixture, Type>::solution> &retSol)
{
    if (T0 < 0)
    {
        FatalErrorInFunction
            << "Negative initial temperature T0: " << T0
            << abort(FatalError);
    }
    /*
        e = -152361.67898158100433647632598876953125;
        rho = 2.7025262888433871921733953058719635009765625;
        T0 = 385.130472898669040660024620592594146728515625;
        p0 = 339463.156064767041243612766265869140625;
        this->X_[0] = 0.99969244178928817756712987829814665019512176513672;
        this->X_[1] = 0.00030755821071176789750084568098031923000235110521317;
    */

    //static int call_n=0;
    //static int inter_n=0;

    scalar T1 = T0;
    scalar p1 = p0;
    auto sol(this->TPN(p1, T1));
    scalar e11 = (this->*Fe)(p1, T1, sol());     //, e21, e12;
    scalar rho11 = (this->*Frho)(p1, T1, sol()); //, rho21, rho12;
    //scalar T2 = T1 + 1, p2 = p1 + 10; //(this->*limit)(T0+1);
    //scalar Ttol = T0 * tol_, ptol = p0 * tol_;
    scalar de = 0, drho = 0; //, dedT = 0, dedp = 0, drhodT = 0, drhodp = 0;
    //scalar de_old = 1e10, drho_old = 1e10;
    //scalar dete;
    int iter = 0;
    scalar maxratio = 0.3;
    scalar dTdE, dPdE, dTdrho, dPdrho;
    //scalar dTdE2, dPdE2, dTdrho2, dPdrho2;
    //typename ISATRhoMultithermo<ThermoMixture, Type>::VLE_D1_data solnn_NIO(this->X_.size());
    //this->TPn_flash_New_TPD_Tudisco_NIO(p1, T1, this->X_, solnn_NIO);
    scalar minf = sqr(e11 - e) + sqr(rho11 - rho);
    scalar t_f;
    std::tie(dTdE, dPdE, dTdrho, dPdrho) = this->dTPdErho(p1, T1, sol());

    //double x1[2], y1[2], x2[2], y2[2], dx[2], dy[2];
    scalar dT = 1e-2;
    scalar dp = 1;
    /*     x1[0] = p0;
    x1[1] = T0;
    x2[0] = p0 + dp;
    x2[1] = T0 + dT;
    y1[0] = (this->*Frho)(x1[0], x1[1]);
    y1[1] = (this->*Fe)(x1[0], x1[1]);
    y2[0] = (this->*Frho)(x2[0], x2[1]);
    y2[1] = (this->*Fe)(x2[0], x2[1]);
    dx[0] = x2[0] - x1[0];
    dx[1] = x2[1] - x1[1];
    dy[0] = y2[0] - y1[0];
    dy[1] = y2[1] - y1[1]; */

    //auto sol(this->TPN(x1[0], x1[1]));
    //double H1[2][2], S[2][2], tempM[2][2], H2[2][2];
    //double dxTdy;
    //std::tie(dTdE, dPdE, dTdrho, dPdrho) = this->dTPdErho(x1[0], x1[1], sol());
    //std::tie(H1[1][1], H1[0][1], H1[1][0], H1[0][0]) = this->dTPdErho(x1[0], x1[1], sol());
    //double e10 = (this->*Fe)(x1[0] + dp, x1[1]);
    //double e01 = (this->*Fe)(x1[0] + dp, x1[1] + dT);
    //double rho10 = (this->*Frho)(x1[0] + dp, x1[1]);
    //double rho01 = (this->*Frho)(x1[0], x1[1] + dT);
    //dedT = (e01 - y1[0]) / dT;
    //drhodT = (rho01 - y1[1]) / dT;
    //dedp = (e10 - y1[0]) / dp;
    //drhodp = (rho10 - y1[1]) / dp;
    ////double ret = dedp * drhodT - dedT * drhodp;
    //H1[0][0] = dedp / ret;
    //H1[0][1] = -drhodp / ret;
    //H1[1][0] = -dedT / ret;
    //H1[1][1] = drhodT / ret;

    //H1[0][0] = 0.01;
    //H1[0][1] = 0;
    //H1[1][0] = 0;
    //H1[1][1] = 0.01;

    //H1[1][1] *= -y1[1]*y1[1];
    //H1[0][1] *= -y1[1]*y1[1];

    do
    {
        /*
        dp = max(mag(dp), 1e-4)*sign(dp);
        dT = max(mag(dT), 1e-5)*sign(dT);
        e12 = (this->*Fe)(p1, T1 + dT);
        e21 = (this->*Fe)(p1 + dp, T1);
        rho12 = (this->*Frho)(p1, T1 + dT);
        rho21 = (this->*Frho)(p1 + dp, T1);
        de_old=de;
        drho_old=drho;
        de = (e11 - e);
        drho = (rho11 - rho);

        if(mag(de)>mag(de_old)||mag(drho)>mag(drho_old))
        {
            relax*=0.5;
        }
        dedt = (e12 - e11) / dT;
        dedp = (e21 - e11) / dp;
        drhodt = (rho12 - rho11) / dT;
        drhodp = (rho21 - rho11) / dp;

        dete = dedt * drhodp - dedp * drhodt;

        dT = -(drhodp * de - dedp * drho) / dete * relax;
        dp = -(dedt * drho - drhodt * de) / dete * relax;
        T1 = T1 + dT;
        p1 = p1 + dp;



        e11 = (this->*Fe)(p1, T1);
        rho11 = (this->*Frho)(p1, T1);
        */
        /*
                dp = max(mag(dp), 1e-4) * sign(dp);
                dT = max(mag(dT), 1e-5) * sign(dT);
                e12 = (this->*Fe)(p1, T1 + dT);
                e21 = (this->*Fe)(p1 + dp, T1);
                rho12 = (this->*Frho)(p1, T1 + dT);
                rho21 = (this->*Frho)(p1 + dp, T1);
                de_old = de;
                drho_old = drho;
                de = (e11 - e);
                drho = (rho11 - rho);

                dedt = (e12 - e11) / dT;
                dedp = (e21 - e11) / dp;
                drhodt = (rho12 - rho11) / dT;
                drhodp = (rho21 - rho11) / dp;
                dete = dedt * drhodp - dedp * drhodt;



                dT = -(drhodp * de - dedp * drho) / dete * relax;
                dp = -(dedt * drho - drhodt * de) / dete * relax;
                T1 = T1 + dT;
                p1 = p1 + dp;
        */

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
        //this->TPn_flash_New_TPD_Tudisco_NIO(p1, T1, this->X_, solnn_NIO);
        retSol.reset(this->TPN(p1, T1).ptr());
        //for (int i = 0; i < this->X_.size(); i++)
        //{
        //    solnn_NIO.Liq.X[i] = retSol().X_liq()[i];
        //    solnn_NIO.Gas.X[i] = retSol().X_gas()[i];
        //}
        //solnn_NIO.vaporFraction = retSol().vaporfra;
        /* if (mag(retSol().vaporfra - solnn_NIO.vaporFraction) > 1e-8)
        {
            Info << retSol().vaporfra << "," << solnn_NIO.vaporFraction << endl;

            Info << endl;
            NotImplemented;
        } */
        e11 = (this->*Fe)(p1, T1, retSol());
        rho11 = (this->*Frho)(p1, T1, retSol());
        t_f = sqr(e11 - e) + sqr(rho11 - rho);

        if (t_f > minf)
        {
            //maxratio *= 0.9;
            if (maxratio > 1e-5)
            {
                maxratio *= 0.9;
            }
        }

        std::tie(dTdE, dPdE, dTdrho, dPdrho) = this->dTPdErho(p1, T1, retSol());
        /*         std::tie(dTdE2, dPdE2, dTdrho2, dPdrho2) = this->dTPdErho_NIO(p1, T1, solnn_NIO);
        if (mag(dTdE - dTdE2) + mag(dPdE - dPdE2) > 1e-8)
        {
            Info << setprecision(50);
            Info << this->X_ << endl;
            Info << retSol().X_gas() << endl;
            Info << solnn_NIO.Gas.X << endl;
            Info << retSol().X_liq() << endl;
            Info << solnn_NIO.Liq.X << endl;
            Info << retSol().vaporfra << endl;
            Info << solnn_NIO.vaporFraction << endl;
            Info << dTdE << "," << dTdE2 << endl;
            Info << dPdE << "," << dPdE2 << endl;
            Info << endl;
            NotImplemented;
        } */
        if (t_f < minf)
        {
            minf = t_f;
        }

        //minf = t_f;

        /*
        dxTdy = dx[0] * dy[0] + dx[1] * dy[1];
        S[0][0] = 1 - dy[0] * dx[0] / dxTdy;
        S[0][1] = -dy[0] * dx[1] / dxTdy;
        S[1][0] = -dy[1] * dx[0] / dxTdy;
        S[1][1] = 1 - dy[1] * dx[1] / dxTdy;

        tempM[0][0] = H1[0][0] * S[0][0] + H1[0][1] * S[1][0];
        tempM[0][1] = H1[0][0] * S[0][1] + H1[0][1] * S[1][1];
        tempM[1][0] = H1[1][0] * S[0][0] + H1[1][1] * S[1][0];
        tempM[1][1] = H1[1][0] * S[0][1] + H1[1][1] * S[1][1];

        H2[0][0] = S[0][0] * tempM[0][0] + S[1][0] * tempM[1][0] + dx[0] * dx[0] / dxTdy;
        H2[0][1] = S[0][0] * tempM[0][1] + S[1][0] * tempM[1][1] + dx[0] * dx[1] / dxTdy;
        H2[1][0] = S[0][1] * tempM[0][0] + S[1][1] * tempM[1][0] + dx[1] * dx[0] / dxTdy;
        H2[1][1] = S[0][1] * tempM[0][1] + S[1][1] * tempM[1][1] + dx[1] * dx[1] / dxTdy;

        x1[0] = x2[0];
        x1[1] = x2[1];

        y1[0] = y2[0];
        y1[1] = y2[1];
        dx[0] = -H2[0][0] * (y1[0] - rho) - H2[0][1] * (y1[1] - e);
        dx[1] = -H2[1][0] * (y1[0] - rho) - H2[1][1] * (y1[1] - e);

        if (mag(dx[0]) > x1[0] * 0.0001)
        {
            dx[0] = dx[0] / mag(dx[0]) * x1[0] * 0.0001;
            dx[1] = dx[1] / mag(dx[0]) * x1[0] * 0.0001;
        }

        if (mag(dx[1]) > x1[1] * 0.0001)
        {
            dx[0] = dx[0] / mag(dx[1]) * x1[1] * 0.0001;
            dx[1] = dx[1] / mag(dx[1]) * x1[1] * 0.0001;
        }

        x2[0] = x1[0] + dx[0];
        x2[1] = x1[1] + dx[1];

        //auto sol(this->TPN(x2[0], x2[1]));

        y2[0] = (this->*Frho)(x2[0], x2[1]);
        y2[1] = (this->*Fe)(x2[0], x2[1]);

        dy[0] = y2[0] - y1[0];
        dy[1] = y2[1] - y1[1];

        H1[0][0] = H2[0][0];
        H1[0][1] = H2[0][1];
        H1[1][0] = H2[1][0];
        H1[1][1] = H2[1][1];
        */
        if (iter++ > 5 * maxIter_)
        {
            WarningInFunction << setprecision(50)
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
            //<< abort(FatalError);
            //continue;
            break;
        }
    } while ((mag(e11 - e) > 1e-1 || mag(rho11 - rho) > 1e-4)); //while ((mag(y2[1] - e) > 1e-2 || mag(y2[0] - rho) > 1e-7)); //while ((mag(e11 - e) > 1e-3 || mag(rho11 - rho) > 1e-5)); // && (mag(dT) > 1e-5 || mag(dp) > 1e-4));

    //return std::make_tuple(x2[0], x2[1]);
    return std::make_tuple(p1, T1);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class ThermoMixture, template <class> class Type>
inline Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::ISATRhoMultithermo(
    const word &name,
    const ISATRhoMultithermo &st)
    : ThermoMixture(name, st)
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class ThermoMixture, template <class> class Type>
inline Foam::word
Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::heName()
{
    return Type<ISATRhoMultithermo<ThermoMixture, Type>>::energyName();
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::Cv(const scalar p, const scalar T) const
{
    autoPtr<typename ISATRhoMultithermo<ThermoMixture, Type>::solution> sol(this->TPN(p, T));
    //return this->Cp(p, T, sol()) - this->CpMCv(p, T, sol());
    return this->Cv_VLE(p, T, sol());
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::Cpv(const scalar p, const scalar T) const
{
    return Type<ISATRhoMultithermo<ThermoMixture, Type>>::Cpv(*this, p, T);
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::gamma(const scalar p, const scalar T) const
{
    const scalar Cp = this->Cp(p, T);
    return Cp / (Cp - this->CpMCv(p, T));
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::CpByCpv(
    const scalar p,
    const scalar T) const
{
    return Type<ISATRhoMultithermo<ThermoMixture, Type>>::CpByCpv(*this, p, T);
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::HE(const scalar p, const scalar T) const
{
    return Type<ISATRhoMultithermo<ThermoMixture, Type>>::HE(*this, p, T);
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::Es(const scalar p, const scalar T) const
{
    /*    scalarList in(this->X().size() + 2), out(this->X().size() + 5);
    for (int i = 0; i < this->X().size(); i++)
    {
        in[i] = this->X_[i];
    }
    in[this->X().size()] = T;
    in[this->X().size() + 1] = p;

    //Info << "psi_HP" << endl;
    autoPtr<typename ISATRhoMultithermo<ThermoMixture, Type>::solution> sol(nullptr);
    Erhovfctree_XTP.call(in, out, sol);
    //psivftree.call(in, out);
    return out[0]; //std::make_tuple(out[0], out[1], out[2], out[3]); */
    return this->Hs(p, T) - p / this->rho(p, T);
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::Es_ac(const scalar p, const scalar T) const
{
    return this->Hs(p, T) - p / this->rho(p, T);
}
template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::THP(const scalar H, const scalar P, const scalar T0) const
{
    return this->T(H, P, T0, &ISATRhoMultithermo<ThermoMixture, Type>::Hs, &ISATRhoMultithermo<ThermoMixture, Type>::Cp, &ISATRhoMultithermo<ThermoMixture, Type>::limit);
}

template <class ThermoMixture, template <class> class Type>
inline std::tuple<Foam::scalar, Foam::scalar, Foam::scalar>
Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::rhovfc_ISAT(const scalar p, const scalar T) const
{
    scalarList in(this->X().size() + 2), out(this->X().size() + 5);
    for (int i = 0; i < this->X().size(); i++)
    {
        in[i] = this->X_[i];
    }
    in[this->X().size()] = T;
    in[this->X().size() + 1] = p;

    //Info << "psi_HP" << endl;
    autoPtr<typename ISATRhoMultithermo<ThermoMixture, Type>::solution> sol(nullptr);
    Erhovfctree_XTP.call(in, out, sol);
    //psivftree.call(in, out);
    //return out[1];//std::make_tuple(out[0], out[1], out[2], out[3]);
    return std::make_tuple(out[1], out[2], out[3]);
    //return this->Hs(p, T) - p / this->rho(p, T);
}

template <class ThermoMixture, template <class> class Type>
inline std::tuple<Foam::scalar, Foam::scalar, Foam::scalar, Foam::scalar>
Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::Erhovfc_ISAT(const scalar p, const scalar T) const
{
    scalarList in(this->X().size() + 2), out(this->X().size() + 5);
    for (int i = 0; i < this->X().size(); i++)
    {
        in[i] = this->X_[i];
    }
    in[this->X().size()] = T;
    in[this->X().size() + 1] = p;

    //Info << "psi_HP" << endl;
    autoPtr<typename ISATRhoMultithermo<ThermoMixture, Type>::solution> sol(nullptr);
    Erhovfctree_XTP.call(in, out, sol);
    //psivftree.call(in, out);
    //return out[1];//std::make_tuple(out[0], out[1], out[2], out[3]);
    return std::make_tuple(out[0], out[1], out[2], out[3]);
    //return this->Hs(p, T) - p / this->rho(p, T);
}

template <class ThermoMixture, template <class> class Type>
inline std::tuple<Foam::scalar, Foam::scalar, Foam::scalar, Foam::scalar>
Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::Erhovfc_XPT(const scalar p, const scalar T) const
{
    scalarList in(this->X().size() + 2), out(4);
    for (int i = 0; i < this->X().size(); i++)
    {
        in[i] = this->X_[i];
    }
    in[this->X().size()] = T;
    in[this->X().size() + 1] = p;

    //autoPtr<typename ISATRhoMultithermo<ThermoMixture, Type>::solution> sol(nullptr);
    typename ISATRhoMultithermo<ThermoMixture, Type>::VLE_D1_data sol_NIO(this->X().size());
    Erhovfctree_XTP.call(in, out, sol_NIO);

    return std::make_tuple(out[0], out[1], out[2], out[3]);
}

template <class ThermoMixture, template <class> class Type>
inline std::tuple<Foam::scalar, Foam::scalar, Foam::scalar, Foam::scalar, Foam::scalar>
Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::Erhovfc_G_rhoY_ISAT(const scalar p, const scalar T, scalarList &Y) const
{
    scalarList in(this->X().size() + 2), out(this->X().size() + 5);
    for (int i = 0; i < this->X().size(); i++)
    {
        in[i] = this->X_[i];
    }
    in[this->X().size()] = T;
    in[this->X().size() + 1] = p;

    //Info << "psi_HP" << endl;
    autoPtr<typename ISATRhoMultithermo<ThermoMixture, Type>::solution> sol(nullptr);
    Erhovfctree_XTP.call(in, out, sol);
    for (int i = 0; i < this->X().size(); i++)
    {
        Y[i] = out[5 + i];
    }
    //psivftree.call(in, out);
    //return out[1];//std::make_tuple(out[0], out[1], out[2], out[3]);
    return std::make_tuple(out[0], out[1], out[2], out[3], out[4]);
    //return this->Hs(p, T) - p / this->rho(p, T);
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::Ea(const scalar p, const scalar T) const
{
    return this->Ha(p, T) - p / this->rho(p, T);
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::G(const scalar p, const scalar T) const
{
    return this->Ha(p, T) - T * this->S(p, T);
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::A(const scalar p, const scalar T) const
{
    return this->Ea(p, T) - T * this->S(p, T);
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::cp(const scalar p, const scalar T) const
{
    return this->Cp(p, T) * this->W();
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::ha(const scalar p, const scalar T) const
{
    return this->Ha(p, T) * this->W();
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::hs(const scalar p, const scalar T) const
{
    return this->Hs(p, T) * this->W();
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::hc() const
{
    return this->Hc() * this->W();
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::s(const scalar p, const scalar T) const
{
    return this->S(p, T) * this->W();
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::he(const scalar p, const scalar T) const
{
    return this->HE(p, T) * this->W();
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::cv(const scalar p, const scalar T) const
{
    return this->Cv(p, T) * this->W();
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::es(const scalar p, const scalar T) const
{
    return this->Es(p, T) * this->W();
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::ea(const scalar p, const scalar T) const
{
    return this->Ea(p, T) * this->W();
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::g(const scalar p, const scalar T) const
{
    return this->G(p, T) * this->W();
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::a(const scalar p, const scalar T) const
{
    return this->A(p, T) * this->W();
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::K(const scalar p, const scalar T) const
{
    scalar arg = -1; //this->Y()*this->G(Pstd, T)/(RR*T);

    if (arg < 600)
    {
        return exp(arg);
    }
    else
    {
        return rootVGreat;
    }
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::Kp(const scalar p, const scalar T) const
{
    return K(p, T);
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::Kc(const scalar p, const scalar T) const
{
    const scalar nm = 1; //this->Y()/this->W();

    if (equal(nm, small))
    {
        return Kp(p, T);
    }
    else
    {
        return Kp(p, T) * pow(Pstd / (RR * T), nm);
    }
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::Kx(
    const scalar p,
    const scalar T) const
{
    const scalar nm = this->Y() / this->W();

    if (equal(nm, small))
    {
        return Kp(p, T);
    }
    else
    {
        return Kp(p, T) * pow(Pstd / p, nm);
    }
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::Kn(
    const scalar p,
    const scalar T,
    const scalar n) const
{
    const scalar nm = this->Y() / this->W();

    if (equal(nm, small))
    {
        return Kp(p, T);
    }
    else
    {
        return Kp(p, T) * pow(n * Pstd / p, nm);
    }
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::THE(
    const scalar he,
    const scalar p,
    const scalar T0) const
{
    return Type<ISATRhoMultithermo<ThermoMixture, Type>>::THE(*this, he, p, T0);
}
//This one !!!

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::THs(
    const scalar hs,
    const scalar p,
    const scalar T0) const
{

    double sum = 0;
    //scalarList in(this->X().size() + 2), out(1);
    scalarList in(this->X().size() + 1), out(1);
    for (int i = 0; i < this->X().size() - 1; i++)
    {
        in[i] = this->X_[i];
        sum += this->X_[i];
    }
    in[this->X().size() - 1] = hs;
    in[this->X().size()] = p;
    //in[this->X().size() + 1] = T0;

    if (1 - sum < 1e-5)
    {
        //Ttree.notCall++;
        return T(hs, p, T0, &ISATRhoMultithermo<ThermoMixture, Type>::Hs, &ISATRhoMultithermo<ThermoMixture, Type>::Cp, &ISATRhoMultithermo<ThermoMixture, Type>::limit);
    }
    //Hnfun.value(in, out);
    //Info << "THS!!" << endl;
    return T(hs, p, T0, &ISATRhoMultithermo<ThermoMixture, Type>::Hs, &ISATRhoMultithermo<ThermoMixture, Type>::Cp, &ISATRhoMultithermo<ThermoMixture, Type>::limit);
    //Ttree.call(in, out);
    //Ttree.call(in, out, T0);
    //return out[0];
}

template <class ThermoMixture, template <class> class Type>
inline std::tuple<Foam::scalar, Foam::scalar> Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::TrhoE(
    const scalar e,
    const scalar rho,
    const scalar p0,
    const scalar T0)
{
    return TrhoEs(e, rho, p0, T0, &ISATRhoMultithermo<ThermoMixture, Type>::Es_ac, &ISATRhoMultithermo<ThermoMixture, Type>::rho);
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::THa(
    const scalar ha,
    const scalar p,
    const scalar T0) const
{
    return T(
        ha,
        p,
        T0,
        &ISATRhoMultithermo<ThermoMixture, Type>::Ha,
        &ISATRhoMultithermo<ThermoMixture, Type>::Cp,
        &ISATRhoMultithermo<ThermoMixture, Type>::limit);
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::TEs(
    const scalar es,
    const scalar p,
    const scalar T0) const
{
    return T(
        es,
        p,
        T0,
        &ISATRhoMultithermo<ThermoMixture, Type>::Es,
        &ISATRhoMultithermo<ThermoMixture, Type>::Cv,
        &ISATRhoMultithermo<ThermoMixture, Type>::limit);
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::TEa(
    const scalar ea,
    const scalar p,
    const scalar T0) const
{
    return T(
        ea,
        p,
        T0,
        &ISATRhoMultithermo<ThermoMixture, Type>::Ea,
        &ISATRhoMultithermo<ThermoMixture, Type>::Cv,
        &ISATRhoMultithermo<ThermoMixture, Type>::limit);
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::dKcdTbyKc(
    const scalar p,
    const scalar T) const
{
    const scalar nm = 1; //this->Y()/this->W();

    if (equal(nm, small))
    {
        return -this->dGdT(Pstd, T); //*this->Y()/RR;
    }
    else
    {
        return -(nm / T + this->dGdT(Pstd, T)); //*this->Y()/RR);
    }
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::dcpdT(const scalar p, const scalar T) const
{
    return this->dCpdT(p, T) * this->W();
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::psi(const scalar p, const scalar T) const
{

    double sum = 0;
    scalarList in(this->X().size() + 1), out(1);
    for (int i = 0; i < this->X().size() - 1; i++)
    {
        in[i] = this->X_[i];
        sum += this->X_[i];
    }
    in[this->X().size() - 1] = p;
    in[this->X().size()] = T;

    if (1 - sum < 1e-5)
    {
        //psitree.notCall++;
        return this->ThermoMixture::psi(p, T);
    }
    // Info<<"psihaha"<<endl;

    //psitree.call(in, out);
    //psivftree.call(in, out);
    return this->ThermoMixture::psi(p, T);
    //return out[0];
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::psi_a(const scalar h, const scalar p, const scalar T0) const
{
    scalar tt = T(h, p, T0, &ISATRhoMultithermo<ThermoMixture, Type>::Hs, &ISATRhoMultithermo<ThermoMixture, Type>::Cp, &ISATRhoMultithermo<ThermoMixture, Type>::limit);
    return this->ThermoMixture::psi(p, tt);
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::psi_HP(const scalar H, const scalar P, const scalar T0) const
{
    scalarList in(this->X().size() + 2), out(1);
    for (int i = 0; i < this->X().size(); i++)
    {
        in[i] = this->X_[i];
    }
    in[this->X().size()] = H;
    in[this->X().size() + 1] = P;

    //Info << "psi_HP" << endl;

    //psitree_HP.call(in, out, T0);
    //psivftree.call(in, out);
    return out[0];
}

template <class ThermoMixture, template <class> class Type>
inline std::tuple<Foam::scalar, Foam::scalar>
Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::Tpsi_HP(const scalar H, const scalar P, const scalar T0) const
{
    scalarList in(this->X().size() + 2), out(2);
    for (int i = 0; i < this->X().size(); i++)
    {
        in[i] = this->X_[i];
    }
    in[this->X().size()] = H;
    in[this->X().size() + 1] = P;

    //Info << "psi_HP" << endl;

    //Tpsitree_HP.call(in, out, T0);
    //psivftree.call(in, out);
    return std::make_tuple(out[0], out[1]);
}
/*
template <class ThermoMixture, template <class> class Type>
inline std::tuple<Foam::scalar, Foam::scalar, Foam::scalar,Foam::scalar>
Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::Tpsivfc_XHP(const scalar H, const scalar P, const scalar T0) const
{
    scalarList in(this->X().size() + 2), out(3);
    for (int i = 0;i < this->X().size();i++)
    {
        in[i] = this->X_[i];
    }
    in[this->X().size()] = H;
    in[this->X().size() + 1] = P;

    //Info << "psi_HP" << endl;

    Tpsivfctree_XHP.call(in, out, T0);
    //psivftree.call(in, out);
    return std::make_tuple(out[0], out[1], out[2], out[3]);
}
*/
template <class ThermoMixture, template <class> class Type>
inline std::tuple<Foam::scalar, Foam::scalar, Foam::scalar, Foam::scalar>
Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::THvfc_XrhoP(const scalar rho, const scalar P, const scalar T0) const
{
    scalarList in(this->X().size() + 2), out(4);
    for (int i = 0; i < this->X().size(); i++)
    {
        in[i] = this->X_[i];
    }
    in[this->X().size()] = rho;
    in[this->X().size() + 1] = P;

    //Info << "psi_HP" << endl;
    //autoPtr<typename ISATRhoMultithermo<ThermoMixture, Type>::solution> sol(nullptr);
    typename ISATRhoMultithermo<ThermoMixture, Type>::VLE_D1_data sol_NIO(this->X().size());
    THvfctree_XrhoP.call(in, out, T0, sol_NIO);
    //psivftree.call(in, out);
    return std::make_tuple(out[0], out[1], out[2], out[3]);
}
template <class ThermoMixture, template <class> class Type>
inline std::tuple<Foam::scalar, Foam::scalar, Foam::scalar, Foam::scalar, Foam::scalar>
Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::THvfc_G_rhoY_XrhoP(const scalar rho, const scalar P, const scalar T0, scalarList &Y) const
{
    scalarList in(this->X().size() + 2), out(this->X().size() + 5);
    for (int i = 0; i < this->X().size(); i++)
    {
        in[i] = this->X_[i];
    }
    in[this->X().size()] = rho;
    in[this->X().size() + 1] = P;

    //Info << "psi_HP" << endl;
    autoPtr<typename ISATRhoMultithermo<ThermoMixture, Type>::solution> sol(nullptr);

    THvfctree_XrhoP.call(in, out, T0, sol);
    for (int i = 0; i < this->X().size(); i++)
    {
        Y[i] = out[5 + i];
    }
    return std::make_tuple(out[0], out[1], out[2], out[3], out[4]);
}
/*
template <class ThermoMixture, template <class> class Type>
inline std::tuple<Foam::scalar, Foam::scalar, Foam::scalar, Foam::scalar>
Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::TPvfc_Erho(const scalar E, const scalar rho, const scalar T0, const scalar P0) const
{
    scalarList in(this->X().size() + 2), out(3);
    for (int i = 0;i < this->X().size();i++)
    {
        in[i] = this->X_[i];
    }
    in[this->X().size()] = E;
    in[this->X().size() + 1] = rho;

    //Info << "psi_HP" << endl;

    TPvfctree_XErho.call(in, out, T0, P0);
    //psivftree.call(in, out);
    return std::make_tuple(out[0], out[1], out[2], out[3]);
}
*/

template <class ThermoMixture, template <class> class Type>
inline std::tuple<Foam::scalar, Foam::scalar, Foam::scalar, Foam::scalar>
Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::TPvfc_XErho(const scalar E, const scalar rho, const scalar T0, const scalar P0) const
{
    scalarList in(this->X().size() + 2), out(4);
    for (int i = 0; i < this->X().size(); i++)
    {
        in[i] = this->X_[i];
    }
    in[this->X().size()] = E;
    in[this->X().size() + 1] = rho;

    //autoPtr<typename ISATRhoMultithermo<ThermoMixture, Type>::solution> sol(nullptr);
    typename ISATRhoMultithermo<ThermoMixture, Type>::VLE_D1_data sol_NIO(this->X().size());
    TPvfctree_XErho.call(in, out, T0, P0, sol_NIO);

    return std::make_tuple(out[0], out[1], out[2], out[3]);
}
/*
template <class ThermoMixture, template <class> class Type>
inline std::tuple<Foam::scalar, Foam::scalar, Foam::scalar, Foam::scalar, Foam::autoPtr<typename ThermoMixture::solution>>
Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::TPvfcsol_Erho(const scalar E, const scalar rho, const scalar T0, const scalar P0) const
{
    scalarList in(this->X().size() + 2), out(3);
    for (int i = 0;i < this->X().size();i++)
    {
        in[i] = this->X_[i];
    }
    in[this->X().size()] = E;
    in[this->X().size() + 1] = rho;

    //Info << "psi_HP" << endl;

    TPvfcsoltree_XErho.call(in, out, T0, P0);
    autoPtr<typename ThermoMixture::solution> pS(new typename thermoType::solution());
    pS().X_liq.reset(new scalarList(this->X().size()));
    pS().X_gas.reset(new scalarList(this->X().size()));
    pS().vaporfra=out[2];
    for(int i=0;i<this->X().size();i++)
    {
        pS().X_gas()[i]=out[4+i];
        pS().X_liq()[i]=out[4+this->X().size()+i];

    }
    //psivftree.call(in, out);
    return std::make_tuple(out[0], out[1], out[2], out[3],pS);
}
*/
/*
template <class ThermoMixture, template <class> class Type>
inline std::tuple<Foam::scalar, Foam::scalar, Foam::scalar, Foam::scalar, Foam::autoPtr<typename ThermoMixture::solution>>
Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::Tpsivfcsol_XHP(const scalar H, const scalar P, const scalar T0) const
{
     scalarList in(this->X().size() + 2), out(3);
    for (int i = 0;i < this->X().size();i++)
    {
        in[i] = this->X_[i];
    }
    in[this->X().size()] = H;
    in[this->X().size() + 1] = P;

    //Info << "psi_HP" << endl;

    Tpsivfcsoltree_XHP.call(in, out, T0);
    autoPtr<typename ThermoMixture::solution> pS(new typename thermoType::solution());
    pS().X_liq.reset(new scalarList(this->X().size()));
    pS().X_gas.reset(new scalarList(this->X().size()));
    pS().vaporfra=out[2];
    for(int i=0;i<this->X().size();i++)
    {
        pS().X_gas()[i]=out[4+i];
        pS().X_liq()[i]=out[4+this->X().size()+i];

    }
    //psivftree.call(in, out);
    return std::make_tuple(out[0], out[1], out[2], out[3],pS);
}
*/

template <class ThermoMixture, template <class> class Type>
inline void Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::newTimeStep()
{
    //TPvfctree_XErho.newTimeStep();
    Erhovfctree_XTP.newTimeStep();
    //Tpsivfctree_XHP.newTimeStep();
    //Tpsivfcsoltree_XHP.newTimeStep();
    THvfctree_XrhoP.newTimeStep();
    TPvfctree_XErho.newTimeStep();
    //kappatree_XTP.newTimeStep();
    //mutree_XTP.newTimeStep();
    //alphahtree_XTP.newTimeStep();
}

template <class ThermoMixture, template <class> class Type>
inline bool Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::newLoop()
{
    return THvfctree_XrhoP.newLoop() || TPvfctree_XErho.newLoop();
}

template <class ThermoMixture, template <class> class Type>
inline void Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::mute_show()
{
    Erhovfctree_XTP.showPerformance();
    Erhovfctree_XTP.mute();
}
/*
template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::mu(const scalar p, const scalar T) const
{

    //double sum = 0;
    scalarList in(this->X().size() + 2), out(1);
    for (int i = 0;i < this->X().size();i++)
    {
        in[i] = this->X_[i];
        //sum += this->X_[i];
    }
    in[this->X().size()] = p;
    in[this->X().size() + 1] = T;

        // Info<<"psihaha"<<endl;

    mutree_XTP.call(in, out);
    return out[0];
}*/
/*
template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::alphah(const scalar p, const scalar T) const
{

    double sum = 0;
    scalarList in(this->X().size() + 2), out(1);
    for (int i = 0;i < this->X().size();i++)
    {
        in[i] = this->X_[i];
        // sum += this->X_[i];
    }
    in[this->X().size()] = p;
    in[this->X().size() + 1] = T;

        // Info<<"psihaha"<<endl;

    alphahtree_XTP.call(in, out);
    return out[0];

}*/
/*
template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::kappa(const scalar p, const scalar T) const
{

    //double sum = 0;
    scalarList in(this->X().size() + 2), out(1);
    for (int i = 0;i < this->X().size();i++)
    {
        in[i] = this->X_[i];
        //sum += this->X_[i];
    }
    in[this->X().size()] = p;
    in[this->X().size() + 1] = T;


    // Info<<"psihaha"<<endl;

    kappatree_XTP.call(in, out);
    return out[0];

}*/

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::vaporfra(const scalar p, const scalar T) const
{

    double sum = 0;
    scalarList in(this->X().size() + 1), out(1);
    for (int i = 0; i < this->X().size() - 1; i++)
    {
        in[i] = this->X_[i];
        sum += this->X_[i];
    }
    in[this->X().size() - 1] = p;
    in[this->X().size()] = T;

    if (1 - sum < 1e-5)
    {
        //vaporfratree.notCall++;
        return this->ThermoMixture::vaporfra(p, T);
    }
    // Info<<"psihaha"<<endl;

    //vaporfratree.call(in, out);
    //psivftree.call(in, out);
    //Info<<in<<endl;
    return this->ThermoMixture::vaporfra(p, T);
    return out[0];
}
template <class ThermoMixture, template <class> class Type>
inline void Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::psi_vf(const scalar p, const scalar T, scalar &psi_out, scalar &vf_out) const
{
    double sum = 0;
    scalarList in(this->X().size() + 1), out(1);
    for (int i = 0; i < this->X().size() - 1; i++)
    {
        in[i] = this->X_[i];
        sum += this->X_[i];
    }
    in[this->X().size() - 1] = p;
    in[this->X().size()] = T;

    if (1 - sum < 1e-5)
    {
        //psivftree.notCall++;
        psi_out = this->ThermoMixture::psi(p, T);
        vf_out = this->ThermoMixture::vaporfra(p, T);
    }

    // Info<<"psihaha"<<endl;

    //psivftree.call(in, out);
    psi_out = this->ThermoMixture::psi(p, T);
    vf_out = this->ThermoMixture::vaporfra(p, T);
    psi_out = out[0];
    vf_out = out[1];
    return;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template <class ThermoMixture, template <class> class Type>
inline void Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::operator+=(
    const ISATRhoMultithermo<ThermoMixture, Type> &st)
{
    ThermoMixture::operator+=(st);
}

template <class ThermoMixture, template <class> class Type>
inline void Foam::species::ISATRhoMultithermo<ThermoMixture, Type>::operator*=(const scalar s)
{
    ThermoMixture::operator*=(s);
}

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template <class ThermoMixture, template <class> class Type>
inline Foam::species::ISATRhoMultithermo<ThermoMixture, Type> Foam::species::operator+(
    const ISATRhoMultithermo<ThermoMixture, Type> &st1,
    const ISATRhoMultithermo<ThermoMixture, Type> &st2)
{
    return ISATRhoMultithermo<ThermoMixture, Type>(
        static_cast<const ThermoMixture &>(st1) + static_cast<const ThermoMixture &>(st2));
}

template <class ThermoMixture, template <class> class Type>
inline Foam::species::ISATRhoMultithermo<ThermoMixture, Type> Foam::species::operator*(
    const scalar s,
    const ISATRhoMultithermo<ThermoMixture, Type> &st)
{
    return ISATRhoMultithermo<ThermoMixture, Type>(
        s * static_cast<const ThermoMixture &>(st));
}

template <class ThermoMixture, template <class> class Type>
inline Foam::species::ISATRhoMultithermo<ThermoMixture, Type> Foam::species::operator==(
    const ISATRhoMultithermo<ThermoMixture, Type> &st1,
    const ISATRhoMultithermo<ThermoMixture, Type> &st2)
{
    return ISATRhoMultithermo<ThermoMixture, Type>(
        static_cast<const ThermoMixture &>(st1) == static_cast<const ThermoMixture &>(st2));
}

// ************************************************************************* //
