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

Class
    Foam::species::ISATRhoMultithermo

Description
    Basic thermodynamics type based on the use of fitting functions for
    cp, h, s obtained from the template argument type ISATRhoMultithermo.  All other
    properties are derived from these primitive functions.

SourceFiles
    ISATRhoMultithermoI.H
    ISATRhoMultithermo.C

\*---------------------------------------------------------------------------*/

#ifndef ISATRhoMultithermo_H
#define ISATRhoMultithermo_H

#include "thermodynamicConstants.H"
#include "ISATmanager.H"
#include <tuple>
using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    //TODO finish this!!!!
/*
    struct ISATfunction
    {
        ISATmanager* p;
        typedef double (ISATmanager::* VLEfunc)(double, double);
        VLEfunc pfunc;
        ISATfunction(ISATmanager p_in, VLEfunc pfunc_in) :p(p_in), pfunc(pfunc_in) {}
        void  value(const scalarList& in, scalarList& out)
        {
            out.resize(1);
            out[0] = (*p).*pfunc(in);//+2000*::Foam::sin(in[1]/10);
        }
        void  derive(const scalarList& in, scalarRectangularMatrix& out)
        {
            scalarList temp = in;
            scalarList tempv;
            value(in, tempv);
            scalarList tempv2;
            scalarList dd;
            scalar dx = 1e-5;
            for (int i = 0;i < in.size();i++)
            {
                temp[i] += dx;
                value(temp, tempv2);
                dd = (tempv2 - tempv) / dx;
                for (int j = 0;j < tempv2.size();j++)
                {
                    out[j][i] = dd[j];
                }
                temp[i] -= dx;
            }
        }
    };
    */


    namespace species
    {

        // Forward declaration of friend functions and operators

        template<class ThermoMixture, template<class> class Type> class ISATRhoMultithermo;

        template<class ThermoMixture, template<class> class Type>
        inline ISATRhoMultithermo<ThermoMixture, Type> operator+
            (
                const ISATRhoMultithermo<ThermoMixture, Type>&,
                const ISATRhoMultithermo<ThermoMixture, Type>&
                );

        template<class ThermoMixture, template<class> class Type>
        inline ISATRhoMultithermo<ThermoMixture, Type> operator*
            (
                const scalar,
                const ISATRhoMultithermo<ThermoMixture, Type>&
                );

        template<class ThermoMixture, template<class> class Type>
        inline ISATRhoMultithermo<ThermoMixture, Type> operator==
            (
                const ISATRhoMultithermo<ThermoMixture, Type>&,
                const ISATRhoMultithermo<ThermoMixture, Type>&
                );


        template<class ThermoMixture, template<class> class Type>
        Ostream& operator<<
            (
                Ostream&,
                const ISATRhoMultithermo<ThermoMixture, Type>&
                );

        template<class Type>
        struct Hnflashfunction
        {
            Type* p;
            Hnflashfunction(Type* pin) :p(pin) {}
            void value(const scalarList& in, scalarList& out, scalar T0) const   //ToDo ADd t0 as sep agru
            {
                scalar H = in[in.size() - 2];
                scalar P = in[in.size() - 1];
                scalarList temp(p->X_.size());
                out.resize(1);
                double sum = 0;
                int n_sp = p->X_.size();
                for (int i = 0;i < n_sp - 1;i++)
                {
                    temp[i] = p->X_[i];
                    p->X_[i] = in[i];
                    sum += in[i];
                }
                temp[n_sp - 1] = p->X_[n_sp - 1];
                p->X_[n_sp - 1] = 1 - sum;
                out[0] = p->T(H, P, T0, &Type::Hs, &Type::Cp, &Type::limit);
                //Info << out << endl;
                for (int i = 0;i < n_sp;i++)
                {
                    p->X_[i] = temp[i];
                }

            }

            void  derive(const scalarList& in, scalarRectangularMatrix& out, scalar T0)
            {
                scalarList temp = in;
                scalarList tempv;
                value(in, tempv, T0);
                scalarList tempv2;
                scalarList dd;
                double sum = 0;
                for (int i = 0;i < in.size() - 2;i++)
                    sum += in[i];
                scalar dx = 1e-6;
                for (int i = 0;i < in.size();i++)
                {
                    temp[i] += dx;
                    value(temp, tempv2, T0);
                    dd = (tempv2 - tempv) / dx;
                    for (int j = 0;j < tempv2.size();j++)
                    {
                        out[i][j] = dd[j];
                    }
                    temp[i] -= dx;
                }
            }
        };

        template<class Type>
        struct thermofunction
        {
            Type* p;
            typedef scalar(Type::* Tfunc) (const scalar p, const scalar T, autoPtr<typename Type::solution>& sol) const;
            Tfunc pfunc;
            thermofunction(Type* pin, Tfunc fin) :p(pin), pfunc(fin) {}
            void value(const scalarList& in, scalarList& out) const
            {
                scalarList temp(p->X_.size());
                out.resize(1);
                double sum = 0;
                int n_sp = p->X_.size();
                for (int i = 0;i < n_sp - 1;i++)
                {
                    temp[i] = p->X_[i];
                    p->X_[i] = in[i];
                    sum += in[i];
                }
                temp[n_sp - 1] = p->X_[n_sp - 1];
                p->X_[n_sp - 1] = 1 - sum;
                autoPtr<typename Type::solution> sol(p->TPN(in[in.size() - 2], in[in.size() - 1]));
                out[0] = (p->*pfunc)(in[in.size() - 2], in[in.size() - 1], sol);
                //Info << out << endl;
                for (int i = 0;i < n_sp;i++)
                {
                    p->X_[i] = temp[i];
                }

            }

            void  derive(const scalarList& in, scalarRectangularMatrix& out)
            {
                scalarList temp = in;
                scalarList tempv;
                value(in, tempv);
                scalarList tempv2;
                scalarList dd;
                double sum = 0;
                for (int i = 0;i < in.size() - 2;i++)
                    sum += in[i];
                scalar dx = 1e-6;
                for (int i = 0;i < in.size();i++)
                {
                    temp[i] += dx;
                    value(temp, tempv2);
                    dd = (tempv2 - tempv) / dx;
                    for (int j = 0;j < tempv2.size();j++)
                    {
                        out[i][j] = dd[j];
                    }
                    temp[i] -= dx;
                }
            }
        };
        template<class Type>
        struct thermofunction_psi
        {
            Type* p;
            typedef scalar(Type::* Tfunc) (const scalar p, const scalar T, autoPtr<typename Type::solution>& sol) const;
            thermofunction_psi(Type* pin) :p(pin) {}
            void value(const scalarList& in, scalarList& out, scalar T0) const
            {
                scalar H = in[in.size() - 2];
                scalar P = in[in.size() - 1];
                scalar T = p->T(H, P, T0, &Type::Hs, &Type::Cp, &Type::limit);
                //scalarList temp(p->X_.size());
                out.resize(1);
                /*
                for (int i = 0;i < p->X_.size();i++)
                {
                    temp[i] = in[i];
                }
                */
                p->setX(in);
                autoPtr<typename Type::solution> sol(p->TPN(P, T));

                out[0] = (p->Type::thermoVLE::psi)(P, T, sol);
            }

            void  derive(const scalarList& in, scalarRectangularMatrix& out, scalar T0)
            {
                //static double maxx = 0;
                scalar H = in[in.size() - 2];
                scalar P = in[in.size() - 1];
                scalar T = p->T(H, P, T0, &Type::Hs, &Type::Cp, &Type::limit);
                /*
                scalarList temp(p->X_.size());
                int n_sp = p->X_.size();
                for (int i = 0;i < n_sp;i++)
                {
                    temp[i] = in[i];
                }
                */
                p->setX(in);
                autoPtr<typename Type::solution> sol(p->TPN(P, T));
                autoPtr<scalarList> grad(p->drhoPdXHP_HsP(P, T, sol()));
                for (int i = 0;i < grad->size();i++)
                {
                    out[i][0] = grad()[i];
                }
                /*
                for (int i = 0;i < p->X_.size();i++)
                {
                    out[i][0] = p->drhoPdXi_HsP(P, T, i, sol());
                }
                out[p->X_.size()][0] = p->drhoPdH_HsP(P, T, sol());
                out[p->X_.size() + 1][0] = p->drhoPdP_HsP(P, T, sol());
                */
                //Info << grad() << endl;
                //Info << out << endl;
                /*
                for (int i = 0;i < 4;i++)
                {
                    if (fabs((out[i][0] - grad()[i]) / grad()[i]) > maxx)
                        maxx = fabs((out[i][0] - grad()[i]) / grad()[i]);
                    //Info << (out[i][0] - grad()[i]) / grad()[i] << ",";
                }
                if (maxx > 1e-06)
                {
                    Info << grad() << endl;
                    Info << out << endl;
                    Info << maxx << endl;
                }*/
                //Info << maxx << endl;
                //FatalErrorInFunction << "vaporf 1" << exit(FatalError);
            }
        };
        template<class Type>
        struct thermofunction_Tpsi
        {
            Type* p;
            thermofunction_Tpsi(Type* pin) :p(pin) {}
            void value(const scalarList& in, scalarList& out, scalar T0) const
            {
                scalar H = in[in.size() - 2];
                scalar P = in[in.size() - 1];
                scalar T = p->T(H, P, T0, &Type::Hs, &Type::Cp, &Type::limit);
                //scalarList temp(p->X_.size());
                out.resize(2);
                /*
                for (int i = 0;i < p->X_.size();i++)
                {
                    temp[i] = in[i];
                }
                */
                p->setX(in);
                autoPtr<typename Type::solution> sol(p->TPN(P, T));
                out[0] = T;
                out[1] = (p->Type::thermoVLE::psi)(P, T, sol);
            }

            void  derive(const scalarList& in, scalarRectangularMatrix& out, scalar T0)
            {
                scalar H = in[in.size() - 2];
                scalar P = in[in.size() - 1];
                scalar T = p->T(H, P, T0, &Type::Hs, &Type::Cp, &Type::limit);

                p->setX(in);
                autoPtr<typename Type::solution> sol(p->TPN(P, T));
                autoPtr<scalarRectangularMatrix> grad(p->dTrhoPdXHP_HsP(P, T, sol()));
                out = grad();
            }
        };

        template<class Type>
        struct thermofunction_Tpsivf
        {
            Type* p;
            thermofunction_Tpsivf(Type* pin) :p(pin) {}
            void value(const scalarList& in, scalarList& out, scalar T0) const
            {
                scalar H = in[in.size() - 2];
                scalar P = in[in.size() - 1];
                scalar T = p->T(H, P, T0, &Type::Hs, &Type::Cp, &Type::limit);
                //scalarList temp(p->X_.size());
                out.resize(3);
                /*
                for (int i = 0;i < p->X_.size();i++)
                {
                    temp[i] = in[i];
                }
                */
                p->setX(in);
                autoPtr<typename Type::solution> sol(p->TPN(P, T));
                out[0] = T;
                out[1] = (p->Type::thermoVLE::psi)(P, T, sol);
                out[2] = sol().vaporfra;
            }

            void  derive(const scalarList& in, scalarRectangularMatrix& out, scalar T0)
            {
                scalar H = in[in.size() - 2];
                scalar P = in[in.size() - 1];
                scalar T = p->T(H, P, T0, &Type::Hs, &Type::Cp, &Type::limit);

                p->setX(in);
                autoPtr<typename Type::solution> sol(p->TPN(P, T));//sol(p->TPn_flash_New_TPD(P, T));//sol(p->TPn_flash(P, T));
                autoPtr<scalarRectangularMatrix> grad(p->dTrhoPvfdXHP_HsP(P, T, sol()));
                out = grad();
            }
        };

        template<class Type>
        struct thermofunction_TPvf
        {
            Type* p;
            thermofunction_TPvf(Type* pin) :p(pin) {}
            void value(const scalarList& in, scalarList& out, scalar T0, scalar P0) const
            {
                scalar E = in[in.size() - 2];
                scalar rho = in[in.size() - 1];
                scalar T, P;
                p->setX(in);
                std::tie(P, T) = p->TrhoE(E, rho, P0, T0);
                autoPtr<typename Type::solution> sol(p->TPN(P, T));
/*
                scalarList test1(2), test2(2);
                scalar e1, e2, rho1, rho2, p1, p2, T1, T2;
                e1 = -115670.06779247121;
                rho1 = 209.08049449382986;
                test1[0] = 3.2509615069070284e-05;
                test1[1] = 0.99996749038493093;
                p->setX(test1);
                std::tie(p1, T1) = p->TrhoE(e1, rho1, P0, T0);

                e2 = -345151.48554157966;
                rho2 = 674.32780349253721;

                test2[0] = 0.99969244178928895;
                test2[1] = 0.0003075582107111563;


                p->setX(test2);
                std::tie(p2, T2) = p->TrhoE(e2, rho2, P0, T0);

                scalar rhomix, emix, al = 0.99, bt = 1 - al, pmix, Tmix;
                scalarList testmix(2);
                rhomix = al * rho1 + bt * rho2;
                emix = (al * rho1 * e1 + bt * rho2 * e2) / rhomix;
                testmix[0] = (al * rho1 * 0.0001 / 86.1754 + bt * rho2 * 0.9999 / 86.1754) / (al * rho1 * 0.0001 / 86.1754 + bt * rho2 * 0.9999 / 86.1754 + al * rho1 * 0.9999 / 28.0134 + bt * rho2 * 0.0001 / 28.0134);
                testmix[1] = 1 - testmix[0];
                p->setX(testmix);
                std::tie(pmix, Tmix) = p->TrhoE(emix, rhomix, P0, T0);
                FatalErrorInFunction
                    << "T1: " << T1
                    << "\nT2: " << T2
                    << "\np1: " << p1
                    << "\np2: " << p2
                    << "\nTmix: " << Tmix
                    << "\npmix: " << pmix
                    << exit(FatalError);
                autoPtr<typename Type::solution> sol1(p->TPN(P, T));
                scalar cc = p->c(P, T, sol1());
                */
                /*
                                autoPtr<typename Type::solution> sol1(p->TPn_flash(P, T));
                                autoPtr<scalarRectangularMatrix> grad(p->dTPvfdXErho_constC(P, T, sol1()));
                                scalar dx = 1e-4;
                                scalarList test(2);
                                test[0] = in[0];
                                test[1] = in[1]+dx;
                                scalar T1, P1;
                                p->setX(test);
                                std::tie(P1, T1) = p->TrhoE(E, rho, P0, T0);
                                autoPtr<typename Type::solution> sol2(p->TPn_flash(P1, T1));
                */

                out.resize(4);

                out[0] = T;
                out[1] = P;
                out[2] = sol().vaporfra;
                out[3] = p->c(P, T, sol());
                /*
                                FatalErrorInFunction
                                    << "T1: " << P
                                    << "\nT2: " << P1
                                    << "\ndEsdx1: " << (P1-P) / dx
                                    << "\ndEsdx1_ AA: " << grad()[1][1]
                                    << exit(FatalError);
                                    */

            }

            void  derive(const scalarList& in, scalarRectangularMatrix& out, scalar T0, scalar P0)
            {
                scalar E = in[in.size() - 2];
                scalar rho = in[in.size() - 1];
                scalar T, P;
                std::tie(P, T) = p->TrhoE(E, rho, P0, T0);

                p->setX(in);
                autoPtr<typename Type::solution> sol(p->TPN(P, T));//sol(p->TPn_flash_New_TPD(P, T));//sol(p->TPn_flash(P, T));
                autoPtr<scalarRectangularMatrix> grad(p->dTPvfdXErho_constC(P, T, sol()));
                //autoPtr<scalarRectangularMatrix> grad(p->dTPvfdXErho(P, T, sol()));
                out = grad();
            }
            bool valid_in(const scalarList& in)
            {
                bool flag = true;
                for (int i = 0;i < in.size() - 2;i++)
                {
                    flag = flag && (in[i] > 0);
                }
                flag = flag && (in[in.size() - 1] > 0);
                return flag;
            }
        };

        template<class Type>
        struct thermofunction_TPX_to_Erhovfc
        {
            Type* p;
            thermofunction_TPX_to_Erhovfc(Type* pin) :p(pin) {}
            void value(const scalarList& in, scalarList& out) const
            {
                scalar T = in[in.size() - 2];
                scalar P = in[in.size() - 1];
                //scalar T, P;
                /*scalarList test(2);
                test[0] = 3.2647302995662537e-05;
                test[1] = 0.99996735269700421;
                T = 310.99984251699999;
                P = 19999979.704999998;
                p->setX(test);
                autoPtr<typename Type::solution> sol1(p->TPN(P, T));
                scalar E1 = p->Hs(P, T, sol1()) - P / p->rho(P, T, sol1());
                scalar dd= p->dEdXi(P,T,0,sol1());
                autoPtr<scalarRectangularMatrix> grad(p->dErhovfcdXTP_simp(P, T, sol1()));
                scalar dx = 1e-4;
                test[0] += dx;
                p->setX(test);
                autoPtr<typename Type::solution> sol2(p->TPN(P, T));
                scalar E2 = p->Hs(P, T, sol2()) - P / p->rho(P, T, sol2());
                */
                /*
                test[0] = 5.881479846724812e-05;//3.2647302995662537e-05;
                test[1] = 0.99994118520153286;//0.99996735269700421;
                T = 310.97695248799999;//310.99984251699999;
                P = 19996764.623;//19999979.704999998;
                */
                //std::tie(P, T) = p->TrhoE(E, rho, P0, T0);
                //scalarList temp(p->X_.size());
                out.resize(4);
                //out.resize(3);
                /*
                for (int i = 0;i < p->X_.size();i++)
                {
                    temp[i] = in[i];
                }
                */
                p->setX(in);

                autoPtr<typename Type::solution> sol(p->TPN(P, T));
                out[0] = p->Hs(P, T, sol()) - P / p->rho(P, T, sol());
                out[1] = p->rho(P, T, sol());
                out[2] = sol().vaporfra;
                out[3] = p->c(P, T, sol());
                /*
                FatalErrorInFunction
                    << "Es1: " << E1
                    << "\nEs2: " << E2
                    << "\ndEsdx1: " << (E2-E1)/dx
                    << "\ndEsdx1_ AA: " << dd
                    << exit(FatalError);
                    */
            }

            void  derive(const scalarList& in, scalarRectangularMatrix& out)
            {
                scalar T = in[in.size() - 2];
                scalar P = in[in.size() - 1];
                //scalar T, P;
                //std::tie(P, T) = p->TrhoE(E, rho, P0, T0);

                p->setX(in);
                autoPtr<typename Type::solution> sol(p->TPN(P, T));//sol(p->TPn_flash_New_TPD(P, T));//sol(p->TPn_flash(P, T));
                autoPtr<scalarRectangularMatrix> grad(p->dErhovfcdXTP_simp(P, T, sol()));
                //autoPtr<scalarRectangularMatrix> grad(p->dTPvfdXErho(P, T, sol()));
                out = grad();
            }
            bool valid_in(const scalarList& in)
            {
                bool flag = true;
                for (int i = 0;i < in.size();i++)
                {
                    flag = flag && (in[i] > 0);
                }
                return flag;
            }
        };
        /*
                template<class Type>
                struct thermofunction_TPvf
                {
                    Type* p;
                    thermofunction_TPvf(Type* pin) :p(pin) {}
                    void value(const scalarList& in, scalarList& out, scalar T0,scalar P0) const
                    {
                        scalar E = in[in.size() - 2];
                        scalar rho = in[in.size() - 1];
                        scalar P,T;
                        std::tie(P,T)=p->TrhoEs(E, rho, P0, T0, &Type::Es, &Type::rho);
                        //scalar T = p->T(H, P, T0, &Type::Hs, &Type::Cp, &Type::limit);
                        //scalarList temp(p->X_.size());
                        out.resize(3);

                        p->setX(in);
                        autoPtr<typename Type::solution> sol(p->TPn_flash(P, T));
                        out[0] = T;
                        out[1] = P;
                        out[2] = sol().vaporfra;
                    }

                    void  derive(const scalarList& in, scalarRectangularMatrix& out, scalar T0,scalar P0)
                    {
                        scalar E = in[in.size() - 2];
                        scalar rho = in[in.size() - 1];
                        //Todo finish these
                        //scalar H = in[in.size() - 2];
                        //scalar P = in[in.size() - 1];
                        //scalar T = p->T(H, P, T0, &Type::Hs, &Type::Cp, &Type::limit);
                        scalar P,T;
                        std::tie(P,T)=p->TrhoEs(E, rho, P0, T0, &Type::Es, &Type::rho);
                        p->setX(in);
                        autoPtr<typename Type::solution> sol(p->TPn_flash(P, T));//sol(p->TPn_flash_New_TPD(P, T));//sol(p->TPn_flash(P, T));
                        //Todo finish these
                        autoPtr<scalarRectangularMatrix> grad(p->dTrhoPvfdXHP_HsP(P, T, sol()));
                        out = grad();
                    }
                };
        */


        template<class Type>
        struct thermofunction2
        {
            Type* p;
            typedef scalar(Type::* Tfunc) (const scalar p, const scalar T) const;
            Tfunc pfunc;
            thermofunction2(Type* pin, Tfunc fin) :p(pin), pfunc(fin) {}
            void value(const scalarList& in, scalarList& out) const
            {
                scalarList temp(p->X_.size());
                out.resize(1);
                double sum = 0;
                int n_sp = p->X_.size();
                for (int i = 0;i < n_sp - 1;i++)
                {
                    temp[i] = p->X_[i];
                    p->X_[i] = in[i];
                    sum += in[i];
                }
                temp[n_sp - 1] = p->X_[n_sp - 1];
                p->X_[n_sp - 1] = 1 - sum;
                out[0] = (p->*pfunc)(in[in.size() - 2], in[in.size() - 1]);
                //Info << out << endl;
                for (int i = 0;i < n_sp;i++)
                {
                    p->X_[i] = temp[i];
                }

            }

            void  derive(const scalarList& in, scalarRectangularMatrix& out)
            {
                scalarList temp = in;
                scalarList tempv;
                value(in, tempv);
                scalarList tempv2;
                scalarList dd;
                double sum = 0;
                for (int i = 0;i < in.size() - 2;i++)
                    sum += in[i];
                scalar dx = 1e-6;
                for (int i = 0;i < in.size();i++)
                {
                    temp[i] += dx;
                    value(temp, tempv2);
                    dd = (tempv2 - tempv) / dx;
                    for (int j = 0;j < tempv2.size();j++)
                    {
                        out[i][j] = dd[j];
                    }
                    temp[i] -= dx;
                }
            }
        };

        template<class Type>
        struct thermofunction3
        {
            Type* p;
            typedef scalar(Type::* Tfunc) (const scalar p, const scalar T) const;
            Tfunc pfunc;
            thermofunction3(Type* pin, Tfunc fin) :p(pin), pfunc(fin) {}
            void value(const scalarList& in, scalarList& out) const
            {
                scalar T = in[in.size() - 2];
                scalar P = in[in.size() - 1];
                scalarList temp(p->X_.size());
                out.resize(1);
                p->setX(in);
                out[0] = (p->*pfunc)(in[in.size() - 2], in[in.size() - 1]);
            }

            void  derive(const scalarList& in, scalarRectangularMatrix& out)
            {
                scalarList temp = in;
                scalarList tempv;
                value(in, tempv);
                scalarList tempv2;
                scalarList dd;
                scalar dx = 1e-5;
                for (int i = 0;i < in.size();i++)
                {
                    temp[i] += dx;
                    value(temp, tempv2);
                    dd = (tempv2 - tempv) / dx;
                    for (int j = 0;j < tempv2.size();j++)
                    {
                        out[i][j] = dd[j];
                    }
                    temp[i] -= dx;
                }
            }
            bool valid_in(const scalarList& in)
            {
                bool flag = true;
                for (int i = 0;i < in.size();i++)
                {
                    flag = flag && (in[i] > 0);
                }
                return flag;
            }
        };

        template<class Type>
        struct psi_vffunction
        {
            Type* p;
            typedef scalar(Type::* Tfunc) (const scalar p, const scalar T, autoPtr<typename Type::solution>& sol) const;
            Tfunc ppsi, pvf;
            psi_vffunction(Type* pin, Tfunc ppsi_in, Tfunc pvf_in) :p(pin), ppsi(ppsi_in), pvf(pvf_in) {}
            void value(const scalarList& in, scalarList& out) const
            {
                scalarList temp(p->X_.size());
                out.resize(2);
                double sum = 0;
                int n_sp = p->X_.size();
                for (int i = 0;i < n_sp - 1;i++)
                {
                    temp[i] = p->X_[i];
                    p->X_[i] = in[i];
                    sum += in[i];
                }
                temp[n_sp - 1] = p->X_[n_sp - 1];
                p->X_[n_sp - 1] = 1 - sum;
                autoPtr<typename Type::solution> sol(p->TPN(in[in.size() - 2], in[in.size() - 1]));
                out[0] = (p->*ppsi)(in[in.size() - 2], in[in.size() - 1], sol);
                out[1] = (p->*pvf)(in[in.size() - 2], in[in.size() - 1], sol);
                //Info << out << endl;
                for (int i = 0;i < n_sp;i++)
                {
                    p->X_[i] = temp[i];
                }

            }

            void  derive(const scalarList& in, scalarRectangularMatrix& out)
            {
                scalarList temp = in;
                scalarList tempv;
                value(in, tempv);
                scalarList tempv2;
                scalarList dd;
                double sum = 0;
                for (int i = 0;i < in.size() - 2;i++)
                    sum += in[i];
                scalar dx = 1e-6;
                for (int i = 0;i < in.size();i++)
                {
                    temp[i] += dx;
                    value(temp, tempv2);
                    dd = (tempv2 - tempv) / dx;
                    for (int j = 0;j < tempv2.size();j++)
                    {
                        out[i][j] = dd[j];
                    }
                    temp[i] -= dx;
                }
            }
        };


        /*---------------------------------------------------------------------------*\
                                   Class ISATRhoMultithermo Declaration
        \*---------------------------------------------------------------------------*/

        template<class ThermoMixture, template<class> class Type>
        class ISATRhoMultithermo
            :
            public ThermoMixture,
            public Type<ISATRhoMultithermo<ThermoMixture, Type>>
        {
            // Private data

                //- Convergence tolerance of energy -> temperature inversion functions
            static const scalar tol_;

            //- Max number of iterations in energy->temperature inversion functions
            static const int maxIter_;
            dictionary ISATDict_;


            ISATRhoMultithermo<ThermoMixture, Type>* temp_p;
            thermofunction_TPvf<ISATRhoMultithermo<ThermoMixture, Type>> TPvffun_XErho;
            thermofunction_TPX_to_Erhovfc<ISATRhoMultithermo<ThermoMixture, Type>>  Erhovfcfun_XTP;
            //thermofunction2<ThermoMixture> HEfun;
            thermofunction3<ThermoMixture> kappafun_XTP;
            thermofunction3<ThermoMixture> mufun_XTP;
            thermofunction3<ThermoMixture> alphahfun_XTP;



            mutable ISATmanager< thermofunction_TPvf<ISATRhoMultithermo<ThermoMixture, Type>>> TPvftree_XErho;
            mutable ISATmanager< thermofunction_TPX_to_Erhovfc<ISATRhoMultithermo<ThermoMixture, Type>>> Erhovfctree_XTP;

            //mutable ISATmanager< thermofunction2<ThermoMixture>> HEtree;
            mutable ISATmanager< thermofunction3<ThermoMixture>> kappatree_XTP;
            mutable ISATmanager< thermofunction3<ThermoMixture>> mutree_XTP;
            mutable ISATmanager< thermofunction3<ThermoMixture>> alphahtree_XTP;





            // Private Member Functions

                //- Return the temperature corresponding to the value of the
                //  thermodynamic property f, given the function f = F(p, T)
                //  and dF(p, T)/dT
            inline scalar T
            (
                scalar f,
                scalar p,
                scalar T0,
                scalar(ISATRhoMultithermo::* F)(const scalar, const scalar) const,
                scalar(ISATRhoMultithermo::* dFdT)(const scalar, const scalar) const,
                scalar(ISATRhoMultithermo::* limit)(const scalar) const
            ) const;
            /*
                    inline scalar T_IQI
                    (
                        scalar f,
                        scalar p,
                        scalar T0,
                        scalar (ISATRhoMultithermo::*F)(const scalar, const scalar) const,
                        scalar (ISATRhoMultithermo::*dFdT)(const scalar, const scalar) const,
                        scalar (ISATRhoMultithermo::*limit)(const scalar) const
                    ) const;
            */
            inline std::tuple<scalar, scalar> TrhoEs
            (
                scalar e,
                scalar rho,
                scalar p0,
                scalar T0,
                scalar(ISATRhoMultithermo<ThermoMixture, Type>::* Fe)(const scalar, const scalar) const,
                scalar(ISATRhoMultithermo<ThermoMixture, Type>::* Frho)(const scalar, const scalar) const
            );

        public:

            //- The thermodynamics of the individual species'
            typedef ThermoMixture thermoVLE;
            typedef ISATRhoMultithermo<ThermoMixture, Type> thermoType;
            typedef typename ThermoMixture::SingleThermoType SingleThermoType;

            // Constructors

                //- Construct from components
            inline ISATRhoMultithermo(const ThermoMixture& sp);
            //- Construct from dictionary
            ISATRhoMultithermo(const dictionary& dict, PtrList<ThermoMixture>& speciesData);

            //- Construct as named copy
            inline ISATRhoMultithermo(const word& name, const ISATRhoMultithermo&);


            inline ISATRhoMultithermo(const word& name, PtrList<SingleThermoType>&, const speciesTable&, const dictionary&, ISATRhoMultithermo<ThermoMixture, Type>* pin);


            // Member Functions

                //- Return the instantiated type name
            static word typeName()
            {
                return ThermoMixture::typeName();
            }

            //- Name of Enthalpy/Internal energy
            static inline word heName();

            inline void newTimeStep();
            // Fundamental properties
            // (These functions must be provided in derived types)

                // Heat capacity at constant pressure [J/(kg K)]
                // inline scalar Cp(const scalar p, const scalar T) const;

                // Sensible enthalpy [J/kg]
                // inline scalar Hs(const scalar p, const scalar T) const;

                // Chemical enthalpy [J/kg]
                // inline scalar Hc() const;

                // Absolute Enthalpy [J/kg]
                // inline scalar Ha(const scalar p, const scalar T) const;


                // Entropy [J/(kg K)]
                // inline scalar S(const scalar p, const scalar T) const;
            inline scalar psi(const scalar p, const scalar T) const;
            inline scalar psi_a(const scalar p, const scalar T, const scalar T0) const;
            inline scalar psi_HP(const scalar p, const scalar T, const scalar T0) const;
            inline std::tuple<scalar, scalar> Tpsi_HP(const scalar p, const scalar T, const scalar T0) const;
            inline std::tuple<scalar, scalar, scalar> Tpsivf_HP(const scalar p, const scalar T, const scalar T0) const;
            inline std::tuple<scalar, scalar, scalar, scalar > TPvf_Erho(const scalar E, const scalar rho, const scalar T0, const scalar P0) const;

            inline scalar mu(const scalar p, const scalar T) const;
            inline scalar alphah(const scalar p, const scalar T) const;
            inline scalar kappa(const scalar p, const scalar T) const;
            inline scalar vaporfra(const scalar p, const scalar T) const;
            inline void psi_vf(const scalar p, const scalar T, scalar& psi_out, scalar& vf_out) const;
            //inline scalar psi(const scalar p, const scalar T) const;

            // Mass specific derived properties

                //- Heat capacity at constant volume [J/(kg K)]
            inline scalar Cv(const scalar p, const scalar T) const;

            //- Heat capacity at constant pressure/volume [J/(kg K)]
            inline scalar Cpv(const scalar p, const scalar T) const;

            //- Gamma = Cp/Cv []
            inline scalar gamma(const scalar p, const scalar T) const;

            //- Ratio of heat capacity at constant pressure to that at
            //  constant pressure/volume []
            inline scalar CpByCpv(const scalar p, const scalar T) const;

            //- Enthalpy/Internal energy [J/kg]
            inline scalar HE(const scalar p, const scalar T) const;

            //- Sensible internal energy [J/kg]
            inline scalar Es(const scalar p, const scalar T) const;
            inline scalar Es_ac(const scalar p, const scalar T) const;

            inline std::tuple<scalar, scalar, scalar> rhovfc_ISAT(const scalar p, const scalar T) const;

            //- Absolute internal energy [J/kg]
            inline scalar Ea(const scalar p, const scalar T) const;

            //- Gibbs free energy [J/kg]
            inline scalar G(const scalar p, const scalar T) const;

            //- Helmholtz free energy [J/kg]
            inline scalar A(const scalar p, const scalar T) const;


            // Mole specific derived properties

                //- Heat capacity at constant pressure [J/(kmol K)]
            inline scalar cp(const scalar p, const scalar T) const;

            //- Absolute Enthalpy [J/kmol]
            inline scalar ha(const scalar p, const scalar T) const;

            //- Sensible enthalpy [J/kmol]
            inline scalar hs(const scalar p, const scalar T) const;

            //- Chemical enthalpy [J/kmol]
            inline scalar hc() const;

            //- Entropy [J/(kmol K)]
            inline scalar s(const scalar p, const scalar T) const;

            //- Enthalpy/Internal energy [J/kmol]
            inline scalar he(const scalar p, const scalar T) const;

            //- Heat capacity at constant volume [J/(kmol K)]
            inline scalar cv(const scalar p, const scalar T) const;

            //- Sensible internal energy [J/kmol]
            inline scalar es(const scalar p, const scalar T) const;

            //- Absolute internal energy [J/kmol]
            inline scalar ea(const scalar p, const scalar T) const;

            //- Gibbs free energy [J/kmol]
            inline scalar g(const scalar p, const scalar T) const;

            //- Helmholtz free energy [J/kmol]
            inline scalar a(const scalar p, const scalar T) const;


            // Equilibrium reaction thermodynamics

                //- Equilibrium constant [] i.t.o fugacities
                //  = PIi(fi/Pstd)^nui
            inline scalar K(const scalar p, const scalar T) const;

            //- Equilibrium constant [] i.t.o. partial pressures
            //  = PIi(pi/Pstd)^nui
            //  For low pressures (where the gas mixture is near perfect) Kp = K
            inline scalar Kp(const scalar p, const scalar T) const;

            //- Equilibrium constant i.t.o. molar concentration
            //  = PIi(ci/cstd)^nui
            //  For low pressures (where the gas mixture is near perfect)
            //  Kc = Kp(pstd/(RR*T))^nu
            inline scalar Kc(const scalar p, const scalar T) const;

            //- Equilibrium constant [] i.t.o. mole-fractions
            //  For low pressures (where the gas mixture is near perfect)
            //  Kx = Kp(pstd/p)^nui
            inline scalar Kx
            (
                const scalar p,
                const scalar T
            ) const;

            //- Equilibrium constant [] i.t.o. number of moles
            //  For low pressures (where the gas mixture is near perfect)
            //  Kn = Kp(n*pstd/p)^nui where n = number of moles in mixture
            inline scalar Kn
            (
                const scalar p,
                const scalar T,
                const scalar n
            ) const;


            // Energy->temperature  inversion functions

                //- Temperature from enthalpy or internal energy
                //  given an initial temperature T0
            inline scalar THE
            (
                const scalar H,
                const scalar p,
                const scalar T0
            ) const;

            //- Temperature from sensible enthalpy given an initial T0
            inline scalar THs
            (
                const scalar Hs,
                const scalar p,
                const scalar T0
            ) const;

            //- Temperature from absolute enthalpy
            //  given an initial temperature T0
            inline scalar THa
            (
                const scalar H,
                const scalar p,
                const scalar T0
            ) const;

            //- Temperature from sensible internal energy
            //  given an initial temperature T0
            inline scalar TEs
            (
                const scalar E,
                const scalar p,
                const scalar T0
            ) const;

            //- Temperature from absolute internal energy
            //  given an initial temperature T0
            inline scalar TEa
            (
                const scalar E,
                const scalar p,
                const scalar T0
            ) const;

            inline std::tuple<scalar, scalar> TrhoE
            (
                const scalar e,
                const scalar rho,
                const scalar p0,
                const scalar T0
            );


            // Derivative term used for Jacobian

                //- Derivative of B (acooding to Niemeyer et al.) w.r.t. temperature
            inline scalar dKcdTbyKc(const scalar p, const scalar T) const;

            //- Derivative of cp w.r.t. temperature
            inline scalar dcpdT(const scalar p, const scalar T) const;


            // I-O

                //- Write to Ostream
            void write(Ostream& os) const;


            // Member operators

            inline void operator+=(const ISATRhoMultithermo&);
            inline void operator*=(const scalar);


            // Friend operators

            friend ISATRhoMultithermo operator+ <ThermoMixture, Type>
                (
                    const ISATRhoMultithermo&,
                    const ISATRhoMultithermo&
                    );

            friend ISATRhoMultithermo operator* <ThermoMixture, Type>
                (
                    const scalar s,
                    const ISATRhoMultithermo&
                    );

            friend ISATRhoMultithermo operator== <ThermoMixture, Type>
                (
                    const ISATRhoMultithermo&,
                    const ISATRhoMultithermo&
                    );



            // Ostream Operator

            friend Ostream& operator<< <ThermoMixture, Type>
                (
                    Ostream&,
                    const ISATRhoMultithermo&
                    );

            template<class TypeName>
            friend struct Hnflashfunction;

            template<class TypeName>
            friend struct thermofunction_psi;
            template<class TypeName>
            friend struct thermofunction_Tpsi;
            template<class TypeName>
            friend struct thermofunction_Tpsivf;
        };




        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    } // End namespace species
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ISATRhoMultithermoI_dev.H"

#ifdef NoRepository
#include "ISATRhoMultithermo_dev.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
