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

Description

SourceFiles
    ISATFunction.H

\*---------------------------------------------------------------------------*/

#ifndef ISATFunction_H
#define ISATFunction_H

namespace Foam
{
    namespace species
    {
        struct basicFunction
        {
            static void read(const dictionary &ISATDict_, scalarList &maxLeafsize_, scalarList &toleranceOut_in, scalarList &initToleranceIn_in, scalarList &scaleIn_in)
            {
                maxLeafsize_ = scalarList(ISATDict_.lookup("maxLeafsize"));
                toleranceOut_in = scalarList(ISATDict_.lookup("toleranceOut"));
                initToleranceIn_in = scalarList(ISATDict_.lookup("initToleranceIn"));
                scaleIn_in = scalarList(ISATDict_.lookup("scaleIn"));
            }
        };
        template <class Type>
        struct Hnflashfunction
        {
            Type *p;
            Hnflashfunction(Type *pin) : p(pin) {}
            void value(const scalarList &in, scalarList &out, scalar T0) const // ToDo ADd t0 as sep agru
            {
                scalar H = in[in.size() - 2];
                scalar P = in[in.size() - 1];
                scalarList temp(p->X_.size());
                out.resize(1);
                double sum = 0;
                int n_sp = p->X_.size();
                for (int i = 0; i < n_sp - 1; i++)
                {
                    temp[i] = p->X_[i];
                    p->X_[i] = in[i];
                    sum += in[i];
                }
                temp[n_sp - 1] = p->X_[n_sp - 1];
                p->X_[n_sp - 1] = 1 - sum;
                out[0] = p->T(H, P, T0, &Type::Hs, &Type::Cp, &Type::limit);
                // Info << out << endl;
                for (int i = 0; i < n_sp; i++)
                {
                    p->X_[i] = temp[i];
                }
            }

            void derive(const scalarList &in, scalarRectangularMatrix &out, scalar T0)
            {
                scalarList temp = in;
                scalarList tempv;
                value(in, tempv, T0);
                scalarList tempv2;
                scalarList dd;
                double sum = 0;
                for (int i = 0; i < in.size() - 2; i++)
                    sum += in[i];
                scalar dx = 1e-6;
                for (int i = 0; i < in.size(); i++)
                {
                    temp[i] += dx;
                    value(temp, tempv2, T0);
                    dd = (tempv2 - tempv) / dx;
                    for (int j = 0; j < tempv2.size(); j++)
                    {
                        out[i][j] = dd[j];
                    }
                    temp[i] -= dx;
                }
            }
        };

        template <class Type>
        struct thermofunction
        {
            Type *p;
            typedef scalar (Type::*Tfunc)(const scalar p, const scalar T, autoPtr<typename Type::solution> &sol) const;
            Tfunc pfunc;
            thermofunction(Type *pin, Tfunc fin) : p(pin), pfunc(fin) {}
            void value(const scalarList &in, scalarList &out) const
            {
                scalarList temp(p->X_.size());
                out.resize(1);
                double sum = 0;
                int n_sp = p->X_.size();
                for (int i = 0; i < n_sp - 1; i++)
                {
                    temp[i] = p->X_[i];
                    p->X_[i] = in[i];
                    sum += in[i];
                }
                temp[n_sp - 1] = p->X_[n_sp - 1];
                p->X_[n_sp - 1] = 1 - sum;
                autoPtr<typename Type::solution> sol(p->TPN(in[in.size() - 2], in[in.size() - 1]));
                out[0] = (p->*pfunc)(in[in.size() - 2], in[in.size() - 1], sol);
                // Info << out << endl;
                for (int i = 0; i < n_sp; i++)
                {
                    p->X_[i] = temp[i];
                }
            }

            void derive(const scalarList &in, scalarRectangularMatrix &out)
            {
                scalarList temp = in;
                scalarList tempv;
                value(in, tempv);
                scalarList tempv2;
                scalarList dd;
                double sum = 0;
                for (int i = 0; i < in.size() - 2; i++)
                    sum += in[i];
                scalar dx = 1e-6;
                for (int i = 0; i < in.size(); i++)
                {
                    temp[i] += dx;
                    value(temp, tempv2);
                    dd = (tempv2 - tempv) / dx;
                    for (int j = 0; j < tempv2.size(); j++)
                    {
                        out[i][j] = dd[j];
                    }
                    temp[i] -= dx;
                }
            }
        };
        template <class Type>
        struct thermofunction_psi
        {
            Type *p;
            typedef scalar (Type::*Tfunc)(const scalar p, const scalar T, autoPtr<typename Type::solution> &sol) const;
            thermofunction_psi(Type *pin) : p(pin) {}
            void value(const scalarList &in, scalarList &out, scalar T0) const
            {
                scalar H = in[in.size() - 2];
                scalar P = in[in.size() - 1];
                scalar T = p->T(H, P, T0, &Type::Hs, &Type::Cp, &Type::limit);
                // scalarList temp(p->X_.size());
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

            void derive(const scalarList &in, scalarRectangularMatrix &out, scalar T0)
            {
                // static double maxx = 0;
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
                for (int i = 0; i < grad->size(); i++)
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
                // Info << grad() << endl;
                // Info << out << endl;
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
                // Info << maxx << endl;
                // FatalErrorInFunction << "vaporf 1" << exit(FatalError);
            }
        };
        template <class Type>
        struct thermofunction_Tpsi
        {
            Type *p;
            thermofunction_Tpsi(Type *pin) : p(pin) {}
            void value(const scalarList &in, scalarList &out, scalar T0) const
            {
                scalar H = in[in.size() - 2];
                scalar P = in[in.size() - 1];
                scalar T = p->T(H, P, T0, &Type::Hs, &Type::Cp, &Type::limit);
                // scalarList temp(p->X_.size());
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

            void derive(const scalarList &in, scalarRectangularMatrix &out, scalar T0)
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

        template <class Type>
        struct thermofunction_Tpsivf
        {
            Type *p;
            thermofunction_Tpsivf(Type *pin) : p(pin) {}
            void value(const scalarList &in, scalarList &out, scalar T0) const
            {
                scalar H = in[in.size() - 2];
                scalar P = in[in.size() - 1];
                p->setX(in);
                scalar T = p->T(H, P, T0, &Type::Hs, &Type::Cp, &Type::limit);
                // scalarList temp(p->X_.size());
                out.resize(4);
                /*
                for (int i = 0;i < p->X_.size();i++)
                {
                    temp[i] = in[i];
                }
                */

                autoPtr<typename Type::solution> sol(p->TPN(P, T));
                out[0] = T;
                out[1] = (p->Type::thermoVLE::psi)(P, T, sol);
                out[2] = sol().vaporfra;
                out[3] = p->c(P, T, sol());
            }

            void derive(const scalarList &in, scalarRectangularMatrix &out, scalar T0)
            {
                scalar H = in[in.size() - 2];
                scalar P = in[in.size() - 1];
                scalar T = p->T(H, P, T0, &Type::Hs, &Type::Cp, &Type::limit);

                p->setX(in);
                autoPtr<typename Type::solution> sol(p->TPN(P, T)); // sol(p->TPn_flash_New_TPD(P, T));//sol(p->TPn_flash(P, T));
                autoPtr<scalarRectangularMatrix> grad(p->dTrhoPvfcdXHP_HsP(P, T, sol()));
                out = grad();
            }
            bool valid_in(const scalarList &in)
            {
                bool flag = true;
                for (int i = 0; i < in.size() - 2; i++)
                {
                    flag = flag && (in[i] > 0);
                }
                flag = flag && (in[in.size() - 1] > 0);
                return flag;
            }
        };

        template <class Type>
        struct thermofunction_TPvf
        {
            Type *p;
            thermofunction_TPvf(Type *pin) : p(pin) {}
            void value(const scalarList &in, scalarList &out, scalar T0, scalar P0) const
            {
                scalar E = in[in.size() - 2];
                scalar rho = in[in.size() - 1];
                scalar T, P;
                p->setX(in);
                std::tie(P, T) = p->TrhoE(E, rho, P0, T0);
                /*
                scalarList test(3);
                test[0]=3.0800157865319561e-05;
                test[1]=0.14733120870077421;
                test[2]=0.85263799114136052;
                p->setX(test);
                P=59984981.409090608;
                T=310.99673481821509;*/
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

            void derive(const scalarList &in, scalarRectangularMatrix &out, scalar T0, scalar P0)
            {
                scalar E = in[in.size() - 2];
                scalar rho = in[in.size() - 1];
                scalar T, P;
                std::tie(P, T) = p->TrhoE(E, rho, P0, T0);

                p->setX(in);
                autoPtr<typename Type::solution> sol(p->TPN(P, T)); // sol(p->TPn_flash_New_TPD(P, T));//sol(p->TPn_flash(P, T));
                autoPtr<scalarRectangularMatrix> grad(p->dTPvfdXErho_constC(P, T, sol()));
                // autoPtr<scalarRectangularMatrix> grad(p->dTPvfdXErho(P, T, sol()));
                out = grad();
            }
            bool valid_in(const scalarList &in)
            {
                bool flag = true;
                for (int i = 0; i < in.size() - 2; i++)
                {
                    flag = flag && (in[i] > 0);
                }
                flag = flag && (in[in.size() - 1] > 0);
                return flag;
            }
        };

        template <class Type>
        struct thermofunction_THvf : basicFunction
        {
            Type *p;
            thermofunction_THvf(Type *pin) : p(pin) {}
            void value(const scalarList &in, scalarList &out, scalar T0, autoPtr<typename Type::solution> &sol) const
            {
                scalar rho = in[in.size() - 2];
                scalar P = in[in.size() - 1];
                p->setX(in);
                //scalar T = p->T(rho, P, T0, &Type::rho, &Type::drhodT, &Type::limit);
                //scalar T = p->T(1 / rho, P, T0, &Type::rrho, &Type::drrhodT, &Type::limit);
                scalar T = p->Tsol(1 / rho, P, T0, &Type::rrho, &Type::drrhodT, &Type::limit, sol);

                //autoPtr<typename Type::solution> sol(p->TPN(P, T));
                //sol.reset(p->TPN(P, T).ptr());

                out.resize(4);

                out[0] = T;
                out[1] = p->Hs(P, T, sol());
                out[2] = sol().vaporfra;
                out[3] = p->c_opt(P, T, rho, sol());
            }

            void derive(const scalarList &in, const scalarList &in2, scalarRectangularMatrix &out, scalar T0, autoPtr<typename Type::solution> &sol)
            {
                scalar rho = in[in.size() - 2];
                scalar P = in[in.size() - 1];
                p->setX(in);
                scalar T = in2[0];
                //scalar T = p->T(rho, P, T0, &Type::rho, &Type::drhodT, &Type::limit);

                //sol.reset(p->TPN(P, T).ptr());
                //autoPtr<typename Type::solution> sol(p->TPN(P, T)); // sol(p->TPn_flash_New_TPD(P, T));//sol(p->TPn_flash(P, T));
                autoPtr<scalarRectangularMatrix> grad(p->dTHvfcdXrhoP(P, T, sol()));
                // autoPtr<scalarRectangularMatrix> grad(p->dTPvfdXErho(P, T, sol()));
                out = grad();
            }
            bool valid_in(const scalarList &in)
            {
                bool flag = true;
                for (int i = 0; i < in.size() - 2; i++)
                {
                    flag = flag && (in[i] > 0);
                }
                flag = flag && (in[in.size() - 1] > 0);
                flag = flag && (in[in.size() - 2] > 0);
                return flag;
            }
            static void read(const dictionary &ISATDict_, scalarList &maxLeafsize_in, scalarList &toleranceOut_in, scalarList &initToleranceIn_in, scalarList &scaleIn_in)
            {
                scalarList maxLeafsize_temp(ISATDict_.lookup("maxLeafsize"));
                if (maxLeafsize_temp.size() != 3 && maxLeafsize_temp.size() != maxLeafsize_in.size())
                {
                    FatalErrorInFunction
                        << "maxLeafsize error"
                        << abort(FatalError);
                }
                if (maxLeafsize_temp.size() == 3)
                {

                    maxLeafsize_in[maxLeafsize_in.size() - 1] = maxLeafsize_temp[maxLeafsize_temp.size() - 1];
                    maxLeafsize_in[maxLeafsize_in.size() - 2] = maxLeafsize_temp[maxLeafsize_temp.size() - 2];
                    for (int i = 0; i < maxLeafsize_in.size() - 2; i++)
                    {
                        maxLeafsize_in[i] = maxLeafsize_temp[0];
                    }
                }
                toleranceOut_in = scalarList(ISATDict_.lookup("toleranceOut"));

                scalarList initToleranceIn_temp(ISATDict_.lookup("initToleranceIn"));
                if (initToleranceIn_temp.size() != 3 && initToleranceIn_temp.size() != initToleranceIn_in.size())
                {
                    FatalErrorInFunction
                        << "initToleranceIn error"
                        << abort(FatalError);
                }
                if (initToleranceIn_temp.size() == 3)
                {

                    initToleranceIn_in[initToleranceIn_in.size() - 1] = initToleranceIn_temp[initToleranceIn_temp.size() - 1];
                    initToleranceIn_in[initToleranceIn_in.size() - 2] = initToleranceIn_temp[initToleranceIn_temp.size() - 2];
                    for (int i = 0; i < initToleranceIn_in.size() - 2; i++)
                    {
                        initToleranceIn_in[i] = initToleranceIn_temp[0];
                    }
                }
                scalarList scaleIn_temp(ISATDict_.lookup("scaleIn"));
                if (scaleIn_temp.size() != 3 && scaleIn_temp.size() != scaleIn_in.size())
                {
                    FatalErrorInFunction
                        << "scaleIn error"
                        << abort(FatalError);
                }
                if (scaleIn_temp.size() == 3)
                {

                    scaleIn_in[scaleIn_in.size() - 1] = scaleIn_temp[scaleIn_temp.size() - 1];
                    scaleIn_in[scaleIn_in.size() - 2] = scaleIn_temp[scaleIn_temp.size() - 2];
                    for (int i = 0; i < scaleIn_in.size() - 2; i++)
                    {
                        scaleIn_in[i] = scaleIn_temp[0];
                    }
                }
            }
        };

        template <class Type>
        struct thermofunction_Tpsivfsol
        {
            Type *p;
            thermofunction_Tpsivfsol(Type *pin) : p(pin) {}
            void value(const scalarList &in, scalarList &out, scalar T0) const
            {
                scalar H = in[in.size() - 2];
                scalar P = in[in.size() - 1];
                p->setX(in);
                scalar T = p->T(H, P, T0, &Type::Hs, &Type::Cp, &Type::limit);

                autoPtr<typename Type::solution> sol(p->TPN(P, T));

                out.resize(2 * in.size());

                out[0] = T;
                out[1] = (p->Type::thermoVLE::psi)(P, T, sol);
                out[2] = sol().vaporfra;
                out[3] = p->c(P, T, sol());
                for (int i = 0; i < in.size() - 2; i++)
                {
                    out[4 + i] = sol().X_gas()[i];
                    out[2 + in.size() + i] = sol().X_liq()[i];
                }
            }

            void derive(const scalarList &in, scalarRectangularMatrix &out, scalar T0)
            {
                scalar H = in[in.size() - 2];
                scalar P = in[in.size() - 1];
                scalar T = p->T(H, P, T0, &Type::Hs, &Type::Cp, &Type::limit); // Todo redundant, removde this

                p->setX(in);
                autoPtr<typename Type::solution> sol(p->TPN(P, T)); // sol(p->TPn_flash_New_TPD(P, T));//sol(p->TPn_flash(P, T));
                autoPtr<scalarRectangularMatrix> grad(p->dTrhoPvfcsoldXHP_HsP(P, T, sol()));
                out = grad();
            }
            bool valid_in(const scalarList &in)
            {
                bool flag = true;
                for (int i = 0; i < in.size() - 2; i++)
                {
                    flag = flag && (in[i] > 0);
                }
                flag = flag && (in[in.size() - 1] > 0);
                return flag;
            }
        };

        template <class Type>
        struct thermofunction_TPX_to_Erhovfc : public basicFunction
        {
            Type *p;
            thermofunction_TPX_to_Erhovfc(Type *pin) : p(pin) {}
            void value(const scalarList &in, scalarList &out, autoPtr<typename Type::solution> &sol) const
            {
                scalar T = in[in.size() - 2];
                scalar P = in[in.size() - 1];
                // scalar T, P;
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
                // std::tie(P, T) = p->TrhoE(E, rho, P0, T0);
                // scalarList temp(p->X_.size());
                out.resize(4);
                // out.resize(3);
                /*
                for (int i = 0;i < p->X_.size();i++)
                {
                    temp[i] = in[i];
                }
                */
                p->setX(in);
                sol.reset(p->TPN(P, T).ptr());
                //autoPtr<typename Type::solution> sol(p->TPN(P, T));
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

            void derive(const scalarList &in, scalarList &in2, scalarRectangularMatrix &out , autoPtr<typename Type::solution> &sol)
            {
                scalar T = in[in.size() - 2];
                scalar P = in[in.size() - 1];
                // scalar T, P;
                // std::tie(P, T) = p->TrhoE(E, rho, P0, T0);

                p->setX(in);
                //autoPtr<typename Type::solution> sol(p->TPN(P, T)); // sol(p->TPn_flash_New_TPD(P, T));//sol(p->TPn_flash(P, T));
                autoPtr<scalarRectangularMatrix> grad(p->dErhovfcdXTP_simp(P, T, sol()));
                // autoPtr<scalarRectangularMatrix> grad(p->dTPvfdXErho(P, T, sol()));
                out = grad();
            }
            bool valid_in(const scalarList &in)
            {
                bool flag = true;
                for (int i = 0; i < in.size(); i++)
                {
                    flag = flag && (in[i] > 0);
                }
                return flag;
            }
            static void read(const dictionary &ISATDict_, scalarList &maxLeafsize_in, scalarList &toleranceOut_in, scalarList &initToleranceIn_in, scalarList &scaleIn_in)
            {
                scalarList maxLeafsize_temp(ISATDict_.lookup("maxLeafsize"));
                if (maxLeafsize_temp.size() != 3 && maxLeafsize_temp.size() != maxLeafsize_in.size())
                {
                    FatalErrorInFunction
                        << "maxLeafsize error"
                        << abort(FatalError);
                }
                if (maxLeafsize_temp.size() == 3)
                {

                    maxLeafsize_in[maxLeafsize_in.size() - 1] = maxLeafsize_temp[maxLeafsize_temp.size() - 1];
                    maxLeafsize_in[maxLeafsize_in.size() - 2] = maxLeafsize_temp[maxLeafsize_temp.size() - 2];
                    for (int i = 0; i < maxLeafsize_in.size() - 2; i++)
                    {
                        maxLeafsize_in[i] = maxLeafsize_temp[0];
                    }
                }
                toleranceOut_in = scalarList(ISATDict_.lookup("toleranceOut"));

                scalarList initToleranceIn_temp(ISATDict_.lookup("initToleranceIn"));
                if (initToleranceIn_temp.size() != 3 && initToleranceIn_temp.size() != initToleranceIn_in.size())
                {
                    FatalErrorInFunction
                        << "initToleranceIn error"
                        << abort(FatalError);
                }
                if (initToleranceIn_temp.size() == 3)
                {

                    initToleranceIn_in[initToleranceIn_in.size() - 1] = initToleranceIn_temp[initToleranceIn_temp.size() - 1];
                    initToleranceIn_in[initToleranceIn_in.size() - 2] = initToleranceIn_temp[initToleranceIn_temp.size() - 2];
                    for (int i = 0; i < initToleranceIn_in.size() - 2; i++)
                    {
                        initToleranceIn_in[i] = initToleranceIn_temp[0];
                    }
                }
                scalarList scaleIn_temp(ISATDict_.lookup("scaleIn"));
                if (scaleIn_temp.size() != 3 && scaleIn_temp.size() != scaleIn_in.size())
                {
                    FatalErrorInFunction
                        << "scaleIn error"
                        << abort(FatalError);
                }
                if (scaleIn_temp.size() == 3)
                {

                    scaleIn_in[scaleIn_in.size() - 1] = scaleIn_temp[scaleIn_temp.size() - 1];
                    scaleIn_in[scaleIn_in.size() - 2] = scaleIn_temp[scaleIn_temp.size() - 2];
                    for (int i = 0; i < scaleIn_in.size() - 2; i++)
                    {
                        scaleIn_in[i] = scaleIn_temp[0];
                    }
                }
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

        template <class Type>
        struct thermofunction2
        {
            Type *p;
            typedef scalar (Type::*Tfunc)(const scalar p, const scalar T) const;
            Tfunc pfunc;
            thermofunction2(Type *pin, Tfunc fin) : p(pin), pfunc(fin) {}
            void value(const scalarList &in, scalarList &out) const
            {
                scalarList temp(p->X_.size());
                out.resize(1);
                double sum = 0;
                int n_sp = p->X_.size();
                for (int i = 0; i < n_sp - 1; i++)
                {
                    temp[i] = p->X_[i];
                    p->X_[i] = in[i];
                    sum += in[i];
                }
                temp[n_sp - 1] = p->X_[n_sp - 1];
                p->X_[n_sp - 1] = 1 - sum;
                out[0] = (p->*pfunc)(in[in.size() - 2], in[in.size() - 1]);
                // Info << out << endl;
                for (int i = 0; i < n_sp; i++)
                {
                    p->X_[i] = temp[i];
                }
            }

            void derive(const scalarList &in, scalarRectangularMatrix &out)
            {
                scalarList temp = in;
                scalarList tempv;
                value(in, tempv);
                scalarList tempv2;
                scalarList dd;
                double sum = 0;
                for (int i = 0; i < in.size() - 2; i++)
                    sum += in[i];
                scalar dx = 1e-6;
                for (int i = 0; i < in.size(); i++)
                {
                    temp[i] += dx;
                    value(temp, tempv2);
                    dd = (tempv2 - tempv) / dx;
                    for (int j = 0; j < tempv2.size(); j++)
                    {
                        out[i][j] = dd[j];
                    }
                    temp[i] -= dx;
                }
            }
            static void read(const dictionary &ISATDict_, scalarList &maxLeafsize_in, scalarList &toleranceOut_in, scalarList &initToleranceIn_in, scalarList &scaleIn_in)
            {
                scalarList maxLeafsize_temp(ISATDict_.lookup("maxLeafsize"));
                if (maxLeafsize_temp.size() == 3)
                {

                    maxLeafsize_in[maxLeafsize_in.size() - 1] = maxLeafsize_temp[maxLeafsize_temp.size() - 1];
                    maxLeafsize_in[maxLeafsize_in.size() - 2] = maxLeafsize_temp[maxLeafsize_temp.size() - 2];
                    for (int i = 0; i < maxLeafsize_in.size() - 2; i++)
                    {
                        maxLeafsize_in[i] = maxLeafsize_temp[0];
                    }
                    FatalErrorInFunction
                        << maxLeafsize_in
                        << abort(FatalError);
                }
                toleranceOut_in = scalarList(ISATDict_.lookup("toleranceOut"));
                scalarList initToleranceIn_temp(ISATDict_.lookup("initToleranceIn"));

                scalarList scaleIn_temp(ISATDict_.lookup("scaleIn"));
            }
        };

        template <class Type>
        struct thermofunction3
        {
            Type *p;
            typedef scalar (Type::*Tfunc)(const scalar p, const scalar T) const;
            Tfunc pfunc;
            thermofunction3(Type *pin, Tfunc fin) : p(pin), pfunc(fin) {}
            void value(const scalarList &in, scalarList &out) const
            {
                scalar T = in[in.size() - 2];
                scalar P = in[in.size() - 1];
                scalarList temp(p->X_.size());
                out.resize(1);
                p->setX(in);
                out[0] = (p->*pfunc)(in[in.size() - 2], in[in.size() - 1]);
            }

            void derive(const scalarList &in, scalarRectangularMatrix &out)
            {
                scalarList temp = in;
                scalarList tempv;
                value(in, tempv);
                scalarList tempv2;
                scalarList dd;
                scalar dx = 1e-5;
                for (int i = 0; i < in.size(); i++)
                {
                    temp[i] += dx;
                    value(temp, tempv2);
                    dd = (tempv2 - tempv) / dx;
                    for (int j = 0; j < tempv2.size(); j++)
                    {
                        out[i][j] = dd[j];
                    }
                    temp[i] -= dx;
                }
            }
            bool valid_in(const scalarList &in)
            {
                bool flag = true;
                for (int i = 0; i < in.size(); i++)
                {
                    flag = flag && (in[i] > 0);
                }
                return flag;
            }
        };

        template <class Type>
        struct psi_vffunction
        {
            Type *p;
            typedef scalar (Type::*Tfunc)(const scalar p, const scalar T, autoPtr<typename Type::solution> &sol) const;
            Tfunc ppsi, pvf;
            psi_vffunction(Type *pin, Tfunc ppsi_in, Tfunc pvf_in) : p(pin), ppsi(ppsi_in), pvf(pvf_in) {}
            void value(const scalarList &in, scalarList &out) const
            {
                scalarList temp(p->X_.size());
                out.resize(2);
                double sum = 0;
                int n_sp = p->X_.size();
                for (int i = 0; i < n_sp - 1; i++)
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
                // Info << out << endl;
                for (int i = 0; i < n_sp; i++)
                {
                    p->X_[i] = temp[i];
                }
            }

            void derive(const scalarList &in, scalarRectangularMatrix &out)
            {
                scalarList temp = in;
                scalarList tempv;
                value(in, tempv);
                scalarList tempv2;
                scalarList dd;
                double sum = 0;
                for (int i = 0; i < in.size() - 2; i++)
                    sum += in[i];
                scalar dx = 1e-6;
                for (int i = 0; i < in.size(); i++)
                {
                    temp[i] += dx;
                    value(temp, tempv2);
                    dd = (tempv2 - tempv) / dx;
                    for (int j = 0; j < tempv2.size(); j++)
                    {
                        out[i][j] = dd[j];
                    }
                    temp[i] -= dx;
                }
            }
        };
    }
}

#endif