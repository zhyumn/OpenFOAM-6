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
        struct thermofunction_TrhoP_to_THvfc : basicFunction
        {
            Type *p;
            thermofunction_TrhoP_to_THvfc(Type *pin) : p(pin) {}
            void value(const scalarList &in, scalarList &out, scalar T0, typename Type::VLE_D1_data &sol_NIO) const
            {
                scalar rho = in[in.size() - 2];
                scalar P = in[in.size() - 1];
                p->setX(in);
                //scalar T = p->T(rho, P, T0, &Type::rho, &Type::drhodT, &Type::limit);
                //scalar T = p->T(1 / rho, P, T0, &Type::rrho, &Type::drrhodT, &Type::limit);

                //scalar T = P * p->W() / (RR * 1.0e-03 * rho);
                //scalar T = p->Tsol(1 / rho, P, T0, &Type::rrho, &Type::drrhodT, &Type::limit, sol);
                scalar T = p->PV_nflash_NIO(P, rho, T0, sol_NIO);
                /*if (mag(T - T2) > 1e-8)
                {
                    Info << "here1" << endl;
                    Info << T << endl;
                    Info << T2 << endl;
                    Info << endl;
                    NotImplemented;
                }*/

                //typename Type::VLE_data *psolution;
                //p->TPn_flash_New_TPD_Tudisco_NIO(P, T, p->X(), psolution);

                //autoPtr<typename Type::solution> sol(p->TPN(P, T));
                /*scalarList in_test(3);
                in_test[0] = 0.301188931652;
                in_test[1] = 0.645773250779;
                in_test[2] = 0.0530378175691;
                p->setX(in_test);
                T = 345.032844024;
                P = 53365896.3196;
                */

                //sol.reset(p->TPN(P, T).ptr());

                out.resize(4);

                /*
                out[0] = T;
                out[1] = p->Hs_ideal(P, T);
                out[2] = 1;
                //out[3] = p->c_opt(P, T, rho, sol());
                out[3] = p->c_ideal(P, T);
                */

                out[0] = T;
                out[1] = p->Hs_NIO(P, T, sol_NIO);
                out[2] = sol_NIO.vaporFraction; //sol().vaporfra;
                //out[3] = p->c_opt(P, T, rho, sol());
                //out[3] = p->c_opt_full(P, T, rho, sol());
                out[3] = p->c_full_NIO(P, T, rho, sol_NIO);

                /*
                scalar xx = p->c_opt(P, T, rho, sol());

                if (mag(xx - out[3]) > 1e-8)
                {
                    FatalErrorInFunction
                        << "mag(xx-out[3])>1e-8"
                        << abort(FatalError);
                }*/
            }

            void derive(const scalarList &in, const scalarList &in2, scalarRectangularMatrix &out, scalar T0, typename Type::VLE_D1_data &sol_NIO)
            {
                //scalar rho = in[in.size() - 2];
                scalar P = in[in.size() - 1];
                p->setX(in);
                scalar T = in2[0];
                //scalar T = p->T(rho, P, T0, &Type::rho, &Type::drhodT, &Type::limit);

                //sol.reset(p->TPN(P, T).ptr());
                //autoPtr<typename Type::solution> sol(p->TPN(P, T)); // sol(p->TPn_flash_New_TPD(P, T));//sol(p->TPn_flash(P, T));
                //autoPtr<scalarRectangularMatrix> grad(p->dTHvfcdXrhoP(P, T, sol()));
                autoPtr<scalarRectangularMatrix> grad(p->dTHvfcdXrhoP_NIO(P, T, sol_NIO));

                /*                 if (mag(grad1()[0][0] - grad()[0][0]) > 1e-8)
                {
                    Info << "here2" << endl;
                    Info << sol().vaporfra << endl;
                    Info << sol_NIO.vaporFraction << endl;
                    Info << sol().X_gas() << endl;
                    Info << sol_NIO.Gas.X << endl;
                    Info << sol().X_liq() << endl;
                    Info << sol_NIO.Liq.X << endl;
                    Info << grad1()[0][0] << endl;
                    Info << grad()[0][0] << endl;
                    Info << endl;
                    NotImplemented;
                } */

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
                else
                {
                    for (int i = 0; i < maxLeafsize_in.size(); i++)
                    {
                        maxLeafsize_in[i] = maxLeafsize_temp[i];
                    }
                }

                toleranceOut_in = scalarList(ISATDict_.lookup("toleranceOut"));
                scalarList toleranceOut_temp(ISATDict_.lookup("toleranceOut"));
                if (toleranceOut_temp.size() != 4)
                {
                    FatalErrorInFunction
                        << "toleranceOut error"
                        << abort(FatalError);
                }

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
                else
                {
                    for (int i = 0; i < initToleranceIn_in.size(); i++)
                    {
                        initToleranceIn_in[i] = initToleranceIn_temp[i];
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
                else
                {
                    for (int i = 0; i < scaleIn_in.size(); i++)
                    {
                        scaleIn_in[i] = scaleIn_temp[i];
                    }
                }
            }
        };

        template <class Type>
        struct thermofunction_TPX_to_Erhovfc : public basicFunction
        {
            Type *p;
            thermofunction_TPX_to_Erhovfc(Type *pin) : p(pin) {}
            void value(const scalarList &in, scalarList &out, typename Type::VLE_D1_data &sol_NIO) const
            {
                scalar T = in[in.size() - 2];
                scalar P = in[in.size() - 1];
                /*
                //scalar T, P;
                scalarList test(2);
                test[0] = 0.997;
                test[1] = 0.003;
                T = 595;
                P = 5e6;
                p->setX(test);
                autoPtr<typename Type::solution> sol1(p->TPN(P, T));
                scalar E1 = p->Hs(P, T, sol1()); //- P / p->rho(P, T, sol1());

                test[0] = 0.999;
                test[1] = 0.001;
                p->setX(test);
                autoPtr<typename Type::solution> sol2(p->TPN(P, T));
                scalar E2 = p->Hs(P, T, sol2());

                FatalErrorInFunction
                    << exit(FatalError);
                    */
                /*
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

                /*
                out[0] = p->Hs_ideal(P, T) - P / p->rho_ideal(P, T);
                out[1] = p->rho_ideal(P, T);
                out[2] = 1;
                out[3] = p->c_ideal(P, T);
                */

                //sol.reset(p->TPN(P, T).ptr());
                p->TPn_flash_New_TPD_Tudisco_NIO(P, T, p->X_, sol_NIO);
                //autoPtr<typename Type::solution> sol(p->TPN(P, T));
                out[0] = p->Hs_NIO(P, T, sol_NIO) - P / p->rho_NIO(P, T, sol_NIO);
                out[1] = p->rho_NIO(P, T, sol_NIO);
                out[2] = sol_NIO.vaporFraction; //sol().vaporfra;
                //out[3] = p->c(P, T, sol());
                out[3] = p->c_full_NIO(P, T, out[1], sol_NIO);
                //scalar c_test = p->c_full_NIO(P, T, out[1], sol_NIO);
                /*                 if (mag(out[3] - c_test) > 1e-8)
                {
                    Info << out[3] << endl;
                    Info << c_test << endl;
                    Info << endl;
                    NotImplemented;
                } */

                /*
                FatalErrorInFunction
                    << "Es1: " << E1
                    << "\nEs2: " << E2
                    << "\ndEsdx1: " << (E2-E1)/dx
                    << "\ndEsdx1_ AA: " << dd
                    << exit(FatalError);
                    */
            }

            void derive(const scalarList &in, scalarList &in2, scalarRectangularMatrix &out, typename Type::VLE_D1_data &sol_NIO)
            {
                scalar T = in[in.size() - 2];
                scalar P = in[in.size() - 1];
                // scalar T, P;
                // std::tie(P, T) = p->TrhoE(E, rho, P0, T0);

                p->setX(in);
                //autoPtr<typename Type::solution> sol(p->TPN(P, T)); // sol(p->TPn_flash_New_TPD(P, T));//sol(p->TPn_flash(P, T));
                //autoPtr<scalarRectangularMatrix> grad(p->dErhovfcdXTP(P, T, sol()));
                autoPtr<scalarRectangularMatrix> grad(p->dErhovfcdXTP_NIO(P, T, sol_NIO));

                /*if (mag(grad()[0][0] - grad2()[0][0]) > 1e-8)
                {
                    Info << grad2()[0][0] << endl;
                    Info << grad()[0][0] << endl;
                    Info << endl;
                    NotImplemented;
                }*/
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
                else
                {
                    for (int i = 0; i < maxLeafsize_in.size(); i++)
                    {
                        maxLeafsize_in[i] = maxLeafsize_temp[i];
                    }
                }

                toleranceOut_in = scalarList(ISATDict_.lookup("toleranceOut"));
                scalarList toleranceOut_temp(ISATDict_.lookup("toleranceOut"));
                if (toleranceOut_temp.size() != 4)
                {
                    FatalErrorInFunction
                        << "toleranceOut error"
                        << abort(FatalError);
                }

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
                else
                {
                    for (int i = 0; i < initToleranceIn_in.size(); i++)
                    {
                        initToleranceIn_in[i] = initToleranceIn_temp[i];
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
                else
                {
                    for (int i = 0; i < scaleIn_in.size(); i++)
                    {
                        scaleIn_in[i] = scaleIn_temp[i];
                    }
                }
            }
        };

        template <class Type>
        struct thermofunction_ErhoX_to_TPvfc : public basicFunction
        {
            Type *p;
            thermofunction_ErhoX_to_TPvfc(Type *pin) : p(pin) {}
            void value(const scalarList &in, scalarList &out, scalar T0, scalar P0, typename Type::VLE_D1_data &sol_NIO) const
            {
                scalar E = in[in.size() - 2];
                scalar rho = in[in.size() - 1];
                scalar T, P;
                //scalar T1, P1;
                p->setX(in);
                //std::tie(P, T) = p->TrhoEs(E, rho, P0, T0, &Type::Es_vle, &Type::rho, sol);
                //typename Type::VLE_D1_data solnn_NIO(p->X_.size());
                std::tie(P, T) = p->UVn_flash_NIO(E, rho, P0, T0, sol_NIO);

                //if (mag(P1 - P) + mag(T1 - T) > 1e-8)
                //{
                //    Info << P1 << "," << P << endl;
                //    Info << T1 << "," << T << endl;
                //    Info << endl;
                //    NotImplemented;
                //}

                //scalarList in_test(2);
                //in_test[0] = 0.126936802377;
                //in_test[1] = 0.873063197623;
                //p->setX(in_test);
                //typename Type::VLE_data sol_NIO(p->X_.size());
                //p->TPn_flash_New_TPD_Tudisco_NIO(P, T, p->X_, sol_NIO);
                //typename Type::VLE_data sol_NIO(p->X_.size());
                /*
                for (int i = 0; i < 1000000; i++)
                {
                    p->TPn_flash_New_TPD_Tudisco_NIO(5000022.16182, 324.483915823, p->X_, sol_NIO);
                }
                NotImplemented;
                */
                //sol.reset(p->TPN(P, T).ptr());
                //autoPtr<typename Type::solution> sol(p->TPN(P, T));

                out.resize(4);

                out[0] = T;
                out[1] = P;
                out[2] = sol_NIO.vaporFraction;
                //out[2] = sol().vaporfra;

                //Info << sol().vaporfra << endl;
                //Info << sol_NIO.vaporFraction << endl;
                //Info << endl;

                /*                 if (mag(sol().vaporfra - sol_NIO.vaporFraction) > 1e-8)
                {
                    Info << sol().vaporfra << endl;
                    Info << sol_NIO.vaporFraction << endl;
                    Info << endl;
                    NotImplemented;
                } */
                /*if ((1 - sol().vaporfra) * sol().vaporfra > 0.1)
                {
                    Info << P << endl;
                    Info << T << endl;
                    Info << in << endl;
                    NotImplemented;
                }*/

                //out[3] = p->c(P, T, sol());
                //scalar c_test = p->c(P, T, sol());
                out[3] = p->c_NIO(P, T, p->rho_NIO(P, T, sol_NIO), sol_NIO);
                //scalar c_test = p->c_opt_full_NIO(P, T, p->rho_NIO(P,T,solnn_NIO), solnn_NIO);

                /*                 if (mag(out[3] - c_test) > 1e-8)
                {
                    Info << out[3] << endl;
                    Info << c_test << endl;
                    Info << endl;
                    NotImplemented;
                } */

                /*
                                FatalErrorInFunction
                                    << "T1: " << P
                                    << "\nT2: " << P1
                                    << "\ndEsdx1: " << (P1-P) / dx
                                    << "\ndEsdx1_ AA: " << grad()[1][1]
                                    << exit(FatalError);
                                    */
            }

            void derive(const scalarList &in, scalarList &in2, scalarRectangularMatrix &out, scalar T0, scalar P0, typename Type::VLE_D1_data &sol_NIO)
            {
                //scalar E = in[in.size() - 2];
                //scalar rho = in[in.size() - 1];
                scalar &T = in2[0];
                scalar &P = in2[1];
                //T = in2[0];
                //P = in2[1];
                //std::tie(P, T) = p->TrhoEs(E, rho, P0, T0, &Type::Es_ac, &Type::rho);

                p->setX(in);
                //autoPtr<typename Type::solution> sol(p->TPN(P, T)); // sol(p->TPn_flash_New_TPD(P, T));//sol(p->TPn_flash(P, T));
                //sol.reset(p->TPN(P, T).ptr());
                //autoPtr<scalarRectangularMatrix> grad1(p->dTPvfc_G_rhoY_dXErho(P, T, sol()));
                autoPtr<scalarRectangularMatrix> grad(p->dTPvfcdXErho_NIO(P, T, sol_NIO));
                /*               autoPtr<scalarRectangularMatrix> grad0(new scalarRectangularMatrix(p->X_.size() + 2, 4));
                for (int i = 0; i < grad0().m(); i++)
                    for (int j = 0; j < grad0().n(); j++)
                    {
                        grad0()(i, j) = grad1()(i, j);
                    }

                if (mag(grad1()[0][0] - grad()[0][0]) > 1e-8)
                {
                    Info << grad1()[0][0] << endl;
                    Info << grad()[0][0] << endl;
                    Info << endl;
                    NotImplemented;
                } */

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
                else
                {
                    for (int i = 0; i < maxLeafsize_in.size(); i++)
                    {
                        maxLeafsize_in[i] = maxLeafsize_temp[i];
                    }
                }

                toleranceOut_in = scalarList(ISATDict_.lookup("toleranceOut"));

                if (toleranceOut_in.size() != 4)
                {
                    FatalErrorInFunction
                        << "toleranceOut error"
                        << abort(FatalError);
                }

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
                else
                {
                    for (int i = 0; i < initToleranceIn_in.size(); i++)
                    {
                        initToleranceIn_in[i] = initToleranceIn_temp[i];
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
                else
                {
                    for (int i = 0; i < scaleIn_in.size(); i++)
                    {
                        scaleIn_in[i] = scaleIn_temp[i];
                    }
                }
            }
        };
    }
}

#endif