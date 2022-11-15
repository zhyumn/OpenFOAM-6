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
        struct thermofunction_TrhoP_to_THvf_G_rhoY : basicFunction
        {
            Type *p;
            thermofunction_TrhoP_to_THvf_G_rhoY(Type *pin) : p(pin) {}
            void value(const scalarList &in, scalarList &out, scalar T0, autoPtr<typename Type::solution> &sol) const
            {
                scalar rho = in[in.size() - 2];
                scalar P = in[in.size() - 1];
                p->setX(in);

                scalar T = p->Tsol(1 / rho, P, T0, &Type::rrho, &Type::drrhodT, &Type::limit, sol);

                out.resize(4 + p->X_.size());

                out[0] = T;
                out[1] = p->Hs(P, T, sol());
                out[2] = sol().vaporfra;

                scalar W_G_ = p->W(sol().X_gas());
                out[3] = P * W_G_ / (p->Type::singlePhaseMixtureThermoType::Z(P, T, sol().X_gas(), 1) * RR * 1.0e-03 * T);
                for (int i = 0; i < sol().X_gas().size(); i++)
                {
                    out[i + 4] = sol().X_gas()[i] * (*p)[i].W() * 1e-3 / W_G_;
                }

                /*
                scalar xx = p->c_opt(P, T, rho, sol());

                if (mag(xx - out[3]) > 1e-8)
                {
                    FatalErrorInFunction
                        << "mag(xx-out[3])>1e-8"
                        << abort(FatalError);
                }*/
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
                autoPtr<scalarRectangularMatrix> grad(p->dTHvf_G_rhoY_dXrhoP(P, T, sol()));
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

                //toleranceOut_in = scalarList(ISATDict_.lookup("toleranceOut"));
                scalarList toleranceOut_temp(ISATDict_.lookup("toleranceOut"));
                if (toleranceOut_temp.size() != 5 && toleranceOut_temp.size() != toleranceOut_in.size())
                {
                    FatalErrorInFunction
                        << "toleranceOut error"
                        << abort(FatalError);
                }
                if (toleranceOut_temp.size() == 5)
                {
                    for (int i = 0; i < 4; i++)
                    {
                        toleranceOut_in[i] = toleranceOut_temp[i];
                    }
                    for (int i = 4; i < toleranceOut_in.size(); i++)
                    {
                        toleranceOut_in[i] = toleranceOut_temp[4];
                    }
                }
                else
                {
                    for (int i = 0; i < toleranceOut_in.size(); i++)
                    {
                        toleranceOut_in[i] = toleranceOut_temp[i];
                    }
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
        struct thermofunction_TPX_to_Erhovfc_G_rhoY : public basicFunction
        {
            Type *p;
            thermofunction_TPX_to_Erhovfc_G_rhoY(Type *pin) : p(pin) {}
            void value(const scalarList &in, scalarList &out, autoPtr<typename Type::solution> &sol) const
            {
                scalar T = in[in.size() - 2];
                scalar P = in[in.size() - 1];

                out.resize(5 + p->X_.size());

                p->setX(in);

                /*
                out[0] = p->Hs_ideal(P, T) - P / p->rho_ideal(P, T);
                out[1] = p->rho_ideal(P, T);
                out[2] = 1;
                out[3] = p->c_ideal(P, T);
                */

                sol.reset(p->TPN(P, T).ptr());
                //autoPtr<typename Type::solution> sol(p->TPN(P, T));
                out[0] = p->Hs(P, T, sol()) - P / p->rho(P, T, sol());
                out[1] = p->rho(P, T, sol());
                out[2] = sol().vaporfra;
                out[3] = p->c_opt_full(P, T,out[1], sol());

                scalar W_G_ = p->W(sol().X_gas());
                out[4] = P * W_G_ / (p->Type::singlePhaseMixtureThermoType::Z(P, T, sol().X_gas(), 1) * RR * 1.0e-03 * T);
                for (int i = 0; i < sol().X_gas().size(); i++)
                {
                    out[i + 5] = sol().X_gas()[i] * (*p)[i].W() * 1e-3 / W_G_;
                }

                /*
                FatalErrorInFunction
                    << "Es1: " << E1
                    << "\nEs2: " << E2
                    << "\ndEsdx1: " << (E2-E1)/dx
                    << "\ndEsdx1_ AA: " << dd
                    << exit(FatalError);
                    */
            }

            void derive(const scalarList &in, scalarList &in2, scalarRectangularMatrix &out, autoPtr<typename Type::solution> &sol)
            {
                scalar T = in[in.size() - 2];
                scalar P = in[in.size() - 1];
                // scalar T, P;
                // std::tie(P, T) = p->TrhoE(E, rho, P0, T0);

                p->setX(in);
                //autoPtr<typename Type::solution> sol(p->TPN(P, T)); // sol(p->TPn_flash_New_TPD(P, T));//sol(p->TPn_flash(P, T));
                autoPtr<scalarRectangularMatrix> grad(p->dErhovfc_G_rhoY_dXTP(P, T, sol()));
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

                //toleranceOut_in = scalarList(ISATDict_.lookup("toleranceOut"));
                scalarList toleranceOut_temp(ISATDict_.lookup("toleranceOut"));
                if (toleranceOut_temp.size() != 6 && toleranceOut_temp.size() != toleranceOut_in.size())
                {
                    FatalErrorInFunction
                        << "toleranceOut error"
                        << abort(FatalError);
                }
                if (toleranceOut_temp.size() == 6)
                {
                    for (int i = 0; i < 5; i++)
                    {
                        toleranceOut_in[i] = toleranceOut_temp[i];
                    }
                    for (int i = 5; i < toleranceOut_in.size(); i++)
                    {
                        toleranceOut_in[i] = toleranceOut_temp[5];
                    }
                }
                else
                {
                    for (int i = 0; i < toleranceOut_in.size(); i++)
                    {
                        toleranceOut_in[i] = toleranceOut_temp[i];
                    }
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
        struct thermofunction_XHsP_to_Trhovfc_G_rhoY : basicFunction
        {
            Type *p;
            thermofunction_XHsP_to_Trhovfc_G_rhoY(Type *pin) : p(pin) {}
            void value(const scalarList &in, scalarList &out, scalar T0, autoPtr<typename Type::solution> &sol) const
            {
                scalar Hs = in[in.size() - 2];
                scalar P = in[in.size() - 1];
                scalar T = p->T(Hs, P, T0, &Type::Hs, &Type::Cp, &Type::limit);
                //std::tie(P,T)=p->TrhoEs(E, rho, P0, T0, &Type::Es, &Type::rho);
                //scalar T = p->T(H, P, T0, &Type::Hs, &Type::Cp, &Type::limit);
                //scalarList temp(p->X_.size());
                out.resize(5 + p->X_.size());

                p->setX(in);
                sol.reset(p->TPN(P, T).ptr());
                //autoPtr<typename Type::solution> sol(p->TPN(P, T).ptr());
                out[0] = T;
                out[1] = p->rho(P, T, sol());
                out[2] = sol().vaporfra;
                out[3] = p->c_opt_full(P, T,out[1], sol());


                scalar W_G_ = p->W(sol().X_gas());
                out[4] = P * W_G_ / (p->Type::singlePhaseMixtureThermoType::Z(P, T, sol().X_gas(), 1) * RR * 1.0e-03 * T);
                for (int i = 0; i < sol().X_gas().size(); i++)
                {
                    out[i + 5] = sol().X_gas()[i] * (*p)[i].W() * 1e-3 / W_G_;
                }
            }

            void derive(const scalarList &in, scalarList &in2, scalarRectangularMatrix &out, scalar T0, autoPtr<typename Type::solution> &sol)
            {
                scalar Hs = in[in.size() - 2];
                scalar P = in[in.size() - 1];

                p->setX(in);
                scalar T = in2[0];
                autoPtr<scalarRectangularMatrix> grad(p->dTrhovfc_G_rhoY_dXHsP(P, T, sol()));
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

                //toleranceOut_in = scalarList(ISATDict_.lookup("toleranceOut"));
                scalarList toleranceOut_temp(ISATDict_.lookup("toleranceOut"));
                if (toleranceOut_temp.size() != 6 && toleranceOut_temp.size() != toleranceOut_in.size())
                {
                    FatalErrorInFunction
                        << "toleranceOut error"
                        << abort(FatalError);
                }
                if (toleranceOut_temp.size() == 6)
                {
                    for (int i = 0; i < 5; i++)
                    {
                        toleranceOut_in[i] = toleranceOut_temp[i];
                    }
                    for (int i = 5; i < toleranceOut_in.size(); i++)
                    {
                        toleranceOut_in[i] = toleranceOut_temp[5];
                    }
                }
                else
                {
                    for (int i = 0; i < toleranceOut_in.size(); i++)
                    {
                        toleranceOut_in[i] = toleranceOut_temp[i];
                    }
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