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
#include "parData.H"

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
                const scalar &rho = in[in.size() - 2];
                const scalar &P = in[in.size() - 1];
                p->setX(in);
                scalar T = p->PV_nflash_NIO(P, rho, T0, sol_NIO);

                out[0] = T;
                out[1] = p->Hs_NIO(P, T, sol_NIO);
                out[2] = sol_NIO.vaporFraction;
                out[3] = p->c_full_NIO(P, T, p->rho_NIO(P, T, sol_NIO), sol_NIO);
            }

            void derive(const scalarList &in, const scalarList &in2, scalarRectangularMatrix &out, scalar T0, typename Type::VLE_D1_data &sol_NIO)
            {
                const scalar &P = in[in.size() - 1];
                p->setX(in);
                const scalar &T = in2[0];
                autoPtr<scalarRectangularMatrix> grad(p->dTHvfcdXrhoP_NIO(P, T, sol_NIO));
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
        struct thermofunction_XEP_to_Trhovfc : basicFunction
        {
            Type *p;
            thermofunction_XEP_to_Trhovfc(Type *pin) : p(pin) {}
            void value(const scalarList &in, scalarList &out, scalar T0, typename Type::VLE_D1_data &sol_NIO) const
            {
                const scalar &E = in[in.size() - 2];
                const scalar &P = in[in.size() - 1];
                p->setX(in);
                scalar T = p->EP_nflash_NIO(E, P, T0, sol_NIO);

                out[0] = T;
                out[1] = p->rho_NIO(P, T, sol_NIO);
                out[2] = sol_NIO.vaporFraction;
                out[3] = p->c_full_NIO(P, T, out[1], sol_NIO);
            }

            void derive(const scalarList &in, const scalarList &in2, Foam::SMatrix<double> &out, scalar T0, typename Type::VLE_D1_data &sol_NIO)
            {
                const scalar &P = in[in.size() - 1];
                p->setX(in);
                const scalar &T = in2[0];
                scalarRectangularMatrix out_tmp(in.size(), 4);
                autoPtr<scalarRectangularMatrix> grad(p->dTrhovfcdXEP_NIO(P, T, sol_NIO));
                out_tmp = grad();
                for (int i = 0; i < in.size(); i++)
                    for (int j = 0; j < 4; j++)
                    {
                        out(i, j) = out_tmp[i][j];
                    }
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
        struct thermofunction_TrhoP_to_TEvfc : basicFunction
        {
            Type *p;
            thermofunction_TrhoP_to_TEvfc(Type *pin) : p(pin) {}
            void value(const scalarList &in, scalarList &out, scalar T0, typename Type::VLE_D1_data &sol_NIO) const
            {
                const scalar &rho = in[in.size() - 2];
                const scalar &P = in[in.size() - 1];
                p->setX(in);
                scalar T = p->PV_nflash_NIO(P, rho, T0, sol_NIO);

                out[0] = T;
                out[1] = p->Es_NIO(P, T, sol_NIO);
                out[2] = sol_NIO.vaporFraction;
                out[3] = p->c_full_NIO(P, T, p->rho_NIO(P, T, sol_NIO), sol_NIO);
            }

            void derive(const scalarList &in, const scalarList &in2, Foam::SMatrix<double> &out, scalar T0, typename Type::VLE_D1_data &sol_NIO)
            {
                const scalar &P = in[in.size() - 1];
                p->setX(in);
                const scalar &T = in2[0];
                scalarRectangularMatrix out_tmp(in.size(), 4);
                p->dTEvfcdXrhoP_NIO(P, T, sol_NIO, out_tmp);
                for (int i = 0; i < in.size(); i++)
                    for (int j = 0; j < 4; j++)
                    {
                        out(i, j) = out_tmp[i][j];
                    }
                //autoPtr<scalarRectangularMatrix> grad(p->dTEvfcdXrhoP_NIO(P, T, sol_NIO,out));
                //out = grad();
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
                const scalar &T = in[in.size() - 2];
                const scalar &P = in[in.size() - 1];
                p->setX(in);
                p->TPn_flash_New_TPD_Tudisco_NIO(P, T, p->X_, sol_NIO);
                out[0] = p->Hs_NIO(P, T, sol_NIO) - P / p->rho_NIO(P, T, sol_NIO);
                out[1] = p->rho_NIO(P, T, sol_NIO);
                out[2] = sol_NIO.vaporFraction;
                out[3] = p->c_full_NIO(P, T, out[1], sol_NIO);
            }

            void derive(const scalarList &in, scalarList &in2, Foam::SMatrix<double> &out, typename Type::VLE_D1_data &sol_NIO)
            {
                const scalar &T = in[in.size() - 2];
                const scalar &P = in[in.size() - 1];
                p->setX(in);
                scalarRectangularMatrix out_tmp(in.size(), 4);
                p->dErhovfcdXTP_NIO(P, T, sol_NIO, out_tmp);

                for (int i = 0; i < in.size(); i++)
                    for (int j = 0; j < 4; j++)
                    {
                        out(i, j) = out_tmp[i][j];
                    }
                //autoPtr<scalarRectangularMatrix> grad(p->dErhovfcdXTP_NIO(P, T, sol_NIO));
                //out = grad();
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
                const scalar &E = in[in.size() - 2];
                const scalar &rho = in[in.size() - 1];
                scalar T, P;
                p->setX(in);
                std::tie(P, T) = p->UVn_flash_NIO(E, rho, P0, T0, sol_NIO);

                out[0] = T;
                out[1] = P;
                out[2] = sol_NIO.vaporFraction;
                out[3] = p->c_NIO(P, T, p->rho_NIO(P, T, sol_NIO), sol_NIO);
            }

            void derive(const scalarList &in, scalarList &in2, Foam::SMatrix<double> &out, scalar T0, scalar P0, typename Type::VLE_D1_data &sol_NIO)
            {
                const scalar &T = in2[0];
                const scalar &P = in2[1];
                p->setX(in);

                scalarRectangularMatrix out_tmp(in.size(), 4);
                p->dTPvfcdXErho_NIO(P, T, sol_NIO, out_tmp);
                for (int i = 0; i < in.size(); i++)
                    for (int j = 0; j < 4; j++)
                    {
                        out(i, j) = out_tmp[i][j];
                    }
                //autoPtr<scalarRectangularMatrix> grad(p->dTPvfcdXErho_NIO(P, T, sol_NIO));
                //out = grad();
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