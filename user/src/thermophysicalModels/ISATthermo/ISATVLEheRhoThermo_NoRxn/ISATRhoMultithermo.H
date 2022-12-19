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
#include "ISATFunction.H"
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

        template <class ThermoMixture, template <class> class Type>
        class ISATRhoMultithermo;

        template <class ThermoMixture, template <class> class Type>
        inline ISATRhoMultithermo<ThermoMixture, Type> operator+(
            const ISATRhoMultithermo<ThermoMixture, Type> &,
            const ISATRhoMultithermo<ThermoMixture, Type> &);

        template <class ThermoMixture, template <class> class Type>
        inline ISATRhoMultithermo<ThermoMixture, Type> operator*(
            const scalar,
            const ISATRhoMultithermo<ThermoMixture, Type> &);

        template <class ThermoMixture, template <class> class Type>
        inline ISATRhoMultithermo<ThermoMixture, Type> operator==(
            const ISATRhoMultithermo<ThermoMixture, Type> &,
            const ISATRhoMultithermo<ThermoMixture, Type> &);

        template <class ThermoMixture, template <class> class Type>
        Ostream &operator<<(
            Ostream &,
            const ISATRhoMultithermo<ThermoMixture, Type> &);

        /*---------------------------------------------------------------------------*\
                                   Class ISATRhoMultithermo Declaration
        \*---------------------------------------------------------------------------*/

        template <class ThermoMixture, template <class> class Type>
        class ISATRhoMultithermo
            : public ThermoMixture,
              public Type<ISATRhoMultithermo<ThermoMixture, Type>>
        {
            // Private data

            //- Convergence tolerance of energy -> temperature inversion functions
            static const scalar tol_;

            //- Max number of iterations in energy->temperature inversion functions
            static const int maxIter_;
            dictionary ISATDict_;

            ISATRhoMultithermo<ThermoMixture, Type> *temp_p;
            thermofunction_TPX_to_Erhovfc<ISATRhoMultithermo<ThermoMixture, Type>> Erhovfcfun_XTP;
            thermofunction_TrhoP_to_THvfc<ISATRhoMultithermo<ThermoMixture, Type>> THvfcfun_XrhoP;
            thermofunction_ErhoX_to_TPvfc<ISATRhoMultithermo<ThermoMixture, Type>> TPvfcfun_XErho;

            mutable ISATmanager<thermofunction_TPX_to_Erhovfc<ISATRhoMultithermo<ThermoMixture, Type>>> Erhovfctree_XTP;
            mutable ISATmanager<thermofunction_TrhoP_to_THvfc<ISATRhoMultithermo<ThermoMixture, Type>>> THvfctree_XrhoP;
            mutable ISATmanager<thermofunction_ErhoX_to_TPvfc<ISATRhoMultithermo<ThermoMixture, Type>>> TPvfctree_XErho;

            // Private Member Functions

            //- Return the temperature corresponding to the value of the
            //  thermodynamic property f, given the function f = F(p, T)
            //  and dF(p, T)/dT
            inline scalar T(
                scalar f,
                scalar p,
                scalar T0,
                scalar (ISATRhoMultithermo::*F)(const scalar, const scalar) const,
                scalar (ISATRhoMultithermo::*dFdT)(const scalar, const scalar) const,
                scalar (ISATRhoMultithermo::*limit)(const scalar) const) const;

            inline scalar Tsol(
                scalar f,
                scalar p,
                scalar T0,
                scalar (ISATRhoMultithermo<ThermoMixture, Type>::*F)(const scalar, const scalar, typename ISATRhoMultithermo<ThermoMixture, Type>::solution &) const,
                scalar (ISATRhoMultithermo<ThermoMixture, Type>::*dFdT)(const scalar, const scalar)
                    const,
                scalar (ISATRhoMultithermo<ThermoMixture, Type>::*limit)(const scalar) const,
                autoPtr<typename ISATRhoMultithermo<ThermoMixture, Type>::solution> &sol) const;
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
            inline std::tuple<scalar, scalar> TrhoEs(
                scalar e,
                scalar rho,
                scalar p0,
                scalar T0,
                scalar (ISATRhoMultithermo<ThermoMixture, Type>::*Fe)(const scalar, const scalar, typename ISATRhoMultithermo<ThermoMixture, Type>::solution &) const,
                scalar (ISATRhoMultithermo<ThermoMixture, Type>::*Frho)(const scalar, const scalar, typename ISATRhoMultithermo<ThermoMixture, Type>::solution &) const,
                autoPtr<typename ISATRhoMultithermo<ThermoMixture, Type>::solution> &retSol);

        public:
            //- The thermodynamics of the individual species'
            typedef ThermoMixture thermoVLE;
            typedef ISATRhoMultithermo<ThermoMixture, Type> thermoType;
            typedef typename ThermoMixture::SingleThermoType SingleThermoType;
            typedef typename ThermoMixture::singlePhaseMixtureThermoType singlePhaseMixtureThermoType;

            // Constructors

            //- Construct from components
            inline ISATRhoMultithermo(const ThermoMixture &sp);
            //- Construct from dictionary
            ISATRhoMultithermo(const dictionary &dict, PtrList<ThermoMixture> &speciesData);

            //- Construct as named copy
            inline ISATRhoMultithermo(const word &name, const ISATRhoMultithermo &);

            inline ISATRhoMultithermo(const word &name, PtrList<SingleThermoType> &, const speciesTable &, const dictionary &, ISATRhoMultithermo<ThermoMixture, Type> *pin);

            // Member Functions

            //- Return the instantiated type name
            static word typeName()
            {
                return ThermoMixture::typeName();
            }

            //- Name of Enthalpy/Internal energy
            static inline word heName();

            inline void newTimeStep();

            inline bool newLoop();
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
            //inline std::tuple<scalar, scalar, scalar, scalar > Tpsivfc_XHP(const scalar p, const scalar T, const scalar T0) const;
            inline std::tuple<scalar, scalar, scalar, scalar> THvfc_XrhoP(const scalar p, const scalar T, const scalar T0) const;
            inline std::tuple<scalar, scalar, scalar, scalar, scalar> THvfc_G_rhoY_XrhoP(const scalar p, const scalar T, const scalar T0, scalarList &Y) const;
            inline std::tuple<scalar, scalar, scalar, scalar, scalar> TPvfc_G_rhoY_XErho(const scalar E, const scalar rho, const scalar T0, const scalar P0, scalarList &Y) const;
            inline std::tuple<scalar, scalar, scalar, scalar> TPvfc_Erho(const scalar E, const scalar rho, const scalar T0, const scalar P0) const;
            inline std::tuple<scalar, scalar, scalar, scalar> TPvfc_XErho(const scalar E, const scalar rho, const scalar T0, const scalar P0) const;
            //inline std::tuple<scalar, scalar, scalar, scalar, autoPtr<typename ThermoMixture::solution>> TPvfcsol_Erho(const scalar E, const scalar rho, const scalar T0, const scalar P0) const;
            inline std::tuple<scalar, scalar, scalar, scalar, autoPtr<typename ThermoMixture::solution>> Tpsivfcsol_XHP(const scalar H, const scalar P, const scalar T0) const;

            inline scalar mu(const scalar p, const scalar T) const;
            inline scalar alphah(const scalar p, const scalar T) const;
            inline scalar kappa(const scalar p, const scalar T) const;
            inline scalar vaporfra(const scalar p, const scalar T) const;
            inline void psi_vf(const scalar p, const scalar T, scalar &psi_out, scalar &vf_out) const;

            inline void mute_show();
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
            inline scalar THP(const scalar H, const scalar P, const scalar T0) const;

            inline std::tuple<scalar, scalar, scalar> rhovfc_ISAT(const scalar p, const scalar T) const;
            inline std::tuple<scalar, scalar, scalar, scalar> Erhovfc_ISAT(const scalar p, const scalar T) const;
            inline std::tuple<scalar, scalar, scalar, scalar> Erhovfc_XPT(const scalar p, const scalar T) const;
            inline std::tuple<scalar, scalar, scalar, scalar, scalar> Erhovfc_G_rhoY_ISAT(const scalar p, const scalar T, scalarList &Y) const;

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
            inline scalar Kx(
                const scalar p,
                const scalar T) const;

            //- Equilibrium constant [] i.t.o. number of moles
            //  For low pressures (where the gas mixture is near perfect)
            //  Kn = Kp(n*pstd/p)^nui where n = number of moles in mixture
            inline scalar Kn(
                const scalar p,
                const scalar T,
                const scalar n) const;

            // Energy->temperature  inversion functions

            //- Temperature from enthalpy or internal energy
            //  given an initial temperature T0
            inline scalar THE(
                const scalar H,
                const scalar p,
                const scalar T0) const;

            //- Temperature from sensible enthalpy given an initial T0
            inline scalar THs(
                const scalar Hs,
                const scalar p,
                const scalar T0) const;

            //- Temperature from absolute enthalpy
            //  given an initial temperature T0
            inline scalar THa(
                const scalar H,
                const scalar p,
                const scalar T0) const;

            //- Temperature from sensible internal energy
            //  given an initial temperature T0
            inline scalar TEs(
                const scalar E,
                const scalar p,
                const scalar T0) const;

            //- Temperature from absolute internal energy
            //  given an initial temperature T0
            inline scalar TEa(
                const scalar E,
                const scalar p,
                const scalar T0) const;

            inline std::tuple<scalar, scalar> TrhoE(
                const scalar e,
                const scalar rho,
                const scalar p0,
                const scalar T0);

            // Derivative term used for Jacobian

            //- Derivative of B (acooding to Niemeyer et al.) w.r.t. temperature
            inline scalar dKcdTbyKc(const scalar p, const scalar T) const;

            //- Derivative of cp w.r.t. temperature
            inline scalar dcpdT(const scalar p, const scalar T) const;

            // I-O

            //- Write to Ostream
            void write(Ostream &os) const;

            // Member operators

            inline void operator+=(const ISATRhoMultithermo &);
            inline void operator*=(const scalar);

            // Friend operators

            friend ISATRhoMultithermo operator+<ThermoMixture, Type>(
                const ISATRhoMultithermo &,
                const ISATRhoMultithermo &);

            friend ISATRhoMultithermo operator*<ThermoMixture, Type>(
                const scalar s,
                const ISATRhoMultithermo &);

            friend ISATRhoMultithermo operator==<ThermoMixture, Type>(
                const ISATRhoMultithermo &,
                const ISATRhoMultithermo &);

            // Ostream Operator

            friend Ostream &operator<<<ThermoMixture, Type>(
                Ostream &,
                const ISATRhoMultithermo &);

            template <class TypeName>
            friend struct Hnflashfunction;

            template <class TypeName>
            friend struct thermofunction_psi;
            template <class TypeName>
            friend struct thermofunction_Tpsi;
            template <class TypeName>
            friend struct thermofunction_Tpsivf;
            template <class TypeName>
            friend struct thermofunction_Tpsivfsol;
            template <class TypeName>
            friend struct thermofunction_TrhoP_to_THvfc;
            template <class TypeName>
            friend struct thermofunction_ErhoX_to_TPvfc;
        };

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    } // End namespace species
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ISATRhoMultithermoI.H"

#ifdef NoRepository
#include "ISATRhoMultithermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
