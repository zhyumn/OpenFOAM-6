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
    Foam::species::multithermo

Description
    Basic thermodynamics type based on the use of fitting functions for
    cp, h, s obtained from the template argument type multithermo.  All other
    properties are derived from these primitive functions.

SourceFiles
    multithermoI.H
    multithermo.C

\*---------------------------------------------------------------------------*/

#ifndef multithermo_H
#define multithermo_H

#include "IOstreams.H"
#include "thermodynamicConstants.H"
#include "speciesTable.H"
#include <tuple>
using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace species
{

// Forward declaration of friend functions and operators

template<class ThermoMixture, template<class> class Type> class multithermo;

template<class ThermoMixture, template<class> class Type>
inline multithermo<ThermoMixture, Type> operator+
(
    const multithermo<ThermoMixture, Type>&,
    const multithermo<ThermoMixture, Type>&
);

template<class ThermoMixture, template<class> class Type>
inline multithermo<ThermoMixture, Type> operator*
(
    const scalar,
    const multithermo<ThermoMixture, Type>&
);

template<class ThermoMixture, template<class> class Type>
inline multithermo<ThermoMixture, Type> operator==
(
    const multithermo<ThermoMixture, Type>&,
    const multithermo<ThermoMixture, Type>&
);


template<class ThermoMixture, template<class> class Type>
Ostream& operator<<
(
    Ostream&,
    const multithermo<ThermoMixture, Type>&
);


/*---------------------------------------------------------------------------*\
                           Class multithermo Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoMixture, template<class> class Type>
class multithermo
:
    public ThermoMixture,
    public Type<multithermo<ThermoMixture, Type>>
{
    // Private data

        //- Convergence tolerance of energy -> temperature inversion functions
        static const scalar tol_;

        //- Max number of iterations in energy->temperature inversion functions
        static const int maxIter_;


    // Private Member Functions

        //- Return the temperature corresponding to the value of the
        //  thermodynamic property f, given the function f = F(p, T)
        //  and dF(p, T)/dT
        inline scalar T
        (
            scalar f,
            scalar p,
            scalar T0,
            scalar (multithermo::*F)(const scalar, const scalar) const,
            scalar (multithermo::*dFdT)(const scalar, const scalar) const,
            scalar (multithermo::*limit)(const scalar) const
        ) const;

        inline std::tuple<scalar,scalar> TrhoEs
        (
            scalar e,
            scalar rho,
            scalar p0,
            scalar T0,
            scalar (multithermo<ThermoMixture, Type>::*Fe)(const scalar, const scalar) const,
            scalar (multithermo<ThermoMixture, Type>::*Frho)(const scalar, const scalar) const
        ) const;

/*
        inline scalar T_IQI
        (
            scalar f,
            scalar p,
            scalar T0,
            scalar (multithermo::*F)(const scalar, const scalar) const,
            scalar (multithermo::*dFdT)(const scalar, const scalar) const,
            scalar (multithermo::*limit)(const scalar) const
        ) const;
*/

public:

    //- The thermodynamics of the individual species'
    typedef multithermo<ThermoMixture, Type> thermoType;
    typedef typename ThermoMixture::SingleThermoType SingleThermoType;

    // Constructors

        //- Construct from components
        inline multithermo(const ThermoMixture& sp);

        //- Construct from dictionary
        multithermo(const dictionary& dict,PtrList<ThermoMixture> &speciesData);

        //- Construct as named copy
        inline multithermo(const word& name, const multithermo&);

        
        inline multithermo(const word& name,PtrList<SingleThermoType>&,const speciesTable& , const dictionary&);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return ThermoMixture::typeName();
        }

        //- Name of Enthalpy/Internal energy
        static inline word heName();


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

            inline std::tuple<scalar,scalar> TrhoE
            (
                const scalar e,
                const scalar rho,
                const scalar p0,
                const scalar T0
            ) const;
        // Derivative term used for Jacobian

            //- Derivative of B (acooding to Niemeyer et al.) w.r.t. temperature
            inline scalar dKcdTbyKc(const scalar p, const scalar T) const;

            //- Derivative of cp w.r.t. temperature
            inline scalar dcpdT(const scalar p, const scalar T) const;


        // I-O

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        inline void operator+=(const multithermo&);
        inline void operator*=(const scalar);


    // Friend operators

        friend multithermo operator+ <ThermoMixture, Type>
        (
            const multithermo&,
            const multithermo&
        );

        friend multithermo operator* <ThermoMixture, Type>
        (
            const scalar s,
            const multithermo&
        );

        friend multithermo operator== <ThermoMixture, Type>
        (
            const multithermo&,
            const multithermo&
        );
        


    // Ostream Operator

        friend Ostream& operator<< <ThermoMixture, Type>
        (
            Ostream&,
            const multithermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace species
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "multithermoI.H"

#ifdef NoRepository
    #include "multithermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
