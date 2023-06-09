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
    Foam::chungTransportMixture

Description
    Constant properties Transport package.
    Templated into a given thermodynamics package (needed for thermal
    conductivity).

SourceFiles
    chungTransportMixtureI.H
    chungTransportMixture.C

\*---------------------------------------------------------------------------*/

#ifndef chungTransportMixture_H
#define chungTransportMixture_H
#include <tuple>
using std::make_tuple;
using std::tie;
using std::tuple;
#include "IOstreams.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    // Forward declaration of friend functions and operators

    template <class ThermoMixture>
    class chungTransportMixture;

    template <class ThermoMixture>
    inline chungTransportMixture<ThermoMixture> operator+(
        const chungTransportMixture<ThermoMixture> &,
        const chungTransportMixture<ThermoMixture> &);

    template <class ThermoMixture>
    inline chungTransportMixture<ThermoMixture> operator*(
        const scalar,
        const chungTransportMixture<ThermoMixture> &);

    template <class ThermoMixture>
    Ostream &operator<<(
        Ostream &,
        const chungTransportMixture<ThermoMixture> &);

    /*---------------------------------------------------------------------------*\
                           Class chungTransportMixture Declaration
\*---------------------------------------------------------------------------*/

    template <class ThermoMixture>
    class chungTransportMixture
        : public ThermoMixture
    {
        // Private data
        scalarListList w_global;
        scalarListList sigma_global;
        scalarListList e_k0_global;
        scalarListList kij_global;
        scalarListList mwij_global;
        scalarListList pow3_sigma_global;
        scalarListList sqr_sigma_global;
        scalarListList sqrt_mwij_global;
        scalarListList mu_M_global;
        scalarListList Mmd_global;
        scalarListList sigmd_global;
        scalarListList sqrt_Mmd_global;
        scalarListList sqr_sigmd_global;
        scalarListList mw_binary_M_global;
        scalarListList e_k_mix_M_global;
        scalarListList w_mix_M_global;

        mutable scalarListList Dimix_opt_list1;
        mutable scalarListList Dimix_opt_list2;
        mutable scalarList Dimix_opt_list3;
        mutable scalarListList kappa_mu_opt_list1;

        // Private Member Functions

        //- Construct from components
        inline chungTransportMixture(
            const ThermoMixture &t);

    public:
        typedef typename ThermoMixture::SingleThermoType SingleThermoType;
        // Constructors

        //- Construct as named copy
        inline chungTransportMixture(const word &, const chungTransportMixture &);

        //- Construct from dictionary
        chungTransportMixture(const dictionary &dict, PtrList<SingleThermoType> &speciesData);

        inline chungTransportMixture(const word &name, PtrList<SingleThermoType> &, const speciesTable &, const dictionary &);

        //- Construct and return a clone
        inline autoPtr<chungTransportMixture> clone() const;

        // Selector from dictionary
        inline static autoPtr<chungTransportMixture> New(const dictionary &dict);

        inline void chung_init();

        // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return ThermoMixture::typeName();
        }

        //- Dynamic viscosity [kg/ms]
        inline scalar mu(const scalar p, const scalar T, const label pf = 2, const scalarList *Xinp = NULL) const;

        inline scalar mu_rho(const scalar p, const scalar T, scalar rho, const scalarList &Xinp) const;
        inline scalar mu_rho(const scalar p, const scalar T, scalar rho) const { return mu_rho(p, T, rho, this->X_); };
        inline scalar mu_opt(const scalar p, const scalar T, scalar rho, const scalarList &Xinp) const;
        inline scalar mu_opt(const scalar p, const scalar T, scalar rho) const { return mu_opt(p, T, rho, this->X_); };
        /*inline scalar mu(const scalar p, const scalar T, const label pf = 2) const
        {
            return mu(p, T, pf, NULL);
        }*/

        inline scalar phi(scalar PRR, scalar TRR) const;
        inline scalar phi_opt(scalar PRR, scalar TRR) const;

        inline scalar Dimix(const scalar p, const scalar T, label speciei) const;
        inline scalar Dimix_opt(const scalar p, const scalar T, label speciei) const;
        
        //- Thermal conductivity [W/mK]
        inline scalar kappa(const scalar p, const scalar T, const label pf = 2, const scalarList *Xinp = NULL) const;

        inline scalar kappa_rho(const scalar p, const scalar T, scalar rho, const scalarList &Xinp) const;
        inline scalar kappa_rho(const scalar p, const scalar T, scalar rho) const { return kappa_rho(p, T, rho, this->X_); };
        inline scalar kappa_opt(const scalar p, const scalar T, scalar rho, const scalarList &Xinp) const;
        inline scalar kappa_opt(const scalar p, const scalar T, scalar rho) const { return kappa_opt(p, T, rho, this->X_); };
        inline tuple<scalar, scalar> kappa_mu_opt(const scalar p, const scalar T, scalar rho) const;
        //- Thermal diffusivity of enthalpy [kg/ms]
        inline scalar alphah(const scalar p, const scalar T) const;

        // Species diffusivity
        // inline scalar D(const scalar p, const scalar T) const;

        //- Write to Ostream
        void write(Ostream &os) const;

        // Member operators

        // inline void operator=(const chungTransportMixture&);

        inline void operator+=(const chungTransportMixture &);

        inline void operator*=(const scalar);

        // Friend operators

        friend chungTransportMixture operator+<ThermoMixture>(
            const chungTransportMixture &,
            const chungTransportMixture &);

        friend chungTransportMixture operator*<ThermoMixture>(
            const scalar,
            const chungTransportMixture &);

        // Ostream Operator

        friend Ostream &operator<<<ThermoMixture>(
            Ostream &,
            const chungTransportMixture &);
    };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "chungTransportMixtureI.H"

#ifdef NoRepository
#include "chungTransportMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
