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

Class
    Foam::VLE

Description
    VLE mixture equation of state.

SourceFiles
    VLEI.H
    VLE.C

\*---------------------------------------------------------------------------*/

#ifndef VLE_H
#define VLE_H

#include <float.h>
#include <tuple>
#include "autoPtr.H"
#include "scalarMatrices.H"
#include "PengRobinsonMixture.H"
using std::make_tuple;
using std::tie;
using std::tuple;
//#include "scalarList.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {

    // Forward declaration of friend functions and operators

    template <class ThermoMixture>
    class VLE;

    template <class ThermoMixture>
    inline VLE<ThermoMixture> operator==(
        const VLE<ThermoMixture>&,
        const VLE<ThermoMixture>&);

    template <class ThermoMixture>
    inline VLE<ThermoMixture> operator*(
        const scalar s,
        const VLE<ThermoMixture>&);

    template <class ThermoMixture>
    Ostream& operator<<(
        Ostream&,
        const VLE<ThermoMixture>&);



    /*---------------------------------------------------------------------------*\
                               Class VLE Declaration
    \*---------------------------------------------------------------------------*/

    template <class ThermoMixture>
    class VLE
        : public ThermoMixture {
    public:
        typedef typename ThermoMixture::SingleThermoType SingleThermoType;
        // Constructors

        //- Construct from components
        inline VLE(
            const ThermoMixture& sp,
            const scalar& Tc,
            const scalar& Vc,
            const scalar& Pc,
            const scalar& omega);

        //- Construct from dictionary
        VLE(const dictionary& dict, PtrList<ThermoMixture>& speciesData);

        //- Construct as named copy
        inline VLE(const word& name, const VLE&);

        inline VLE(const word& name, PtrList<SingleThermoType>&, const speciesTable&, const dictionary&);
        //- Construct and return a clone
        inline autoPtr<VLE> clone() const;

        // Selector from dictionary
        inline static autoPtr<VLE> New(
            const dictionary& dict);

        // Member functions

        //- Return the instantiated type name
        static word typeName() {
            return "VLE" + word(ThermoMixture::typeName());
        }
        struct solution {
            autoPtr<scalarList> X_liq, X_gas, equalconstant;
            scalar vaporfra;
            solution() : X_liq(), X_gas(), equalconstant(), vaporfra(0) {};
            solution(const solution& ref) : X_liq(ref.X_liq), X_gas(ref.X_gas), equalconstant(ref.equalconstant), vaporfra(ref.vaporfra) {};
            void reset(const solution& ref) { X_liq = ref.X_liq; X_gas = ref.X_gas; equalconstant = ref.equalconstant; vaporfra = ref.vaporfra; }
        };
        // Fundamental properties

        //- Is the equation of state is incompressible i.e. rho != f(p)
        static const bool incompressible = false;

        //- Is the equation of state is isochoric i.e. rho = const
        static const bool isochoric = false;

        static bool noVLE;
        bool inputK;
        scalarList Kinit;

        //- Return density [kg/m^3]
        inline scalar rho(scalar p, scalar T) const;
        inline scalar rho(scalar p, scalar T, solution& sol) const;
        inline scalar drhodT(scalar p, scalar T, solution& sol) const;
        inline scalar drhodT(scalar p, scalar T) const;
        inline scalar drhodP(scalar p, scalar T, solution& sol) const;
        inline scalar drhodXi(scalar p, scalar T, label di, solution& sol) const;

        inline scalar rho_noVLE(scalar p, scalar T, solution& sol) const;
        inline scalar drhodT_noVLE(scalar p, scalar T, solution& sol) const;
        inline scalar drhodP_noVLE(scalar p, scalar T, solution& sol) const;

        inline scalar dEdT(scalar p, scalar T, solution& sol) const;
        inline scalar dEdP(scalar p, scalar T, solution& sol) const;
        inline scalar dEdXi(scalar p, scalar T, label di, solution& sol) const;

        inline scalar dTdE_rhoX(scalar p, scalar T, solution& sol) const;
        inline scalar dTdrho_EX(scalar p, scalar T, solution& sol) const;
        inline scalar dPdE_rhoX(scalar p, scalar T, solution& sol) const;
        inline scalar dPdrho_EX(scalar p, scalar T, solution& sol) const;

        inline scalar dTdXi_Erho(scalar p, scalar T, label di, solution& sol) const;
        inline scalar dPdXi_Erho(scalar p, scalar T, label di, solution& sol) const;

        inline scalar dvfdE_rhoX(scalar p, scalar T, solution& sol) const;
        inline scalar dvfdrho_EX(scalar p, scalar T, solution& sol) const;
        inline scalar dvfdXi_Erho(scalar p, scalar T, label di, solution& sol) const;
        //- Return enthalpy departure [J/kg]
        inline scalar H(const scalar p, const scalar T) const;

        //- Return Cp departure [J/(kg K]
        inline scalar Cp_num(scalar p, scalar T) const;
        inline scalar Cp(scalar p, scalar T) const;
        inline scalar Cp(scalar p, scalar T, solution& sol) const;
        inline scalar Cp_Hs(scalar p, scalar T, solution& sol) const;
        inline scalar dHadXi(scalar p, scalar T, label di, solution& sol) const;
        inline scalar dHadP(scalar p, scalar T, solution& sol) const;
        inline scalar dHsdXi(scalar p, scalar T, label di, solution& sol) const;
        inline scalar dHsdP(scalar p, scalar T, solution& sol) const;
        inline scalar dSdP(scalar p, scalar T, solution& sol) const;
        inline scalar dSdT(scalar p, scalar T, solution& sol) const;

        //- Return entropy [J/(kg K)]
        inline scalar S(const scalar p, const scalar T) const;

        //- Return compressibility rho/p [s^2/m^2]
        inline scalar psi(scalar p, scalar T) const;
        inline scalar psi(scalar p, scalar T, autoPtr<solution>& sol) const;

        //- Return compression factor []
        inline scalar Z(scalar p, scalar T) const;
        inline scalar Z(scalar p, scalar T, solution& sol) const;

        inline scalar kappa_dev(scalar p, scalar T) const;
        inline scalar mu_dev(scalar p, scalar T) const;
        inline scalar alphah_dev(scalar p, scalar T) const;
        //inline scalar Dimix_dev(scalar p, scalar T,label di) const;

        inline scalar dZdT(scalar p, scalar T, solution& sol) const;
        inline scalar dZdP(scalar p, scalar T, solution& sol) const;
        inline scalar dZdXi(scalar p, scalar T, label di, solution& sol) const;

        //- Return (Cp - Cv) [J/(kg K]
        inline scalar CpMCv(scalar p, scalar T) const;
        inline scalar CpMCv(scalar p, scalar T, solution& sol) const;
        inline scalar Cv_VLE(scalar p, scalar T, solution& sol) const;
        inline scalar Cv_gas(scalar p, scalar T, solution& sol) const;
        inline scalar A(scalar p, scalar T) const;
        //inline scalar B(scalar p, scalar T) const;
        inline scalar Ha(scalar p, scalar T) const;
        inline scalar Ha(scalar p, scalar T, solution& sol) const;
        inline scalar Hc() const;
        inline scalar Hs(scalar p, scalar T) const;
        inline scalar Hs(scalar p, scalar T, solution& sol) const;
        inline scalar G(scalar p, scalar T) const;
        inline scalar dGdT(scalar p, scalar T) const;
        inline scalar dCpdT(scalar p, scalar T) const;
        inline scalar limit(scalar T) const;
        //solve cubic equatoin
        static inline scalar cubicSolve(scalar a0, scalar a1, scalar a2, label flag);  //flag 0 smallest>0 ,flag 1 largest

        inline scalar vaporfra(scalar p, scalar T) const;
        inline scalar vaporfra(scalar p, scalar T, autoPtr<solution>& sol) const;

        inline scalar alpha(scalar p, scalar T, solution& sol) const;
        inline label solve_vaporf(scalar p, scalar T, solution& sol) const;
        inline label noramlize(scalarList& X_in) const;
        inline autoPtr<scalarList> wilsoneq(scalar p, scalar T) const;
        inline autoPtr<scalarList> dvidT(scalar p, scalar T, solution& sol) const;
        inline scalar dvfdT(scalar p, scalar T, solution& sol) const;
        inline scalar dvfdP(scalar p, scalar T, solution& sol) const;
        inline scalar dvfdXi(scalar p, scalar T, label di, solution& sol) const;
        inline autoPtr<scalarList> dXdT_gas(scalar p, scalar T, solution& sol) const;
        inline autoPtr<scalarList> dXdP_gas(scalar p, scalar T, solution& sol) const;
        inline autoPtr<scalarList> dXdT_liq(scalar p, scalar T, solution& sol) const;
        inline autoPtr<scalarList> dXdP_liq(scalar p, scalar T, solution& sol) const;
        inline scalar dWdT_gas(scalar p, scalar T, solution& sol) const;
        inline scalar dWdP_gas(scalar p, scalar T, solution& sol) const;
        inline scalar dWdXi_gas(scalar p, scalar T, label di, solution& sol) const;
        inline scalar dWdT_liq(scalar p, scalar T, solution& sol) const;
        inline scalar dWdXi_liq(scalar p, scalar T, label di, solution& sol) const;
        inline scalar dWdXi(label di) const;

        inline autoPtr<scalarList> dvidP(scalar p, scalar T, solution& sol) const;

        inline autoPtr<scalarList> dvidXi(scalar p, scalar T, label di, solution& sol) const;
        inline autoPtr<scalarList> dXdXi_gas(scalar p, scalar T, label di, solution& sol) const;
        inline autoPtr<scalarList> dXdXi_liq(scalar p, scalar T, label di, solution& sol) const;
        inline autoPtr<scalarList> myinitial(scalar p, scalar T) const;
        autoPtr<solution> TPn_flash(scalar p, scalar T) const;
        autoPtr<solution> TPn_flash_New(scalar p, scalar T) const;
        autoPtr<solution> TPn_flash_New2(scalar p, scalar T, bool stable, bool isVapor, Foam::scalarList& K_init) const;
        autoPtr<solution> TPn_flash_Matheis(scalar p, scalar T, const scalarList& Xinp, bool stable, bool isVapor, Foam::scalarList& K_init, bool& success) const;
        autoPtr<solution> TPn_flash_New_TPD_Tudisco(scalar p, scalar T, const scalarList& Xinp) const;
        autoPtr<solution> TPn_flash_New_TPD_Tudisco(scalar p, scalar T) const { return TPn_flash_New_TPD_Tudisco(p, T, this->X_); };
        autoPtr<solution> TPn_flash_New_TPD(scalar p, scalar T) const;
        inline tuple<bool, autoPtr<scalarList>, bool> solveTPD_BFGS(scalar p, scalar T) const;
        inline tuple<bool, Foam::autoPtr<List<scalarList>>, bool, autoPtr<scalarList>> solveTPD_BFGS_v2(scalar p, scalar T, const scalarList& Xinp) const;
        inline tuple<bool, Foam::autoPtr<List<scalarList>>, bool, autoPtr<scalarList>> solveTPD_BFGS_v2(scalar p, scalar T) const { return solveTPD_BFGS_v2(p, T, this->X_); };
        void TPn_validation(scalar p, scalar T) const;

        //HPX system
        inline scalar dPdP() const { return 1.0; }
        inline scalar dPdXi(label) const { return 0; }
        inline scalar dPdH() const { return 0; }
        inline scalar dXidP(label) const { return 0; }
        inline scalar dXidXi(label ai, label bi) const { return ai == bi ? 1.0 : 0.0; }
        inline scalar dXidH(label) const { return 0; }
        inline scalar dTdP_HP(scalar p, scalar T, solution& sol) const;
        inline scalar dTdH_HP(scalar p, scalar T, solution& sol) const;
        inline scalar dTdXi_HP(scalar p, scalar T, label di, solution& sol) const;

        inline scalar drhodH_HP(scalar p, scalar T, solution& sol) const;
        inline scalar drhodP_HP(scalar p, scalar T, solution& sol) const;
        inline scalar drhodXi_HP(scalar p, scalar T, label di, solution& sol) const;
        inline scalar drhoPdH_HP(scalar p, scalar T, solution& sol) const;
        inline scalar drhoPdP_HP(scalar p, scalar T, solution& sol) const;
        inline scalar drhoPdXi_HP(scalar p, scalar T, label di, solution& sol) const;

        inline scalar dvfdH_HP(scalar p, scalar T, solution& sol) const;
        inline scalar dvfdP_HP(scalar p, scalar T, solution& sol) const;
        inline scalar dvfdXi_HP(scalar p, scalar T, label di, solution& sol) const;

        inline scalar dTdP_HsP(scalar p, scalar T, solution& sol) const;
        inline scalar dTdH_HsP(scalar p, scalar T, solution& sol) const;
        inline scalar dTdXi_HsP(scalar p, scalar T, label di, solution& sol) const;

        inline scalar drhodH_HsP(scalar p, scalar T, solution& sol) const;
        inline scalar drhodP_HsP(scalar p, scalar T, solution& sol) const;
        inline scalar drhodXi_HsP(scalar p, scalar T, label di, solution& sol) const;
        inline scalar drhoPdH_HsP(scalar p, scalar T, solution& sol) const;
        inline scalar drhoPdP_HsP(scalar p, scalar T, solution& sol) const;
        inline scalar drhoPdXi_HsP(scalar p, scalar T, label di, solution& sol) const;
        inline autoPtr<scalarList> drhoPdXHP_HsP(scalar p, scalar T, solution& sol) const;
        inline autoPtr<scalarRectangularMatrix> dTrhoPdXHP_HsP(scalar p, scalar T, solution& sol) const;
        inline autoPtr<scalarRectangularMatrix> dTrhoPvfdXHP_HsP(scalar p, scalar T, solution& sol) const;
        inline autoPtr<scalarRectangularMatrix> dTrhoPvfcdXHP_HsP(scalar p, scalar T, solution& sol) const;
        inline autoPtr<scalarRectangularMatrix> dTrhoPvfcsoldXHP_HsP(scalar p, scalar T, solution& sol) const;
        inline autoPtr<scalarRectangularMatrix> dTPvfdXErho(scalar p, scalar T, solution& sol) const;
        inline autoPtr<scalarRectangularMatrix> dTPvfdXErho_constC(scalar p, scalar T, solution& sol) const;
        //inline autoPtr<scalarRectangularMatrix> dTPvfsoldXErho_constC(scalar p, scalar T, solution& sol) const;
        inline autoPtr<scalarRectangularMatrix> dTPvfdXErho_simp(scalar p, scalar T, solution& sol) const;
        inline autoPtr<scalarRectangularMatrix> dTHvfcdXrhoP(scalar p, scalar T, solution& sol)const;

        inline autoPtr<scalarRectangularMatrix> dErhovfcdXTP_simp(scalar p, scalar T, solution& sol) const;

        inline tuple<scalar, scalar, scalar, scalar> dTPdErho(scalar p, scalar T, solution& sol) const;


        inline scalar dvfdH_HsP(scalar p, scalar T, solution& sol) const;
        inline scalar dvfdP_HsP(scalar p, scalar T, solution& sol) const;
        inline scalar dvfdXi_HsP(scalar p, scalar T, label di, solution& sol) const;

        inline scalar kappaT(scalar p, scalar T, solution& sol) const;
        inline scalar alphaP(scalar p, scalar T, solution& sol) const;
        inline scalar kappaS(scalar p, scalar T, solution& sol) const;
        inline scalar c(scalar p, scalar T, solution& sol) const;
        inline scalar c(scalar p, scalar T) const;

        inline scalar kappaT_noVLE(scalar p, scalar T, solution& sol) const;
        inline scalar alphaP_noVLE(scalar p, scalar T, solution& sol) const;
        inline scalar kappaS_noVLE(scalar p, scalar T, solution& sol) const;
        inline scalar c_noVLE(scalar p, scalar T, solution& sol) const;
        inline scalar c_noVLE(scalar p, scalar T) const;


        //mutable solution temp_sol;
        //mutable scalar temp_p,temp_T;
        /* friend VLE operator==<ThermoMixture>(
               const VLE &,
               const VLE &);

           friend VLE operator*<ThermoMixture>(
               const scalar &,
               const VLE &);
               */
               // IO

               //- Write to Ostream
        void write(Ostream& os) const;

        // Ostream Operator

        friend Ostream& operator<<<ThermoMixture>(
            Ostream&,
            const VLE&);
    };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    template<class T>
    void autoPtrReset(autoPtr<T>& ap, const autoPtr<T>& bp)
    {

        ap.reset(const_cast<autoPtr<T>&>(bp).ptr());
        const_cast<autoPtr<T>&>(bp).reset();
        //pa.reset(pb.operator->());
        //autoPtr<T>(bp).reset();
    }
}  // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //




#include "VLEI.H"
#include "VLE_AF.H"

#ifdef NoRepository
#include "VLE.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //