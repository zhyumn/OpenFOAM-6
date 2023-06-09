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
    Foam::PengRobinsonMixture

Description
    PengRobinsonMixture mixture equation of state.

SourceFiles
    PengRobinsonMixtureI.H
    PengRobinsonMixture.C

\*---------------------------------------------------------------------------*/

#ifndef PengRobinsonMixture_H
#define PengRobinsonMixture_H
#include "autoPtr.H"
#include "speciesTable.H"
#include <tuple>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    // Forward declaration of friend functions and operators

    template <class ThermoMixture>
    class PengRobinsonMixture;

    template <class ThermoMixture>
    inline PengRobinsonMixture<ThermoMixture> operator+(
        const PengRobinsonMixture<ThermoMixture> &,
        const PengRobinsonMixture<ThermoMixture> &);

    template <class ThermoMixture>
    inline PengRobinsonMixture<ThermoMixture> operator*(
        const scalar,
        const PengRobinsonMixture<ThermoMixture> &);

    template <class ThermoMixture>
    inline PengRobinsonMixture<ThermoMixture> operator==(
        const PengRobinsonMixture<ThermoMixture> &,
        const PengRobinsonMixture<ThermoMixture> &);

    template <class ThermoMixture>
    Ostream &operator<<(
        Ostream &,
        const PengRobinsonMixture<ThermoMixture> &);

    /*---------------------------------------------------------------------------*\
                           Class PengRobinsonMixture Declaration
\*---------------------------------------------------------------------------*/

    template <class ThermoMixture>
    class PengRobinsonMixture
        : public ThermoMixture
    {
    public:
        scalarListList bico;
        scalarListList oneMbico;
        //labelList IJlistI;
        //labelList IJlistJ;
        label phaseFlag; //gas 1,liq 0

    public:
        typedef typename ThermoMixture::SingleThermoType SingleThermoType;
        typedef typename SingleThermoType::EOS EOS;
        // Constructors

        //- Construct from components

        //- Construct from dictionary
        PengRobinsonMixture(const dictionary &dict, PtrList<SingleThermoType> &);

        //- Construct as named copy
        inline PengRobinsonMixture(const word &name, const PengRobinsonMixture &);

        inline PengRobinsonMixture(const word &name, PtrList<SingleThermoType> &, const speciesTable &, const dictionary &);

        //- Construct and return a clone
        inline autoPtr<PengRobinsonMixture> clone() const;

        // Selector from dictionary
        inline static autoPtr<PengRobinsonMixture> New(
            const dictionary &dict);

        // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return word(ThermoMixture::typeName());
        }
        struct Basic_data
        {
            scalar sqrtT;
            scalarList B, Tsqrt_AbyP;
            Basic_data(label nspecie) : sqrtT(0), B(nspecie), Tsqrt_AbyP(nspecie){};
        };
        struct Phase_data
        {
            Basic_data &basic;
            scalarList X;
            scalar Z, A, B;
            Phase_data(label nspecie, Basic_data &basic_in) : basic(basic_in), X(nspecie), Z(0){};
        };

        // Fundamental properties
        inline scalar A(scalar p, scalar T, const scalarList &Xinp) const;
        inline scalar A(scalar p, scalar T) const { return A(p, T, this->X_); }

        inline void A_A2(scalar p, scalar T, const scalarList &Xinp, scalar &A, scalarList &A2) const;
        inline scalar dAdT(scalar p, scalar T, const scalarList &Xinp) const;
        inline scalar dAdP(scalar p, scalar T, const scalarList &Xinp) const;
        inline scalar dAdxi(scalar p, scalar T, label di, const scalarList &Xinp) const;
        inline scalar d2AdTdxi(scalar p, scalar T, label di, const scalarList &Xinp) const;
        inline scalar B(scalar p, scalar T, const scalarList &Xinp) const;
        inline scalar B(scalar p, scalar T) const { return B(p, T, this->X_); }
        inline scalar dBdT(scalar p, scalar T, const scalarList &Xinp) const;
        inline scalar dBdP(scalar p, scalar T, const scalarList &Xinp) const;
        inline scalar dBdxi(scalar p, scalar T, label di, const scalarList &Xinp) const;
        inline scalar a(scalar p, scalar T, const scalarList &Xinp) const;
        inline scalar b(scalar p, scalar T, const scalarList &Xinp) const;
        inline scalar dadt(scalar p, scalar T, const scalarList &Xinp) const;
        inline scalar dadp(scalar p, scalar T, const scalarList *Xinp = NULL) const;
        inline scalar d2adtdxi(scalar p, scalar T, label di, const scalarList &Xinp) const;
        inline scalar dadxi(scalar p, scalar T, label di, const scalarList &Xinp) const;
        inline scalar ddT_dadt(scalar p, scalar T, const scalarList *Xinp = NULL) const;
        inline scalar d2adt2(scalar p, scalar T, const scalarList &Xinp) const;
        inline autoPtr<scalarList> A2(scalar p, scalar T, const scalarList &Xinp) const;
        inline autoPtr<scalarList> A2(scalar p, scalar T) const { return A2(p, T, this->X_); };
        inline autoPtr<scalarList> dA2dT(scalar p, scalar T, const scalarList &Xinp) const;
        inline autoPtr<scalarList> dA2dP(scalar p, scalar T, const scalarList &Xinp) const;
        inline autoPtr<scalarList> dA2dxi(scalar p, scalar T, label di, const scalarList &Xinp) const;

        inline static scalar cubicSolve(scalar a2, scalar a1, scalar a0, label flag);
        inline static std::tuple<int, double, double, double> cubicSolve(scalar a2, scalar a1, scalar a0);
        inline static std::tuple<int, double, double, double> cubicSolve2(scalar a2, scalar a1, scalar a0);
        //- Return density [kg/m^3]
        inline scalar rho(scalar p, scalar T, const scalarList &Xinp, const label pf = 2) const;
        inline scalar rho(scalar p, scalar T, scalar Zin, const scalarList &Xinp, const label pf = 2) const;
        inline scalar drhodT(scalar p, scalar T, const scalarList &Xinp, const label pf = 2) const;
        inline scalar drhodP(scalar p, scalar T, const scalarList &Xinp, const label pf = 2) const;
        inline scalar drhodT(scalar p, scalar T, scalar Zin, const scalarList &Xinp, const label pf = 2) const;
        inline scalar drhodP(scalar p, scalar T, scalar Zin, const scalarList &Xinp, const label pf = 2) const;
        //inline scalar V(scalar p, scalar T, const label *pf = NULL, const scalarList *Xinp = NULL) const;
        inline scalar Ha(scalar p, scalar T, const scalarList &Xinp, const label pf = 2) const;
        inline scalar Ha(scalar p, scalar T, const scalar Zin, const scalarList &Xinp, const label pf = 2) const;
        inline scalar Hs(scalar p, scalar T, const scalarList &Xinp, const label pf = 2) const;
        inline scalar Hs(scalar p, scalar T, const scalar Zin, const scalarList &Xinp, const label pf = 2) const;
        inline scalar dHadT(scalar p, scalar T, const scalar Zin, const scalarList &Xinp, const label pf = 2) const;
        inline scalar dHadP(scalar p, scalar T, const scalar Zin, const scalarList &Xinp, const label pf = 2) const;
        inline scalar dHadxi(scalar p, scalar T, label di, const scalar Zin, const scalarList &Xinp, const label pf = 2) const;
        inline scalar dHadT(scalar p, scalar T, const scalarList &Xinp, const label pf = 2) const;
        inline scalar dHadP(scalar p, scalar T, const scalarList &Xinp, const label pf = 2) const;
        inline scalar dHadxi(scalar p, scalar T, label di, const scalarList &Xinp, const label pf = 2) const;
        inline scalar dHadxi_opt(scalar p, scalar T, label di, const scalarList &Xinp, const label pf = 2) const;
        inline scalar dHsdT(scalar p, scalar T, const scalar Zin, const scalarList &Xinp, const label pf = 2) const;
        inline scalar dHsdP(scalar p, scalar T, const scalar Zin, const scalarList &Xinp, const label pf = 2) const;
        inline scalar dHsdxi(scalar p, scalar T, label di, const scalar Zin, const scalarList &Xinp, const label pf = 2) const;
        inline scalar dHsdT(scalar p, scalar T, const scalarList &Xinp, const label pf = 2) const;
        inline scalar dHsdP(scalar p, scalar T, const scalarList &Xinp, const label pf = 2) const;
        inline scalar dHsdxi(scalar p, scalar T, label di, const scalarList &Xinp, const label pf = 2) const;
        inline scalar dSdT(scalar p, scalar T, const scalarList &Xinp, const label pf = 2) const;
        inline scalar dSdP(scalar p, scalar T, const scalarList &Xinp, const label pf = 2) const;
        inline scalar dSdxi(scalar p, scalar T, label di, const scalarList &Xinp, const label pf = 2) const;
        inline scalar Hc(const scalarList &Xinp) const;
        inline scalar Hideal(scalar p, scalar T, const scalarList &Xinp) const;
        inline scalar Hideal_Mole(scalar p, scalar T, const scalarList &Xinp) const;
        inline scalar dHidealdT(scalar p, scalar T, const scalarList &Xinp) const;
        inline scalar dHidealdxi(scalar p, scalar T, label di, const scalarList &Xinp) const;
        inline scalar dHidealdP(scalar p, scalar T, const scalarList *Xinp = NULL) const;
        inline scalar dSidealdT(scalar p, scalar T, const scalarList &Xinp) const;
        inline scalar dSidealdP(scalar p, scalar T, const scalarList &Xinp) const;
        inline scalar dSidealdxi(scalar p, scalar T, label di, const scalarList &Xinp) const;
        inline scalar Hig_phase(scalar T, const scalarList *Xinp = NULL) const;
        inline scalar Hig2_phase(scalar T, const scalarList *Xinp = NULL) const;
        inline scalar V(scalar p, scalar T, const label pf = 2, const scalarList *Xinp = NULL) const;
        inline scalar CpMCv(scalar p, scalar T, const label pf = 2, const scalarList *Xinp = NULL) const;
        inline scalar Cp(scalar p, scalar T, const label pf = 2, const scalarList *Xinp = NULL) const;
        inline scalar S(scalar p, scalar T, const scalarList &Xinp, const label pf = 2) const;
        inline scalar Sideal_Mole(scalar p, scalar T, const scalarList &Xinp) const;
        inline scalar Sideal(scalar p, scalar T, const scalarList &Xinp) const;
        inline scalar Gideal_Mole(scalar p, scalar T, const scalarList &Xinp) const;
        inline scalar Gideal(scalar p, scalar T, const scalarList &Xinp) const;
        inline scalar muideal_Mole(scalar p, scalar T, label di, const scalarList *Xinp = NULL) const;
        inline scalar Aideal(scalar p, scalar T, const scalarList &Xinp) const;
        inline scalar Aideal_Mole(scalar p, scalar T, const scalarList *Xinp = NULL) const;
        inline scalar G(scalar p, scalar T, const label pf = 2, const scalarList *Xinp = NULL) const;
        inline scalar G_Mole(scalar p, scalar T, const scalarList &Xinp) const;
        inline scalar G_departure_Mole(scalar p, scalar T, const scalarList &Xinp) const;
        inline scalar G_TPD(scalar p, scalar T, const scalarList &Xinp) const;
        inline scalar A_TPD(scalar p, scalar T, const scalarList &Xinp) const;
        //- Return enthalpy departure [J/kg]
        inline scalar H(const scalar p, const scalar T) const;

        //- Return Cp departure [J/(kg K)]
        //inline scalar Cp(scalar p, scalar T) const;

        //- Return entropy [J/(kg K)]
        inline scalar S(const scalar p, const scalar T) const;

        //- Return compressibility rho/p [s^2/m^2]
        inline scalar psi(scalar p, scalar T) const;

        //- Return compression factor []
        inline scalar Z(scalar p, scalar T, const scalarList &Xinp, const label pf = 2) const;
        inline scalar Z(scalar p, scalar T, scalar A_, scalar B_, const scalarList &Xinp, const label pf = 2) const;
        inline scalar Z(scalar p, scalar T, const label pf = 2) const { return Z(p, T, this->X_, pf); }
        inline std::tuple<scalar, scalar, scalar> Z_coe(scalar p, scalar T, const scalarList &Xinp) const;
        inline std::tuple<double, double> Z_opt(scalar p, scalar T, scalar A, const scalarList &Xinp) const;
        inline std::tuple<double, double> Z_opt2(scalar p, scalar T, scalar A, const scalarList &Xinp) const;
        inline std::tuple<double, double> Z_opt2(scalar p, scalar T, scalar A, scalar B, const scalarList &Xinp) const;
        inline std::tuple<double, double> Z_opt2(scalar p, scalar T, const scalarList &Xinp) const;
        inline std::tuple<int, double, double, double> Z_full(scalar p, scalar T, const scalarList &Xinp) const;
        inline scalar Z_gibbs(scalar p, scalar T, const scalarList &Xinp) const;
        inline scalar Z_gibbs_mix(scalar p, scalar T, const scalarList &Xinp) const;
        inline scalar Z_gibbs_mix(scalar p, scalar T, scalar A_, scalar B_, const scalarList &Xinp) const;
        inline std::tuple<double, bool> Z_gibbs2(scalar p, scalar T, const scalarList &Xinp) const;
        inline scalar vol_gibbs(scalar p, scalar T, const scalarList *Xinp) const;

        inline scalar dZdT(scalar p, scalar T, const scalarList &Xinp, const label pf = 2) const;
        inline scalar dZdT(scalar p, scalar T, scalar Zin, const scalarList &Xinp, const label pf = 2) const;

        inline scalar dZdP(scalar p, scalar T, const scalarList &Xinp, const label pf = 2) const;
        inline scalar dZdP(scalar p, scalar T, scalar Zin, const scalarList &Xinp, const label pf = 2) const;

        inline scalar dZdxi(scalar p, scalar T, label di, const scalarList &Xinp, const label pf = 2) const;
        inline scalar dZdxi(scalar p, scalar T, label di, scalar Zin, const scalarList &Xinp, const label pf = 2) const;

        //- Return molecular weight [kg/mol]
        inline scalar W(const scalarList &Xinp) const;
        inline scalar W() const { return W(this->X_); };

        inline autoPtr<scalarList> fugacityCoefficient(scalar p, scalar T, const scalarList &Xinp, const label pf = 2) const;
        inline autoPtr<scalarList> fugacityCoefficient(scalar p, scalar T, const label pf = 2) const { return fugacityCoefficient(p, T, this->X_, pf); };
        inline autoPtr<scalarList> Ln_fugacityCoefficient(scalar p, scalar T, const scalarList &Xinp, const label pf = 2) const;
        inline bool Ln_fugacityCoefficient_opt(scalar p, scalar T, const scalarList &Xinp, scalarList &fugcoef, scalarList &B_, scalarList &Tsqrt_AbyP_, scalar &z, scalar &g, const label pf) const;
        //inline void Ln_fugacityCoefficient_opt_NIO(scalar p, scalar T, const scalarList &Xinp, scalarList &fugcoef, Phase_data &phaseData) const;
        inline bool Ln_fugacityCoefficient_opt(scalar p, scalar T, const scalarList &Xinp, scalarList &fugcoef, scalarList &B_, scalarList &Tsqrt_AbyP_, const label pf) const;
        inline autoPtr<scalarList> Ln_fugacityCoefficient(scalar p, scalar T, const label pf = 2) const { return Ln_fugacityCoefficient(p, T, this->X_, pf); }
        inline autoPtr<scalarList> ln_fugacityCoefficient(scalar p, scalar T, scalar z, const scalarList &Xinp) const;

        inline autoPtr<scalarList> ddT_Ln_fugacityCoefficient(scalar p, scalar T, const scalarList &Xinp, const label pf = 2) const;

        inline autoPtr<scalarList> ddP_Ln_fugacityCoefficient(scalar p, scalar T, const scalarList &Xinp, const label pf = 2) const;

        inline autoPtr<scalarList> ddxi_Ln_fugacityCoefficient(scalar p, scalar T, label di, const scalarList &Xinp, const label pf = 2) const;

        //- Return species compressibility coefficient in each phase []
        inline void fugcoef_phase(int phasei, scalar p, scalar T, scalar *comp_inp, scalar *fugcoef) const;

        /*
        //- Return mass fraction of 0 species, temporary
        inline void initialization(scalar p, scalar T) const;

        //- Return mass fraction of 0 species, temporary
        inline void initPTXdiagram(scalar press, scalar temp, scalar *comp, scalar *comp_liq, scalar *comp_gas) const;

        //- Generate Thermo table
        inline void genTabulation(scalar press, scalar temp, scalar *comp) const;

        //- Generate Diffusion table
        inline void genDiffusion(scalar press, scalar temp, scalar *comp) const;

        //- Generate Hx flash
        inline void Hn_flash(scalar press, scalar temp, scalar *comp) const;

        //- Generate Miture Critical Point
        inline void initMix_critical(scalar *comp) const;

        //- Return mass fraction of 0 species, temporary
        inline Foam::scalar Evaluate_alpha(scalar press, scalar temp, scalar vaporfra, scalar *comp_liq, scalar *comp_gas, scalar *comp) const;

        //- Return phase compositions and phase state
        inline void initTPD(scalar press, scalar temp, scalar *comp, scalar *comp_liq, scalar *comp_gas, scalar &tpdtest) const;

        //- Return phase mixture ideal internal energy [J/kg]
        inline Foam::scalar sieig_phase(scalar T, scalar *comp_phase) const;

        //- Return phase mixture ideal entropy [J/(kg K)]
        inline Foam::scalar Sig_phase(scalar T, scalar *comp) const;

        //- Return phase mixture ideal enthalpy [J/kg]
        inline Foam::scalar Hig_phase(scalar T, scalar *comp_phase) const;

        //- Return phase mixture ideal enthalpy [J/kg]
        inline Foam::scalar Hig2_phase(scalar T, scalar *comp_phase) const;

        //- Return phase mixture ideal enthalpy [J/kg]
        inline void cvcpig_phase(scalar T, scalar *comp_phase, scalar &cpideal, scalar &cvideal) const;

        //- Return composition []
        inline void calccomp(scalar *comp) const;

        //- Return molecular weight [kg/mol]
        inline scalar mwmix(scalar *comp_inp) const;

        //- Return phase compostions [kg/mol]
        inline scalar yTox_phase(scalar yGasInp, scalar *comp_gas, scalar *comp_liq) const;
*/
        //- Is the equation of state is incompressible i.e. rho != f(p)
        static const bool incompressible = false;

        //- Is the equation of state is isochoric i.e. rho = const
        static const bool isochoric = false;
        /*
        //- Return density [kg/m^3]
        inline scalar rhomix_phase(int phasei, scalar p, scalar T, scalar *comp_inp) const;

        //- Return specific volume m^3/mol]
        inline scalar volmmix_phase(int phasei, scalar p, scalar T, scalar *comp_inp) const;

        //- Return speices enthalpy departure [J/kg]
        inline scalar H_phase(int phasei, const scalar p, const scalar T, scalar *comp) const;

        //- Return speices entropy [J/(kg K)]
        inline scalar S_phase(int phasei, const scalar p, const scalar T, scalar *comp) const;

        //- Return compressibility rho/p [s^2/m^2]
        inline scalar psi_phase(int phasei, scalar p, scalar T, scalar *comp) const;
*/
        //- Return (Cp - Cv) [J/(kg K]
        inline scalar CpMCv(scalar p, scalar T) const;
        /*
        //- Return ABZmix: Amixz, Bmixz, Amix2, Bz []
        inline void ABZmix(scalar p, scalar T, scalar *comp_inp, scalar &Amixz, scalar &Bmixz, scalar *Amix2, scalar *Bz) const;

        //- Return abmix: amix, bmix, dadt, dadt2
        inline void abmix(scalar T, scalar *comp_inp, scalar &amix, scalar &bmix, scalar &dadt, scalar &dadt2) const;

        //- Return compression factor []
        inline scalar Zmix(int phasei, scalar Amixz, scalar Bmixz) const;

        //- Return compression factor root []
        inline void Root_Z(int phasei, scalar Amixz, scalar Bmixz, scalar *rootz) const;

        //- Return species compressibility coefficient in each phase []
        inline void fugcoef_phase(int phasei, scalar p, scalar T, scalar *comp_inp, scalar *fugcoef) const;

        //- Return species compressibility coefficient in mixture []
        inline void fugcoef_mixture(scalar p, scalar T, scalar *comp_inp, scalar *fugcoef_mix, scalar &gibbsenergy) const;

        //- Return phase gibbs energy [J/kg]
        inline void gibbs_phase(int phasei, scalar p, scalar T, scalar *comp_inp, scalar &gibbsenergy) const;

        //- Return equl_constant [J/kg]
        inline void wilsoneq(scalar p, scalar T, scalar *equl_constant) const;

        //- Return density of each species [kg/m3]
        inline Foam::scalar rho_spd(int phasei, label spid, scalar p, scalar T) const;

        //- Return density of each species [kg/m3]
        inline Foam::scalar vol_spd(int phasei, label spid, scalar p, scalar T) const;

        //- Return specific heat capibility for mixture [J/kgK]
        inline Foam::scalar cvmix_phase(int phasei, scalar p, scalar T, scalar *comp_phase) const;
        inline Foam::scalar cpmix_phase(int phasei, scalar p, scalar T, scalar *comp_phase) const;

        //- Return internal energy for mixture [kg/m3]
        inline Foam::scalar siemix_phase(int phasei, scalar p, scalar T, scalar *comp_phase) const;

        //- Return enthalpy for mixture [J/kg]
        inline Foam::scalar Hmix(int phasei, scalar p, scalar T, scalar *comp_phase) const;

        //- Return enthalpy for mixture [J/kg]
        inline Foam::scalar Hmix2(int phasei, scalar p, scalar T, scalar *comp_phase) const;

        //- Return pure component EOS parameters
        inline void ab_spd(scalar p, scalar T, scalar Tcsp, scalar Pcsp, scalar omegasp, scalar &alphaz, scalar &az, scalar &bz, scalar &dadtz) const;

        //- Return enthalpy for spd mixture [J/kg]
        inline Foam::scalar Hspd_mix(scalar p, scalar T, scalar *comp_inp, int spflag) const;

        //- Return soundspeed for mixture [kg/m3]
        inline void soundspeedmix(int phasei, scalar p, scalar T, scalar *comp, scalar &ct_phase, scalar &cs_phase) const;

        //- Return OAD217_function []
        inline Foam::scalar OAD217_function(scalar XX, scalar YY) const;

        //- Return OAD013_function []
        inline void OAD013_function(scalar *comp_in, scalar *comp_out, int &ierr) const;

        //- Return OAD012_function []
        inline void OAD012_function(scalar p, scalar T, scalar *comp_inp, scalar *equalconstant, scalar *comp_liq, scalar *comp_gas, scalar &vaporf, int &ierr) const;

        //- Return OAD012_function []
        inline void OAD012_singlephase(scalar p, scalar T, scalar *comp_inp, scalar *equalconstant, scalar tpdvalue, scalar &vaporfra, scalar *comp_liq, scalar *comp_gas) const;

        //- acceleration
        inline void GDEM_SS(scalar *vt, scalar *vt1, scalar *v2, scalar *vout) const;

        //- TPn_flash
        inline void TPn_flash(scalar p, scalar T, scalar *comp, scalar *comp_liq, scalar *comp_gas, scalar &vaporfra, scalar *equalconstant) const;

        //- TPn_flash
        inline void TPn_simple(scalar p, scalar T, scalar *comp, scalar *comp_liq, scalar *comp_gas, scalar &vaporfra, scalar *equalconstant) const;

        //- TPn_output
        inline void TPn_output(scalar p, scalar T, scalar *comp, scalar *comp_liq, scalar *comp_gas, scalar *svphase, scalar *rhophase, scalar *zphase, scalar *cpphase, scalar *cvphase, scalar *ctphase, scalar *csphase, scalar *siephase, scalar *gamphase, scalar *dpdtphase, scalar *dvdtphase) const;

        //- Checkphasesplit
        inline void Checkphasesplit(scalar p, scalar T, scalar *comp_check, scalar *comp, scalar &tpd0, scalar &trifile, scalar *comp_normal) const;

        //- TPD
        inline void TPD(scalar p, scalar T, scalar *comp, scalar *comp_liq, scalar *comp_gas, scalar &tpdtest) const;

        //- Thermal conductivity
        inline void kappa_phase(int phasei, scalar p, scalar T, scalar *comp_phase, scalar &thcond_phase, scalar &viscos) const;

        //- Mass diffusivity under low pressure, based on Ping's C.P.
        inline scalar Dij_lowP(scalar p, scalar T, scalar *comp_inp) const;

        //- Mass diffusivity corrections for high pressure, based on Ping's C.P.
        inline scalar Dij_correct(scalar p, scalar T, scalar *comp_inp) const;

        //- Mass diffusivity under high pressure, based on Ping's C.P.
        inline scalar Dij_highP(scalar p, scalar T, scalar *comp_inp) const;

        //- Mass diffusivity under low pressure, based on Bellan's binary crtical point
        inline scalar Dij_lowP_pesudo(scalar p, scalar T, scalar *comp_inp) const;

        //- Mass diffusivity corrections for high pressure, based on Bellan's C.P.
        inline scalar Dij_correct_pesudo(scalar p, scalar T, scalar *comp_inp) const;

        //- Mass diffusivity under high pressure, based on Bellan's binary crtical point
        inline scalar Dij_highP_pesudo(scalar p, scalar T, scalar *comp_inp) const;

        //- Binary critical Temperature
        inline scalar TcBi(scalar *comp_inp) const;

        //- Binary critical Pressure
        inline scalar PcBi(scalar *comp_inp) const;

        //- Thermal conductivity
        inline void Mixcritic_Bellan(scalar &Tcmix, scalar &Pcmix, scalar &Vcmix, scalar &Zcmix) const;

        //- Thermal conductivity
        inline void from_sie_ro_to_PT(scalar sie_inp, scalar ro_inp, scalar *comp_inp, scalar &pres, scalar &temp) const;
*/
        // IO

        //- Write to Ostream
        void write(Ostream &os) const;

        // Ostream Operator

        friend Ostream &operator<<<ThermoMixture>(
            Ostream &,
            const PengRobinsonMixture &);
    };
    // * * * * * * * * * * * * * * * * Global data  * * * * * * * * * * * * * * //

    // num_spe for every element in the periodic table
    //extern PengRobinsonMixture num_spe;
    //extern PengRobinsonMixture pr_spe, Tr_spe;

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PengRobinsonMixtureI.H"
#include "PengRobinsonMixture_NIO.H"

#ifdef NoRepository
#include "PengRobinsonMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
