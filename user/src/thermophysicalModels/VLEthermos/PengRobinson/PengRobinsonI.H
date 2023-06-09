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

\*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* PengRobinson EOS                                                          */
/* Thermal solver                                                            */
/* Developed by   Ping Yi,  03/01/2019,   UMN                                */
/* mix means each phase mixture; ave means two-phase mixture                 */
/*---------------------------------------------------------------------------*/

#include "PengRobinson.H"
#include "specie.H"
#include "mathematicalConstants.H"
#include "OFstream.H"
#include "IFstream.H"
#include "string.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template <class Specie>
inline Foam::PengRobinson<Specie>::PengRobinson(
    const Specie &sp,
    const scalar &Tc,
    const scalar &Vc,
    const scalar &Pc,
    const scalar &omega)
    : Specie(sp),
      Tc_(Tc),
      Vc_(Vc),
      Pc_(Pc),
      omega_(omega),
      phaseFlag(1)
{
    Zc_ = Pc_ * Vc_ / (RR * Tc_);
    sqrt_rPc_ = sqrt(1 / Pc_);

    if (omega_ > 0.49)
    {
        kappa_ = omega_ * (omega_ * (omega_ * (0.016666) - 0.164423) + 1.48503) + 0.379642;
    }
    else
    {
        kappa_ = omega_ * (omega_ * (-0.26992) + 1.54226) + 0.37464;
    }

    sqrt_rTc_ = sqrt(1 / Tc_);
    coef_ = sqrt(0.45724) * sqrt_rPc_ * Tc_;
    Aa_ = coef_ * (1.0 + kappa_);
    Ab_ = -coef_ * kappa_ * sqrt_rTc_;

    dAdTa_ = -0.45724 * Tc_ * sqr(kappa_) / Pc_;
    dAdTb_ = 3 * 0.45724 * kappa_ * (1 + kappa_) / Pc_ * Tc_ * sqrt(Tc_);
    dAdTc_ = -2 * 0.45724 * sqr(Tc_) * sqr(1 + kappa_) / Pc_;

    dAdPa_ = 0.45724 * Tc_ * sqr(kappa_) / Pc_;
    dAdPb_ = -2 * 0.45724 * kappa_ * (1 + kappa_) / Pc_ * Tc_ * sqrt(Tc_);
    dAdPc_ = 0.45724 * sqr(Tc_) * sqr(1 + kappa_) / Pc_;

    aa_ = RR * 1.0e-03 * Aa_;
    ab_ = RR * 1.0e-03 * Ab_;

    dadta_ = 0.45724 * sqr(RR * 1.0e-03) * Tc_ * sqr(kappa_) / Pc_;
    dadtb_ = -0.45724 * sqr(RR * 1.0e-03) * kappa_ * (1 + kappa_) / Pc_ * Tc_ * sqrt(Tc_);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class Specie>
inline Foam::PengRobinson<Specie>::PengRobinson(
    const word &name,
    const PengRobinson &pg)
    : Specie(name, pg),
      Tc_(pg.Tc_),
      Vc_(pg.Vc_),
      Zc_(pg.Zc_),
      Pc_(pg.Pc_),
      omega_(pg.omega_)
{
    sqrt_rPc_ = pg.sqrt_rPc_;
    sqrt_rTc_ = pg.sqrt_rTc_;
    kappa_ = pg.kappa_;
    coef_ = pg.coef_;

    Aa_ = pg.Aa_;
    Ab_ = pg.Ab_;

    dAdTa_ = pg.dAdTa_;
    dAdTb_ = pg.dAdTb_;
    dAdTc_ = pg.dAdTc_;

    dAdPa_ = pg.dAdPa_;
    dAdPb_ = pg.dAdPb_;
    dAdPc_ = pg.dAdPc_;

    aa_ = pg.aa_;
    ab_ = pg.ab_;

    dadta_ = pg.dadta_;
    dadtb_ = pg.dadtb_;
}

template <class Specie>
inline Foam::autoPtr<Foam::PengRobinson<Specie>>
Foam::PengRobinson<Specie>::clone() const
{
    return autoPtr<PengRobinson<Specie>>(
        new PengRobinson<Specie>(*this));
}

template <class Specie>
inline Foam::autoPtr<Foam::PengRobinson<Specie>>
Foam::PengRobinson<Specie>::New(
    const dictionary &dict)
{
    return autoPtr<PengRobinson<Specie>>(
        new PengRobinson<Specie>(dict));
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template <class Specie>
inline Foam::scalar Foam::PengRobinson<Specie>::rho(
    scalar p,
    scalar T) const
{
    const scalar Z = this->Z(p, T);
    return p / (Z * this->R() * T);
}
template <class Specie>
inline Foam::scalar Foam::PengRobinson<Specie>::A(
    scalar P,
    scalar T) const
{
    scalar kappa;
    scalar alpha;
    if (omega_ > 0.49)
    {
        kappa = omega_ * (omega_ * (omega_ * (0.016666) - 0.164423) + 1.48503) + 0.379642;
    }
    else
    {
        kappa = omega_ * (omega_ * (-0.26992) + 1.54226) + 0.37464;
    }
    alpha = sqr(1.0 + kappa * (1.0 - sqrt(T / Tc_)));
    return 0.45724 * alpha * P / Pc_ / sqr(T / Tc_);
}
template <class Specie>
inline Foam::scalar Foam::PengRobinson<Specie>::Tsqrt_AbyP(
    scalar sqrtT) const
{
    //scalar kappa;
    //scalar sqrtalpha;
    /*if (omega_ > 0.49)
    {
        kappa = omega_ * (omega_ * (omega_ * (0.016666) - 0.164423) + 1.48503) + 0.379642;
    }
    else
    {
        kappa = omega_ * (omega_ * (-0.26992) + 1.54226) + 0.37464;
    }
    */
    //sqrtalpha = 1.0 + kappa_ * (1.0 - sqrtT * sqrt_rTc_);

    //return coef_ * sqrtalpha;
    return Aa_ + Ab_ * sqrtT;
}
template <class Specie>
inline Foam::scalar Foam::PengRobinson<Specie>::sqrt_ATbyP_M(
    scalar T) const
{
    //scalar kappa;
    scalar sqrtalpha;
    /*
    if (omega_ > 0.49)
    {
        kappa = omega_ * (omega_ * (omega_ * (0.016666) - 0.164423) + 1.48503) + 0.379642;
    }
    else
    {
        kappa = omega_ * (omega_ * (-0.26992) + 1.54226) + 0.37464;
    }
    */
    sqrtalpha = 1.0 + kappa_ * (1.0 - sqrt(T) / sqrt(Tc_));
    //scalar ce=sqrt(0.45724)*sqrt_rPc_*Tc_;
    return sqrt_rPc_ * Tc_ * sqrtalpha;
}

template <class Specie>
inline Foam::scalar Foam::PengRobinson<Specie>::dAdT(
    scalar P,
    scalar T) const
{
    scalar kappa;
    scalar alpha;
    if (omega_ > 0.49)
    {
        kappa = omega_ * (omega_ * (omega_ * (0.016666) - 0.164423) + 1.48503) + 0.379642;
    }
    else
    {
        kappa = omega_ * (omega_ * (-0.26992) + 1.54226) + 0.37464;
    }
    alpha = sqr(1.0 + kappa * (1.0 - sqrt(T / Tc_)));
    scalar dalphadT = -kappa * (1.0 + kappa * (1.0 - sqrt(T / Tc_))) / sqrt(T * Tc_);
    return 0.45724 * dalphadT * P / Pc_ / sqr(T / Tc_) - 2 * 0.45724 * alpha * P / Pc_ / sqr(T / Tc_) / T;
}

template <class Specie>
inline Foam::scalar Foam::PengRobinson<Specie>::sqrTdAdTbyP(
    scalar rSqrtT,
    scalar rT) const
{
    return dAdTa_ + dAdTb_ * rSqrtT + dAdTc_ * rT;
}

template <class Specie>
inline Foam::scalar Foam::PengRobinson<Specie>::dAdP(
    scalar P,
    scalar T) const
{
    scalar kappa;
    scalar alpha;
    if (omega_ > 0.49)
    {
        kappa = omega_ * (omega_ * (omega_ * (0.016666) - 0.164423) + 1.48503) + 0.379642;
    }
    else
    {
        kappa = omega_ * (omega_ * (-0.26992) + 1.54226) + 0.37464;
    }
    alpha = sqr(1.0 + kappa * (1.0 - sqrt(T / Tc_)));
    return 0.45724 * alpha / Pc_ / sqr(T / Tc_);
}

template <class Specie>
inline Foam::scalar Foam::PengRobinson<Specie>::TdAdP(
    scalar rSqrtT,
    scalar rT) const
{

    return dAdPa_ + dAdPb_ * rSqrtT + dAdPc_ * rT;
}

template <class Specie>
inline Foam::scalar Foam::PengRobinson<Specie>::B(
    scalar P,
    scalar T) const
{
    return 0.07780 * P / Pc_ * Tc_ / T;
}

template <class Specie>
inline Foam::scalar Foam::PengRobinson<Specie>::dBdT(
    scalar P,
    scalar T) const
{
    return -0.07780 * P / Pc_ * Tc_ / (T * T);
}

template <class Specie>
inline Foam::scalar Foam::PengRobinson<Specie>::dBdP(
    scalar P,
    scalar T) const
{
    return 0.07780 / Pc_ * Tc_ / T;
}

template <class Specie>
inline Foam::scalar Foam::PengRobinson<Specie>::EOSalpha(
    scalar P,
    scalar T) const
{
    scalar kappa;
    if (omega_ > 0.49)
    {
        kappa = omega_ * (omega_ * (omega_ * (0.016666) - 0.164423) + 1.48503) + 0.379642;
    }
    else
    {
        kappa = omega_ * (omega_ * (-0.26992) + 1.54226) + 0.37464;
    }
    return sqr(1.0 + kappa * (1.0 - sqrt(T / Tc_)));
}
template <class Specie>
inline Foam::scalar Foam::PengRobinson<Specie>::a(
    scalar P,
    scalar T) const
{
    return 0.45724 * EOSalpha(P, T) * sqr(RR * 1.0e-03 * Tc_) / Pc_;
}

template <class Specie>
inline Foam::scalar Foam::PengRobinson<Specie>::sqrta(
    scalar sqrtT) const
{
    return aa_ + ab_ * sqrtT;
}
template <class Specie>
inline Foam::scalar Foam::PengRobinson<Specie>::dadt(
    scalar P,
    scalar T) const
{
    scalar kappa;
    if (omega_ > 0.49)
    {
        kappa = omega_ * (omega_ * (omega_ * (0.016666) - 0.164423) + 1.48503) + 0.379642;
    }
    else
    {
        kappa = omega_ * (omega_ * (-0.26992) + 1.54226) + 0.37464;
    }
    return 0.45724 * (1.0 + kappa * (1.0 - sqrt(T / Tc_))) * sqr(RR * 1.0e-03 * Tc_) / Pc_ * (-kappa / sqrt(Tc_ * T));
}

template <class Specie>
inline Foam::scalar Foam::PengRobinson<Specie>::dadt_opt(
    scalar rSqrtT) const
{
    return dadta_ + dadtb_ * rSqrtT;
}
template <class Specie>
inline Foam::scalar Foam::PengRobinson<Specie>::dadp(
    scalar P,
    scalar T) const
{
    return 0.0;
}
template <class Specie>
inline Foam::scalar Foam::PengRobinson<Specie>::d2adt2(
    scalar P,
    scalar T) const
{
    scalar kappa;
    if (omega_ > 0.49)
    {
        kappa = omega_ * (omega_ * (omega_ * (0.016666) - 0.164423) + 1.48503) + 0.379642;
    }
    else
    {
        kappa = omega_ * (omega_ * (-0.26992) + 1.54226) + 0.37464;
    }
    return 0.45724 * ((1.0 + kappa * (1.0 - sqrt(T / Tc_))) * (kappa / (2 * T * sqrt(Tc_ * T))) + kappa * kappa / (2 * Tc_ * T)) * sqr(RR * 1.0e-03 * Tc_) / Pc_;
}
template <class Specie>
inline Foam::scalar Foam::PengRobinson<Specie>::b(
    scalar P,
    scalar T) const
{
    return 0.07780 * RR * 1.0e-03 * Tc_ / Pc_;
}
template <class Specie>
inline Foam::scalar Foam::PengRobinson<Specie>::dbdT(
    scalar P,
    scalar T) const
{
    return 0;
}
template <class Specie>
inline Foam::scalar Foam::PengRobinson<Specie>::dbdP(
    scalar P,
    scalar T) const
{
    return 0;
}
/*
inline Foam::scalar Foam::PengRobinson<Specie>::c(
    scalar P,
    scalar T) const
{
    scalar kappa;
    scalar alpha;
    if (omega_ > 0.49)
    {
        kappa = omega_ * (omega_ * (omega_ * (0.016666) - 0.164423) + 1.48503) + 0.379642;
    }
    else
    {
        kappa = omega_ * (omega_ * (-0.26992) + 1.54226) + 0.37464;
    }
    return kappa / sqrt(Tc_);
}
*/

template <class Specie>
inline Foam::scalar Foam::PengRobinson<Specie>::cubicSolve(
    scalar a2,
    scalar a1,
    scalar a0,
    label flag)
{

    scalar Q = (3.0 * a1 - a2 * a2) / 9.0;
    scalar Rl = (9.0 * a2 * a1 - 27.0 * a0 - 2.0 * a2 * a2 * a2) / 54.0;
    scalar Q3 = Q * Q * Q;
    scalar D = Q3 + Rl * Rl;
    scalar root = -1;
    if (D <= 0)
    {
        scalar th = ::acos(Rl / sqrt(-Q3));
        scalar qm = 2 * sqrt(-Q);
        scalar r1 = qm * cos(th / 3) - a2 / 3;
        scalar r2 =
            qm * cos((th + 2 * constant::mathematical::pi) / 3) - a2 / 3;
        scalar r3 =
            qm * cos((th + 4 * constant::mathematical::pi) / 3) - a2 / 3;

        // three root;
        if (flag == 1) //gas
        {
            root = max(r1, max(r2, r3));
        }
        else //liquid
        {
            root = max(r1, max(r2, r3));
            root = min(r1 > 1e-3 ? r1 : root, root);
            root = min(r2 > 1e-3 ? r2 : root, root);
            root = min(r3 > 1e-3 ? r3 : root, root);
        }
    }
    else
    {
        // One root is real: in supercritical zone
        scalar D05 = sqrt(D);
        scalar S = sign(Rl + D05) * pow(mag(Rl + D05), 1.0 / 3.0);
        scalar Tl = sign(Rl - D05) * pow(mag(Rl - D05), 1.0 / 3.0);
        /*
        scalar RlD = Rl + D05;
        scalar S = 0.0;
        if (RlD > 0.0)
        {
            S = sign(Rl + D05) pow(mag(Rl + D05), 1.0 / 3.0);
        }
        else
        {
            S = -pow(-RlD, 1.0 / 3.0);
        }
        scalar Tl = 0;
        if (D05 > Rl)
        {
            Tl = -pow(mag(Rl - D05), 1.0 / 3.0);
        }
        else
        {
            Tl = pow(Rl - D05, 1.0 / 3.0);
        }
        */
        root = S + Tl - a2 / 3;
    }

    if (root <= 0.0)
        FatalErrorInFunction << "Z is negative! Z = " << root << exit(FatalError);
    return root;
}

//can be used for phases or mixtures, depending on the Azin and Bzin, only one root
template <class Specie>
inline Foam::scalar Foam::PengRobinson<Specie>::Z(
    scalar p,
    scalar T) const
{
    scalar A = this->A(p, T), B = this->B(p, T);
    scalar a2 = B - 1.0;
    scalar a1 = A - 2.0 * B - 3.0 * B * B;
    scalar a0 = -A * B + B * B + B * B * B;
    return cubicSolve(a2, a1, a0, phaseFlag);
}

template <class Specie>
inline Foam::scalar Foam::PengRobinson<Specie>::H(scalar p, scalar T) const // ideal gas  //TODO impreve code structure
{
    return 0;
}

template <class Specie>
inline Foam::scalar Foam::PengRobinson<Specie>::Cp(scalar p, scalar T) const // ideal gas
{
    return 0;
}

template <class Specie>
inline Foam::scalar Foam::PengRobinson<Specie>::S(scalar p, scalar T) const // ideal gas
{
    return -this->R() * log(p / Pstd);
}

template <class Specie>
inline Foam::scalar Foam::PengRobinson<Specie>::dSdP(scalar p, scalar T) const // ideal gas
{
    return -this->R() / p;
}

template <class Specie>
inline Foam::scalar Foam::PengRobinson<Specie>::dSdT(scalar p, scalar T) const // ideal gas
{
    return 0;
}

template <class Specie>
inline Foam::scalar Foam::PengRobinson<Specie>::psi(
    scalar p,
    scalar T) const
{
    const scalar Z = this->Z(p, T);

    return 1.0 / (Z * this->R() * T);
}

template <class Specie>
inline Foam::scalar Foam::PengRobinson<Specie>::CpMCv //Todo
    (
        scalar p,
        scalar T) const
{
    const scalar Tr = T / Tc_;
    const scalar a = 0.45724 * sqr(RR * Tc_) / Pc_;
    const scalar b = 0.07780 * RR * Tc_ / Pc_;
    const scalar kappa = 0.37464 + 1.54226 * omega_ - 0.26992 * sqr(omega_);
    const scalar alpha = sqr(1 + kappa * (1 - sqrt(Tr)));

    const scalar A = alpha * a * p / sqr(RR * T);
    const scalar B = b * p / (RR * T);

    const scalar Z = this->Z(p, T);

    const scalar ap = kappa * a * (kappa / Tc_ - (1 + kappa) / sqrt(T * Tc_));
    const scalar M = (sqr(Z) + 2 * B * Z - sqr(B)) / (Z - B);
    const scalar N = ap * B / (b * RR);

    return this->R() * sqr(M - N) / (sqr(M) - 2 * A * (Z + B));
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template <class Specie>
inline void Foam::PengRobinson<Specie>::operator+=(
    const PengRobinson<Specie> &pg)
{
    scalar Y1 = this->Y();  //YP: Y() is mole number of component in mixture
    Specie::operator+=(pg); //YP: is it mean that sum all Y?  pg means all species Y, so after this, Y() means sumY

    if (mag(this->Y()) > small)
    {
        Y1 /= this->Y();                      // Y1 means mole fraction
        const scalar Y2 = pg.Y() / this->Y(); // pg.Y() means all species corresponding value; Y2 means others mole fraction

        Tc_ = Y1 * Tc_ + Y2 * pg.Tc_;
        Vc_ = Y1 * Vc_ + Y2 * pg.Vc_;
        Zc_ = Y1 * Zc_ + Y2 * pg.Zc_;
        Pc_ = RR * Zc_ * Tc_ / Vc_;
        omega_ = Y1 * omega_ + Y2 * pg.omega_;
    }
}

template <class Specie>
inline void Foam::PengRobinson<Specie>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template <class Specie>
Foam::PengRobinson<Specie> Foam::operator+(
    const PengRobinson<Specie> &pg1,
    const PengRobinson<Specie> &pg2)
{
    Specie sp(
        static_cast<const Specie &>(pg1) + static_cast<const Specie &>(pg2));

    if (mag(sp.Y()) < small)
    {
        return PengRobinson<Specie>(
            sp,
            pg1.Tc_,
            pg1.Vc_,
            pg1.Zc_,
            pg1.Pc_,
            pg1.omega_

        );
    }
    else
    {
        const scalar Y1 = pg1.Y() / sp.Y();
        const scalar Y2 = pg2.Y() / sp.Y();

        const scalar Tc = Y1 * pg1.Tc_ + Y2 * pg2.Tc_;
        const scalar Vc = Y1 * pg1.Vc_ + Y2 * pg2.Vc_;
        const scalar Zc = Y1 * pg1.Zc_ + Y2 * pg2.Zc_;

        return PengRobinson<Specie>(
            sp,
            Tc,
            Vc,
            Zc,
            RR * Zc * Tc / Vc,
            Y1 * pg1.omega_ + Y2 * pg2.omega_);
    }
}

template <class Specie>
Foam::PengRobinson<Specie> Foam::operator*(
    const scalar s,
    const PengRobinson<Specie> &pg)
{
    return PengRobinson<Specie>(
        s * static_cast<const Specie &>(pg),
        pg.Tc_,
        pg.Vc_,
        pg.Zc_,
        pg.Pc_,
        pg.omega_);
}

template <class Specie>
Foam::PengRobinson<Specie> Foam::operator==(
    const PengRobinson<Specie> &pg1,
    const PengRobinson<Specie> &pg2)
{
    Specie sp(
        static_cast<const Specie &>(pg1) == static_cast<const Specie &>(pg2));

    const scalar Y1 = pg1.Y() / sp.Y();
    const scalar Y2 = pg2.Y() / sp.Y();

    const scalar Tc = Y2 * pg2.Tc_ - Y1 * pg1.Tc_;
    const scalar Vc = Y2 * pg2.Vc_ - Y1 * pg1.Vc_;
    const scalar Zc = Y2 * pg2.Zc_ - Y1 * pg1.Zc_;

    return PengRobinson<Specie>(
        sp,
        Tc,
        Vc,
        Zc,
        Y2 * pg2.omega_ - Y1 * pg1.omega_);
}

//temp E
/*
template <class Specie>
inline Foam::scalar Foam::PengRobinson<Specie>::Hig_phase(
    scalar T) const
{
    return RR * 1.0e-03 * T * (T * (T * (T * (T * Hig_phase_[4] / 5.0 + Hig_phase_[3] / 4.0) + Hig_phase_[2] / 3.0) + Hig_phase_[1] / 2.0) + Hig_phase_[0]);
}

template <class Specie>
inline Foam::scalar Foam::PengRobinson<Specie>::Hig2_phase(
    scalar T) const
{
    //return   Hig2_phase_[0]*T +Hig2_phase_[1]*Hig2_phase_[2]*(2.0/(exp(2.0*Hig2_phase_[2]/T)-1.0))+ 2.0*Hig2_phase_[3]*Hig2_phase_[4]*(1.0/(exp(2.0*Hig2_phase_[4]/T)+1.0));
    return Hig2_phase_[0] + Hig2_phase_[1] * T + Hig2_phase_[2] * Hig2_phase_[3] * ((exp(2.0 * Hig2_phase_[3] / T) + 1) / (exp(2.0 * Hig2_phase_[3] / T) - 1.0)) + Hig2_phase_[4] * Hig2_phase_[5] * ((exp(2.0 * Hig2_phase_[5] / T) - 1.0) / (exp(2.0 * Hig2_phase_[5] / T) + 1.0));
}
*/
// ************************************************************************* //
