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

\*---------------------------------------------------------------------------*/

#include "ISATmultithermo.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template <class ThermoMixture, template <class> class Type>
inline Foam::species::ISATmultithermo<ThermoMixture, Type>::ISATmultithermo(
    const ThermoMixture& sp)
    : ThermoMixture(sp)
{
}
template <class ThermoMixture, template <class> class Type>
inline Foam::species::ISATmultithermo<ThermoMixture, Type>::ISATmultithermo(
    const word& name,
    PtrList<SingleThermoType>& speciesData,
    const speciesTable& specieNames,
    const dictionary& thermoDict)
    :
    ThermoMixture(name, speciesData, specieNames, thermoDict),
    ISATDict_
    (
        IFstream
        (
            fileName(thermoDict.lookup("ISATDictFile")).expand()
        )()
    ),

    Tpsivffun_HP(this),
    mufun(this, &ThermoMixture::mu),
    alphahfun(this, &ThermoMixture::alphah),

    Tpsivftree_HP(speciesData.size() + 2, 3, Tpsivffun_HP, "Tpsivftree_HP",ISATDict_),
    mutree(speciesData.size() + 1, 1, mufun, "mutree",ISATDict_),
    alphahtree(speciesData.size() + 1, 1, alphahfun, "alphahtree",ISATDict_)

{
    mutree.relepsilon() = 1;//2e-3;
    alphahtree.relepsilon() = 0;// 2e-3;
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar Foam::species::ISATmultithermo<ThermoMixture, Type>::T(
    scalar f,
    scalar p,
    scalar T0,
    scalar(ISATmultithermo<ThermoMixture, Type>::* F)(const scalar, const scalar) const,
    scalar(ISATmultithermo<ThermoMixture, Type>::* dFdT)(const scalar, const scalar)
    const,
    scalar(ISATmultithermo<ThermoMixture, Type>::* limit)(const scalar) const) const
{

    if (T0 < 0)
    {
        FatalErrorInFunction
            << "Negative initial temperature T0: " << T0
            << abort(FatalError);
    }
    //static int call_n=0; 
    //static int inter_n=0;
    scalar T1 = T0;
    scalar F1 = (this->*F)(p, T1), F2 = F1;
    scalar T2 = (this->*limit)(T0 + 1);
    // if(F1>f)
    //     T2=(this->*limit)(T0-30);
     //scalar F1=(this->*F)(p, T1);
     //scalar F2=(this->*F)(p, T2);
     //scalar Tnew = T0;
    scalar Ttol = T0 * tol_;
    scalar dF = 0, dfdt = 0;
    int iter = 0;
    scalar relax = 1;
    do
    {
        if (fabs(F2) > fabs(F1))
            relax *= 0.5;
        F1 = F2;
        F2 = (this->*F)(p, T2);

        dF = (F2 - f);
        dfdt = (F2 - F1) / (T2 - T1);//(this->*dFdT)(p, Test);
        T1 = T2;
        if (mag(relax * dF / dfdt) > 30)
            T2 = (this->*limit)(T1 - sign(relax * dF / dfdt) * 30);
        else
            T2 = (this->*limit)(T1 - relax * dF / dfdt);
        //Tnew =
         //   (this->*limit)(Test - ((this->*F)(p, Test) - f) / (this->*dFdT)(p, Test));

        if (T2 == 0)
        {
            FatalErrorInFunction << setprecision(50)
                << "Maximum number of iterations exceeded: " << maxIter_
                << " f: " << f
                << " p: " << p
                << " T0: " << T0
                << "  Tg: " << T2
                << " X: " << this->X_
                << abort(FatalError);
        }
        if (iter++ > maxIter_)
        {
            FatalErrorInFunction << setprecision(50)
                << "Maximum number of iterations exceeded: " << maxIter_
                << " f: " << f
                << " p: " << p
                << " T0: " << T0
                << "  Tg: " << T2
                << " X: " << this->X_
                << abort(FatalError);
        }

    } while (mag(T2 - T1) > Ttol && mag(f - F2) > 0.0001);
    //call_n++;
    //inter_n+=iter;
    //if(call_n%10000==0)
    //Info<<"-------------------call N="<<call_n<<",average iter="<<inter_n*1.0/call_n<<endl;
    return T2;

}

template <class ThermoMixture, template <class> class Type>
inline std::tuple<Foam::scalar, Foam::scalar> Foam::species::ISATmultithermo<ThermoMixture, Type>::TrhoEs(
    scalar e,
    scalar rho,
    scalar p0,
    scalar T0,
    scalar (ISATmultithermo<ThermoMixture, Type>::*Fe)(const scalar, const scalar) const,
    scalar (ISATmultithermo<ThermoMixture, Type>::*Frho)(const scalar, const scalar)
        const) const
{
    if (T0 < 0)
    {
        FatalErrorInFunction
            << "Negative initial temperature T0: " << T0
            << abort(FatalError);
    }
    //static int call_n=0;
    //static int inter_n=0;
    scalar dT = 1e-2;
    scalar dp = 1;
    scalar T1 = T0;
    scalar p1 = p0;
    scalar e11 = (this->*Fe)(p1, T1), e21, e12;
    scalar rho11 = (this->*Frho)(p1, T1), rho21, rho12;
    scalar T2 = T1 + 1, p2 = p1 + 10;  //(this->*limit)(T0+1);
    scalar Ttol = T0 * tol_, ptol = p0 * tol_;
    scalar de = 0, drho = 0, dedt = 0, dedp = 0, drhodt = 0, drhodp = 0;
    scalar de_old=1e10,drho_old=1e10;
    scalar dete;
    int iter = 0;
    scalar relax = 1;
    //Info<<endl;
    do
    {
        //Info<<T1<<","<<p1<<endl;
        dp = max(mag(dp), 1e-4)*sign(dp);
        dT = max(mag(dT), 1e-5)*sign(dT);
        e12 = (this->*Fe)(p1, T1 + dT);
        e21 = (this->*Fe)(p1 + dp, T1);
        rho12 = (this->*Frho)(p1, T1 + dT);
        rho21 = (this->*Frho)(p1 + dp, T1);
        de_old=de;
        drho_old=drho;
        de = (e11 - e);
        drho = (rho11 - rho);

        if(mag(de)>mag(de_old)||mag(drho)>mag(drho_old))
        {
            relax*=0.5;
        }
        dedt = (e12 - e11) / dT;  //(this->*dFdT)(p, Test);
        dedp = (e21 - e11) / dp;
        drhodt = (rho12 - rho11) / dT;
        drhodp = (rho21 - rho11) / dp;

        dete = dedt * drhodp - dedp * drhodt;
        /*if(fabs(dete)<1e-8)
        {
            Info<<"dedp="<<dedp<<",drhodp="<<drhodp<<",dp="<<dp<<endl;
            dT = -drho  / drhodt * relax;
            dp = -de  / dedp * relax;

        }
        else{
            */
        dT = -(drhodp * de - dedp * drho) / dete * relax;
        dp = -(dedt * drho - drhodt * de) / dete * relax;
        T1 = T1 + dT;
        p1 = p1 + dp;

        //}


        e11 = (this->*Fe)(p1, T1);
        rho11 = (this->*Frho)(p1, T1);

        if (iter++ > 100*maxIter_)
        {
            FatalErrorInFunction << setprecision(50)
                                 << "Maximum number of iterations exceeded: " << maxIter_
                                 << " e: " << e
                                 << " rho: " << rho
                                 << " T0: " << T0
                                 << "  p0: " << p0
                                 << " T1: " << T1
                                 << "  p1: " << p1
                                 << " X: " << this->X_
                                 << abort(FatalError);
        }

    } while (mag(dT) > Ttol || (mag(e11 - e) > 1e-5 || mag(rho11 - rho) > 1e-7|| mag(dT) > 1e-5|| mag(dp) > 1e-4 ));
    //call_n++;
    //inter_n+=iter;
    //if(call_n%10000==0)
    //Info<<"-------------------call N="<<call_n<<",average iter="<<inter_n*1.0/call_n<<endl;
    return std::make_tuple(p1, T1);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class ThermoMixture, template <class> class Type>
inline Foam::species::ISATmultithermo<ThermoMixture, Type>::ISATmultithermo(
    const word& name,
    const ISATmultithermo& st)
    : ThermoMixture(name, st)
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class ThermoMixture, template <class> class Type>
inline Foam::word
Foam::species::ISATmultithermo<ThermoMixture, Type>::heName()
{
    return Type<ISATmultithermo<ThermoMixture, Type>>::energyName();
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::Cv(const scalar p, const scalar T) const
{
    return this->Cp(p, T) - this->CpMCv(p, T);
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::Cpv(const scalar p, const scalar T) const
{
    return Type<ISATmultithermo<ThermoMixture, Type>>::Cpv(*this, p, T);
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::gamma(const scalar p, const scalar T) const
{
    const scalar Cp = this->Cp(p, T);
    return Cp / (Cp - this->CpMCv(p, T));
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::CpByCpv(
    const scalar p,
    const scalar T) const
{
    return Type<ISATmultithermo<ThermoMixture, Type>>::CpByCpv(*this, p, T);
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::HE(const scalar p, const scalar T) const
{
    return Type<ISATmultithermo<ThermoMixture, Type>>::HE(*this, p, T);
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::Es(const scalar p, const scalar T) const
{
    return this->Hs(p, T) - p / this->rho(p, T);
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::Ea(const scalar p, const scalar T) const
{
    return this->Ha(p, T) - p / this->rho(p, T);
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::G(const scalar p, const scalar T) const
{
    return this->Ha(p, T) - T * this->S(p, T);
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::A(const scalar p, const scalar T) const
{
    return this->Ea(p, T) - T * this->S(p, T);
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::cp(const scalar p, const scalar T) const
{
    return this->Cp(p, T) * this->W();
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::ha(const scalar p, const scalar T) const
{
    return this->Ha(p, T) * this->W();
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::hs(const scalar p, const scalar T) const
{
    return this->Hs(p, T) * this->W();
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::hc() const
{
    return this->Hc() * this->W();
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::s(const scalar p, const scalar T) const
{
    return this->S(p, T) * this->W();
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::he(const scalar p, const scalar T) const
{
    return this->HE(p, T) * this->W();
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::cv(const scalar p, const scalar T) const
{
    return this->Cv(p, T) * this->W();
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::es(const scalar p, const scalar T) const
{
    return this->Es(p, T) * this->W();
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::ea(const scalar p, const scalar T) const
{
    return this->Ea(p, T) * this->W();
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::g(const scalar p, const scalar T) const
{
    return this->G(p, T) * this->W();
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::a(const scalar p, const scalar T) const
{
    return this->A(p, T) * this->W();
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::K(const scalar p, const scalar T) const
{
    scalar arg = -1; //this->Y()*this->G(Pstd, T)/(RR*T);

    if (arg < 600)
    {
        return exp(arg);
    }
    else
    {
        return rootVGreat;
    }
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::Kp(const scalar p, const scalar T) const
{
    return K(p, T);
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::Kc(const scalar p, const scalar T) const
{
    const scalar nm = 1; //this->Y()/this->W();

    if (equal(nm, small))
    {
        return Kp(p, T);
    }
    else
    {
        return Kp(p, T) * pow(Pstd / (RR * T), nm);
    }
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar Foam::species::ISATmultithermo<ThermoMixture, Type>::Kx(
    const scalar p,
    const scalar T) const
{
    const scalar nm = this->Y() / this->W();

    if (equal(nm, small))
    {
        return Kp(p, T);
    }
    else
    {
        return Kp(p, T) * pow(Pstd / p, nm);
    }
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar Foam::species::ISATmultithermo<ThermoMixture, Type>::Kn(
    const scalar p,
    const scalar T,
    const scalar n) const
{
    const scalar nm = this->Y() / this->W();

    if (equal(nm, small))
    {
        return Kp(p, T);
    }
    else
    {
        return Kp(p, T) * pow(n * Pstd / p, nm);
    }
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar Foam::species::ISATmultithermo<ThermoMixture, Type>::THE(
    const scalar he,
    const scalar p,
    const scalar T0) const
{
    return Type<ISATmultithermo<ThermoMixture, Type>>::THE(*this, he, p, T0);
}
//This one !!!

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar Foam::species::ISATmultithermo<ThermoMixture, Type>::THs(
    const scalar hs,
    const scalar p,
    const scalar T0) const
{

    double sum = 0;
    //scalarList in(this->X().size() + 2), out(1);
    scalarList in(this->X().size() + 1), out(1);
    for (int i = 0;i < this->X().size() - 1;i++)
    {
        in[i] = this->X_[i];
        sum += this->X_[i];
    }
    in[this->X().size() - 1] = hs;
    in[this->X().size()] = p;
    //in[this->X().size() + 1] = T0;

    if (1 - sum < 1e-5)
    {
        //Ttree.notCall++;
        return T(hs, p, T0, &ISATmultithermo<ThermoMixture, Type>::Hs, &ISATmultithermo<ThermoMixture, Type>::Cp, &ISATmultithermo<ThermoMixture, Type>::limit);
    }
    //Hnfun.value(in, out);
    //Info << "THS!!" << endl;
    return T(hs, p, T0, &ISATmultithermo<ThermoMixture, Type>::Hs, &ISATmultithermo<ThermoMixture, Type>::Cp, &ISATmultithermo<ThermoMixture, Type>::limit);
    //Ttree.call(in, out);
    //Ttree.call(in, out, T0);
    //return out[0];
}

template <class ThermoMixture, template <class> class Type>
inline std::tuple<Foam::scalar, Foam::scalar> Foam::species::ISATmultithermo<ThermoMixture, Type>::TrhoE(
    const scalar e,
    const scalar rho,
    const scalar p0,
    const scalar T0) const
{
    return TrhoEs(e, rho, p0, T0, &ISATmultithermo<ThermoMixture, Type>::Es, &ISATmultithermo<ThermoMixture, Type>::rho);
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar Foam::species::ISATmultithermo<ThermoMixture, Type>::THa(
    const scalar ha,
    const scalar p,
    const scalar T0) const
{
    return T(
        ha,
        p,
        T0,
        &ISATmultithermo<ThermoMixture, Type>::Ha,
        &ISATmultithermo<ThermoMixture, Type>::Cp,
        &ISATmultithermo<ThermoMixture, Type>::limit);
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar Foam::species::ISATmultithermo<ThermoMixture, Type>::TEs(
    const scalar es,
    const scalar p,
    const scalar T0) const
{
    return T(
        es,
        p,
        T0,
        &ISATmultithermo<ThermoMixture, Type>::Es,
        &ISATmultithermo<ThermoMixture, Type>::Cv,
        &ISATmultithermo<ThermoMixture, Type>::limit);
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar Foam::species::ISATmultithermo<ThermoMixture, Type>::TEa(
    const scalar ea,
    const scalar p,
    const scalar T0) const
{
    return T(
        ea,
        p,
        T0,
        &ISATmultithermo<ThermoMixture, Type>::Ea,
        &ISATmultithermo<ThermoMixture, Type>::Cv,
        &ISATmultithermo<ThermoMixture, Type>::limit);
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::dKcdTbyKc(
    const scalar p,
    const scalar T) const
{
    const scalar nm = 1; //this->Y()/this->W();

    if (equal(nm, small))
    {
        return -this->dGdT(Pstd, T); //*this->Y()/RR;
    }
    else
    {
        return -(nm / T + this->dGdT(Pstd, T)); //*this->Y()/RR);
    }
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::dcpdT(const scalar p, const scalar T) const
{
    return this->dCpdT(p, T) * this->W();
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::psi(const scalar p, const scalar T) const
{

    double sum = 0;
    scalarList in(this->X().size() + 1), out(1);
    for (int i = 0;i < this->X().size() - 1;i++)
    {
        in[i] = this->X_[i];
        sum += this->X_[i];
    }
    in[this->X().size() - 1] = p;
    in[this->X().size()] = T;

    if (1 - sum < 1e-5)
    {
        //psitree.notCall++;
        return this->ThermoMixture::psi(p, T);
    }
    // Info<<"psihaha"<<endl;

    //psitree.call(in, out);
    //psivftree.call(in, out);
    return this->ThermoMixture::psi(p, T);
    //return out[0];

}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::psi_a(const scalar h, const scalar p, const scalar T0) const
{
    scalar tt = T(h, p, T0, &ISATmultithermo<ThermoMixture, Type>::Hs, &ISATmultithermo<ThermoMixture, Type>::Cp, &ISATmultithermo<ThermoMixture, Type>::limit);
    return this->ThermoMixture::psi(p, tt);
}


template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::psi_HP(const scalar H, const scalar P, const scalar T0) const
{
    scalarList in(this->X().size() + 2), out(1);
    for (int i = 0;i < this->X().size();i++)
    {
        in[i] = this->X_[i];
    }
    in[this->X().size()] = H;
    in[this->X().size() + 1] = P;

    //Info << "psi_HP" << endl;

    //psitree_HP.call(in, out, T0);
    //psivftree.call(in, out);
    return out[0];
}


template <class ThermoMixture, template <class> class Type>
inline std::tuple<Foam::scalar, Foam::scalar>
Foam::species::ISATmultithermo<ThermoMixture, Type>::Tpsi_HP(const scalar H, const scalar P, const scalar T0) const
{
    scalarList in(this->X().size() + 2), out(2);
    for (int i = 0;i < this->X().size();i++)
    {
        in[i] = this->X_[i];
    }
    in[this->X().size()] = H;
    in[this->X().size() + 1] = P;

    //Info << "psi_HP" << endl;

    //Tpsitree_HP.call(in, out, T0);
    //psivftree.call(in, out);
    return std::make_tuple(out[0], out[1]);

}

template <class ThermoMixture, template <class> class Type>
inline std::tuple<Foam::scalar, Foam::scalar, Foam::scalar>
Foam::species::ISATmultithermo<ThermoMixture, Type>::Tpsivf_HP(const scalar H, const scalar P, const scalar T0) const
{
    scalarList in(this->X().size() + 2), out(3);
    for (int i = 0;i < this->X().size();i++)
    {
        in[i] = this->X_[i];
    }
    in[this->X().size()] = H;
    in[this->X().size() + 1] = P;

    //Info << "psi_HP" << endl;

    Tpsivftree_HP.call(in, out, T0);
    //psivftree.call(in, out);
    return std::make_tuple(out[0], out[1], out[2]);
}

template <class ThermoMixture, template <class> class Type>
inline std::tuple<Foam::scalar, Foam::scalar, Foam::scalar>
Foam::species::ISATmultithermo<ThermoMixture, Type>::TPvf_Erho(const scalar E, const scalar rho, const scalar T0,const scalar P0) const
{
    scalarList in(this->X().size() + 2), out(3);
    for (int i = 0;i < this->X().size();i++)
    {
        in[i] = this->X_[i];
    }
    in[this->X().size()] = E;
    in[this->X().size() + 1] = rho;

    //Info << "psi_HP" << endl;

    Tpsivftree_HP.call(in, out, T0);
    //psivftree.call(in, out);
    return std::make_tuple(out[0], out[1], out[2]);
}


template <class ThermoMixture, template <class> class Type>
inline void Foam::species::ISATmultithermo<ThermoMixture, Type>::newTimeStep()
{
    Tpsivftree_HP.newTimeStep();
}
template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::mu(const scalar p, const scalar T) const
{

    double sum = 0;
    scalarList in(this->X().size() + 1), out(1);
    for (int i = 0;i < this->X().size() - 1;i++)
    {
        in[i] = this->X_[i];
        sum += this->X_[i];
    }
    in[this->X().size() - 1] = p;
    in[this->X().size()] = T;

    if (1 - sum < 1e-5)
    {
        mutree.notCall++;
        return this->ThermoMixture::mu(p, T);
    }
    // Info<<"psihaha"<<endl;

    mutree.call(in, out);
    return out[0];
}
template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::alphah(const scalar p, const scalar T) const
{

    double sum = 0;
    scalarList in(this->X().size() + 1), out(1);
    for (int i = 0;i < this->X().size() - 1;i++)
    {
        in[i] = this->X_[i];
        sum += this->X_[i];
    }
    in[this->X().size() - 1] = p;
    in[this->X().size()] = T;

    if (1 - sum < 1e-5)
    {
        alphahtree.notCall++;
        return this->ThermoMixture::alphah(p, T);
    }
    // Info<<"psihaha"<<endl;

    alphahtree.call(in, out);
    return out[0];

}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::vaporfra(const scalar p, const scalar T) const
{

    double sum = 0;
    scalarList in(this->X().size() + 1), out(1);
    for (int i = 0;i < this->X().size() - 1;i++)
    {
        in[i] = this->X_[i];
        sum += this->X_[i];
    }
    in[this->X().size() - 1] = p;
    in[this->X().size()] = T;

    if (1 - sum < 1e-5)
    {
        //vaporfratree.notCall++;
        return this->ThermoMixture::vaporfra(p, T);
    }
    // Info<<"psihaha"<<endl;

    //vaporfratree.call(in, out);
    //psivftree.call(in, out);
    //Info<<in<<endl;
    return this->ThermoMixture::vaporfra(p, T);
    return out[0];

}
template <class ThermoMixture, template <class> class Type>
inline void Foam::species::ISATmultithermo<ThermoMixture, Type>::psi_vf(const scalar p, const scalar T, scalar& psi_out, scalar& vf_out) const
{
    double sum = 0;
    scalarList in(this->X().size() + 1), out(1);
    for (int i = 0;i < this->X().size() - 1;i++)
    {
        in[i] = this->X_[i];
        sum += this->X_[i];
    }
    in[this->X().size() - 1] = p;
    in[this->X().size()] = T;

    if (1 - sum < 1e-5)
    {
        //psivftree.notCall++;
        psi_out = this->ThermoMixture::psi(p, T);
        vf_out = this->ThermoMixture::vaporfra(p, T);
    }

    // Info<<"psihaha"<<endl;

    //psivftree.call(in, out);
    psi_out = this->ThermoMixture::psi(p, T);
    vf_out = this->ThermoMixture::vaporfra(p, T);
    psi_out = out[0];
    vf_out = out[1];
    return;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template <class ThermoMixture, template <class> class Type>
inline void Foam::species::ISATmultithermo<ThermoMixture, Type>::operator+=(
    const ISATmultithermo<ThermoMixture, Type>& st)
{
    ThermoMixture::operator+=(st);
}

template <class ThermoMixture, template <class> class Type>
inline void Foam::species::ISATmultithermo<ThermoMixture, Type>::operator*=(const scalar s)
{
    ThermoMixture::operator*=(s);
}

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template <class ThermoMixture, template <class> class Type>
inline Foam::species::ISATmultithermo<ThermoMixture, Type> Foam::species::operator+(
    const ISATmultithermo<ThermoMixture, Type>& st1,
    const ISATmultithermo<ThermoMixture, Type>& st2)
{
    return ISATmultithermo<ThermoMixture, Type>(
        static_cast<const ThermoMixture&>(st1) + static_cast<const ThermoMixture&>(st2));
}

template <class ThermoMixture, template <class> class Type>
inline Foam::species::ISATmultithermo<ThermoMixture, Type> Foam::species::operator*(
    const scalar s,
    const ISATmultithermo<ThermoMixture, Type>& st)
{
    return ISATmultithermo<ThermoMixture, Type>(
        s * static_cast<const ThermoMixture&>(st));
}

template <class ThermoMixture, template <class> class Type>
inline Foam::species::ISATmultithermo<ThermoMixture, Type> Foam::species::operator==(
    const ISATmultithermo<ThermoMixture, Type>& st1,
    const ISATmultithermo<ThermoMixture, Type>& st2)
{
    return ISATmultithermo<ThermoMixture, Type>(
        static_cast<const ThermoMixture&>(st1) == static_cast<const ThermoMixture&>(st2));
}

// ************************************************************************* //
