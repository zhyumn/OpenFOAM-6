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

#include "VLE.H"
#include <tuple>
#include "IOstreams.H"

// template <class ThermoMixture>
// inline tuple<scalar,scalar,scalar,scalar> Foam::VLE<ThermoMixture>::PF_fields_lv(scalar p, scalar T , scalarList &Y_liq, scalarList &Y_gas) const{
//     autoPtr<solution> sol(TPn_flash_New_TPD_Tudisco(p, T));
//     scalarList &X_liq = sol().X_liq();
//     scalarList &X_gas = sol().X_gas();
//     scalar Wliq = this->W(sol().X_liq());
//     scalar Wgas = this->W(sol().X_gas());
//     scalar rhoLiq = this->ThermoMixture::rho(p, T, sol().X_liq(), 0);
//     scalar rhoGas = this->ThermoMixture::rho(p, T, sol().X_gas(), 1);
//     scalar ha_liq,ha_gas;

//     scalar temp;

//     forAll(Y_liq,i){
//         if(rhoLiq => rhoGas){
//             Y_liq[i] = this->operator[](i).W()*1E-3*X_liq[i]/Wliq;
//             Y_gas[i] = this->operator[](i).W()*1E-3*X_gas[i]/Wgas;
//         }
//         else if(rhoLiq < rhoGas){
//             Y_liq[i] = this->operator[](i).W()*1E-3*X_gas[i]/Wgas;
//             Y_gas[i] = this->operator[](i).W()*1E-3*X_liq[i]/Wliq;
//         }
//     }

//     if(rhoLiq >= rhoGas){
//         ha_liq = this->ThermoMixture::Ha(p, T,sol().X_liq(), 0);
//         ha_gas = this->ThermoMixture::Ha(p, T,sol().X_gas(), 1);
//     }
//     else{
//         ha_liq = this->ThermoMixture::Ha(p, T,sol().X_gas(), 0);
//         ha_gas = this->ThermoMixture::Ha(p, T,sol().X_liq(), 1);
//     }

//     if(rhoLiq < rhoGas){
//         temp = rhoLiq;
//         rhoLiq = rhoGas;
//         rhoGas = temp;
//     }

//     return make_tuple(rhoLiq,rhoGas,ha_liq,ha_gas);
// }

template <class ThermoMixture>
inline tuple<scalar,scalar> Foam::VLE<ThermoMixture>::rho_lv(scalar p, scalar T) const{
    autoPtr<solution> sol(TPn_flash_New_TPD_Tudisco(p, T));
    scalar rhoLiq = this->ThermoMixture::rho(p, T, sol().X_liq(), 0);
    scalar rhoGas = this->ThermoMixture::rho(p, T, sol().X_gas(), 1);
    scalar temp;
    if(rhoLiq < rhoGas){
        temp = rhoLiq;
        rhoLiq = rhoGas;
        rhoGas = temp;
    }
    return make_tuple(rhoLiq,rhoGas);
}

template <class ThermoMixture>
inline void Foam::VLE<ThermoMixture>::Y_lv(scalar p, scalar T, scalarList &Y_liq, scalarList &Y_gas) const{
    autoPtr<solution> sol(TPn_flash_New_TPD_Tudisco(p, T));
    scalarList &X_liq = sol().X_liq();
    scalarList &X_gas = sol().X_gas();
    scalar Wliq = this->W(sol().X_liq());
    scalar Wgas = this->W(sol().X_gas());
    scalar rhoLiq = this->ThermoMixture::rho(p, T, sol().X_liq(), 0);
    scalar rhoGas = this->ThermoMixture::rho(p, T, sol().X_gas(), 1);

    // autoPtr<List<scalarList>> Yliq(new List<scalarList>(this->X_.size()));
    // autoPtr<List<scalarList>> Ygas(new List<scalarList>(this->X_.size()));
    // List<scalarList> &Y_liq = Yliq();
    // List<scalarList> &Y_gas = Ygas();
    // if(X_gas[1]<X_gas[0]){
    //     std::cout<<"y line"<<std::endl;
    //     std::cout<<rhoLiq<<"\t"<<rhoGas<<"\t"<<X_gas[0]<<"\t"<<X_gas[1]<<std::endl;
    // }

    forAll(Y_liq,i){
        if(rhoLiq >= rhoGas){
            Y_liq[i] = this->operator[](i).W()*1E-3*X_liq[i]/Wliq;
            Y_gas[i] = this->operator[](i).W()*1E-3*X_gas[i]/Wgas;
        }
        else if(rhoLiq < rhoGas){
            Y_liq[i] = this->operator[](i).W()*1E-3*X_gas[i]/Wgas;
            Y_gas[i] = this->operator[](i).W()*1E-3*X_liq[i]/Wliq;
        }
    }
}

template <class ThermoMixture>
inline tuple<scalar,scalar> Foam::VLE<ThermoMixture>::h_lv(scalar p, scalar T) const{
    autoPtr<solution> sol(TPn_flash_New_TPD_Tudisco(p, T));
    scalar rhoLiq = this->ThermoMixture::rho(p, T, sol().X_liq(), 0);
    scalar rhoGas = this->ThermoMixture::rho(p, T, sol().X_gas(), 1);
    scalar ha_liq,ha_gas;
    if(rhoLiq >= rhoGas){
        ha_liq = this->ThermoMixture::Hs(p, T,sol().X_liq(), 0);
        ha_gas = this->ThermoMixture::Hs(p, T,sol().X_gas(), 1);
    }
    else{
        ha_liq = this->ThermoMixture::Hs(p, T,sol().X_gas(), 0);
        ha_gas = this->ThermoMixture::Hs(p, T,sol().X_liq(), 1);
    }
    return make_tuple(ha_liq,ha_gas);
}