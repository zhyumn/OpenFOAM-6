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

#include "ISATVLEheRhoThermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template <class BasicPsiThermo, class MixtureType>
void Foam::ISATVLEheRhoThermo<BasicPsiThermo, MixtureType>::calculate_init()
{
    scalarField &hCells = this->he_;
    const scalarField &pCells = this->p_;

    scalarField &TCells = this->T_.primitiveFieldRef();
    scalarField &psiCells = this->psi_.primitiveFieldRef();
    scalarField &muCells = this->mu_.primitiveFieldRef();
    scalarField &alphaCells = this->alpha_.primitiveFieldRef();
    scalarField &vaporfracCells = this->vaporfrac_.primitiveFieldRef();
    scalarField &soundspeedCells = this->soundspeed_.primitiveFieldRef();
    scalarField &rhoCells = this->rho_.primitiveFieldRef();
    scalarField &kappaCells = this->kappa_.primitiveFieldRef();
    scalarField &rho_G_Cells = this->rho_G_.primitiveFieldRef();
    //scalarField &ZCells = this->Z_.primitiveFieldRef();
    scalar p_temp, T_temp, he_temp;
    scalarList Y_temp(MixtureType::Y().size());

    //PF fields
    scalarField &rho_lCells = this->rho_l_.primitiveFieldRef();
    scalarField &rho_vCells = this->rho_v_.primitiveFieldRef();
    scalarField &h_lCells = this->h_l_.primitiveFieldRef();
    scalarField &h_vCells = this->h_v_.primitiveFieldRef();
    scalarList Y_l_temp(MixtureType::Y().size());
    scalarList Y_v_temp(MixtureType::Y().size());
    

    forAll(TCells, celli)
    {
        const typename MixtureType::thermoType &mixture_ =
            this->cellMixture(celli);

        //std::tie(rhoCells[celli], vaporfracCells[celli], soundspeedCells[celli]) = mixture_.rhovfc_ISAT(pCells[celli], TCells[celli]);
        //std::tie(he_temp, rhoCells[celli], vaporfracCells[celli], soundspeedCells[celli], rho_G_Cells[celli]) = mixture_.Erhovfc_G_rhoY_ISAT(pCells[celli], TCells[celli], Y_temp);
        std::tie(hCells[celli], rhoCells[celli], vaporfracCells[celli], soundspeedCells[celli], rho_G_Cells[celli]) = mixture_.Erhovfc_G_rhoY_ISAT(pCells[celli], TCells[celli], Y_temp);
        //ZCells[celli] = mixture_.Z(pCells[celli], TCells[celli]);
        for (int i = 0; i < Y_temp.size(); i++)
        {
            Y_G_List_[i][celli] = Y_temp[i];
        }
        psiCells[celli] = rhoCells[celli] / pCells[celli];

        //PF terms update
        std::tie(rho_lCells[celli],rho_vCells[celli]) = mixture_.rho_lv(pCells[celli], TCells[celli]);
        std::tie(h_lCells[celli],h_vCells[celli]) = mixture_.h_lv(pCells[celli], TCells[celli]);
        mixture_.Y_lv(pCells[celli], TCells[celli], Y_l_temp, Y_v_temp);
        for (int i = 0; i < Y_temp.size(); i++)
        {
            Y_l_[i][celli] = Y_l_temp[i];
            Y_v_[i][celli] = Y_v_temp[i];
        }
    }

    if (!inviscid_)
    {
        forAll(TCells, celli)
        {
            const typename MixtureType::thermoType &mixture_ =
                this->cellMixture(celli);
            tie(kappaCells[celli], muCells[celli]) = mixture_.kappa_mu_opt(pCells[celli], TCells[celli], rhoCells[celli]);
            //muCells[celli] = mixture_.mu_opt(pCells[celli], TCells[celli], rhoCells[celli]);
            //kappaCells[celli] = mixture_.kappa_opt(pCells[celli], TCells[celli], rhoCells[celli]);
        }

        forAll(heList_, i)
        {
            forAll(TCells, celli)
            {
                const typename MixtureType::thermoType &mixture_ =
                    this->cellMixture(celli);
                heList_[i].primitiveFieldRef()[celli] = this->speciesData()[i].Hs(pCells[celli], TCells[celli]);

                Dimix_[i].primitiveFieldRef()[celli] = mixture_.Dimix_opt(pCells[celli], TCells[celli], i);
            }
        }
    }

    volScalarField::Boundary &pBf =
        this->p_.boundaryFieldRef();

    volScalarField::Boundary &TBf =
        this->T_.boundaryFieldRef();

    volScalarField::Boundary &psiBf =
        this->psi_.boundaryFieldRef();

    volScalarField::Boundary &heBf =
        this->he().boundaryFieldRef();

    volScalarField::Boundary &muBf =
        this->mu_.boundaryFieldRef();

    volScalarField::Boundary &alphaBf =
        this->alpha_.boundaryFieldRef();

    volScalarField::Boundary &vaporfracBf =
        this->vaporfrac_.boundaryFieldRef();

    volScalarField::Boundary &soundspeedBf =
        this->soundspeed_.boundaryFieldRef();

    volScalarField::Boundary &rhoBf =
        this->rho_.boundaryFieldRef();

    volScalarField::Boundary &kappaBf =
        this->kappa_.boundaryFieldRef();

    volScalarField::Boundary &rho_G_Bf =
        this->rho_G_.boundaryFieldRef();
    
    //PF boundary terms
    volScalarField::Boundary &rho_l_Bf =
        this->rho_l_.boundaryFieldRef();
    
    volScalarField::Boundary &rho_v_Bf =
        this->rho_v_.boundaryFieldRef();

    volScalarField::Boundary &h_l_Bf =
        this->h_l_.boundaryFieldRef();

    volScalarField::Boundary &h_v_Bf =
        this->h_v_.boundaryFieldRef();

    //volScalarField::Boundary &ZBf =
    //    this->Z_.boundaryFieldRef();

    // volScalarField::Boundary& entropyBf =
    // this->entropy_.boundaryFieldRef();

    forAll(this->T_.boundaryField(), patchi)
    {
        fvPatchScalarField &pp = pBf[patchi];
        fvPatchScalarField &pT = TBf[patchi];
        fvPatchScalarField &ppsi = psiBf[patchi];
        fvPatchScalarField &phe = heBf[patchi];
        fvPatchScalarField &pmu = muBf[patchi];
        fvPatchScalarField &palpha = alphaBf[patchi];
        fvPatchScalarField &pvaporfrac = vaporfracBf[patchi];
        fvPatchScalarField &psoundspeed = soundspeedBf[patchi];
        fvPatchScalarField &prho = rhoBf[patchi];
        fvPatchScalarField &pkappa = kappaBf[patchi];
        fvPatchScalarField &prho_G = rho_G_Bf[patchi];
        //fvPatchScalarField &pZ = ZBf[patchi];

        //PF fields
        fvPatchScalarField &prho_l = rho_l_Bf[patchi];
        fvPatchScalarField &prho_v = rho_v_Bf[patchi];
        fvPatchScalarField &ph_l = h_l_Bf[patchi];
        fvPatchScalarField &ph_v = h_v_Bf[patchi];

        // fvPatchScalarField& pentropy = entropyBf[patchi];

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType &mixture_ =
                    this->patchFaceMixture(patchi, facei);
                //std::tie(prho[facei], pvaporfrac[facei], psoundspeed[facei]) = mixture_.rhovfc_ISAT(pp[facei], pT[facei]);
                std::tie(he_temp, prho[facei], pvaporfrac[facei], psoundspeed[facei], prho_G[facei]) = mixture_.Erhovfc_G_rhoY_ISAT(pp[facei], pT[facei], Y_temp);
                //pZ[facei] = mixture_.Z(pp[facei], pT[facei]);
                for (int i = 0; i < Y_temp.size(); i++)
                {
                    Y_G_List_[i].boundaryFieldRef()[patchi][facei] = Y_temp[i];
                }
                ppsi[facei] = prho[facei] / pp[facei];

                //PF terms update
                std::tie(prho_l[facei],prho_v[facei]) = mixture_.rho_lv(pp[facei], pT[facei]);
                std::tie(ph_l[facei],ph_v[facei]) = mixture_.h_lv(pp[facei], pT[facei]);
                mixture_.Y_lv(pp[facei], pT[facei], Y_l_temp, Y_v_temp);
                for (int i = 0; i < Y_temp.size(); i++)
                {
                    Y_l_[i].boundaryFieldRef()[patchi][facei] = Y_l_temp[i];
                    Y_v_[i].boundaryFieldRef()[patchi][facei] = Y_v_temp[i];
                }
            }
        }
        else
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType &mixture_ =
                    this->patchFaceMixture(patchi, facei);
                std::tie(he_temp, prho[facei], pvaporfrac[facei], psoundspeed[facei], prho_G[facei]) = mixture_.Erhovfc_G_rhoY_ISAT(pp[facei], pT[facei], Y_temp);
                //pZ[facei] = mixture_.Z(pp[facei], pT[facei]);
                for (int i = 0; i < Y_temp.size(); i++)
                {
                    Y_G_List_[i].boundaryFieldRef()[patchi][facei] = Y_temp[i];
                }
                ppsi[facei] = prho[facei] / pp[facei];

                //PF terms update
                std::tie(prho_l[facei],prho_v[facei]) = mixture_.rho_lv(pp[facei], pT[facei]);
                std::tie(ph_l[facei],ph_v[facei]) = mixture_.h_lv(pp[facei], pT[facei]);
                mixture_.Y_lv(pp[facei], pT[facei], Y_l_temp, Y_v_temp);
                for (int i = 0; i < Y_temp.size(); i++)
                {
                    Y_l_[i].boundaryFieldRef()[patchi][facei] = Y_l_temp[i];
                    Y_v_[i].boundaryFieldRef()[patchi][facei] = Y_v_temp[i];
                }
            }
        }

        if (!inviscid_)
        {
            if (pT.fixesValue())
            {
                forAll(pT, facei)
                {
                    const typename MixtureType::thermoType &mixture_ =
                        this->patchFaceMixture(patchi, facei);
                    tie(pkappa[facei], pmu[facei]) = mixture_.kappa_mu_opt(pp[facei], pT[facei], prho[facei]);
                    //pmu[facei] = mixture_.mu_opt(pp[facei], pT[facei], prho[facei]);
                    //pkappa[facei] = mixture_.kappa_opt(pp[facei], pT[facei], prho[facei]);

                    // pmu[facei] = mixture_.mu_rho(pp[facei], pT[facei], prho[facei]);
                    // pkappa[facei] = mixture_.kappa_rho(pp[facei], pT[facei]);
                    // palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
                }
            }
            else
            {
                forAll(pT, facei)
                {
                    const typename MixtureType::thermoType &mixture_ =
                        this->patchFaceMixture(patchi, facei);
                    tie(pkappa[facei], pmu[facei]) = mixture_.kappa_mu_opt(pp[facei], pT[facei], prho[facei]);
                    //pmu[facei] = mixture_.mu_opt(pp[facei], pT[facei], prho[facei]);
                    //pkappa[facei] = mixture_.kappa_opt(pp[facei], pT[facei], prho[facei]);
                    // pmu[facei] = mixture_.mu_rho(pp[facei], pT[facei], prho[facei]);
                    // pkappa[facei] = mixture_.kappa_rho(pp[facei], pT[facei]);
                    // palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
                }
            }
        }
    }
    // this->rho_= this->psi_* this->p_;

    if (!inviscid_)
    {
        forAll(heList_, i)
        {
            volScalarField::Boundary &DimixBf =
                this->Dimix_[i].boundaryFieldRef();

            volScalarField::Boundary &heListBf =
                this->heList_[i].boundaryFieldRef();

            forAll(this->T_.boundaryField(), patchi)
            {
                fvPatchScalarField &pDimix = DimixBf[patchi];
                fvPatchScalarField &pheList = heListBf[patchi];
                fvPatchScalarField &pT = TBf[patchi];
                fvPatchScalarField &pp = pBf[patchi];
                if (pT.fixesValue())
                {
                    forAll(pT, facei)
                    {
                        const typename MixtureType::thermoType &mixture_ =
                            this->patchFaceMixture(patchi, facei);

                        pDimix[facei] = mixture_.Dimix_opt(pp[facei], pT[facei], i);

                        pheList[facei] = this->speciesData()[i].Hs(pp[facei], pT[facei]);
                    }
                }
                else
                {
                    forAll(pT, facei)
                    {
                        const typename MixtureType::thermoType &mixture_ =
                            this->patchFaceMixture(patchi, facei);
                        pDimix[facei] = mixture_.Dimix_opt(pp[facei], pT[facei], i);
                        pheList[facei] = this->speciesData()[i].Hs(pp[facei], pT[facei]); // TODO check Hs
                    }
                }
            }
        }
    }
    this->mute_show();
}

template <class BasicPsiThermo, class MixtureType>
void Foam::ISATVLEheRhoThermo<BasicPsiThermo, MixtureType>::calculate()
{
    // static scalar maxvlll = 0;
    static bool boundary_flag = false;
    scalarField &hCells = this->he_;
    scalarField &pCells = this->p_;
    static int debug_flag = 0;
    debug_flag++;
    scalarField &TCells = this->T_.primitiveFieldRef();
    scalarField &psiCells = this->psi_.primitiveFieldRef();
    scalarField &muCells = this->mu_.primitiveFieldRef();
    scalarField &alphaCells = this->alpha_.primitiveFieldRef();
    scalarField &vaporfracCells = this->vaporfrac_.primitiveFieldRef();
    scalarField &soundspeedCells = this->soundspeed_.primitiveFieldRef();
    scalarField &rhoCells = this->rho_.primitiveFieldRef();
    scalarField &kappaCells = this->kappa_.primitiveFieldRef();
    scalarField &rho_G_Cells = this->rho_G_.primitiveFieldRef();
    //scalarField &ZCells = this->Z_.primitiveFieldRef();
    scalar tempT, tempP, maxdT = 0, maxdP = 0, tempmu, temppsi, tempHe;
    static scalar maxdmu = 0;
    static int timeflag = 0;
    scalarList Y_temp(MixtureType::Y().size());
    scalar VLEtime = 0;

    //PF fields
    scalarField &rho_lCells = this->rho_l_.primitiveFieldRef();
    scalarField &rho_vCells = this->rho_v_.primitiveFieldRef();
    scalarField &h_lCells = this->h_l_.primitiveFieldRef();
    scalarField &h_vCells = this->h_v_.primitiveFieldRef();
    scalarList Y_l_temp(MixtureType::Y().size());
    scalarList Y_v_temp(MixtureType::Y().size());

    if (!boundary_flag)
    {
        this->newTimeStep();
        clockTime_.timeIncrement();
        if (DF_)
        {
            do
            {
                forAll(TCells, celli)
                {
                    const typename MixtureType::thermoType &mixture_ = this->cellMixture(celli);

                    // std::tie(TCells[celli], psiCells[celli], vaporfracCells[celli], soundspeedCells[celli]) = mixture_.Tpsivfc_XHP(hCells[celli] + pCells[celli] / rhoCells[celli], pCells[celli], TCells[celli]);
                    std::tie(TCells[celli], hCells[celli], vaporfracCells[celli], soundspeedCells[celli], rho_G_[celli]) = mixture_.THvfc_G_rhoY_XrhoP(rhoCells[celli], pCells[celli], TCells[celli], Y_temp);
                    //std::tie(hCells[celli], rhoCells[celli], vaporfracCells[celli], soundspeedCells[celli], rho_G_Cells[celli]) = mixture_.Erhovfc_G_rhoY_ISAT(pCells[celli], TCells[celli], Y_temp);

                    hCells[celli] -= pCells[celli] / rhoCells[celli];

                    for (int i = 0; i < Y_temp.size(); i++)
                    {
                        Y_G_List_[i][celli] = Y_temp[i];
                    }
                    // autoPtr<typename MixtureType::thermoType::solution> sol;
                    // std::tie(TCells[celli], temppsi, vaporfracCells[celli], soundspeedCells[celli], sol) = mixture_.Tpsivfcsol_XHP(hCells[celli] + pCells[celli] / rhoCells[celli], pCells[celli], TCells[celli]);
                    // rhoCells[celli] = psiCells[celli] * pCells[celli];
                    psiCells[celli] = rhoCells[celli] / pCells[celli];

                    // PF terms
                    std::tie(rho_lCells[celli],rho_vCells[celli]) = mixture_.rho_lv(pCells[celli], TCells[celli]);
                    std::tie(h_lCells[celli],h_vCells[celli]) = mixture_.h_lv(pCells[celli], TCells[celli]);
                    mixture_.Y_lv(pCells[celli], TCells[celli], Y_l_temp, Y_v_temp);
                    for (int i = 0; i < Y_temp.size(); i++)
                    {
                        Y_l_[i][celli] = Y_l_temp[i];
                        Y_v_[i][celli] = Y_v_temp[i];
                    }
                }
            } while (this->newLoop());
        }
        else
        {
            do
            {
                forAll(TCells, celli)
                {
                    const typename MixtureType::thermoType &mixture_ = this->cellMixture(celli);

                    std::tie(TCells[celli], pCells[celli], vaporfracCells[celli], soundspeedCells[celli], rho_G_[celli]) = mixture_.TPvfc_G_rhoY_XErho(hCells[celli], rhoCells[celli], TCells[celli], pCells[celli], Y_temp);

                    for (int i = 0; i < Y_temp.size(); i++)
                    {
                        Y_G_List_[i][celli] = Y_temp[i];
                    }
                    psiCells[celli] = rhoCells[celli] / pCells[celli];

                    // PF terms
                    std::tie(rho_lCells[celli],rho_vCells[celli]) = mixture_.rho_lv(pCells[celli], TCells[celli]);
                    std::tie(h_lCells[celli],h_vCells[celli]) = mixture_.h_lv(pCells[celli], TCells[celli]);
                    mixture_.Y_lv(pCells[celli], TCells[celli], Y_l_temp, Y_v_temp);
                    for (int i = 0; i < Y_temp.size(); i++)
                    {
                        Y_l_[i][celli] = Y_l_temp[i];
                        Y_v_[i][celli] = Y_v_temp[i];
                    }
                }
            } while (this->newLoop());
        }
        
        if (ISATlog_)
        {
            VLEtime += clockTime_.timeIncrement();
            cpuISAT_VLE_()
                << this->time().timeOutputValue()
                << ",    " << VLEtime << endl;
        }

        if (!inviscid_)
        {
            forAll(TCells, celli)
            {
                const typename MixtureType::thermoType &mixture_ =
                    this->cellMixture(celli);
                tie(kappaCells[celli], muCells[celli]) = mixture_.kappa_mu_opt(pCells[celli], TCells[celli], rhoCells[celli]);
                //muCells[celli] = mixture_.mu_opt(pCells[celli], TCells[celli], rhoCells[celli]);
                //kappaCells[celli] = mixture_.kappa_opt(pCells[celli], TCells[celli], rhoCells[celli]);
                // alphaCells[celli] = mixture_.alphah(pCells[celli], TCells[celli]);//kappaCells[celli] / mixture_.Cp(pCells[celli], TCells[celli]);
            }

            forAll(heList_, i)
            {
                forAll(TCells, celli)
                {
                    const typename MixtureType::thermoType &mixture_ =
                        this->cellMixture(celli);
                    heList_[i].primitiveFieldRef()[celli] = this->speciesData()[i].Hs(pCells[celli], TCells[celli]);

                    Dimix_[i].primitiveFieldRef()[celli] = mixture_.Dimix_opt(pCells[celli], TCells[celli], i);
                }
            }
        }
    }
    else
    {

        volScalarField::Boundary &pBf =
            this->p_.boundaryFieldRef();

        volScalarField::Boundary &rhoBf =
            this->rho_.boundaryFieldRef();

        volScalarField::Boundary &TBf =
            this->T_.boundaryFieldRef();

        volScalarField::Boundary &psiBf =
            this->psi_.boundaryFieldRef();

        volScalarField::Boundary &heBf =
            this->he().boundaryFieldRef();

        volScalarField::Boundary &muBf =
            this->mu_.boundaryFieldRef();

        volScalarField::Boundary &alphaBf =
            this->alpha_.boundaryFieldRef();

        volScalarField::Boundary &vaporfracBf =
            this->vaporfrac_.boundaryFieldRef();

        volScalarField::Boundary &soundspeedBf =
            this->soundspeed_.boundaryFieldRef();

        volScalarField::Boundary &kappaBf =
            this->kappa_.boundaryFieldRef();

        volScalarField::Boundary &rho_G_Bf =
            this->rho_G_.boundaryFieldRef();

        //PF boundary terms
        volScalarField::Boundary &rho_l_Bf =
            this->rho_l_.boundaryFieldRef();
        
        volScalarField::Boundary &rho_v_Bf =
            this->rho_v_.boundaryFieldRef();

        volScalarField::Boundary &h_l_Bf =
            this->h_l_.boundaryFieldRef();

        volScalarField::Boundary &h_v_Bf =
            this->h_v_.boundaryFieldRef();

        //volScalarField::Boundary &ZBf =
        //    this->Z_.boundaryFieldRef();

        forAll(this->T_.boundaryField(), patchi)
        {
            fvPatchScalarField &pp = pBf[patchi];
            fvPatchScalarField &prho = rhoBf[patchi];
            fvPatchScalarField &pT = TBf[patchi];
            fvPatchScalarField &ppsi = psiBf[patchi];
            fvPatchScalarField &phe = heBf[patchi];
            fvPatchScalarField &pmu = muBf[patchi];
            fvPatchScalarField &palpha = alphaBf[patchi];
            fvPatchScalarField &pvaporfrac = vaporfracBf[patchi];
            fvPatchScalarField &psoundspeed = soundspeedBf[patchi];
            fvPatchScalarField &pkappa = kappaBf[patchi];
            fvPatchScalarField &prho_G = rho_G_Bf[patchi];
            //fvPatchScalarField &pZ = ZBf[patchi];

            //PF fields
            fvPatchScalarField &prho_l = rho_l_Bf[patchi];
            fvPatchScalarField &prho_v = rho_v_Bf[patchi];
            fvPatchScalarField &ph_l = h_l_Bf[patchi];
            fvPatchScalarField &ph_v = h_v_Bf[patchi];

            if (pT.fixesValue())
            {
                forAll(pT, facei)
                {
                    const typename MixtureType::thermoType &mixture_ =
                        this->patchFaceMixture(patchi, facei);
                    //std::tie(prho[facei], pvaporfrac[facei], psoundspeed[facei]) = mixture_.rhovfc_ISAT(pp[facei], pT[facei]);
                    std::tie(phe[facei], prho[facei], pvaporfrac[facei], psoundspeed[facei], prho_G[facei]) = mixture_.Erhovfc_G_rhoY_ISAT(pp[facei], pT[facei], Y_temp);
                    //pZ[facei] = mixture_.Z(pp[facei], pT[facei]);
                    for (int i = 0; i < Y_temp.size(); i++)
                    {
                        Y_G_List_[i].boundaryFieldRef()[patchi][facei] = Y_temp[i];
                    }

                    //phe[facei] = mixture_.HE(pp[facei], pT[facei]);
                    ppsi[facei] = prho[facei] / pp[facei];

                    //PF terms update
                    std::tie(prho_l[facei],prho_v[facei]) = mixture_.rho_lv(pp[facei], pT[facei]);
                    std::tie(ph_l[facei],ph_v[facei]) = mixture_.h_lv(pp[facei], pT[facei]);
                    mixture_.Y_lv(pp[facei], pT[facei], Y_l_temp, Y_v_temp);
                    for (int i = 0; i < Y_temp.size(); i++)
                    {
                        Y_l_[i].boundaryFieldRef()[patchi][facei] = Y_l_temp[i];
                        Y_v_[i].boundaryFieldRef()[patchi][facei] = Y_v_temp[i];
                    }
                }
            }
            else
            {
                forAll(pT, facei)
                {
                    const typename MixtureType::thermoType &mixture_ =
                        this->patchFaceMixture(patchi, facei);

                    //std::tie(phe[facei], prho[facei], pvaporfrac[facei], psoundspeed[facei]) = mixture_.Erhovfc_ISAT(pp[facei], pT[facei]);
                    std::tie(phe[facei], prho[facei], pvaporfrac[facei], psoundspeed[facei], prho_G[facei]) = mixture_.Erhovfc_G_rhoY_ISAT(pp[facei], pT[facei], Y_temp);
                    //pZ[facei] = mixture_.Z(pp[facei], pT[facei]);
                    for (int i = 0; i < Y_temp.size(); i++)
                    {
                        Y_G_List_[i].boundaryFieldRef()[patchi][facei] = Y_temp[i];
                    }

                    ppsi[facei] = prho[facei] / pp[facei];

                    //PF terms update
                    std::tie(prho_l[facei],prho_v[facei]) = mixture_.rho_lv(pp[facei], pT[facei]);
                    std::tie(ph_l[facei],ph_v[facei]) = mixture_.h_lv(pp[facei], pT[facei]);
                    mixture_.Y_lv(pp[facei], pT[facei], Y_l_temp, Y_v_temp);
                    for (int i = 0; i < Y_temp.size(); i++)
                    {
                        Y_l_[i].boundaryFieldRef()[patchi][facei] = Y_l_temp[i];
                        Y_v_[i].boundaryFieldRef()[patchi][facei] = Y_v_temp[i];
                    }
                }
            }

            if (!inviscid_)
            {
                forAll(pT, facei)
                {
                    const typename MixtureType::thermoType &mixture_ =
                        this->patchFaceMixture(patchi, facei);
                    tie(pkappa[facei], pmu[facei]) = mixture_.kappa_mu_opt(pp[facei], pT[facei], prho[facei]);
                    //pmu[facei] = mixture_.mu_opt(pp[facei], pT[facei], prho[facei]);
                    //pkappa[facei] = mixture_.kappa_opt(pp[facei], pT[facei], prho[facei]);
                }
            }
        }

        if (!inviscid_)
        {
            forAll(heList_, i)
            {
                volScalarField::Boundary &DimixBf =
                    this->Dimix_[i].boundaryFieldRef();

                volScalarField::Boundary &heListBf =
                    this->heList_[i].boundaryFieldRef();
                forAll(this->T_.boundaryField(), patchi)
                {
                    fvPatchScalarField &pDimix = DimixBf[patchi];
                    fvPatchScalarField &pheList = heListBf[patchi];
                    fvPatchScalarField &pT = TBf[patchi];
                    fvPatchScalarField &pp = pBf[patchi];
                    if (pT.fixesValue())
                    {
                        forAll(pT, facei)
                        {
                            const typename MixtureType::thermoType &mixture_ =
                                this->patchFaceMixture(patchi, facei);

                            pDimix[facei] = mixture_.Dimix_opt(pp[facei], pT[facei], i);
                            pheList[facei] = this->speciesData()[i].Hs(pp[facei], pT[facei]);
                        }
                    }
                    else
                    {
                        forAll(pT, facei)
                        {
                            const typename MixtureType::thermoType &mixture_ =
                                this->patchFaceMixture(patchi, facei);
                            pDimix[facei] = mixture_.Dimix_opt(pp[facei], pT[facei], i);

                            pheList[facei] = this->speciesData()[i].Hs(pp[facei], pT[facei]);
                        }
                    }
                }
            }
        }
    }
    boundary_flag = !boundary_flag;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class BasicPsiThermo, class MixtureType>
Foam::ISATVLEheRhoThermo<BasicPsiThermo, MixtureType>::ISATVLEheRhoThermo(
    const fvMesh &mesh,
    const word &phaseName)
    : heThermo<BasicPsiThermo, MixtureType>(mesh, phaseName),
      vaporfrac_(
          IOobject(
              "vaporfrac",
              mesh.time().timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::AUTO_WRITE),
          mesh,
          dimless),
      soundspeed_(
          IOobject(
              "soundspeed",
              mesh.time().timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::AUTO_WRITE),
          mesh,
          dimVelocity),
      kappa_(
          IOobject(
              "thermo:kappa",
              mesh.time().timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          mesh,
          dimEnergy / (dimTime * dimLength * dimTemperature)),
      rho_G_(
          IOobject(
              "thermo:rho_G",
              mesh.time().timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          mesh,
          dimDensity),
      rho_l_(
          IOobject(
              "thermo:rho_l",
              mesh.time().timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          mesh,
          dimDensity),
      rho_v_(
          IOobject(
              "thermo:rho_v",
              mesh.time().timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          mesh,
          dimDensity),
      h_l_(
          IOobject(
              "thermo:h_l",
              mesh.time().timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          mesh,
          this->he_.dimensions()),
      h_v_(
          IOobject(
              "thermo:h_v",
              mesh.time().timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          mesh,
          this->he_.dimensions()),
      /*       Z_(
          IOobject(
              "thermo:Z",
              mesh.time().timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          mesh,
          dimless), */
      Dimix_(MixtureType::Y().size()),
      heList_(MixtureType::Y().size()),
      WList_(MixtureType::Y().size()),
      Y_G_List_(MixtureType::Y().size()),
      Y_l_(MixtureType::Y().size()),
      Y_v_(MixtureType::Y().size()),
      inviscid_(false)
{
    IOdictionary thermoDict(
        IOobject(
            this->phasePropertyName(this->dictName, phaseName),
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false));
    inviscid_ = thermoDict.lookupOrDefault<bool>("inviscid", false);
    nloop = thermoDict.lookupOrDefault<label>("nloop", 1);
    DF_ = thermoDict.lookupOrDefault<bool>("doubleFlux", true);
    ISATlog_ = thermoDict.lookupOrDefault<bool>("ISATlog", false);
    //noVLE_ = thermoDict.lookupOrDefault<bool>("noVLE", false);
    //MixtureType::thermoType::noVLE = noVLE_;
    // FatalErrorInFunction
    //     << "inviscid_:" <<inviscid_
    //     << exit(FatalError);
    if (ISATlog_)
    {
        cpuISAT_VLE_ = logFile("VLEtime", mesh);
    }

    forAll(Dimix_, i)
    {
        Dimix_.set(
            i,
            new volScalarField(
                IOobject(
                    this->phasePropertyName("thermo:Dimix:") + MixtureType::species()[i],
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE),
                mesh,
                dimensionSet(0, 2, -1, 0, 0)));
        heList_.set(
            i,
            new volScalarField(

                IOobject(
                    "hei:" + MixtureType::species()[i],
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE),
                mesh,
                this->he_.dimensions()));
        WList_.set(
            i,
            new volScalarField(
                IOobject(
                    "W:" + MixtureType::species()[i],
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false),
                mesh,
                dimMass / dimMoles));
        Y_G_List_.set(
            i,
            new volScalarField(
                IOobject(
                    "Y_G:" + MixtureType::species()[i],
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE),
                mesh,
                dimless));
        Y_l_.set(
            i,
            new volScalarField(
                IOobject(
                    "Y_l:" + MixtureType::species()[i],
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE),
                mesh,
                dimless));
        Y_v_.set(
            i,
            new volScalarField(
                IOobject(
                    "Y_v:" + MixtureType::species()[i],
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE),
                mesh,
                dimless));
    }

    forAll(WList_, i)
    {
        forAll(WList_[i].primitiveFieldRef(), celli)
        {

            WList_[i].primitiveFieldRef()[celli] = this->speciesData()[i].W();
        }
    }
    calculate_init();

    // Switch on saving old time
    this->psi_.oldTime();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class BasicPsiThermo, class MixtureType>
Foam::ISATVLEheRhoThermo<BasicPsiThermo, MixtureType>::~ISATVLEheRhoThermo()
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class BasicPsiThermo, class MixtureType>
void Foam::ISATVLEheRhoThermo<BasicPsiThermo, MixtureType>::correct()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    // force the saving of the old-time values
    this->psi_.oldTime();

    calculate();

    if (debug)
    {
        Info << "    Finished" << endl;
    }
}

template <class BasicPsiThermo, class MixtureType>
inline Foam::autoPtr<Foam::OFstream>
Foam::ISATVLEheRhoThermo<BasicPsiThermo, MixtureType>::logFile(const word &name, const fvMesh &mesh) const
{
    mkDir(mesh.time().path() / "ISAT_VLE");
    return autoPtr<OFstream>(
        new OFstream(
            mesh.time().path() / "ISAT_VLE" / name));
}

// ************************************************************************* //
