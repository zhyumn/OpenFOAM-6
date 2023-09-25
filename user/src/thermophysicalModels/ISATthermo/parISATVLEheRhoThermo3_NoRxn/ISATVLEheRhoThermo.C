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
#include "SUPstream.H"
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template <class BasicPsiThermo, class MixtureType>
void Foam::ISATVLEheRhoThermo<BasicPsiThermo, MixtureType>::calculate_init()
{
    //const scalarField &hCells = this->he_;
    scalarField &hCells = this->he_;
    const scalarField &pCells = this->p_;

    scalarField &TCells = this->T_.primitiveFieldRef();
    scalarField &psiCells = this->psi_.primitiveFieldRef();
    scalarField &muCells = this->mu_.primitiveFieldRef();
    //scalarField &alphaCells = this->alpha_.primitiveFieldRef();
    scalarField &vaporfracCells = this->vaporfrac_.primitiveFieldRef();
    scalarField &soundspeedCells = this->soundspeed_.primitiveFieldRef();
    scalarField &rhoCells = this->rho_.primitiveFieldRef();
    scalarField &kappaCells = this->kappa_.primitiveFieldRef();

    //scalar p_temp, T_temp,
    //scalar he_temp;
    scalarList Y_temp(MixtureType::Y().size());
    //Pout << "!!!!!!!!test-1" << endl;
    forAll(TCells, celli)
    {
        const typename MixtureType::thermoType &mixture_ =
            this->cellMixture(celli);

        //std::tie(he_temp, rhoCells[celli], vaporfracCells[celli], soundspeedCells[celli]) = mixture_.Erhovfc_XPT(pCells[celli], TCells[celli]);
        std::tie(hCells[celli], rhoCells[celli], vaporfracCells[celli], soundspeedCells[celli]) = mixture_.Erhovfc_XPT(pCells[celli], TCells[celli]);
        psiCells[celli] = rhoCells[celli] / pCells[celli];
    }
    //Pout << "!!!!!!!!test-2" << endl;
    if (!inviscid_)
    {
        forAll(TCells, celli)
        {
            const typename MixtureType::thermoType &mixture_ =
                this->cellMixture(celli);
            tie(kappaCells[celli], muCells[celli]) = mixture_.kappa_mu_opt(pCells[celli], TCells[celli], rhoCells[celli]);
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

    //volScalarField::Boundary &alphaBf =
    //    this->alpha_.boundaryFieldRef();

    volScalarField::Boundary &vaporfracBf =
        this->vaporfrac_.boundaryFieldRef();

    volScalarField::Boundary &soundspeedBf =
        this->soundspeed_.boundaryFieldRef();

    volScalarField::Boundary &rhoBf =
        this->rho_.boundaryFieldRef();

    volScalarField::Boundary &kappaBf =
        this->kappa_.boundaryFieldRef();

    //volScalarField::Boundary &rho_G_Bf =
    //this->rho_G_.boundaryFieldRef();

    //volScalarField::Boundary &ZBf =
    //    this->Z_.boundaryFieldRef();

    // volScalarField::Boundary& entropyBf =
    // this->entropy_.boundaryFieldRef();
    //Pout << "!!!!!!!!test-2-1" << endl;
    forAll(this->T_.boundaryField(), patchi)
    {
        fvPatchScalarField &pp = pBf[patchi];
        fvPatchScalarField &pT = TBf[patchi];
        fvPatchScalarField &ppsi = psiBf[patchi];
        fvPatchScalarField &phe = heBf[patchi];
        fvPatchScalarField &pmu = muBf[patchi];
        //fvPatchScalarField &palpha = alphaBf[patchi];
        fvPatchScalarField &pvaporfrac = vaporfracBf[patchi];
        fvPatchScalarField &psoundspeed = soundspeedBf[patchi];
        fvPatchScalarField &prho = rhoBf[patchi];
        fvPatchScalarField &pkappa = kappaBf[patchi];
        //fvPatchScalarField &prho_G = rho_G_Bf[patchi];
        //fvPatchScalarField &pZ = ZBf[patchi];

        // fvPatchScalarField& pentropy = entropyBf[patchi];

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType &mixture_ =
                    this->patchFaceMixture(patchi, facei);

                //std::tie(he_temp, prho[facei], pvaporfrac[facei], psoundspeed[facei]) = mixture_.Erhovfc_XPT(pp[facei], pT[facei]);
                std::tie(phe[facei], prho[facei], pvaporfrac[facei], psoundspeed[facei]) = mixture_.Erhovfc_XPT(pp[facei], pT[facei]);

                ppsi[facei] = prho[facei] / pp[facei];
            }
        }
        else
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType &mixture_ =
                    this->patchFaceMixture(patchi, facei);
                //std::tie(he_temp, prho[facei], pvaporfrac[facei], psoundspeed[facei]) = mixture_.Erhovfc_XPT(pp[facei], pT[facei]);
                std::tie(phe[facei], prho[facei], pvaporfrac[facei], psoundspeed[facei]) = mixture_.Erhovfc_XPT(pp[facei], pT[facei]);

                ppsi[facei] = prho[facei] / pp[facei];
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
                }
            }
            else
            {
                forAll(pT, facei)
                {
                    const typename MixtureType::thermoType &mixture_ =
                        this->patchFaceMixture(patchi, facei);
                    tie(pkappa[facei], pmu[facei]) = mixture_.kappa_mu_opt(pp[facei], pT[facei], prho[facei]);
                }
            }
        }
    }
    // this->rho_= this->psi_* this->p_;
    //Pout << "!!!!!!!!test-2-2" << endl;
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
    //Pout << "!!!!!!!!test-2-3" << endl;
    this->mute_show();
    //Pout << "!!!!!!!!test-2-4" << endl;
}

void my_quicksort(labelField &arrays1, labelField &arrays2, int left, int right)
{
    if (left >= right)
    {
        return;
    }
    int i = left, j = right;
    label init1 = arrays1[i];
    label init2 = arrays2[i];

    while (j > i)
    {
        for (; arrays1[j] > init1 && j > i; j--)
            ;
        if (i < j)
        {
            arrays1[i] = arrays1[j];
            arrays2[i] = arrays2[j];
            i++;
        }
        for (; arrays1[i] < init1 && j > i; i++)
            ;
        if (i < j)
        {
            arrays1[j] = arrays1[i];
            arrays2[j] = arrays2[i];
            j--;
        }
    }
    arrays1[i] = init1;
    arrays2[i] = init2;
    my_quicksort(arrays1, arrays2, left, i - 1);
    my_quicksort(arrays1, arrays2, i + 1, right);
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
    //scalarField &alphaCells = this->alpha_.primitiveFieldRef();
    scalarField &vaporfracCells = this->vaporfrac_.primitiveFieldRef();
    scalarField &soundspeedCells = this->soundspeed_.primitiveFieldRef();
    scalarField &rhoCells = this->rho_.primitiveFieldRef();
    scalarField &rhodCells = this->rho_d.primitiveFieldRef();
    scalarField &kappaCells = this->kappa_.primitiveFieldRef();

    bool retreived_tmp;
    label N_add_tmp;
    //scalarField &rho_G_Cells = this->rho_G_.primitiveFieldRef();
    //scalarField &ZCells = this->Z_.primitiveFieldRef();
    //scalar tempT, tempP, maxdT = 0, maxdP = 0, tempmu, temppsi, tempHe;
    //static scalar maxdmu = 0;
    //static int timeflag = 0;
    scalarList Y_temp(MixtureType::Y().size());
    scalar VLEtime = 0;
    //Pout << "!!!!!!!!test-3" << endl;
    if (!boundary_flag)
    {
        this->newTimeStep();
        SUPstream::Sync();
        clockTime_.timeIncrement();
        //scalar starttime = MPI_Wtime();
        if (scheme_ == "doubleFlux")
        {
            do
            {
                forAll(TCells, celli)
                {
                    const typename MixtureType::thermoType &mixture_ = this->cellMixture(celli);

                    //std::tie(TCells[celli], hCells[celli], vaporfracCells[celli], soundspeedCells[celli]) = mixture_.THvfc_XrhoP(rhoCells[celli], pCells[celli], TCells[celli]);

                    //hCells[celli] -= pCells[celli] / rhoCells[celli];
                    std::tie(TCells[celli], hCells[celli], vaporfracCells[celli], soundspeedCells[celli]) = mixture_.TEvfc_XrhoP(rhoCells[celli], pCells[celli], TCells[celli]);
                    psiCells[celli] = rhoCells[celli] / pCells[celli];
                }
            } while (this->newLoop());
        }
        else if (scheme_ == "doubleFlux++")
        {
            do
            {
                forAll(TCells, celli)
                {
                    const typename MixtureType::thermoType &mixture_ = this->cellMixture(celli);

                    //std::tie(TCells[celli], hCells[celli], vaporfracCells[celli], soundspeedCells[celli]) = mixture_.THvfc_XrhoP(rhoCells[celli], pCells[celli], TCells[celli]);

                    //hCells[celli] -= pCells[celli] / rhoCells[celli];
                    std::tie(TCells[celli], rhodCells[celli], vaporfracCells[celli], soundspeedCells[celli]) = mixture_.Trhovfc_XEP(hCells[celli], pCells[celli], TCells[celli]);
                    psiCells[celli] = rhoCells[celli] / pCells[celli];
                }
            } while (this->newLoop());
        }
        else if (scheme_ == "conservativeFlux")
        {
            if (!loadBalance_)
            {
                do
                {
                    forAll(TCells, celli)
                    {
                        const typename MixtureType::thermoType &mixture_ = this->cellMixture(celli);

                        std::tie(TCells[celli], pCells[celli], vaporfracCells[celli], soundspeedCells[celli], retreived_tmp) = mixture_.TPvfc_XErho(hCells[celli], rhoCells[celli], TCells[celli], pCells[celli]);

                        psiCells[celli] = rhoCells[celli] / pCells[celli];
                    }
                } while (this->newLoop());
            }
            else
            {
                Nsend = 0;
                Nreceive = 0;
                spare = false;
                spare_cpu--;
                nfinished_block = 0;

                int batch_size = Batch_Size;
                int nloop = TCells.size() / batch_size;
                int last_size = TCells.size() % batch_size;
                if (last_size != 0)
                {
                    nloop++;
                }
                else
                {
                    last_size = batch_size;
                }

                int iter_link = head_link;

                for (int i = 0; i < nloop;)
                {
                    int size_i = batch_size;
                    if (iter_link == nloop - 1)
                    {
                        size_i = last_size;
                    }
                    int first_i = iter_link * batch_size;

                    N_add_tmp = size_i;
                    for (int j = 0; j < size_i; j++)
                    {
                        label celli = first_i + j;
                        const typename MixtureType::thermoType &mixture_ = this->cellMixture(celli);
                        std::tie(TCells[celli], pCells[celli], vaporfracCells[celli], soundspeedCells[celli], retreived_tmp) = mixture_.TPvfc_XErho(hCells[celli], rhoCells[celli], TCells[celli], pCells[celli]);

                        psiCells[celli] = rhoCells[celli] / pCells[celli];
                        N_add_tmp -= retreived_tmp;
                    }
                    N_add_[iter_link] = N_add_tmp;
                    nfinished_block++;
                    iter_link = link_[iter_link];
                    i++;

                    if (finished_head[rank] != -1)
                    {
                        int iter_start;
                        int iter;
                        int tail;

                        pl_finished[rank].lock();
                        iter_start = finished_head[rank];
                        iter = iter_start;
                        finished_head[rank] = -1;
                        pl_finished[rank].unlock();

                        tail = iter;
                        while (iter != -1)
                        {
                            for (int j = 0; j < ja[iter].size_i; j++)
                            {
                                label celli = ja[iter].first_i + j;
                                TCells[celli] = ja[iter].jobs[j].out[0];
                                pCells[celli] = ja[iter].jobs[j].out[1];
                                vaporfracCells[celli] = ja[iter].jobs[j].out[2];
                                soundspeedCells[celli] = ja[iter].jobs[j].out[3];
                                psiCells[celli] = rhoCells[celli] / pCells[celli];
                            }
                            tail = iter;
                            N_add_[ja[iter].N_batch] = ja[iter].N_add;
                            iter = ja[iter].next;
                            nfinished_block++;
                        }
                        l_empty.lock();
                        if (empty_tail == -1)
                        {
                            empty_tail = tail;
                            empty_head = iter_start;
                        }
                        else
                        {
                            ja[empty_tail].next = iter_start;
                            empty_tail = tail;
                        }
                        l_empty.unlock();
                    }

                    if (i < nloop && Nsend < Nplan && l_empty.try_lock())
                    {
                        int size_i = batch_size;
                        if (iter_link == nloop - 1)
                        {
                            size_i = last_size;
                        }
                        int first_i = iter_link * batch_size;

                        if (empty_head != -1)
                        {
                            int iter;
                            iter = empty_head;
                            empty_head = ja[empty_head].next;
                            if (empty_head == -1)
                            {
                                empty_tail = -1;
                            }
                            l_empty.unlock();

                            ja[iter].first_i = first_i;
                            ja[iter].size_i = size_i;

                            for (int j = 0; j < ja[iter].size_i; j++)
                            {
                                label celli = ja[iter].first_i + j;
                                const typename MixtureType::thermoType &mixture_ = this->cellMixture(celli);
                                ja[iter].jobs[j].p0 = pCells[celli];
                                ja[iter].jobs[j].T0 = TCells[celli];
                                for (int k = 0; k < mixture_.X().size(); k++)
                                    ja[iter].jobs[j].in[k] = mixture_.X()[k];
                                ja[iter].jobs[j].in[mixture_.X().size()] = hCells[celli];
                                ja[iter].jobs[j].in[mixture_.X().size() + 1] = rhoCells[celli];
                            }
                            Nsend++;
                            ja[iter].rank = rank;
                            ja[iter].N_batch = iter_link;

                            l_filled.lock();
                            if (filled_tail == -1)
                            {
                                ja[iter].next = -1;
                                filled_tail = iter;
                                filled_head = iter;
                            }
                            else
                            {
                                ja[iter].next = -1;
                                ja[filled_tail].next = iter;
                                filled_tail = iter;
                            }
                            l_filled.unlock();
                            iter_link = link_[iter_link];
                            i++;
                        }
                        else
                        {
                            l_empty.unlock();
                        }
                    }

                    if (i < nloop && spare_cpu > 0 && l_sender.try_lock())
                    {
                        for (int ii = 0; ii <= spare_cpu; ii++)
                        {
                            int size_i = batch_size;
                            if (iter_link == nloop - 1)
                            {
                                size_i = last_size;
                            }
                            int first_i = iter_link * batch_size;

                            if (empty_head != -1 && i < nloop)
                            {
                                int iter;
                                l_empty.lock();
                                if (empty_head != -1)
                                {
                                    iter = empty_head;
                                    empty_head = ja[empty_head].next;
                                    if (empty_head == -1)
                                    {
                                        empty_tail = -1;
                                    }
                                    l_empty.unlock();

                                    ja[iter].first_i = first_i;
                                    ja[iter].size_i = size_i;

                                    for (int j = 0; j < ja[iter].size_i; j++)
                                    {
                                        label celli = ja[iter].first_i + j;
                                        const typename MixtureType::thermoType &mixture_ = this->cellMixture(celli);
                                        ja[iter].jobs[j].p0 = pCells[celli];
                                        ja[iter].jobs[j].T0 = TCells[celli];
                                        for (int k = 0; k < mixture_.X().size(); k++)
                                            ja[iter].jobs[j].in[k] = mixture_.X()[k];
                                        ja[iter].jobs[j].in[mixture_.X().size()] = hCells[celli];
                                        ja[iter].jobs[j].in[mixture_.X().size() + 1] = rhoCells[celli];
                                    }
                                    Nsend++;
                                    ja[iter].rank = rank;
                                    ja[iter].N_batch = iter_link;

                                    l_filled.lock();
                                    if (filled_tail == -1)
                                    {
                                        ja[iter].next = -1;
                                        filled_tail = iter;
                                        filled_head = iter;
                                    }
                                    else
                                    {
                                        ja[iter].next = -1;
                                        ja[filled_tail].next = iter;
                                        filled_tail = iter;
                                    }
                                    l_filled.unlock();
                                    iter_link = link_[iter_link];
                                    i++;
                                }
                                else
                                {
                                    l_empty.unlock();
                                }
                            }
                            else
                                break;
                        }
                        l_sender.unlock();
                    }
                }
                //Pout << "2nd loop" << endl;
                while (1)
                {

                    if (l_receiver.try_lock())
                    {
                        int iter = -1;
                        if (filled_head != -1)
                        {
                            l_filled.lock();
                            iter = filled_head;

                            filled_head = ja[filled_head].next;
                            if (filled_head == -1)
                            {
                                filled_tail = -1;
                            }
                            l_filled.unlock();
                        }
                        l_receiver.unlock();
                        if (iter != -1)
                        {
                            if (spare)
                            {
                                spare = false;
                                spare_cpu--;
                            }
                            int N_add = ja[iter].size_i;
                            for (int j = 0; j < ja[iter].size_i; j++)
                            {
                                const typename MixtureType::thermoType &mixture_ = this->cellMixture(0);
                                std::tie(ja[iter].jobs[j].out[0], ja[iter].jobs[j].out[1], ja[iter].jobs[j].out[2], ja[iter].jobs[j].out[3], retreived_tmp) = mixture_.TPvfc_XErho(ja[iter].jobs[j].in, ja[iter].jobs[j].T0, ja[iter].jobs[j].p0);
                                ja[iter].jobs[j].retrieved = retreived_tmp;
                                N_add -= retreived_tmp;
                            }
                            ja[iter].N_add = N_add;
                            Nreceive++;
                            pl_finished[ja[iter].rank].lock();
                            ja[iter].next = finished_head[ja[iter].rank];
                            finished_head[ja[iter].rank] = iter;
                            pl_finished[ja[iter].rank].unlock();
                        }
                    }
                    if (finished_head[rank] != -1)
                    {
                        int iter_start;
                        int iter;
                        int tail;

                        pl_finished[rank].lock();
                        iter_start = finished_head[rank];
                        iter = iter_start;
                        finished_head[rank] = -1;
                        pl_finished[rank].unlock();

                        tail = iter;
                        while (iter != -1)
                        {
                            for (int j = 0; j < ja[iter].size_i; j++)
                            {
                                label celli = ja[iter].first_i + j;
                                TCells[celli] = ja[iter].jobs[j].out[0];
                                pCells[celli] = ja[iter].jobs[j].out[1];
                                vaporfracCells[celli] = ja[iter].jobs[j].out[2];
                                soundspeedCells[celli] = ja[iter].jobs[j].out[3];
                            }
                            N_add_[ja[iter].N_batch] = ja[iter].N_add;
                            tail = iter;
                            iter = ja[iter].next;
                            nfinished_block++;
                        }
                        l_empty.lock();

                        if (empty_tail == -1)
                        {
                            empty_tail = tail;
                            empty_head = iter_start;
                        }
                        else
                        {
                            ja[empty_tail].next = iter_start;
                            empty_tail = tail;
                        }
                        l_empty.unlock();
                    }
                    if (!spare)
                    {
                        spare = true;
                        spare_cpu++;
                    }
                    if (spare_cpu.load() == n_cpu && nfinished_block == nloop && filled_head == -1)
                    {
                        break;
                    }
                    /*                     if (filled_head != -1)
                    {
                        Pout << "1" << endl;
                    }
                    if (spare_cpu.load() != n_cpu)
                    {
                        Pout << "2" << endl;
                    }
                    if (nfinished_block != nloop)
                    {
                        Pout << "3" << endl;
                    } */
                }

                for (int iter = 0; iter < Nbatch; iter++)
                {
                    index_[iter] = iter;
                }
                //Pout << N_add_ << endl;
                my_quicksort(N_add_, index_, 0, Nbatch - 1);
                head_link = index_[0];
                label iter_ = head_link;
                for (int iter = 0; iter < Nbatch - 1; iter++)
                {
                    link_[iter_] = index_[iter + 1];
                    iter_ = link_[iter_];
                }
                link_[iter_] = -1;

                Nplan = Nsend - Nreceive;
                //Pout << link_ << endl;
                SUPstream::Sync();
            }
        }
        else if (scheme_ == "mix")
        {
            do
            {
                forAll(TCells, celli)
                {
                    const typename MixtureType::thermoType &mixture_ = this->cellMixture(celli);
                    if ((*FCcell)[celli] == 0)
                    {
                        //std::tie(TCells[celli], hCells[celli], vaporfracCells[celli], soundspeedCells[celli]) = mixture_.THvfc_XrhoP(rhoCells[celli], pCells[celli], TCells[celli]);
                        //hCells[celli] -= pCells[celli] / rhoCells[celli];

                        std::tie(TCells[celli], hCells[celli], vaporfracCells[celli], soundspeedCells[celli]) = mixture_.TEvfc_XrhoP(rhoCells[celli], pCells[celli], TCells[celli]);
                    }
                    else if ((*FCcell)[celli] == 1)
                    {
                        std::tie(TCells[celli], pCells[celli], vaporfracCells[celli], soundspeedCells[celli], retreived_tmp) = mixture_.TPvfc_XErho(hCells[celli], rhoCells[celli], TCells[celli], pCells[celli]);
                    }

                    psiCells[celli] = rhoCells[celli] / pCells[celli];
                }
            } while (this->newLoop());
        }

        //scalar endtime = MPI_Wtime();

        if (ISATlog_)
        {
            VLEtime += clockTime_.timeIncrement();
            cpuISAT_VLE_()
                << this->time().timeOutputValue()
                << ",    " << VLEtime << endl;
            tree_size_()
                << this->time().timeOutputValue()
                << ",    " << this->cellMixture(0).TPvfc_XErho_treesize() << ",    " << this->cellMixture(0).TPvfc_XErho_treesize_loc() << endl;

            /*             cpuISAT_VLE_()
                << this->time().timeOutputValue()
                << ",    " << endtime - starttime << endl; */
        }
        SUPstream::Sync();
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

        //volScalarField::Boundary &alphaBf =
        //    this->alpha_.boundaryFieldRef();

        volScalarField::Boundary &vaporfracBf =
            this->vaporfrac_.boundaryFieldRef();

        volScalarField::Boundary &soundspeedBf =
            this->soundspeed_.boundaryFieldRef();

        volScalarField::Boundary &kappaBf =
            this->kappa_.boundaryFieldRef();

        //volScalarField::Boundary &rho_G_Bf =
        //    this->rho_G_.boundaryFieldRef();

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
            //fvPatchScalarField &palpha = alphaBf[patchi];
            fvPatchScalarField &pvaporfrac = vaporfracBf[patchi];
            fvPatchScalarField &psoundspeed = soundspeedBf[patchi];
            fvPatchScalarField &pkappa = kappaBf[patchi];
            //fvPatchScalarField &prho_G = rho_G_Bf[patchi];
            //fvPatchScalarField &pZ = ZBf[patchi];

            if (pT.fixesValue())
            {
                forAll(pT, facei)
                {
                    const typename MixtureType::thermoType &mixture_ =
                        this->patchFaceMixture(patchi, facei);

                    std::tie(phe[facei], prho[facei], pvaporfrac[facei], psoundspeed[facei]) = mixture_.Erhovfc_XPT(pp[facei], pT[facei]);

                    ppsi[facei] = prho[facei] / pp[facei];
                }
            }
            else
            {
                forAll(pT, facei)
                {
                    const typename MixtureType::thermoType &mixture_ =
                        this->patchFaceMixture(patchi, facei);

                    std::tie(phe[facei], prho[facei], pvaporfrac[facei], psoundspeed[facei]) = mixture_.Erhovfc_XPT(pp[facei], pT[facei]);

                    ppsi[facei] = prho[facei] / pp[facei];
                }
            }

            if (!inviscid_)
            {
                forAll(pT, facei)
                {
                    const typename MixtureType::thermoType &mixture_ =
                        this->patchFaceMixture(patchi, facei);
                    tie(pkappa[facei], pmu[facei]) = mixture_.kappa_mu_opt(pp[facei], pT[facei], prho[facei]);
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

      rho_d(
          IOobject(
              "rho_d",
              mesh.time().timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          mesh,
          dimDensity),
      /*       rho_G_(
          IOobject(
              "thermo:rho_G",
              mesh.time().timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          mesh,
          dimDensity), */
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
      //Y_G_List_(MixtureType::Y().size()),
      inviscid_(false),
      n_block(SUPstream::node_manager.size * 2),
      Batch_Size(this->lookupOrDefault("batchSize", 100)),
      Nbatch(vaporfrac_.size() / Batch_Size + (vaporfrac_.size() % Batch_Size > 0 ? 1 : 0)),
      link_(Nbatch, 0),
      link_rev(Nbatch, 0),
      N_add_(Nbatch, 0),
      index_(Nbatch, 0),
      ja(SUPstream::node_manager, n_block, sizeof(jobArray) + sizeof(jobInput) * (Batch_Size - 1)),
      sspare_cpu(SUPstream::node_manager), spare_cpu(sspare_cpu().var),
      sfinished_head(SUPstream::node_manager, SUPstream::node_manager.size), finished_head(&sfinished_head()),
      l_empty(SUPstream::node_manager, SUPstream::Sync),
      l_filled(SUPstream::node_manager, SUPstream::Sync),
      l_finished(SUPstream::node_manager, SUPstream::Sync),
      sl_finished(SUPstream::node_manager, SUPstream::node_manager.size),
      pl_finished(sl_finished.ptr()),
      l_receiver(SUPstream::node_manager, SUPstream::Sync), l_sender(SUPstream::node_manager, SUPstream::Sync),
      sempty_head(SUPstream::node_manager), sempty_tail(SUPstream::node_manager),
      empty_head(sempty_head()), empty_tail(sempty_tail()),
      sfilled_head(SUPstream::node_manager), sfilled_tail(SUPstream::node_manager),
      filled_head(sfilled_head()), filled_tail(sfilled_tail()),
      rank(SUPstream::node_manager.rank)
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
    //DF_ = thermoDict.lookupOrDefault<bool>("doubleFlux", true);
    ISATlog_ = thermoDict.lookupOrDefault<bool>("ISATlog", false);
    scheme_ = word(thermoDict.lookup("scheme"));
    n_cpu = SUPstream::node_manager.size;

    loadBalance_ = thermoDict.lookupOrDefault<bool>("loadBalance", false);
    //noVLE_ = thermoDict.lookupOrDefault<bool>("noVLE", false);
    //MixtureType::thermoType::noVLE = noVLE_;
    // FatalErrorInFunction
    //     << "inviscid_:" <<inviscid_
    //     << exit(FatalError);
    FCcell = &mesh.objectRegistry::lookupObjectRef<volScalarField>("FCcell");
    if (ISATlog_)
    {
        cpuISAT_VLE_ = logFile("VLEtime", mesh);
        tree_size_ = logFile("TreeSize", mesh);
    }
    SUPstream::Sync();
    if (SUPstream::node_manager.rank == 0)
    {
        spare_cpu = n_cpu;
        for (int i = 0; i < n_block; i++)
        {
            ja[i].next = i + 1;
        }
        ja[n_block - 1].next = -1;
        empty_head = 0;
        empty_tail = n_block - 1;
        filled_head = -1;
        filled_tail = -1;
    }
    finished_head[rank] = -1;
    SUPstream::Sync();

    head_link = 0;
    //head_add = -1;
    for (int i = 0; i < Nbatch; i++)
    {
        link_[i] = i + 1;
        N_add_[i] = 0;
    }
    link_[link_.size() - 1] = -1;

    Nplan = 0;

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
