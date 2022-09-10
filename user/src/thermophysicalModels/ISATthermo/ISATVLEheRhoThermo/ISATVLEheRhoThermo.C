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
    const scalarField &hCells = this->he_;
    const scalarField &pCells = this->p_;

    scalarField &TCells = this->T_.primitiveFieldRef();
    scalarField &psiCells = this->psi_.primitiveFieldRef();
    scalarField &muCells = this->mu_.primitiveFieldRef();
    scalarField &alphaCells = this->alpha_.primitiveFieldRef();
    scalarField &vaporfracCells = this->vaporfrac_.primitiveFieldRef();
    scalarField &soundspeedCells = this->soundspeed_.primitiveFieldRef();
    scalarField &rhoCells = this->rho_.primitiveFieldRef();
    scalarField &kappaCells = this->kappa_.primitiveFieldRef();
    scalar p_temp, T_temp;
    if (noVLE_)
    {
        forAll(TCells, celli)
        {
            const typename MixtureType::thermoType &mixture_ =
                this->cellMixture(celli);
            auto sol(mixture_.TPN(pCells[celli], TCells[celli]));
            rhoCells[celli] = mixture_.rho_noVLE(pCells[celli], TCells[celli], sol());
            soundspeedCells[celli] = mixture_.c_noVLE(pCells[celli], TCells[celli], sol());
            psiCells[celli] = rhoCells[celli] / pCells[celli];
            vaporfracCells[celli] = 1;
        }
    }
    else
    {
        forAll(TCells, celli)
        {
            const typename MixtureType::thermoType &mixture_ =
                this->cellMixture(celli);

            std::tie(rhoCells[celli], vaporfracCells[celli], soundspeedCells[celli]) = mixture_.rhovfc_ISAT(pCells[celli], TCells[celli]);
            psiCells[celli] = rhoCells[celli] / pCells[celli];
        }
    }
    if (!inviscid_)
    {
        forAll(TCells, celli)
        {
            const typename MixtureType::thermoType &mixture_ =
                this->cellMixture(celli);
            muCells[celli] = mixture_.mu_rho(pCells[celli], TCells[celli], rhoCells[celli]);
            kappaCells[celli] = mixture_.kappa_rho(pCells[celli], TCells[celli], rhoCells[celli]);
        }

        forAll(heList_, i)
        {
            forAll(TCells, celli)
            {
                const typename MixtureType::thermoType &mixture_ =
                    this->cellMixture(celli);
                heList_[i].primitiveFieldRef()[celli] = this->speciesData()[i].Hs(pCells[celli], TCells[celli]);

                Dimix_[i].primitiveFieldRef()[celli] = mixture_.Dimix(pCells[celli], TCells[celli], i);
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
        // fvPatchScalarField& pentropy = entropyBf[patchi];
        if (noVLE_)
        {
            if (pT.fixesValue())
            {
                forAll(pT, facei)
                {
                    const typename MixtureType::thermoType &mixture_ =
                        this->patchFaceMixture(patchi, facei);

                    auto sol(mixture_.TPN(pp[facei], pT[facei]));
                    prho[facei] = mixture_.rho_noVLE(pp[facei], pT[facei], sol());
                    psoundspeed[facei] = mixture_.c_noVLE(pp[facei], pT[facei], sol());
                    ppsi[facei] = prho[facei] / pp[facei];
                    pvaporfrac[facei] = 1;
                }
            }
            else
            {
                forAll(pT, facei)
                {
                    const typename MixtureType::thermoType &mixture_ =
                        this->patchFaceMixture(patchi, facei);
                    auto sol(mixture_.TPN(pp[facei], pT[facei]));
                    prho[facei] = mixture_.rho_noVLE(pp[facei], pT[facei], sol());
                    psoundspeed[facei] = mixture_.c_noVLE(pp[facei], pT[facei], sol());
                    // std::tie(rhoCells[celli], vaporfracCells[celli], soundspeedCells[celli]) = mixture_.rhovfc_ISAT(pCells[celli], TCells[celli]);
                    ppsi[facei] = prho[facei] / pp[facei];
                    pvaporfrac[facei] = 1;
                }
            }
        }
        else
        {
            if (pT.fixesValue())
            {
                forAll(pT, facei)
                {
                    const typename MixtureType::thermoType &mixture_ =
                        this->patchFaceMixture(patchi, facei);
                    std::tie(prho[facei], pvaporfrac[facei], psoundspeed[facei]) = mixture_.rhovfc_ISAT(pp[facei], pT[facei]);
                    ppsi[facei] = prho[facei] / pp[facei];
                }
            }
            else
            {
                forAll(pT, facei)
                {
                    const typename MixtureType::thermoType &mixture_ =
                        this->patchFaceMixture(patchi, facei);
                    std::tie(prho[facei], pvaporfrac[facei], psoundspeed[facei]) = mixture_.rhovfc_ISAT(pp[facei], pT[facei]);
                    ppsi[facei] = prho[facei] / pp[facei];
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
                    pmu[facei] = mixture_.mu_rho(pp[facei], pT[facei], prho[facei]);
                    pkappa[facei] = mixture_.kappa_rho(pp[facei], pT[facei], prho[facei]);
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
                    pmu[facei] = mixture_.mu_rho(pp[facei], pT[facei], prho[facei]);
                    pkappa[facei] = mixture_.kappa_rho(pp[facei], pT[facei], prho[facei]);
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

                        pDimix[facei] = mixture_.Dimix(pp[facei], pT[facei], i);

                        pheList[facei] = this->speciesData()[i].Hs(pp[facei], pT[facei]);
                    }
                }
                else
                {
                    forAll(pT, facei)
                    {
                        const typename MixtureType::thermoType &mixture_ =
                            this->patchFaceMixture(patchi, facei);
                        pDimix[facei] = mixture_.Dimix(pp[facei], pT[facei], i);
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
    scalar tempT, tempP, maxdT = 0, maxdP = 0, tempmu, temppsi, tempHe;
    static scalar maxdmu = 0;
    static int timeflag = 0;
    this->newTimeStep();
    if (noVLE_)
    {
        forAll(TCells, celli)
        {
            const typename MixtureType::thermoType &mixture_ =
                this->cellMixture(celli);
            TCells[celli] = mixture_.THP(hCells[celli] + pCells[celli] / rhoCells[celli], pCells[celli], TCells[celli]);
            vaporfracCells[celli] = 1;
            soundspeedCells[celli] = mixture_.c_noVLE(pCells[celli], TCells[celli]); //, sol());
            // std::tie(TCells[celli], temppsi, vaporfracCells[celli], soundspeedCells[celli]) = mixture_.Tpsivfc_XHP(hCells[celli] + pCells[celli] / rhoCells[celli], pCells[celli], TCells[celli]);
            psiCells[celli] = rhoCells[celli] / pCells[celli];
        }
    }
    else
    {
        forAll(TCells, celli)
        {
            const typename MixtureType::thermoType &mixture_ = this->cellMixture(celli);

            // std::tie(TCells[celli], psiCells[celli], vaporfracCells[celli], soundspeedCells[celli]) = mixture_.Tpsivfc_XHP(hCells[celli] + pCells[celli] / rhoCells[celli], pCells[celli], TCells[celli]);
            std::tie(TCells[celli], tempHe, vaporfracCells[celli], soundspeedCells[celli]) = mixture_.THvfc_XrhoP(rhoCells[celli], pCells[celli], TCells[celli]);
            hCells[celli] = mixture_.HE(pCells[celli], TCells[celli]);
            // hCells[celli] -= pCells[celli] / rhoCells[celli];
            // autoPtr<typename MixtureType::thermoType::solution> sol;
            // std::tie(TCells[celli], temppsi, vaporfracCells[celli], soundspeedCells[celli], sol) = mixture_.Tpsivfcsol_XHP(hCells[celli] + pCells[celli] / rhoCells[celli], pCells[celli], TCells[celli]);
            // rhoCells[celli] = psiCells[celli] * pCells[celli];
            psiCells[celli] = rhoCells[celli] / pCells[celli];
        }
    }
    if (!inviscid_)
    {
        forAll(TCells, celli)
        {
            const typename MixtureType::thermoType &mixture_ =
                this->cellMixture(celli);
            muCells[celli] = mixture_.mu_rho(pCells[celli], TCells[celli], rhoCells[celli]);
            kappaCells[celli] = mixture_.kappa_rho(pCells[celli], TCells[celli], rhoCells[celli]);
            // alphaCells[celli] = mixture_.alphah(pCells[celli], TCells[celli]);//kappaCells[celli] / mixture_.Cp(pCells[celli], TCells[celli]);
        }

        forAll(heList_, i)
        {
            forAll(TCells, celli)
            {
                const typename MixtureType::thermoType &mixture_ =
                    this->cellMixture(celli);
                heList_[i].primitiveFieldRef()[celli] = this->speciesData()[i].Hs(pCells[celli], TCells[celli]);

                Dimix_[i].primitiveFieldRef()[celli] = mixture_.Dimix(pCells[celli], TCells[celli], i);
            }
        }
    }

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
        if (noVLE_) // todo add transport for no VLE
        {
            if (pT.fixesValue())
            {
                forAll(pT, facei)
                {
                    const typename MixtureType::thermoType &mixture_ =
                        this->patchFaceMixture(patchi, facei);

                    phe[facei] = mixture_.HE(pp[facei], pT[facei]);
                    // ppsi[facei] = prho[facei] / pp[facei];

                    auto sol(mixture_.TPN(pp[facei], pT[facei]));
                    prho[facei] = mixture_.rho_noVLE(pp[facei], pT[facei], sol());
                    psoundspeed[facei] = mixture_.c_noVLE(pp[facei], pT[facei], sol());
                    // hCells[celli] =mixture_.HE(pCells[celli], TCells[celli]);
                    // std::tie(rhoCells[celli], vaporfracCells[celli], soundspeedCells[celli]) = mixture_.rhovfc_ISAT(pCells[celli], TCells[celli]);
                    ppsi[facei] = prho[facei] / pp[facei];
                    pvaporfrac[facei] = 1;
                }
            }
            else
            {
                forAll(pT, facei)
                {
                    const typename MixtureType::thermoType &mixture_ =
                        this->patchFaceMixture(patchi, facei);
                    pT[facei] = mixture_.THP(phe[facei] + pp[facei] / prho[facei], pp[facei], pT[facei]);
                    pvaporfrac[facei] = 1;
                    psoundspeed[facei] = mixture_.c_noVLE(pp[facei], pT[facei]);
                    // std::tie(pT[facei], temppsi, pvaporfrac[facei], psoundspeed[facei]) = mixture_.Tpsivfc_XHP(phe[facei] - pp[facei] / prho[facei], pp[facei], pT[facei]);

                    phe[facei] = mixture_.HE(pp[facei], pT[facei]);
                    ppsi[facei] = prho[facei] / pp[facei];
                }
            }
        }
        else
        {
            if (pT.fixesValue())
            {
                forAll(pT, facei)
                {
                    const typename MixtureType::thermoType &mixture_ =
                        this->patchFaceMixture(patchi, facei);
                    std::tie(prho[facei], pvaporfrac[facei], psoundspeed[facei]) = mixture_.rhovfc_ISAT(pp[facei], pT[facei]);

                    phe[facei] = mixture_.HE(pp[facei], pT[facei]);
                    ppsi[facei] = prho[facei] / pp[facei];
                }
            }
            else
            {
                forAll(pT, facei)
                {
                    const typename MixtureType::thermoType &mixture_ =
                        this->patchFaceMixture(patchi, facei);

                    // std::tie(pT[facei], ppsi[facei], pvaporfrac[facei]) = mixture_.Tpsivf_HP(phe[facei], pp[facei], pT[facei]);
                    // std::tie(pT[facei], pp[facei], pvaporfrac[facei], psoundspeed[facei]) = mixture_.TPvf_Erho(phe[facei], prho[facei], pT[facei], pp[facei]);

                    // std::tie(pT[facei], ppsi[facei], pvaporfrac[facei], psoundspeed[facei]) = mixture_.Tpsivfc_XHP(phe[facei] + pp[facei] / prho[facei], pp[facei], pT[facei]);
                    std::tie(pT[facei], tempHe, pvaporfrac[facei], psoundspeed[facei]) = mixture_.THvfc_XrhoP(prho[facei], pp[facei], pT[facei]);

                    ppsi[facei] = prho[facei] / pp[facei];
                    // autoPtr<typename MixtureType::thermoType::solution> sol;
                    // std::tie(pT[facei], temppsi, pvaporfrac[facei], psoundspeed[facei],sol) = mixture_.Tpsivfcsol_XHP(phe[facei] + pp[facei] / prho[facei], pp[facei], pT[facei]);
                    // ppsi[facei] = prho[facei] / pp[facei];
                }
            }
        }
        if (!inviscid_)
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType &mixture_ =
                    this->patchFaceMixture(patchi, facei);
                pmu[facei] = mixture_.mu_rho(pp[facei], pT[facei], prho[facei]);
                pkappa[facei] = mixture_.kappa_rho(pp[facei], pT[facei], prho[facei]);
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

                        pDimix[facei] = mixture_.Dimix(pp[facei], pT[facei], i);
                        pheList[facei] = this->speciesData()[i].Hs(pp[facei], pT[facei]);
                    }
                }
                else
                {
                    forAll(pT, facei)
                    {
                        const typename MixtureType::thermoType &mixture_ =
                            this->patchFaceMixture(patchi, facei);
                        pDimix[facei] = mixture_.Dimix(pp[facei], pT[facei], i);

                        pheList[facei] = this->speciesData()[i].Hs(pp[facei], pT[facei]);
                    }
                }
            }
        }
    }
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
              IOobject::NO_WRITE),
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
      Dimix_(MixtureType::Y().size()),
      heList_(MixtureType::Y().size()),
      WList_(MixtureType::Y().size()),
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
    noVLE_ = thermoDict.lookupOrDefault<bool>("noVLE", false);
    MixtureType::thermoType::noVLE = noVLE_;
    // FatalErrorInFunction
    //     << "inviscid_:" <<inviscid_
    //     << exit(FatalError);

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

// ************************************************************************* //
