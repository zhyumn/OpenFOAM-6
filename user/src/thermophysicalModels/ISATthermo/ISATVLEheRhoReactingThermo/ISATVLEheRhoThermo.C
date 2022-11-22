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
    scalarField &rho_G_Cells = this->rho_G_.primitiveFieldRef();
    //scalarField &ZCells = this->Z_.primitiveFieldRef();
    scalar p_temp, T_temp, he_temp;
    scalarList Y_temp(MixtureType::Y().size());

    forAll(TCells, celli)
    {
        const typename MixtureType::thermoType &mixture_ =
            this->cellMixture(celli);

        //std::tie(rhoCells[celli], vaporfracCells[celli], soundspeedCells[celli]) = mixture_.rhovfc_ISAT(pCells[celli], TCells[celli]);
        std::tie(he_temp, rhoCells[celli], vaporfracCells[celli], soundspeedCells[celli], rho_G_Cells[celli]) = mixture_.Erhovfc_G_rhoY_ISAT(pCells[celli], TCells[celli], Y_temp);
        //ZCells[celli] = mixture_.Z(pCells[celli], TCells[celli]);
        for (int i = 0; i < Y_temp.size(); i++)
        {
            Y_G_List_[i][celli] = Y_temp[i];
        }
        psiCells[celli] = rhoCells[celli] / pCells[celli];
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
                heList_[i].primitiveFieldRef()[celli] = this->speciesData()[i].Ha(pCells[celli], TCells[celli]);

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

    volScalarField::Boundary &rho_G_Bf =
        this->rho_G_.boundaryFieldRef();

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

                        pheList[facei] = this->speciesData()[i].Ha(pp[facei], pT[facei]);
                    }
                }
                else
                {
                    forAll(pT, facei)
                    {
                        const typename MixtureType::thermoType &mixture_ =
                            this->patchFaceMixture(patchi, facei);
                        pDimix[facei] = mixture_.Dimix(pp[facei], pT[facei], i);
                        pheList[facei] = this->speciesData()[i].Ha(pp[facei], pT[facei]); // TODO check Hs
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

    if (!boundary_flag)
    {
        this->newTimeStep();

        forAll(TCells, celli)
        {
            const typename MixtureType::thermoType &mixture_ = this->cellMixture(celli);

            // std::tie(TCells[celli], psiCells[celli], vaporfracCells[celli], soundspeedCells[celli]) = mixture_.Tpsivfc_XHP(hCells[celli] + pCells[celli] / rhoCells[celli], pCells[celli], TCells[celli]);
            std::tie(TCells[celli], hCells[celli], vaporfracCells[celli], soundspeedCells[celli], rho_G_[celli]) = mixture_.THvfc_G_rhoY_XrhoP(rhoCells[celli], pCells[celli], TCells[celli], Y_temp);
            //hCells[celli] = mixture_.HE(pCells[celli], TCells[celli]);
            hCells[celli] -= pCells[celli] / rhoCells[celli];
            //ZCells[celli] = mixture_.Z(pCells[celli], TCells[celli]);
            for (int i = 0; i < Y_temp.size(); i++)
            {
                Y_G_List_[i][celli] = Y_temp[i];
            }
            // autoPtr<typename MixtureType::thermoType::solution> sol;
            // std::tie(TCells[celli], temppsi, vaporfracCells[celli], soundspeedCells[celli], sol) = mixture_.Tpsivfcsol_XHP(hCells[celli] + pCells[celli] / rhoCells[celli], pCells[celli], TCells[celli]);
            // rhoCells[celli] = psiCells[celli] * pCells[celli];
            psiCells[celli] = rhoCells[celli] / pCells[celli];
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
                    heList_[i].primitiveFieldRef()[celli] = this->speciesData()[i].Ha(pCells[celli], TCells[celli]);

                    Dimix_[i].primitiveFieldRef()[celli] = mixture_.Dimix(pCells[celli], TCells[celli], i);
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
                            pheList[facei] = this->speciesData()[i].Ha(pp[facei], pT[facei]);
                        }
                    }
                    else
                    {
                        forAll(pT, facei)
                        {
                            const typename MixtureType::thermoType &mixture_ =
                                this->patchFaceMixture(patchi, facei);
                            pDimix[facei] = mixture_.Dimix(pp[facei], pT[facei], i);

                            pheList[facei] = this->speciesData()[i].Ha(pp[facei], pT[facei]);
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
      rho_G_(
          IOobject(
              "thermo:rho_G",
              mesh.time().timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          mesh,
          dimDensity),
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
    //noVLE_ = thermoDict.lookupOrDefault<bool>("noVLE", false);
    //MixtureType::thermoType::noVLE = noVLE_;
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

template <class BasicPsiThermo, class MixtureType>
void Foam::ISATVLEheRhoThermo<BasicPsiThermo, MixtureType>::customChemistry(scalar temperature, scalar pressure, scalar density)
{

    // Mw kg/mol
    const typename MixtureType::thermoType &mixture_ =
        this->cellMixture(0);

    W_r[0] = 0.170338;//(mixture_[0].W())*1E-3; // c12h26
    W_r[1] = 0.032;//(mixture_[1].W())*1E-3; // O2
    W_r[2] = 0.04401;//(mixture_[2].W())*1E-3; // co2
    W_r[3] =  0.0280101;//(mixture_[3].W())*1E-3; // co
    W_r[4] = 0.018;//(mixture_[4].W())*1E-3; // h2o
    W_r[5] = 0.0280134;//(mixture_[5].W())*1E-3; // N2

    Hc[0] = -290.9; // c12h26   kJ/mol
    Hc[1] = 0.0; // O2
    Hc[2] = -393.51; // co2
    Hc[3] = -110.53; // co
    Hc[4] = -241.826; // h2o
    Hc[5] = 0.0;// N2

    for(int i = 0;i<6;i++)
    {
        Hc[i] = Hc[i] * 1E3 / W_r[i];  //J/kg
        //Info<<Hc[i]<<endl;
    }
    //Info<<temperature<<"\t"<<pressure<<"\t"<<"\t"<<density;

    setThermo(temperature, pressure, density);
}

template <class BasicPsiThermo, class MixtureType>
void Foam::ISATVLEheRhoThermo<BasicPsiThermo, MixtureType>::setThermo(scalar temperature, scalar pressure, scalar density)
{
    
    T_r = temperature;
    P_r = pressure;
    rho_r = density;
}

template <class BasicPsiThermo, class MixtureType>
void Foam::ISATVLEheRhoThermo<BasicPsiThermo, MixtureType>::setY(double Y, int speciesLabel)
{
    Ygas_r[speciesLabel] = Y;

    //Info<<"y \t"<<speciesLabel<<"\t"<<Ygas_r[speciesLabel]<<endl;
}

template <class BasicPsiThermo, class MixtureType>
void Foam::ISATVLEheRhoThermo<BasicPsiThermo, MixtureType>::setConc()
{

    for (int i = 0; i < 6; i++)
    {
        c_r[i] = rho_r * Ygas_r[i] *1E-6 / W_r[i];
    }
    // for(int i=0;i<6;i++){
    //     Info<<"Species \t"<<i<<"\t conc = "<<c_r[i]<<"\t"<<Ygas_r[i]<<endl;
    // }
}

template <class BasicPsiThermo, class MixtureType>
void Foam::ISATVLEheRhoThermo<BasicPsiThermo, MixtureType>::setRate()
{
    double A1;
    double R1, R2f, R2b, R2;
    double t0, t1, t2, t3, t4, t5, t6;
    double phi_st;
    double x_f, x_o, Mavg, Temp, err ;
    double T_ini,Tf,To;
    double Z,yf,yo,yff,yoo,s,S;

    // reactions 2 step
    // c12h26 + 12.5O2 -> 12CO + 13H2O   R1 = A1 exp(-31944/(RT)) [c12h26]^0.25 [O2]^1.25
    // CO + 0.5O2 <=> CO2 R2 = 3.98x10^14 exp(-40000/RT)[CO]^1 [H2O]^0.5 [O2]^0.25 - 5x10^8 exp(-40000/RT)[CO2]

    setConc();

    // coefficients for reaction from Hickel, 2022 LES reacting and non reacting transcritical fuel sprays ...
    t0 = 27.38;
    t1 = -2.13;
    t2 = -2.05;
    t3 = 1.89;
    t4 = -0.01;
    t5 = 2.87E-4;
    t6 = 8.43;

    Mavg = 0.0;

    for(int i=0;i<6;i++){
        Mavg += Ygas_r[i]/W_r[i];
    }

    Mavg = 1.0/Mavg;
    Tf = 363.0;
    To = 900.0;

    //mixture fraction estimation
    s = 32.0*(12.0 + 26.0/4.0)/(12.0*12.0+26.0);
    yff = 1.0;
    yoo = 0.15;
    yf = Ygas_r[0]/yff;
    yo = Ygas_r[1]/yoo;
    S = s*yff/(yoo+1E-30);
    Z = (S*yf - yo + 1)/(S+1);

    T_ini = Z*Tf + (1-Z)*To;

    //if (err>5) Info<<"Gas phase temperature "<<T_r<<"\t"<<Temp<<endl; 

    x_f = Ygas_r[0]*Mavg/W_r[0];
    x_o = Ygas_r[1]*Mavg/W_r[1];

    phi_st = (x_f/(x_o+1E-30))/(1.0/18.5);
    
    // 0 = c12h26. 1 = o2, 2 = co2, 3 = co, 4 = h2o, 5 = inert
    A1 = exp(t0 + t1 * exp(t2 * phi_st) + t3 * tanh((t4 + t5 * phi_st) * T_ini + t6));
    R1 = A1 * exp(-31944.0 / (R_r * T_r)) * pow(c_r[0], 0.25) * pow(c_r[1], 1.25);
    R2f = 3.98 * 1E14 * exp(-40000.0 / (R_r * T_r)) * pow(c_r[3], 1.0) * pow(c_r[4], 0.5) * pow(c_r[1], 0.25);
    R2b = 5.0 * 1E8 * exp(-40000.0 / (R_r * T_r)) * pow(c_r[2], 1.0);
    
    if(Ygas_r[0] < 0.000001 || Ygas_r[1] < 0.000001){   
        R1 = 0.0; 
    }

    if(Ygas_r[1] < 0.000001 || Ygas_r[3] < 0.000001 || Ygas_r[4] < 0.000001){   
        R2f = 0.0; 
    }

    if(Ygas_r[2] < 0.000001){   
        R2b = 0.0; 
    } 

    R2 = R2f - R2b;   
    
    reactRate_r[0] = -R1;
    reactRate_r[1] = -12.5 * R1 - 0.5 * R2;
    reactRate_r[2] = R2;
    reactRate_r[3] = 12.0 * R1 - R2;
    reactRate_r[4] = 13.0 * R1;
    reactRate_r[5] = 0.0; // inert

    for(int i=0;i<6;i++){
        reactRate_r[i] = reactRate_r[i]*W_r[i]*1E6;
    }
}

template <class BasicPsiThermo, class MixtureType>
double Foam::ISATVLEheRhoThermo<BasicPsiThermo, MixtureType>::RR(int speciesLabel)
{
    //Info<<"react rate \t"<<speciesLabel<<" = "<<reactRate_r[speciesLabel]<<"\t"<<c_r[speciesLabel]<<endl;
    return reactRate_r[speciesLabel];
}

template <class BasicPsiThermo, class MixtureType>
double Foam::ISATVLEheRhoThermo<BasicPsiThermo, MixtureType>::Qdot()
{
    double q = 0.0;
    const typename MixtureType::thermoType &mixture_ =
        this->cellMixture(0);

    for (int i = 0; i < 6; i++)
    {
        //q = q - mixture_[0].Hc(i) * reactRate_r[i];
        q -= Hc[i]*reactRate_r[i];
    }
    // for(int i=0;i<6;i++){
    //     Info<<"Species \t"<<i<<" temperature \t"<<T_r<<"\t formation enthalpy = "<<mixture_[i].Hc()<<endl;
    // }
    //Info<<"q \t"<<q<<endl;

    return q;
}

template <class BasicPsiThermo, class MixtureType>
double Foam::ISATVLEheRhoThermo<BasicPsiThermo, MixtureType>::gasPhaseRho(scalar pressure, scalar temperature)
{
    const typename MixtureType::thermoType &mixture_ = this->cellMixture(0);
    return mixture_.rho(pressure, temperature);
}

template <class BasicPsiThermo, class MixtureType>
double Foam::ISATVLEheRhoThermo<BasicPsiThermo, MixtureType>::gasPhaseY(label celli, label speciesLabel)
{
    const typename MixtureType::thermoType &mixture_ = this->cellMixture(celli);
    return 0;
}


// template <class BasicPsiThermo, class MixtureType>
// void Foam::ISATVLEheRhoThermo<BasicPsiThermo, MixtureType>:: customCHemistry(volScalarField temperature, volScalarField pressure, volScalarField density, PtrList<volScalarField> Y)
// {
//     Temperature = temperature;
//     Pressure = pressure;
//     rho = density;
//     forAll(Y,i)
//     {
//         Yrxn[i] = Y[i];
//     }

//     W_r[0] = 0.170338;//(mixture_[0].W())*1E-3; // c12h26
//     W_r[1] = 0.032;//(mixture_[1].W())*1E-3; // O2
//     W_r[2] = 0.04401;//(mixture_[2].W())*1E-3; // co2
//     W_r[3] =  0.0280101;//(mixture_[3].W())*1E-3; // co
//     W_r[4] = 0.018;//(mixture_[4].W())*1E-3; // h2o
//     W_r[5] = 0.0280134;//(mixture_[5].W())*1E-3; // N2

//     Hc[0] = -290.9; // c12h26   kJ/mol
//     Hc[1] = 0; // O2
//     Hc[2] = -393.51; // co2
//     Hc[3] = -110.53; // co
//     Hc[4] = -241.826; // h2o
//     Hc[5] = 0;// N2

//     for(int i = 0;i<6;i++)
//     {
//         Hc[i] = Hc[i] * 1E3 / W_r[i];  //J/kg
//     }


//     calculateReactions();
// }


// template <class BasicPsiThermo, class MixtureType>
// void Foam::ISATVLEheRhoThermo<BasicPsiThermo, MixtureType>:: calculateReactions()
// {
//     forAll(Temperature,celli)
//     {
//         forAll(Yrxn,i)
//         {
//             conc[i][celli] =  rho[celli] * Yrxn[i][celli] * 1E-6 / W_r[i];
//         }

//         scalar R1 = R1(celli);
//         scalar R2 = R2(celli);

//         rate[0][celli] = -R1;
//         rate[1][celli] = -12.5 * R1 - 0.5 * R2;
//         rate[2][celli] = R2;
//         rate[3][celli] = 12.0 * R1 - R2;
//         rate[4][celli] = 13.0 * R1;
//         rate[5][celli] = 0.0; // inert

//         Qdot[celli] = 0;
//         forAll(Yrxn,i)
//         {
//             Qdot[celli] = Qdot[celli] - rate[i][celli]*Hc[i];
//         }
//     }
// }

// template <class BasicPsiThermo, class MixtureType>
// volScalarField Foam::ISATVLEheRhoThermo<BasicPsiThermo, MixtureType>:: RR(label i)
// {
//     return rate[i];
// }

// template <class BasicPsiThermo, class MixtureType>
// volScalarField Foam::ISATVLEheRhoThermo<BasicPsiThermo, MixtureType>:: Qdot()
// {
//      return Qdot;
// }