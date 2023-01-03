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
    scalar he_temp;
    scalarList Y_temp(MixtureType::Y().size());

    forAll(TCells, celli)
    {
        const typename MixtureType::thermoType &mixture_ =
            this->cellMixture(celli);

        //std::tie(he_temp, rhoCells[celli], vaporfracCells[celli], soundspeedCells[celli]) = mixture_.Erhovfc_XPT(pCells[celli], TCells[celli]);
        std::tie(hCells[celli], rhoCells[celli], vaporfracCells[celli], soundspeedCells[celli]) = mixture_.Erhovfc_XPT(pCells[celli], TCells[celli]);
        psiCells[celli] = rhoCells[celli] / pCells[celli];
    }

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

    //volScalarField::Boundary &heBf =
    //    this->he().boundaryFieldRef();

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

    forAll(this->T_.boundaryField(), patchi)
    {
        fvPatchScalarField &pp = pBf[patchi];
        fvPatchScalarField &pT = TBf[patchi];
        fvPatchScalarField &ppsi = psiBf[patchi];
        //fvPatchScalarField &phe = heBf[patchi];
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

                std::tie(he_temp, prho[facei], pvaporfrac[facei], psoundspeed[facei]) = mixture_.Erhovfc_XPT(pp[facei], pT[facei]);

                ppsi[facei] = prho[facei] / pp[facei];
            }
        }
        else
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType &mixture_ =
                    this->patchFaceMixture(patchi, facei);
                std::tie(he_temp, prho[facei], pvaporfrac[facei], psoundspeed[facei]) = mixture_.Erhovfc_XPT(pp[facei], pT[facei]);

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
    //scalarField &alphaCells = this->alpha_.primitiveFieldRef();
    scalarField &vaporfracCells = this->vaporfrac_.primitiveFieldRef();
    scalarField &soundspeedCells = this->soundspeed_.primitiveFieldRef();
    scalarField &rhoCells = this->rho_.primitiveFieldRef();
    scalarField &kappaCells = this->kappa_.primitiveFieldRef();
    //scalarField &rho_G_Cells = this->rho_G_.primitiveFieldRef();
    //scalarField &ZCells = this->Z_.primitiveFieldRef();
    //scalar tempT, tempP, maxdT = 0, maxdP = 0, tempmu, temppsi, tempHe;
    //static scalar maxdmu = 0;
    //static int timeflag = 0;
    scalarList Y_temp(MixtureType::Y().size());
    scalar VLEtime = 0;

    if (!boundary_flag)
    {
        this->newTimeStep();
        clockTime_.timeIncrement();
        if (scheme_ == "doubleFlux")
        {
            do
            {
                forAll(TCells, celli)
                {
                    const typename MixtureType::thermoType &mixture_ = this->cellMixture(celli);

                    std::tie(TCells[celli], hCells[celli], vaporfracCells[celli], soundspeedCells[celli]) = mixture_.THvfc_XrhoP(rhoCells[celli], pCells[celli], TCells[celli]);

                    hCells[celli] -= pCells[celli] / rhoCells[celli];

                    psiCells[celli] = rhoCells[celli] / pCells[celli];
                }
            } while (this->newLoop());
        }
        else if (scheme_ == "conservativeFlux")
        {
            do
            {
                forAll(TCells, celli)
                {
                    const typename MixtureType::thermoType &mixture_ = this->cellMixture(celli);

                    std::tie(TCells[celli], pCells[celli], vaporfracCells[celli], soundspeedCells[celli]) = mixture_.TPvfc_XErho(hCells[celli], rhoCells[celli], TCells[celli], pCells[celli]);

                    psiCells[celli] = rhoCells[celli] / pCells[celli];
                }
            } while (this->newLoop());
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
                        std::tie(TCells[celli], hCells[celli], vaporfracCells[celli], soundspeedCells[celli]) = mixture_.THvfc_XrhoP(rhoCells[celli], pCells[celli], TCells[celli]);
                        hCells[celli] -= pCells[celli] / rhoCells[celli];
                    }
                    else if ((*FCcell)[celli] == 1)
                    {
                        std::tie(TCells[celli], pCells[celli], vaporfracCells[celli], soundspeedCells[celli]) = mixture_.TPvfc_XErho(hCells[celli], rhoCells[celli], TCells[celli], pCells[celli]);
                    }

                    psiCells[celli] = rhoCells[celli] / pCells[celli];
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
    //DF_ = thermoDict.lookupOrDefault<bool>("doubleFlux", true);
    ISATlog_ = thermoDict.lookupOrDefault<bool>("ISATlog", false);
    scheme_ = word(thermoDict.lookup("scheme"));
    //noVLE_ = thermoDict.lookupOrDefault<bool>("noVLE", false);
    //MixtureType::thermoType::noVLE = noVLE_;
    // FatalErrorInFunction
    //     << "inviscid_:" <<inviscid_
    //     << exit(FatalError);
    FCcell = &mesh.objectRegistry::lookupObjectRef<volScalarField>("FCcell");
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
