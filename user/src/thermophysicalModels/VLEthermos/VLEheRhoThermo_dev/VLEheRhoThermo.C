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

#include "VLEheRhoThermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
void Foam::VLEheRhoThermo<BasicPsiThermo, MixtureType>::calculate()
{

    //double tempT;
    const scalarField& hCells = this->he_;
    scalarField& pCells = this->p_;

    scalarField& TCells = this->T_.primitiveFieldRef();
    scalarField& psiCells = this->psi_.primitiveFieldRef();
    scalarField& muCells = this->mu_.primitiveFieldRef();
    scalarField& alphaCells = this->alpha_.primitiveFieldRef();
    scalarField& vaporfracCells = this->vaporfrac_.primitiveFieldRef();
    scalarField& soundspeedCells = this->soundspeed_.primitiveFieldRef();
    scalarField& entropyCells = this->entropy_.primitiveFieldRef();
    const scalarField& rhoCells = this->rho_.primitiveFieldRef();
    scalarField& kappaCells = this->kappa_.primitiveFieldRef();
    //scalar p_temp, T_temp;

    forAll(TCells, celli)
    {
        const typename MixtureType::thermoType& mixture_ =
            this->cellMixture(celli);
        /*
        TCells[celli] = mixture_.THE
        (
            hCells[celli],
            pCells[celli],
            TCells[celli]
        );*/
        std::tie(pCells[celli], TCells[celli]) = mixture_.TrhoE
        (
            hCells[celli],
            rhoCells[celli],
            pCells[celli],
            TCells[celli]
        );

        //psiCells[celli] = mixture_.psi(pCells[celli], TCells[celli]);
        psiCells[celli] = rhoCells[celli] / pCells[celli];
        muCells[celli] = mixture_.mu_dev(pCells[celli], TCells[celli]);
        kappaCells[celli] = mixture_.kappa_dev(pCells[celli], TCells[celli]);
        alphaCells[celli] = kappaCells[celli] / mixture_.Cp(pCells[celli], TCells[celli]);
        //alphaCells[celli] = mixture_.alphah(pCells[celli], TCells[celli]);
        vaporfracCells[celli] = mixture_.vaporfra(pCells[celli], TCells[celli]);
        soundspeedCells[celli] = mixture_.c(pCells[celli], TCells[celli]);
        entropyCells[celli] = mixture_.S(pCells[celli], TCells[celli]);
    }

    forAll(heList_, i)
    {
        forAll(TCells, celli)
        {
            const typename MixtureType::thermoType& mixture_ =
                this->cellMixture(celli);
            heList_[i].primitiveFieldRef()[celli]
                = this->speciesData()[i].HE(pCells[celli], TCells[celli]);

            Dimix_[i].primitiveFieldRef()[celli]
                = mixture_.Dimix(pCells[celli], TCells[celli], i);

        }
    }


    volScalarField::Boundary& pBf =
        this->p_.boundaryFieldRef();

    volScalarField::Boundary& TBf =
        this->T_.boundaryFieldRef();

    volScalarField::Boundary& psiBf =
        this->psi_.boundaryFieldRef();

    volScalarField::Boundary& heBf =
        this->he().boundaryFieldRef();

    volScalarField::Boundary& muBf =
        this->mu_.boundaryFieldRef();

    volScalarField::Boundary& alphaBf =
        this->alpha_.boundaryFieldRef();

    volScalarField::Boundary& vaporfracBf =
        this->vaporfrac_.boundaryFieldRef();

    volScalarField::Boundary& soundspeedBf =
        this->soundspeed_.boundaryFieldRef();

    volScalarField::Boundary& entropyBf =
        this->entropy_.boundaryFieldRef();

    volScalarField::Boundary& kappaBf =
        this->kappa_.boundaryFieldRef();
    //scalarField& kappaCells = this->kappa_.primitiveFieldRef();

    forAll(this->T_.boundaryField(), patchi)
    {
        fvPatchScalarField& pp = pBf[patchi];
        fvPatchScalarField& pT = TBf[patchi];
        fvPatchScalarField& ppsi = psiBf[patchi];
        fvPatchScalarField& phe = heBf[patchi];
        fvPatchScalarField& pmu = muBf[patchi];
        fvPatchScalarField& palpha = alphaBf[patchi];
        fvPatchScalarField& pvaporfrac = vaporfracBf[patchi];
        fvPatchScalarField& psoundspeed = soundspeedBf[patchi];
        fvPatchScalarField& pentropy = entropyBf[patchi];
        fvPatchScalarField& pkappa = kappaBf[patchi];

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                phe[facei] = mixture_.HE(pp[facei], pT[facei]);

                ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                pmu[facei] = mixture_.mu_dev(pp[facei], pT[facei]);
                pkappa[facei] = mixture_.kappa_dev(pp[facei], pT[facei]);
                palpha[facei] = pkappa[facei] / mixture_.Cp(pp[facei], pT[facei]);
                //palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
                pvaporfrac[facei] = mixture_.vaporfra(pp[facei], pT[facei]);
                psoundspeed[facei] = mixture_.c(pp[facei], pT[facei]);
                pentropy[facei] = mixture_.S(pp[facei], pT[facei]);
            }
        }
        else
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                pT[facei] = mixture_.THE(phe[facei], pp[facei], pT[facei]);

                ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                pmu[facei] = mixture_.mu_dev(pp[facei], pT[facei]);
                pkappa[facei] = mixture_.kappa_dev(pp[facei], pT[facei]);
                palpha[facei] = pkappa[facei] / mixture_.Cp(pp[facei], pT[facei]);
                //palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
                pvaporfrac[facei] = mixture_.vaporfra(pp[facei], pT[facei]);
                psoundspeed[facei] = mixture_.c(pp[facei], pT[facei]);
                pentropy[facei] = mixture_.S(pp[facei], pT[facei]);
            }
        }
    }
    forAll(heList_, i)
    {
        volScalarField::Boundary& DimixBf =
            this->Dimix_[i].boundaryFieldRef();


        volScalarField::Boundary& heListBf =
            this->heList_[i].boundaryFieldRef();
        forAll(this->T_.boundaryField(), patchi)
        {
            fvPatchScalarField& pDimix = DimixBf[patchi];
            fvPatchScalarField& pheList = heListBf[patchi];
            fvPatchScalarField& pT = TBf[patchi];
            fvPatchScalarField& pp = pBf[patchi];
            if (pT.fixesValue())
            {
                forAll(pT, facei)
                {
                    const typename MixtureType::thermoType& mixture_ =
                        this->patchFaceMixture(patchi, facei);

                    pDimix[facei]
                        = mixture_.Dimix(pp[facei], pT[facei], i);
                    pheList[facei]
                        = this->speciesData()[i].HE(pp[facei], pT[facei]);
                }
            }
            else
            {
                forAll(pT, facei)
                {
                    const typename MixtureType::thermoType& mixture_ =
                        this->patchFaceMixture(patchi, facei);
                    pDimix[facei]
                        = mixture_.Dimix(pp[facei], pT[facei], i);

                    pheList[facei]
                        = this->speciesData()[i].HE(pp[facei], pT[facei]);

                }
            }


        }
    }
}

template<class BasicPsiThermo, class MixtureType>
void Foam::VLEheRhoThermo<BasicPsiThermo, MixtureType>::calculate_init()
{

    //double tempT;
    const scalarField& hCells = this->he_;
    const scalarField& pCells = this->p_;

    scalarField& TCells = this->T_.primitiveFieldRef();
    scalarField& psiCells = this->psi_.primitiveFieldRef();
    scalarField& muCells = this->mu_.primitiveFieldRef();
    scalarField& alphaCells = this->alpha_.primitiveFieldRef();
    scalarField& vaporfracCells = this->vaporfrac_.primitiveFieldRef();
    scalarField& soundspeedCells = this->soundspeed_.primitiveFieldRef();
    scalarField& entropyCells = this->entropy_.primitiveFieldRef();
    scalarField& rhoCells = this->rho_.primitiveFieldRef();
    scalarField& kappaCells = this->kappa_.primitiveFieldRef();
    //scalar p_temp, T_temp;

    forAll(TCells, celli)
    {
        const typename MixtureType::thermoType& mixture_ =
            this->cellMixture(celli);
        //tempT = mixture_.THE(hCells[celli], pCells[celli], TCells[celli]);

        TCells[celli] = mixture_.THE
        (
            hCells[celli],
            pCells[celli],
            TCells[celli]
        );
        for (int i = 180;i <= 600;i += 10)
        {
            Info << i << "," << mixture_.Dimix(1e5, i,0)  << endl;
        }
        FatalErrorInFunction
            << "DEBUG: T:" << TCells[celli] << ", p:" << pCells[celli]
            << " mu=" << mixture_.Dimix(1e5, 250, 0)
            << exit(FatalError);

        psiCells[celli] = mixture_.psi(pCells[celli], TCells[celli]);
        rhoCells[celli] = psiCells[celli] * pCells[celli];
        muCells[celli] = mixture_.mu_dev(pCells[celli], TCells[celli]);
        kappaCells[celli] = mixture_.kappa_dev(pCells[celli], TCells[celli]);
        alphaCells[celli] = kappaCells[celli] / mixture_.Cp(pCells[celli], TCells[celli]);
        //alphaCells[celli] = mixture_.alphah(pCells[celli], TCells[celli]);
        vaporfracCells[celli] = mixture_.vaporfra(pCells[celli], TCells[celli]);
        //muCells[celli] = mixture_.mu(100000, 293);
        //alphaCells[celli] = mixture_.alphah(100000, 293);
        //vaporfracCells[celli] = mixture_.vaporfra(100000, 293);
        //vaporfracCells[celli] = mixture_.vaporfra(100000, 293);
        //vaporfracCells[celli] = mixture_.vaporfra(100000, 293);
        //FatalErrorInFunction
        //    << "DEBUG: T:" <<TCells[celli]<<", p:"<<pCells[celli]
        //    << exit(FatalError);
        //vaporfracCells[celli] = mixture_.vaporfra(pCells[celli], TCells[celli]);
        //FatalErrorInFunction
        //    << "DEBUG: T:" <<TCells[celli]<<", p:"<<pCells[celli]
        //    << exit(FatalError);
        soundspeedCells[celli] = mixture_.c(pCells[celli], TCells[celli]);
        entropyCells[celli] = mixture_.S(pCells[celli], TCells[celli]);


    }
    forAll(heList_, i)
    {
        forAll(TCells, celli)
        {
            const typename MixtureType::thermoType& mixture_ =
                this->cellMixture(celli);
            heList_[i].primitiveFieldRef()[celli]
                = this->speciesData()[i].HE(pCells[celli], TCells[celli]);

            Dimix_[i].primitiveFieldRef()[celli]
                = mixture_.Dimix(pCells[celli], TCells[celli], i);
        }
    }

    volScalarField::Boundary& pBf =
        this->p_.boundaryFieldRef();

    volScalarField::Boundary& TBf =
        this->T_.boundaryFieldRef();

    volScalarField::Boundary& psiBf =
        this->psi_.boundaryFieldRef();

    volScalarField::Boundary& heBf =
        this->he().boundaryFieldRef();

    volScalarField::Boundary& muBf =
        this->mu_.boundaryFieldRef();

    volScalarField::Boundary& alphaBf =
        this->alpha_.boundaryFieldRef();

    volScalarField::Boundary& vaporfracBf =
        this->vaporfrac_.boundaryFieldRef();

    volScalarField::Boundary& soundspeedBf =
        this->soundspeed_.boundaryFieldRef();

    volScalarField::Boundary& entropyBf =
        this->entropy_.boundaryFieldRef();

    volScalarField::Boundary& kappaBf =
        this->kappa_.boundaryFieldRef();

    forAll(this->T_.boundaryField(), patchi)
    {
        fvPatchScalarField& pp = pBf[patchi];
        fvPatchScalarField& pT = TBf[patchi];
        fvPatchScalarField& ppsi = psiBf[patchi];
        fvPatchScalarField& phe = heBf[patchi];
        fvPatchScalarField& pmu = muBf[patchi];
        fvPatchScalarField& palpha = alphaBf[patchi];
        fvPatchScalarField& pvaporfrac = vaporfracBf[patchi];
        fvPatchScalarField& psoundspeed = soundspeedBf[patchi];
        fvPatchScalarField& pentropy = entropyBf[patchi];
        fvPatchScalarField& pkappa = kappaBf[patchi];

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                phe[facei] = mixture_.HE(pp[facei], pT[facei]);

                ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                pmu[facei] = mixture_.mu_dev(pp[facei], pT[facei]);
                pkappa[facei] = mixture_.kappa_dev(pp[facei], pT[facei]);
                palpha[facei] = pkappa[facei] / mixture_.Cp(pp[facei], pT[facei]);
                //palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
                pvaporfrac[facei] = mixture_.vaporfra(pp[facei], pT[facei]);
                psoundspeed[facei] = mixture_.c(pp[facei], pT[facei]);
                pentropy[facei] = mixture_.S(pp[facei], pT[facei]);
            }
        }
        else
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                pT[facei] = mixture_.THE(phe[facei], pp[facei], pT[facei]);

                ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                pmu[facei] = mixture_.mu_dev(pp[facei], pT[facei]);
                pkappa[facei] = mixture_.kappa_dev(pp[facei], pT[facei]);
                palpha[facei] = pkappa[facei] / mixture_.Cp(pp[facei], pT[facei]);
                //palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
                pvaporfrac[facei] = mixture_.vaporfra(pp[facei], pT[facei]);
                psoundspeed[facei] = mixture_.c(pp[facei], pT[facei]);
                pentropy[facei] = mixture_.S(pp[facei], pT[facei]);
            }
        }
    }
    forAll(heList_, i)
    {
        volScalarField::Boundary& DimixBf =
            this->Dimix_[i].boundaryFieldRef();

        volScalarField::Boundary& heListBf =
            this->heList_[i].boundaryFieldRef();

        forAll(this->T_.boundaryField(), patchi)
        {
            fvPatchScalarField& pDimix = DimixBf[patchi];
            fvPatchScalarField& pheList = heListBf[patchi];
            fvPatchScalarField& pT = TBf[patchi];
            fvPatchScalarField& pp = pBf[patchi];
            if (pT.fixesValue())
            {
                forAll(pT, facei)
                {
                    const typename MixtureType::thermoType& mixture_ =
                        this->patchFaceMixture(patchi, facei);

                    pDimix[facei]
                        = mixture_.Dimix(pp[facei], pT[facei], i);
                    pheList[facei]
                        = this->speciesData()[i].HE(pp[facei], pT[facei]);
                }
            }
            else
            {
                forAll(pT, facei)
                {
                    const typename MixtureType::thermoType& mixture_ =
                        this->patchFaceMixture(patchi, facei);
                    pDimix[facei]
                        = mixture_.Dimix(pp[facei], pT[facei], i);
                    pheList[facei]
                        = this->speciesData()[i].HE(pp[facei], pT[facei]);

                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::VLEheRhoThermo<BasicPsiThermo, MixtureType>::VLEheRhoThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
    :
    heThermo<BasicPsiThermo, MixtureType>(mesh, phaseName),
    vaporfrac_
    (
        IOobject
        (
            "vaporfrac",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimless
    ),
    soundspeed_
    (
        IOobject
        (
            "soundspeed",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimVelocity
    ),
    entropy_
    (
        IOobject
        (
            "entropy",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimEnergy / dimTime
    ),
    kappa_
    (
        IOobject
        (
            "thermo:kappa",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimEnergy / (dimTime * dimLength * dimTemperature)
    ),
    Dimix_(MixtureType::Y().size()),
    heList_(MixtureType::Y().size()),
    WList_(MixtureType::Y().size())
    /*,
    rho_
    (
        IOobject
        (
            "rho",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimMass/dimVolume
    )*/
{
    forAll(Dimix_, i)
    {
        Dimix_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    this->phasePropertyName("thermo:Dimix:") + MixtureType::species()[i],
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionSet(0, 2, -1, 0, 0)
            )
        );
        heList_.set
        (
            i,
            new volScalarField
            (

                IOobject
                (
                    "hei:" + MixtureType::species()[i],
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                this->he_.dimensions()
            )
        );
        WList_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "W:" + MixtureType::species()[i],
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh,
                dimMass / dimMoles
            )
        );
    }


    forAll(WList_, i)
    {
        forAll(WList_[i].primitiveFieldRef(), celli)
        {

            WList_[i].primitiveFieldRef()[celli]
                = this->speciesData()[i].W();
        }
    }
    calculate_init();

    // Switch on saving old time
    this->psi_.oldTime();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::VLEheRhoThermo<BasicPsiThermo, MixtureType>::~VLEheRhoThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
void Foam::VLEheRhoThermo<BasicPsiThermo, MixtureType>::correct()
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
