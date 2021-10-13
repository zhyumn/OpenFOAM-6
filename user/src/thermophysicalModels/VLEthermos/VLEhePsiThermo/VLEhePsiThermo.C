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

#include "VLEhePsiThermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
void Foam::VLEhePsiThermo<BasicPsiThermo, MixtureType>::calculate()
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
    scalar p_temp,T_temp;

    forAll(TCells, celli)
    {
        const typename MixtureType::thermoType& mixture_ =
            this->cellMixture(celli);
        //tempT = mixture_.THE(hCells[celli], pCells[celli], TCells[celli]);
        //if (tempT < 452)
        //{
            //Info << hCells[900] << "," << hCells[100] << endl;
        //    FatalErrorInFunction << hCells[900] << "," << hCells[100] << "T=" << tempT << ",h=" << hCells[celli] << ",p=" << pCells[celli] << abort(FatalError);
        //}
        /*
        TCells[celli] = mixture_.THE
        (
            hCells[celli],
            pCells[celli],
            TCells[celli]
        );*/
        std::tie(p_temp,T_temp) =mixture_.TrhoE
        (
            hCells[celli],
            rhoCells[celli],
            pCells[celli],
            TCells[celli]
        );
        pCells[celli]=p_temp;
        TCells[celli]=T_temp;

        psiCells[celli] = mixture_.psi(pCells[celli], TCells[celli]);

        muCells[celli] = mixture_.mu(pCells[celli], TCells[celli]);
        alphaCells[celli] = mixture_.alphah(pCells[celli], TCells[celli]);
        vaporfracCells[celli] = mixture_.vaporfra(pCells[celli], TCells[celli]);
        soundspeedCells[celli] = mixture_.c(pCells[celli], TCells[celli]);
        entropyCells[celli] = mixture_.S(pCells[celli], TCells[celli]);
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

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                phe[facei] = mixture_.HE(pp[facei], pT[facei]);

                ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                pmu[facei] = mixture_.mu(pp[facei], pT[facei]);
                palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
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
                pmu[facei] = mixture_.mu(pp[facei], pT[facei]);
                palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
                pvaporfrac[facei] = mixture_.vaporfra(pp[facei], pT[facei]);
                psoundspeed[facei] = mixture_.c(pp[facei], pT[facei]);
                pentropy[facei] = mixture_.S(pp[facei], pT[facei]);
            }
        }
    }
}

template<class BasicPsiThermo, class MixtureType>
void Foam::VLEhePsiThermo<BasicPsiThermo, MixtureType>::calculate_init()
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
    //scalarField& rhoCells = this->rho_.primitiveFieldRef();
    scalar p_temp,T_temp;

    forAll(TCells, celli)
    {
        const typename MixtureType::thermoType& mixture_ =
            this->cellMixture(celli);
        //tempT = mixture_.THE(hCells[celli], pCells[celli], TCells[celli]);
        //if (tempT < 452)
        //{
            //Info << hCells[900] << "," << hCells[100] << endl;
        //    FatalErrorInFunction << hCells[900] << "," << hCells[100] << "T=" << tempT << ",h=" << hCells[celli] << ",p=" << pCells[celli] << abort(FatalError);
        //}

        TCells[celli] = mixture_.THE
        (
            hCells[celli],
            pCells[celli],
            TCells[celli]
        );
        

        psiCells[celli] = mixture_.psi(pCells[celli], TCells[celli]);

        muCells[celli] = mixture_.mu(pCells[celli], TCells[celli]);
        alphaCells[celli] = mixture_.alphah(pCells[celli], TCells[celli]);
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

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                phe[facei] = mixture_.HE(pp[facei], pT[facei]);

                ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                pmu[facei] = mixture_.mu(pp[facei], pT[facei]);
                palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
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
                pmu[facei] = mixture_.mu(pp[facei], pT[facei]);
                palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
                pvaporfrac[facei] = mixture_.vaporfra(pp[facei], pT[facei]);
                psoundspeed[facei] = mixture_.c(pp[facei], pT[facei]);
                pentropy[facei] = mixture_.S(pp[facei], pT[facei]);
            }
        }
    }
    //this->rho_= this->psi_* this->p_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::VLEhePsiThermo<BasicPsiThermo, MixtureType>::VLEhePsiThermo
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
        dimEnergy/dimTime
    ),
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
    )
{
    calculate_init();

    // Switch on saving old time
    this->psi_.oldTime();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::VLEhePsiThermo<BasicPsiThermo, MixtureType>::~VLEhePsiThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
void Foam::VLEhePsiThermo<BasicPsiThermo, MixtureType>::correct()
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
