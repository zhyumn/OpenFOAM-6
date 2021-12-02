/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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

#include "parallelISAT_chem.H"
namespace Foam
{
    template<class CompType, class ThermoType>
    parallelISAT_chem<CompType, ThermoType>* ISAT_chem<CompType, ThermoType>::pISAT = nullptr;
    template<class CompType, class ThermoType>
    std::ostream& operator<<(std::ostream& out, typename ISAT_chem<CompType, ThermoType>::leafData& A)
    {
        //out << A.phi_[0] << ", " << A.Rv;
        return out;
    }
    template<class CompType, class ThermoType>
    std::ostream& operator<<(std::ostream& out, typename ISAT_chem<CompType, ThermoType>::nodeData& A)
    {
        //out << A.v;
        return out;
    }
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

//const dataType Foam::parallelISAT_chem::staticData();


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //



template<class CompType, class ThermoType>
Foam::parallelISAT_chem<CompType, ThermoType>::parallelISAT_chem(SUPstream::mpi_manager& manager_in, label nmem, SUPstream::mpi_sync& sync_in)
    :
    parallelISAT<ISAT_chem<CompType, ThermoType>, chemistryTabulationMethod<CompType, ThermoType>>(manager_in, nmem, sync_in, dict, chemistry)
{
    Foam::ISAT_chem<CompType, ThermoType>::pISAT = this;
}






// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //




// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
template<class CompType, class ThermoType>
Foam::parallelISAT_chem<CompType, ThermoType>::~parallelISAT_chem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
template<class CompType, class ThermoType>
typename Foam::ISAT_chem<CompType, ThermoType>::outputType Foam::ISAT_chem<CompType, ThermoType>::leafData::func(const inputType& x)
{
    return x;
}
template<class CompType, class ThermoType>
void Foam::ISAT_chem<CompType, ThermoType>::leafData::gradFunc(const inputType& x, typename Foam::ISAT_chem<CompType, ThermoType>::gradientType& dxdy)
{
    return;
}
/*
template<class CompType, class ThermoType>
void Foam::ISAT_chem<CompType, ThermoType>::leafData::set(const inputType& x)
{
    phi_ = x; Rphi_ = func(x); computeA(phi_, Rphi_, A_); EOA = Foam::ISAT_chem<CompType, ThermoType>::pISAT->tolerance() / max(fabs(A_), 1);
}
*/
/*
template<typename ...Args>
template<class CompType, class ThermoType>
void Foam::ISAT_chem<CompType, ThermoType>::leafData::set(const inputType& x, const outputType& y, Args ...args)
{
    phi_ = x; Rphi_ = y; computeA(phi_, Rphi_, A_, args...); EOA = Foam::ISAT_chem<CompType, ThermoType>::pISAT->tolerance() / max(fabs(A_), 1);
}
*/
/*
template<class CompType, class ThermoType>
void Foam::ISAT_chem<CompType, ThermoType>::leafData::computeA
(
    scalarSquareMatrix& A,
    const scalarField& Rphiq,
    const scalar rhoi,
    const scalar dt
)
{
    bool mechRedActive = pISAT->chemistry_.mechRed()->active();
    label speciesNumber = pISAT->chemistry_.nSpecie();
    scalarField Rcq(pISAT->chemistry_.nEqns() + nAdditionalEqns_ - 2);
    for (label i = 0; i < speciesNumber; i++)
    {
        label s2c = i;
        if (mechRedActive)
        {
            s2c = pISAT->chemistry_.simplifiedToCompleteIndex()[i];
        }
        Rcq[i] = rhoi * Rphiq[s2c] / pISAT->chemistry_.specieThermo()[s2c].W();
    }
    Rcq[speciesNumber] = Rphiq[Rphiq.size() - nAdditionalEqns_];
    Rcq[speciesNumber + 1] = Rphiq[Rphiq.size() - nAdditionalEqns_ + 1];
    if (pISAT->variableTimeStep())
    {
        Rcq[speciesNumber + 2] = Rphiq[Rphiq.size() - nAdditionalEqns_ + 2];
    }

    // Aaa is computed implicitly,
    // A is given by A = C(psi0, t0+dt), where C is obtained through solving
    // d/dt C(psi0, t) = J(psi(t))C(psi0, t)
    // If we solve it implicitly:
    // (C(psi0, t0+dt) - C(psi0, t0))/dt = J(psi(t0+dt))C(psi0, t0+dt)
    // The Jacobian is thus computed according to the mapping
    // C(psi0,t0+dt)*(I-dt*J(psi(t0+dt))) = C(psi0, t0)
    // A = C(psi0,t0)/(I-dt*J(psi(t0+dt)))
    // where C(psi0,t0) = I
    scalarField dcdt(speciesNumber + 2, Zero);
    pISAT->chemistry_.jacobian(runTime_.value(), Rcq, dcdt, A);

    // The jacobian is computed according to the molar concentration
    // the following conversion allows the code to use A with mass fraction
    for (label i = 0; i < speciesNumber; i++)
    {
        label si = i;

        if (mechRedActive)
        {
            si = pISAT->chemistry_.simplifiedToCompleteIndex()[i];
        }

        for (label j = 0; j < speciesNumber; j++)
        {
            label sj = j;
            if (mechRedActive)
            {
                sj = pISAT->chemistry_.simplifiedToCompleteIndex()[j];
            }
            A(i, j) *=
                -dt * pISAT->chemistry_.specieThermo()[si].W()
                / pISAT->chemistry_.specieThermo()[sj].W();
        }

        A(i, i) += 1;
        // Columns for pressure and temperature
        A(i, speciesNumber) *=
            -dt * pISAT->chemistry_.specieThermo()[si].W() / rhoi;
        A(i, speciesNumber + 1) *=
            -dt * pISAT->chemistry_.specieThermo()[si].W() / rhoi;
    }

    // For the temperature and pressure lines, ddc(dTdt)
    // should be converted in ddY(dTdt)
    for (label i = 0; i < speciesNumber; i++)
    {
        label si = i;
        if (mechRedActive)
        {
            si = pISAT->chemistry_.simplifiedToCompleteIndex()[i];
        }

        A(speciesNumber, i) *=
            -dt * rhoi / pISAT->chemistry_.specieThermo()[si].W();
        A(speciesNumber + 1, i) *=
            -dt * rhoi / pISAT->chemistry_.specieThermo()[si].W();
    }

    A(speciesNumber, speciesNumber) = -dt * A(speciesNumber, speciesNumber) + 1;

    A(speciesNumber + 1, speciesNumber + 1) =
        -dt * A(speciesNumber + 1, speciesNumber + 1) + 1;

    if (pISAT->variableTimeStep())
    {
        A(speciesNumber + 2, speciesNumber + 2) = 1;
    }

    // Inverse of (I-dt*J(psi(t0+dt)))
    LUscalarMatrix LUA(A);
    LUA.inv(A);

    // After inversion, lines of p and T are set to 0 except diagonal.  This
    // avoid skewness of the ellipsoid of accuracy and potential issues in the
    // binary tree.
    for (label i = 0; i < speciesNumber; i++)
    {
        A(speciesNumber, i) = 0;
        A(speciesNumber + 1, i) = 0;
    }
}


template<class CompType, class ThermoType>
bool Foam::ISAT_chem<CompType, ThermoType>::leafData::grow(const scalarField& phiq)
{
    scalarField dphi(phiq - phi());
    label dim = completeSpaceSize();
    label initNActiveSpecies(nActiveSpecies_);
    bool isMechRedActive = chemistry_.mechRed()->active();

    if (isMechRedActive)
    {
        label activeAdded(0);
        DynamicList<label> dimToAdd(0);

        // check if the difference of active species is lower than the maximum
        // number of new dimensions allowed
        for (label i = 0; i < completeSpaceSize() - nAdditionalEqns_; i++)
        {
            // first test if the current chemPoint has an inactive species
            // corresponding to an active one in the query point
            if
                (
                    completeToSimplifiedIndex_[i] == -1
                    && chemistry_.completeToSimplifiedIndex()[i] != -1
                    )
            {
                activeAdded++;
                dimToAdd.append(i);
            }
            // then test if an active species in the current chemPoint
            // corresponds to an inactive on the query side
            if
                (
                    completeToSimplifiedIndex_[i] != -1
                    && chemistry_.completeToSimplifiedIndex()[i] == -1
                    )
            {
                activeAdded++;
                // we don't need to add a new dimension but we count it to have
                // control on the difference through maxNumNewDim
            }
            // finally test if both points have inactive species but
            // with a dphi!=0
            if
                (
                    completeToSimplifiedIndex_[i] == -1
                    && chemistry_.completeToSimplifiedIndex()[i] == -1
                    && dphi[i] != 0
                    )
            {
                activeAdded++;
                dimToAdd.append(i);
            }
        }

        // if the number of added dimension is too large, growth fail
        if (activeAdded > maxNumNewDim_)
        {
            return false;
        }

        // the number of added dimension to the current chemPoint
        nActiveSpecies_ += dimToAdd.size();
        simplifiedToCompleteIndex_.setSize(nActiveSpecies_);
        forAll(dimToAdd, i)
        {
            label si = nActiveSpecies_ - dimToAdd.size() + i;
            // add the new active species
            simplifiedToCompleteIndex_[si] = dimToAdd[i];
            completeToSimplifiedIndex_[dimToAdd[i]] = si;
        }

        // update LT and A :
        //-add new column and line for the new active species
        //-transfer last two lines of the previous matrix (p and T) to the end
        //  (change the diagonal position)
        //-set all element of the new lines and columns to zero except diagonal
        //  (=1/(tolerance*scaleFactor))
        if (nActiveSpecies_ > initNActiveSpecies)
        {
            scalarSquareMatrix LTvar = LT_; // take a copy of LT_
            scalarSquareMatrix Avar = A_; // take a copy of A_
            LT_ = scalarSquareMatrix(nActiveSpecies_ + nAdditionalEqns_, Zero);
            A_ = scalarSquareMatrix(nActiveSpecies_ + nAdditionalEqns_, Zero);

            // write the initial active species
            for (label i = 0; i < initNActiveSpecies; i++)
            {
                for (label j = 0; j < initNActiveSpecies; j++)
                {
                    LT_(i, j) = LTvar(i, j);
                    A_(i, j) = Avar(i, j);
                }
            }

            // write the columns for temperature and pressure
            for (label i = 0; i < initNActiveSpecies; i++)
            {
                for (label j = 1; j >= 0; j--)
                {
                    LT_(i, nActiveSpecies_ + j) = LTvar(i, initNActiveSpecies + j);
                    A_(i, nActiveSpecies_ + j) = Avar(i, initNActiveSpecies + j);
                    LT_(nActiveSpecies_ + j, i) = LTvar(initNActiveSpecies + j, i);
                    A_(nActiveSpecies_ + j, i) = Avar(initNActiveSpecies + j, i);
                }
            }
            // end with the diagonal elements for temperature and pressure
            LT_(nActiveSpecies_, nActiveSpecies_) =
                LTvar(initNActiveSpecies, initNActiveSpecies);
            A_(nActiveSpecies_, nActiveSpecies_) =
                Avar(initNActiveSpecies, initNActiveSpecies);
            LT_(nActiveSpecies_ + 1, nActiveSpecies_ + 1) =
                LTvar(initNActiveSpecies + 1, initNActiveSpecies + 1);
            A_(nActiveSpecies_ + 1, nActiveSpecies_ + 1) =
                Avar(initNActiveSpecies + 1, initNActiveSpecies + 1);

            if (variableTimeStep())
            {
                LT_(nActiveSpecies_ + 2, nActiveSpecies_ + 2) =
                    LTvar(initNActiveSpecies + 2, initNActiveSpecies + 2);
                A_(nActiveSpecies_ + 2, nActiveSpecies_ + 2) =
                    Avar(initNActiveSpecies + 2, initNActiveSpecies + 2);
            }

            for (label i = initNActiveSpecies; i < nActiveSpecies_;i++)
            {
                LT_(i, i) =
                    1.0
                    / (tolerance_ * scaleFactor_[simplifiedToCompleteIndex_[i]]);
                A_(i, i) = 1;
            }
        }

        dim = nActiveSpecies_ + nAdditionalEqns_;
    }

    // beginning of grow algorithm
    scalarField phiTilde(dim, 0);
    scalar normPhiTilde = 0;
    // p' = L^T.(p-phi)

    for (label i = 0; i < dim; i++)
    {
        for (label j = i; j < dim - nAdditionalEqns_; j++)// LT is upper triangular
        {
            label sj = j;
            if (isMechRedActive)
            {
                sj = simplifiedToCompleteIndex_[j];
            }
            phiTilde[i] += LT_(i, j) * dphi[sj];
        }

        phiTilde[i] += LT_(i, dim - nAdditionalEqns_) * dphi[idT_];
        phiTilde[i] += LT_(i, dim - nAdditionalEqns_ + 1) * dphi[idp_];

        if (variableTimeStep())
        {
            phiTilde[i] += LT_(i, dim - nAdditionalEqns_ + 2) * dphi[iddeltaT_];
        }
        normPhiTilde += sqr(phiTilde[i]);
    }

    scalar invSqrNormPhiTilde = 1.0 / normPhiTilde;
    normPhiTilde = sqrt(normPhiTilde);

    // gamma = (1/|p'| - 1)/|p'|^2
    scalar gamma = (1 / normPhiTilde - 1) * invSqrNormPhiTilde;
    scalarField u(gamma * phiTilde);
    scalarField v(dim, 0);

    for (label i = 0; i < dim; i++)
    {
        for (label j = 0; j <= i;j++)
        {
            v[i] += phiTilde[j] * LT_(j, i);
        }
    }

    qrUpdate(LT_, dim, u, v);
    nGrowth_++;

    return true;
}


template<class CompType, class ThermoType>
void Foam::ISAT_chem<CompType, ThermoType>::leafData::retrieve
(
    const inputType& x, outputType& y
)
{
    label nEqns = pISAT->chemistry_.nEqns(); // Species, T, p
    bool mechRedActive = pISAT->chemistry_.mechRed()->active();
    Rphiq = Rphi_;
    scalarField dphi(x - phi_);
    const scalarSquareMatrix& gradientsMatrix = A_;
    List<label>& completeToSimplified(phi0->completeToSimplifiedIndex());

    // Rphiq[i]=Rphi0[i]+A(i, j)dphi[j]
    // where Aij is dRi/dphi_j
    for (label i = 0; i < nEqns - nAdditionalEqns_; i++)
    {
        if (mechRedActive)
        {
            label si = completeToSimplified[i];
            // The species is active
            if (si != -1)
            {
                for (label j = 0; j < nEqns - 2; j++)
                {
                    label sj = completeToSimplified[j];
                    if (sj != -1)
                    {
                        Rphiq[i] += gradientsMatrix(si, sj) * dphi[j];
                    }
                }
                Rphiq[i] +=
                    gradientsMatrix(si, phi0->nActiveSpecies()) * dphi[nEqns - 2];
                Rphiq[i] +=
                    gradientsMatrix(si, phi0->nActiveSpecies() + 1)
                    * dphi[nEqns - 1];

                if (pISAT->variableTimeStep())
                {
                    Rphiq[i] +=
                        gradientsMatrix(si, phi0->nActiveSpecies() + 2)
                        * dphi[nEqns];
                }

                // As we use an approximation of A, Rphiq should be checked for
                // negative values
                Rphiq[i] = max(0, Rphiq[i]);
            }
            // The species is not active A(i, j) = I(i, j)
            else
            {
                Rphiq[i] += dphi[i];
                Rphiq[i] = max(0, Rphiq[i]);
            }
        }
        else // Mechanism reduction is not active
        {
            for (label j = 0; j < nEqns; j++)
            {
                Rphiq[i] += gradientsMatrix(i, j) * dphi[j];
            }
            // As we use a first order gradient matrix, Rphiq should be checked
            // for negative values
            Rphiq[i] = max(0, Rphiq[i]);
        }
    }
}


template<class CompType, class ThermoType>
bool Foam::ISAT_chem<CompType, ThermoType>::leafData::inEOA(const scalarField& phiq)
{
    scalarField dphi(phiq - phi());
    bool isMechRedActive = chemistry_.mechRed()->active();
    label dim(0);
    if (isMechRedActive)
    {
        dim = nActiveSpecies_;
    }
    else
    {
        dim = completeSpaceSize() - nAdditionalEqns_;
    }

    scalar epsTemp = 0;
    List<scalar> propEps(completeSpaceSize(), scalar(0));

    for (label i = 0; i < completeSpaceSize() - nAdditionalEqns_; i++)
    {
        scalar temp = 0;

        // When mechanism reduction is inactive OR on active species multiply L
        // by dphi to get the distance in the active species direction else (for
        // inactive species), just multiply the diagonal element and dphi
        if
            (
                !(isMechRedActive)
                || (isMechRedActive && completeToSimplifiedIndex_[i] != -1)
                )
        {
            label si = (isMechRedActive) ? completeToSimplifiedIndex_[i] : i;

            for (label j = si; j < dim; j++)// LT is upper triangular
            {
                label sj = (isMechRedActive) ? simplifiedToCompleteIndex_[j] : j;
                temp += LT_(si, j) * dphi[sj];
            }

            temp += LT_(si, dim) * dphi[idT_];
            temp += LT_(si, dim + 1) * dphi[idp_];
            if (variableTimeStep())
            {
                temp += LT_(si, dim + 2) * dphi[iddeltaT_];
            }
        }
        else
        {
            temp = dphi[i] / (tolerance_ * scaleFactor_[i]);
        }

        epsTemp += sqr(temp);

        if (printProportion_)
        {
            propEps[i] = temp;
        }
    }

    // Temperature
    if (variableTimeStep())
    {
        epsTemp +=
            sqr
            (
                LT_(dim, dim) * dphi[idT_]
                + LT_(dim, dim + 1) * dphi[idp_]
                + LT_(dim, dim + 2) * dphi[iddeltaT_]
            );
    }
    else
    {
        epsTemp +=
            sqr
            (
                LT_(dim, dim) * dphi[idT_]
                + LT_(dim, dim + 1) * dphi[idp_]
            );
    }

    // Pressure
    if (variableTimeStep())
    {
        epsTemp +=
            sqr
            (
                LT_(dim + 1, dim + 1) * dphi[idp_]
                + LT_(dim + 1, dim + 2) * dphi[iddeltaT_]
            );
    }
    else
    {
        epsTemp += sqr(LT_(dim + 1, dim + 1) * dphi[idp_]);
    }

    if (variableTimeStep())
    {
        epsTemp += sqr(LT_[dim + 2][dim + 2] * dphi[iddeltaT_]);
    }

    if (printProportion_)
    {
        propEps[idT_] = sqr
        (
            LT_(dim, dim) * dphi[idT_]
            + LT_(dim, dim + 1) * dphi[idp_]
        );

        propEps[idp_] =
            sqr(LT_(dim + 1, dim + 1) * dphi[idp_]);

        if (variableTimeStep())
        {
            propEps[iddeltaT_] =
                sqr(LT_[dim + 2][dim + 2] * dphi[iddeltaT_]);
        }

    }

    if (sqrt(epsTemp) > 1 + tolerance_)
    {
        if (printProportion_)
        {
            scalar max = -1;
            label maxIndex = -1;
            for (label i = 0; i < completeSpaceSize(); i++)
            {
                if (max < propEps[i])
                {
                    max = propEps[i];
                    maxIndex = i;
                }
            }
            word propName;
            if (maxIndex >= completeSpaceSize() - nAdditionalEqns_)
            {
                if (maxIndex == idT_)
                {
                    propName = "T";
                }
                else if (maxIndex == idp_)
                {
                    propName = "p";
                }
                else if (maxIndex == iddeltaT_)
                {
                    propName = "deltaT";
                }
            }
            else
            {
                propName = chemistry_.Y()[maxIndex].member();
            }
            Info << "Direction maximum impact to error in ellipsoid: "
                << propName << endl;
            Info << "Proportion to the total error on the retrieve: "
                << max / (epsTemp + small) << endl;
        }
        return false;
    }
    else
    {
        return true;
    }
}

template<class CompType, class ThermoType>
void Foam::ISAT_chem<CompType, ThermoType>::nodeData::calcV
(
    const leafData& elementLeft,
    const leafData& elementRight,
    scalarField& v
)
{
    // LT is the transpose of the L matrix
    scalarSquareMatrix& LT = elementLeft.LT_;
    bool mechReductionActive = elementLeft.pISAT->chemistry().mechRed()->active();

    // Difference of composition in the full species domain
    scalarField phiDif(elementRight.phi_ - elementLeft.phi_);
    const scalarField& scaleFactor(elementLeft.scaleFactor_);
    scalar epsTol = elementLeft.tolerance_;

    // v = LT.T()*LT*phiDif
    for (label i = 0; i < elementLeft->completeSpaceSize(); i++)
    {
        label si = i;
        bool outOfIndexI = true;
        if (mechReductionActive)
        {
            if (i < elementLeft->completeSpaceSize() - nAdditionalEqns_)
            {
                si = elementLeft->completeToSimplifiedIndex()[i];
                outOfIndexI = (si == -1);
            }
            else // temperature and pressure
            {
                outOfIndexI = false;
                const label dif =
                    i - (elementLeft->completeSpaceSize() - nAdditionalEqns_);
                si = elementLeft->nActiveSpecies() + dif;
            }
        }
        if (!mechReductionActive || (mechReductionActive && !(outOfIndexI)))
        {
            v[i] = 0;
            for (label j = 0; j < elementLeft->completeSpaceSize(); j++)
            {
                label sj = j;
                bool outOfIndexJ = true;
                if (mechReductionActive)
                {
                    if (j < elementLeft->completeSpaceSize() - nAdditionalEqns_)
                    {
                        sj = elementLeft->completeToSimplifiedIndex()[j];
                        outOfIndexJ = (sj == -1);
                    }
                    else
                    {
                        outOfIndexJ = false;
                        const label dif =
                            j
                            - (
                                elementLeft->completeSpaceSize()
                                - nAdditionalEqns_
                                );
                        sj = elementLeft->nActiveSpecies() + dif;
                    }
                }
                if
                    (
                        !mechReductionActive
                        || (mechReductionActive && !(outOfIndexJ))
                        )
                {
                    // Since L is a lower triangular matrix k=0->min(i, j)
                    for (label k = 0; k <= min(si, sj); k++)
                    {
                        v[i] += LT(k, si) * LT(k, sj) * phiDif[j];
                    }
                }
            }
        }
        else
        {
            // When it is an inactive species the diagonal element of LT is
            // 1/(scaleFactor*epsTol)
            v[i] = phiDif[i] / sqr(scaleFactor[i] * epsTol);
        }
    }
}*/


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


// ************************************************************************* //
