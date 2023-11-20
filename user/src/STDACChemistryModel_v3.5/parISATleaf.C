/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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

#include "parISATleaf.H"
#include "SVD.H"
#include "IOstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
struct static_data;
void Foam::chemistryTabulationMethodSs::leafData::qrDecompose(
    const label nCols,
    scalarSquareMatrix &R)
{
    scalarField c(nCols);
    scalarField d(nCols);
    scalar scale, sigma, sum;

    for (label k = 0; k < nCols - 1; k++)
    {
        scale = 0;
        for (label i = k; i < nCols; i++)
        {
            scale = max(scale, mag(R(i, k)));
        }
        if (scale == 0)
        {
            c[k] = d[k] = 0;
        }
        else
        {
            for (label i = k; i < nCols; i++)
            {
                R(i, k) /= scale;
            }
            sum = 0;
            for (label i = k; i < nCols; i++)
            {
                sum += sqr(R(i, k));
            }
            sigma = sign(R(k, k)) * sqrt(sum);
            R(k, k) += sigma;
            c[k] = sigma * R(k, k);
            d[k] = -scale * sigma;
            for (label j = k + 1; j < nCols; j++)
            {
                sum = 0;
                for (label i = k; i < nCols; i++)
                {
                    sum += R(i, k) * R(i, j);
                }
                scalar tau = sum / c[k];
                for (label i = k; i < nCols; i++)
                {
                    R(i, j) -= tau * R(i, k);
                }
            }
        }
    }

    d[nCols - 1] = R(nCols - 1, nCols - 1);

    // form R
    for (label i = 0; i < nCols; i++)
    {
        R(i, i) = d[i];
        for (label j = 0; j < i; j++)
        {
            R(i, j) = 0;
        }
    }
}

//template <class CompType, class ThermoType>
Foam::label Foam::chemistryTabulationMethodSs::leafData::
    simplifiedToCompleteIndex(
        const label i)
{
    //const label &nAdditionalEqns_ = pISAT->nAdditionalEqns_;
    //const label &completeSpaceSize_ = pISAT->completeSpaceSize_;

    const label &nAdditionalEqns_ = static_data::nAdditionalEqns_;
    const label &completeSpaceSize_ = static_data::completeSpaceSize_;

    if (i < nActiveSpecies_)
    {
        return simplifiedToCompleteIndex_[i];
    }
    else if (i == nActiveSpecies_)
    {
        return completeSpaceSize_ - nAdditionalEqns_;
    }
    else if (i == nActiveSpecies_ + 1)
    {
        return completeSpaceSize_ - nAdditionalEqns_ + 1;
    }
    else if (static_data::variableTimeStep && (i == nActiveSpecies_ + 2))
    {
        return completeSpaceSize_ - nAdditionalEqns_ + 2;
    }
    else
    {
        return -1;
    }
}

//template <class CompType, class ThermoType>
bool Foam::chemistryTabulationMethodSs::leafData::inEOA(const scalarField &phiq)
{
    //const label &completeSpaceSize = pISAT->completeSpaceSize_;
    //const label &nAdditionalEqns_ = pISAT->nAdditionalEqns_;
    //const scalar &tolerance_ = pISAT->tolerance();
    //const label &idT_ = pISAT->idT_;
    //const label &idp_ = pISAT->idp_;
    //const label &iddeltaT_ = pISAT->iddeltaT_;
    //const bool printProportion_ = pISAT->printProportion_;

    const label &completeSpaceSize = static_data::completeSpaceSize_;
    const label &nAdditionalEqns_ = static_data::nAdditionalEqns_;
    const scalar &tolerance_ = static_data::tolerance_;
    const label &idT_ = static_data::idT_;
    const label &idp_ = static_data::idp_;
    const label &iddeltaT_ = static_data::iddeltaT_;

    scalarField dphi(completeSpaceSize);
    for (label i = 0; i < completeSpaceSize; i++)
    {
        dphi[i] = phiq[i] - phi_[i];
    }
    //scalarField dphi(phiq-phi());
    //bool isMechRedActive = pISAT->chemistry_.mechRed()->active();
    bool isMechRedActive = static_data::mechRedActive;
    label dim(0);
    if (isMechRedActive)
    {
        dim = nActiveSpecies_;
    }
    else
    {
        dim = completeSpaceSize - nAdditionalEqns_;
    }

    scalar epsTemp = 0;
    //List<scalar> propEps(completeSpaceSize, scalar(0));
    //SSquareMatrix<scalar> &LT_ref = LT_;
    FSSquareMatrix<scalar> LT_ref(LT_);

    for (label i = 0; i < completeSpaceSize - nAdditionalEqns_; i++)
    {
        scalar temp = 0;

        // When mechanism reduction is inactive OR on active species multiply L
        // by dphi to get the distance in the active species direction else (for
        // inactive species), just multiply the diagonal element and dphi
        if (
            !(isMechRedActive) || (isMechRedActive && completeToSimplifiedIndex_[i] != -1))
        {
            label si = (isMechRedActive) ? completeToSimplifiedIndex_[i] : i;

            for (label j = si; j < dim; j++) // LT is upper triangular
            {
                label sj = (isMechRedActive) ? simplifiedToCompleteIndex_[j] : j;
                temp += LT_ref(si, j) * dphi[sj];
            }

            temp += LT_ref(si, dim) * dphi[idT_];
            temp += LT_ref(si, dim + 1) * dphi[idp_];
            if (static_data::variableTimeStep)
            {
                temp += LT_ref(si, dim + 2) * dphi[iddeltaT_];
            }
        }
        else
        {
            temp = dphi[i] / (tolerance_ * scaleFactor_[i]);
        }
        //epsTemp0 += sqr(temp);
        epsTemp += sqr(temp);

        /*         if (printProportion_)
        {
            propEps[i] = temp;
        } */
    }

    // Temperature
    if (static_data::variableTimeStep)
    {
        //epsTemp1 +=
        //    sqr(
        //        LT_(dim, dim) * dphi[idT_]
        //        + LT_(dim, dim + 1) * dphi[idp_]
        //        + LT_(dim, dim + 2) * dphi[iddeltaT_]
        //    );
        epsTemp +=
            sqr(
                LT_ref(dim, dim) * dphi[idT_] + LT_ref(dim, dim + 1) * dphi[idp_] + LT_ref(dim, dim + 2) * dphi[iddeltaT_]);
    }
    else
    {
        epsTemp +=
            sqr(
                LT_ref(dim, dim) * dphi[idT_] + LT_ref(dim, dim + 1) * dphi[idp_]);
    }

    // Pressure
    if (static_data::variableTimeStep)
    {
        //epsTemp2 += dphi[idp_];
        epsTemp +=
            sqr(
                LT_ref(dim + 1, dim + 1) * dphi[idp_] + LT_ref(dim + 1, dim + 2) * dphi[iddeltaT_]);
    }
    else
    {
        epsTemp += sqr(LT_ref(dim + 1, dim + 1) * dphi[idp_]);
    }

    if (static_data::variableTimeStep)
    {
        //epsTemp3 += sqr(LT_(dim + 2, dim + 2) * dphi[iddeltaT_]);
        epsTemp += sqr(LT_ref(dim + 2, dim + 2) * dphi[iddeltaT_]);
    }

    /*     if (printProportion_)
    {
        propEps[idT_] = sqr(
            LT_ref(dim, dim) * dphi[idT_] + LT_ref(dim, dim + 1) * dphi[idp_]);

        propEps[idp_] =
            sqr(LT_ref(dim + 1, dim + 1) * dphi[idp_]);

        if (pISAT->variableTimeStep())
        {
            propEps[iddeltaT_] =
                sqr(LT_ref(dim + 2, dim + 2) * dphi[iddeltaT_]);
        }
    } */
    //bool zz = sqrt(epsTemp) > 1 + tolerance_;

    if (sqrt(epsTemp) > 1 + tolerance_)
    {
        /*         if (printProportion_)
        {
            scalar max = -1;
            label maxIndex = -1;
            for (label i = 0; i < completeSpaceSize; i++)
            {
                if (max < propEps[i])
                {
                    max = propEps[i];
                    maxIndex = i;
                }
            }
            word propName;
            if (maxIndex >= completeSpaceSize - nAdditionalEqns_)
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
                propName = pISAT->chemistry_.Y()[maxIndex].member();
            }
            Info << "Direction maximum impact to error in ellipsoid: "
                 << propName << endl;
            Info << "Proportion to the total error on the retrieve: "
                 << max / (epsTemp + small) << endl;
        } */

        return false;
    }
    else
    {
        return true;
    }
}


//template <class CompType, class ThermoType>
//template <class M>
void Foam::chemistryTabulationMethodSs::leafData::qrUpdate(
    FSSquareMatrix<scalar> &R,
    const label n,
    const Foam::scalarField &u,
    const Foam::scalarField &v)
{
    label k;

    scalarField w(u);
    for (k = n - 1; k >= 0; k--)
    {
        if (w[k] != 0)
        {
            break;
        }
    }

    if (k < 0)
    {
        k = 0;
    }

    for (label i = k - 1; i >= 0; i--)
    {
        rotate(R, i, w[i], -w[i + 1], n);
        if (w[i] == 0)
        {
            w[i] = mag(w[i + 1]);
        }
        else if (mag(w[i]) > mag(w[i + 1]))
        {
            w[i] = mag(w[i]) * sqrt(1.0 + sqr(w[i + 1] / w[i]));
        }
        else
        {
            w[i] = mag(w[i + 1]) * sqrt(1.0 + sqr(w[i] / w[i + 1]));
        }
    }

    for (label i = 0; i < n; i++)
    {
        R(0, i) += w[0] * v[i];
    }

    for (label i = 0; i < k; i++)
    {
        rotate(R, i, R(i, i), -R(i + 1, i), n);
    }
}

//template <class CompType, class ThermoType>
template <class M>
void Foam::chemistryTabulationMethodSs::leafData::rotate(
    M &R,
    const label i,
    const scalar a,
    const scalar b,
    label n)
{
    scalar c, fact, s, w, y;
    if (a == 0)
    {
        c = 0;
        s = (b >= 0 ? 1 : -1);
    }
    else if (mag(a) > mag(b))
    {
        fact = b / a;
        c = sign(a) / sqrt(1.0 + sqr(fact));
        s = fact * c;
    }
    else
    {
        fact = a / b;
        s = sign(b) / sqrt(1.0 + sqr(fact));
        c = fact * s;
    }
    for (label j = i; j < n; j++)
    {
        y = R(i, j);
        w = R(i + 1, j);
        R(i, j) = c * y - s * w;
        R(i + 1, j) = s * y + c * w;
    }
}

//template <class CompType, class ThermoType>
void Foam::chemistryTabulationMethodSs::leafData::retrieve(
    const scalarField &x,
    scalarField &Rphiq)
{
    //timeTag_ = pISAT->chemistry_.timeSteps();
    //const label &nAdditionalEqns_ = pISAT->nAdditionalEqns_;
    const label &nAdditionalEqns_ = static_data::nAdditionalEqns_;
    //label nEqns = pISAT->chemistry_.nEqns(); // Species, T, p
    label nEqns = static_data::nEqns; // Species, T, p
    //bool mechRedActive = pISAT->chemistry_.mechRed()->active();
    bool mechRedActive = static_data::mechRedActive;
    scalarField dphi(x.size());
    for (label i = 0; i < x.size(); i++)
    {
        Rphiq[i] = Rphi_[i];
        dphi[i] = x[i] - phi_[i];
    }

    //const SSquareMatrix<scalar> &gradientsMatrix = A_;
    FSSquareMatrix<scalar> gradientsMatrix(A_);

    SList<label> &completeToSimplified = completeToSimplifiedIndex_;

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
                    gradientsMatrix(si, nActiveSpecies_) * dphi[nEqns - 2];
                Rphiq[i] +=
                    gradientsMatrix(si, nActiveSpecies_ + 1) * dphi[nEqns - 1];

                if (static_data::variableTimeStep) // (pISAT->variableTimeStep()) //(static_data::variableTimeStep)
                {
                    Rphiq[i] +=
                        gradientsMatrix(si, nActiveSpecies_ + 2) * dphi[nEqns];
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

//template <class CompType, class ThermoType>
bool Foam::chemistryTabulationMethodSs::leafData::checkSolution(
    const scalarField &phiq,
    const scalarField &Rphiq)
{
    scalar eps2 = 0;
    //const label &iddeltaT_ = pISAT->iddeltaT_;
    //const label &idT_ = pISAT->idT_;
    //const label &idp_ = pISAT->idp_;
    //const scalar &tolerance_ = pISAT->tolerance();
    //const label &completeSpaceSize = pISAT->completeSpaceSize_;
    //const label &nAdditionalEqns_ = pISAT->nAdditionalEqns_;

    const label &iddeltaT_ = static_data::iddeltaT_;
    const label &idT_ = static_data::idT_;
    const label &idp_ = static_data::idp_;
    const scalar &tolerance_ = static_data::tolerance_;
    const label &completeSpaceSize = static_data::completeSpaceSize_;
    const label &nAdditionalEqns_ = static_data::nAdditionalEqns_;
    scalarField dR(phiq.size());
    scalarField dphi(phiq.size());

    for (label i = 0; i < phiq.size(); i++)
    {
        dR[i] = Rphiq[i] - Rphi_[i];
        dphi[i] = phiq[i] - phi_[i];
    }

    const SList<scalar> &scaleFactorV = scaleFactor_;
    //const SSquareMatrix<scalar> &Avar = A_;
    FSSquareMatrix<scalar> Avar(A_);
    //bool isMechRedActive = pISAT->chemistry_.mechRed()->active();
    bool isMechRedActive = static_data::mechRedActive;
    scalar dRl = 0;
    label dim = completeSpaceSize - 2;
    if (isMechRedActive)
    {
        dim = nActiveSpecies_;
    }

    // Since we build only the solution for the species, T and p are not
    // included
    for (label i = 0; i < completeSpaceSize - nAdditionalEqns_; i++)
    {
        dRl = 0;
        if (isMechRedActive)
        {
            label si = completeToSimplifiedIndex_[i];

            // If this species is active
            if (si != -1)
            {
                for (label j = 0; j < dim; j++)
                {
                    label sj = simplifiedToCompleteIndex_[j];
                    dRl += Avar(si, j) * dphi[sj];
                }
                dRl += Avar(si, nActiveSpecies_) * dphi[idT_];
                dRl += Avar(si, nActiveSpecies_ + 1) * dphi[idp_];
                if (static_data::variableTimeStep)
                {
                    dRl += Avar(si, nActiveSpecies_ + 2) * dphi[iddeltaT_];
                }
            }
            else
            {
                dRl = dphi[i];
            }
        }
        else
        {
            for (label j = 0; j < completeSpaceSize; j++)
            {
                dRl += Avar(i, j) * dphi[j];
            }
        }
        eps2 += sqr((dR[i] - dRl) / scaleFactorV[i]);
    }

    eps2 = sqrt(eps2);

    if (eps2 > tolerance_)
    {
        return false;
    }
    else
    {
        // if the solution is in the ellipsoid of accuracy
        return true;
    }
}

void Foam::chemistryTabulationMethodSs::nodeData::calcV(
    leafData &elementLeft,
    leafData &elementRight,
    SList<scalar> &v, SUPstream::mpi_mutex &mem_mutex)
{

    //const label &nAdditionalEqns_ = pISAT->nAdditionalEqns_;
    const label &nAdditionalEqns_ = static_data::nAdditionalEqns_;
    // LT is the transpose of the L matrix
    //gradientType_old &LT = elementLeft.LT_old;
    //SSquareMatrix<scalar> &LT = elementLeft.LT_;
    mem_mutex.lock();
    FSSquareMatrix<scalar> LT(elementLeft.LT_);
    mem_mutex.unlock();
    //bool mechReductionActive = pISAT->chemistry_.mechRed()->active();
    bool mechReductionActive = static_data::mechRedActive;

    // Difference of composition in the full species domain
    scalarField phiDif(elementRight.phi_.size_);

    for (label i = 0; i < phiDif.size(); i++)
    {
        phiDif[i] = elementRight.phi_[i] - elementLeft.phi_[i];
    }

    const SList<scalar> &scaleFactor = elementLeft.scaleFactor_;
    //const scalar &epsTol = pISAT->tolerance();
    const scalar &epsTol = static_data::tolerance_;

    // v = LT.T()*LT*phiDif
    for (label i = 0; i < static_data::completeSpaceSize_; i++)
    {
        label si = i;
        bool outOfIndexI = true;
        if (mechReductionActive)
        {
            if (i < static_data::completeSpaceSize_ - nAdditionalEqns_)
            {
                si = elementLeft.completeToSimplifiedIndex_[i];
                outOfIndexI = (si == -1);
            }
            else // temperature and pressure
            {
                outOfIndexI = false;
                const label dif =
                    i - (static_data::completeSpaceSize_ - nAdditionalEqns_);
                si = elementLeft.nActiveSpecies_ + dif;
            }
        }
        if (!mechReductionActive || (mechReductionActive && !(outOfIndexI)))
        {
            v[i] = 0;
            for (label j = 0; j < static_data::completeSpaceSize_; j++)
            {
                label sj = j;
                bool outOfIndexJ = true;
                if (mechReductionActive)
                {
                    if (j < static_data::completeSpaceSize_ - nAdditionalEqns_)
                    {
                        sj = elementLeft.completeToSimplifiedIndex_[j];
                        outOfIndexJ = (sj == -1);
                    }
                    else
                    {
                        outOfIndexJ = false;
                        const label dif =
                            j - (static_data::completeSpaceSize_ - nAdditionalEqns_);
                        sj = elementLeft.nActiveSpecies_ + dif;
                    }
                }
                if (
                    !mechReductionActive || (mechReductionActive && !(outOfIndexJ)))
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
}

//template <class CompType, class ThermoType>
void Foam::chemistryTabulationMethodSs::nodeData::calcV(
    leafData &elementLeft,
    leafData &elementRight,
    SList<scalar> &v)
{

    const label &nAdditionalEqns_ = static_data::nAdditionalEqns_;
    // LT is the transpose of the L matrix
    //gradientType_old &LT = elementLeft.LT_old;
    //SSquareMatrix<scalar> &LT = elementLeft.LT_;
    FSSquareMatrix<scalar> LT(elementLeft.LT_);
    bool mechReductionActive = static_data::mechRedActive;

    // Difference of composition in the full species domain
    scalarField phiDif(elementRight.phi_.size_);

    for (label i = 0; i < phiDif.size(); i++)
    {
        phiDif[i] = elementRight.phi_[i] - elementLeft.phi_[i];
    }

    const SList<scalar> &scaleFactor = elementLeft.scaleFactor_;
    const scalar &epsTol = static_data::tolerance_;

    // v = LT.T()*LT*phiDif
    for (label i = 0; i < static_data::completeSpaceSize_; i++)
    {
        label si = i;
        bool outOfIndexI = true;
        if (mechReductionActive)
        {
            if (i < static_data::completeSpaceSize_ - nAdditionalEqns_)
            {
                si = elementLeft.completeToSimplifiedIndex_[i];
                outOfIndexI = (si == -1);
            }
            else // temperature and pressure
            {
                outOfIndexI = false;
                const label dif =
                    i - (static_data::completeSpaceSize_ - nAdditionalEqns_);
                si = elementLeft.nActiveSpecies_ + dif;
            }
        }
        if (!mechReductionActive || (mechReductionActive && !(outOfIndexI)))
        {
            v[i] = 0;
            for (label j = 0; j < static_data::completeSpaceSize_; j++)
            {
                label sj = j;
                bool outOfIndexJ = true;
                if (mechReductionActive)
                {
                    if (j < static_data::completeSpaceSize_ - nAdditionalEqns_)
                    {
                        sj = elementLeft.completeToSimplifiedIndex_[j];
                        outOfIndexJ = (sj == -1);
                    }
                    else
                    {
                        outOfIndexJ = false;
                        const label dif =
                            j - (static_data::completeSpaceSize_ - nAdditionalEqns_);
                        sj = elementLeft.nActiveSpecies_ + dif;
                    }
                }
                if (
                    !mechReductionActive || (mechReductionActive && !(outOfIndexJ)))
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
}
//template <class CompType, class ThermoType>
Foam::scalar Foam::chemistryTabulationMethodSs::nodeData::calcA(
    const leafData &elementLeft,
    const leafData &elementRight)
{
    scalarField phih(elementLeft.phi_.size_);
    //scalarField phih((elementLeft->phi() + elementRight->phi())/2);
    for (label i = 0; i < phih.size(); i++)
    {
        phih[i] = (elementLeft.phi_[i] + elementRight.phi_[i]) / 2;
    }
    scalar a = 0;
    for (label i = 0; i < phih.size(); i++)
    {
        a += v_[i] * phih[i];
    }

    return a;
}
// ************************************************************************* //
