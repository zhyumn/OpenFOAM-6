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

#include "ISAT_manager.H"
#include <iomanip>
#include "LUscalarMatrix.H"
#include "SVD.H"
#include "SortableList.H"
//#include<iostream>

namespace Foam
{
    namespace chemistryTabulationMethodSs
    {
        //template<>
        //label ISAT_manager<psiReactionThermo, gasHThermoPhysics>::out = 0;

        //template<>
        //label ISAT_chem<psiReactionThermo, gasHThermoPhysics>::out = 0;

        bool static_data::mechRedActive = true;
        label static_data::nAdditionalEqns_ = 1;
        bool static_data::variableTimeStep = true;
        label static_data::nEqns = 0;
        label static_data::completeSpaceSize_ = 0;
        scalar static_data::tolerance_ = 0.0;
        label static_data::idT_ = 0;
        label static_data::idp_ = 0;
        label static_data::iddeltaT_ = 0;
        //template <>
        //ISAT_manager<psiReactionThermo, gasHThermoPhysics> *ISAT_chem<psiReactionThermo, gasHThermoPhysics>::pISAT = nullptr;
        //template <class CompType, class ThermoType>
        /*         std::ostream &operator<<(std::ostream &out, leafData &A)
        {
            out << A.v << ", " << A.Rv;
            return out;
        } */
        //template <class CompType, class ThermoType>
        /*         std::ostream &operator<<(std::ostream &out, nodeData &A)
        {
            out << A.v;
            return out;
        } */
    }
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

//const dataType Foam::ISAT_manager::staticData();

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class CompType, class ThermoType>
Foam::chemistryTabulationMethodSs::ISAT_manager<CompType, ThermoType>::ISAT_manager(
    const dictionary &dict,
    STDACChemistryModel<CompType, ThermoType> &chemistry)
    : chemistryTabulationMethodS<CompType, ThermoType>(
          dict,
          chemistry),
      tree(*this,
           SUPstream::node_manager,
           label(readLabel(this->coeffsDict_.lookup("maxNLeafs")) * 1.2),
           readLabel(this->coeffsDict_.lookup("maxNLeafs")),
           SUPstream::Sync),
      runTime_(chemistry.time()),
      completeSpaceSize_(chemistry.nEqns() + ((this->variableTimeStep()) ? 1 : 0)),
      scaleFactor_(completeSpaceSize_, 1),
      maxNumNewDim_(this->coeffsDict_.lookupOrDefault("maxNumNewDim", 0)),
      chPMaxLifeTime_(
          this->coeffsDict_.lookupOrDefault("chPMaxLifeTime", INT_MAX)),
      maxGrowth_(this->coeffsDict_.lookupOrDefault("maxGrowth", INT_MAX)),
      minBalanceThreshold_(
          this->coeffsDict_.lookupOrDefault(
              "minBalanceThreshold", 0.1 * tree.maxNLeafs_)),
      maxDepthFactor_(
          this->coeffsDict_.lookupOrDefault(
              "maxDepthFactor",
              (tree.maxNLeafs_ - 1) / (log(scalar(tree.maxNLeafs_)) / log(2.0))))
{
    //this->maxNLeafs_ = readLabel(this->coeffsDict_.lookup("maxNLeafs"));

    //Foam::chemistryTabulationMethodSs::ISAT_chem<CompType, ThermoType>::pISAT = this;
    completeSpaceSize_ = chemistry.nEqns() + ((this->variableTimeStep()) ? 1 : 0);

    if (this->active_)
    {
        dictionary scaleDict(this->coeffsDict_.subDict("scaleFactor"));
        label Ysize = this->chemistry_.Y().size();
        scalar otherScaleFactor = readScalar(scaleDict.lookup("otherSpecies"));
        for (label i = 0; i < Ysize; i++)
        {
            if (!scaleDict.found(this->chemistry_.Y()[i].member()))
            {
                scaleFactor_[i] = otherScaleFactor;
            }
            else
            {
                scaleFactor_[i] =
                    readScalar(
                        scaleDict.lookup(this->chemistry_.Y()[i].member()));
            }
        }
        scaleFactor_[Ysize] = readScalar(scaleDict.lookup("Temperature"));
        scaleFactor_[Ysize + 1] = readScalar(scaleDict.lookup("Pressure"));
        if (this->variableTimeStep())
        {
            scaleFactor_[Ysize + 2] = readScalar(scaleDict.lookup("deltaT"));
        }
    }
    if (this->variableTimeStep())
    {
        nAdditionalEqns_ = 3;
        iddeltaT_ = completeSpaceSize_ - 1;
    }
    else
    {
        nAdditionalEqns_ = 2;
        iddeltaT_ = completeSpaceSize_;
    }

    idT_ = completeSpaceSize_ - nAdditionalEqns_;
    idp_ = completeSpaceSize_ - nAdditionalEqns_ + 1;

    printProportion_ = false;
    //maxNumNewDim_ = 0;

    if (this->log())
    {
        nRetrievedFile_ = chemistry.logFile("found_isat.out");
        nGrowthFile_ = chemistry.logFile("growth_isat.out");
        nAddFile_ = chemistry.logFile("add_isat.out");
        sizeFile_ = chemistry.logFile("size_isat.out");
        totalGrowthFile_ = chemistry.logFile("total_growth_isat.out");
    }
    static_data::mechRedActive = this->chemistry_.mechRed()->active();
    static_data::nAdditionalEqns_ = nAdditionalEqns_;
    static_data::variableTimeStep = this->variableTimeStep();
    static_data::nEqns = chemistry.nEqns();
    static_data::completeSpaceSize_ = completeSpaceSize_;
    static_data::tolerance_ = this->tolerance();
    static_data::idT_ = idT_;
    static_data::idp_ = idp_;
    static_data::iddeltaT_ = iddeltaT_;
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
template <class CompType, class ThermoType>
Foam::chemistryTabulationMethodSs::ISAT_manager<CompType, ThermoType>::~ISAT_manager()
{
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

/* template <class CompType, class ThermoType>
void Foam::chemistryTabulationMethodSs::ISAT_chem<CompType, ThermoType>::leafData::set(const scalarField &x, const scalarField &y, SUPstream::mpi_mutex &mem_mutex, const scalar rhoi, const scalar dt)
{
    scalarSquareMatrix A(0); //TODO
    phi_ = x;
    Rphi_ = y;

    computeA(x, y, A_, A, rhoi, dt, mem_mutex);
    //pISAT->computeA(x, y, A_, A, rhoi, dt, mem_mutex);
    computeLT(A, mem_mutex);
    //pISAT->computeLT(*this,A, mem_mutex);
    timeTag_ = pISAT->chemistry_.timeSteps();
    nGrowth_ = 0;
} */

template <class CompType, class ThermoType>
void Foam::chemistryTabulationMethodSs::ISAT_manager<CompType, ThermoType>::setLeaf(
    leafData &leaf_input,
    const scalarField &x,
    const scalarField &y,
    SUPstream::mpi_mutex &mem_mutex,
    const scalar rhoi,
    const scalar dt)
{
    scalarSquareMatrix A(0); //TODO
    leaf_input.phi_ = x;
    leaf_input.Rphi_ = y;

    //computeA(x, y, A_, A, rhoi, dt, mem_mutex);
    this->computeA(x, y, leaf_input.A_, A, rhoi, dt, mem_mutex);
    //computeLT(A, mem_mutex);
    this->computeLT(leaf_input, A, mem_mutex);
    leaf_input.timeTag_ = this->chemistry_.timeSteps();
    leaf_input.nGrowth_ = 0;
}

/* template <class CompType, class ThermoType>
void Foam::chemistryTabulationMethodSs::ISAT_chem<CompType, ThermoType>::leafData::computeA(
    const scalarField &phiq,
    const scalarField &Rphiq,
    SSquareMatrix<scalar> &A_out,
    scalarSquareMatrix &A,
    const scalar rhoi,
    const scalar dt, SUPstream::mpi_mutex &mem_mutex)
{

    bool mechRedActive = pISAT->chemistry_.mechRed()->active();
    label speciesNumber = pISAT->chemistry_.nSpecie();

    const label &nAdditionalEqns_ = pISAT->nAdditionalEqns_;
    const Time &runTime_ = pISAT->runTime_;

    label Asize = pISAT->chemistry_.nEqns() + nAdditionalEqns_ - 2;
    A.setSize(Asize);
    mem_mutex.lock();
    A_out.init(Asize);
    mem_mutex.unlock();

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
                -dt * pISAT->chemistry_.specieThermo()[si].W() / pISAT->chemistry_.specieThermo()[sj].W();
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

    for (label i = 0; i < Asize; i++)
        for (label j = 0; j < Asize; j++)
        {
            A_out(i, j) = A(i, j);
        }
} */

template <class CompType, class ThermoType>
void Foam::chemistryTabulationMethodSs::ISAT_manager<CompType, ThermoType>::computeA(
    const scalarField &phiq,
    const scalarField &Rphiq,
    SSquareMatrix<scalar> &A_out,
    scalarSquareMatrix &A,
    const scalar rhoi,
    const scalar dt, SUPstream::mpi_mutex &mem_mutex)
{

    bool mechRedActive = this->chemistry_.mechRed()->active();
    label speciesNumber = this->chemistry_.nSpecie();

    const label &nAdditionalEqns_ = this->nAdditionalEqns_;
    const Time &runTime_ = this->runTime_;

    label Asize = this->chemistry_.nEqns() + nAdditionalEqns_ - 2;
    A.setSize(Asize);
    mem_mutex.lock();
    A_out.init(Asize);
    mem_mutex.unlock();

    scalarField Rcq(this->chemistry_.nEqns() + nAdditionalEqns_ - 2);
    for (label i = 0; i < speciesNumber; i++)
    {
        label s2c = i;
        if (mechRedActive)
        {
            s2c = this->chemistry_.simplifiedToCompleteIndex()[i];
        }
        Rcq[i] = rhoi * Rphiq[s2c] / this->chemistry_.specieThermo()[s2c].W();
    }
    Rcq[speciesNumber] = Rphiq[Rphiq.size() - nAdditionalEqns_];
    Rcq[speciesNumber + 1] = Rphiq[Rphiq.size() - nAdditionalEqns_ + 1];
    if (this->variableTimeStep())
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
    this->chemistry_.jacobian(runTime_.value(), Rcq, dcdt, A);

    // The jacobian is computed according to the molar concentration
    // the following conversion allows the code to use A with mass fraction
    for (label i = 0; i < speciesNumber; i++)
    {
        label si = i;

        if (mechRedActive)
        {
            si = this->chemistry_.simplifiedToCompleteIndex()[i];
        }

        for (label j = 0; j < speciesNumber; j++)
        {
            label sj = j;
            if (mechRedActive)
            {
                sj = this->chemistry_.simplifiedToCompleteIndex()[j];
            }
            A(i, j) *=
                -dt * this->chemistry_.specieThermo()[si].W() / this->chemistry_.specieThermo()[sj].W();
        }

        A(i, i) += 1;
        // Columns for pressure and temperature
        A(i, speciesNumber) *=
            -dt * this->chemistry_.specieThermo()[si].W() / rhoi;
        A(i, speciesNumber + 1) *=
            -dt * this->chemistry_.specieThermo()[si].W() / rhoi;
    }

    // For the temperature and pressure lines, ddc(dTdt)
    // should be converted in ddY(dTdt)
    for (label i = 0; i < speciesNumber; i++)
    {
        label si = i;
        if (mechRedActive)
        {
            si = this->chemistry_.simplifiedToCompleteIndex()[i];
        }

        A(speciesNumber, i) *=
            -dt * rhoi / this->chemistry_.specieThermo()[si].W();
        A(speciesNumber + 1, i) *=
            -dt * rhoi / this->chemistry_.specieThermo()[si].W();
    }

    A(speciesNumber, speciesNumber) = -dt * A(speciesNumber, speciesNumber) + 1;

    A(speciesNumber + 1, speciesNumber + 1) = -dt * A(speciesNumber + 1, speciesNumber + 1) + 1;

    if (this->variableTimeStep())
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

    for (label i = 0; i < Asize; i++)
        for (label j = 0; j < Asize; j++)
        {
            A_out(i, j) = A(i, j);
        }
}

/* template <class CompType, class ThermoType>
void Foam::chemistryTabulationMethodSs::ISAT_chem<CompType, ThermoType>::leafData::computeLT(
    scalarSquareMatrix &A, SUPstream::mpi_mutex &mem_mutex)
{
    const scalar &tolerance_ = pISAT->tolerance();
    const label &nAdditionalEqns_ = pISAT->nAdditionalEqns_;
    const label &completeSpaceSize = pISAT->completeSpaceSize_;
    const label &iddeltaT_ = pISAT->iddeltaT_;
    STDACChemistryModel<CompType, ThermoType> &chemistry = pISAT->chemistry_;

    for (label i = 0; i < completeSpaceSize; i++)
    {
        scaleFactor_[i] = pISAT->scaleFactor_[i];
    }
    if (pISAT->variableTimeStep())
    {
        scaleFactor_[iddeltaT_] *= phi_[iddeltaT_] / tolerance_;
    }
    // const label& idT_ = pISAT->idT_;
    //const label& idp_ = pISAT->idp_;
    nActiveSpecies_ = chemistry.mechRed()->NsSimp();
    bool isMechRedActive = chemistry.mechRed()->active();
    if (isMechRedActive)
    {
        for (label i = 0; i < completeSpaceSize - nAdditionalEqns_; i++)
        {
            completeToSimplifiedIndex_[i] =
                chemistry.completeToSimplifiedIndex()[i];
        }
        for (label i = 0; i < nActiveSpecies_; i++)
        {
            simplifiedToCompleteIndex_[i] =
                chemistry.simplifiedToCompleteIndex()[i];
        }
    }

    label reduOrCompDim = completeSpaceSize;
    if (isMechRedActive)
    {
        reduOrCompDim = nActiveSpecies_ + nAdditionalEqns_;
    }

    // SVD decomposition A = U*D*V^T
    SVD svdA(A);

    scalarDiagonalMatrix D(reduOrCompDim);
    const scalarDiagonalMatrix &S = svdA.S();

    // Replace the value of vector D by max(D, 1/2), first ISAT paper
    for (label i = 0; i < reduOrCompDim; i++)
    {
        D[i] = max(S[i], 0.5);
    }

    // Rebuild A with max length, tol and scale factor before QR decomposition
    scalarRectangularMatrix Atilde(reduOrCompDim);

    // Result stored in Atilde
    multiply(Atilde, svdA.U(), D, svdA.V().T());

    for (label i = 0; i < reduOrCompDim; i++)
    {
        for (label j = 0; j < reduOrCompDim; j++)
        {
            label compi = i;

            if (isMechRedActive)
            {
                compi = simplifiedToCompleteIndex(i);
            }

            // SF*A/tolerance
            // (where SF is diagonal with inverse of scale factors)
            // SF*A is the same as dividing each line by the scale factor
            // corresponding to the species of this line
            //Pout << "hahahhahaha=" << (tolerance_) << endl;
            Atilde(i, j) /= (tolerance_ * pISAT->scaleFactor_[compi]);
        }
    }

    // The object LT_ (the transpose of the Q) describe the EOA, since we have
    // A^T B^T B A that should be factorized into L Q^T Q L^T and is set in the
    // qrDecompose function
    scalarSquareMatrix LT = scalarSquareMatrix(Atilde);

    qrDecompose(reduOrCompDim, LT);
    mem_mutex.lock();
    LT_.init(reduOrCompDim);
    mem_mutex.unlock();
    for (label i = 0; i < reduOrCompDim; i++)
    {
        for (label j = 0; j < reduOrCompDim; j++)
        {
            LT_(i, j) = LT(i, j);
        }
    }
} */

template <class CompType, class ThermoType>
void Foam::chemistryTabulationMethodSs::ISAT_manager<CompType, ThermoType>::computeLT(
    leafData &leaf_input, scalarSquareMatrix &A, SUPstream::mpi_mutex &mem_mutex)
{
    const scalar &tolerance_ = this->tolerance();
    const label &nAdditionalEqns_ = this->nAdditionalEqns_;
    const label &completeSpaceSize = this->completeSpaceSize_;
    const label &iddeltaT_ = this->iddeltaT_;
    STDACChemistryModel<CompType, ThermoType> &chemistry = this->chemistry_;

    for (label i = 0; i < completeSpaceSize; i++)
    {
        leaf_input.scaleFactor_[i] = this->scaleFactor_[i];
    }
    if (this->variableTimeStep())
    {
        leaf_input.scaleFactor_[iddeltaT_] *= leaf_input.phi_[iddeltaT_] / tolerance_;
    }
    // const label& idT_ = this->idT_;
    //const label& idp_ = this->idp_;
    leaf_input.nActiveSpecies_ = chemistry.mechRed()->NsSimp();
    bool isMechRedActive = chemistry.mechRed()->active();
    if (isMechRedActive)
    {
        for (label i = 0; i < completeSpaceSize - nAdditionalEqns_; i++)
        {
            leaf_input.completeToSimplifiedIndex_[i] = chemistry.completeToSimplifiedIndex()[i];
        }
        for (label i = 0; i < leaf_input.nActiveSpecies_; i++)
        {
            leaf_input.simplifiedToCompleteIndex_[i] = chemistry.simplifiedToCompleteIndex()[i];
        }
    }

    label reduOrCompDim = completeSpaceSize;
    if (isMechRedActive)
    {
        reduOrCompDim = leaf_input.nActiveSpecies_ + nAdditionalEqns_;
    }

    // SVD decomposition A = U*D*V^T
    SVD svdA(A);

    scalarDiagonalMatrix D(reduOrCompDim);
    const scalarDiagonalMatrix &S = svdA.S();

    // Replace the value of vector D by max(D, 1/2), first ISAT paper
    for (label i = 0; i < reduOrCompDim; i++)
    {
        D[i] = max(S[i], 0.5);
    }

    // Rebuild A with max length, tol and scale factor before QR decomposition
    scalarRectangularMatrix Atilde(reduOrCompDim);

    // Result stored in Atilde
    multiply(Atilde, svdA.U(), D, svdA.V().T());

    for (label i = 0; i < reduOrCompDim; i++)
    {
        for (label j = 0; j < reduOrCompDim; j++)
        {
            label compi = i;

            if (isMechRedActive)
            {
                compi = leaf_input.simplifiedToCompleteIndex(i);
            }

            // SF*A/tolerance
            // (where SF is diagonal with inverse of scale factors)
            // SF*A is the same as dividing each line by the scale factor
            // corresponding to the species of this line
            //Pout << "hahahhahaha=" << (tolerance_) << endl;
            Atilde(i, j) /= (tolerance_ * this->scaleFactor_[compi]);
        }
    }

    // The object LT_ (the transpose of the Q) describe the EOA, since we have
    // A^T B^T B A that should be factorized into L Q^T Q L^T and is set in the
    // qrDecompose function
    scalarSquareMatrix LT = scalarSquareMatrix(Atilde);

    qrDecompose(reduOrCompDim, LT);
    mem_mutex.lock();
    leaf_input.LT_.init(reduOrCompDim);
    mem_mutex.unlock();
    for (label i = 0; i < reduOrCompDim; i++)
    {
        for (label j = 0; j < reduOrCompDim; j++)
        {
            leaf_input.LT_(i, j) = LT(i, j);
        }
    }
}

template <class CompType, class ThermoType>
void Foam::chemistryTabulationMethodSs::ISAT_manager<CompType, ThermoType>::qrDecompose(
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


/* template <class CompType, class ThermoType>
bool Foam::chemistryTabulationMethodSs::ISAT_chem<CompType, ThermoType>::leafData::grow(const scalarField &phiq, SUPstream::mpi_mutex &mem_mutex)
{

    const label &completeSpaceSize = pISAT->completeSpaceSize_;
    const label &nAdditionalEqns_ = pISAT->nAdditionalEqns_;
    STDACChemistryModel<CompType, ThermoType> &chemistry_ = pISAT->chemistry_;

    const label &idT_ = pISAT->idT_;
    const label &idp_ = pISAT->idp_;
    const label &iddeltaT_ = pISAT->iddeltaT_;
    const scalar &tolerance_ = pISAT->tolerance();

    scalarField dphi(completeSpaceSize);
    for (label i = 0; i < completeSpaceSize; i++)
    {
        dphi[i] = phiq[i] - phi_[i];
    }

    //scalarField dphi(phiq - phi());
    label dim = pISAT->completeSpaceSize_;

    label initNActiveSpecies(nActiveSpecies_);
    bool isMechRedActive = pISAT->chemistry_.mechRed()->active();

    //SSquareMatrix<scalar> &LT_ref = LT_;
    //SSquareMatrix<scalar> &A_ref = A_;
    mem_mutex.lock();
    FSSquareMatrix<scalar> LT_ref(LT_);
    FSSquareMatrix<scalar> A_ref(A_);
    mem_mutex.unlock();

    if (isMechRedActive)
    {
        label activeAdded(0);
        DynamicList<label> dimToAdd(0);

        // check if the difference of active species is lower than the maximum
        // number of new dimensions allowed
        for (label i = 0; i < completeSpaceSize - nAdditionalEqns_; i++)
        {
            // first test if the current chemPoint has an inactive species
            // corresponding to an active one in the query point
            if (
                completeToSimplifiedIndex_[i] == -1 && chemistry_.completeToSimplifiedIndex()[i] != -1)
            {
                activeAdded++;
                dimToAdd.append(i);
            }
            // then test if an active species in the current chemPoint
            // corresponds to an inactive on the query side
            if (
                completeToSimplifiedIndex_[i] != -1 && chemistry_.completeToSimplifiedIndex()[i] == -1)
            {
                activeAdded++;
                // we don't need to add a new dimension but we count it to have
                // control on the difference through maxNumNewDim
            }
            // finally test if both points have inactive species but
            // with a dphi!=0
            if (
                completeToSimplifiedIndex_[i] == -1 && chemistry_.completeToSimplifiedIndex()[i] == -1 && dphi[i] != 0)
            {
                activeAdded++;
                dimToAdd.append(i);
            }
        }

        // if the number of added dimension is too large, growth fail
        if (activeAdded > pISAT->maxNumNewDim_)
        {
            return false;
        }

        // the number of added dimension to the current chemPoint
        nActiveSpecies_ += dimToAdd.size();
        mem_mutex.lock();
        simplifiedToCompleteIndex_.setSize(nActiveSpecies_);
        mem_mutex.unlock();
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
            //label initSize = initNActiveSpecies + nAdditionalEqns_;

            scalarSquareMatrix LTvar(LT_.size_);
            scalarSquareMatrix Avar(A_.size_);

            for (label i = 0; i < A_.size_; i++)
            {
                for (label j = 0; j < A_.size_; j++)
                {
                    Avar(i, j) = A_ref(i, j);
                    A_ref(i, j) = 0;
                }
            }
            for (label i = 0; i < LT_.size_; i++)
            {
                for (label j = 0; j < LT_.size_; j++)
                {
                    LTvar(i, j) = LT_ref(i, j);
                    LT_ref(i, j) = 0;
                }
            }
            mem_mutex.lock();
            LT_ref.setSize(nActiveSpecies_ + nAdditionalEqns_);
            A_ref.setSize(nActiveSpecies_ + nAdditionalEqns_);
            mem_mutex.unlock();
            for (label i = 0; i < LT_.size_; i++)
            {
                for (label j = 0; j < LT_.size_; j++)
                {
                    A_ref(i, j) = 0;
                    LT_ref(i, j) = 0;
                }
            }

            //scalarSquareMatrix LTvar = LT_; // take a copy of LT_  //Todo
            //scalarSquareMatrix Avar = A_; // take a copy of A_
            //LT_ = scalarSquareMatrix(nActiveSpecies_ + nAdditionalEqns_, Zero);
            //A_ = scalarSquareMatrix(nActiveSpecies_ + nAdditionalEqns_, Zero);

            // write the initial active species

            for (label i = 0; i < initNActiveSpecies; i++)
            {
                for (label j = 0; j < initNActiveSpecies; j++)
                {
                    LT_ref(i, j) = LTvar(i, j);
                    A_ref(i, j) = Avar(i, j);
                }
            }

            // write the columns for temperature and pressure
            for (label i = 0; i < initNActiveSpecies; i++)
            {
                for (label j = 1; j >= 0; j--)
                {
                    LT_ref(i, nActiveSpecies_ + j) = LTvar(i, initNActiveSpecies + j);
                    A_ref(i, nActiveSpecies_ + j) = Avar(i, initNActiveSpecies + j);
                    LT_ref(nActiveSpecies_ + j, i) = LTvar(initNActiveSpecies + j, i);
                    A_ref(nActiveSpecies_ + j, i) = Avar(initNActiveSpecies + j, i);
                }
            }

            // end with the diagonal elements for temperature and pressure
            LT_ref(nActiveSpecies_, nActiveSpecies_) =
                LTvar(initNActiveSpecies, initNActiveSpecies);
            A_ref(nActiveSpecies_, nActiveSpecies_) =
                Avar(initNActiveSpecies, initNActiveSpecies);
            LT_ref(nActiveSpecies_ + 1, nActiveSpecies_ + 1) =
                LTvar(initNActiveSpecies + 1, initNActiveSpecies + 1);
            A_ref(nActiveSpecies_ + 1, nActiveSpecies_ + 1) =
                Avar(initNActiveSpecies + 1, initNActiveSpecies + 1);

            if (pISAT->variableTimeStep())
            {
                LT_ref(nActiveSpecies_ + 2, nActiveSpecies_ + 2) =
                    LTvar(initNActiveSpecies + 2, initNActiveSpecies + 2);
                A_ref(nActiveSpecies_ + 2, nActiveSpecies_ + 2) =
                    Avar(initNActiveSpecies + 2, initNActiveSpecies + 2);
            }

            for (label i = initNActiveSpecies; i < nActiveSpecies_; i++)
            {
                LT_ref(i, i) =
                    1.0 / (tolerance_ * scaleFactor_[simplifiedToCompleteIndex_[i]]);
                A_ref(i, i) = 1;
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
        for (label j = i; j < dim - nAdditionalEqns_; j++) // LT is upper triangular
        {
            label sj = j;
            if (isMechRedActive)
            {
                sj = simplifiedToCompleteIndex_[j];
            }
            phiTilde[i] += LT_ref(i, j) * dphi[sj];
        }

        phiTilde[i] += LT_ref(i, dim - nAdditionalEqns_) * dphi[idT_];
        phiTilde[i] += LT_ref(i, dim - nAdditionalEqns_ + 1) * dphi[idp_];

        if (pISAT->variableTimeStep())
        {
            phiTilde[i] += LT_ref(i, dim - nAdditionalEqns_ + 2) * dphi[iddeltaT_];
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
        for (label j = 0; j <= i; j++)
        {
            v[i] += phiTilde[j] * LT_ref(j, i);
        }
    }

    qrUpdate(LT_ref, dim, u, v);
    nGrowth_++;

    return true;
} */

template <class CompType, class ThermoType>
bool Foam::chemistryTabulationMethodSs::ISAT_manager<CompType, ThermoType>::grow(
    leafData &leaf_input, const scalarField &phiq, SUPstream::mpi_mutex &mem_mutex)
{

    const label &completeSpaceSize = this->completeSpaceSize_;
    const label &nAdditionalEqns_ = this->nAdditionalEqns_;
    STDACChemistryModel<CompType, ThermoType> &chemistry_ = this->chemistry_;

    const label &idT_ = this->idT_;
    const label &idp_ = this->idp_;
    const label &iddeltaT_ = this->iddeltaT_;
    const scalar &tolerance_ = this->tolerance();
    label &nActiveSpecies_ = leaf_input.nActiveSpecies_;
    scalarField dphi(completeSpaceSize);
    for (label i = 0; i < completeSpaceSize; i++)
    {
        dphi[i] = phiq[i] - leaf_input.phi_[i];
    }

    //scalarField dphi(phiq - phi());
    label dim = this->completeSpaceSize_;

    label initNActiveSpecies(leaf_input.nActiveSpecies_);
    bool isMechRedActive = this->chemistry_.mechRed()->active();

    //SSquareMatrix<scalar> &LT_ref = LT_;
    //SSquareMatrix<scalar> &A_ref = A_;
    mem_mutex.lock();
    FSSquareMatrix<scalar> LT_ref(leaf_input.LT_);
    FSSquareMatrix<scalar> A_ref(leaf_input.A_);
    mem_mutex.unlock();

    if (isMechRedActive)
    {
        label activeAdded(0);
        DynamicList<label> dimToAdd(0);

        // check if the difference of active species is lower than the maximum
        // number of new dimensions allowed
        for (label i = 0; i < completeSpaceSize - nAdditionalEqns_; i++)
        {
            // first test if the current chemPoint has an inactive species
            // corresponding to an active one in the query point
            if (
                leaf_input.completeToSimplifiedIndex_[i] == -1 && chemistry_.completeToSimplifiedIndex()[i] != -1)
            {
                activeAdded++;
                dimToAdd.append(i);
            }
            // then test if an active species in the current chemPoint
            // corresponds to an inactive on the query side
            if (
                leaf_input.completeToSimplifiedIndex_[i] != -1 && chemistry_.completeToSimplifiedIndex()[i] == -1)
            {
                activeAdded++;
                // we don't need to add a new dimension but we count it to have
                // control on the difference through maxNumNewDim
            }
            // finally test if both points have inactive species but
            // with a dphi!=0
            if (
                leaf_input.completeToSimplifiedIndex_[i] == -1 && chemistry_.completeToSimplifiedIndex()[i] == -1 && dphi[i] != 0)
            {
                activeAdded++;
                dimToAdd.append(i);
            }
        }

        // if the number of added dimension is too large, growth fail
        if (activeAdded > this->maxNumNewDim_)
        {
            return false;
        }

        // the number of added dimension to the current chemPoint
        nActiveSpecies_ += dimToAdd.size();
        mem_mutex.lock();
        leaf_input.simplifiedToCompleteIndex_.setSize(nActiveSpecies_);
        mem_mutex.unlock();
        forAll(dimToAdd, i)
        {
            label si = nActiveSpecies_ - dimToAdd.size() + i;
            // add the new active species
            leaf_input.simplifiedToCompleteIndex_[si] = dimToAdd[i];
            leaf_input.completeToSimplifiedIndex_[dimToAdd[i]] = si;
        }

        // update LT and A :
        //-add new column and line for the new active species
        //-transfer last two lines of the previous matrix (p and T) to the end
        //  (change the diagonal position)
        //-set all element of the new lines and columns to zero except diagonal
        //  (=1/(tolerance*scaleFactor))
        if (nActiveSpecies_ > initNActiveSpecies)
        {
            //label initSize = initNActiveSpecies + nAdditionalEqns_;

            scalarSquareMatrix LTvar(leaf_input.LT_.size_);
            scalarSquareMatrix Avar(leaf_input.A_.size_);

            for (label i = 0; i < leaf_input.A_.size_; i++)
            {
                for (label j = 0; j < leaf_input.A_.size_; j++)
                {
                    Avar(i, j) = A_ref(i, j);
                    A_ref(i, j) = 0;
                }
            }
            for (label i = 0; i < leaf_input.LT_.size_; i++)
            {
                for (label j = 0; j < leaf_input.LT_.size_; j++)
                {
                    LTvar(i, j) = LT_ref(i, j);
                    LT_ref(i, j) = 0;
                }
            }
            mem_mutex.lock();
            LT_ref.setSize(nActiveSpecies_ + nAdditionalEqns_);
            A_ref.setSize(nActiveSpecies_ + nAdditionalEqns_);
            mem_mutex.unlock();
            for (label i = 0; i < leaf_input.LT_.size_; i++)
            {
                for (label j = 0; j < leaf_input.LT_.size_; j++)
                {
                    A_ref(i, j) = 0;
                    LT_ref(i, j) = 0;
                }
            }

            //scalarSquareMatrix LTvar = LT_; // take a copy of LT_  //Todo
            //scalarSquareMatrix Avar = A_; // take a copy of A_
            //LT_ = scalarSquareMatrix(nActiveSpecies_ + nAdditionalEqns_, Zero);
            //A_ = scalarSquareMatrix(nActiveSpecies_ + nAdditionalEqns_, Zero);

            // write the initial active species

            for (label i = 0; i < initNActiveSpecies; i++)
            {
                for (label j = 0; j < initNActiveSpecies; j++)
                {
                    LT_ref(i, j) = LTvar(i, j);
                    A_ref(i, j) = Avar(i, j);
                }
            }

            // write the columns for temperature and pressure
            for (label i = 0; i < initNActiveSpecies; i++)
            {
                for (label j = 1; j >= 0; j--)
                {
                    LT_ref(i, nActiveSpecies_ + j) = LTvar(i, initNActiveSpecies + j);
                    A_ref(i, nActiveSpecies_ + j) = Avar(i, initNActiveSpecies + j);
                    LT_ref(nActiveSpecies_ + j, i) = LTvar(initNActiveSpecies + j, i);
                    A_ref(nActiveSpecies_ + j, i) = Avar(initNActiveSpecies + j, i);
                }
            }

            // end with the diagonal elements for temperature and pressure
            LT_ref(nActiveSpecies_, nActiveSpecies_) =
                LTvar(initNActiveSpecies, initNActiveSpecies);
            A_ref(nActiveSpecies_, nActiveSpecies_) =
                Avar(initNActiveSpecies, initNActiveSpecies);
            LT_ref(nActiveSpecies_ + 1, nActiveSpecies_ + 1) =
                LTvar(initNActiveSpecies + 1, initNActiveSpecies + 1);
            A_ref(nActiveSpecies_ + 1, nActiveSpecies_ + 1) =
                Avar(initNActiveSpecies + 1, initNActiveSpecies + 1);

            if (static_data::variableTimeStep)
            {
                LT_ref(nActiveSpecies_ + 2, nActiveSpecies_ + 2) =
                    LTvar(initNActiveSpecies + 2, initNActiveSpecies + 2);
                A_ref(nActiveSpecies_ + 2, nActiveSpecies_ + 2) =
                    Avar(initNActiveSpecies + 2, initNActiveSpecies + 2);
            }

            for (label i = initNActiveSpecies; i < nActiveSpecies_; i++)
            {
                LT_ref(i, i) =
                    1.0 / (tolerance_ * scaleFactor_[leaf_input.simplifiedToCompleteIndex_[i]]);
                A_ref(i, i) = 1;
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
        for (label j = i; j < dim - nAdditionalEqns_; j++) // LT is upper triangular
        {
            label sj = j;
            if (isMechRedActive)
            {
                sj = leaf_input.simplifiedToCompleteIndex_[j];
            }
            phiTilde[i] += LT_ref(i, j) * dphi[sj];
        }

        phiTilde[i] += LT_ref(i, dim - nAdditionalEqns_) * dphi[idT_];
        phiTilde[i] += LT_ref(i, dim - nAdditionalEqns_ + 1) * dphi[idp_];

        if (static_data::variableTimeStep)
        {
            phiTilde[i] += LT_ref(i, dim - nAdditionalEqns_ + 2) * dphi[iddeltaT_];
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
        for (label j = 0; j <= i; j++)
        {
            v[i] += phiTilde[j] * LT_ref(j, i);
        }
    }

    leaf_input.qrUpdate(LT_ref, dim, u, v);
    leaf_input.nGrowth_++;

    return true;
}


//template <class CompType, class ThermoType>


template <class CompType, class ThermoType>
void Foam::chemistryTabulationMethodSs::ISAT_manager<CompType, ThermoType>::balance()
{
    //this->valid();
    scalarField mean(completeSpaceSize_, 0.0);

    //1) walk through the entire tree by starting with the tree's most left
    // chemPoint
    SharedPointer<typename parallelISAT<ISAT_manager<CompType, ThermoType>>::Leaf> x = tree.treeMin();
    List<SharedPointer<typename parallelISAT<ISAT_manager<CompType, ThermoType>>::Leaf>> chemPoints(tree.size_leaf_);
    label chPi = 0;

    //2) compute the mean composition
    while (x.notNULL())
    {
        const SList<scalar> &phij = x->phi();
        for (int i = 0; i < phij.size_; i++)
        {
            mean[i] += phij[i];
        }

        chemPoints[chPi++] = x;
        x = tree.treeSuccessor(x);
    }
    x = chemPoints[0];
    for (int i = 0; i < x->phi().size_; i++)
    {
        mean[i] /= tree.size_leaf_;
    }
    //mean /= size_;
    //3) compute the variance for each space direction
    List<scalar> variance(completeSpaceSize_, 0.0);
    forAll(chemPoints, j)
    {
        const SList<scalar> &phij = chemPoints[j]->phi();
        forAll(variance, vi)
        {
            variance[vi] += sqr(phij[vi] - mean[vi]);
        }
    }

    //4) analyze what is the direction of the maximal variance
    scalar maxVariance(-1.0);
    label maxDir(-1);
    forAll(variance, vi)
    {
        if (maxVariance < variance[vi])
        {
            maxVariance = variance[vi];
            maxDir = vi;
        }
    }
    // maxDir indicates the direction of maximum variance
    // we create the new root node by taking the two extreme points
    // in this direction if these extreme points were not deleted in the
    // cleaning that come before the balance function they are still important
    // and the tree should therefore take them into account
    SortableList<scalar> phiMaxDir(chemPoints.size(), 0.0);
    forAll(chemPoints, j)
    {
        phiMaxDir[j] = chemPoints[j]->phi()[maxDir];
    }

    phiMaxDir.sort();
    // delete reference to all node since the tree is reshaped
    tree.deleteAllNode();
    tree.root_.setNULL();

    // add the node for the two extremum
    SharedPointer<typename parallelISAT<ISAT_manager<CompType, ThermoType>>::Node> newNode = tree.node_manager.New(completeSpaceSize_);
    if (newNode.isNULL())
    {
        FatalErrorInFunction
            << "run out of memory"
            << exit(FatalError);
    }
    newNode->set(*chemPoints[phiMaxDir.indices()[0]], *chemPoints[phiMaxDir.indices()[phiMaxDir.size() - 1]]);
    //newNode->parent_.offset = -1;
    newNode->parent_.setNULL();
    newNode->leafLeft_ = chemPoints[phiMaxDir.indices()[0]];
    newNode->leafRight_ = chemPoints[phiMaxDir.indices()[phiMaxDir.size() - 1]];
    newNode->nodeLeft_.setNULL();
    newNode->nodeRight_.setNULL();

    /*label newNode = new bn
    (
        chemPoints[phiMaxDir.indices()[0]],
        chemPoints[phiMaxDir.indices()[phiMaxDir.size() - 1]],
        nullptr
    );*/
    tree.root_ = newNode;

    chemPoints[phiMaxDir.indices()[0]]->node_2 = newNode;
    chemPoints[phiMaxDir.indices()[phiMaxDir.size() - 1]]->node_2 = newNode;

    for (label cpi = 1; cpi < chemPoints.size() - 1; cpi++)
    {
        SharedPointer<typename parallelISAT<ISAT_manager<CompType, ThermoType>>::Leaf> phi0;

        phi0 = tree.binaryTreeSearch(chemPoints[phiMaxDir.indices()[cpi]]->phi(), tree.root_);

        // add the chemPoint
        SharedPointer<typename parallelISAT<ISAT_manager<CompType, ThermoType>>::Node> nodeToAdd = tree.node_manager.New(completeSpaceSize_);

        if (nodeToAdd.isNULL())
        {
            FatalErrorInFunction
                << "run out of memory"
                << exit(FatalError);
        }
        nodeToAdd->set(*phi0, *chemPoints[phiMaxDir.indices()[cpi]]);
        nodeToAdd->parent_ = phi0->node_2;
        nodeToAdd->leafLeft_ = phi0;
        nodeToAdd->leafRight_ = chemPoints[phiMaxDir.indices()[cpi]];
        nodeToAdd->nodeLeft_.setNULL();
        nodeToAdd->nodeRight_.setNULL();
        //new bn(phi0, chemPoints[phiMaxDir.indices()[cpi]], phi0->node());
        // make the parent of phi0 point to the newly created node
        tree.insertNode(phi0, nodeToAdd);
        phi0->node_2 = nodeToAdd;
        chemPoints[phiMaxDir.indices()[cpi]]->node_2 = nodeToAdd;
    }
}

template <class CompType, class ThermoType>
bool Foam::chemistryTabulationMethodSs::ISAT_manager<CompType, ThermoType>::cleanAndBalance()
{
    if (tree.size_leaf_ < 1)
    {
        return false;
    }
    bool treeModified(false);

    // Check all chemPoints to see if we need to delete some of the chemPoints
    // according to the ellapsed time and number of growths
    SharedPointer<typename parallelISAT<ISAT_manager<CompType, ThermoType>>::Leaf> x = tree.treeMin();

    while (x.notNULL())
    {

        SharedPointer<typename parallelISAT<ISAT_manager<CompType, ThermoType>>::Leaf> xtmp = tree.treeSuccessor(x);

        scalar elapsedTimeSteps = this->chemistry_.timeSteps() - x->timeTag();

        //(elapsedTimeSteps > chPMaxLifeTime_);

        if ((elapsedTimeSteps > chPMaxLifeTime_) || (x->nGrowth() > maxGrowth_))
        {
            tree.deleteLeaf(x);
            treeModified = true;
        }
        x = xtmp;
    }

    //MRUList_.clear();

    // Check if the tree should be balanced according to criterion:
    //  -the depth of the tree bigger than a*log2(size), log2(size) being the
    //      ideal depth (e.g. 4 leafs can be stored in a tree of depth 2)
    if (
        size() > minBalanceThreshold_ && tree.depth() >
                                             maxDepthFactor_ * log2(scalar(size())))
    {
        balance();
        treeModified = true;
    }

    // Return a bool to specify if the tree structure has been modified and is
    // now below the user specified limit (true if not full)
    return (treeModified && !tree.isFull());
}
// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

// ************************************************************************* //
