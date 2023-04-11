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

#include "STDACChemistryModel.H"
#include "UniformField.H"
#include "localEulerDdtScheme.H"
#include "clockTime.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class ReactionThermo, class ThermoType>
Foam::STDACChemistryModel<ReactionThermo, ThermoType>::STDACChemistryModel(
    ReactionThermo &thermo)
    : StandardChemistryModel<ReactionThermo, ThermoType>(thermo),
      variableTimeStep_(
          this->mesh().time().controlDict().lookupOrDefault(
              "adjustTimeStep",
              false) ||
          fv::localEulerDdt::enabled(this->mesh())),
      timeSteps_(0),
      NsDAC_(this->nSpecie_),
      completeC_(this->nSpecie_, 0),
      reactionsDisabled_(this->reactions_.size(), false),
      specieComp_(this->nSpecie_),
      completeToSimplifiedIndex_(this->nSpecie_, -1),
      simplifiedToCompleteIndex_(this->nSpecie_),
      tabulationResults_(
          IOobject(
              thermo.phasePropertyName("TabulationResults"),
              this->time().timeName(),
              this->mesh(),
              IOobject::NO_READ,
              IOobject::AUTO_WRITE),
          this->mesh(),
          scalar(0)),
      n_block(SUPstream::node_manager.size * 2),
      //sja_(SUPstream::node_manager, n_block,sizeof(jobInput)*Batch_Size),
      //ja(&sja_()),
      Batch_Size(this->subDict("tabulation").lookupOrDefault("batchSize", 100)),
      ja(SUPstream::node_manager, n_block, sizeof(jobArray) + sizeof(jobInput) * (Batch_Size - 1)),
      ja_loc(Batch_Size),
      sspare_cpu(SUPstream::node_manager), spare_cpu(sspare_cpu()),
      sempty_head(SUPstream::node_manager), sempty_tail(SUPstream::node_manager),
      empty_head(sempty_head()), empty_tail(sempty_tail()),
      sfilled_head(SUPstream::node_manager), sfilled_tail(SUPstream::node_manager),
      filled_head(sfilled_head()), filled_tail(sfilled_tail()),
      sfinished_head(SUPstream::node_manager, SUPstream::node_manager.size), finished_head(&sfinished_head()),
      rank(SUPstream::node_manager.rank),
      l_empty(SUPstream::node_manager, SUPstream::Sync),
      l_filled(SUPstream::node_manager, SUPstream::Sync),
      l_finished(SUPstream::node_manager, SUPstream::Sync),
      l_receiver(SUPstream::node_manager, SUPstream::Sync), l_sender(SUPstream::node_manager, SUPstream::Sync)
{
    SUPstream::Sync();
    n_cpu = SUPstream::node_manager.size;
    if (SUPstream::node_manager.rank == 0)
    {
        spare_cpu = n_cpu;
        for (int i = 0; i < n_block; i++)
        {
            ja[i].next = i + 1;
            //ja[i].filled = false;
            //ja[i].finished = false;
        }
        ja[n_block - 1].next = -1;
        empty_head = 0;
        empty_tail = n_block - 1;
        filled_head = -1;
        filled_tail = -1;

        //l_receiver.lock();
    }
    finished_head[rank] = -1;
    SUPstream::Sync();
    loadBalance_ = this->subDict("tabulation").lookupOrDefault("loadBalance", true);
    Info << "STDACChemistryModel!!!! variableTimeStep_=" << variableTimeStep_ << endl;
    basicSpecieMixture &composition = this->thermo().composition();

    // Store the species composition according to the species index
    speciesTable speciesTab = composition.species();

    const HashTable<List<specieElement>> &specComp =
        dynamicCast<const reactingMixture<ThermoType> &>(this->thermo())
            .specieComposition();

    forAll(specieComp_, i)
    {
        specieComp_[i] = specComp[this->Y()[i].member()];
    }

    mechRed_ = chemistryReductionMethodS<ReactionThermo, ThermoType>::New(
        *this,
        *this);

    // When the mechanism reduction method is used, the 'active' flag for every
    // species should be initialized (by default 'active' is true)
    if (mechRed_->active())
    {
        forAll(this->Y(), i)
        {
            IOobject header(
                this->Y()[i].name(),
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ);

            // Check if the species file is provided, if not set inactive
            // and NO_WRITE
            if (!header.typeHeaderOk<volScalarField>(true))
            {
                composition.setInactive(i);
                this->Y()[i].writeOpt() = IOobject::NO_WRITE;
            }
        }
    }

    tabulation_ = chemistryTabulationMethodS<ReactionThermo, ThermoType>::New(
        *this,
        *this);

    if (mechRed_->log())
    {
        cpuReduceFile_ = logFile("cpu_reduce.out");
        nActiveSpeciesFile_ = logFile("nActiveSpecies.out");
    }

    if (tabulation_->log())
    {
        cpuAddFile_ = logFile("cpu_add.out");
        cpuGrowFile_ = logFile("cpu_grow.out");
        cpuRetrieveFile_ = logFile("cpu_retrieve.out");
        cpuTotalFile_ = logFile("cpu_total.out");
    }

    if (mechRed_->log() || tabulation_->log())
    {
        cpuSolveFile_ = logFile("cpu_solve.out");
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class ReactionThermo, class ThermoType>
Foam::STDACChemistryModel<ReactionThermo, ThermoType>::~STDACChemistryModel()
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class ReactionThermo, class ThermoType>
void Foam::STDACChemistryModel<ReactionThermo, ThermoType>::omega(
    const scalarField &c, // Contains all species even when mechRed is active
    const scalar T,
    const scalar p,
    scalarField &dcdt) const
{
    const bool reduced = mechRed_->active();

    scalar pf, cf, pr, cr;
    label lRef, rRef;

    dcdt = Zero;

    forAll(this->reactions_, i)
    {
        if (!reactionsDisabled_[i])
        {
            const Reaction<ThermoType> &R = this->reactions_[i];

            scalar omegai = R.omega(
                p, T, c, pf, cf, lRef, pr, cr, rRef);

            forAll(R.lhs(), s)
            {
                label si = R.lhs()[s].index;
                if (reduced)
                {
                    si = completeToSimplifiedIndex_[si];
                }

                const scalar sl = R.lhs()[s].stoichCoeff;
                dcdt[si] -= sl * omegai;
            }
            forAll(R.rhs(), s)
            {
                label si = R.rhs()[s].index;
                if (reduced)
                {
                    si = completeToSimplifiedIndex_[si];
                }

                const scalar sr = R.rhs()[s].stoichCoeff;
                dcdt[si] += sr * omegai;
            }
        }
    }
}

template <class ReactionThermo, class ThermoType>
Foam::scalar Foam::STDACChemistryModel<ReactionThermo, ThermoType>::omega(
    const Reaction<ThermoType> &R,
    const scalarField &c, // Contains all species even when mechRed is active
    const scalar T,
    const scalar p,
    scalar &pf,
    scalar &cf,
    label &lRef,
    scalar &pr,
    scalar &cr,
    label &rRef) const
{
    const scalar kf = R.kf(p, T, c);
    const scalar kr = R.kr(kf, p, T, c);

    const label Nl = R.lhs().size();
    const label Nr = R.rhs().size();

    label slRef = 0;
    lRef = R.lhs()[slRef].index;

    pf = kf;
    for (label s = 1; s < Nl; s++)
    {
        const label si = R.lhs()[s].index;

        if (c[si] < c[lRef])
        {
            const scalar exp = R.lhs()[slRef].exponent;
            pf *= pow(max(c[lRef], 0), exp);
            lRef = si;
            slRef = s;
        }
        else
        {
            const scalar exp = R.lhs()[s].exponent;
            pf *= pow(max(c[si], 0), exp);
        }
    }
    cf = max(c[lRef], 0);

    {
        const scalar exp = R.lhs()[slRef].exponent;
        if (exp < 1)
        {
            if (cf > small)
            {
                pf *= pow(cf, exp - 1);
            }
            else
            {
                pf = 0;
            }
        }
        else
        {
            pf *= pow(cf, exp - 1);
        }
    }

    label srRef = 0;
    rRef = R.rhs()[srRef].index;

    // Find the matrix element and element position for the rhs
    pr = kr;
    for (label s = 1; s < Nr; s++)
    {
        const label si = R.rhs()[s].index;
        if (c[si] < c[rRef])
        {
            const scalar exp = R.rhs()[srRef].exponent;
            pr *= pow(max(c[rRef], 0), exp);
            rRef = si;
            srRef = s;
        }
        else
        {
            const scalar exp = R.rhs()[s].exponent;
            pr *= pow(max(c[si], 0), exp);
        }
    }
    cr = max(c[rRef], 0);

    {
        const scalar exp = R.rhs()[srRef].exponent;
        if (exp < 1)
        {
            if (cr > small)
            {
                pr *= pow(cr, exp - 1);
            }
            else
            {
                pr = 0;
            }
        }
        else
        {
            pr *= pow(cr, exp - 1);
        }
    }

    return pf * cf - pr * cr;
}

template <class ReactionThermo, class ThermoType>
void Foam::STDACChemistryModel<ReactionThermo, ThermoType>::derivatives(
    const scalar time,
    const scalarField &c,
    scalarField &dcdt) const
{
    const bool reduced = mechRed_->active();

    const scalar T = c[this->nSpecie_];
    const scalar p = c[this->nSpecie_ + 1];

    if (reduced)
    {
        // When using DAC, the ODE solver submit a reduced set of species the
        // complete set is used and only the species in the simplified mechanism
        // are updated
        this->c_ = completeC_;

        // Update the concentration of the species in the simplified mechanism
        // the other species remain the same and are used only for third-body
        // efficiencies
        for (label i = 0; i < NsDAC_; i++)
        {
            this->c_[simplifiedToCompleteIndex_[i]] = max(c[i], 0);
        }
    }
    else
    {
        for (label i = 0; i < this->nSpecie(); i++)
        {
            this->c_[i] = max(c[i], 0);
        }
    }

    omega(this->c_, T, p, dcdt);

    // Constant pressure
    // dT/dt = ...
    scalar rho = 0;
    for (label i = 0; i < this->c_.size(); i++)
    {
        const scalar W = this->specieThermo_[i].W();
        rho += W * this->c_[i];
    }

    scalar cp = 0;
    for (label i = 0; i < this->c_.size(); i++)
    {
        // cp function returns [J/(kmol K)]
        cp += this->c_[i] * this->specieThermo_[i].cp(p, T);
    }
    cp /= rho;

    // When mechanism reduction is active
    // dT is computed on the reduced set since dcdt is null
    // for species not involved in the simplified mechanism
    scalar dT = 0;
    for (label i = 0; i < this->nSpecie_; i++)
    {
        label si;
        if (reduced)
        {
            si = simplifiedToCompleteIndex_[i];
        }
        else
        {
            si = i;
        }

        // ha function returns [J/kmol]
        const scalar hi = this->specieThermo_[si].ha(p, T);
        dT += hi * dcdt[i];
    }
    dT /= rho * cp;

    dcdt[this->nSpecie_] = -dT;

    // dp/dt = ...
    dcdt[this->nSpecie_ + 1] = 0;
}

template <class ReactionThermo, class ThermoType>
void Foam::STDACChemistryModel<ReactionThermo, ThermoType>::jacobian(
    const scalar t,
    const scalarField &c,
    scalarField &dcdt,
    scalarSquareMatrix &J) const
{
    const bool reduced = mechRed_->active();

    // If the mechanism reduction is active, the computed Jacobian
    // is compact (size of the reduced set of species)
    // but according to the information of the complete set
    // (i.e. for the third-body efficiencies)

    const scalar T = c[this->nSpecie_];
    const scalar p = c[this->nSpecie_ + 1];

    if (reduced)
    {
        this->c_ = completeC_;
        for (label i = 0; i < NsDAC_; i++)
        {
            this->c_[simplifiedToCompleteIndex_[i]] = max(c[i], 0);
        }
    }
    else
    {
        forAll(this->c_, i)
        {
            this->c_[i] = max(c[i], 0);
        }
    }

    J = Zero;
    dcdt = Zero;
    scalarField hi(this->c_.size());
    scalarField cpi(this->c_.size());
    forAll(hi, i)
    {
        hi[i] = this->specieThermo_[i].ha(p, T);
        cpi[i] = this->specieThermo_[i].cp(p, T);
    }

    scalar omegaI = 0;

    forAll(this->reactions_, ri)
    {
        if (!reactionsDisabled_[ri])
        {
            const Reaction<ThermoType> &R = this->reactions_[ri];
            scalar kfwd, kbwd;
            R.dwdc(
                p,
                T,
                this->c_,
                J,
                dcdt,
                omegaI,
                kfwd,
                kbwd,
                reduced,
                completeToSimplifiedIndex_);
            R.dwdT(
                p,
                T,
                this->c_,
                omegaI,
                kfwd,
                kbwd,
                J,
                reduced,
                completeToSimplifiedIndex_,
                this->nSpecie_);
        }
    }

    // The species derivatives of the temperature term are partially computed
    // while computing dwdc, they are completed hereunder:
    scalar cpMean = 0;
    scalar dcpdTMean = 0;
    forAll(this->c_, i)
    {
        cpMean += this->c_[i] * cpi[i]; // J/(m3.K)
        // Already multiplied by rho
        dcpdTMean += this->c_[i] * this->specieThermo_[i].dcpdT(p, T);
    }

    scalar dTdt = 0;
    forAll(hi, i)
    {
        if (reduced)
        {
            const label si = completeToSimplifiedIndex_[i];
            if (si != -1)
            {
                dTdt += hi[i] * dcdt[si]; // J/(m3.s)
            }
        }
        else
        {
            dTdt += hi[i] * dcdt[i]; // J/(m3.s)
        }
    }
    dTdt /= -cpMean; // K/s
    dcdt[this->nSpecie_] = dTdt;

    for (label i = 0; i < this->nSpecie_; i++)
    {
        J(this->nSpecie_, i) = 0;
        for (label j = 0; j < this->nSpecie_; j++)
        {
            const label sj = reduced ? simplifiedToCompleteIndex_[j] : j;
            J(this->nSpecie_, i) += hi[sj] * J(j, i);
        }
        const label si = reduced ? simplifiedToCompleteIndex_[i] : i;
        J(this->nSpecie_, i) += cpi[si] * dTdt; // J/(mol.s)
        J(this->nSpecie_, i) /= -cpMean;        // K/s / (mol/m3)
    }

    // ddT of dTdt
    J(this->nSpecie_, this->nSpecie_) = 0;
    for (label i = 0; i < this->nSpecie_; i++)
    {
        const label si = reduced ? simplifiedToCompleteIndex_[i] : i;
        J(this->nSpecie_, this->nSpecie_) +=
            cpi[si] * dcdt[i] + hi[si] * J(i, this->nSpecie_);
    }
    J(this->nSpecie_, this->nSpecie_) += dTdt * dcpdTMean;
    J(this->nSpecie_, this->nSpecie_) /= -cpMean;
    J(this->nSpecie_, this->nSpecie_) += dTdt / T;
}
template <class ReactionThermo, class ThermoType>
Foam::scalar Foam::STDACChemistryModel<ReactionThermo, ThermoType>::solve(
    const scalar rho, scalar p, scalar T, scalarField &c, scalarField &c0, scalarField &phiq, scalarField &Rphiq, scalar &deltaTChem_, const scalar &deltaT, label &growOrAdd, label &celli, scalarField &RR)
{
    const bool reduced = mechRed_->active();
    scalar timeLeft = deltaT;
    scalar deltaTMin = great;
    // When tabulation is active (short-circuit evaluation for retrieve)
    // It first tries to retrieve the solution of the system with the
    // information stored through the tabulation method
    if (tabulation_->retrieve(phiq, Rphiq))
    {
        // Retrieved solution stored in Rphiq
        for (label i = 0; i < this->nSpecie(); i++)
        {
            c[i] = rho * Rphiq[i] / this->specieThermo_[i].W();
        }

        //searchISATCpuTime_ += clockTime_.timeIncrement();
    }
    // This position is reached when tabulation is not used OR
    // if the solution is not retrieved.
    // In the latter case, it adds the information to the tabulation
    // (it will either expand the current data or add a new stored point).
    else
    {
        // Store total time waiting to attribute to add or grow
        //scalar timeTmp = clockTime_.timeIncrement();

        if (reduced)
        {
            // Reduce mechanism change the number of species (only active)
            mechRed_->reduceMechanism(c, T, p);
            //nActiveSpecies += mechRed_->NsSimp();
            //nAvg++;
            //scalar timeIncr = clockTime_.timeIncrement();
            //reduceMechCpuTime_ += timeIncr;
            //timeTmp += timeIncr;
        }

        // Calculate the chemical source terms
        while (timeLeft > small)
        {
            scalar dt = timeLeft;
            if (reduced)
            {
                // completeC_ used in the overridden ODE methods
                // to update only the active species
                completeC_ = c;

                // Solve the reduced set of ODE
                this->solve(
                    simplifiedC_, T, p, dt, deltaTChem_);

                for (label i = 0; i < NsDAC_; i++)
                {
                    c[simplifiedToCompleteIndex_[i]] = simplifiedC_[i];
                }
            }
            else
            {
                this->solve(c, T, p, dt, deltaTChem_);
            }
            timeLeft -= dt;
        }

        {
            //scalar timeIncr = clockTime_.timeIncrement();
            //solveChemistryCpuTime_ += timeIncr;
            //timeTmp += timeIncr;
        }

        // If tabulation is used, we add the information computed here to
        // the stored points (either expand or add)

        forAll(c, i)
        {
            Rphiq[i] = c[i] / rho * this->specieThermo_[i].W();
        }
        if (tabulation_->variableTimeStep())
        {
            Rphiq[Rphiq.size() - 3] = T;
            Rphiq[Rphiq.size() - 2] = p;
            Rphiq[Rphiq.size() - 1] = deltaT;
        }
        else
        {
            Rphiq[Rphiq.size() - 2] = T;
            Rphiq[Rphiq.size() - 1] = p;
        }
        growOrAdd =
            tabulation_->add(phiq, Rphiq, rho, deltaT);
        //if (growOrAdd)
        //{
        //    this->setTabulationResultsAdd(celli);
        //addNewLeafCpuTime_ += clockTime_.timeIncrement() + timeTmp;
        // }
        //else
        //{
        //    this->setTabulationResultsGrow(celli);
        //growCpuTime_ += clockTime_.timeIncrement() + timeTmp;
        //}

        // When operations are done and if mechanism reduction is active,
        // the number of species (which also affects nEqns) is set back
        // to the total number of species (stored in the mechRed object)
        if (reduced)
        {
            this->nSpecie_ = mechRed_->nSpecie();
        }
        deltaTMin = deltaTChem_;
        //deltaTMin = min(deltaTChem_, deltaTMin);

        deltaTChem_ = min(deltaTChem_, this->deltaTChemMax_);
    }

    // Set the RR vector (used in the solver)
    for (label i = 0; i < this->nSpecie_; i++)
    {
        RR[i] =
            (c[i] - c0[i]) * this->specieThermo_[i].W() / deltaT;
    }
    return deltaTMin;
}

template <class ReactionThermo, class ThermoType>
Foam::scalar Foam::STDACChemistryModel<ReactionThermo, ThermoType>::solve(
    jobdata &in)
{
    scalar &rho = in.rho;
    scalar &p = in.p;
    scalar &T = in.T;
    //FixedList<scalar, DateSize1> &c = in.c;
    FixedList<scalar, DateSize1> &c0 = in.c0;
    //FixedList<scalar, DateSize2> &phiq = in.phiq;
    //FixedList<scalar, DateSize2> &Rphiq = in.Rphiq;

    scalarField phiq(in.phiq.size()), Rphiq(in.phiq.size());
    scalarField c(in.c.size());
    for (label i = 0; i < phiq.size(); i++)
    {
        phiq[i] = in.phiq[i];
        Rphiq[i] = in.Rphiq[i];
    }
    for (label i = 0; i < c.size(); i++)
    {
        c[i] = in.c[i];
    }

    scalar &deltaTChem_ = in.deltaTChem_;
    scalar &deltaT = in.deltaT;
    label &growOrAdd = in.growOrAdd;
    label &celli = in.celli;
    FixedList<scalar, DateSize1> &RR = in.RR;

    const bool reduced = mechRed_->active();
    scalar timeLeft = deltaT;
    scalar deltaTMin = great;
    // When tabulation is active (short-circuit evaluation for retrieve)
    // It first tries to retrieve the solution of the system with the
    // information stored through the tabulation method
    if (tabulation_->retrieve(phiq, Rphiq))
    {
        // Retrieved solution stored in Rphiq
        for (label i = 0; i < this->nSpecie(); i++)
        {
            c[i] = rho * Rphiq[i] / this->specieThermo_[i].W();
        }

        //searchISATCpuTime_ += clockTime_.timeIncrement();
    }
    // This position is reached when tabulation is not used OR
    // if the solution is not retrieved.
    // In the latter case, it adds the information to the tabulation
    // (it will either expand the current data or add a new stored point).
    else
    {
        // Store total time waiting to attribute to add or grow
        //scalar timeTmp = clockTime_.timeIncrement();

        if (reduced)
        {
            // Reduce mechanism change the number of species (only active)
            mechRed_->reduceMechanism(c, T, p);
            //nActiveSpecies += mechRed_->NsSimp();
            //nAvg++;
            //scalar timeIncr = clockTime_.timeIncrement();
            //reduceMechCpuTime_ += timeIncr;
            //timeTmp += timeIncr;
        }

        // Calculate the chemical source terms
        while (timeLeft > small)
        {
            scalar dt = timeLeft;
            if (reduced)
            {
                // completeC_ used in the overridden ODE methods
                // to update only the active species
                completeC_ = c;

                // Solve the reduced set of ODE
                this->solve(
                    simplifiedC_, T, p, dt, deltaTChem_);

                for (label i = 0; i < NsDAC_; i++)
                {
                    c[simplifiedToCompleteIndex_[i]] = simplifiedC_[i];
                }
            }
            else
            {
                this->solve(c, T, p, dt, deltaTChem_);
            }
            timeLeft -= dt;
        }

        {
            //scalar timeIncr = clockTime_.timeIncrement();
            //solveChemistryCpuTime_ += timeIncr;
            //timeTmp += timeIncr;
        }

        // If tabulation is used, we add the information computed here to
        // the stored points (either expand or add)

        forAll(c, i)
        {
            Rphiq[i] = c[i] / rho * this->specieThermo_[i].W();
        }
        if (tabulation_->variableTimeStep())
        {
            Rphiq[Rphiq.size() - 3] = T;
            Rphiq[Rphiq.size() - 2] = p;
            Rphiq[Rphiq.size() - 1] = deltaT;
        }
        else
        {
            Rphiq[Rphiq.size() - 2] = T;
            Rphiq[Rphiq.size() - 1] = p;
        }
        growOrAdd =
            tabulation_->add(phiq, Rphiq, rho, deltaT);
        //if (growOrAdd)
        //{
        //    this->setTabulationResultsAdd(celli);
        //addNewLeafCpuTime_ += clockTime_.timeIncrement() + timeTmp;
        // }
        //else
        //{
        //    this->setTabulationResultsGrow(celli);
        //growCpuTime_ += clockTime_.timeIncrement() + timeTmp;
        //}

        // When operations are done and if mechanism reduction is active,
        // the number of species (which also affects nEqns) is set back
        // to the total number of species (stored in the mechRed object)
        if (reduced)
        {
            this->nSpecie_ = mechRed_->nSpecie();
        }
        deltaTMin = deltaTChem_;
        //deltaTMin = min(deltaTChem_, deltaTMin);

        deltaTChem_ = min(deltaTChem_, this->deltaTChemMax_);
    }

    // Set the RR vector (used in the solver)
    for (label i = 0; i < this->nSpecie_; i++)
    {
        RR[i] =
            (c[i] - c0[i]) * this->specieThermo_[i].W() / deltaT;
    }
    return deltaTMin;
}

template <class ReactionThermo, class ThermoType>
Foam::scalar Foam::STDACChemistryModel<ReactionThermo, ThermoType>::solve(
    jobInput &in)
{
    scalar &rho = in.rho;
    scalar &p = in.p;
    scalar &T = in.T;
    //FixedList<scalar, DateSize1> &c = in.c;
    FixedList<scalar, DateSize1> &c0 = in.c0;
    //FixedList<scalar, DateSize2> &phiq = in.phiq;
    //FixedList<scalar, DateSize2> &Rphiq = in.Rphiq;

    scalarField phiq(in.phiq.size()), Rphiq(in.phiq.size());
    scalarField c(in.c.size());
    for (label i = 0; i < phiq.size(); i++)
    {
        phiq[i] = in.phiq[i];
        Rphiq[i] = in.Rphiq[i];
    }
    for (label i = 0; i < c.size(); i++)
    {
        c[i] = in.c[i];
    }

    scalar &deltaTChem_ = in.deltaTChem_;
    scalar &deltaT = in.deltaT;
    label &growOrAdd = in.growOrAdd;
    //label &celli = in.celli;
    FixedList<scalar, DateSize1> &RR = in.RR;

    const bool reduced = mechRed_->active();
    scalar timeLeft = deltaT;
    scalar deltaTMin = great;
    // When tabulation is active (short-circuit evaluation for retrieve)
    // It first tries to retrieve the solution of the system with the
    // information stored through the tabulation method
    if (tabulation_->retrieve(phiq, Rphiq))
    {
        // Retrieved solution stored in Rphiq
        for (label i = 0; i < this->nSpecie(); i++)
        {
            c[i] = rho * Rphiq[i] / this->specieThermo_[i].W();
        }

        //searchISATCpuTime_ += clockTime_.timeIncrement();
    }
    // This position is reached when tabulation is not used OR
    // if the solution is not retrieved.
    // In the latter case, it adds the information to the tabulation
    // (it will either expand the current data or add a new stored point).
    else
    {
        // Store total time waiting to attribute to add or grow
        //scalar timeTmp = clockTime_.timeIncrement();

        if (reduced)
        {
            // Reduce mechanism change the number of species (only active)
            mechRed_->reduceMechanism(c, T, p);
            //nActiveSpecies += mechRed_->NsSimp();
            //nAvg++;
            //scalar timeIncr = clockTime_.timeIncrement();
            //reduceMechCpuTime_ += timeIncr;
            //timeTmp += timeIncr;
        }

        // Calculate the chemical source terms
        while (timeLeft > small)
        {
            scalar dt = timeLeft;
            if (reduced)
            {
                // completeC_ used in the overridden ODE methods
                // to update only the active species
                completeC_ = c;

                // Solve the reduced set of ODE
                this->solve(
                    simplifiedC_, T, p, dt, deltaTChem_);

                for (label i = 0; i < NsDAC_; i++)
                {
                    c[simplifiedToCompleteIndex_[i]] = simplifiedC_[i];
                }
            }
            else
            {
                this->solve(c, T, p, dt, deltaTChem_);
            }
            timeLeft -= dt;
        }

        {
            //scalar timeIncr = clockTime_.timeIncrement();
            //solveChemistryCpuTime_ += timeIncr;
            //timeTmp += timeIncr;
        }

        // If tabulation is used, we add the information computed here to
        // the stored points (either expand or add)

        forAll(c, i)
        {
            Rphiq[i] = c[i] / rho * this->specieThermo_[i].W();
        }
        if (tabulation_->variableTimeStep())
        {
            Rphiq[Rphiq.size() - 3] = T;
            Rphiq[Rphiq.size() - 2] = p;
            Rphiq[Rphiq.size() - 1] = deltaT;
        }
        else
        {
            Rphiq[Rphiq.size() - 2] = T;
            Rphiq[Rphiq.size() - 1] = p;
        }
        //Info<< SUPstream::node_manager.rank << "," << phiq << endl;
        growOrAdd =
            tabulation_->add(phiq, Rphiq, rho, deltaT);
        //if (growOrAdd)
        //{
        //    this->setTabulationResultsAdd(celli);
        //addNewLeafCpuTime_ += clockTime_.timeIncrement() + timeTmp;
        // }
        //else
        //{
        //    this->setTabulationResultsGrow(celli);
        //growCpuTime_ += clockTime_.timeIncrement() + timeTmp;
        //}

        // When operations are done and if mechanism reduction is active,
        // the number of species (which also affects nEqns) is set back
        // to the total number of species (stored in the mechRed object)
        if (reduced)
        {
            this->nSpecie_ = mechRed_->nSpecie();
        }
        deltaTMin = deltaTChem_;
        //deltaTMin = min(deltaTChem_, deltaTMin);

        deltaTChem_ = min(deltaTChem_, this->deltaTChemMax_);
    }

    // Set the RR vector (used in the solver)
    for (label i = 0; i < this->nSpecie_; i++)
    {
        RR[i] =
            (c[i] - c0[i]) * this->specieThermo_[i].W() / deltaT;
    }
    return deltaTMin;
}

template <class ReactionThermo, class ThermoType>
template <class DeltaTType>
Foam::scalar Foam::STDACChemistryModel<ReactionThermo, ThermoType>::solve(
    const DeltaTType &deltaT)
{

    // Increment counter of time-step
    timeSteps_++;

    const bool reduced = mechRed_->active();

    label nAdditionalEqn = (tabulation_->variableTimeStep() ? 1 : 0);

    basicSpecieMixture &composition = this->thermo().composition();

    // CPU time analysis
    const clockTime clockTime_ = clockTime();
    //clockTime_.timeIncrement();
    scalar reduceMechCpuTime_ = 0;
    scalar addNewLeafCpuTime_ = 0;
    scalar growCpuTime_ = 0;
    scalar solveChemistryCpuTime_ = 0;
    scalar searchISATCpuTime_ = 0;

    this->resetTabulationResults();

    // Average number of active species
    scalar nActiveSpecies = 0;
    scalar nAvg = 0;

    BasicChemistryModel<ReactionThermo>::correct();

    scalar deltaTMin = great;

    if (!this->chemistry_)
    {
        return deltaTMin;
    }

    const volScalarField rho(
        IOobject(
            "rho",
            this->time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false),
        this->thermo().rho());

    const scalarField &T = this->thermo().T();
    const scalarField &p = this->thermo().p();

    scalarField c(this->nSpecie_);
    scalarField c0(this->nSpecie_);
    scalarField RR(this->nSpecie_);

    // Composition vector (Yi, T, p)
    scalarField phiq(this->nEqns() + nAdditionalEqn);

    scalarField Rphiq(this->nEqns() + nAdditionalEqn);
    SUPstream::Sync();
    scalar solveCpuTime_ = 0;
    clockTime_.timeIncrement();
    if (tabulation_->active())
    {
        if (loadBalance_)
        {

            spare = false;
            spare_cpu--;
            nfinished_block = 0;

            int batch_size = Batch_Size;
            int nloop = rho.size() / batch_size;
            int last_size = rho.size() % batch_size;
            if (last_size != 0)
            {
                nloop++;
            }
            else
            {
                last_size = batch_size;
            }

            for (int i = 0; i < nloop;)
            {
                int size_i = batch_size;
                if (i == nloop - 1)
                {
                    size_i = last_size;
                }
                int first_i = i * batch_size;

                ja_loc.first_i = first_i;
                ja_loc.size_i = size_i;
                for (int j = 0; j < ja_loc.size_i; j++)
                {

                    ja_loc.jobs[j].rho = rho[j + ja_loc.first_i];
                    ja_loc.jobs[j].p = p[j + ja_loc.first_i];
                    ja_loc.jobs[j].T = T[j + ja_loc.first_i];
                    ja_loc.jobs[j].deltaT = deltaT[j + ja_loc.first_i];

                    for (label k = 0; k < this->nSpecie_; k++)
                    {
                        ja_loc.jobs[j].c[k] = ja_loc.jobs[j].rho * this->Y_[k][j + ja_loc.first_i] / this->specieThermo_[k].W();
                        ja_loc.jobs[j].c0[k] = ja_loc.jobs[j].c[k];
                        ja_loc.jobs[j].phiq[k] = this->Y()[k][j + ja_loc.first_i];
                    }

                    ja_loc.jobs[j].phiq[this->nSpecie()] = ja_loc.jobs[j].T;
                    ja_loc.jobs[j].phiq[this->nSpecie() + 1] = ja_loc.jobs[j].p;
                    if (tabulation_->variableTimeStep())
                    {
                        ja_loc.jobs[j].phiq[this->nSpecie() + 2] = deltaT[j + ja_loc.first_i];
                    }

                    // Initialise time progress

                    // Not sure if this is necessary
                    ja_loc.jobs[j].Rphiq = Zero;

                    //clockTime_.timeIncrement();

                    ja_loc.jobs[j].deltaTChem_ = this->deltaTChem_[j + ja_loc.first_i];
                }

                for (int j = 0; j < ja_loc.size_i; j++)
                {
                    deltaTMin = min(solve(ja_loc.jobs[j]), deltaTMin);
                    this->deltaTChem_[j + ja_loc.first_i] = ja_loc.jobs[j].deltaTChem_;
                    if (ja_loc.jobs[j].growOrAdd)
                    {
                        this->setTabulationResultsAdd(j + ja_loc.first_i);
                    }
                    else
                    {
                        this->setTabulationResultsGrow(j + ja_loc.first_i);
                    }

                    // Set the RR vector (used in the solver)
                    for (label k = 0; k < this->nSpecie_; k++)
                    {
                        this->RR_[k][j + ja_loc.first_i] = ja_loc.jobs[j].RR[k];
                    }
                }

                nfinished_block++;
                i++;

                if (finished_head[rank] != -1)
                {
                    int iter_start;
                    int iter;
                    int tail;

                    l_finished.lock();
                    iter_start = finished_head[rank];
                    iter = iter_start;
                    finished_head[rank] = -1;
                    l_finished.unlock();

                    tail = iter;
                    while (iter != -1)
                    {

                        for (int j = 0; j < ja[iter].size_i; j++)
                        {
                            this->deltaTChem_[j + ja[iter].first_i] = ja[iter].jobs[j].deltaTChem_;
                            if (ja[iter].jobs[j].growOrAdd)
                            {
                                this->setTabulationResultsAdd(j + ja[iter].first_i);
                            }
                            else
                            {
                                this->setTabulationResultsGrow(j + ja[iter].first_i);
                            }

                            // Set the RR vector (used in the solver)
                            for (label k = 0; k < this->nSpecie_; k++)
                            {
                                this->RR_[k][j + ja[iter].first_i] = ja[iter].jobs[j].RR[k];
                            }
                        }
                        tail = iter;
                        iter = ja[iter].next;
                        nfinished_block++;
                    }

                    //ja[0].finished = false;
                    //ja[0].filled = false;

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

                if (i < nloop && spare_cpu > 0 && l_sender.try_lock())
                {
                    for (int ii = 0; ii <= spare_cpu; ii++)
                    {
                        int size_i = batch_size;
                        if (i == nloop - 1)
                        {
                            size_i = last_size;
                        }
                        int first_i = i * batch_size;

                        if (empty_head != -1 && i < nloop)
                        {
                            int iter;

                            l_empty.lock();
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
                                //ja[iter].finished = false;
                                ja[iter].jobs[j].rho = rho[j + ja[iter].first_i];
                                ja[iter].jobs[j].p = p[j + ja[iter].first_i];
                                ja[iter].jobs[j].T = T[j + ja[iter].first_i];
                                ja[iter].jobs[j].deltaT = deltaT[j + ja[iter].first_i];

                                for (label k = 0; k < this->nSpecie_; k++)
                                {
                                    ja[iter].jobs[j].c[k] = ja[iter].jobs[j].rho * this->Y_[k][j + ja[iter].first_i] / this->specieThermo_[k].W();
                                    ja[iter].jobs[j].c0[k] = ja[iter].jobs[j].c[k];
                                    ja[iter].jobs[j].phiq[k] = this->Y()[k][j + ja[iter].first_i];
                                }

                                ja[iter].jobs[j].phiq[this->nSpecie()] = ja[iter].jobs[j].T;
                                ja[iter].jobs[j].phiq[this->nSpecie() + 1] = ja[iter].jobs[j].p;
                                if (tabulation_->variableTimeStep())
                                {
                                    ja[iter].jobs[j].phiq[this->nSpecie() + 2] = deltaT[j + ja[iter].first_i];
                                }

                                // Initialise time progress

                                // Not sure if this is necessary
                                ja[iter].jobs[j].Rphiq = Zero;

                                //clockTime_.timeIncrement();

                                //label growOrAdd;
                                ja[iter].jobs[j].deltaTChem_ = this->deltaTChem_[j + ja[iter].first_i];
                            }

                            ja[iter].rank = rank;
                            //ja[0].finished = false;
                            //ja[iter].filled = true;

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
                            i++;
                        }
                        else
                            break;
                        l_sender.unlock();
                    }
                }
            }
            //int xx = 0;

            //std::cout << SUPstream::node_manager.rank << "!!!!!!!!!!!!!!!!!!!0 " << std::endl;
            while (1)
            {
                //l_receiver.lock();
                if (l_receiver.try_lock())
                {
                    int iter = -1;
                    if (filled_head != -1)
                    {

                        l_filled.lock();
                        iter = filled_head;
                        //std::cout << SUPstream::node_manager.rank << "," << filled_head << "," << ja[iter].first_i << std::endl;
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
                        for (int j = 0; j < ja[iter].size_i; j++)
                        {
                            deltaTMin = min(solve(ja[iter].jobs[j]), deltaTMin);
                        }

                        l_finished.lock();
                        ja[iter].next = finished_head[ja[iter].rank];
                        finished_head[ja[iter].rank] = iter;
                        //ja[iter].finished = true;
                        l_finished.unlock();
                    }

                    //
                }
                if (finished_head[rank] != -1)
                {
                    int iter_start;
                    int iter;
                    int tail;

                    l_finished.lock();
                    iter_start = finished_head[rank];
                    iter = iter_start;
                    finished_head[rank] = -1;
                    l_finished.unlock();

                    tail = iter;
                    while (iter != -1)
                    {

                        for (int j = 0; j < ja[iter].size_i; j++)
                        {
                            this->deltaTChem_[j + ja[iter].first_i] = ja[iter].jobs[j].deltaTChem_;
                            if (ja[iter].jobs[j].growOrAdd)
                            {
                                this->setTabulationResultsAdd(j + ja[iter].first_i);
                            }
                            else
                            {
                                this->setTabulationResultsGrow(j + ja[iter].first_i);
                            }

                            //                     this->deltaTChem_[celli] =
                            //min(this->deltaTChem_[celli], this->deltaTChemMax_);

                            // Set the RR vector (used in the solver)
                            for (label k = 0; k < this->nSpecie_; k++)
                            {
                                this->RR_[k][j + ja[iter].first_i] = ja[iter].jobs[j].RR[k];
                            }
                        }
                        tail = iter;
                        iter = ja[iter].next;
                        nfinished_block++;
                    }

                    //ja[0].finished = false;
                    //ja[0].filled = false;

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
                //std::cout << SUPstream::node_manager.rank << "!!!!!!!!!!!!!!!!!!!nfinished_block= "<< nfinished_block << std::endl;
                if (spare_cpu == n_cpu && nfinished_block == nloop)
                {
                    break;
                }

                //l_sender.unlock();
            }
            //std::cout << SUPstream::node_manager.rank << "!!!!!!!!!!!!!!!!!!!1 " << std::endl;
            SUPstream::Sync();
        }
        //l_receiver.lock();

        else
        {

            int batch_size = Batch_Size;
            int nloop = rho.size() / batch_size;
            int last_size = rho.size() % batch_size;
            if (last_size != 0)
            {
                nloop++;
            }
            else
            {
                last_size = batch_size;
            }
            for (int i = 0; i < nloop; i++)
            {
                int size_i = batch_size;
                //int last_i = (i + 1) * batch_size;
                if (i == nloop - 1)
                {
                    //last_i = i * batch_size + last_size;
                    size_i = last_size;
                }
                int first_i = i * batch_size;

                ja_loc.first_i = first_i;
                ja_loc.size_i = size_i;
                //ja.last_i = last_i;

                for (int j = 0; j < ja_loc.size_i; j++)
                {

                    ja_loc.jobs[j].rho = rho[j + ja_loc.first_i];
                    ja_loc.jobs[j].p = p[j + ja_loc.first_i];
                    ja_loc.jobs[j].T = T[j + ja_loc.first_i];
                    ja_loc.jobs[j].deltaT = deltaT[j + ja_loc.first_i];

                    //tmp.celli = celli;
                    //tmp.rho = rho[celli];
                    //tmp.p = p[celli];
                    //tmp.T = T[celli];
                    //tmp.deltaT = deltaT[celli];

                    //const scalar rhoi = rho[celli];
                    //scalar pi = p[celli];
                    //scalar Ti = T[celli];

                    for (label k = 0; k < this->nSpecie_; k++)
                    {
                        //c[i] = rhoi * this->Y_[i][celli] / this->specieThermo_[i].W();

                        //tmp.c[i] = tmp.rho * this->Y_[i][celli] / this->specieThermo_[i].W();
                        //tmp.c0[i] = tmp.c[i];
                        //tmp.phiq[i] = this->Y()[i][celli];

                        ja_loc.jobs[j].c[k] = ja_loc.jobs[j].rho * this->Y_[k][j + ja_loc.first_i] / this->specieThermo_[k].W();
                        ja_loc.jobs[j].c0[k] = ja_loc.jobs[j].c[k];
                        ja_loc.jobs[j].phiq[k] = this->Y()[k][j + ja_loc.first_i];
                    }
                    //phiq[this->nSpecie()] = Ti;
                    //phiq[this->nSpecie() + 1] = pi;

                    //tmp.phiq[this->nSpecie()] = tmp.T;
                    //tmp.phiq[this->nSpecie() + 1] = tmp.p;
                    //if (tabulation_->variableTimeStep())
                    //{
                    //    tmp.phiq[this->nSpecie() + 2] = deltaT[celli];
                    //}

                    ja_loc.jobs[j].phiq[this->nSpecie()] = ja_loc.jobs[j].T;
                    ja_loc.jobs[j].phiq[this->nSpecie() + 1] = ja_loc.jobs[j].p;
                    if (tabulation_->variableTimeStep())
                    {
                        ja_loc.jobs[j].phiq[this->nSpecie() + 2] = deltaT[j + ja_loc.first_i];
                    }

                    // Initialise time progress
                    //scalar timeLeft = deltaT[j + ja.first_i];

                    // Not sure if this is necessary
                    //tmp.Rphiq = Zero;
                    ja_loc.jobs[j].Rphiq = Zero;

                    //clockTime_.timeIncrement();

                    //label growOrAdd;
                    ja_loc.jobs[j].deltaTChem_ = this->deltaTChem_[j + ja_loc.first_i];
                }

                for (int j = 0; j < ja_loc.size_i; j++)
                {
                    //deltaTMin = min(solve(rhoi, pi, Ti, c, c0, phiq, Rphiq, this->deltaTChem_[celli], deltaT[celli], growOrAdd, celli, RR), deltaTMin);
                    deltaTMin = min(solve(ja_loc.jobs[j]), deltaTMin);
                    this->deltaTChem_[j + ja_loc.first_i] = ja_loc.jobs[j].deltaTChem_;
                    if (ja_loc.jobs[j].growOrAdd)
                    {
                        this->setTabulationResultsAdd(j + ja_loc.first_i);
                    }
                    else
                    {
                        this->setTabulationResultsGrow(j + ja_loc.first_i);
                    }

                    /*                     this->deltaTChem_[celli] =
                        min(this->deltaTChem_[celli], this->deltaTChemMax_); */

                    // Set the RR vector (used in the solver)
                    for (label k = 0; k < this->nSpecie_; k++)
                    {
                        this->RR_[k][j + ja_loc.first_i] = ja_loc.jobs[j].RR[k];
                    }
                }
            }
        }
    }
    else
    {
        //std::cout << "xzxhi!!!!!!!!!!!!!!!!!!!!!!!!!!!" << rho.size() << std::endl;
        forAll(rho, celli)
        {
            const scalar rhoi = rho[celli];
            scalar pi = p[celli];
            scalar Ti = T[celli];

            for (label i = 0; i < this->nSpecie_; i++)
            {
                c[i] = rhoi * this->Y_[i][celli] / this->specieThermo_[i].W();
                c0[i] = c[i];
                phiq[i] = this->Y()[i][celli];
            }
            phiq[this->nSpecie()] = Ti;
            phiq[this->nSpecie() + 1] = pi;
            if (tabulation_->variableTimeStep())
            {
                phiq[this->nSpecie() + 2] = deltaT[celli];
            }

            // Initialise time progress
            scalar timeLeft = deltaT[celli];

            // Not sure if this is necessary
            Rphiq = Zero;

            clockTime_.timeIncrement();

            // When tabulation is active (short-circuit evaluation for retrieve)
            // It first tries to retrieve the solution of the system with the
            // information stored through the tabulation method
            if (tabulation_->active() && tabulation_->retrieve(phiq, Rphiq))
            {
                // Retrieved solution stored in Rphiq
                for (label i = 0; i < this->nSpecie(); i++)
                {
                    c[i] = rhoi * Rphiq[i] / this->specieThermo_[i].W();
                }

                searchISATCpuTime_ += clockTime_.timeIncrement();
            }
            // This position is reached when tabulation is not used OR
            // if the solution is not retrieved.
            // In the latter case, it adds the information to the tabulation
            // (it will either expand the current data or add a new stored point).
            else
            {
                // Store total time waiting to attribute to add or grow
                scalar timeTmp = clockTime_.timeIncrement();

                if (reduced)
                {
                    // Reduce mechanism change the number of species (only active)
                    mechRed_->reduceMechanism(c, Ti, pi);
                    nActiveSpecies += mechRed_->NsSimp();
                    nAvg++;
                    scalar timeIncr = clockTime_.timeIncrement();
                    reduceMechCpuTime_ += timeIncr;
                    timeTmp += timeIncr;
                }

                // Calculate the chemical source terms
                while (timeLeft > small)
                {
                    scalar dt = timeLeft;
                    if (reduced)
                    {
                        // completeC_ used in the overridden ODE methods
                        // to update only the active species
                        completeC_ = c;

                        // Solve the reduced set of ODE
                        this->solve(
                            simplifiedC_, Ti, pi, dt, this->deltaTChem_[celli]);

                        for (label i = 0; i < NsDAC_; i++)
                        {
                            c[simplifiedToCompleteIndex_[i]] = simplifiedC_[i];
                        }
                    }
                    else
                    {
                        this->solve(c, Ti, pi, dt, this->deltaTChem_[celli]);
                    }
                    timeLeft -= dt;
                }

                {
                    scalar timeIncr = clockTime_.timeIncrement();
                    solveChemistryCpuTime_ += timeIncr;
                    timeTmp += timeIncr;
                }

                // If tabulation is used, we add the information computed here to
                // the stored points (either expand or add)
                if (tabulation_->active())
                {
                    forAll(c, i)
                    {
                        Rphiq[i] = c[i] / rhoi * this->specieThermo_[i].W();
                    }
                    if (tabulation_->variableTimeStep())
                    {
                        Rphiq[Rphiq.size() - 3] = Ti;
                        Rphiq[Rphiq.size() - 2] = pi;
                        Rphiq[Rphiq.size() - 1] = deltaT[celli];
                    }
                    else
                    {
                        Rphiq[Rphiq.size() - 2] = Ti;
                        Rphiq[Rphiq.size() - 1] = pi;
                    }
                    label growOrAdd =
                        tabulation_->add(phiq, Rphiq, rhoi, deltaT[celli]);
                    if (growOrAdd)
                    {
                        this->setTabulationResultsAdd(celli);
                        addNewLeafCpuTime_ += clockTime_.timeIncrement() + timeTmp;
                    }
                    else
                    {
                        this->setTabulationResultsGrow(celli);
                        growCpuTime_ += clockTime_.timeIncrement() + timeTmp;
                    }
                }

                // When operations are done and if mechanism reduction is active,
                // the number of species (which also affects nEqns) is set back
                // to the total number of species (stored in the mechRed object)
                if (reduced)
                {
                    this->nSpecie_ = mechRed_->nSpecie();
                }
                deltaTMin = min(this->deltaTChem_[celli], deltaTMin);

                this->deltaTChem_[celli] =
                    min(this->deltaTChem_[celli], this->deltaTChemMax_);
            }

            // Set the RR vector (used in the solver)
            for (label i = 0; i < this->nSpecie_; i++)
            {
                this->RR_[i][celli] =
                    (c[i] - c0[i]) * this->specieThermo_[i].W() / deltaT[celli];
            }
        }
    }
    solveCpuTime_ += clockTime_.timeIncrement();

    if (mechRed_->log() || tabulation_->log())
    {
        cpuSolveFile_()
            << this->time().timeOutputValue()
            << "    " << solveChemistryCpuTime_ << endl;
    }

    if (mechRed_->log())
    {
        cpuReduceFile_()
            << this->time().timeOutputValue()
            << "    " << reduceMechCpuTime_ << endl;
    }

    if (tabulation_->active())
    {
        // Every time-step, look if the tabulation should be updated
        tabulation_->update();

        // Write the performance of the tabulation
        tabulation_->writePerformance();

        if (tabulation_->log())
        {
            cpuRetrieveFile_()
                << this->time().timeOutputValue()
                << "    " << searchISATCpuTime_ << endl;

            cpuGrowFile_()
                << this->time().timeOutputValue()
                << "    " << growCpuTime_ << endl;

            cpuAddFile_()
                << this->time().timeOutputValue()
                << "    " << addNewLeafCpuTime_ << endl;

            cpuTotalFile_()
                << this->time().timeOutputValue()
                << "    " << solveCpuTime_ << endl;
        }
    }

    if (reduced && nAvg && mechRed_->log())
    {
        // Write average number of species
        nActiveSpeciesFile_()
            << this->time().timeOutputValue()
            << "    " << nActiveSpecies / nAvg << endl;
    }

    if (Pstream::parRun())
    {
        List<bool> active(composition.active());
        Pstream::listCombineGather(active, orEqOp<bool>());
        Pstream::listCombineScatter(active);

        forAll(active, i)
        {
            if (active[i])
            {
                composition.setActive(i);
            }
        }
    }

    forAll(this->Y(), i)
    {
        if (composition.active(i))
        {
            this->Y()[i].writeOpt() = IOobject::AUTO_WRITE;
        }
    }

    return deltaTMin;
}

template <class ReactionThermo, class ThermoType>
template <class DeltaTType>
Foam::scalar Foam::STDACChemistryModel<ReactionThermo, ThermoType>::solve2(
    const DeltaTType &deltaT)
{
    //std::cout << "xxxxxxzzhi!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
    // Increment counter of time-step
    timeSteps_++;

    const bool reduced = mechRed_->active();

    label nAdditionalEqn = (tabulation_->variableTimeStep() ? 1 : 0);

    basicSpecieMixture &composition = this->thermo().composition();

    // CPU time analysis
    const clockTime clockTime_ = clockTime();
    clockTime_.timeIncrement();
    scalar reduceMechCpuTime_ = 0;
    scalar addNewLeafCpuTime_ = 0;
    scalar growCpuTime_ = 0;
    scalar solveChemistryCpuTime_ = 0;
    scalar searchISATCpuTime_ = 0;

    this->resetTabulationResults();

    // Average number of active species
    scalar nActiveSpecies = 0;
    scalar nAvg = 0;

    BasicChemistryModel<ReactionThermo>::correct();

    scalar deltaTMin = great;

    if (!this->chemistry_)
    {
        return deltaTMin;
    }

    const volScalarField rho(
        IOobject(
            "rho",
            this->time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false),
        this->thermo().rho());

    const scalarField &T = this->thermo().T();
    const scalarField &p = this->thermo().p();

    scalarField c(this->nSpecie_);
    scalarField c0(this->nSpecie_);
    scalarField RR(this->nSpecie_);

    //jobdata tmp;
    //tmp.c.setSize(this->nSpecie_);
    //tmp.c0.setSize(this->nSpecie_);
    //tmp.RR.setSize(this->nSpecie_);

    // Composition vector (Yi, T, p)
    scalarField phiq(this->nEqns() + nAdditionalEqn);

    scalarField Rphiq(this->nEqns() + nAdditionalEqn);

    //tmp.phiq.setSize(this->nEqns() + nAdditionalEqn);
    //tmp.Rphiq.setSize(this->nEqns() + nAdditionalEqn);
    /*     FatalErrorInFunction << "nSpecie_=" << this->nSpecie_
                         << " ,this->nEqns() + nAdditionalEqn= " << this->nEqns() + nAdditionalEqn
                         << abort(FatalError); */
    //std::cout << "xxxxxxhi!!!!!!!!!!!!!!!!!!!!!!!!!!!" << rho.size() << std::endl;
    if (tabulation_->active())
    {
        if (SUPstream::node_manager.rank != 2)
        {
            //spare_cpu = 0;

            //std::cout << "xxhi!!!!!!!!!!!!!!!!!!!!!!!!!!!" << rho.size() << std::endl;

            int batch_size = Batch_Size;
            int nloop = rho.size() / batch_size;
            int last_size = rho.size() % batch_size;
            if (last_size != 0)
            {
                nloop++;
            }
            else
            {
                last_size = batch_size;
            }
            for (int i = 0; i < nloop; i++)
            {
                int size_i = batch_size;
                //int last_i = (i + 1) * batch_size;
                if (i == nloop - 1)
                {
                    //last_i = i * batch_size + last_size;
                    size_i = last_size;
                }
                int first_i = i * batch_size;

                ja_loc.first_i = first_i;
                ja_loc.size_i = size_i;
                //ja.last_i = last_i;

                for (int j = 0; j < ja_loc.size_i; j++)
                {

                    ja_loc.jobs[j].rho = rho[j + ja_loc.first_i];
                    ja_loc.jobs[j].p = p[j + ja_loc.first_i];
                    ja_loc.jobs[j].T = T[j + ja_loc.first_i];
                    ja_loc.jobs[j].deltaT = deltaT[j + ja_loc.first_i];

                    //tmp.celli = celli;
                    //tmp.rho = rho[celli];
                    //tmp.p = p[celli];
                    //tmp.T = T[celli];
                    //tmp.deltaT = deltaT[celli];

                    //const scalar rhoi = rho[celli];
                    //scalar pi = p[celli];
                    //scalar Ti = T[celli];

                    for (label k = 0; k < this->nSpecie_; k++)
                    {
                        //c[i] = rhoi * this->Y_[i][celli] / this->specieThermo_[i].W();

                        //tmp.c[i] = tmp.rho * this->Y_[i][celli] / this->specieThermo_[i].W();
                        //tmp.c0[i] = tmp.c[i];
                        //tmp.phiq[i] = this->Y()[i][celli];

                        ja_loc.jobs[j].c[k] = ja_loc.jobs[j].rho * this->Y_[k][j + ja_loc.first_i] / this->specieThermo_[k].W();
                        ja_loc.jobs[j].c0[k] = ja_loc.jobs[j].c[k];
                        ja_loc.jobs[j].phiq[k] = this->Y()[k][j + ja_loc.first_i];
                    }
                    //phiq[this->nSpecie()] = Ti;
                    //phiq[this->nSpecie() + 1] = pi;

                    //tmp.phiq[this->nSpecie()] = tmp.T;
                    //tmp.phiq[this->nSpecie() + 1] = tmp.p;
                    //if (tabulation_->variableTimeStep())
                    //{
                    //    tmp.phiq[this->nSpecie() + 2] = deltaT[celli];
                    //}

                    ja_loc.jobs[j].phiq[this->nSpecie()] = ja_loc.jobs[j].T;
                    ja_loc.jobs[j].phiq[this->nSpecie() + 1] = ja_loc.jobs[j].p;
                    if (tabulation_->variableTimeStep())
                    {
                        ja_loc.jobs[j].phiq[this->nSpecie() + 2] = deltaT[j + ja_loc.first_i];
                    }

                    // Initialise time progress
                    //scalar timeLeft = deltaT[j + ja.first_i];

                    // Not sure if this is necessary
                    //tmp.Rphiq = Zero;
                    ja_loc.jobs[j].Rphiq = Zero;

                    clockTime_.timeIncrement();

                    //label growOrAdd;
                    ja_loc.jobs[j].deltaTChem_ = this->deltaTChem_[j + ja_loc.first_i];
                }

                for (int j = 0; j < ja_loc.size_i; j++)
                {
                    //deltaTMin = min(solve(rhoi, pi, Ti, c, c0, phiq, Rphiq, this->deltaTChem_[celli], deltaT[celli], growOrAdd, celli, RR), deltaTMin);
                    deltaTMin = min(solve(ja_loc.jobs[j]), deltaTMin);
                    this->deltaTChem_[j + ja_loc.first_i] = ja_loc.jobs[j].deltaTChem_;
                    if (ja_loc.jobs[j].growOrAdd)
                    {
                        this->setTabulationResultsAdd(j + ja_loc.first_i);
                    }
                    else
                    {
                        this->setTabulationResultsGrow(j + ja_loc.first_i);
                    }

                    /*                     this->deltaTChem_[celli] =
                        min(this->deltaTChem_[celli], this->deltaTChemMax_); */

                    // Set the RR vector (used in the solver)
                    for (label k = 0; k < this->nSpecie_; k++)
                    {
                        this->RR_[k][j + ja_loc.first_i] = ja_loc.jobs[j].RR[k];
                    }
                }
            }
            int xx = 0;
            spare = false;
            spare_cpu--;
            SUPstream::Sync();
            //std::cout << SUPstream::node_manager.rank << "!!!!!!!!!!!!!!!!!!!0 " << std::endl;
            while (1)
            {
                //l_receiver.lock();
                if (l_receiver.try_lock())
                {
                    int iter = -1;
                    if (filled_head != -1)
                    {

                        l_filled.lock();
                        iter = filled_head;
                        //std::cout << SUPstream::node_manager.rank << "," << filled_head << "," << ja[iter].first_i << std::endl;
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
                        //std::cout << SUPstream::node_manager.rank << ", working for" << ja[iter].rank << "," << ja[iter].first_i << std::endl;
                        if (spare)
                        {
                            spare = false;
                            spare_cpu--;
                        }
                        xx++;
                        for (int j = 0; j < ja[iter].size_i; j++)
                        {
                            //deltaTMin = min(solve(rhoi, pi, Ti, c, c0, phiq, Rphiq, this->deltaTChem_[celli], deltaT[celli], growOrAdd, celli, RR), deltaTMin);
                            deltaTMin = min(solve(ja[iter].jobs[j]), deltaTMin);
                        }

                        //std::cout << SUPstream::node_manager.rank << ", work finished," << ja[iter].rank << "," << ja[iter].first_i << std::endl;

                        l_finished.lock();
                        ja[iter].next = finished_head[ja[iter].rank];
                        finished_head[ja[iter].rank] = iter;
                        ja[iter].finished = true;
                        l_finished.unlock();
                    }

                    //
                }
                //std::cout << SUPstream::node_manager.rank << "," << spare_cpu << "," << n_cpu << ",spare " << spare << std::endl;
                if (!spare)
                {
                    spare = true;
                    spare_cpu++;
                }
                //std::cout << spare_cpu << std::endl;
                if (spare_cpu == n_cpu)
                {
                    //std::cout << SUPstream::node_manager.rank << "," << nloop << "," << xx << std::endl;
                    break;
                }

                //l_sender.unlock();
            }
            SUPstream::Sync();
        }
        //l_receiver.lock();

        else
        {

            spare = false;
            spare_cpu--;
            SUPstream::Sync();
            //std::cout << SUPstream::node_manager.rank << "!!!!!!!!!!!!!!!!!!!1 " << std::endl;
            int batch_size = Batch_Size;
            int nloop = rho.size() / batch_size;
            int last_size = rho.size() % batch_size;
            if (last_size != 0)
            {
                nloop++;
            }
            else
            {
                last_size = batch_size;
            }

            int nfinished_block = 0;
            for (int i = 0;;)
            {
                if (l_sender.try_lock())
                {
                    int size_i = batch_size;
                    //int last_i = (i + 1) * batch_size;
                    if (i == nloop - 1)
                    {
                        //last_i = i * batch_size + last_size;
                        size_i = last_size;
                    }
                    int first_i = i * batch_size;

                    if (finished_head[rank] != -1)
                    {
                        //std::cout << "hi!!!!!!!!!!!!!!!!!!!!!!!!!!!" << rho.size() << std::endl;

                        int iter_start;
                        int iter;
                        int tail;

                        l_finished.lock();
                        iter_start = finished_head[rank];
                        iter = iter_start;
                        finished_head[rank] = -1;
                        l_finished.unlock();

                        tail = iter;
                        //std::cout << SUPstream::node_manager.rank << ",finished " << iter_start << "," << ja[iter].first_i << std::endl;
                        while (iter != -1)
                        {

                            for (int j = 0; j < ja[iter].size_i; j++)
                            {
                                //deltaTMin = min(solve(rhoi, pi, Ti, c, c0, phiq, Rphiq, this->deltaTChem_[celli], deltaT[celli], growOrAdd, celli, RR), deltaTMin);
                                //deltaTMin = min(solve(ja.jobs[j]), deltaTMin);
                                this->deltaTChem_[j + ja[iter].first_i] = ja[iter].jobs[j].deltaTChem_;
                                if (ja[iter].jobs[j].growOrAdd)
                                {
                                    this->setTabulationResultsAdd(j + ja[iter].first_i);
                                }
                                else
                                {
                                    this->setTabulationResultsGrow(j + ja[iter].first_i);
                                }

                                /*                     this->deltaTChem_[celli] =
                        min(this->deltaTChem_[celli], this->deltaTChemMax_); */

                                // Set the RR vector (used in the solver)
                                for (label k = 0; k < this->nSpecie_; k++)
                                {
                                    this->RR_[k][j + ja[iter].first_i] = ja[iter].jobs[j].RR[k];
                                }
                            }
                            tail = iter;
                            iter = ja[iter].next;
                            nfinished_block++;
                        }

                        ja[0].finished = false;
                        ja[0].filled = false;

                        l_empty.lock();

                        //std::cout << SUPstream::node_manager.rank << "," << empty_tail << std::endl;
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
                        //????
                    }

                    if (empty_head != -1 && i < nloop)
                    {
                        //std::cout << SUPstream::node_manager.rank << "," << empty_head << "," << first_i << "," << rho.size() << std::endl;
                        int iter;

                        l_empty.lock();
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
                            ja[iter].finished = false;
                            ja[iter].jobs[j].rho = rho[j + ja[iter].first_i];
                            ja[iter].jobs[j].p = p[j + ja[iter].first_i];
                            ja[iter].jobs[j].T = T[j + ja[iter].first_i];
                            ja[iter].jobs[j].deltaT = deltaT[j + ja[iter].first_i];

                            //tmp.celli = celli;
                            //tmp.rho = rho[celli];
                            //tmp.p = p[celli];
                            //tmp.T = T[celli];
                            //tmp.deltaT = deltaT[celli];

                            //const scalar rhoi = rho[celli];
                            //scalar pi = p[celli];
                            //scalar Ti = T[celli];

                            for (label k = 0; k < this->nSpecie_; k++)
                            {
                                //c[i] = rhoi * this->Y_[i][celli] / this->specieThermo_[i].W();

                                //tmp.c[i] = tmp.rho * this->Y_[i][celli] / this->specieThermo_[i].W();
                                //tmp.c0[i] = tmp.c[i];
                                //tmp.phiq[i] = this->Y()[i][celli];

                                ja[iter].jobs[j].c[k] = ja[iter].jobs[j].rho * this->Y_[k][j + ja[iter].first_i] / this->specieThermo_[k].W();
                                ja[iter].jobs[j].c0[k] = ja[iter].jobs[j].c[k];
                                ja[iter].jobs[j].phiq[k] = this->Y()[k][j + ja[iter].first_i];
                            }
                            //phiq[this->nSpecie()] = Ti;
                            //phiq[this->nSpecie() + 1] = pi;

                            //tmp.phiq[this->nSpecie()] = tmp.T;
                            //tmp.phiq[this->nSpecie() + 1] = tmp.p;
                            //if (tabulation_->variableTimeStep())
                            //{
                            //    tmp.phiq[this->nSpecie() + 2] = deltaT[celli];
                            //}

                            ja[iter].jobs[j].phiq[this->nSpecie()] = ja[iter].jobs[j].T;
                            ja[iter].jobs[j].phiq[this->nSpecie() + 1] = ja[iter].jobs[j].p;
                            if (tabulation_->variableTimeStep())
                            {
                                ja[iter].jobs[j].phiq[this->nSpecie() + 2] = deltaT[j + ja[iter].first_i];
                            }

                            // Initialise time progress
                            //scalar timeLeft = deltaT[j + ja.first_i];

                            // Not sure if this is necessary
                            //tmp.Rphiq = Zero;
                            ja[iter].jobs[j].Rphiq = Zero;

                            clockTime_.timeIncrement();

                            //label growOrAdd;
                            ja[iter].jobs[j].deltaTChem_ = this->deltaTChem_[j + ja[iter].first_i];
                        }

                        ja[iter].rank = rank;
                        //ja[0].finished = false;
                        ja[iter].filled = true;

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
                        i++;
                    }
                    l_sender.unlock();
                }

                if (nfinished_block == nloop)
                {
                    break;
                }
            }
            while (1)
            {
                if (!spare)
                {
                    spare = true;
                    spare_cpu++;
                    //std::cout << << spare_cpu << std::endl;
                }
                //spare_cpu << std::endl;
                //std::cout << SUPstream::node_manager.rank << "," << spare_cpu << "," << n_cpu << std::endl;
                if (spare_cpu == n_cpu)
                {
                    break;
                }
            }
            SUPstream::Sync();
        }
    }
    else
    {
        //std::cout << "xzxhi!!!!!!!!!!!!!!!!!!!!!!!!!!!" << rho.size() << std::endl;
        forAll(rho, celli)
        {
            const scalar rhoi = rho[celli];
            scalar pi = p[celli];
            scalar Ti = T[celli];

            for (label i = 0; i < this->nSpecie_; i++)
            {
                c[i] = rhoi * this->Y_[i][celli] / this->specieThermo_[i].W();
                c0[i] = c[i];
                phiq[i] = this->Y()[i][celli];
            }
            phiq[this->nSpecie()] = Ti;
            phiq[this->nSpecie() + 1] = pi;
            if (tabulation_->variableTimeStep())
            {
                phiq[this->nSpecie() + 2] = deltaT[celli];
            }

            // Initialise time progress
            scalar timeLeft = deltaT[celli];

            // Not sure if this is necessary
            Rphiq = Zero;

            clockTime_.timeIncrement();

            // When tabulation is active (short-circuit evaluation for retrieve)
            // It first tries to retrieve the solution of the system with the
            // information stored through the tabulation method
            if (tabulation_->active() && tabulation_->retrieve(phiq, Rphiq))
            {
                // Retrieved solution stored in Rphiq
                for (label i = 0; i < this->nSpecie(); i++)
                {
                    c[i] = rhoi * Rphiq[i] / this->specieThermo_[i].W();
                }

                searchISATCpuTime_ += clockTime_.timeIncrement();
            }
            // This position is reached when tabulation is not used OR
            // if the solution is not retrieved.
            // In the latter case, it adds the information to the tabulation
            // (it will either expand the current data or add a new stored point).
            else
            {
                // Store total time waiting to attribute to add or grow
                scalar timeTmp = clockTime_.timeIncrement();

                if (reduced)
                {
                    // Reduce mechanism change the number of species (only active)
                    mechRed_->reduceMechanism(c, Ti, pi);
                    nActiveSpecies += mechRed_->NsSimp();
                    nAvg++;
                    scalar timeIncr = clockTime_.timeIncrement();
                    reduceMechCpuTime_ += timeIncr;
                    timeTmp += timeIncr;
                }

                // Calculate the chemical source terms
                while (timeLeft > small)
                {
                    scalar dt = timeLeft;
                    if (reduced)
                    {
                        // completeC_ used in the overridden ODE methods
                        // to update only the active species
                        completeC_ = c;

                        // Solve the reduced set of ODE
                        this->solve(
                            simplifiedC_, Ti, pi, dt, this->deltaTChem_[celli]);

                        for (label i = 0; i < NsDAC_; i++)
                        {
                            c[simplifiedToCompleteIndex_[i]] = simplifiedC_[i];
                        }
                    }
                    else
                    {
                        this->solve(c, Ti, pi, dt, this->deltaTChem_[celli]);
                    }
                    timeLeft -= dt;
                }

                {
                    scalar timeIncr = clockTime_.timeIncrement();
                    solveChemistryCpuTime_ += timeIncr;
                    timeTmp += timeIncr;
                }

                // If tabulation is used, we add the information computed here to
                // the stored points (either expand or add)
                if (tabulation_->active())
                {
                    forAll(c, i)
                    {
                        Rphiq[i] = c[i] / rhoi * this->specieThermo_[i].W();
                    }
                    if (tabulation_->variableTimeStep())
                    {
                        Rphiq[Rphiq.size() - 3] = Ti;
                        Rphiq[Rphiq.size() - 2] = pi;
                        Rphiq[Rphiq.size() - 1] = deltaT[celli];
                    }
                    else
                    {
                        Rphiq[Rphiq.size() - 2] = Ti;
                        Rphiq[Rphiq.size() - 1] = pi;
                    }
                    label growOrAdd =
                        tabulation_->add(phiq, Rphiq, rhoi, deltaT[celli]);
                    if (growOrAdd)
                    {
                        this->setTabulationResultsAdd(celli);
                        addNewLeafCpuTime_ += clockTime_.timeIncrement() + timeTmp;
                    }
                    else
                    {
                        this->setTabulationResultsGrow(celli);
                        growCpuTime_ += clockTime_.timeIncrement() + timeTmp;
                    }
                }

                // When operations are done and if mechanism reduction is active,
                // the number of species (which also affects nEqns) is set back
                // to the total number of species (stored in the mechRed object)
                if (reduced)
                {
                    this->nSpecie_ = mechRed_->nSpecie();
                }
                deltaTMin = min(this->deltaTChem_[celli], deltaTMin);

                this->deltaTChem_[celli] =
                    min(this->deltaTChem_[celli], this->deltaTChemMax_);
            }

            // Set the RR vector (used in the solver)
            for (label i = 0; i < this->nSpecie_; i++)
            {
                this->RR_[i][celli] =
                    (c[i] - c0[i]) * this->specieThermo_[i].W() / deltaT[celli];
            }
        }
    }

    if (mechRed_->log() || tabulation_->log())
    {
        cpuSolveFile_()
            << this->time().timeOutputValue()
            << "    " << solveChemistryCpuTime_ << endl;
    }

    if (mechRed_->log())
    {
        cpuReduceFile_()
            << this->time().timeOutputValue()
            << "    " << reduceMechCpuTime_ << endl;
    }

    if (tabulation_->active())
    {
        // Every time-step, look if the tabulation should be updated
        tabulation_->update();

        // Write the performance of the tabulation
        tabulation_->writePerformance();

        if (tabulation_->log())
        {
            cpuRetrieveFile_()
                << this->time().timeOutputValue()
                << "    " << searchISATCpuTime_ << endl;

            cpuGrowFile_()
                << this->time().timeOutputValue()
                << "    " << growCpuTime_ << endl;

            cpuAddFile_()
                << this->time().timeOutputValue()
                << "    " << addNewLeafCpuTime_ << endl;
        }
    }

    if (reduced && nAvg && mechRed_->log())
    {
        // Write average number of species
        nActiveSpeciesFile_()
            << this->time().timeOutputValue()
            << "    " << nActiveSpecies / nAvg << endl;
    }

    if (Pstream::parRun())
    {
        List<bool> active(composition.active());
        Pstream::listCombineGather(active, orEqOp<bool>());
        Pstream::listCombineScatter(active);

        forAll(active, i)
        {
            if (active[i])
            {
                composition.setActive(i);
            }
        }
    }

    forAll(this->Y(), i)
    {
        if (composition.active(i))
        {
            this->Y()[i].writeOpt() = IOobject::AUTO_WRITE;
        }
    }

    return deltaTMin;
}

template <class ReactionThermo, class ThermoType>
Foam::scalar Foam::STDACChemistryModel<ReactionThermo, ThermoType>::solve(
    const scalar deltaT)
{
    // Don't allow the time-step to change more than a factor of 2
    return min(
        this->solve<UniformField<scalar>>(UniformField<scalar>(deltaT)),
        2 * deltaT);
}

template <class ReactionThermo, class ThermoType>
Foam::scalar Foam::STDACChemistryModel<ReactionThermo, ThermoType>::solve(
    const scalarField &deltaT)
{
    return this->solve<scalarField>(deltaT);
}

template <class ReactionThermo, class ThermoType>
void Foam::STDACChemistryModel<ReactionThermo, ThermoType>::
    setTabulationResultsAdd(
        const label celli)
{
    tabulationResults_[celli] = 0;
}

template <class ReactionThermo, class ThermoType>
void Foam::STDACChemistryModel<ReactionThermo, ThermoType>::
    setTabulationResultsGrow(const label celli)
{
    tabulationResults_[celli] = 1;
}

template <class ReactionThermo, class ThermoType>
void Foam::STDACChemistryModel<ReactionThermo, ThermoType>::
    setTabulationResultsRetrieve(
        const label celli)
{
    tabulationResults_[celli] = 2;
}

// ************************************************************************* //
