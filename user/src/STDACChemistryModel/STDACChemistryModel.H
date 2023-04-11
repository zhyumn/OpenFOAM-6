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

Class
    Foam::STDACChemistryModel

Description
    Extends StandardChemistryModel by adding the STDAC method.

    References:
    \verbatim
        Contino, F., Jeanmart, H., Lucchini, T., & D’Errico, G. (2011).
        Coupling of in situ adaptive tabulation and dynamic adaptive chemistry:
        An effective method for solving combustion in engine simulations.
        Proceedings of the Combustion Institute, 33(2), 3057-3064.

        Contino, F., Lucchini, T., D'Errico, G., Duynslaegher, C.,
        Dias, V., & Jeanmart, H. (2012).
        Simulations of advanced combustion modes using detailed chemistry
        combined with tabulation and mechanism reduction techniques.
        SAE International Journal of Engines,
        5(2012-01-0145), 185-196.

        Contino, F., Foucher, F., Dagaut, P., Lucchini, T., D’Errico, G., &
        Mounaïm-Rousselle, C. (2013).
        Experimental and numerical analysis of nitric oxide effect on the
        ignition of iso-octane in a single cylinder HCCI engine.
        Combustion and Flame, 160(8), 1476-1483.

        Contino, F., Masurier, J. B., Foucher, F., Lucchini, T., D’Errico, G., &
        Dagaut, P. (2014).
        CFD simulations using the STDAC method to model iso-octane combustion
        for a large range of ozone seeding and temperature conditions
        in a single cylinder HCCI engine.
        Fuel, 137, 179-184.
    \endverbatim

SourceFiles
    STDACChemistryModelI.H
    STDACChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef STDACChemistryModel_H
#define STDACChemistryModel_H

#include "StandardChemistryModel.H"
#include "chemistryReductionMethodS.H"
#include "chemistryTabulationMethodS.H"
#include "OFstream.H"
#include "SUPstream.H"

#define DateSize1 53
#define DateSize2 56
//#define Batch_Size 100 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class STDACChemistryModel Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionThermo, class ThermoType>
class STDACChemistryModel
:
    public StandardChemistryModel<ReactionThermo, ThermoType>
{
    // Private member data

        bool variableTimeStep_;

        label timeSteps_;

        // Mechanism reduction
        label NsDAC_;
        scalarField completeC_;
        scalarField simplifiedC_;
        Field<bool> reactionsDisabled_;
        List<List<specieElement>> specieComp_;
        Field<label> completeToSimplifiedIndex_;
        DynamicList<label> simplifiedToCompleteIndex_;
        autoPtr<chemistryReductionMethodS<ReactionThermo, ThermoType>>
            mechRed_;

        // Tabulation
        autoPtr<chemistryTabulationMethodS<ReactionThermo, ThermoType>>
            tabulation_;

        // Log file for the average time spent reducing the chemistry
        autoPtr<OFstream> cpuReduceFile_;

        // Write average number of species
        autoPtr<OFstream> nActiveSpeciesFile_;

        //- Log file for the average time spent adding tabulated data
        autoPtr<OFstream> cpuAddFile_;

        //- Log file for the average time spent growing tabulated data
        autoPtr<OFstream> cpuGrowFile_;

        //- Log file for the average time spent retrieving tabulated data
        autoPtr<OFstream> cpuRetrieveFile_;

        autoPtr<OFstream> cpuTotalFile_;

        //- Log file for average time spent solving the chemistry
        autoPtr<OFstream> cpuSolveFile_;

        // Field containing information about tabulation:
        // 0 -> add (direct integration)
        // 1 -> grow
        // 2 -> retrieve
        volScalarField tabulationResults_;


    // Private Member Functions

        //- Disallow copy constructor
        STDACChemistryModel(const STDACChemistryModel&);

        //- Disallow default bitwise assignment
        void operator=(const STDACChemistryModel&);

        //- Solve the reaction system for the given time step
        //  of given type and return the characteristic time
        //  Variable number of species added
        template<class DeltaTType>
        scalar solve(const DeltaTType& deltaT);

        template<class DeltaTType>
        scalar solve2(const DeltaTType& deltaT);


public:

    //- Runtime type information
    TypeName("STDAC");


    // Constructors

        //- Construct from thermo
        STDACChemistryModel(ReactionThermo& thermo);


    //- Destructor
    virtual ~STDACChemistryModel();


    // Member Functions

        //- Return true if the time-step is variable and/or non-uniform
        inline bool variableTimeStep() const;

        //- Return the number of chemistry evaluations (used by ISAT)
        inline label timeSteps() const;

        //- Create and return a STDAC log file of the given name
        inline autoPtr<OFstream> logFile(const word& name) const;

        inline PtrList<volScalarField>& Y();

        //- dc/dt = omega, rate of change in concentration, for each species
        virtual void omega
        (
            const scalarField& c,
            const scalar T,
            const scalar p,
            scalarField& dcdt
        ) const;

        //- Return the reaction rate for reaction r and the reference
        //  species and charateristic times
        virtual scalar omega
        (
            const Reaction<ThermoType>& r,
            const scalarField& c,
            const scalar T,
            const scalar p,
            scalar& pf,
            scalar& cf,
            label& lRef,
            scalar& pr,
            scalar& cr,
            label& rRef
        ) const;


        // Chemistry model functions (overriding functions in
        // StandardChemistryModel to use the private solve function)

            //- Solve the reaction system for the given time step
            //  and return the characteristic time
            virtual scalar solve(const scalar deltaT);

            //- Solve the reaction system for the given time step
            //  and return the characteristic time
            virtual scalar solve(const scalarField& deltaT);


        // ODE functions (overriding functions in StandardChemistryModel to take
        // into account the variable number of species)

            virtual void derivatives
            (
                const scalar t,
                const scalarField& c,
                scalarField& dcdt
            ) const;

            virtual void jacobian
            (
                const scalar t,
                const scalarField& c,
                scalarField& dcdt,
                scalarSquareMatrix& J
            ) const;

            virtual void solve
            (
                scalarField& c,
                scalar& T,
                scalar& p,
                scalar& deltaT,
                scalar& subDeltaT
            ) const = 0;

            struct jobdata
            {
                scalar rho;
                scalar p;
                scalar T;
                FixedList<scalar, DateSize1> c;
                FixedList<scalar, DateSize1> c0; 
                FixedList<scalar, DateSize2> phiq;
                FixedList<scalar, DateSize2> Rphiq;
                scalar deltaTChem_; 
                scalar deltaT;
                label growOrAdd; 
                label celli;
                FixedList<scalar, DateSize1> RR;
            };

            struct jobInput
            {
                scalar rho;
                scalar p;
                scalar T;
                FixedList<scalar, DateSize1> c;
                FixedList<scalar, DateSize1> c0; 
                FixedList<scalar, DateSize2> phiq;
                FixedList<scalar, DateSize2> Rphiq;
                scalar deltaTChem_; 
                scalar deltaT;
                label growOrAdd; 
                //label celli;
                FixedList<scalar, DateSize1> RR;
            };

            struct jobArray
            {
                //bool finished;
                //bool filled;
                label rank;
                label next;
                label first_i;
                //label last_i;
                label size_i;
                //jobInput jobs[Batch_Size];
                jobInput jobs[1];
                //label celli;
            };

            struct jobArray2
            {
                //bool finished;
                //bool filled;
                label rank;
                label next;
                label first_i;
                //label last_i;
                label size_i;
                //jobInput jobs[Batch_Size];
                jobInput *jobs;

                 jobArray2(int size){jobs=new jobInput[size];}
                 ~jobArray2(){delete[] jobs;}
                //label celli;
            };
            label n_block;
            int  Batch_Size;
            //SUPstream::shared_data<jobArray> sja_;
            SUPstream::shared_data2<jobArray> ja;
            //jobArray *ja;
            label ja_size;
            jobArray2 ja_loc;
            SUPstream::shared_data<std::atomic<int>> sspare_cpu;
            std::atomic<int> &spare_cpu;
            SUPstream::shared_data<int> sempty_head,sempty_tail;
            int & empty_head, & empty_tail;
            SUPstream::shared_data<int> sfilled_head, sfilled_tail;
            int &filled_head, &filled_tail;
            SUPstream::shared_data<int>  sfinished_head;
            int *finished_head;
            int & rank;
            int nfinished_block ;
            bool loadBalance_;

            SUPstream::mpi_mutex l_empty;
            SUPstream::mpi_mutex l_filled;
            SUPstream::mpi_mutex l_finished;

            bool spare;
            label n_cpu;
            //SUPstream::shared_data<jobArray> sja_;
            //jobArray ja;
            SUPstream::mpi_mutex l_receiver;
            SUPstream::mpi_mutex l_sender;
            scalar solve
            (
                const scalar rho,
                scalar p, 
                scalar T,
                scalarField &c,
                scalarField &c0, 
                scalarField &phiq,
                scalarField &Rphiq,
                scalar &deltaTChem_, 
                const scalar &deltaT, 
                label &growOrAdd, 
                label &celli, 
                scalarField &RR
            );
            scalar solve
            (
                jobdata& in
            );

            scalar solve
            (
                jobInput& in
            );


        // Mechanism reduction access functions

            inline void setNsDAC(const label newNsDAC);

            inline void setNSpecie(const label newNs);

            inline scalarField& completeC();

            inline scalarField& simplifiedC();

            inline Field<bool>& reactionsDisabled();

            inline bool active(const label i) const;

            inline void setActive(const label i);

            inline DynamicList<label>& simplifiedToCompleteIndex();

            inline Field<label>& completeToSimplifiedIndex();

            inline const Field<label>& completeToSimplifiedIndex() const;

            inline List<List<specieElement>>& specieComp();

            inline
                autoPtr<chemistryReductionMethodS<ReactionThermo, ThermoType>>&
                mechRed();

            tmp<volScalarField> tabulationResults() const
            {
                return tabulationResults_;
            }

            void setTabulationResultsAdd(const label celli);

            void setTabulationResultsGrow(const label celli);

            void setTabulationResultsRetrieve(const label celli);

            inline void resetTabulationResults();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "STDACChemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "STDACChemistryModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
