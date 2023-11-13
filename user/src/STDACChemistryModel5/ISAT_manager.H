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

Class
    Foam::ISAT_manager

Description

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef ISAT_manager_H
#define ISAT_manager_H

#include "parallelISAT.H"
#include <iostream>
#include "parData.H"
#include "slab.H"
#include "parISATleaf.H"
#include "binaryTreeS2.H"
#include "chemPointISATS.H"
#define ISATcache
//#include "FixedSquareMatrix.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//#define Datasize 56

//#define DateSize1 53
//#define DateSize2 56
#define LIKELY(exp) __builtin_expect(exp, 1)
#define UNLIKELY(exp) __builtin_expect(exp, 0)
namespace Foam
{
    namespace chemistryTabulationMethodSs
    {

        template <class CompType, class ThermoType>
        class ISAT_manager;

        /*---------------------------------------------------------------------------*\
                                 Class ISAT_manager Declaration
        \*---------------------------------------------------------------------------*/
        //template <class CompType, class ThermoType>

        template <class CompType, class ThermoType>
        class ISAT_manager
            : public chemistryTabulationMethodS<CompType, ThermoType>
        {
        public:
            parallelISAT<ISAT_manager<CompType, ThermoType>> tree;
#ifdef ISATcache
            binaryTreeS<CompType, ThermoType> chemisTree_;
            bool cleaningRequired_;
#endif

            //- Number of equations in addition to the species eqs.
            label nAdditionalEqns_;

            const Time &runTime_;

            //scalar tolerance_;

            label completeSpaceSize_;

            label iddeltaT_;

            scalarField scaleFactor_;

            bool printProportion_;

            label maxNumNewDim_;

            label idT_;
            label idp_;

            label chPMaxLifeTime_;
            label maxGrowth_;

            scalar minBalanceThreshold_;

            scalar maxDepthFactor_;

            autoPtr<OFstream> nRetrievedFile_;
            autoPtr<OFstream> nGrowthFile_;
            autoPtr<OFstream> nAddFile_;
            autoPtr<OFstream> sizeFile_;
            autoPtr<OFstream> totalGrowthFile_;

        public:
            TypeName("parallelISAT");

            //- Construct from components
            ISAT_manager(const dictionary &dict, STDACChemistryModel<CompType, ThermoType> &chemistry);

            //- Destructor
            ~ISAT_manager();

            // Write
            void balance();

            bool cleanAndBalance();

            bool update()
            {
                bool ret = false;
                tree.sync_.sync();
                if (tree.manager_.rank == 0)
                {
                    //parallelISAT<ISAT_chem<CompType, ThermoType>, emptyClass>::update();
                    if (this->chemistry_.timeSteps() % 20 == 0)
                        ret = cleanAndBalance();
                    if (tree.isFull())
                    {
                        tree.clear();
                        ret = true;
                    }
                    //Pout << "\n???here!"<< ","<< "\n"<< endl;
                    //parallelISAT<ISAT_chem<CompType, ThermoType>, emptyClass>::update();
                    tree.update();
                }
                tree.sync_.sync();
                return ret;
            }
            void setLeaf(
                Leaf &leaf_input,
                const scalarField &x,
                const scalarField &y,
                SUPstream::mpi_mutex &mem_mutex,
                const scalar rhoi,
                const scalar dt);

            void computeA(
                const scalarField &phiq,
                const scalarField &Rphiq,
                SSquareMatrix<scalar> &A_out,
                scalarSquareMatrix &A,
                const scalar rhoi,
                const scalar dt, SUPstream::mpi_mutex &mem_mutex);

            void computeLT(Leaf &leaf_input, scalarSquareMatrix &A, SUPstream::mpi_mutex &mem_mutex);
            void qrDecompose(
                const label nCols,
                scalarSquareMatrix &R);

            bool grow(Leaf &leaf_input, const scalarField &phiq, SUPstream::mpi_mutex &mem_mutex);
            bool grow(chemPointISATS<CompType, ThermoType> *phi0, const scalarField &phiq, const scalarField &Rphiq);
            bool sync()
            {
                bool ret = false;
                tree.sync_.sync();
                return ret;
            }
            label size() { return tree.size_leaf(); };
            void writePerformance()
            {

                if (this->log())
                {
                    nRetrievedFile_()
                        << runTime_.timeOutputValue() << "    " << tree.nRetrieved_ << endl;
                    tree.nRetrieved_ = 0;

                    nGrowthFile_()
                        << runTime_.timeOutputValue() << "    " << tree.nGrowth_ << endl;
                    tree.nGrowth_ = 0;

                    nAddFile_()
                        << runTime_.timeOutputValue() << "    " << tree.nAdd_ << endl;
                    tree.nAdd_ = 0;

                    sizeFile_()
                        << runTime_.timeOutputValue() << "    " << size() << endl;

                    totalGrowthFile_()
                        << runTime_.timeOutputValue() << "    " << tree.totalGrowth_ << endl;
                    tree.sync_.sync();
                    tree.totalGrowth_ = 0;
                }
            }
            bool retrieve(
                const scalarField &phiq,
                scalarField &Rphiq)
            {
#ifdef ISATcache
                return chemisTree_.retrieve(phiq, Rphiq);
#endif
                return tree.retrieve(phiq, Rphiq);
            };

            void computeA(
                scalarSquareMatrix &A,
                const scalarField &Rphiq,
                const scalar rho,
                const scalar dt);
            label add(
                const scalarField &phiq,
                const scalarField &Rphiq,
                const scalar rho,
                const scalar deltaT)
            {

#ifdef ISATcache
                chemisTree_.add(phiq, Rphiq, rho, deltaT);

#endif

                return tree.add(phiq, Rphiq, rho, deltaT);
            };

            //friend ISAT_chem<CompType, ThermoType>;
        };

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    } // End namespace chemistryTabulationMethodSs
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#ifdef NoRepository
//#include "ISAT_manager.C"
//#endif

#include "ISAT_managerI.H"
//#include "ISAT_managerI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
