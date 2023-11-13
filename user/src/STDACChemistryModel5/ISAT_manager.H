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
#include "binaryTreeS.H"
#include "chemPointISAT.H"
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
            binaryTree<CompType, ThermoType> chemisTree_;
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
            bool grow(chemPointISAT<CompType, ThermoType> *phi0, const scalarField &phiq, const scalarField &Rphiq);
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

                chemisTree_.retrieve(phiq, Rphiq);
                bool retrieved(false);
                chemPointISAT<CompType, ThermoType> *phi0;
                // If the tree is not empty
                if (chemisTree_.size())
                {
                    chemisTree_.binaryTreeSearch(phiq, chemisTree_.root(), phi0);

                    // lastSearch keeps track of the chemPoint we obtain by the regular
                    // binary tree search
                    //lastSearch_ = phi0;
                    if (phi0->inEOA(phiq))
                    {
                        retrieved = true;
                    }
                    // After a successful secondarySearch, phi0 store a pointer to the
                    // found chemPoint
                    else if (chemisTree_.secondaryBTSearch(phiq, phi0))
                    {
                        retrieved = true;
                    }
                    //else if (MRURetrieve_)
                    //{
                    //    typename SLList<
                    //        chemPointISAT<CompType, ThermoType> *>::iterator iter = MRUList_.begin();
                    //    for (; iter != MRUList_.end(); ++iter)
                    //    {
                    //        phi0 = iter();
                    //        if (phi0->inEOA(phiq))
                    //        {
                    //            retrieved = true;
                    //            break;
                    //        }
                    //    }
                    //}
                }
                // The tree is empty, retrieved is still false
                //else
                //{
                // There is no chempoints that we can try to grow
                //lastSearch_ = nullptr;
                //}

                if (retrieved)
                {
                    phi0->increaseNumRetrieve();
                    scalar elapsedTimeSteps =
                        this->chemistry_.timeSteps() - phi0->timeTag();

                    // Raise a flag when the chemPoint has been used more than the allowed
                    // number of time steps
                    if (elapsedTimeSteps > chPMaxLifeTime_ && !phi0->toRemove())
                    {
                        cleaningRequired_ = true;
                        phi0->toRemove() = true;
                    }
                    //lastSearch_->lastTimeUsed() = this->chemistry_.timeSteps();
                    //addToMRU(phi0);
                    //calcNewC(phi0, phiq, Rphiq);
                    phi0->retrieve(phiq, Rphiq);
                    chemisTree_.nRetrieved()++;
                    return true;
                }
                else
                {
                    // This point is reached when every retrieve trials have failed
                    // or if the tree is empty
                    return false;
                }

#endif

                //bool flag = parallelISAT<ISAT_chem<CompType, ThermoType>, emptyClass>::retrieve(phiQ, RphiQ);
                //return parallelISAT<ISAT_chem<CompType, ThermoType>, emptyClass>::retrieve(phiQ, RphiQ);
                return tree.retrieve(phiq, Rphiq);
            };

            /*             bool retrieve(
                const FixedList<scalar, DateSize2> &phiQ,
                FixedList<scalar, DateSize2> &RphiQ)
            {
                scalarField phiQ_2(phiQ.size()), RphiQ_2(phiQ.size());
                for (label i = 0; i < phiQ.size(); i++)
                {
                    phiQ_2[i] = phiQ[i];
                    RphiQ_2[i] = RphiQ[i];
                }
                //bool ret = parallelISAT<ISAT_chem<CompType, ThermoType>, emptyClass>::retrieve(phiQ_2, RphiQ_2);
                bool ret = tree.retrieve(phiQ_2, RphiQ_2);

                for (label i = 0; i < phiQ.size(); i++)
                {
                    phiQ[i] = phiQ_2[i];
                    RphiQ[i] = RphiQ_2[i];
                }
                //bool flag = parallelISAT<ISAT_chem<CompType, ThermoType>, emptyClass>::retrieve(phiQ, RphiQ);
                return ret;
            }; */
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

                label growthOrAddFlag = 1;
                // If lastSearch_ holds a valid pointer to a chemPoint AND the growPoints_
                // option is on, the code first tries to grow the point hold by lastSearch_
                //if (lastSearch_ && growPoints_)
                //{
                chemPointISAT<CompType, ThermoType> *phi0 = nullptr;
                chemisTree_.binaryTreeSearch(phiq, chemisTree_.root(), phi0);
                if (grow(phi0, phiq, Rphiq))
                {
                    chemisTree_.nGrowth()++;
                    growthOrAddFlag = 0;
                    //addToMRU(lastSearch_);
                    //the structure of the tree is not modified, return false
                    return growthOrAddFlag;
                }
                //}

                // If the code reach this point, it is either because lastSearch_ is not
                // valid, OR because growPoints_ is not on, OR because the grow operation
                // has failed. In the three cases, a new point is added to the tree.
                /* if (chemisTree_.isFull())
                {
                    // If cleanAndBalance operation do not result in a reduction of the tree
                    // size, the last possibility is to delete completely the tree.
                    // It can be partially rebuild with the MRU list if this is used.
                    if (!cleanAndBalance())
                    {
                        DynamicList<chemPointISAT<CompType, ThermoType> *> tempList;
                        if (maxMRUSize_ > 0)
                        {
                            // Create a copy of each chemPointISAT of the MRUList_ before
                            // they are deleted
                            typename SLList<
                                chemPointISAT<CompType, ThermoType> *>::iterator iter = MRUList_.begin();
                            for (; iter != MRUList_.end(); ++iter)
                            {
                                tempList.append(
                                    new chemPointISAT<CompType, ThermoType>(*iter()));
                            }
                        }
                        chemisTree_.clear();

                        // Pointers to chemPoint are not valid anymore, clear the list
                        //MRUList_.clear();

                        // Construct the tree without giving a reference to attach to it
                        // since the structure has been completely discarded
                        chemPointISAT<CompType, ThermoType> *nulPhi = 0;
                        forAll(tempList, i)
                        {
                            chemisTree_.insertNewLeaf(
                                tempList[i]->phi(),
                                tempList[i]->Rphi(),
                                tempList[i]->A(),
                                scaleFactor_,
                                this->tolerance(),
                                scaleFactor_.size(),
                                nulPhi);
                            deleteDemandDrivenData(tempList[i]);
                        }
                    }

                    // The structure has been changed, it will force the binary tree to
                    // perform a new search and find the most appropriate point still stored
                    //lastSearch_ = nullptr;
                } */

                // Compute the A matrix needed to store the chemPoint.
                label ASize = this->chemistry_.nEqns() + nAdditionalEqns_ - 2;
                scalarSquareMatrix A(ASize, Zero);
                computeA(A, Rphiq, rho, deltaT);
                chemPointISAT<CompType, ThermoType> *NP = nullptr;
                chemisTree_.insertNewLeaf(
                    phiq,
                    Rphiq,
                    A,
                    scaleFactor_,
                    this->tolerance(),
                    scaleFactor_.size(),
                    NP // lastSearch_ may be nullptr (handled by binaryTree)
                );
                /*                 if (lastSearch_ != nullptr)
                {
                    addToMRU(lastSearch_);
                } */
                chemisTree_.nAdd()++;

                return growthOrAddFlag;
#endif
                //label ret = parallelISAT<ISAT_chem<CompType, ThermoType>, emptyClass>::add(phiQ, RphiQ, rho, deltaT);
                //return parallelISAT<ISAT_chem<CompType, ThermoType>, emptyClass>::add(phiQ, RphiQ, rho, deltaT);
                return tree.add(phiq, Rphiq, rho, deltaT);
            };

            /*             label add(
                const FixedList<scalar, DateSize2> &phiQ,
                const FixedList<scalar, DateSize2> &RphiQ,
                const scalar rho,
                const scalar deltaT)
            {
                //label ret = parallelISAT<ISAT_chem<CompType, ThermoType>, emptyClass>::add
                scalarField phiQ_2(phiQ.size()), RphiQ_2(phiQ.size());
                for (label i = 0; i < phiQ.size(); i++)
                {
                    phiQ_2[i] = phiQ[i];
                    RphiQ_2[i] = RphiQ[i];
                }
                //bool ret = parallelISAT<ISAT_chem<CompType, ThermoType>, emptyClass>::add(phiQ_2, RphiQ_2, rho, deltaT);
                bool ret = tree.add(phiQ_2, RphiQ_2, rho, deltaT);
                for (label i = 0; i < phiQ.size(); i++)
                {
                    phiQ[i] = phiQ_2[i];
                    RphiQ[i] = RphiQ_2[i];
                }

                return ret;
            }; */

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
