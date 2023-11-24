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
    Foam::chemistryTabulationMethods::ISAT

Description
    Implementation of the ISAT (In-situ adaptive tabulation), for chemistry
    calculation.

    Reference:
    \verbatim
        Pope, S. B. (1997).
        Computationally efficient implementation of combustion chemistry using
        in situ adaptive tabulation.
        Combustion Theory and Modelling, 1, 41-63.
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef ISATMANAGER_H
#define ISATMANAGER_H

#include "ISATbinaryTree.H"
#include "parISATbinaryTree.H"
#include <fstream>
//#include "ISATleaf.H"
//#include "ISATNode.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    /*---------------------------------------------------------------------------*\
                                Class ISAT Declaration
    \*---------------------------------------------------------------------------*/

    template <class FuncType>
    class ISATmanager
    {
        // Private data

        //- List of the stored 'points' organized in a binary tree

        label Ninput, Noutput;
        dictionary ISATDict_;
        parISATbinaryTree tableTree_;
        ISATbinaryTree tableTreeL_;

        FuncType *pfunc;
        double epsilon_;
        double relepsilon_;
        scalarRectangularMatrix scaleFactor_;

        scalarRectangularMatrix scaleIn_;
        scalarRectangularMatrix toleranceOut2_;
        scalarRectangularMatrix initToleranceIn2_;

        scalarRectangularMatrix init_elp_;

        bool modified_;
        label nloop_;

        label timeSteps_;
        bool noISAT_;
        label checkInterval_;
        scalar maxDepthFactor_;

        scalarList maxLeafsize_;
        mutable label nRetrieved_;
        mutable label nGrowth_;
        mutable label nAdd_;
        mutable label nCall_;
        mutable label nRRetrieved_;
        mutable label nRRetrieved2_;
        mutable label nRGrowth_;
        mutable label nRGrowth2_;
        mutable label nRGrowth3_;
        mutable label nRAdd_;
        mutable label nRAdd2_;
        mutable label nRAdd3_;
        mutable label nRCall_;

        word treename_;
        bool muted_;

    public:
        label notCall = 0;

        //- Runtime type information
        //TypeName("ISAT");

        // Constructors

        //- Construct from dictionary
        ISATmanager(label in_n, label out_n, FuncType &func, const word &name_in, const dictionary &ISATDict);

        // Destructor
        virtual ~ISATmanager();

        // Member Functions
        inline parISATbinaryTree &tableTree()
        {
            return tableTree_;
        }
        template <class... Args>
        void add(const scalarList &value, scalarList &out, SharedPointer<parISATleaf> &pleaf, bool growflag, bool &growsuc, Args &... arg);

        template <class... Args>
        void addL(const scalarList &value, scalarList &out, Args &... arg);

        void addL_in(const scalarList &value, scalarList &out, SharedPointer<parISATleaf> &pleaf);

        inline void mute()
        {
            muted_ = true;
        }
        inline void unmute()
        {
            muted_ = false;
        }

        inline scalar &epsilon()
        {
            return epsilon_;
        }

        inline const scalar &epsilon() const
        {
            return epsilon_;
        }

        inline scalar &relepsilon()
        {
            return relepsilon_;
        }

        inline const scalar &relepsilon() const
        {
            return relepsilon_;
        }

        SharedPointer<parISATleaf> search(const scalarList &value);
        template <class... Args>
        bool call(const Foam::scalarList &value, scalarList &out, Args &... arg);

        void tablevalue(const Foam::scalarList &value, scalarList &out);

        bool retrieve(const Foam::scalarList &value, scalarList &out);

        //bool grow(SharedPointer<parISATleaf>plf, const scalarList &dvalue, scalarList &data1, const scalarList &data2);
        bool grow2(SharedPointer<parISATleaf> plf, const scalarList &dvalue, scalarList &data1, SharedPointer<parISATleaf> &plf_out);
        bool growL(ISATleaf *plf, const scalarList &dvalue, scalarList &data1);

        double distance(const scalarList &l, const scalarList &r);
        double normalized_distance(const scalarList &l, const scalarList &r);

        double norm(const scalarList &l);

        bool newLoop()
        {
            nloop_++;
            if (nloop_ == 1)
                return false;
            bool ret = modified_;
            modified_ = false;
            return ret;
        }
        void test(SharedPointer<parISATNode> subroot, int level, const string es)
        {
            bool flag = true;
            if (subroot->nodeLeft_.isNULL())
            {
                flag = flag && subroot->leafLeft_.notNULL();
                if (flag)
                {
                    flag = flag && subroot->leafLeft_->value_.size() == Ninput;
                }
                if (flag == false)
                {
                    FatalErrorInFunction << "broken tree at " << level << es << exit(FatalError);
                }
            }
            else
            {
                flag = flag && subroot->leafLeft_.isNULL();
                if (flag == false)
                {
                    FatalErrorInFunction << "broken tree at " << level << es << "\n"
                                         << tableTree_.size() << exit(FatalError);
                }
            }
            if (flag == false)
            {
                FatalErrorInFunction << "broken tree at " << level << es << exit(FatalError);
            }
            if (subroot->nodeRight_.isNULL())
            {
                flag = flag && subroot->leafRight_.notNULL();
                if (flag)
                {
                    flag = flag && subroot->leafRight_->value_.size() == Ninput;
                }
                if (flag == false)
                {
                    FatalErrorInFunction << "broken tree at " << level << es << exit(FatalError);
                }
            }
            else
            {
                flag = flag && subroot->leafRight_.isNULL();
                if (flag == false)
                {
                    FatalErrorInFunction << "broken tree at " << level << es << exit(FatalError);
                }
            }
            if (flag == false)
            {
                FatalErrorInFunction << "broken tree at " << level << es << exit(FatalError);
            }
            if (subroot->nodeLeft_.notNULL())
            {
                test(subroot->nodeLeft_, level + 1, es);
            }
            if (subroot->nodeRight_.notNULL())
            {
                test(subroot->nodeRight_, level + 1, es);
            }
        }
        void test(const string es)
        {
            bool flag = true;
            if (tableTree_.size() == 0)
            {
                flag = flag && tableTree_.root_.isNULL();
                if (flag == false)
                {
                    FatalErrorInFunction << "broken tree at 0" << es << exit(FatalError);
                }
            }
            else if (tableTree_.size() == 1)
            {
                flag = flag && tableTree_.root_->leafRight_.isNULL();
                flag = flag && tableTree_.root_->nodeLeft_.isNULL();
                flag = flag && tableTree_.root_->nodeRight_.isNULL();
                flag = flag && tableTree_.root_->parent_.isNULL();
                flag = flag && tableTree_.root_->leafLeft_->node_ == tableTree_.root_;
                flag = flag && tableTree_.root_->leafLeft_->value_.size() == Ninput;

                if (flag == false)
                {
                    FatalErrorInFunction << "broken tree at 0" << es << exit(FatalError);
                }
            }
            else
            {
                test(tableTree_.root_, 1, es);
            }
        }

        void print(std::ofstream &out, SharedPointer<parISATNode> subroot)
        {
            out << "{";
            if (subroot->nodeLeft_.isNULL())
            {
                out << "(";
                out << subroot->leafLeft_.offset << "," << subroot->leafLeft_->value_.size();
                out << ")";
            }
            else
            {
                print(out, subroot->nodeLeft_);
            }

            if (subroot->nodeRight_.isNULL())
            {
                out << "(";
                out << subroot->leafRight_.offset << "," << subroot->leafRight_->value_.size();
                out << ")";
            }
            else
            {
                print(out, subroot->nodeRight_);
            }

            out << "}";
        }

        void print(std::ofstream &out)
        {
            if (tableTree_.size() <= 1)
            {
                return;
            }
            print(out, tableTree_.root_);
            out << std::endl;
        }

        void newTimeStep()
        {
            nloop_ = 0;
            modified_ = false;
            timeSteps_++;
            tableTreeL_.updateMaxNLeaf();
            //tableTree_.updateMaxNLeaf();

            tableTreeL_.timeTagList().newTimeStep();

            if (tableTreeL_.maxTimeTag_ != 0)
            {
                int n = tableTreeL_.timeTagList().sizeSum(tableTreeL_.maxTimeTag_ + 1) - tableTreeL_.timeTagList().sizeSum(tableTreeL_.maxTimeTag_);

                for (int i = 0; i < n; i++)
                {
                    tableTreeL_.deleteLeaf(tableTreeL_.timeTagList().pop());
                }
            }
            if (tableTreeL_.size() > 1)
            {
                if (timeSteps_ % checkInterval_ == 0 && tableTreeL_.depth() * 1.0 / ceil(log2(tableTreeL_.size())) > maxDepthFactor_)
                {
                    tableTreeL_.balance(scaleIn_);
                    //showPerformance();
                }
            }
            tableTree_.sync_.sync();
            if (tableTree_.manager_.rank == 0)
            {

                if (!(tableTree_.size() == 0 || tableTree_.size() == 1))
                {

                    if (timeSteps_ % checkInterval_ == 0)
                    {
                        tableTree_.clean(timeSteps_);
                    }
                    if (timeSteps_ % checkInterval_ == 0 && tableTree_.depth() * 1.0 / ceil(log2(tableTree_.size())) > maxDepthFactor_)
                    {

                        tableTree_.balance(scaleIn_);
                        showPerformance();
                    }
                }
            }
            if (timeSteps_ % checkInterval_ == 0)
            {
                showPerformance();
            }
            tableTree_.sync_.sync();

            tableTree_.mem_lock.lock();
            tableTree_.leaf_manager.reuse();
            tableTree_.node_manager.reuse();

            tableTree_.leaf_manager.clean();
            tableTree_.node_manager.clean();
            tableTree_.mem_lock.unlock();

            tableTreeL_.updateMaxNLeaf(tableTree_.size());
        }

        label treesize()
        {
            return tableTree_.size();
        }

        label treesize_loc()
        {
            return tableTreeL_.size();
        }

        inline const scalarRectangularMatrix &scaleFactor() const
        {
            return scaleFactor_;
        }

        inline scalarRectangularMatrix &scaleFactor()
        {
            return scaleFactor_;
        }

        inline const scalarRectangularMatrix &scaleIn() const
        {
            return scaleIn_;
        }

        inline scalarRectangularMatrix &scaleIn()
        {
            return scaleIn_;
        }

        void showPerformance() const;
        const word &treename() const
        {
            return treename_;
        }
        inline const scalarRectangularMatrix &init_elp() const
        {
            return init_elp_;
        }
        inline scalarRectangularMatrix &init_elp()
        {
            return init_elp_;
        }
    };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "parISATmanager.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
