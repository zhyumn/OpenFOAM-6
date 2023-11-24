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
    Foam::parallelISAT

Description

SourceFiles
    ISATI.H
    ISAT.C
    ISATIO.C

\*---------------------------------------------------------------------------*/

#ifndef PARALLELISAT_H
#define PARALLELISAT_H

#include "SUPstream.H"
#include <fstream>
#include "IOstreams.H"
#include "slab.H"
#include "parData.H"
#include "parISATleafNode.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace chemistryTabulationMethodSs
    {
        // Forward declaration of classes
        class Istream;
        class Ostream;

        // Forward declaration of friend functions and operators
        template <typename T>
        class parallelISAT;

        /*---------------------------------------------------------------------------*\
                             Class parallelISAT Declaration
    \*---------------------------------------------------------------------------*/
        template <typename T>
        class parallelISAT
        {

        public:
            //typedef typename T::ISAT_tmp::SList<scalar> SList<scalar>;
            //typedef typename T::ISAT_tmp::outputType outputType;

            T &ISATmanager;
            //protected:
            // Private data

            //- Description of data_
            SUPstream::mpi_manager &manager_;
            SUPstream::mpi_sync &sync_;
            label maxNLeafs_;
            //label nmem_leaf_;
            //label nmem_node_;

            scalar tolerance_;

            // Private Member Functions

            //- Disallow default bitwise copy construct
            parallelISAT(const parallelISAT &);

            //- Disallow default bitwise assignment
            void operator=(const parallelISAT &);

            template <class Type>
            class Mem_manager
            {
                label pemptylist_;
                //SUPstream::shared_data<label> &emptylist_;
                label &tail_;
                label &size_;
                //label &tmp_size_;

            public:
                Mem_manager(label &tail_in, label &size_in /*, label &tmp_size_in*/) : pemptylist_(-1), tail_(tail_in), size_(size_in) /*, tmp_size_(tmp_size_in)*/ { tail_ = 0; }
                inline SharedPointer<Type> New(int size_in)
                {
                    size_++;
                    auto ret = pslab->alloc(sizeof(Type));
                    if (ret == sptr_NULL)
                    {
                        FatalErrorInFunction << "Run out of memory\n"
                                             << exit(FatalError);
                    }
                    (*SharedPointer<Type>(ret)).init(size_in);
                    return SharedPointer<Type>(ret);
                }
                inline SharedPointer<Type> New(const Type &x)
                {
                    size_++;
                    auto ret = pslab->alloc(sizeof(Type));
                    if (ret == sptr_NULL)
                    {
                        FatalErrorInFunction << "Run out of memory\n"
                                             << exit(FatalError);
                    }
                    (*SharedPointer<Type>(ret)).init(x);
                    return SharedPointer<Type>(ret);
                }
                inline void Delete(SharedPointer<Type> &in)
                {
                    size_--;
                    //tmp_size_++;
                    ///emptylist_[tail_] = in.offset;
                    //tail_++;
                    in->setNext(pemptylist_);
                    pemptylist_ = in.offset;
                }
                inline void update()
                {

                    SharedPointer<Type> ptr, tmp;
                    ptr.offset = pemptylist_;
                    while (ptr.offset != -1)
                    {
                        tmp = ptr;
                        ptr.offset = ptr->next();
                        tmp->Delete();
                        pslab->free(sizeof(Type), tmp.offset);
                    }

                    //Pout << "tail_ =  " << tail_<<endl;

                    //tail_ = 0;
                    pemptylist_ = -1;
                }
                /*                 inline void report()
                {
                    std::cout << "Mem_manager" << std::endl;
                    std::cout << "Total:" << (tmp_size_ + size_) * sizeof(Type) << " Used:" << size_ * sizeof(Type) << " Temp:" << (tmp_size_) * sizeof(Type) << std::endl;
                } */
            };
            struct Private_variable
            {
                label size_leaf;
                //label head_leaf;
                label tail_leaf;
                label tmp_size_leaf;
                //label temp_tail_leaf;
                label m_size_leaf;

                label size_node;
                //label heaf_node;
                label tail_node;
                label tmp_size_node;
                //label temp_tail_node;
                label m_size_node;
                std::atomic<int> leafN;

                //label root;
                SharedPointer<Node> root_;
                //std::atomic<size_t> root_;

                //label totalGrowth;
            };

            SUPstream::shared_data<Private_variable> pv_;
            label &size_leaf_;
            //label &head_leaf;
            label &tail_leaf;
            //label &temp_tail_leaf;

            std::atomic<int> &leafN;

            label &size_node;
            //label &head_node;
            label &tail_node;
            //label &temp_tail_node;

            SharedPointer<Node> &root_;
            //std::atomic<size_t> *root_atomic;
            //label &totalGrowth_;

            SUPstream::mpi_mutex write_lock;
            SUPstream::mpi_mutex mem_lock;

            Mem_manager<Leaf> leaf_manager;
            Mem_manager<Node> node_manager;

            label nRetrieved_;
            label nGrowth_;
            label nAdd_;

        public:
            // Static data members

            //- Static data staticData
            //static const dataType staticData;

            // Constructors

            //- Construct null
            //parallelISAT();

            //- Construct from components
            template <typename... Args>
            parallelISAT(T &ISATmanager_in, SUPstream::mpi_manager &manager_in, /*label nmem,*/ label maxNleafs, SUPstream::mpi_sync &sync_in, Args... args);

            //- Construct from Istream
            //ISAT(Istream&);

            //- Construct as copy
            //ISAT(const ISAT&);

            // Selectors

            //- Select null constructed
            //static autoPtr<ISAT> New();

            //- Destructor
            ~parallelISAT();

            // Member Functions

            // Access
            bool isFull();
            const label &size_leaf() { return size_leaf_; }
            SharedPointer<Leaf> binaryTreeSearch(const SList<scalar> &value, SharedPointer<Node> node_index);
            SharedPointer<Leaf> binaryTreeSearch(const scalarField &value);
            SharedPointer<Leaf> binaryTreeSearch(const scalarField &value, SharedPointer<Node> node_index);
            SharedPointer<Leaf> binaryTreeSearch(const scalarField &value, SharedPointer<Node> node_index, SharedPointer<Node> &node_ret, int &side);
            //void writeDot(string name);
            //void writeDot(label nodeindex, std::ofstream &fout);
            /*             scalar leafvalue(SharedPointer<Leaf> index)
            {
                if (index.isNULL())
                {
                    return -1;
                };
                return index->v;
            }; */

            bool retrieve(const scalarField &x, scalarField &y);

            bool retrieve(const scalarField &x, scalarField &y, SharedPointer<Leaf> &pleaf);

            // Check

            // Edit
            bool grow(const scalarField &x, SharedPointer<Leaf> &pleaf);
            bool grow(SharedPointer<Leaf> &pleaf_out, const scalarField &x, SharedPointer<Leaf> &pleaf);
            //bool grow(const SList<scalar> &value);

            template <typename... Args>
            label add(const scalarField &x, const scalarField &y, Args... args);
            template <typename... Args>
            label add(SharedPointer<Leaf> &pleaf_out, bool &growsuc, const scalarField &x, const scalarField &y, Args... args);
            //label add(const SList<scalar> &x);
/*             void report()
            {
                std::cout << "----------------------- report " << size_leaf_ << std::endl;
                pslab->report();
                leaf_manager.report();
                node_manager.report();
                std::cout << "-----------------------  " << std::endl;
            } */
            bool update()
            {
                //report();

                leaf_manager.update();
                node_manager.update();

                return true;
            }
            void clear()
            {
                deleteAll();

                size_leaf_ = 0;

                size_node = 0;
                root_.setNULL();
            }

            inline void valid();
            //serial function
            SharedPointer<Leaf> treeMin(SharedPointer<Node> subTreeRoot);

            inline SharedPointer<Leaf> treeMin()
            {
                return treeMin(root_);
            }

            SharedPointer<Leaf> treeSuccessor(SharedPointer<Leaf> x);

            SharedPointer<Leaf> leafPSibling(SharedPointer<Leaf> x);

            SharedPointer<Node> nodeSibling_leaf(SharedPointer<Leaf> y);
            SharedPointer<Node> nodeSibling_node(SharedPointer<Node> y);

            void transplant(SharedPointer<Node> u, SharedPointer<Node> v);

            label depth(SharedPointer<Node> subTreeRoot);

            inline label depth()
            {
                return depth(root_);
            }

            void deleteLeaf(SharedPointer<Leaf> index_phi0);

            void deleteAllNode(SharedPointer<Node> subTreeRoot);

            inline void deleteAllNode()
            {
                deleteAllNode(root_);
                //node_manager.update();
            }

            void deleteAll(SharedPointer<Node> subTreeRoot);

            inline void deleteAll()
            {
                deleteAll(root_);
                //node_manager.update();
            }

            void insertNode(
                SharedPointer<Leaf> phi0,
                SharedPointer<Node> newNode);
            bool treeCheck(SharedPointer<Node> index);
            inline bool treeCheck()
            {
                if (size_leaf_ == 0)
                {
                    if (root_.isNULL())
                        return true;
                    else
                    {
                        FatalErrorInFunction
                            << "p0";
                        //<< abort(FatalError);
                        return false;
                    }
                }
                if (size_leaf_ == 1)
                {
                    bool ret;
                    if (root_.isNULL())
                    {
                        FatalErrorInFunction
                            << "p1 ";
                        // << abort(FatalError);
                        return false;
                    }
                    ret = (!root_->leafLeft_.isNULL());
                    ret = ret && (root_->nodeLeft_.isNULL());
                    ret = ret && (root_->nodeRight_.isNULL());
                    ret = ret && (root_->leafRight_.isNULL());
                    ret = ret && (root_->parent_.isNULL());
                    ret = ret && (root_->leafLeft_->node_ == root_);
                    if (ret == false)
                    {
                        FatalErrorInFunction
                            << "p2 ";
                        //<< abort(FatalError);
                    }
                    return ret;
                }
                if (root_->parent_.notNULL())
                {
                    FatalErrorInFunction
                        << "p3 ";
                    //<< abort(FatalError);
                    return false;
                }

                return treeCheck(root_);
            }
            /* virtual label add
         (
             const scalarField& phiq,
             const scalarField& Rphiq,
             const scalar rho,
             const scalar deltaT
         );
         */
            // Write

            // Member Operators

            //void operator=(const ISAT&);

            // Friend Functions

            // Friend Operators

            // IOstream Operators

            //friend Istream& operator>>(Istream&, parallelISAT&);
            //friend Ostream& operator<<(Ostream&, const parallelISAT&);
        };
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "parallelISATI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
