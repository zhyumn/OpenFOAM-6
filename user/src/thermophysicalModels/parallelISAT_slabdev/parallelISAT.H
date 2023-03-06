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
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    // Forward declaration of classes
    class Istream;
    class Ostream;

    // Forward declaration of friend functions and operators
    template <typename T, typename T2>
    class parallelISAT;
    //Istream& operator>>(Istream&, parallelISAT&);
    //Ostream& operator<<(Ostream&, const parallelISAT&);

    /*---------------------------------------------------------------------------*\
                             Class parallelISAT Declaration
    \*---------------------------------------------------------------------------*/
    template <typename T, typename T2>
    class parallelISAT : public T2
    {

    public:
        typedef typename T::inputType inputType;
        typedef typename T::outputType outputType;
        typedef typename T::inputType_old inputType_old;
        typedef typename T::outputType_old outputType_old;

    protected:
        // Private data

        //- Description of data_
        //dataType data_;
        SUPstream::mpi_manager &manager_;
        SUPstream::mpi_sync &sync_;
        label maxNLeafs_;
        label nmem_leaf_;
        label nmem_node_;

        scalar tolerance_;

        // Private Member Functions

        //- Disallow default bitwise copy construct
        parallelISAT(const parallelISAT &);

        //- Disallow default bitwise assignment
        void operator=(const parallelISAT &);
        struct Node;
        struct Leaf : public T::leafData
        {

            //Node* node_;
            //label v;
            //scalar v;

            //- Pointer to the chemistryModel object
            //TDACChemistryModel<CompType, ThermoType>& chemistry_;

            //- Vector storing the composition, temperature and pressure
            //  and deltaT if a variable time step is set on
            //scalarField phi_;

            //- Vector storing the mapping of the composition phi
            //scalarField Rphi_;

            //- LT the transpose of the L matrix describing the Ellipsoid Of
            //  Accuracy use List of Lists to be able to change size if DAC is used
            //scalarSquareMatrix LT_;

            //- A the mapping gradient matrix
            //scalarSquareMatrix A_;

            //- Vector storing the scale factor
            //scalarField scaleFactor_;

            //- Reference to the node in the binary tree holding this chemPoint
            // Only for write
            //label node_;
            SharedPointer<Node> node_2;
            void init(const Leaf &x)
            {
                T::leafData::init(x);
                node_2.offset = x.node_2.offset;
            }
            void init(int x)
            {
                T::leafData::init(x);
                //node_2.offset = x.offset;
            }

            //- The size of the composition space (size of the vector phi)
            //label completeSpaceSize_;

            //- Number of times the element has been grown
            //label nGrowth_;

            //- Tolerance for the Ellipsoid of accuracy
            //static scalar tolerance_;

            //- Number of active species stored in the chemPoint
            //label nActiveSpecies_;

            //- Vectors that store the index conversion between the simplified
            //  and the complete chemical mechanism
            //List<label> simplifiedToCompleteIndex_;

            //label timeTag_;
            //label lastTimeUsed_;

            //bool toRemove_;

            //label maxNumNewDim_;

            //Switch printProportion_;

            //- Variable to store the number of retrieves the chemPoint
            //  will generate at each time step
            //label numRetrieve_;

            //- Variable to store the number of time steps the chempoint is allowed
            //   to still live according to the maxChPLifeTime_ parameter
            //label nLifeTime_;

            //List<label> completeToSimplifiedIndex_;

            //- Number of equations in addition to the species eqs.
            //label nAdditionalEqns_;

            //label idT_;
            //label idp_;
            //label iddeltaT_;
        };
        struct Node : public T::nodeData
        {
            //Leaf* leaf_;
            //label v;
            //scalar v;

            SharedPointer<Leaf> leafLeft_2;
            //label leafLeft_;

            //- Element on the right
            SharedPointer<Leaf> leafRight_2;
            //label leafRight_;

            //- Node which follows on the left

            //first check node
            SharedPointer<Node> nodeLeft_2;
            //label nodeLeft_;

            //- Node which follows on the right
            SharedPointer<Node> nodeRight_2;
            //label nodeRight_;

            //- Parent node  Only for write
            SharedPointer<Node> parent_2;
            //label parent_;

            //- Number of equations in addition to the species eqs.
            //label nAdditionalEqns_;

            //scalarField v_;

            //scalar a_;
        };
        //SUPstream::shared_data<Leaf> leafmem;
        SUPstream::shared_data<label> emptylist_leaf;
        //SUPstream::shared_data<Node> nodemem;
        //SUPstream::shared_data<label> emptylist_node;
        SUPstream::shared_data<label> emptylist_node;

        template <class Type>
        class Mem_manager
        {

            label &size_;
            label &tmp_size_;
            SUPstream::shared_data<label> &emptylist_;
            label &tail_;

        public:
            Mem_manager(SUPstream::shared_data<label> &emptylist_in, label &tail_in, label &size_in, label &tmp_size_in) : emptylist_(emptylist_in), tail_(tail_in), size_(size_in), tmp_size_(tmp_size_in) { tail_ = 0; }
            inline SharedPointer<Type> New(int size_in)
            {
                size_++;
                auto ret = pslab->alloc(sizeof(Type));
                (*SharedPointer<Type>(ret)).init(size_in);
                /*
                std::cout << "----------------------- new " << ret << std::endl;
                pslab->report();
                report();
                std::cout << "----------------------- " << std::endl;
                if (size_ == 2)
                    FatalErrorInFunction
                        << "!!!!!!!!"
                        << abort(FatalError);
                        */
                //std::cout << pslab->manager_.rank << "AAAAAAAAAAAAAAAAAA " << ret <<" "<<sizeof(Type)  << std::endl;
                return SharedPointer<Type>(ret);
            }
            inline SharedPointer<Type> New(const Type &x)
            {
                size_++;
                auto ret = pslab->alloc(sizeof(Type));
                (*SharedPointer<Type>(ret)).init(x);
                /*
                std::cout << "----------------------- new " << ret << std::endl;
                pslab->report();
                report();
                std::cout << "----------------------- " << std::endl;
                if (size_ == 2)
                    FatalErrorInFunction
                        << "!!!!!!!!"
                        << abort(FatalError);
                        */
                //std::cout << pslab->manager_.rank << "AAAAAAAAAAAAAAAAAA " << ret <<" "<<sizeof(Type)  << std::endl;
                return SharedPointer<Type>(ret);
            }
            inline void Delete(SharedPointer<Type> &in)
            {
                size_--;
                tmp_size_++;
                //std::cout << "######################## delete" << std::endl;
                /*
                std::cout << "----------------------- delete" << std::endl;
                pslab->report();
                report();
                std::cout << "----------------------- " << std::endl;
                if (size_ == 2)
                    FatalErrorInFunction
                        << "!!!!!!!!"
                        << abort(FatalError);
                        */
                //std::cout << pslab->manager_.rank << "AAAAAAAAAAAAAAAAAA- " << in.offset << std::endl;
                emptylist_[tail_] = in.offset;
                tail_++;
                //pslab->free(sizeof(Type), in.offset);
            }
            inline void update()
            {
                tmp_size_ = 0;
                //std::cout << pslab->manager_.rank << sizeof(Type) << " AAAAAAAAAAAAAAAAAA1 " << tail_ << std::endl;
                while (tail_ > 0)
                {
                    tail_--;
                    //std::cout << "######################## update " <<  <<std::endl;
                    (*SharedPointer<Type>(emptylist_[tail_])).Delete();
                    pslab->free(sizeof(Type), emptylist_[tail_]);
                }

                //std::cout << pslab->manager_.rank << sizeof(Type) << " AAAAAAAAAAAAAAAAAA " << tail_ << std::endl;
            }
            inline void report()
            {
                std::cout << "Mem_manager" << std::endl;
                std::cout << "Total:" << (tmp_size_ + size_) * sizeof(Type) << " Used:" << size_ * sizeof(Type) << " Temp:" << (tmp_size_) * sizeof(Type) << std::endl;
            }
            /*void clear()
            {
                head_ = 0;
                temp_tail_ = nmem_ - 1;
                tail_ = nmem_ - 1;
                for (label i = 0;i < nmem_;i++)
                {
                    emptylist_[i] = i;
                }
            }*/
        };
        struct Private_variable
        {
            //         head      temp_tail       tail
            //         |         |               |
            //    | | |X|X|X|X|X|X|Y|Y|Y|Y|Y|Y|Y|Y| | | |
            label size_leaf;
            label head_leaf;
            label tail_leaf;
            label tmp_size_leaf;
            label temp_tail_leaf;
            label m_size_leaf;

            label size_node;
            label heaf_node;
            label tail_node;
            label tmp_size_node;
            label temp_tail_node;
            label m_size_node;

            label root;
            SharedPointer<Node> root_2;

            label totalGrowth;
            //label zzzz;
            //label zzzz2;
        };

        SUPstream::shared_data<Private_variable> pv_;
        SUPstream::shared_data<label> r_;
        label &size_leaf_;
        label &head_leaf;
        label &tail_leaf;
        label &temp_tail_leaf;

        label &size_node;
        label &head_node;
        label &tail_node;
        label &temp_tail_node;

        label &root_;
        SharedPointer<Node> &root_2;
        label &totalGrowth_;

        //label& zzzz_;
        //label& zzzz2_;
        SUPstream::mpi_mutex write_lock;

        Mem_manager<Node> node_manager;
        Mem_manager<Leaf> leaf_manager;

        label nRetrieved_;
        label nGrowth_;
        label nAdd_;

    public:
        //Node tnode;
        //Leaf tleaf;
        // Static data members

        //- Static data staticData
        //static const dataType staticData;

        // Constructors

        //- Construct null
        //parallelISAT();

        //- Construct from components
        template <typename... Args>
        parallelISAT(SUPstream::mpi_manager &manager_in, label nmem, label maxNleafs, SUPstream::mpi_sync &sync_in, Args... args);

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
        //inline scalar& tolerance() { return tolerance_; }
        SharedPointer<Leaf> binaryTreeSearch(const inputType_old &value, SharedPointer<Node> node_index);
        SharedPointer<Leaf> binaryTreeSearch(const inputType &value, SharedPointer<Node> node_index);
        SharedPointer<Leaf> binaryTreeSearch(const inputType_old &value);
        SharedPointer<Leaf> binaryTreeSearch(const scalarField &value);
        SharedPointer<Leaf> binaryTreeSearch(const scalarField &value, SharedPointer<Node> node_index);
        void writeDot(string name);
        void writeDot(label nodeindex, std::ofstream &fout);
        scalar leafvalue(SharedPointer<Leaf> index)
        {
            if (index.isNULL())
            {
                return -1;
            };
            return index->v;
        };

        bool retrieve(const scalarField &x, scalarField &y);
        bool retrieve(const inputType_old &x, outputType_old &y);
        // Check

        // Edit
        bool grow(const scalarField &x, SharedPointer<Leaf> &pleaf);
        bool grow(const inputType &value);
        //void insert(const inputType& value);

        template <typename... Args>
        label add(const scalarField &x, const scalarField &y, Args... args);
        template <typename... Args>
        label add(const inputType_old &x, const outputType_old &y, Args... args);
        label add(const inputType &x);
        void report()
        {
            std::cout << "----------------------- report " << size_leaf_ << std::endl;
            pslab->d_report();
            leaf_manager.report();
            node_manager.report();
            std::cout << "-----------------------  " << std::endl;
        }
        bool update()
        {
            report();
            /*             std::cout << "----------------------- update1 " << size_leaf_ << std::endl;
            pslab->report();
            leaf_manager.report();
            node_manager.report();
            std::cout << "-----------------------  " << std::endl; */
            leaf_manager.update();
            node_manager.update();
            /*             std::cout << "----------------------- update2 " << size_leaf_ << std::endl;
            pslab->report();
            leaf_manager.report();
            node_manager.report();
            std::cout << "-----------------------  " << std::endl; */
            return true;
        }
        void clear()
        {
            deleteAll();
            //leaf_manager.clear();
            size_leaf_ = 0;
            //node_manager.clear();
            size_node = 0;
            root_2.offset = -1;
        }
        //serial function
        SharedPointer<Leaf> treeMin(SharedPointer<Node> subTreeRoot);

        inline SharedPointer<Leaf> treeMin()
        {
            return treeMin(root_2);
        }

        SharedPointer<Leaf> treeSuccessor(SharedPointer<Leaf> x);

        SharedPointer<Leaf> leafPSibling(SharedPointer<Leaf> x);

        SharedPointer<Node> nodeSibling_leaf(SharedPointer<Leaf> y);
        SharedPointer<Node> nodeSibling_node(SharedPointer<Node> y);

        void transplant(SharedPointer<Node> u, SharedPointer<Node> v);

        label depth(SharedPointer<Node> subTreeRoot);

        inline label depth()
        {
            return depth(root_2);
        }

        void deleteLeaf(SharedPointer<Leaf> index_phi0);

        void deleteAllNode(SharedPointer<Node> subTreeRoot);

        inline void deleteAllNode()
        {
            deleteAllNode(root_2);
            //node_manager.update();
        }

        void deleteAll(SharedPointer<Node> subTreeRoot);

        inline void deleteAll()
        {
            deleteAll(root_2);
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
                if (root_2.isNULL())
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
                if (root_2.isNULL())
                {
                    FatalErrorInFunction
                        << "p1 ";
                    // << abort(FatalError);
                    return false;
                }
                ret = (!root_2->leafLeft_2.isNUL());
                ret = ret && (root_2->nodeLeft_.isNUL());
                ret = ret && (root_2->nodeRight_.isNUL());
                ret = ret && (root_2->leafRight_.isNUL());
                ret = ret && (root_2->parent_.isNUL());
                ret = ret && (root_2->leafLeft_2->node_2.offset == root_2.offset);
                if (ret == false)
                {
                    FatalErrorInFunction
                        << "p2 ";
                    //<< abort(FatalError);
                }
                return ret;
            }
            if (!root_2->parent_2.isNULL())
            {
                FatalErrorInFunction
                    << "p3 ";
                //<< abort(FatalError);
                return false;
            }

            return treeCheck(root_2);
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

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "parallelISATI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
