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
    Foam::parISATbinaryTree

Description

    Data storage of the data according to a binary
    tree structure.

            0 (root node)
         /     \
        0       0
      /   \   /   \
     L     R L     0
    / \
   L   R

    L: leafLeft_
    R: leafRight_

\*---------------------------------------------------------------------------*/

#ifndef parISATbinaryTree_H
#define parISATbinaryTree_H

#include "dictionary.H"
//#include "ISATleaf.H"
#include "parISATNode.H"
#include "List2D.H"
#include "SUPstream.H"
//#include "binaryNode.H"
//#include "chemPointISAT.H"

namespace Foam
{

    class parISATleaf;

    class parISATNode;

    //template<class CompType, class ThermoType>

    class parISATbinaryTree
    {

    public:
    protected:
        SUPstream::mpi_manager &manager_;
        SUPstream::mpi_sync &sync_;

    public:
        bool constNLeafs_;
        label NLeafBase_;
        scalar NLeafRatio_;
        label NLeafRange_;
        label NLeafMax_;

        label maxTimeTag_;

        //labelList 2dlist_s;

        SUPstream::shared_data<label> emptylist_leaf;
        SUPstream::shared_data<label> emptylist_node;

        template <class Type>
        class Mem_manager
        {
            SUPstream::shared_data<label> &emptylist_;
            int n_in_, n_out_;
            label &tail_;
            label &size_;
            label &tmp_size_;

        public:
            Mem_manager(SUPstream::shared_data<label> &emptylist_in, int n_in, int n_out, label &tail_in, label &size_in, label &tmp_size_in)
                : emptylist_(emptylist_in), n_in_(n_in), n_out_(n_out), tail_(tail_in),
                  size_(size_in), tmp_size_(tmp_size_in) { tail_ = 0; }
            inline SharedPointer<Type> New()
            {
                size_++;
                auto ret = pslab->alloc(sizeof(Type));
                if (ret == sptr_NULL)
                {
                    FatalErrorInFunction << "Run out of memory\n"
                                         << exit(FatalError);
                }
                (*SharedPointer<Type>(ret)).init(n_in_, n_out_);
                return SharedPointer<Type>(ret);
            }
            /*             inline SharedPointer<Type> New(const Type &x)
            {
                size_++;
                auto ret = pslab->alloc(sizeof(Type));
                if (ret == sptr_NULL)
                {
                    FatalErrorInFunction << "Run out of memory\n"
                                         << exit(FatalError);
                }
                (*SharedPointer<Type>(ret)).init(n_in_, n_out_);
                return SharedPointer<Type>(ret);
            } */
            inline void Delete(SharedPointer<Type> &in)
            {
                size_--;
                tmp_size_++;
                emptylist_[tail_] = in.offset;
                tail_++;
            }
            inline void update()
            {
                tmp_size_ = 0;
                while (tail_ > 0)
                {
                    tail_--;
                    (*SharedPointer<Type>(emptylist_[tail_])).Delete();
                    pslab->free(sizeof(Type), emptylist_[tail_]);
                }
            }
            inline void report()
            {
                std::cout << "Mem_manager" << std::endl;
                std::cout << "Total:" << (tmp_size_ + size_) * sizeof(Type) << " Used:" << size_ * sizeof(Type) << " Temp:" << (tmp_size_) * sizeof(Type) << std::endl;
            }
        };

        struct Private_variable
        {
            std::atomic<label> size_leaf;
            label head_leaf;
            label tail_leaf;
            label tmp_size_leaf;
            //label temp_tail_leaf;
            label m_size_leaf;

            std::atomic<label> size_node;
            label heaf_node;
            label tail_node;
            label tmp_size_node;
            //label temp_tail_node;
            label m_size_node;

            std::atomic<label> NRetrieved;
            std::atomic<label> NGrowth;
            std::atomic<label> NAdd;
            std::atomic<label> NCall;

            std::atomic<label> NF1;
            std::atomic<label> NF2;
            std::atomic<label> NF3;
            std::atomic<label> NF4;
            std::atomic<label> NF5;
            std::atomic<label> NF6;
            std::atomic<label> NF7;
            //label root;
            SharedPointer<parISATNode> root_;
            //std::atomic<size_t> root_;

            label totalGrowth;
        };
        SUPstream::shared_data<Private_variable> pv_;
        std::atomic<label> &size_leaf_;
        label &head_leaf;
        label &tail_leaf;

        std::atomic<label> &size_node;
        label &head_node;
        label &tail_node;

        std::atomic<label> &NRetrieved;
        std::atomic<label> &NGrowth;
        std::atomic<label> &NAdd;
        std::atomic<label> &NCall;

        label NRetrievedL;
        label NGrowthL;
        label NAddL;
        label NCallL;

        std::atomic<label> &NF1;
        std::atomic<label> &NF2;
        std::atomic<label> &NF3;
        std::atomic<label> &NF4;
        std::atomic<label> &NF5;
        std::atomic<label> &NF6;
        std::atomic<label> &NF7;

        SUPstream::mpi_mutex write_lock;
        SUPstream::mpi_mutex mem_lock;

        Mem_manager<parISATleaf> leaf_manager;
        Mem_manager<parISATNode> node_manager;

        label n_in_;
        label n_out_;
        //- Reference to the chemistryModel
        // TDACChemistryModel<CompType, ThermoType>& chemistry_;

        //- Root node of the binary tree
        SharedPointer<parISATNode> &root_;

        label NtimeTag_;

        labelList start_s;
        labelList N_s;

        //- Maximum number of elements in the binary tree
        label maxNLeafs_;

        //- Size of the BST (= number of chemPoint stored)
        //label size_;

        //List2D<ISATleaf> timeTagList_;

        void insertNode(
            SharedPointer<parISATleaf> &phi0,
            SharedPointer<parISATNode> &newNode);

        SharedPointer<parISATleaf> Sibling(SharedPointer<parISATleaf> &x);
        SharedPointer<parISATNode> nodeSibling(SharedPointer<parISATleaf> &x);
        //- Replace the binaryNode u with v
        void transplant(SharedPointer<parISATNode> &u, SharedPointer<parISATNode> &v);

    public:
        //- Constructors

        //- Construct from dictionary and chemistryOnLineLibrary
        parISATbinaryTree(
            label in_n,
            label out_n,
            dictionary ISATDict)
            : manager_(SUPstream::node_manager), sync_(SUPstream::Sync),
              emptylist_leaf(manager_, 100000), emptylist_node(manager_, 100000),
              pv_(manager_),
              size_leaf_(pv_().size_leaf), head_leaf(pv_().head_leaf),
              tail_leaf(pv_().tail_leaf),
              size_node(pv_().size_node), head_node(pv_().heaf_node),
              tail_node(pv_().tail_node),
              NRetrieved(pv_().NRetrieved), NGrowth(pv_().NGrowth), NAdd(pv_().NAdd), NCall(pv_().NCall),
              NRetrievedL(0), NGrowthL(0), NAddL(0), NCallL(0),
              NF1(pv_().NF1), NF2(pv_().NF2), NF3(pv_().NF3),
              NF4(pv_().NF4), NF5(pv_().NF5), NF6(pv_().NF6),
              NF7(pv_().NF7),
              write_lock(manager_, sync_), mem_lock(manager_, sync_),
              leaf_manager(emptylist_leaf, in_n, out_n, tail_leaf, pv_().m_size_leaf, pv_().tmp_size_leaf),
              node_manager(emptylist_node, in_n, out_n, tail_node, pv_().m_size_node, pv_().tmp_size_node),
              n_in_(in_n), n_out_(out_n), root_(pv_().root_), NtimeTag_(readLabel(ISATDict.lookup("NtimeTag"))),
              start_s(NtimeTag_), N_s(NtimeTag_) //, size_(0) //, timeTagList_(readLabel(ISATDict.lookup("NtimeTag")))
        {
            sync_.sync();
            if (manager_.rank == 0)
            {
                root_.setNULL();
                NRetrieved = 0;
                NGrowth = 0;
                NAdd = 0;
                NCall = 0;

                NF1 = 0;
                NF2 = 0;
                NF3 = 0;
                NF4 = 0;
                NF5 = 0;
                NF6 = 0;
                NF7 = 0;
            }
            sync_.sync();

            label timeTagList = readLabel(ISATDict.lookup("NtimeTag"));
            dictionary &NleafDcit = ISATDict.subDict("maxNLeafs");
            if (word(NleafDcit.lookup("type")) == word("const"))
            {
                constNLeafs_ = true;
                maxNLeafs_ = readLabel(NleafDcit.lookup("value"));
                maxTimeTag_ = 0;
            }
            else if (word(NleafDcit.lookup("type")) == word("adaptive"))
            {
                constNLeafs_ = false;
                NLeafBase_ = readLabel(NleafDcit.lookup("base"));
                NLeafRatio_ = readScalar(NleafDcit.lookup("ratio"));
                NLeafRange_ = readLabel(NleafDcit.lookup("range"));
                NLeafRange_ = min(NLeafRange_, timeTagList - 1);
                NLeafMax_ = readLabel(NleafDcit.lookup("max"));
                maxNLeafs_ = NLeafMax_;
                maxTimeTag_ = 0;
            }
            else if (word(NleafDcit.lookup("type")) == word("maxTimeTag"))
            {
                constNLeafs_ = false;
                NLeafBase_ = readLabel(NleafDcit.lookup("base"));
                NLeafRatio_ = readScalar(NleafDcit.lookup("ratio"));
                NLeafRange_ = readLabel(NleafDcit.lookup("range"));
                NLeafRange_ = min(NLeafRange_, timeTagList - 1);
                NLeafMax_ = readLabel(NleafDcit.lookup("max"));
                maxNLeafs_ = NLeafMax_;
                maxTimeTag_ = readLabel(NleafDcit.lookup("maxTimeTag"));
            }
            else
            {
                FatalErrorInFunction
                    << "maxNLeafs type can only be const and adaptive."
                    << exit(FatalError);
            }
        };

        inline label maxNLeaf()
        {
            if (constNLeafs_ == false)
            {
                label sum = 0;
                for (int i = 0; i < NLeafRange_; i++)
                    sum += N_s[i];
                //label sum = timeTagList_.sizeSum(NLeafRange_);
                if (maxTimeTag_ != 0)
                {
                    label sum2 = 0;
                    for (int i = 0; i < maxTimeTag_; i++)
                        sum2 += N_s[i];
                    return min(min(max(NLeafBase_, sum * NLeafRatio_), NLeafMax_), sum2);
                }
                else
                {
                    return min(max(NLeafBase_, sum * NLeafRatio_), NLeafMax_);
                }
            }
            return maxNLeafs_;
        }

        inline SharedPointer<parISATNode> root()
        {
            return root_;
        }

        inline const SharedPointer<parISATNode> root() const
        {
            return root_;
        }

        inline label n_in() const
        {
            return n_in_;
        }
        inline label n_out() const
        {
            return n_out_;
        }

        /*         inline List2D<ISATleaf> &timeTagList()
        {
            return timeTagList_;
        }

        inline const List2D<ISATleaf> &timeTagList() const
        {
            return timeTagList_;
        } */
        void deleteSubTree(SharedPointer<parISATNode> subTreeRoot);

        inline void deleteSubTree()
        {
            deleteSubTree(root_);
        }

        inline void deleteSubTree_Node()
        {
            deleteSubTree_Node(root_);
        }

        void deleteSubTree_Node(SharedPointer<parISATNode> subTreeRoot);

        inline label size() const
        {
            return size_leaf_;
        }
        //- Computes iteratively the depth of the subTree
        label depth(SharedPointer<parISATNode> &subTreeRoot) const;

        inline label depth() const
        {
            return depth(root_);
        }
        void deleteLeaf(SharedPointer<parISATleaf> pleaf);

        //- Member functions

        inline label maxNLeafs()
        {
            return maxNLeafs_;
        }
        inline label maxNLeafs() const
        {
            return maxNLeafs_;
        }
        inline void setmaxNLeafs(label n)
        {
            maxNLeafs_ = n;
        }

        // Insert a new leaf starting from the parent node of phi0
        // Parameters: phi0 the leaf to replace by a node
        // phiq the new composition to store
        // Rphiq the mapping of the new composition point
        // A the mapping gradient matrix
        // B the matrix used to initialize the EOA
        // nCols the size of the matrix
        // Returns: void
        // Description :
        //1) Create a new leaf with the data to initialize the EOA and to
        // retrieve the mapping by linear interpolation (the EOA is
        // initialize in the chemPoint constructor)
        //2) Get the parent node of phi0 and connect a new node in place of the
        // leaf of phi0. This new node is constructed with phi0 on the left
        // and phiq on the right (the hyperplane is computed inside the
        // binaryNode constructor)

        //SharedPointer<parISATleaf> insertNewLeaf(
        //    const scalarList &value,
        //    const scalarList &data);

        // Search the parISATbinaryTree until the nearest leaf of a specified
        // leaf is found.
        void binaryTreeSearch(
            const scalarList &value,
            SharedPointer<parISATNode> &node,
            SharedPointer<parISATleaf> &nearest);

        void binaryTreeSearch(
            const scalarList &value,
            SharedPointer<parISATNode> node,
            SharedPointer<parISATleaf> &nearest,
            SharedPointer<parISATNode> &node_ret, int &side);

        void eval(const scalarList &value, scalarList &ret);
        /*void print
        (
            Ostream& OFout,
            ISATNode* node = nullptr,
            int a = 0
        );
        */

        //- Removes every entries of the tree and delete the associated objects
        void clear();

        ~parISATbinaryTree();

        SharedPointer<parISATNode> balance_build(scalarRectangularMatrix &scaleIn, const SharedPointer<parISATNode> &pNode, int start, int end, List<SharedPointer<parISATleaf>> &arrays);
        void balance(scalarRectangularMatrix &scaleIn);
        //void balance();
        void getMiddle(List<SharedPointer<parISATleaf>> &arrays, int ic, int left, int right, int posleft, int posright);
        //inline void balance(scalarRectangularMatrix & scaleIn);
        void getleaflist(SharedPointer<parISATNode> node, int &tail, List<SharedPointer<parISATleaf>> &list);
        void Nleaflist(SharedPointer<parISATNode> node, int &ret);
        //- ListFull
        bool isFull();

        void clean(scalar timestep);

        void sort(List<SharedPointer<parISATleaf>> &arrays, int left, int right);
        void quicksort(List<SharedPointer<parISATleaf>> &arrays, int left, int right);
        template <typename T>
        friend class ISATmanager;
    };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#ifdef NoRepository
//#include "parISATbinaryTree.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
