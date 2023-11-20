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

#include "parISATbinaryTree.H"
#include "SortableList.H"
#include "demandDrivenData.H"
#include "IOstream.H"
#include <fstream>
#define LIKELY(exp) __builtin_expect(exp, 1)
#define UNLIKELY(exp) __builtin_expect(exp, 0)
// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::parISATbinaryTree::insertNode(
    SharedPointer<parISATleaf> &phi0,
    SharedPointer<parISATNode> &newNode)
{
    if (phi0 == phi0->node_->leafRight_) // phi0 is on the right
    {
        phi0->node_->leafRight_ = sptr_NULL;
        phi0->node_->nodeRight_ = newNode;
        return;
    }
    else if (phi0 == phi0->node_->leafLeft_) // phi0 is on the left
    {
        phi0->node_->leafLeft_ = sptr_NULL;
        phi0->node_->nodeLeft_ = newNode;
        return;
    }

    // if we reach this point, there is an issue with the addressing
    FatalErrorInFunction
        << "trying to insert a node with a wrong pointer to a chemPoint"
        << exit(FatalError);
}

Foam::SharedPointer<Foam::parISATleaf> Foam::parISATbinaryTree::Sibling(SharedPointer<parISATleaf> &x)
{
    if (size() > 1)
    {
        if (x == x->node_->leafLeft_)
        {
            // x is on the left, return right side
            // might return nullptr if the right side is a node
            return x->node_->leafRight_;
        }
        else if (x == x->node_->leafRight_)
        {
            // x is on the right, return left side
            return x->node_->leafLeft_;
        }
        else
        {
            FatalErrorInFunction
                << "wrong addressing of the initial leaf"
                << exit(FatalError);
            return SharedPointer<parISATleaf>(sptr_NULL);
        }
    }
    // there is only one leaf attached to the root_, no sibling
    return SharedPointer<parISATleaf>(sptr_NULL);
}

Foam::SharedPointer<Foam::parISATNode> Foam::parISATbinaryTree::nodeSibling(SharedPointer<parISATleaf> &x)
{
    if (size() > 1)
    {
        if (x == x->node_->leafLeft_)
        {
            // x is on the left, return right side
            return x->node_->nodeRight_;
        }
        else if (x == x->node_->leafRight_)
        {
            // x is on the right, return left side
            return x->node_->nodeLeft_;
        }
        else
        {
            FatalErrorInFunction
                << "wrong addressing of the initial leaf"
                << exit(FatalError);
            return SharedPointer<parISATNode>(sptr_NULL);
        }
    }
    return SharedPointer<parISATNode>(sptr_NULL);
}

void Foam::parISATbinaryTree::transplant(SharedPointer<parISATNode> &u, SharedPointer<parISATNode> &v)
{
    if (v != sptr_NULL)
    {
        // u is root_
        if (u->parent_.isNULL())
        {
            root_ = v;
        }
        // u is on the left of its parent
        else if (u == u->parent_->nodeLeft_)
        {
            u->parent_->nodeLeft_ = v;
        }
        // u is ont the right of its parent
        else if (u == u->parent_->nodeRight_)
        {
            u->parent_->nodeRight_ = v;
        }
        else
        {
            FatalErrorInFunction
                << "wrong addressing of the initial node"
                << exit(FatalError);
        }
        v->parent_ = u->parent_;
    }
    else
    {
        FatalErrorInFunction
            << "trying to transplant a nullptr node"
            << exit(FatalError);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::parISATbinaryTree::depth(SharedPointer<parISATNode> &subTreeRoot) const
{
    // when we reach the leaf, we return 0
    if (subTreeRoot.isNULL())
    {
        return 0;
    }
    else
    {
        return 1 + max(
                       depth(subTreeRoot->nodeLeft_),
                       depth(subTreeRoot->nodeRight_));
    }
}

/* Foam::SharedPointer<Foam::parISATleaf>  Foam::parISATbinaryTree::insertNewLeaf(
    const scalarList &value,
    const scalarList &data)
{
    SharedPointer<parISATleaf> newleaf;
    if (size() == 0) // no points are stored
    {
        // create an empty binary node and point root_ to it
        //root_ = new ISATNode();
        root_ = node_manager.New(); //new ISATNode();
        // create the new chemPoint which holds the composition point
        // phiq and the data to initialize the EOA
        newleaf = leaf_manager.New(); //new ISATleaf(n_in_, n_out_, value, root_, data);
        newleaf.init(n_in_, n_out_, value, root_, data);
        root_->leafLeft_ = newleaf;
    }
    else // at least one point stored
    {
        // no reference chemPoint, a BT search is required

        SharedPointer<parISATleaf> pleaf;
        binaryTreeSearch(value, root_, pleaf);

        // access to the parent node of the chemPoint
        SharedPointer<parISATNode> parentNode = pleaf->node_;

        // create the new chemPoint which holds the composition point
        // phiq and the data to initialize the EOA
        newleaf = leaf_manager.New(); //new ISATleaf(n_in_, n_out_, value, nullptr, data);
        newleaf.init(n_in_, n_out_, value, root_, data);
        // insert new node on the parent node in the position of the
        // previously stored leaf (phi0)
        // the new node contains phi0 on the left and phiq on the right
        // the hyper plane is computed in the binaryNode constructor
        SharedPointer<parISATNode> newNode;
        if (size() > 1)
        {
            newNode = node_manager.New(); // new ISATNode(pleaf, newleaf, parentNode);
            newNode.init(pleaf, newleaf, parentNode);
            // make the parent of phi0 point to the newly created node
            insertNode(pleaf, newNode);
        }
        else // size() == 1 (because not equal to 0)
        {
            // when size is 1, the binaryNode is without hyperplane
            deleteDemandDrivenData(root_);
            newNode = node_manager.New(); //new ISATNode(pleaf, newleaf, nullptr);
            newNode.init(pleaf, newleaf, nullptr);
            root_ = newNode;
        }

        pleaf->node_ = newNode;
        newleaf->node_ = newNode;
    }
    size_leaf_++;
    //newleaf->pTimeTagList_ = timeTagList_.insert(newleaf);
    //heap_.insert(newleaf);
    return newleaf;
} */

void Foam::parISATbinaryTree::binaryTreeSearch(
    const scalarList &value,
    SharedPointer<parISATNode> &node,
    SharedPointer<parISATleaf> &nearest)
{
    if (LIKELY(sizeLargerThan1()))
    {
        scalar vPhi = 0.0;
        const SList<scalar> &v = node->v_;
        const scalar &a = node->a_;
        // compute v*phi
        forAll(v, i)
            vPhi += value[i] * v[i];

        if (vPhi > a) // on right side (side of the newly added point)
        {
            if (node->nodeRight_.notNULL())
            {
                binaryTreeSearch(value, node->nodeRight_, nearest);
            }
            else // the terminal node is reached, store leaf on the right
            {
                nearest = node->leafRight_;
                return;
            }
        }
        else // on left side (side of the previously stored point)
        {
            if (node->nodeLeft_.notNULL())
            {
                binaryTreeSearch(value, node->nodeLeft_, nearest);
            }
            else // the terminal node is reached, return element on right
            {
                nearest = node->leafLeft_;
                return;
            }
        }
    }
    // only one point stored (left element of the root)
    else if (size() == 1)
    {
        nearest = root_->leafLeft_;
    }
    else // no point stored
    {
        nearest = sptr_NULL;
    }
}

void Foam::parISATbinaryTree::binaryTreeSearch(
    const scalarList &value,
    SharedPointer<parISATNode> node,
    SharedPointer<parISATleaf> &nearest,
    SharedPointer<parISATNode> &node_ret, int &side)
{
    if (LIKELY(sizeLargerThan1()))
    {
        scalar vPhi = 0.0;
        //const SList<scalar> &v = node->v_;
        //const scalar &a = node->a_;
        // compute v*phi
        forAll(node->v_, i)
            vPhi += value[i] * node->v_[i];

        if (vPhi > node->a_) // on right side (side of the newly added point)
        {
            if (node->nodeRight_.notNULL())
            {
                binaryTreeSearch(value, node->nodeRight_, nearest, node_ret, side);
            }
            else // the terminal node is reached, store leaf on the right
            {
                nearest = node->leafRight_;
                node_ret = node;
                side = 1; //0 left, 1 right
                return;
            }
        }
        else // on left side (side of the previously stored point)
        {
            if (node->nodeLeft_.notNULL())
            {
                binaryTreeSearch(value, node->nodeLeft_, nearest, node_ret, side);
            }
            else // the terminal node is reached, return element on right
            {
                nearest = node->leafLeft_;
                node_ret = node;
                side = 0; //0 left, 1 right
                return;
            }
        }
    }
    // only one point stored (left element of the root)
    else if (size() == 1)
    {
        node_ret = root_;
        side = 0; //0 left, 1 right
        nearest = root_->leafLeft_;
    }
    else // no point stored
    {
        node_ret = sptr_NULL;
        side = 0;
        nearest = sptr_NULL;
    }
}

void Foam::parISATbinaryTree::eval(const scalarList &value, scalarList &ret)
{
    SharedPointer<parISATleaf> pleaf;
    binaryTreeSearch(value, root_, pleaf);
    if (pleaf.isNULL())
        return;
    for (int i = 0; i < pleaf->data_.size(); i++)
        ret[i] = pleaf->data_[i];
    for (int i = 0; i < n_out_; i++)
    {
        for (int j = 0; j < n_in_; j++)
        {
            ret[i] += pleaf->A_(j, i) * (value[j] - pleaf->value_[j]);
        }
    }
}

bool Foam::parISATbinaryTree::isFull()
{
    return size_leaf_ >= maxNLeafs_;
}

void Foam::parISATbinaryTree::deleteSubTree(SharedPointer<parISATNode> subTreeRoot)
{
    if (subTreeRoot.notNULL())
    {
        leaf_manager.Delete(subTreeRoot->leafLeft_);
        leaf_manager.Delete(subTreeRoot->leafRight_);
        deleteSubTree(subTreeRoot->nodeLeft_);
        deleteSubTree(subTreeRoot->nodeRight_);
        node_manager.Delete(subTreeRoot);

        //deleteDemandDrivenData(subTreeRoot->leafLeft_);
        //deleteDemandDrivenData(subTreeRoot->leafRight_);

        //deleteDemandDrivenData(subTreeRoot);
    }
}

void Foam::parISATbinaryTree::deleteSubTree_Node(SharedPointer<parISATNode> subTreeRoot)
{
    if (subTreeRoot.notNULL())
    {
        deleteSubTree_Node(subTreeRoot->nodeLeft_);
        deleteSubTree_Node(subTreeRoot->nodeRight_);
        node_manager.Delete(subTreeRoot);
        //deleteDemandDrivenData(subTreeRoot);
    }
}

void Foam::parISATbinaryTree::clear()
{
    // Recursively delete the element in the subTree
    if (size() > 0)
        deleteSubTree();
    // Reset root node (should already be nullptr)
    root_ = sptr_NULL;

    // Reset size_leaf
    size_leaf_ = 0;
}

Foam::parISATbinaryTree::~parISATbinaryTree()
{
    if (manager_.rank == 0)
        clear();
}

void Foam::parISATbinaryTree::deleteLeaf(SharedPointer<parISATleaf> pleaf)
{
    if (size() == 1) // only one point is stored
    {
        leaf_manager.Delete(pleaf);
        node_manager.Delete(root_);
        //deleteDemandDrivenData(pleaf);
        //deleteDemandDrivenData(root_);
    }
    else if (size() > 1)
    {
        SharedPointer<parISATNode> z = pleaf->node_;
        SharedPointer<parISATNode> x;
        SharedPointer<parISATleaf> siblingpleaf = Sibling(pleaf);

        if (siblingpleaf.notNULL()) // the sibling of phi0 is a chemPoint
        {
            // z was root (only two chemPoints in the tree)
            if (z->parent_.isNULL())
            {
                root_ = node_manager.New(); //new ISATNode();
                root_->leafLeft_ = siblingpleaf;
                siblingpleaf->node_ = root_;
            }
            else if (z == z->parent_->nodeLeft_)
            {
                z->parent_->leafLeft_ = siblingpleaf;
                z->parent_->nodeLeft_ = sptr_NULL;
                siblingpleaf->node_ = z->parent_;
            }
            else if (z == z->parent_->nodeRight_)
            {
                z->parent_->leafRight_ = siblingpleaf;
                z->parent_->nodeRight_ = sptr_NULL;
                siblingpleaf->node_ = z->parent_;
            }
            else
            {
                FatalErrorInFunction
                    << "wrong addressing of the initial leaf"
                    << exit(FatalError);
            }
        }
        else
        {
            x = nodeSibling(pleaf);
            if (x.notNULL())
            {
                transplant(z, x);
            }
            else
            {
                FatalErrorInFunction
                    << "inconsistent structure of the tree, no leaf and no node"
                    << exit(FatalError);
            }
        }
        //deleteDemandDrivenData(pleaf);
        leaf_manager.Delete(pleaf);
        node_manager.Delete(z);
        //deleteDemandDrivenData(z);
    }
    size_leaf_--;
}
void Foam::parISATbinaryTree::getMiddle(List<SharedPointer<parISATleaf>> &arrays, int ic, int left, int right, int posleft, int posright)
{
    if (left == right)
    {
        return;
    }
    int i = left, j = right;
    SharedPointer<parISATleaf> ptemp = arrays[i];
    while (j > i)
    {
        for (; arrays[j]->value_[ic] > ptemp->value_[ic] && j > i; j--)
            ;
        if (i < j)
            arrays[i++] = arrays[j];
        for (; arrays[i]->value_[ic] < ptemp->value_[ic] && j > i; i++)
            ;
        if (i < j)
            arrays[j--] = arrays[i];
    }
    arrays[i] = ptemp;
    if (i < posleft)
    {
        getMiddle(arrays, ic, i + 1, right, posleft, posright);
    }
    else if (i > posright)
    {
        getMiddle(arrays, ic, left, i - 1, posleft, posright);
    }
    else
    {
        if (i > posleft)
        {
            getMiddle(arrays, ic, left, i - 1, posleft, posleft);
        }
        if (i < posright)
        {
            getMiddle(arrays, ic, i + 1, right, posright, posright);
        }
    }
}

void Foam::parISATbinaryTree::quicksort(List<SharedPointer<parISATleaf>> &arrays, int left, int right)
{
    if (left >= right)
    {
        return;
    }
    int i = left, j = right;
    SharedPointer<parISATleaf> init = arrays[i];
    while (j > i)
    {
        for (; arrays[j]->lastUsed < init->lastUsed && j > i; j--)
            ;
        if (i < j)
            arrays[i++] = arrays[j];
        for (; arrays[i]->lastUsed > init->lastUsed && j > i; i++)
            ;
        if (i < j)
            arrays[j--] = arrays[i];
    }
    arrays[i] = init;
    quicksort(arrays, left, i - 1);
    quicksort(arrays, i + 1, right);
}

Foam::SharedPointer<Foam::parISATNode> Foam::parISATbinaryTree::balance_build(scalarRectangularMatrix &scaleIn, const SharedPointer<parISATNode> &pNode, int start, int end, List<SharedPointer<parISATleaf>> &arrays)
{
    scalarList mean(n_in_), variance(n_in_);
    label maxDir(0);
    scalar maxvariance = 0;
    int len = end - start + 1;
    for (int i = 0; i < n_in_; i++)
    {
        mean[i] = 0;
        for (int j = start; j <= end; j++)
        {
            mean[i] += arrays[j]->value_[i];
        }
        mean[i] /= len;
    }
    for (int i = 0; i < n_in_; i++)
    {
        variance[i] = 0;
        for (int j = start; j <= end; j++)
        {
            variance[i] += sqr(arrays[j]->value_[i] - mean[i]);
        }
        variance[i] /= len * sqr(scaleIn[i][i]);
        if (maxvariance < variance[i])
        {
            maxvariance = variance[i];
            maxDir = i;
        }
    }
    int mid = (start + end) / 2;
    getMiddle(arrays, maxDir, start, end, mid, mid + 1);
    SharedPointer<parISATNode> pret = node_manager.New(); //ISATNode(arrays[mid], arrays[mid + 1], pNode, maxDir);
    pret->set(arrays[mid], arrays[mid + 1], pNode, maxDir);
    if (start == mid)
    {
        pret->leafLeft_ = arrays[mid];
        arrays[mid]->node_ = pret;
    }
    else
    {
        pret->nodeLeft_ = balance_build(scaleIn, pret, start, mid, arrays);
    }
    if (end == mid + 1)
    {
        pret->leafRight_ = arrays[end];
        arrays[end]->node_ = pret;
    }
    else
    {
        pret->nodeRight_ = balance_build(scaleIn, pret, mid + 1, end, arrays);
    }
    return pret;
}
void Foam::parISATbinaryTree::Nleaflist(SharedPointer<parISATNode> node, int &ret)
{
    if (node->nodeLeft_.notNULL())
    {
        Nleaflist(node->nodeLeft_, ret);
    }
    else // the terminal node is reached, return element on right
    {
        ret++;
    }
    if (node->nodeRight_.notNULL())
    {
        Nleaflist(node->nodeRight_, ret);
    }
    else // the terminal node is reached, store leaf on the right
    {
        ret++;
    }
}
void Foam::parISATbinaryTree::getleaflist(SharedPointer<parISATNode> node, int &tail, List<SharedPointer<parISATleaf>> &list)
{
    if (node->nodeLeft_.notNULL())
    {
        getleaflist(node->nodeLeft_, tail, list);
    }
    else // the terminal node is reached, return element on right
    {
        list[tail] = node->leafLeft_;
        tail++;
    }
    if (node->nodeRight_.notNULL())
    {
        getleaflist(node->nodeRight_, tail, list);
    }
    else // the terminal node is reached, store leaf on the right
    {
        list[tail] = node->leafRight_;
        tail++;
    }
}

void Foam::parISATbinaryTree::clean(scalar timestep)
{
    if (size() < 2)
        return;

    List<SharedPointer<parISATleaf>> pleaf(size_leaf_);
    int tail = 0;
    getleaflist(root_, tail, pleaf);
    quicksort(pleaf, 0, tail - 1);

    if (pleaf[0]->lastUsed > timestep - 1)
    {
        FatalErrorInFunction
            << "lastUsed is wrong"
            << exit(FatalError);
    }
    int timeN = timestep - 1;

    start_s[0] = 0;

    int j = 0;
    int i = 0;
    for (; i < NtimeTag_; i++)
    {
        start_s[i] = j;
        for (; j < size_leaf_ && pleaf[j]->lastUsed == timestep - i - 1; j++)
            ;
        N_s[i] = j - start_s[i];
    }
    if (j != size_leaf_)
    {
        N_s[NtimeTag_ - 1] = size_leaf_ - start_s[NtimeTag_ - 1];
    }
    Info << "table info:" << endl;
    Info << "size:" << size_leaf_ << endl;

    for (int iter = 0; iter < NtimeTag_; iter++)
    {
        Info << iter << ":" << N_s[iter] << endl;
    }
    label newmax = maxNLeaf();
    /*    i = 0;
     for (; i < NtimeTag_; i++)
    {
        Info << start_s[i] << "," << N_s[i] << endl;
    }

    for (int i = 0; i < size_leaf_; i++)
    {

        Info << pleaf[i]->lastUsed << ",";
    } */
    /*     FatalErrorInFunction
        << "test " << newmax
        << exit(FatalError); */

    for (int i = newmax; i < size_leaf_; i++)
    {
        deleteLeaf(pleaf[i]);
    }
    /*     for (int i = 0; i < size_leaf_; i++)
    {
        if (pleaf[i]->lastUsed < timeN)
        {
            timeN--;
            start_s[j + 1] = i;
            N_s[j] = start_s[j + 1] - start_s[j];
        }
        Info << pleaf[i]->lastUsed << ",";
    } */
}

void Foam::parISATbinaryTree::balance(scalarRectangularMatrix &scaleIn)
{
    if (size() <= 2)
        return;

    /*     int nleaf1 = 0;
    int nleaf2 = 0;
    Nleaflist(root_, nleaf1);

    if (nleaf1 != size_leaf_)
    {
        FatalErrorInFunction
            << "wrong leaf listXX"
            << size_leaf_ << "\n"
            << nleaf1
            << exit(FatalError);
    } */

    List<SharedPointer<parISATleaf>> pleaf(size_leaf_);
    int tail = 0;
    getleaflist(root_, tail, pleaf);

    if (tail != size_leaf_)
    {
        FatalErrorInFunction
            << "wrong leaf list"
            << size_leaf_ << "\n"
            << tail
            << exit(FatalError);
    }
    /*     for (int i = 0; i < size_leaf_; i++)
    {
        out << "(" << pleaf[i].offset << "," << pleaf[i]->value_.size() << ")";
        out.flush();
    }
    out << std::endl; */

    //timeTagList_.leaflist(pleaf.data());
    deleteSubTree_Node();
    root_ = balance_build(scaleIn, SharedPointer<parISATNode>(sptr_NULL), 0, size_leaf_ - 1, pleaf);
    /*     Nleaflist(root_, nleaf2);

    if (nleaf2 != size_leaf_)
    {
        FatalErrorInFunction
            << "wrong leaf listYY"
            << size_leaf_ << "\n"
            << nleaf2
            << exit(FatalError);
    } */
}

// ************************************************************************* //
