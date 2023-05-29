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

#include "parallelISAT.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

//const dataType Foam::ISAT::staticData();

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//Foam::ISAT::ISAT()
//:
//baseClassName(),
//data_()
//{}
template <typename T, typename T2>
template <typename... Args>
Foam::parallelISAT<T, T2>::parallelISAT(SUPstream::mpi_manager &manager_in, label nmem, label maxNleafs, SUPstream::mpi_sync &sync_in, Args... args)
    : T2(args...), manager_(manager_in), sync_(sync_in), maxNLeafs_(maxNleafs), nmem_leaf_(nmem), nmem_node_(nmem),
      pv_(manager_in), emptylist_leaf(manager_in, 100000), emptylist_node(manager_in, 100000),
      size_leaf_(pv_().size_leaf), head_leaf(pv_().head_leaf),
      tail_leaf(pv_().tail_leaf), temp_tail_leaf(pv_().temp_tail_leaf),
      size_node(pv_().size_node), head_node(pv_().heaf_node),
      tail_node(pv_().tail_node), temp_tail_node(pv_().temp_tail_node),
      root_2(pv_().root_2), totalGrowth_(pv_().totalGrowth),
      write_lock(manager_in, sync_in),
      leaf_manager(emptylist_leaf, tail_leaf, pv_().m_size_leaf, pv_().tmp_size_leaf),
      node_manager(emptylist_node, tail_node, pv_().m_size_node, pv_().tmp_size_node),
      nRetrieved_(0),
      nGrowth_(0),
      nAdd_(0)
{
    size_leaf_ = 0;
    tail_leaf = 0;
    tail_node = 0;

    pv_().tmp_size_leaf = 0;
    pv_().tmp_size_node = 0;
    pv_().m_size_leaf = 0;
    pv_().m_size_node = 0;
    size_node = 0;

    root_2.offset = -1;
    totalGrowth_ = 0;
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
template <typename T, typename T2>
Foam::parallelISAT<T, T2>::~parallelISAT()
{
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
template <typename T, typename T2>
Foam::SharedPointer<typename Foam::parallelISAT<T, T2>::Leaf> Foam::parallelISAT<T, T2>::binaryTreeSearch(const scalarField &value)
{
    auto ret = binaryTreeSearch(value, root_2);
    return ret;
}

template <typename T, typename T2>
Foam::SharedPointer<typename Foam::parallelISAT<T, T2>::Leaf> Foam::parallelISAT<T, T2>::binaryTreeSearch(
    const scalarField &value,
    SharedPointer<Node> node_ptr)
{
    if (size_leaf_ > 1)
    {

        if (node_ptr->goLeft(value))
        {
            if (!node_ptr->nodeLeft_2.isNULL())
            {
                return binaryTreeSearch(value, node_ptr->nodeLeft_2);
            }
            else
            {
                return node_ptr->leafLeft_2;
            }
        }
        else
        {
            if (!node_ptr->nodeRight_2.isNULL())
            {
                return binaryTreeSearch(value, node_ptr->nodeRight_2);
            }
            else
            {
                return node_ptr->leafRight_2;
            }
        }
    }
    // only one point stored (left element of the root)
    else if (size_leaf_ == 1)
    {
        return root_2->leafLeft_2;
    }
    else // no point stored
    {
        SharedPointer<Leaf> ret;
        ret.offset = -1;
        return ret;
    }
}

template <typename T, typename T2>
Foam::SharedPointer<typename Foam::parallelISAT<T, T2>::Leaf> Foam::parallelISAT<T, T2>::binaryTreeSearch(
    const inputType &value,
    SharedPointer<Node> node_ptr)
{

    if (size_leaf_ > 1)
    {
        if (node_ptr->goLeft(value))
        {
            if (!node_ptr->nodeLeft_2.isNULL())
            {
                return binaryTreeSearch(value, node_ptr->nodeLeft_2);
            }
            else
            {
                return node_ptr->leafLeft_2;
            }
        }
        else
        {
            if (!node_ptr->nodeRight_2.isNULL())
            {
                return binaryTreeSearch(value, node_ptr->nodeRight_2);
            }
            else
            {
                return node_ptr->leafRight_2;
            }
        }
    }
    // only one point stored (left element of the root)
    else if (size_leaf_ == 1)
    {
        return root_2->leafLeft_2;
    }
    else // no point stored
    {
        SharedPointer<Leaf> ret;
        ret.offset = -1;
        return ret;
    }
}
template <typename T, typename T2>
bool Foam::parallelISAT<T, T2>::grow(const inputType &x)
{
    bool success = false;
    SharedPointer<Leaf> &pleaf = binaryTreeSearch(x);
    outputType y_approx;
    pleaf->retrieve(x, y_approx);
    if (!pleaf->inEOA(x))
    {

        if (pleaf->checkSolution(x, y_approx))
        {
            SharedPointer<Leaf> &pleaf_new = leaf_manager.New();

            if (!pleaf_new.isNull())
            {
                *pleaf_new = *pleaf;
                pleaf_new->grow(x);

                if (pleaf->node_2->leafLeft_2.offset == pleaf.offset)
                {
                    pleaf->node_2->leafLeft_2.offset = pleaf_new.offset;
                }
                else
                {
                    pleaf->node_2->leafRight_2.offset = pleaf_new.offset;
                }
                leaf_manager.Delete(pleaf);
                success = true;
            }
        }
    }
    return success;
}
template <typename T, typename T2>
bool Foam::parallelISAT<T, T2>::grow(const scalarField &x, SharedPointer<Leaf> &pleaf)
{

    bool success = false;

    SharedPointer<Leaf> pleaf_new = leaf_manager.New(*pleaf);

    if (!pleaf_new.isNULL())
    {

        pleaf_new->node_2.offset = pleaf->node_2.offset;
        success = pleaf_new->grow(x);
        if (success)
        {
            if (pleaf->node_2->leafLeft_2.offset == pleaf.offset)
            {
                pleaf->node_2->leafLeft_2.offset = pleaf_new.offset;
            }
            else
            {
                pleaf->node_2->leafRight_2.offset = pleaf_new.offset;
            }
            leaf_manager.Delete(pleaf);
        }
        else
        {
            leaf_manager.Delete(pleaf_new);
        }
    }
    return success;
}

template <typename T, typename T2>
bool Foam::parallelISAT<T, T2>::retrieve(const scalarField &x, scalarField &y)
{

    SharedPointer<Leaf> pleaf = binaryTreeSearch(x);
    if (pleaf.isNULL())
    {
        return false;
    }
    if (pleaf->inEOA(x))
    {
        pleaf->retrieve(x, y);
        nRetrieved_++;
        return true;
    }
    else
    {
        return false;
    }
}

template <typename T, typename T2>
bool Foam::parallelISAT<T, T2>::isFull()
{
    return size_leaf_ >= maxNLeafs_;
}

template <typename T, typename T2>
template <typename... Args>
Foam::label Foam::parallelISAT<T, T2>::add(
    const scalarField &x, const scalarField &y, Args... args)
{
    write_lock.lock();

    label growthOrAddFlag = 1;

    if (isFull())
    {
        write_lock.unlock();
        return growthOrAddFlag;
    }

    if (size_leaf_ == 0)
    {
        SharedPointer<Node> pnode = node_manager.New(x.size());

        if (pnode.isNULL())
        {
            growthOrAddFlag = 0;

            write_lock.unlock();
            return growthOrAddFlag;
        }
        SharedPointer<Leaf> pleaf = leaf_manager.New(x.size());

        if (pleaf.isNULL())
        {
            growthOrAddFlag = 0;
            node_manager.Delete(pnode);

            write_lock.unlock();
            return growthOrAddFlag;
        }

        pnode->parent_2.offset = -1;
        pnode->nodeLeft_2.offset = -1;
        pnode->nodeRight_2.offset = -1;
        pnode->leafLeft_2.offset = pleaf.offset;
        pnode->leafRight_2.offset = -1;
        pnode->set();

        pleaf->node_2.offset = pnode.offset;
        pleaf->set(x, y, args...);
        root_2.offset = pnode.offset;
        size_leaf_++;
        size_node++;
        nAdd_++;

        write_lock.unlock();
        return growthOrAddFlag;
    }
    SharedPointer<Leaf> pleaf = binaryTreeSearch(x);

    if (pleaf->inEOA(x))
    {
        growthOrAddFlag = 0;
        write_lock.unlock();

        return growthOrAddFlag;
    }
    if (pleaf->checkSolution(x, y))
    {
        if (grow(x, pleaf))
        {
            nGrowth_++;
            totalGrowth_++;

            growthOrAddFlag = 0;

            write_lock.unlock();
            return growthOrAddFlag;
        }
    }

    //Todo  leaf_index can be obtain by arg
    SharedPointer<Leaf> pleaf2 = binaryTreeSearch(x, root_2);

    if (pleaf2->inEOA(x))
    {
        growthOrAddFlag = 0;
        write_lock.unlock();
        return growthOrAddFlag;
    }
    SharedPointer<Node> parentNode = pleaf2->node_2;

    SharedPointer<Node> pnode_new = node_manager.New(x.size());

    if (pnode_new.isNULL())
    {
        growthOrAddFlag = 0;
        write_lock.unlock();
        return growthOrAddFlag;
    }
    SharedPointer<Leaf> pleaf_new = leaf_manager.New(x.size());

    if (pleaf_new.isNULL())
    {
        growthOrAddFlag = 0;

        node_manager.Delete(pnode_new);

        write_lock.unlock();
        return growthOrAddFlag;
    }

    pleaf_new->set(x, y, args...);
    pleaf_new->node_2.offset = pnode_new.offset;
    pleaf2->node_2.offset = pnode_new.offset;

    pnode_new->parent_2.offset = -1;
    pnode_new->nodeLeft_2.offset = -1;
    pnode_new->nodeRight_2.offset = -1;
    pnode_new->set(*pleaf_new, *pleaf2); //v = (leafmem[new_leaf_index].v + leafmem[leaf_index].v) / 2;

    if (pnode_new->goLeft(*pleaf_new))
    {
        pnode_new->leafLeft_2.offset = pleaf_new.offset;
        pnode_new->leafRight_2.offset = pleaf2.offset;
    }
    else
    {
        pnode_new->leafLeft_2.offset = pleaf2.offset;
        pnode_new->leafRight_2.offset = pleaf_new.offset;
    }

    if (size_leaf_ > 1)
    {
        pnode_new->parent_2.offset = parentNode.offset;
        if (parentNode->leafLeft_2.offset == pleaf2.offset)
        {
            parentNode->nodeLeft_2.offset = pnode_new.offset;
            parentNode->leafLeft_2.offset = -1;
        }
        else
        {
            parentNode->nodeRight_2.offset = pnode_new.offset;
            parentNode->leafRight_2.offset = -1;
        }
        size_leaf_++;
        size_node++;
        nAdd_++;

        write_lock.unlock();
        return growthOrAddFlag;
    }
    else
    {
        root_2.offset = pnode_new.offset;

        node_manager.Delete(parentNode);

        size_leaf_++;
        nAdd_++;

        write_lock.unlock();

        return growthOrAddFlag;
    }
}

template <typename T, typename T2>
Foam::SharedPointer<typename Foam::parallelISAT<T, T2>::Leaf> Foam::parallelISAT<T, T2>::treeMin(SharedPointer<Node> subTreeRoot)
{

    if (!subTreeRoot.isNULL())
    {
        SharedPointer<Node> ret;
        ret.offset = subTreeRoot.offset;
        while (!ret->nodeLeft_2.isNULL())
        {
            ret.offset = ret->nodeLeft_2.offset;
        }
        return ret->leafLeft_2;
    }
    else
    {
        SharedPointer<Leaf> ret;
        ret.offset = -1;
        return ret;
    }
}

template <typename T, typename T2>
Foam::SharedPointer<typename Foam::parallelISAT<T, T2>::Leaf> Foam::parallelISAT<T, T2>::treeSuccessor(SharedPointer<Leaf> x)
{
    SharedPointer<Leaf> ret;
    ret.offset = -1;
    if (size_leaf_ > 1)
    {

        if (x.offset == x->node_2->leafLeft_2.offset)
        {
            if (x->node_2->nodeRight_2.offset == -1)
            {
                return x->node_2->leafRight_2;
            }
            else
            {
                return treeMin(x->node_2->nodeRight_2);
            }
        }
        else if (x.offset == x->node_2->leafRight_2.offset)
        {
            SharedPointer<Node> y;
            y.offset = x->node_2.offset;
            while (!y->parent_2.isNULL())
            {
                if (y.offset == y->parent_2->nodeLeft_2.offset)
                {
                    if (y->parent_2->nodeRight_2.isNULL())
                    {
                        return y->parent_2->leafRight_2;
                    }
                    else
                    {
                        return treeMin(y->parent_2->nodeRight_2);
                    }
                }
                y.offset = y->parent_2.offset;
            }
            // when we reach this point, y points to the root and
            // never entered in the if loop (coming from the right)
            // so we are at the tree maximum and there is no successor
            return ret;
        }
        else
        {
            FatalErrorInFunction
                << "inconsistent structure of the tree, no leaf and no node"
                << exit(FatalError);
            return ret;
        }
    }

    return ret;
}
template <typename T, typename T2>
Foam::SharedPointer<typename Foam::parallelISAT<T, T2>::Leaf> Foam::parallelISAT<T, T2>::leafPSibling(SharedPointer<Leaf> x)
{
    if (size_leaf_ > 1)
    {
        if (x.offset == x->node_2->leafLeft_2.offset)
        {
            // x is on the left, return right side
            // might return nullptr if the right side is a node
            return x->node_2->leafRight_2;
        }
        else if (x.offset == x->node_2->leafRight_2.offset)
        {
            // x is on the right, return left side
            return x->node_2->leafLeft_2;
        }
        else
        {
            FatalErrorInFunction
                << "wrong addressing of the initial leaf"
                << exit(FatalError);
            SharedPointer<Leaf> ret;
            ret.offset = -1;
            return ret;
        }
    }
    // there is only one leaf attached to the root_, no sibling
    SharedPointer<Leaf> ret;
    ret.offset = -1;
    return ret;
}

template <typename T, typename T2>
Foam::SharedPointer<typename Foam::parallelISAT<T, T2>::Node> Foam::parallelISAT<T, T2>::nodeSibling_node(SharedPointer<Node> y)
{
    if (!y->parent_2.isNULL())
    {
        if (y.offset == y->parent_2->nodeLeft_2.offset)
        {
            // y is on the left, return right side
            return y->parent_2->nodeRight_2;
        }
        else if (y.offset == y->parent_2->nodeRight_2.offset)
        {
            return y->parent_2->nodeLeft_2;
        }
        else
        {
            FatalErrorInFunction
                << "wrong addressing of the initial node"
                << exit(FatalError);
            SharedPointer<Node> ret;
            ret.offset = -1;
            return ret;
        }
    }
    SharedPointer<Node> ret;
    ret.offset = -1;
    return ret;
}

template <typename T, typename T2>
Foam::SharedPointer<typename Foam::parallelISAT<T, T2>::Node> Foam::parallelISAT<T, T2>::nodeSibling_leaf(SharedPointer<Leaf> x)
{
    if (size_leaf_ > 1)
    {
        if (x.offset == x->node_2->leafLeft_2.offset)
        {
            // x is on the left, return right side
            return x->node_2->nodeRight_2;
        }
        else if (x.offset == x->node_2->leafRight_2.offset)
        {
            // x is on the right, return left side
            return x->node_2->nodeLeft_2;
        }
        else
        {
            FatalErrorInFunction
                << "wrong addressing of the initial leaf"
                << exit(FatalError);

            SharedPointer<Node> ret;
            ret.offset = -1;
            return ret;
        }
    }
    SharedPointer<Node> ret;
    ret.offset = -1;
    return ret;
}

template <typename T, typename T2>
void Foam::parallelISAT<T, T2>::transplant(SharedPointer<Node> u, SharedPointer<Node> v)
{
    if (v != -1)
    {
        // u is root_
        if (u->parent_2.isNULL())
        {
            root_2.offset = v.offset;
        }
        // u is on the left of its parent
        else if (u.offset == u->parent_2->nodeLeft_2.offset)
        {
            u->parent_2->nodeLeft_2.offset = v.offset;
        }
        // u is ont the right of its parent
        else if (u.offset == u->parent_2->nodeRight_2.offset)
        {
            u->parent_2->nodeRight_2.offset = v.offset;
        }
        else
        {
            FatalErrorInFunction
                << "wrong addressing of the initial node"
                << exit(FatalError);
        }
        v->parent_2.offset = u->parent_2.offset;
    }
    else
    {
        FatalErrorInFunction
            << "trying to transplant a nullptr node"
            << exit(FatalError);
    }
}

template <typename T, typename T2>
void Foam::parallelISAT<T, T2>::deleteLeaf(SharedPointer<Leaf> index_phi0)
{
    if (size_leaf_ == 1) // only one point is stored
    {
        leaf_manager.Delete(index_phi0);
        node_manager.Delete(root_2);
        root_2.offset = -1;
        size_node--;
        //deleteDemandDrivenData(phi0);
        //deleteDemandDrivenData(root_);
    }
    else if (size_leaf_ > 1)
    {
        SharedPointer<Node> z = index_phi0->node_2;
        label x;
        SharedPointer<Leaf> siblingPhi0 = leafPSibling(index_phi0);

        if (!siblingPhi0.isNULL()) // the sibling of phi0 is a chemPoint
        {
            // z was root (only two chemPoints in the tree)
            if (z->parent_2.isNULL())
            {
                root_2 = node_manager.New(siblingPhi0->phi().size_);

                if (root_2.isNULL())
                {
                    FatalErrorInFunction
                        << "run out of memory"
                        << exit(FatalError);
                }
                size_node++;
                root_2->nodeRight_2.offset = -1;
                root_2->leafRight_2.offset = -1;
                root_2->nodeLeft_2.offset = -1;
                root_2->leafLeft_2.offset = siblingPhi0.offset;
                root_2->parent_2.offset = -1;
                siblingPhi0->node_2.offset = root_2.offset;
            }
            else if (z.offset == z->parent_2->nodeLeft_2.offset)
            {
                z->parent_2->leafLeft_2.offset = siblingPhi0.offset;
                z->parent_2->nodeLeft_2.offset = -1;
                siblingPhi0->node_2.offset = z->parent_2.offset;

            }
            else if (z.offset == z->parent_2->nodeRight_2.offset)
            {
                z->parent_2->leafRight_2.offset = siblingPhi0.offset;
                z->parent_2->nodeRight_2.offset = -1;
                siblingPhi0->node_2.offset = z->parent_2.offset;

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
            SharedPointer<Node> x = nodeSibling_leaf(index_phi0);
            if (!x.isNULL())
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

        leaf_manager.Delete(index_phi0);
        node_manager.Delete(z);
        size_node--;
        //deleteDemandDrivenData(phi0);
        //deleteDemandDrivenData(z);
    }
    size_leaf_--;
}

template <typename T, typename T2>
void Foam::parallelISAT<T, T2>::deleteAllNode(SharedPointer<Node> subTreeRoot)
{
    if (!subTreeRoot.isNULL())
    {
        deleteAllNode(subTreeRoot->nodeLeft_2);
        deleteAllNode(subTreeRoot->nodeRight_2);
        node_manager.Delete(subTreeRoot);
        size_node--;
        //deleteDemandDrivenData(subTreeRoot);
    }
}

template <typename T, typename T2>
void Foam::parallelISAT<T, T2>::deleteAll(SharedPointer<Node> subTreeRoot)
{
    if (!subTreeRoot.isNULL())
    {
        deleteAll(subTreeRoot->nodeLeft_2);
        deleteAll(subTreeRoot->nodeRight_2);
        if (!subTreeRoot->leafLeft_2.isNULL())
        {
            leaf_manager.Delete(subTreeRoot->leafLeft_2);
            size_leaf_--;
        }

        if (!subTreeRoot->leafRight_2.isNULL())
        {
            leaf_manager.Delete(subTreeRoot->leafRight_2);
            size_leaf_--;
        }
        node_manager.Delete(subTreeRoot);
        size_node--;
        //deleteDemandDrivenData(subTreeRoot);
    }
}

template <typename T, typename T2>
Foam::label Foam::parallelISAT<T, T2>::depth(SharedPointer<Node> subTreeRoot)
{
    // when we reach the leaf, we return 0
    if (subTreeRoot.isNULL())
    {
        return 0;
    }
    else
    {
        return 1 + max(
                       depth(subTreeRoot->nodeLeft_2),
                       depth(subTreeRoot->nodeRight_2));
    }
}
template <typename T, typename T2>
bool Foam::parallelISAT<T, T2>::treeCheck(SharedPointer<Node> index)
{
    bool leftnode = (!index->nodeLeft_2.isNULL());
    bool leftleaf = (!index->leafLeft_2.isNULL());
    if ((leftnode && leftleaf) || ((!leftnode) && (!leftleaf)))
    {
        FatalErrorInFunction
            << "p4 ";
        //<< abort(FatalError);
        return false;
    }
    bool rightnode = (!index->nodeRight_2.isNULL());
    bool rightleaf = (!index->leafRight_2.isNULL());
    if ((rightnode && rightleaf) || ((!rightnode) && (!rightleaf)))
    {
        FatalErrorInFunction
            << "p5 ";
        //<< abort(FatalError);
        return false;
    }
    if (leftnode)
    {
        if (index->nodeLeft_2->parent_2.offset != index.offset)
        {
            FatalErrorInFunction
                << "p6 ";
            //<< abort(FatalError);

            return false;
        }
    }
    else
    {
        if (index->leafLeft_2->node_2.offset != index.offset)

        {
            FatalErrorInFunction
                << "p7 ";
            //<< abort(FatalError);
            return false;
        }
    }
    if (rightnode)
    {
        if (index->nodeRight_2->parent_2.offset != index.offset)
        {
            FatalErrorInFunction
                << "p8 ";
            //<< abort(FatalError);
            return false;
        }
    }
    else
    {
        if (index->leafRight_2->node_2.offset != index.offset)
        {
            FatalErrorInFunction
                << "p9 ";
            //<< abort(FatalError);
            return false;
        }
    }
    if (leftnode)
    {
        if (treeCheck(index->nodeLeft_) == false)
        {
            FatalErrorInFunction
                << "p10 ";
            //<< abort(FatalError);
            return false;
        }
    }
    if (rightnode)
    {
        if (treeCheck(index->nodeRight_) == false)
        {
            FatalErrorInFunction
                << "p11 ";
            //<< abort(FatalError);
            return false;
        }
    }
    return true;
}

template <typename T, typename T2>
void Foam::parallelISAT<T, T2>::insertNode(
    SharedPointer<Leaf> phi0,
    SharedPointer<Node> newNode)
{
    if (phi0.offset == phi0->node_2->leafRight_2.offset) // phi0 is on the right
    {
        phi0->node_2->leafRight_2.offset = -1;
        phi0->node_2->nodeRight_2.offset = newNode.offset;

        return;
    }
    else if (phi0.offset == phi0->node_2->leafLeft_2.offset) // phi0 is on the left
    {
        phi0->node_2->leafLeft_2.offset = -1;
        phi0->node_2->nodeLeft_2.offset = newNode.offset;

        return;
    }

    // if we reach this point, there is an issue with the addressing
    FatalErrorInFunction
        << "trying to insert a node with a wrong pointer to a chemPoint"
        << exit(FatalError);
}

template <typename T, typename T2>
void Foam::parallelISAT<T, T2>::valid()
{
    SharedPointer<Leaf> x = this->treeMin();
    List<SharedPointer<Leaf>> chemPoints(this->size_leaf_);

    label chPi = 0;

    while (x != -1)
    {
        if (x->phi().size_ < 1)
        {
            FatalErrorInFunction
                << ",size=" << x->phi().size_
                << abort(FatalError);
        }
        chemPoints[chPi++].offset = x.offset;
        //const typename DataType::inputType &phij = x->phi();
        x = this->treeSuccessor(x);
    }
    if (chPi != this->size_leaf_)
        FatalErrorInFunction
            << "size=" << chPi << " " << this->size_leaf_
            << abort(FatalError);

}
/*
template <typename T, typename T2>
void Foam::parallelISAT<T, T2>::writeDot(string name)
{
    std::ofstream fout(name);
    Pout << "writeDot,root_=" << root_ << endl;
    fout << "digraph G{" << std::endl;
    if (root_ != -1)
    {
        fout << "node" << root_ << "[label=\"" << nodemem[root_] << "\"];" << std::endl;
        writeDot(root_, fout);
    }
    fout << "}" << std::endl;
}

template <typename T, typename T2>
void Foam::parallelISAT<T, T2>::writeDot(label nodeindex, std::ofstream &fout)
{
    bool flag_left = false, flag_right = false;
    if (nodemem[nodeindex].nodeLeft_ != -1)
    {
        fout << "node" << nodemem[nodeindex].nodeLeft_ << "[label=\"" << nodemem[nodemem[nodeindex].nodeLeft_] << "\"];" << std::endl;
        fout << "node" << nodeindex << " -> "
             << "node" << nodemem[nodeindex].nodeLeft_ << ";" << std::endl;
        flag_left = true;
    }
    if (nodemem[nodeindex].leafLeft_ != -1)
    {
        fout << "leaf" << nodemem[nodeindex].leafLeft_ << "[label=\"" << leafmem[nodemem[nodeindex].leafLeft_] << "\"];" << std::endl;
        fout << "node" << nodeindex << " -> "
             << "leaf" << nodemem[nodeindex].leafLeft_ << ";" << std::endl;
    }
    fout << "M" << nodeindex << "[style=invis];" << std::endl;
    fout << "node" << nodeindex << " -> "
         << "M" << nodeindex << "[style=invis];" << std::endl;
    if (nodemem[nodeindex].nodeRight_ != -1)
    {
        fout << "node" << nodemem[nodeindex].nodeRight_ << "[label=\"" << nodemem[nodemem[nodeindex].nodeRight_] << "\"];" << std::endl;
        fout << "node" << nodeindex << " -> "
             << "node" << nodemem[nodeindex].nodeRight_ << ";" << std::endl;
        flag_right = true;
    }
    if (nodemem[nodeindex].leafRight_ != -1)
    {
        fout << "leaf" << nodemem[nodeindex].leafRight_ << "[label=\"" << leafmem[nodemem[nodeindex].leafRight_] << "\"];" << std::endl;
        fout << "node" << nodeindex << " -> "
             << "leaf" << nodemem[nodeindex].leafRight_ << ";" << std::endl;
    }
    if (flag_left)
        writeDot(nodemem[nodeindex].nodeLeft_, fout);
    if (flag_right)
        writeDot(nodemem[nodeindex].nodeRight_, fout);
    return;
}
*/
// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //
/*
void Foam::ISAT::operator=(const ISAT& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorInFunction
            << "Attempted assignment to self"
            << abort(FatalError);
    }
}
*/
// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

// ************************************************************************* //
