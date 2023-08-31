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
#define LIKELY(exp) __builtin_expect(exp, 1)
#define UNLIKELY(exp) __builtin_expect(exp, 0)
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
      emptylist_leaf(manager_in, 100000), emptylist_node(manager_in, 100000), pv_(manager_in),
      size_leaf_(pv_().size_leaf), head_leaf(pv_().head_leaf),
      tail_leaf(pv_().tail_leaf), temp_tail_leaf(pv_().temp_tail_leaf),
      size_node(pv_().size_node), head_node(pv_().heaf_node),
      tail_node(pv_().tail_node), temp_tail_node(pv_().temp_tail_node),
      root_(pv_().root_), totalGrowth_(pv_().totalGrowth), //,root_atomic((std::atomic<long unsigned int>*)&pv_().root_.offset)
      write_lock(manager_in, sync_in), mem_lock(manager_in, sync_in),
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

    root_.setNULL();
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
    auto ret = binaryTreeSearch(value, root_);
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
            if (node_ptr->nodeLeft_.notNULL())
            {
                return binaryTreeSearch(value, node_ptr->nodeLeft_);
            }
            else
            {
                return node_ptr->leafLeft_;
            }
        }
        else
        {
            if (node_ptr->nodeRight_.notNULL())
            {
                return binaryTreeSearch(value, node_ptr->nodeRight_);
            }
            else
            {
                return node_ptr->leafRight_;
            }
        }
    }
    // only one point stored (left element of the root)
    else if (size_leaf_ == 1)
    {
        return root_->leafLeft_;
    }
    else // no point stored
    {
        SharedPointer<Leaf> ret;
        ret.setNULL();
        return ret;
    }
}

template <typename T, typename T2>
Foam::SharedPointer<typename Foam::parallelISAT<T, T2>::Leaf> Foam::parallelISAT<T, T2>::binaryTreeSearch(
    const scalarField &value,
    SharedPointer<Node> node_ptr, SharedPointer<Node> &node_ret, int &side)
{
    if (size_leaf_ > 1)
    {

        if (node_ptr->goLeft(value))
        {
            if (node_ptr->nodeLeft_.notNULL())
            {
                return binaryTreeSearch(value, node_ptr->nodeLeft_, node_ret, side);
            }
            else
            {
                node_ret = node_ptr;
                side = 0; //0 left, 1 right

                return node_ptr->leafLeft_;
            }
        }
        else
        {
            if (node_ptr->nodeRight_.notNULL())
            {
                return binaryTreeSearch(value, node_ptr->nodeRight_, node_ret, side);
            }
            else
            {
                node_ret = node_ptr;
                side = 1; //0 left, 1 right
                return node_ptr->leafRight_;
            }
        }
    }
    // only one point stored (left element of the root)
    else if (size_leaf_ == 1)
    {
        node_ret = root_;
        side = 0; //0 left, 1 right
        return root_->leafLeft_;
    }
    else // no point stored
    {
        SharedPointer<Leaf> ret;
        ret.setNULL();
        node_ret.setNULL();
        side = 0;
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
            if (node_ptr->nodeLeft_.notNULL())
            {
                return binaryTreeSearch(value, node_ptr->nodeLeft_);
            }
            else
            {
                return node_ptr->leafLeft_;
            }
        }
        else
        {
            if (node_ptr->nodeRight_.notNULL())
            {
                return binaryTreeSearch(value, node_ptr->nodeRight_);
            }
            else
            {
                return node_ptr->leafRight_;
            }
        }
    }
    // only one point stored (left element of the root)
    else if (size_leaf_ == 1)
    {
        return root_->leafLeft_;
    }
    else // no point stored
    {
        SharedPointer<Leaf> ret;
        ret.setNULL();
        return ret;
    }
}
/* template <typename T, typename T2>
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

                if (pleaf->node_2->leafLeft_ == pleaf)
                {
                    pleaf->node_2->leafLeft_ = pleaf_new;
                }
                else
                {
                    pleaf->node_2->leafRight_ = pleaf_new;
                }
                leaf_manager.Delete(pleaf);
                success = true;
            }
        }
    }
    return success;
} */
template <typename T, typename T2>
bool Foam::parallelISAT<T, T2>::grow(const scalarField &x, SharedPointer<Leaf> &pleaf)
{

    bool success = false;
    mem_lock.lock();
    SharedPointer<Leaf> pleaf_new = leaf_manager.New(*pleaf);
    mem_lock.unlock();

    if (pleaf_new.notNULL())
    {

        pleaf_new->node_2 = pleaf->node_2;
        success = pleaf_new->grow(x, mem_lock);
        if (success)
        {
            SharedPointer<Leaf> *leaf;
            if (pleaf->node_2->leafLeft_ == pleaf)
            {
                leaf = &pleaf->node_2->leafLeft_;
                //pleaf->node_2->leafLeft_ = pleaf_new;
            }
            else
            {
                leaf = &pleaf->node_2->leafRight_;
                //pleaf->node_2->leafRight_ = pleaf_new;
            }

            std::atomic<size_t> *leaf_atomic = reinterpret_cast<std::atomic<long unsigned int> *>(&leaf->offset);
            //std::atomic<size_t> *leaf_atomic = (std::atomic<long unsigned int> *)(&leaf->offset);
            leaf_atomic->store(pleaf_new.offset, std::memory_order_release);
            mem_lock.lock();
            leaf_manager.Delete(pleaf);
            mem_lock.unlock();
        }
        else
        {
            mem_lock.lock();
            leaf_manager.Delete(pleaf_new);
            mem_lock.unlock();
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
    //write_lock.lock();
    std::atomic<size_t> *root_atomic = reinterpret_cast<std::atomic<long unsigned int> *>(&root_.offset);
    label growthOrAddFlag = 1;

    if (isFull())
    {
        //write_lock.unlock();
        return growthOrAddFlag;
    }

    if (size_leaf_ == 0)
    {
        write_lock.lock();
        if (size_leaf_ == 0)
        {
            mem_lock.lock();
            SharedPointer<Node> pnode = node_manager.New(x.size());
            SharedPointer<Leaf> pleaf = leaf_manager.New(x.size());
            mem_lock.unlock();
            /*             if (UNLIKELY(pnode.isNULL()))
            {
                growthOrAddFlag = 0;

                write_lock.unlock();
                return growthOrAddFlag;
            }
            
            if (UNLIKELY(pleaf.isNULL()))
            {
                growthOrAddFlag = 0;
                node_manager.Delete(pnode);

                write_lock.unlock();
                return growthOrAddFlag;
            } */
            //Pout << "\n???here!"<< "," << size_leaf_ << "\n"<< endl;
            pnode->parent_.setNULL();
            pnode->nodeLeft_.setNULL();
            pnode->nodeRight_.setNULL();
            pnode->leafLeft_ = pleaf;
            pnode->leafRight_.setNULL();
            pnode->set();

            pleaf->node_2 = pnode;
            pleaf->set(x, y, mem_lock, args...);

            //root_ = pnode;
            //((std::atomic<size_t>)root_.offset).store(pnode.offset, std::memory_order_release);

            root_atomic->store(pnode.offset, std::memory_order_release);
            size_leaf_++;
            size_node++;
            nAdd_++;

            write_lock.unlock();
            return growthOrAddFlag;
        }
        else
        {
            write_lock.unlock();
        }
    }

    /*     SharedPointer<Leaf> pleaf = binaryTreeSearch(x);

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
    }    */

    //Todo  leaf_index can be obtain by arg
    SharedPointer<Node> parentNode;
    int ret_side; // 0 left, 1 right
    SharedPointer<Leaf> pleaf2 = binaryTreeSearch(x, SharedPointer<Node>(root_atomic->load(std::memory_order_acquire)), parentNode, ret_side);
    //Pout << "\n!!!!!!!!!!!" << pleaf2->phi().size_ << "," << size_leaf_ << "\n"<< endl;

    if (pleaf2.isNULL() || pleaf2->inEOA(x))
    {
        growthOrAddFlag = 0;
        //write_lock.unlock();
        return growthOrAddFlag;
    }
    //SharedPointer<Node> parentNode = pleaf2->node_2;
    if (pleaf2->checkSolution(x, y))
    {
        if ((*parentNode).mutex.try_lock())
        {
            bool tmp_flag = false;
            if (ret_side == 0)
            {
                tmp_flag = (parentNode->leafLeft_ == pleaf2 && pleaf2->node_2 == parentNode);
            }
            else
            {
                tmp_flag = (parentNode->leafRight_ == pleaf2 && pleaf2->node_2 == parentNode);
            }
            if (tmp_flag)
            {
                if (grow(x, pleaf2))
                {
                    nGrowth_++;
                    totalGrowth_++;

                    growthOrAddFlag = 0;
                    (*parentNode).mutex.unlock();
                    //write_lock.unlock();
                    return growthOrAddFlag;
                }
            }
            (*parentNode).mutex.unlock();
        }
    }

    //

    mem_lock.lock();
    SharedPointer<Node> pnode_new = node_manager.New(x.size());
    SharedPointer<Leaf> pleaf_new = leaf_manager.New(x.size());
    mem_lock.unlock();

    //write_lock.lock();

    /*     if (pnode_new.isNULL())
    {
        growthOrAddFlag = 0;
        write_lock.unlock();
        return growthOrAddFlag;
    }
    
    if (pleaf_new.isNULL())
    {
        growthOrAddFlag = 0;

        node_manager.Delete(pnode_new);

        write_lock.unlock();
        return growthOrAddFlag;
    } */

    pleaf_new->set(x, y, mem_lock, args...);
    pleaf_new->node_2 = pnode_new;

    pnode_new->parent_.setNULL();
    pnode_new->nodeLeft_.setNULL();
    pnode_new->nodeRight_.setNULL();
    pnode_new->set(*pleaf_new, *pleaf2, mem_lock); //v = (leafmem[new_leaf_index].v + leafmem[leaf_index].v) / 2;

    if (!(*parentNode).mutex.try_lock())
    {
        mem_lock.lock();
        node_manager.Delete(pnode_new);
        leaf_manager.Delete(pleaf_new);
        mem_lock.unlock();
        growthOrAddFlag = 0;
        return growthOrAddFlag;
    }
    bool tmp_flag = false;
    if (ret_side == 0)
    {
        tmp_flag = (parentNode->leafLeft_ == pleaf2 && pleaf2->node_2 == parentNode);
    }
    else
    {
        tmp_flag = (parentNode->leafRight_ == pleaf2 && pleaf2->node_2 == parentNode);
    }
    if (!tmp_flag)
    {
        (*parentNode).mutex.unlock();
        growthOrAddFlag = 0;
        return growthOrAddFlag;
    }

    pleaf2->node_2 = pnode_new;

    if (pnode_new->goLeft(*pleaf_new))
    {
        pnode_new->leafLeft_ = pleaf_new;
        pnode_new->leafRight_ = pleaf2;
    }
    else
    {
        pnode_new->leafLeft_ = pleaf2;
        pnode_new->leafRight_ = pleaf_new;
    }

    if (size_leaf_ > 1)
    {
        pnode_new->parent_ = parentNode;
        SharedPointer<Leaf> *leaf;
        SharedPointer<Node> *node;
        if (parentNode->leafLeft_ == pleaf2)
        {
            //parentNode->nodeLeft_ = pnode_new;
            // parentNode->leafLeft_.setNULL();
            leaf = &parentNode->leafLeft_;
            node = &parentNode->nodeLeft_;
        }
        else
        {
            //parentNode->nodeRight_ = pnode_new;
            //parentNode->leafRight_.setNULL();
            leaf = &parentNode->leafRight_;
            node = &parentNode->nodeRight_;
        }
        //node->offset = pnode_new.offset;
        //leaf->offset = sptr_NULL;

        std::atomic<size_t> *node_atomic = reinterpret_cast<std::atomic<long unsigned int> *>(&node->offset);
        std::atomic<size_t> *leaf_atomic = reinterpret_cast<std::atomic<long unsigned int> *>(&leaf->offset);

        //std::atomic<size_t> *node_atomic = (std::atomic<long unsigned int> *)(&node->offset);
        //std::atomic<size_t> *leaf_atomic = (std::atomic<long unsigned int> *)(&leaf->offset);

        node_atomic->store(pnode_new.offset, std::memory_order_release);
        leaf_atomic->store(sptr_NULL, std::memory_order_release);

        //((std::atomic<size_t>)node->offset).store(pnode_new.offset, std::memory_order_release);
        //((std::atomic<size_t>)leaf->offset).store(sptr_NULL, std::memory_order_release);

        size_leaf_++;
        size_node++;
        nAdd_++;
        (*parentNode).mutex.unlock();
        //write_lock.unlock();
        return growthOrAddFlag;
    }
    else
    {

        std::atomic<size_t> *root_atomic = reinterpret_cast<std::atomic<long unsigned int> *>(&root_.offset);
        //std::atomic<size_t> *root_atomic = (std::atomic<long unsigned int> *)(&root_.offset);

        root_atomic->store(pnode_new.offset, std::memory_order_release);
        //root_ = pnode_new;
        //((std::atomic<size_t>)root_.offset).store(pnode_new.offset, std::memory_order_release);
        mem_lock.lock();
        node_manager.Delete(parentNode);
        mem_lock.unlock();
        size_leaf_++;
        nAdd_++;

        //write_lock.unlock();
        //(*parentNode).mutex.unlock();
        return growthOrAddFlag;
    }
}

template <typename T, typename T2>
Foam::SharedPointer<typename Foam::parallelISAT<T, T2>::Leaf> Foam::parallelISAT<T, T2>::treeMin(SharedPointer<Node> subTreeRoot)
{

    if (subTreeRoot.notNULL())
    {
        SharedPointer<Node> ret;
        ret = subTreeRoot;
        while (ret->nodeLeft_.notNULL())
        {
            ret = ret->nodeLeft_;
        }
        return ret->leafLeft_;
    }
    else
    {
        SharedPointer<Leaf> ret;
        ret.setNULL();
        return ret;
    }
}

template <typename T, typename T2>
Foam::SharedPointer<typename Foam::parallelISAT<T, T2>::Leaf> Foam::parallelISAT<T, T2>::treeSuccessor(SharedPointer<Leaf> x)
{
    SharedPointer<Leaf> ret;
    ret.setNULL();
    if (size_leaf_ > 1)
    {

        if (x == x->node_2->leafLeft_)
        {
            if (x->node_2->nodeRight_.isNULL())
            {
                return x->node_2->leafRight_;
            }
            else
            {
                return treeMin(x->node_2->nodeRight_);
            }
        }
        else if (x == x->node_2->leafRight_)
        {
            SharedPointer<Node> y;
            y = x->node_2;
            while (y->parent_.notNULL())
            {
                if (y == y->parent_->nodeLeft_)
                {
                    if (y->parent_->nodeRight_.isNULL())
                    {
                        return y->parent_->leafRight_;
                    }
                    else
                    {
                        return treeMin(y->parent_->nodeRight_);
                    }
                }
                y = y->parent_;
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
        if (x == x->node_2->leafLeft_)
        {
            // x is on the left, return right side
            // might return nullptr if the right side is a node
            return x->node_2->leafRight_;
        }
        else if (x == x->node_2->leafRight_)
        {
            // x is on the right, return left side
            return x->node_2->leafLeft_;
        }
        else
        {
            FatalErrorInFunction
                << "wrong addressing of the initial leaf"
                << exit(FatalError);
            SharedPointer<Leaf> ret;
            ret.setNULL();
            return ret;
        }
    }
    // there is only one leaf attached to the root_, no sibling
    SharedPointer<Leaf> ret;
    ret.setNULL();
    return ret;
}

template <typename T, typename T2>
Foam::SharedPointer<typename Foam::parallelISAT<T, T2>::Node> Foam::parallelISAT<T, T2>::nodeSibling_node(SharedPointer<Node> y)
{
    if (y->parent_.notNULL())
    {
        if (y == y->parent_->nodeLeft_)
        {
            // y is on the left, return right side
            return y->parent_->nodeRight_;
        }
        else if (y == y->parent_->nodeRight_)
        {
            return y->parent_->nodeLeft_;
        }
        else
        {
            FatalErrorInFunction
                << "wrong addressing of the initial node"
                << exit(FatalError);
            SharedPointer<Node> ret;
            ret.setNULL();
            return ret;
        }
    }
    SharedPointer<Node> ret;
    ret.setNULL();
    return ret;
}

template <typename T, typename T2>
Foam::SharedPointer<typename Foam::parallelISAT<T, T2>::Node> Foam::parallelISAT<T, T2>::nodeSibling_leaf(SharedPointer<Leaf> x)
{
    if (size_leaf_ > 1)
    {
        if (x == x->node_2->leafLeft_)
        {
            // x is on the left, return right side
            return x->node_2->nodeRight_;
        }
        else if (x == x->node_2->leafRight_)
        {
            // x is on the right, return left side
            return x->node_2->nodeLeft_;
        }
        else
        {
            FatalErrorInFunction
                << "wrong addressing of the initial leaf"
                << exit(FatalError);

            SharedPointer<Node> ret;
            ret.setNULL();
            return ret;
        }
    }
    SharedPointer<Node> ret;
    ret.setNULL();
    return ret;
}

template <typename T, typename T2>
void Foam::parallelISAT<T, T2>::transplant(SharedPointer<Node> u, SharedPointer<Node> v)
{
    if (v.notNULL())
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

template <typename T, typename T2>
void Foam::parallelISAT<T, T2>::deleteLeaf(SharedPointer<Leaf> index_phi0)
{
    if (size_leaf_ == 1) // only one point is stored
    {
        leaf_manager.Delete(index_phi0);
        node_manager.Delete(root_);
        root_.setNULL();
        size_node--;
        //deleteDemandDrivenData(phi0);
        //deleteDemandDrivenData(root_);
    }
    else if (size_leaf_ > 1)
    {
        SharedPointer<Node> z = index_phi0->node_2;
        //label x;
        SharedPointer<Leaf> siblingPhi0 = leafPSibling(index_phi0);

        if (siblingPhi0.notNULL()) // the sibling of phi0 is a chemPoint
        {
            // z was root (only two chemPoints in the tree)
            if (z->parent_.isNULL())
            {
                root_ = node_manager.New(siblingPhi0->phi().size_);
                (*root_).mutex.unlock();

                if (root_.isNULL())
                {
                    FatalErrorInFunction
                        << "run out of memory"
                        << exit(FatalError);
                }
                size_node++;
                root_->nodeRight_.setNULL();
                root_->leafRight_.setNULL();
                root_->nodeLeft_.setNULL();
                root_->leafLeft_ = siblingPhi0;
                root_->parent_.setNULL();
                siblingPhi0->node_2 = root_;
            }
            else if (z == z->parent_->nodeLeft_)
            {
                z->parent_->leafLeft_ = siblingPhi0;
                z->parent_->nodeLeft_.setNULL();
                siblingPhi0->node_2 = z->parent_;
            }
            else if (z == z->parent_->nodeRight_)
            {
                z->parent_->leafRight_ = siblingPhi0;
                z->parent_->nodeRight_.setNULL();
                siblingPhi0->node_2 = z->parent_;
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
    if (subTreeRoot.notNULL())
    {
        deleteAllNode(subTreeRoot->nodeLeft_);
        deleteAllNode(subTreeRoot->nodeRight_);
        node_manager.Delete(subTreeRoot);
        size_node--;
        //deleteDemandDrivenData(subTreeRoot);
    }
}

template <typename T, typename T2>
void Foam::parallelISAT<T, T2>::deleteAll(SharedPointer<Node> subTreeRoot)
{
    if (subTreeRoot.notNULL())
    {
        deleteAll(subTreeRoot->nodeLeft_);
        deleteAll(subTreeRoot->nodeRight_);
        if (subTreeRoot->leafLeft_.notNULL())
        {
            leaf_manager.Delete(subTreeRoot->leafLeft_);
            size_leaf_--;
        }

        if (subTreeRoot->leafRight_.notNULL())
        {
            leaf_manager.Delete(subTreeRoot->leafRight_);
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
                       depth(subTreeRoot->nodeLeft_),
                       depth(subTreeRoot->nodeRight_));
    }
}
template <typename T, typename T2>
bool Foam::parallelISAT<T, T2>::treeCheck(SharedPointer<Node> index)
{
    bool leftnode = (index->nodeLeft_.notNULL());
    bool leftleaf = (index->leafLeft_.notNULL());
    if ((leftnode && leftleaf) || ((!leftnode) && (!leftleaf)))
    {
        FatalErrorInFunction
            << "p4 ";
        //<< abort(FatalError);
        return false;
    }
    bool rightnode = (index->nodeRight_.notNULL());
    bool rightleaf = (index->leafRight_.notNULL());
    if ((rightnode && rightleaf) || ((!rightnode) && (!rightleaf)))
    {
        FatalErrorInFunction
            << "p5 ";
        //<< abort(FatalError);
        return false;
    }
    if (leftnode)
    {
        if (index->nodeLeft_->parent_ != index)
        {
            FatalErrorInFunction
                << "p6 ";
            //<< abort(FatalError);

            return false;
        }
    }
    else
    {
        if (index->leafLeft_->node_2 != index)

        {
            FatalErrorInFunction
                << "p7 ";
            //<< abort(FatalError);
            return false;
        }
    }
    if (rightnode)
    {
        if (index->nodeRight_->parent_ != index)
        {
            FatalErrorInFunction
                << "p8 ";
            //<< abort(FatalError);
            return false;
        }
    }
    else
    {
        if (index->leafRight_->node_2 != index)
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
    if (phi0 == phi0->node_2->leafRight_) // phi0 is on the right
    {
        phi0->node_2->leafRight_.setNULL();
        phi0->node_2->nodeRight_ = newNode;

        return;
    }
    else if (phi0 == phi0->node_2->leafLeft_) // phi0 is on the left
    {
        phi0->node_2->leafLeft_.setNULL();
        phi0->node_2->nodeLeft_ = newNode;

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

    while (x.notNULL())
    {
        if (x->phi().size_ < 1)
        {
            FatalErrorInFunction
                << ",size=" << x->phi().size_
                << abort(FatalError);
        }
        chemPoints[chPi++] = x;
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
