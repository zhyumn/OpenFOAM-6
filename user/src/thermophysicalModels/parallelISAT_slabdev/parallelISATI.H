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
      pv_(manager_in), r_(manager_in, 12), emptylist_leaf(manager_in, 100000), emptylist_node(manager_in, 100000),
      size_leaf_(pv_().size_leaf), head_leaf(pv_().head_leaf),
      tail_leaf(pv_().tail_leaf), temp_tail_leaf(pv_().temp_tail_leaf),
      size_node(pv_().size_node), head_node(pv_().heaf_node),
      tail_node(pv_().tail_node), temp_tail_node(pv_().temp_tail_node),
      root_2(pv_().root_2), totalGrowth_(pv_().totalGrowth), //zzzz_(pv_().zzzz), zzzz2_(pv_().zzzz2),
      write_lock(manager_in, sync_in),
      leaf_manager(emptylist_leaf, tail_leaf, pv_().m_size_leaf, pv_().tmp_size_leaf),
      node_manager(emptylist_node, tail_node, pv_().m_size_node, pv_().tmp_size_node),
      nRetrieved_(0),
      nGrowth_(0),
      nAdd_(0)
{
    /*
    for (label i = 0;i < nmem_leaf_;i++)
    {
        emptylist_leaf[i] = i;
    }
    for (label i = 0;i < nmem_node_;i++)
    {
        emptylist_node[i] = i;
    }*/
    size_leaf_ = 0;
    tail_leaf = 0;
    tail_node = 0;
    //head_leaf = 0;
    //tail_leaf = nmem_leaf_ - 1;
    //temp_tail_leaf = nmem_leaf_ - 1;
    pv_().tmp_size_leaf = 0;
    pv_().tmp_size_node = 0;
    pv_().m_size_leaf = 0;
    pv_().m_size_node = 0;
    size_node = 0;
    //heaf_node = 0;
    //tail_node = nmem_node_ - 1;
    //temp_tail_node = nmem_node_ - 1;

    //root_ = -1;
    root_2.offset = -1;
    totalGrowth_ = 0;
    if (manager_.rank == 0)
        for (int i = 0; i < 12; i++)
            r_[i] = 0;

    //tleaf.node_ = &tnode;
    //tleaf.v = 1;
    //tnode.leaf_ = &tleaf;
    //tnode.v = 2;
    //zzzz_ = 0;
    //zzzz2_ = 0;
}
/*
Foam::ISAT::ISAT(const dataType& data)
:
    baseClassName(),
    data_(data)
{}


Foam::ISAT::ISAT(const ISAT&)
:
    baseClassName(),
    data_()
{}
*/

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //
/*
Foam::autoPtr<Foam::ISAT>
Foam::ISAT::New()
{
    return autoPtr<ISAT>(new ISAT);
}
*/

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
template <typename T, typename T2>
Foam::parallelISAT<T, T2>::~parallelISAT()
{
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
template <typename T, typename T2>
Foam::SharedPointer<typename Foam::parallelISAT<T, T2>::Leaf> Foam::parallelISAT<T, T2>::binaryTreeSearch(const scalarField &value)
{
    //std::cout << "Here:126 " << root_2.offset << std::endl;
    auto ret = binaryTreeSearch(value, root_2);
    //std::cout << "Here:128" << std::endl;
    return ret;
}
template <typename T, typename T2>
Foam::SharedPointer<typename Foam::parallelISAT<T, T2>::Leaf> Foam::parallelISAT<T, T2>::binaryTreeSearch(const inputType_old &value)
{
    //std::cout << "Here:126 " << root_2.offset << std::endl;
    auto ret = binaryTreeSearch(value, root_2);
    //std::cout << "Here:128" << std::endl;
    return ret;
}
template <typename T, typename T2>
Foam::SharedPointer<typename Foam::parallelISAT<T, T2>::Leaf> Foam::parallelISAT<T, T2>::binaryTreeSearch(
    const inputType_old &value,
    SharedPointer<Node> node_ptr)
{
    //std::cout << "Here:134 " << node_ptr.offset << std::endl;
    if (size_leaf_ > 1)
    {
        //scalar vPhi=0.0;
        //const scalarField& v = node->v();
        //const scalar& a = node->a();
        // compute v*phi
        //for (label i=0; i<phiq.size(); i++) vPhi += phiq[i]*v[i];
        //std::cout << "Here:145 " << std::endl;
        if (node_ptr->goLeft(value))
        {
            //std::cout << "Here:148 " << std::endl;
            if (!node_ptr->nodeLeft_2.isNULL())
            {
                //std::cout << "Here:151 " << std::endl;
                return binaryTreeSearch(value, node_ptr->nodeLeft_2);
            }
            else
            {
                //std::cout << "Here:156 " << std::endl;
                return node_ptr->leafLeft_2;
            }
        }
        else
        {
            //std::cout << "Here:162 " << std::endl;
            if (!node_ptr->nodeRight_2.isNULL())
            {
                //std::cout << "Here:165 " << node_ptr.offset << " " << node_ptr->nodeRight_2.offset << " " << node_ptr->nodeLeft_2.offset << " " << node_ptr->leafRight_2.offset << " " << node_ptr->leafLeft_2.offset << std::endl;
                return binaryTreeSearch(value, node_ptr->nodeRight_2);
            }
            else
            {
                //std::cout << "Here:170 " << std::endl;
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
    const scalarField &value,
    SharedPointer<Node> node_ptr)
{
    //std::cout << "Here:134 " << node_ptr.offset << std::endl;
    if (size_leaf_ > 1)
    {
        //scalar vPhi=0.0;
        //const scalarField& v = node->v();
        //const scalar& a = node->a();
        // compute v*phi
        //for (label i=0; i<phiq.size(); i++) vPhi += phiq[i]*v[i];
        //std::cout << "Here:145 " << std::endl;
        if (node_ptr->goLeft(value))
        {
            //std::cout << "Here:148 " << std::endl;
            if (!node_ptr->nodeLeft_2.isNULL())
            {
                //std::cout << "Here:151 " << std::endl;
                return binaryTreeSearch(value, node_ptr->nodeLeft_2);
            }
            else
            {
                //std::cout << "Here:156 " << std::endl;
                return node_ptr->leafLeft_2;
            }
        }
        else
        {
            //std::cout << "Here:162 " << std::endl;
            if (!node_ptr->nodeRight_2.isNULL())
            {
                //std::cout << "Here:165 " << node_ptr.offset << " " << node_ptr->nodeRight_2.offset << " " << node_ptr->nodeLeft_2.offset << " " << node_ptr->leafRight_2.offset << " " << node_ptr->leafLeft_2.offset << std::endl;
                return binaryTreeSearch(value, node_ptr->nodeRight_2);
            }
            else
            {
                //std::cout << "Here:170 " << std::endl;
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
    //std::cout << "Here:134 " << node_ptr.offset << std::endl;
    if (size_leaf_ > 1)
    {
        //scalar vPhi=0.0;
        //const scalarField& v = node->v();
        //const scalar& a = node->a();
        // compute v*phi
        //for (label i=0; i<phiq.size(); i++) vPhi += phiq[i]*v[i];
        //std::cout << "Here:145 " << std::endl;
        if (node_ptr->goLeft(value))
        {
            //std::cout << "Here:148 " << std::endl;
            if (!node_ptr->nodeLeft_2.isNULL())
            {
                //std::cout << "Here:151 " << std::endl;
                return binaryTreeSearch(value, node_ptr->nodeLeft_2);
            }
            else
            {
                //std::cout << "Here:156 " << std::endl;
                return node_ptr->leafLeft_2;
            }
        }
        else
        {
            //std::cout << "Here:162 " << std::endl;
            if (!node_ptr->nodeRight_2.isNULL())
            {
                //std::cout << "Here:165 " << node_ptr.offset << " " << node_ptr->nodeRight_2.offset << " " << node_ptr->nodeLeft_2.offset << " " << node_ptr->leafRight_2.offset << " " << node_ptr->leafLeft_2.offset << std::endl;
                return binaryTreeSearch(value, node_ptr->nodeRight_2);
            }
            else
            {
                //std::cout << "Here:170 " << std::endl;
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
    //write_lock.lock();
    bool success = false;
    SharedPointer<Leaf> &pleaf = binaryTreeSearch(x);
    outputType y_approx;
    pleaf->retrieve(x, y_approx);
    if (!pleaf->inEOA(x))
    {
        //outputType y_exact = Leaf::func(x);
        if (pleaf->checkSolution(x, y_approx))
        {
            SharedPointer<Leaf> &pleaf_new = leaf_manager.New();

            if (!pleaf_new.isNull())
            {
                *pleaf_new = *pleaf;
                pleaf_new->grow(x); //.EOA = fabs(x - leafmem[leaf_index].v);
                //label node_index = leafmem[index].node_;

                if (pleaf->node_2->leafLeft_2.offset == pleaf.offset)
                {
                    pleaf->node_2->leafLeft_2.offset = pleaf_new.offset;
                }
                else
                {
                    pleaf->node_2->leafRight_2.offset = pleaf_new.offset;
                }
                //std::cout << "215 CCCCCCCCCCCCCCCCCC " << pleaf.offset << std::endl;
                leaf_manager.Delete(pleaf);
                success = true;
            }
        }
    }
    //write_lock.unlock();
    return success;
}
template <typename T, typename T2>
bool Foam::parallelISAT<T, T2>::grow(const scalarField &x, SharedPointer<Leaf> &pleaf)
{
    //write_lock.lock();
    bool success = false;
    //SharedPointer<Leaf> &pleaf = binaryTreeSearch(x);
    //outputType y_approx;
    //pleaf->retrieve(x, y_approx);

    SharedPointer<Leaf> pleaf_new = leaf_manager.New(*pleaf);
    //SharedPointer<Leaf> pleaf_new = leaf_manager.New(x.size());

    if (!pleaf_new.isNULL())
    {
        //pleaf_new->init(*pleaf);
        pleaf_new->node_2.offset = pleaf->node_2.offset;
        success = pleaf_new->grow(x); //.EOA = fabs(x - leafmem[leaf_index].v);
                                      //label node_index = leafmem[index].node_;
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
            //std::cout << "215 CCCCCCCCCCCCCCCCCC " << pleaf.offset << std::endl;
            leaf_manager.Delete(pleaf);
        }
        else
        {
            leaf_manager.Delete(pleaf_new);
        }
        //success = true;
    }

    //write_lock.unlock();
    return success;
}

template <typename T, typename T2>
bool Foam::parallelISAT<T, T2>::retrieve(const inputType_old &x, outputType_old &y)
{
    //std::cout << "Here:215" << std::endl;
    //static label zzz = 0;
    /*if (!this->treeCheck())
    {
        FatalErrorInFunction
            << "broken tree retr in func 0\n"
            << abort(FatalError);
    }*/
    SharedPointer<Leaf> pleaf = binaryTreeSearch(x);
    //std::cout << "Here:224" << std::endl;
    if (pleaf.isNULL())
    {
        //y = Leaf::func(x);
        //Pout << "add" << endl;
        //add(x);
        /*if (!this->treeCheck())
        {
            FatalErrorInFunction
                << "broken tree retr in func 1\n"
                << abort(FatalError);
        }*/
        return false;
    }
    //outputType y_approx;
    //leafmem[index].retrieve(x, y_approx);
    /*if (T::out == 1)
    {
        zzz++;
        if (zzz == 3)
        {
            for (int i = 0;i < 13;i++)
            {
                for (int j = 0;j < 13;j++)
                {
                    Info << leafmem[index].LT_(i, j) << " ";
                }
                Info << endl;
            }

            ErrorInFunction
                << "out=" << T::out
                << "\nnActiveSpecies_=" << leafmem[index].nActiveSpecies_
                << abort(FatalError);
        }
    }*/
    //std::cout << "Here:260" << std::endl;
    if (pleaf->inEOA(x))
    {
        //Pout << "inEOA" << endl;
        /*if (!this->treeCheck())
        {
            FatalErrorInFunction
                << "broken tree retr in func 1.5\n"
                << abort(FatalError);
        }*/
        pleaf->retrieve(x, y);

        nRetrieved_++;
        /*if (!this->treeCheck())
        {
            FatalErrorInFunction
                << "broken tree retr in func 2\n"
                << abort(FatalError);
        }*/
        return true;
    }
    else
    {
        //zzzz2_++;
        /*if (zzzz2_ == 4)
        {
            Pout << "x=" << x << endl;
            //Pout << "y=" << y << endl;
            //FatalErrorInFunction
            //    << abort(FatalError);
        }*/

        //Pout << "not inEOA" << endl;
        //y = Leaf::func(x);
        //if (fabs(y_approx - y) > tolerance_ || !grow(x)) //ToDo
        //{
        //Pout << "add" << endl;
        //    add(x);
        //}
        /*if (!this->treeCheck())
        {
            FatalErrorInFunction
                << "broken tree retr in func 3\n"
                << abort(FatalError);
        }*/
        return false;
    }
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
        //std::cout << this->manager_.rank << "here:502" << std::endl;
        pleaf->retrieve(x, y);
        //std::cout << this->manager_.rank << "here:504" << std::endl;
        nRetrieved_++;
        return true;
    }
    else
    {
        return false;
    }
}
/*
template<typename T, typename T2>
void Foam::parallelISAT<T,T2>::insert(const inputType& value)
{
    write_lock.lock();
    //Pout << "haha" << endl;
    //Pout << "root_=" << root_ << endl;
    if (size_leaf == 0)
    {
        //Pout << "haha1" << endl;
        label node_index = node_manager.New();
        if (node_index == -1)
        {
            write_lock.unlock();
            return;
        }
        label leaf_index = leaf_manager.New();
        if (leaf_index == -1)
        {
            node_manager.Delete(node_index);
            write_lock.unlock();
            return;
        }
        //Pout << "haha2,leaf_index=" << leaf_index << endl;
        nodemem[node_index].parent_ = -1;
        nodemem[node_index].nodeLeft_ = -1;
        nodemem[node_index].nodeRight_ = -1;
        nodemem[node_index].leafLeft_ = leaf_index;
        nodemem[node_index].leafRight_ = -1;
        nodemem[node_index].set();
        //Pout << "haha3,node_index=" << node_index << endl;
        leafmem[leaf_index].node_ = node_index;
        leafmem[leaf_index].set(value);
        root_ = node_index;
        size_leaf++;
        size_node++;
        //Pout << value << "inserted!" << endl;
        write_lock.unlock();
        return;
    }
    else
    {
        //Pout << "heihei2,value=" << value << endl;
        //Todo  leaf_index can be obtain by arg
        label leaf_index = binaryTreeSearch(value, root_);

        if (leafmem[leaf_index].collide(value))
        {
            write_lock.unlock();
            return;
        }
        label parentNode = leafmem[leaf_index].node_;

        label new_node_index = node_manager.New();
        if (new_node_index == -1)
        {
            write_lock.unlock();
            return;
        }
        label new_leaf_index = leaf_manager.New();
        if (new_leaf_index == -1)
        {
            node_manager.Delete(new_node_index);
            write_lock.unlock();
            return;
        }
        //Pout << "heihei3,new_leaf_index=" << new_leaf_index << endl;
        leafmem[new_leaf_index].set(value);
        leafmem[new_leaf_index].node_ = new_node_index;
        leafmem[leaf_index].node_ = new_node_index;

        nodemem[new_node_index].parent_ = -1;
        nodemem[new_node_index].nodeLeft_ = -1;
        nodemem[new_node_index].nodeRight_ = -1;
        nodemem[new_node_index].set(leafmem[new_leaf_index], leafmem[leaf_index]);//v = (leafmem[new_leaf_index].v + leafmem[leaf_index].v) / 2;
        if (nodemem[new_node_index].goLeft(leafmem[new_leaf_index]))
        {
            nodemem[new_node_index].leafLeft_ = new_leaf_index;
            nodemem[new_node_index].leafRight_ = leaf_index;
        }
        else
        {
            nodemem[new_node_index].leafLeft_ = leaf_index;
            nodemem[new_node_index].leafRight_ = new_leaf_index;
        }
        if (size_leaf > 1)
        {
            nodemem[new_node_index].parent_ = parentNode;
            if (nodemem[parentNode].leafLeft_ == leaf_index)
            {
                nodemem[parentNode].nodeLeft_ = new_node_index;
                nodemem[parentNode].leafLeft_ = -1;
            }
            else
            {
                nodemem[parentNode].nodeRight_ = new_node_index;
                nodemem[parentNode].leafRight_ = -1;
            }
            size_leaf++;
            size_node++;
            //Pout << value << "inserted!" << endl;
            write_lock.unlock();
            return;
        }
        else
        {
            root_ = new_node_index;
            node_manager.Delete(parentNode);
            size_leaf++;
            //Pout << value << "inserted!" << endl;
            write_lock.unlock();
            return;
        }
    }
}
*/
template <typename T, typename T2>
Foam::label Foam::parallelISAT<T, T2>::add(
    const inputType &x)
{
    return add(x, T::leafData::func(x));
}

template <typename T, typename T2>
bool Foam::parallelISAT<T, T2>::isFull()
{
    return size_leaf_ >= maxNLeafs_;
}
template <typename T, typename T2>
template <typename... Args>
Foam::label Foam::parallelISAT<T, T2>::add(
    const inputType_old &x, const outputType_old &y, Args... args)
{
    write_lock.lock();
    //if (!treeCheck())
    //{
    //    FatalErrorInFunction
    //        << "pos0"
    //        << abort(FatalError);
    //}
    //static label out0 = 0, out1 = 0, out2 = 0, out3 = 0, out4 = 0, out5 = 0, zz0 = 0, zz1 = 0;
    /*if (out0 == 100)
    {
        FatalErrorInFunction
            << "out0=" << out0
            << "\nout1=" << out1
            << "\nout2=" << out2
            << "\nout3=" << out3
            << "\nout4=" << out4
            << "\nout5=" << out5
            << abort(FatalError);
    }*/

    //out0++;
    label growthOrAddFlag = 1;

    if (isFull())
    {
        write_lock.unlock();
        return growthOrAddFlag;
    }
    //Pout << "haha" << endl;
    //Pout << "root_=" << root_ << endl;
    if (size_leaf_ == 0)
    {
        //if (T::out == 1)
        //    Pout << "haha1" << endl;
        SharedPointer<Node> pnode = node_manager.New(x.size());
        //std::cout << "Here:472" << std::endl;
        //label node_index = node_manager.New();
        if (pnode.isNULL())
        {
            growthOrAddFlag = 0;
            //if (!treeCheck())
            //{
            //    FatalErrorInFunction
            //        << "pos01"
            //         << abort(FatalError);
            //}
            write_lock.unlock();
            return growthOrAddFlag;
        }
        SharedPointer<Leaf> pleaf = leaf_manager.New(x.size());
        //std::cout << "Here:488" << std::endl;
        //label leaf_index = leaf_manager.New();
        if (pleaf.isNULL())
        {
            growthOrAddFlag = 0;
            node_manager.Delete(pnode);
            //std::cout << "Here:492" << std::endl;
            //if (!treeCheck())
            //{
            //    FatalErrorInFunction
            //        << "pos02"
            //        << abort(FatalError);
            //}
            write_lock.unlock();
            return growthOrAddFlag;
        }
        //Pout << "haha2,leaf_index=" << leaf_index << endl;
        pnode->parent_2.offset = -1;
        pnode->nodeLeft_2.offset = -1;
        pnode->nodeRight_2.offset = -1;
        pnode->leafLeft_2.offset = pleaf.offset;
        pnode->leafRight_2.offset = -1;
        pnode->set();
        //Pout << "haha3,node_index=" << node_index << endl;
        pleaf->node_2.offset = pnode.offset;
        pleaf->set(x, y, args...);
        root_2.offset = pnode.offset;
        size_leaf_++;
        size_node++;
        nAdd_++;

        //out1++;
        //if (!treeCheck())
        //{
        //    FatalErrorInFunction
        //        << "pos1"
        //        << abort(FatalError);
        //}
        //Pout << value << "inserted!" << endl;
        write_lock.unlock();
        return growthOrAddFlag;
    }
    SharedPointer<Leaf> pleaf = binaryTreeSearch(x);
    //outputType y_approx = leafmem[index].retrieve(value);
    //Pout << "aa!" << leafmem[index].inEOA(value) << endl;
    if (pleaf->inEOA(x))
    {
        //if (T::out == 1)
        //    Pout << "haha2" << endl;
        growthOrAddFlag = 0;
        //if (!treeCheck())
        //{
        //    FatalErrorInFunction
        //        << "pos1,1"
        //        << abort(FatalError);
        //}
        write_lock.unlock();
        //out2++;

        return growthOrAddFlag;
    }
    if (pleaf->checkSolution(x, y))
    {
        /*if (zzzz2_ == 4 && x[0] < 1e-10)
        {
            Pout << "try grow"
                << x;
            //<< abort(FatalError);
        }*/

        //zz0++;
        if (pleaf->grow(x))
        {
            /*if (zzzz2_ == 4 && x[0] < 1e-10)
            {
                FatalErrorInFunction
                    << "grow ss"
                    << abort(FatalError);
            }*/

            //zzzz_++;
            //zz1++;
            //if (T::out == 1)
            //    Pout << "haha3" << endl;
            nGrowth_++;
            totalGrowth_++;
            //out3++;
            growthOrAddFlag = 0;
            //if (!treeCheck())
            //{
            //    FatalErrorInFunction
            //        << "pos2"
            //        << abort(FatalError);
            //}
            write_lock.unlock();
            return growthOrAddFlag;
        }
        //if (zzzz2_ == 4 && x[0] < 1e-10)
        //{
        //    FatalErrorInFunction
        //        << "grow f"
        //        << abort(FatalError);
        //}
    }

    //Pout << "!!!!!!!!!!zz0=" << zz0 << ",zz1=" << zz1 << endl;
    //Pout << "heihei2,value=" << value << endl;
    //Todo  leaf_index can be obtain by arg
    SharedPointer<Leaf> pleaf2 = binaryTreeSearch(x, root_2);

    if (pleaf2->inEOA(x))
    {
        //if (T::out == 1)
        //    Pout << "haha2" << endl;
        growthOrAddFlag = 0;
        //if (!treeCheck())
        //{
        //    FatalErrorInFunction
        //        << "pos2,1"
        //        << abort(FatalError);
        //}
        write_lock.unlock();
        return growthOrAddFlag;
    }
    SharedPointer<Node> parentNode = pleaf2->node_2;

    SharedPointer<Node> pnode_new = node_manager.New(x.size());

    //std::cout << "BBBBBBBBBBBBBBBB " << parentNode.offset << std::endl;
    //std::cout << "BBBBBBBBBBBBBBBB " << pnode_new.offset << std::endl;
    //std::cout << "Here:612" << std::endl;
    if (pnode_new.isNULL())
    {
        growthOrAddFlag = 0;
        //if (!treeCheck())
        //{
        //   FatalErrorInFunction
        //        << "pos2,2"
        //        << abort(FatalError);
        //}
        write_lock.unlock();
        return growthOrAddFlag;
    }
    SharedPointer<Leaf> pleaf_new = leaf_manager.New(x.size());
    //std::cout << "Here:629" << std::endl;
    if (pleaf_new.isNULL())
    {
        growthOrAddFlag = 0;
        //std::cout << "651 CCCCCCCCCCCCCCCCCC " << pnode_new.offset << std::endl;
        node_manager.Delete(pnode_new);
        //std::cout << "Here:631" << std::endl;
        //if (!treeCheck())
        //{
        //    FatalErrorInFunction
        //        << "pos2,3"
        //        << abort(FatalError);
        //}
        write_lock.unlock();
        return growthOrAddFlag;
    }
    //std::cout << "BBBBBBBBBBBBBBBBC " << parentNode.offset << std::endl;
    //std::cout << "Here:644" << std::endl;
    //Pout << "heihei3,new_leaf_index=" << new_leaf_index << endl;
    pleaf_new->set(x, y, args...);
    pleaf_new->node_2.offset = pnode_new.offset;
    //std::cout << "BBBBBBBBBBBBBBBBF " << parentNode.offset << " " << pleaf2->node_2.offset << " " << pnode_new.offset << std::endl;
    pleaf2->node_2.offset = pnode_new.offset;
    //std::cout << "Here:649" << std::endl;
    //std::cout << "BBBBBBBBBBBBBBBBE " << parentNode.offset << " " << pleaf2->node_2.offset << " " << pnode_new.offset << std::endl;
    pnode_new->parent_2.offset = -1;
    pnode_new->nodeLeft_2.offset = -1;
    pnode_new->nodeRight_2.offset = -1;
    pnode_new->set(*pleaf_new, *pleaf2); //v = (leafmem[new_leaf_index].v + leafmem[leaf_index].v) / 2;
    //std::cout << "Here:654" << std::endl;
    //std::cout << "BBBBBBBBBBBBBBBBD " << parentNode.offset << std::endl;
    if (pnode_new->goLeft(*pleaf_new))
    {
        //std::cout << "Here:657" << std::endl;
        pnode_new->leafLeft_2.offset = pleaf_new.offset;
        pnode_new->leafRight_2.offset = pleaf2.offset;
    }
    else
    {
        //std::cout << "Here:662" << std::endl;
        pnode_new->leafLeft_2.offset = pleaf2.offset;
        pnode_new->leafRight_2.offset = pleaf_new.offset;
    }
    //std::cout << "BBBBBBBBBBBBBBBBA " << parentNode.offset << std::endl;
    if (size_leaf_ > 1)
    {
        //std::cout << "Here:668" << std::endl;
        pnode_new->parent_2.offset = parentNode.offset;
        if (parentNode->leafLeft_2.offset == pleaf2.offset)
        {
            parentNode->nodeLeft_2.offset = pnode_new.offset;
            parentNode->leafLeft_2.offset = -1;
        }
        else
        {
            //std::cout << "BBBBBBBBBBBBBBBBZ " << parentNode.offset << std::endl;
            //std::cout << "BBBBBBBBBBBBBBBBZ " << pnode_new.offset << std::endl;
            parentNode->nodeRight_2.offset = pnode_new.offset;
            /*if (parentNode->nodeRight_2.offset == parentNode.offset)
            {
                FatalErrorInFunction
                    << "696 "
                    << parentNode.offset
                    << " "
                    << parentNode->nodeRight_2.offset
                    << exit(FatalError);
            }*/
            parentNode->leafRight_2.offset = -1;
        }
        size_leaf_++;
        size_node++;
        nAdd_++;
        //out4++;
        //if (T::out == 1)
        //    Pout << "haha3" << endl;
        //if (!treeCheck())
        //{
        //    FatalErrorInFunction
        //        << "pos3"
        //        << abort(FatalError);
        //}
        //Pout << value << "inserted!" << endl;
        write_lock.unlock();
        //std::cout << "Here:694" << std::endl;
        return growthOrAddFlag;
    }
    else
    {
        root_2.offset = pnode_new.offset;
        //std::cout << "728 CCCCCCCCCCCCCCCCCC " << parentNode.offset << std::endl;
        node_manager.Delete(parentNode);
        //std::cout << "Here:701" << std::endl;
        size_leaf_++;
        nAdd_++;
        //out5++;
        //if (T::out == 1)
        //    Pout << "haha4" << endl;
        //Pout << value << "inserted!" << endl;
        //Pout << "pos4" << endl;
        //if (!treeCheck())
        //{
        //    FatalErrorInFunction
        //        << "pos4"
        //        << abort(FatalError);
        //}
        write_lock.unlock();

        return growthOrAddFlag;
    }
}

template <typename T, typename T2>
template <typename... Args>
Foam::label Foam::parallelISAT<T, T2>::add(
    const scalarField &x, const scalarField &y, Args... args)
{
    write_lock.lock();
    //std::cout << "add=========";
    //for (int i = 0; i < 12; i++)
    //{
    //    std::cout << "  " << r_[i];
    //}
    //std::cout << std::endl;

    r_[0]++;
    //if (!treeCheck())
    //{
    //    FatalErrorInFunction
    //        << "pos0"
    //        << abort(FatalError);
    //}
    //static label out0 = 0, out1 = 0, out2 = 0, out3 = 0, out4 = 0, out5 = 0, zz0 = 0, zz1 = 0;
    /*if (out0 == 100)
    {
        FatalErrorInFunction
            << "out0=" << out0
            << "\nout1=" << out1
            << "\nout2=" << out2
            << "\nout3=" << out3
            << "\nout4=" << out4
            << "\nout5=" << out5
            << abort(FatalError);
    }*/

    //out0++;
    label growthOrAddFlag = 1;

    if (isFull())
    {
        r_[1]++;
        write_lock.unlock();
        return growthOrAddFlag;
    }
    //Pout << "haha" << endl;
    //Pout << "root_=" << root_ << endl;
    if (size_leaf_ == 0)
    {
        //if (T::out == 1)
        //    Pout << "haha1" << endl;
        SharedPointer<Node> pnode = node_manager.New(x.size());
        //std::cout << "Here:472" << std::endl;
        //label node_index = node_manager.New();
        if (pnode.isNULL())
        {
            growthOrAddFlag = 0;
            //if (!treeCheck())
            //{
            //    FatalErrorInFunction
            //        << "pos01"
            //         << abort(FatalError);
            //}
            r_[2]++;
            write_lock.unlock();
            return growthOrAddFlag;
        }
        SharedPointer<Leaf> pleaf = leaf_manager.New(x.size());
        //std::cout << "Here:488" << std::endl;
        //label leaf_index = leaf_manager.New();
        if (pleaf.isNULL())
        {
            growthOrAddFlag = 0;
            node_manager.Delete(pnode);
            //std::cout << "Here:492" << std::endl;
            //if (!treeCheck())
            //{
            //    FatalErrorInFunction
            //        << "pos02"
            //        << abort(FatalError);
            //}
            r_[3]++;
            write_lock.unlock();
            return growthOrAddFlag;
        }
        //Pout << "haha2,leaf_index=" << leaf_index << endl;
        pnode->parent_2.offset = -1;
        pnode->nodeLeft_2.offset = -1;
        pnode->nodeRight_2.offset = -1;
        pnode->leafLeft_2.offset = pleaf.offset;
        pnode->leafRight_2.offset = -1;
        pnode->set();
        //Pout << "haha3,node_index=" << node_index << endl;
        pleaf->node_2.offset = pnode.offset;
        pleaf->set(x, y, args...);
        root_2.offset = pnode.offset;
        size_leaf_++;
        size_node++;
        nAdd_++;

        //out1++;
        //if (!treeCheck())
        //{
        //    FatalErrorInFunction
        //        << "pos1"
        //        << abort(FatalError);
        //}
        //Pout << value << "inserted!" << endl;
        r_[4]++;
        write_lock.unlock();
        return growthOrAddFlag;
    }
    SharedPointer<Leaf> pleaf = binaryTreeSearch(x);
    //outputType y_approx = leafmem[index].retrieve(value);
    //Pout << "aa!" << leafmem[index].inEOA(value) << endl;
    if (pleaf->inEOA(x))
    {
        //if (T::out == 1)
        //    Pout << "haha2" << endl;
        growthOrAddFlag = 0;
        //if (!treeCheck())
        //{
        //    FatalErrorInFunction
        //        << "pos1,1"
        //        << abort(FatalError);
        //}
        r_[5]++;
        write_lock.unlock();
        //out2++;

        return growthOrAddFlag;
    }
    if (pleaf->checkSolution(x, y))
    {
        /*if (zzzz2_ == 4 && x[0] < 1e-10)
        {
            Pout << "try grow"
                << x;
            //<< abort(FatalError);
        }*/

        //zz0++;
        if (grow(x, pleaf))
        {
            /*if (zzzz2_ == 4 && x[0] < 1e-10)
            {
                FatalErrorInFunction
                    << "grow ss"
                    << abort(FatalError);
            }*/

            //zzzz_++;
            //zz1++;
            //if (T::out == 1)
            //    Pout << "haha3" << endl;
            nGrowth_++;
            totalGrowth_++;
            //out3++;
            growthOrAddFlag = 0;
            //if (!treeCheck())
            //{
            //    FatalErrorInFunction
            //        << "pos2"
            //        << abort(FatalError);
            //}
            r_[6]++;
            write_lock.unlock();
            return growthOrAddFlag;
        }
        //if (zzzz2_ == 4 && x[0] < 1e-10)
        //{
        //    FatalErrorInFunction
        //        << "grow f"
        //        << abort(FatalError);
        //}
    }

    //Pout << "!!!!!!!!!!zz0=" << zz0 << ",zz1=" << zz1 << endl;
    //Pout << "heihei2,value=" << value << endl;
    //Todo  leaf_index can be obtain by arg
    SharedPointer<Leaf> pleaf2 = binaryTreeSearch(x, root_2);

    if (pleaf2->inEOA(x))
    {
        //if (T::out == 1)
        //    Pout << "haha2" << endl;
        growthOrAddFlag = 0;
        //if (!treeCheck())
        //{
        //    FatalErrorInFunction
        //        << "pos2,1"
        //        << abort(FatalError);
        //}
        r_[7]++;
        write_lock.unlock();
        return growthOrAddFlag;
    }
    SharedPointer<Node> parentNode = pleaf2->node_2;

    SharedPointer<Node> pnode_new = node_manager.New(x.size());

    //std::cout << "BBBBBBBBBBBBBBBB " << parentNode.offset << std::endl;
    //std::cout << "BBBBBBBBBBBBBBBB " << pnode_new.offset << std::endl;
    //std::cout << "Here:612" << std::endl;
    if (pnode_new.isNULL())
    {
        growthOrAddFlag = 0;
        //if (!treeCheck())
        //{
        //   FatalErrorInFunction
        //        << "pos2,2"
        //        << abort(FatalError);
        //}
        r_[8]++;
        write_lock.unlock();
        return growthOrAddFlag;
    }
    SharedPointer<Leaf> pleaf_new = leaf_manager.New(x.size());
    //std::cout << "Here:629" << std::endl;
    if (pleaf_new.isNULL())
    {
        growthOrAddFlag = 0;
        //std::cout << "651 CCCCCCCCCCCCCCCCCC " << pnode_new.offset << std::endl;
        node_manager.Delete(pnode_new);
        //std::cout << "Here:631" << std::endl;
        //if (!treeCheck())
        //{
        //    FatalErrorInFunction
        //        << "pos2,3"
        //        << abort(FatalError);
        //}
        r_[9]++;
        write_lock.unlock();
        return growthOrAddFlag;
    }
    //std::cout << "BBBBBBBBBBBBBBBBC " << parentNode.offset << std::endl;
    //std::cout << "Here:644" << std::endl;
    //Pout << "heihei3,new_leaf_index=" << new_leaf_index << endl;
    pleaf_new->set(x, y, args...);
    pleaf_new->node_2.offset = pnode_new.offset;
    //std::cout << "BBBBBBBBBBBBBBBBF " << parentNode.offset << " " << pleaf2->node_2.offset << " " << pnode_new.offset << std::endl;
    pleaf2->node_2.offset = pnode_new.offset;
    //std::cout << "Here:649" << std::endl;
    //std::cout << "BBBBBBBBBBBBBBBBE " << parentNode.offset << " " << pleaf2->node_2.offset << " " << pnode_new.offset << std::endl;
    pnode_new->parent_2.offset = -1;
    pnode_new->nodeLeft_2.offset = -1;
    pnode_new->nodeRight_2.offset = -1;
    pnode_new->set(*pleaf_new, *pleaf2); //v = (leafmem[new_leaf_index].v + leafmem[leaf_index].v) / 2;
    //std::cout << "Here:654" << std::endl;
    //std::cout << "BBBBBBBBBBBBBBBBD " << parentNode.offset << std::endl;
    if (pnode_new->goLeft(*pleaf_new))
    {
        //std::cout << "Here:657" << std::endl;
        pnode_new->leafLeft_2.offset = pleaf_new.offset;
        pnode_new->leafRight_2.offset = pleaf2.offset;
    }
    else
    {
        //std::cout << "Here:662" << std::endl;
        pnode_new->leafLeft_2.offset = pleaf2.offset;
        pnode_new->leafRight_2.offset = pleaf_new.offset;
    }
    //std::cout << "BBBBBBBBBBBBBBBBA " << parentNode.offset << std::endl;
    if (size_leaf_ > 1)
    {
        //std::cout << "Here:668" << std::endl;
        pnode_new->parent_2.offset = parentNode.offset;
        if (parentNode->leafLeft_2.offset == pleaf2.offset)
        {
            parentNode->nodeLeft_2.offset = pnode_new.offset;
            parentNode->leafLeft_2.offset = -1;
        }
        else
        {
            //std::cout << "BBBBBBBBBBBBBBBBZ " << parentNode.offset << std::endl;
            //std::cout << "BBBBBBBBBBBBBBBBZ " << pnode_new.offset << std::endl;
            parentNode->nodeRight_2.offset = pnode_new.offset;
            /*if (parentNode->nodeRight_2.offset == parentNode.offset)
            {
                FatalErrorInFunction
                    << "696 "
                    << parentNode.offset
                    << " "
                    << parentNode->nodeRight_2.offset
                    << exit(FatalError);
            }*/
            parentNode->leafRight_2.offset = -1;
        }
        size_leaf_++;
        size_node++;
        nAdd_++;
        //out4++;
        //if (T::out == 1)
        //    Pout << "haha3" << endl;
        //if (!treeCheck())
        //{
        //    FatalErrorInFunction
        //        << "pos3"
        //        << abort(FatalError);
        //}
        //Pout << value << "inserted!" << endl;
        r_[10]++;
        write_lock.unlock();
        //std::cout << "Here:694" << std::endl;
        return growthOrAddFlag;
    }
    else
    {
        root_2.offset = pnode_new.offset;
        //std::cout << "728 CCCCCCCCCCCCCCCCCC " << parentNode.offset << std::endl;
        node_manager.Delete(parentNode);
        //std::cout << "Here:701" << std::endl;
        size_leaf_++;
        nAdd_++;
        //out5++;
        //if (T::out == 1)
        //    Pout << "haha4" << endl;
        //Pout << value << "inserted!" << endl;
        //Pout << "pos4" << endl;
        //if (!treeCheck())
        //{
        //    FatalErrorInFunction
        //        << "pos4"
        //        << abort(FatalError);
        //}
        r_[11]++;
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
        //std::cout << "956 CCCCCCCCCCCCCCCCCC " << index_phi0.offset << std::endl;
        leaf_manager.Delete(index_phi0);
        //std::cout << "Here:927" << std::endl;
        //std::cout << "962 CCCCCCCCCCCCCCCCCC " << root_2.offset << std::endl;
        node_manager.Delete(root_2);
        //std::cout << "Here:925" << std::endl;
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
                //std::cout << "Here:943" << std::endl;
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
                //root_->leafLeft() = siblingPhi0;
                siblingPhi0->node_2.offset = root_2.offset;
            }
            else if (z.offset == z->parent_2->nodeLeft_2.offset)
            {
                z->parent_2->leafLeft_2.offset = siblingPhi0.offset;
                //z->parent()->leafLeft() = siblingPhi0;
                z->parent_2->nodeLeft_2.offset = -1;
                //z->parent()->nodeLeft() = nullptr;
                siblingPhi0->node_2.offset = z->parent_2.offset;
                //siblingPhi0->node() = z->parent();
            }
            else if (z.offset == z->parent_2->nodeRight_2.offset)
            {
                z->parent_2->leafRight_2.offset = siblingPhi0.offset;
                //z->parent()->leafRight() = siblingPhi0;
                z->parent_2->nodeRight_2.offset = -1;
                //z->parent()->nodeRight() = nullptr;
                siblingPhi0->node_2.offset = z->parent_2.offset;
                //siblingPhi0->node() = z->parent();
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
        //std::cout << "1033 CCCCCCCCCCCCCCCCCC " << index_phi0.offset << std::endl;

        leaf_manager.Delete(index_phi0);
        //std::cout << "Here:1004" << std::endl;
        //std::cout << "1041 CCCCCCCCCCCCCCCCCC " << z.offset << std::endl;
        node_manager.Delete(z);
        //std::cout << "Here:1006" << std::endl;
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
        //std::cout << "1058 CCCCCCCCCCCCCCCCCC " << subTreeRoot.offset << std::endl;
        node_manager.Delete(subTreeRoot);
        //std::cout << "Here:1016" << std::endl;
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
        //phi0->node()->leafRight() = nullptr;
        //phi0->node()->nodeRight() = newNode;
        return;
    }
    else if (phi0.offset == phi0->node_2->leafLeft_2.offset) // phi0 is on the left
    {
        phi0->node_2->leafLeft_2.offset = -1;
        phi0->node_2->nodeLeft_2.offset = newNode.offset;
        //phi0->node()->leafLeft() = nullptr;
        //phi0->node()->nodeLeft() = newNode;
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
    /*     forAll(chemPoints, j)
    {
        const typename DataType::inputType &phij = chemPoints[j]->phi();
        std::cout << chPi << "ccx =" << j << " " << phij.size_ << std::endl;
    } */
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
