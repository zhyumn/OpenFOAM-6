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
template<typename T, typename T2>
template<typename ...Args>
Foam::parallelISAT<T, T2>::parallelISAT(SUPstream::mpi_manager& manager_in, label nmem, label maxNleafs, SUPstream::mpi_sync& sync_in, Args ...args)
    :T2(args...), manager_(manager_in),sync_(sync_in), maxNLeafs_(maxNleafs), nmem_leaf_(nmem), nmem_node_(nmem), leafmem(manager_in, nmem_leaf_),
    emptylist_leaf(manager_in, nmem_leaf_), nodemem(manager_in, nmem_node_),
    emptylist_node(manager_in, nmem_node_), pv_(manager_in),
    size_leaf_(pv_().size_leaf), head_leaf(pv_().head_leaf),
    tail_leaf(pv_().tail_leaf), temp_tail_leaf(pv_().temp_tail_leaf),
    size_node(pv_().size_node), head_node(pv_().heaf_node),
    tail_node(pv_().tail_node), temp_tail_node(pv_().temp_tail_node),
    root_(pv_().root),totalGrowth_(pv_().totalGrowth), //zzzz_(pv_().zzzz), zzzz2_(pv_().zzzz2),
    write_lock(manager_in, sync_in),
    leaf_manager(emptylist_leaf, head_leaf, tail_leaf, temp_tail_leaf, nmem_leaf_),
    node_manager(emptylist_node, head_node, tail_node, temp_tail_node, nmem_node_),
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
    //head_leaf = 0;
    //tail_leaf = nmem_leaf_ - 1;
    //temp_tail_leaf = nmem_leaf_ - 1;

    size_node = 0;
    //heaf_node = 0;
    //tail_node = nmem_node_ - 1;
    //temp_tail_node = nmem_node_ - 1;

    root_ = -1;
    totalGrowth_=0;

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
template<typename T, typename T2>
Foam::parallelISAT<T, T2>::~parallelISAT()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
template<typename T, typename T2>
Foam::label Foam::parallelISAT<T, T2>::binaryTreeSearch(inputType value)
{
    return binaryTreeSearch(value, root_);
}
template<typename T, typename T2>
Foam::label Foam::parallelISAT<T, T2>::binaryTreeSearch
(
    inputType value,
    label node_index
)
{
    if (size_leaf_ > 1)
    {
        //scalar vPhi=0.0;
        //const scalarField& v = node->v();
        //const scalar& a = node->a();
        // compute v*phi
        //for (label i=0; i<phiq.size(); i++) vPhi += phiq[i]*v[i];
        if (nodemem[node_index].goLeft(value))
        {
            if (nodemem[node_index].nodeLeft_ != -1)
            {
                return binaryTreeSearch(value, nodemem[node_index].nodeLeft_);
            }
            else
            {
                return nodemem[node_index].leafLeft_;
            }
        }
        else
        {
            if (nodemem[node_index].nodeRight_ != -1)
            {
                return binaryTreeSearch(value, nodemem[node_index].nodeRight_);
            }
            else
            {
                return nodemem[node_index].leafRight_;
            }

        }
    }
    // only one point stored (left element of the root)
    else if (size_leaf_ == 1)
    {
        return nodemem[root_].leafLeft_;
    }
    else // no point stored
    {
        return -1;
    }
}
template<typename T, typename T2>
bool Foam::parallelISAT<T, T2>::grow(const inputType& x)
{
    //write_lock.lock();
    bool success = false;
    label index = binaryTreeSearch(x);
    outputType y_approx;
    leafmem[index].retrieve(x, y_approx);
    if (!leafmem[index].inEOA(x))
    {
        //outputType y_exact = Leaf::func(x);
        if (leafmem[index].checkSolution(x, y_approx))
        {
            label leaf_index = leaf_manager.New();
            if (leaf_index != -1)
            {
                leafmem[leaf_index] = leafmem[index];
                leafmem[leaf_index].grow(x);//.EOA = fabs(x - leafmem[leaf_index].v);
                label node_index = leafmem[index].node_;

                if (nodemem[node_index].leafLeft_ == index)
                {
                    nodemem[node_index].leafLeft_ = leaf_index;
                }
                else
                {
                    nodemem[node_index].leafRight_ = leaf_index;
                }
                leaf_manager.Delete(index);
                success = true;
            }
        }
    }
    //write_lock.unlock();
    return success;
}
template<typename T, typename T2>
bool Foam::parallelISAT<T, T2>::retrieve(const inputType& x, outputType& y)
{
    //static label zzz = 0;
    /*if (!this->treeCheck())
    {
        FatalErrorInFunction
            << "broken tree retr in func 0\n"
            << abort(FatalError);
    }*/
    label index = binaryTreeSearch(x);
    if (index == -1)
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
    if (leafmem[index].inEOA(x))
    {
        //Pout << "inEOA" << endl;
        /*if (!this->treeCheck())
        {
            FatalErrorInFunction
                << "broken tree retr in func 1.5\n"
                << abort(FatalError);
        }*/
        leafmem[index].retrieve(x, y);


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
template<typename T, typename T2>
Foam::label Foam::parallelISAT<T, T2>::add
(
    const inputType& x
)
{
    return add(x, T::leafData::func(x));
}


template<typename T, typename T2>
bool Foam::parallelISAT<T, T2>::isFull()
{
    return size_leaf_ >= maxNLeafs_;
}
template<typename T, typename T2>
template<typename ...Args>
Foam::label Foam::parallelISAT<T, T2>::add
(
    const inputType& x, const outputType& y, Args ...args
)
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
        label node_index = node_manager.New();
        if (node_index == -1)
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
        label leaf_index = leaf_manager.New();
        if (leaf_index == -1)
        {
            growthOrAddFlag = 0;
            node_manager.Delete(node_index);
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
        nodemem[node_index].parent_ = -1;
        nodemem[node_index].nodeLeft_ = -1;
        nodemem[node_index].nodeRight_ = -1;
        nodemem[node_index].leafLeft_ = leaf_index;
        nodemem[node_index].leafRight_ = -1;
        nodemem[node_index].set();
        //Pout << "haha3,node_index=" << node_index << endl;
        leafmem[leaf_index].node_ = node_index;
        leafmem[leaf_index].set(x, y, args...);
        root_ = node_index;
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
    label index = binaryTreeSearch(x);
    //outputType y_approx = leafmem[index].retrieve(value);
    //Pout << "aa!" << leafmem[index].inEOA(value) << endl;
    if (leafmem[index].inEOA(x))
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
    if (leafmem[index].checkSolution(x, y))
    {
        /*if (zzzz2_ == 4 && x[0] < 1e-10)
        {
            Pout << "try grow"
                << x;
            //<< abort(FatalError);
        }*/

        //zz0++;
        if (leafmem[index].grow(x))
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
    label leaf_index = binaryTreeSearch(x, root_);

    if (leafmem[leaf_index].inEOA(x))
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
    label parentNode = leafmem[leaf_index].node_;

    label new_node_index = node_manager.New();
    if (new_node_index == -1)
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
    label new_leaf_index = leaf_manager.New();
    if (new_leaf_index == -1)
    {
        growthOrAddFlag = 0;
        node_manager.Delete(new_node_index);
        //if (!treeCheck())
        //{
        //    FatalErrorInFunction
        //        << "pos2,3"
        //        << abort(FatalError);
        //}
        write_lock.unlock();
        return growthOrAddFlag;
    }
    //Pout << "heihei3,new_leaf_index=" << new_leaf_index << endl;
    leafmem[new_leaf_index].set(x, y, args...);
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
    if (size_leaf_ > 1)
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
        return growthOrAddFlag;
    }
    else
    {
        root_ = new_node_index;
        node_manager.Delete(parentNode);
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


template<typename T, typename T2>
Foam::label Foam::parallelISAT<T, T2>::treeMin(label subTreeRoot)
{
    if (subTreeRoot != -1)
    {
        while (nodemem[subTreeRoot].nodeLeft_ != -1)
        {
            subTreeRoot = nodemem[subTreeRoot].nodeLeft_;
        }
        return nodemem[subTreeRoot].leafLeft_;
    }
    else
    {
        return -1;
    }
}

template<typename T, typename T2>
Foam::label Foam::parallelISAT<T, T2>::treeSuccessor(label x)
{
    if (size_leaf_ > 1)
    {
        if (x == nodemem[leafmem[x].node_].leafLeft_)
        {
            if (nodemem[leafmem[x].node_].nodeRight_ == -1)
            {
                return nodemem[leafmem[x].node_].leafRight_;
            }
            else
            {
                return treeMin(nodemem[leafmem[x].node_].nodeRight_);
            }
        }
        else if (x == nodemem[leafmem[x].node_].leafRight_)
        {
            label y = leafmem[x].node_;
            while ((nodemem[y].parent_ != -1))
            {
                if (y == nodemem[nodemem[y].parent_].nodeLeft_)
                {
                    if (nodemem[nodemem[y].parent_].nodeRight_ == -1)
                    {
                        return nodemem[nodemem[y].parent_].leafRight_;
                    }
                    else
                    {
                        return treeMin(nodemem[nodemem[y].parent_].nodeRight_);
                    }
                }
                y = nodemem[y].parent_;
            }
            // when we reach this point, y points to the root and
            // never entered in the if loop (coming from the right)
            // so we are at the tree maximum and there is no successor
            return -1;
        }
        else
        {
            FatalErrorInFunction
                << "inconsistent structure of the tree, no leaf and no node"
                << exit(FatalError);
            return -1;
        }
    }

    return -1;
}
template<typename T, typename T2>
Foam::label Foam::parallelISAT<T, T2>::leafPSibling(label x)
{
    if (size_leaf_ > 1)
    {
        if (x == nodemem[leafmem[x].node_].leafLeft_)
        {
            // x is on the left, return right side
            // might return nullptr if the right side is a node
            return nodemem[leafmem[x].node_].leafRight_;
        }
        else if (x == nodemem[leafmem[x].node_].leafRight_)
        {
            // x is on the right, return left side
            return nodemem[leafmem[x].node_].leafLeft_;
        }
        else
        {
            FatalErrorInFunction
                << "wrong addressing of the initial leaf"
                << exit(FatalError);
            return -1;
        }
    }
    // there is only one leaf attached to the root_, no sibling
    return -1;
}

template<typename T, typename T2>
Foam::label Foam::parallelISAT<T, T2>::nodeSibling_node(label y)
{
    if (nodemem[y].parent_ != -1)
    {
        if (y == nodemem[nodemem[y].parent_].nodeLeft_)
        {
            // y is on the left, return right side
            return nodemem[nodemem[y].parent_].nodeRight_;
        }
        else if (y == nodemem[nodemem[y].parent_].nodeRight_)
        {
            return nodemem[nodemem[y].parent_].nodeLeft_;
        }
        else
        {
            FatalErrorInFunction
                << "wrong addressing of the initial node"
                << exit(FatalError);
            return -1;
        }
    }
    return -1;
}

template<typename T, typename T2>
Foam::label Foam::parallelISAT<T, T2>::nodeSibling_leaf(label x)
{
    if (size_leaf_ > 1)
    {
        if (x == nodemem[leafmem[x].node_].leafLeft_)
        {
            // x is on the left, return right side
            return nodemem[leafmem[x].node_].nodeRight_;
        }
        else if (x == nodemem[leafmem[x].node_].leafRight_)
        {
            // x is on the right, return left side
            return nodemem[leafmem[x].node_].nodeLeft_;
        }
        else
        {
            FatalErrorInFunction
                << "wrong addressing of the initial leaf"
                << exit(FatalError);
            return -1;
        }
    }
    return -1;
}

template<typename T, typename T2>
void Foam::parallelISAT<T, T2>::transplant(label u, label v)
{
    if (v != -1)
    {
        // u is root_
        if (nodemem[u].parent_ == -1)
        {
            root_ = v;
        }
        // u is on the left of its parent
        else if (u == nodemem[nodemem[u].parent_].nodeLeft_)
        {
            nodemem[nodemem[u].parent_].nodeLeft_ = v;
        }
        // u is ont the right of its parent
        else if (u == nodemem[nodemem[u].parent_].nodeRight_)
        {
            nodemem[nodemem[u].parent_].nodeRight_ = v;
        }
        else
        {
            FatalErrorInFunction
                << "wrong addressing of the initial node"
                << exit(FatalError);
        }
        nodemem[v].parent_ = nodemem[u].parent_;
    }
    else
    {
        FatalErrorInFunction
            << "trying to transplant a nullptr node"
            << exit(FatalError);
    }
}

template<typename T, typename T2>
void Foam::parallelISAT<T, T2>::deleteLeaf(label index_phi0)
{
    if (size_leaf_ == 1) // only one point is stored
    {
        leaf_manager.Delete(index_phi0);
        node_manager.Delete(root_);
        root_ = -1;
        size_node--;
        //deleteDemandDrivenData(phi0);
        //deleteDemandDrivenData(root_);
    }
    else if (size_leaf_ > 1)
    {
        label z = leafmem[index_phi0].node_;
        label x;
        label siblingPhi0 = leafPSibling(index_phi0);

        if (siblingPhi0 != -1)// the sibling of phi0 is a chemPoint
        {
            // z was root (only two chemPoints in the tree)
            if (nodemem[z].parent_ == -1)
            {
                root_ = node_manager.New();
                if (root_ == -1)
                {
                    FatalErrorInFunction
                        << "run out of memory"
                        << exit(FatalError);
                }
                size_node++;
                nodemem[root_].nodeRight_ = -1;
                nodemem[root_].leafRight_ = -1;
                nodemem[root_].nodeLeft_ = -1;
                nodemem[root_].leafLeft_ = siblingPhi0;
                nodemem[root_].parent_ = -1;
                //root_->leafLeft() = siblingPhi0;
                leafmem[siblingPhi0].node_ = root_;
            }
            else if (z == nodemem[nodemem[z].parent_].nodeLeft_)
            {
                nodemem[nodemem[z].parent_].leafLeft_ = siblingPhi0;
                //z->parent()->leafLeft() = siblingPhi0;
                nodemem[nodemem[z].parent_].nodeLeft_ = -1;
                //z->parent()->nodeLeft() = nullptr;
                leafmem[siblingPhi0].node_ = nodemem[z].parent_;
                //siblingPhi0->node() = z->parent();
            }
            else if (z == nodemem[nodemem[z].parent_].nodeRight_)
            {
                nodemem[nodemem[z].parent_].leafRight_ = siblingPhi0;
                //z->parent()->leafRight() = siblingPhi0;
                nodemem[nodemem[z].parent_].nodeRight_ = -1;
                //z->parent()->nodeRight() = nullptr;
                leafmem[siblingPhi0].node_ = nodemem[z].parent_;
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
            x = nodeSibling_leaf(index_phi0);
            if (x != -1)
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


template<typename T, typename T2>
void Foam::parallelISAT<T, T2>::deleteAllNode(label subTreeRoot)
{
    if (subTreeRoot != -1)
    {
        deleteAllNode(nodemem[subTreeRoot].nodeLeft_);
        deleteAllNode(nodemem[subTreeRoot].nodeRight_);
        node_manager.Delete(subTreeRoot);
        size_node--;
        //deleteDemandDrivenData(subTreeRoot);
    }
}

template<typename T, typename T2>
Foam::label Foam::parallelISAT<T, T2>::depth(label subTreeRoot)
{
    // when we reach the leaf, we return 0
    if (subTreeRoot == -1)
    {
        return 0;
    }
    else
    {
        return
            1
            + max
            (
                depth(nodemem[subTreeRoot].nodeLeft_),
                depth(nodemem[subTreeRoot].nodeRight_)
            );
    }
}
template<typename T, typename T2>
bool Foam::parallelISAT<T, T2>::treeCheck(label index)
{
    bool leftnode = (nodemem[index].nodeLeft_ != -1);
    bool leftleaf = (nodemem[index].leafLeft_ != -1);
    if ((leftnode && leftleaf) || ((!leftnode) && (!leftleaf)))
    {
        FatalErrorInFunction
            << "p4 ";
        //<< abort(FatalError);
        return false;
    }
    bool rightnode = (nodemem[index].nodeRight_ != -1);
    bool rightleaf = (nodemem[index].leafRight_ != -1);
    if ((rightnode && rightleaf) || ((!rightnode) && (!rightleaf)))
    {
        FatalErrorInFunction
            << "p5 ";
        //<< abort(FatalError);
        return false;
    }
    if (leftnode)
    {
        if (nodemem[nodemem[index].nodeLeft_].parent_ != index)
        {
            FatalErrorInFunction
                << "p6 ";
            //<< abort(FatalError);

            return false;
        }
    }
    else
    {
        if (leafmem[nodemem[index].leafLeft_].node_ != index)

        {
            FatalErrorInFunction
                << "p7 ";
            //<< abort(FatalError);
            return false;
        }
    }
    if (rightnode)
    {
        if (nodemem[nodemem[index].nodeRight_].parent_ != index)
        {
            FatalErrorInFunction
                << "p8 ";
            //<< abort(FatalError);
            return false;
        }
    }
    else
    {
        if (leafmem[nodemem[index].leafRight_].node_ != index)
        {
            FatalErrorInFunction
                << "p9 ";
            //<< abort(FatalError);
            return false;
        }
    }
    if (leftnode)
    {
        if (treeCheck(nodemem[index].nodeLeft_) == false)
        {
            FatalErrorInFunction
                << "p10 ";
            //<< abort(FatalError);
            return false;
        }
    }
    if (rightnode)
    {
        if (treeCheck(nodemem[index].nodeRight_) == false)
        {
            FatalErrorInFunction
                << "p11 ";
            //<< abort(FatalError);
            return false;
        }
    }
    return true;
}

template<typename T, typename T2>
void Foam::parallelISAT<T, T2>::insertNode
(
    label phi0,
    label newNode
)
{
    if (phi0 == nodemem[leafmem[phi0].node_].leafRight_)// phi0 is on the right
    {
        nodemem[leafmem[phi0].node_].leafRight_ = -1;
        nodemem[leafmem[phi0].node_].nodeRight_ = newNode;
        //phi0->node()->leafRight() = nullptr;
        //phi0->node()->nodeRight() = newNode;
        return;
    }
    else if (phi0 == nodemem[leafmem[phi0].node_].leafLeft_)// phi0 is on the left
    {
        nodemem[leafmem[phi0].node_].leafLeft_ = -1;
        nodemem[leafmem[phi0].node_].nodeLeft_ = newNode;
        //phi0->node()->leafLeft() = nullptr;
        //phi0->node()->nodeLeft() = newNode;
        return;

    }

    // if we reach this point, there is an issue with the addressing
    FatalErrorInFunction
        << "trying to insert a node with a wrong pointer to a chemPoint"
        << exit(FatalError);
}


template<typename T, typename T2>
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

template<typename T, typename T2>
void Foam::parallelISAT<T, T2>::writeDot(label nodeindex, std::ofstream& fout)
{
    bool flag_left = false, flag_right = false;
    if (nodemem[nodeindex].nodeLeft_ != -1)
    {
        fout << "node" << nodemem[nodeindex].nodeLeft_ << "[label=\"" << nodemem[nodemem[nodeindex].nodeLeft_] << "\"];" << std::endl;
        fout << "node" << nodeindex << " -> " << "node" << nodemem[nodeindex].nodeLeft_ << ";" << std::endl;
        flag_left = true;
    }
    if (nodemem[nodeindex].leafLeft_ != -1)
    {
        fout << "leaf" << nodemem[nodeindex].leafLeft_ << "[label=\"" << leafmem[nodemem[nodeindex].leafLeft_] << "\"];" << std::endl;
        fout << "node" << nodeindex << " -> " << "leaf" << nodemem[nodeindex].leafLeft_ << ";" << std::endl;
    }
    fout << "M" << nodeindex << "[style=invis];" << std::endl;
    fout << "node" << nodeindex << " -> " << "M" << nodeindex << "[style=invis];" << std::endl;
    if (nodemem[nodeindex].nodeRight_ != -1)
    {
        fout << "node" << nodemem[nodeindex].nodeRight_ << "[label=\"" << nodemem[nodemem[nodeindex].nodeRight_] << "\"];" << std::endl;
        fout << "node" << nodeindex << " -> " << "node" << nodemem[nodeindex].nodeRight_ << ";" << std::endl;
        flag_right = true;
    }
    if (nodemem[nodeindex].leafRight_ != -1)
    {
        fout << "leaf" << nodemem[nodeindex].leafRight_ << "[label=\"" << leafmem[nodemem[nodeindex].leafRight_] << "\"];" << std::endl;
        fout << "node" << nodeindex << " -> " << "leaf" << nodemem[nodeindex].leafRight_ << ";" << std::endl;
    }
    if (flag_left)
        writeDot(nodemem[nodeindex].nodeLeft_, fout);
    if (flag_right)
        writeDot(nodemem[nodeindex].nodeRight_, fout);
    return;
}

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
