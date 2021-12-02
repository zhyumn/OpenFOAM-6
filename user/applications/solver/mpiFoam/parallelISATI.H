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
Foam::parallelISAT<T, T2>::parallelISAT(SUPstream::mpi_manager& manager_in, label nmem, SUPstream::mpi_sync& sync_in, Args ...args)
    :T2(args...),nmem_leaf_(nmem), nmem_node_(nmem), tolerance_(1e-1), leafmem(manager_in, nmem_leaf_),
    emptylist_leaf(manager_in, nmem_leaf_), nodemem(manager_in, nmem_node_),
    emptylist_node(manager_in, nmem_node_), pv_(manager_in),
    size_leaf(pv_().size_leaf), head_leaf(pv_().head_leaf),
    tail_leaf(pv_().tail_leaf), temp_tail_leaf(pv_().temp_tail_leaf),
    size_node(pv_().size_node), head_node(pv_().heaf_node),
    tail_node(pv_().tail_node), temp_tail_node(pv_().temp_tail_node),
    root_(pv_().root), write_lock(manager_in, sync_in),
    leaf_manager(emptylist_leaf, head_leaf, tail_leaf, temp_tail_leaf, nmem_leaf_),
    node_manager(emptylist_node, head_node, tail_node, temp_tail_node, nmem_node_)
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
    size_leaf = 0;
    //head_leaf = 0;
    //tail_leaf = nmem_leaf_ - 1;
    //temp_tail_leaf = nmem_leaf_ - 1;

    size_node = 0;
    //heaf_node = 0;
    //tail_node = nmem_node_ - 1;
    //temp_tail_node = nmem_node_ - 1;

    root_ = -1;

    //tleaf.node_ = &tnode;
    //tleaf.v = 1;
    //tnode.leaf_ = &tleaf;
    //tnode.v = 2;

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
    if (size_leaf > 1)
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
    else if (size_leaf == 1)
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
    write_lock.lock();
    bool success = false;
    label index = binaryTreeSearch(x);
    outputType y_approx;
    leafmem[index].retrieve(x, y_approx);
    if (!leafmem[index].inEOA(x))
    {
        outputType y_exact = Leaf::func(x);
        if (y_approx - y_exact <= tolerance_)
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
    write_lock.unlock();
    return success;
}
template<typename T, typename T2>
bool Foam::parallelISAT<T, T2>::retrieve(const inputType& x, outputType& y)
{
    label index = binaryTreeSearch(x);
    if (index == -1)
    {
        y = Leaf::func(x);
        //Pout << "add" << endl;
        add(x);
        return true;
    }
    outputType y_approx;
    leafmem[index].retrieve(x,y_approx);
    if (leafmem[index].inEOA(x))
    {
        //Pout << "inEOA" << endl;
        leafmem[index].retrieve(x,y);
        return true;
    }
    else
    {
        //Pout << "not inEOA" << endl;
        y = Leaf::func(x);
        if (fabs(y_approx - y) > tolerance_ || !grow(x))
        {
            //Pout << "add" << endl;
            add(x);
        }
        return true;
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
template<typename ...Args>
Foam::label Foam::parallelISAT<T, T2>::add
(
    const inputType& x, const outputType& y, Args ...args
)
{
    write_lock.lock();

    label growthOrAddFlag = 1;
    //Pout << "haha" << endl;
    //Pout << "root_=" << root_ << endl;
    if (size_leaf == 0)
    {
        //Pout << "haha1" << endl;
        label node_index = node_manager.New();
        if (node_index == -1)
        {
            write_lock.unlock();
            return growthOrAddFlag;
        }
        label leaf_index = leaf_manager.New();
        if (leaf_index == -1)
        {
            node_manager.Delete(node_index);
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
        leafmem[leaf_index].set(x, y);
        root_ = node_index;
        size_leaf++;
        size_node++;
        //Pout << value << "inserted!" << endl;
        write_lock.unlock();
        return growthOrAddFlag;
    }
    label index = binaryTreeSearch(x);
    //outputType y_approx = leafmem[index].retrieve(value);
    //Pout << "aa!" << leafmem[index].inEOA(value) << endl;
    if (leafmem[index].inEOA(x))
    {
        write_lock.unlock();
        return growthOrAddFlag;
    }

    //Pout << "heihei2,value=" << value << endl;
    //Todo  leaf_index can be obtain by arg
    label leaf_index = binaryTreeSearch(x, root_);

    if (leafmem[leaf_index].inEOA(x))
    {
        write_lock.unlock();
        return growthOrAddFlag;
    }
    label parentNode = leafmem[leaf_index].node_;

    label new_node_index = node_manager.New();
    if (new_node_index == -1)
    {
        write_lock.unlock();
        return growthOrAddFlag;
    }
    label new_leaf_index = leaf_manager.New();
    if (new_leaf_index == -1)
    {
        node_manager.Delete(new_node_index);
        write_lock.unlock();
        return growthOrAddFlag;
    }
    //Pout << "heihei3,new_leaf_index=" << new_leaf_index << endl;
    leafmem[new_leaf_index].set(x, y);
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
        return growthOrAddFlag;
    }
    else
    {
        root_ = new_node_index;
        node_manager.Delete(parentNode);
        size_leaf++;
        //Pout << value << "inserted!" << endl;
        write_lock.unlock();
        return growthOrAddFlag;
    }
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
