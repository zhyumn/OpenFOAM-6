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

#include "binaryTreeS2.H"
#include "SortableList.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template <class CompType, class ThermoType>
void Foam::chemistryTabulationMethodSs::binaryTreeS<CompType, ThermoType>::insertNode(
    chP *&phi0,
    bn *&newNode)
{
    if (phi0 == phi0->node()->leafRight()) // phi0 is on the right
    {
        phi0->node()->leafRight() = nullptr;
        phi0->node()->nodeRight() = newNode;
        return;
    }
    else if (phi0 == phi0->node()->leafLeft()) // phi0 is on the left
    {
        phi0->node()->leafLeft() = nullptr;
        phi0->node()->nodeLeft() = newNode;
        return;
    }

    // if we reach this point, there is an issue with the addressing
    FatalErrorInFunction
        << "trying to insert a node with a wrong pointer to a chemPoint"
        << exit(FatalError);
}

template <class CompType, class ThermoType>
bool Foam::chemistryTabulationMethodSs::binaryTreeS<CompType, ThermoType>::inSubTree(
    const scalarField &phiq,
    bn *y,
    chP *x)
{
    if ((n2ndSearch_ < max2ndSearch_) && (y != nullptr))
    {
        scalar vPhi = 0;
        const scalarField &v = y->v();
        const scalar a = y->a();
        // compute v*phi
        for (label i = 0; i < phiq.size(); i++)
        {
            vPhi += phiq[i] * v[i];
        }
        if (vPhi <= a) // on the left side of the node
        {
            if (y->nodeLeft() == nullptr) // left is a chemPoint
            {
                n2ndSearch_++;
                if (y->leafLeft()->inEOA(phiq))
                {
                    x = y->leafLeft();
                    return true;
                }
            }
            else // the left side is a node
            {
                if (inSubTree(phiq, y->nodeLeft(), x))
                {
                    return true;
                }
            }

            // not on the left side, try the right side
            if ((n2ndSearch_ < max2ndSearch_) && y->nodeRight() == nullptr)
            {
                n2ndSearch_++;
                // we reach the end of the subTree we can return the result
                if (y->leafRight()->inEOA(phiq))
                {
                    x = y->leafRight();
                    return true;
                }
                else
                {
                    x = nullptr;
                    return false;
                }
            }
            else // test for n2ndSearch is done in the call of inSubTree
            {
                return inSubTree(phiq, y->nodeRight(), x);
            }
        }
        else // on right side (symmetric of above)
        {
            if (y->nodeRight() == nullptr)
            {
                n2ndSearch_++;
                if (y->leafRight()->inEOA(phiq))
                {
                    return true;
                }
            }
            else // the right side is a node
            {
                if (inSubTree(phiq, y->nodeRight(), x))
                {
                    x = y->leafRight();
                    return true;
                }
            }
            // if we reach this point, the retrieve has
            // failed on the right side, explore the left side
            if ((n2ndSearch_ < max2ndSearch_) && y->nodeLeft() == nullptr)
            {
                n2ndSearch_++;
                if (y->leafLeft()->inEOA(phiq))
                {
                    x = y->leafLeft();
                    return true;
                }
                else
                {
                    x = nullptr;
                    return false;
                }
            }
            else
            {
                return inSubTree(phiq, y->nodeLeft(), x);
            }
        }
    }
    else
    {
        return false;
    }
}

template <class CompType, class ThermoType>
void Foam::chemistryTabulationMethodSs::binaryTreeS<CompType, ThermoType>::deleteSubTree(bn *subTreeRoot)
{
    if (subTreeRoot != nullptr)
    {
        deleteDemandDrivenData(subTreeRoot->leafLeft());
        deleteDemandDrivenData(subTreeRoot->leafRight());
        deleteSubTree(subTreeRoot->nodeLeft());
        deleteSubTree(subTreeRoot->nodeRight());
        deleteDemandDrivenData(subTreeRoot);
    }
}

template <class CompType, class ThermoType>
void Foam::chemistryTabulationMethodSs::binaryTreeS<CompType, ThermoType>::transplant(bn *u, bn *v)
{
    if (v != nullptr)
    {
        // u is root_
        if (u->parent() == nullptr)
        {
            root_ = v;
        }
        // u is on the left of its parent
        else if (u == u->parent()->nodeLeft())
        {
            u->parent()->nodeLeft() = v;
        }
        // u is ont the right of its parent
        else if (u == u->parent()->nodeRight())
        {
            u->parent()->nodeRight() = v;
        }
        else
        {
            FatalErrorInFunction
                << "wrong addressing of the initial node"
                << exit(FatalError);
        }
        v->parent() = u->parent();
    }
    else
    {
        FatalErrorInFunction
            << "trying to transplant a nullptr node"
            << exit(FatalError);
    }
}

template <class CompType, class ThermoType>
Foam::chemPointISATS<CompType, ThermoType> *
Foam::chemistryTabulationMethodSs::binaryTreeS<CompType, ThermoType>::chemPSibling(bn *y)
{
    if (y->parent() != nullptr)
    {
        if (y == y->parent()->nodeLeft()) // y is on the left, return right side
        {
            // might return nullptr if the right side is a node
            return y->parent()->leafRight();
        }
        else if (y == y->parent()->nodeRight())
        {
            return y->parent()->leafLeft();
        }
        else
        {
            FatalErrorInFunction
                << "wrong addressing of the initial node"
                << exit(FatalError);
            return nullptr;
        }
    }

    // the binaryNodeS is root_ and has no sibling
    return nullptr;
}

template <class CompType, class ThermoType>
Foam::chemPointISATS<CompType, ThermoType> *
Foam::chemistryTabulationMethodSs::binaryTreeS<CompType, ThermoType>::chemPSibling(chP *x)
{
    if (size_ > 1)
    {
        if (x == x->node()->leafLeft())
        {
            // x is on the left, return right side
            // might return nullptr if the right side is a node
            return x->node()->leafRight();
        }
        else if (x == x->node()->leafRight())
        {
            // x is on the right, return left side
            return x->node()->leafLeft();
        }
        else
        {
            FatalErrorInFunction
                << "wrong addressing of the initial leaf"
                << exit(FatalError);
            return nullptr;
        }
    }
    // there is only one leaf attached to the root_, no sibling
    return nullptr;
}

template <class CompType, class ThermoType>
Foam::binaryNodeS<CompType, ThermoType> *
Foam::chemistryTabulationMethodSs::binaryTreeS<CompType, ThermoType>::nodeSibling(bn *y)
{
    if (y->parent() != nullptr)
    {
        if (y == y->parent()->nodeLeft())
        {
            // y is on the left, return right side
            return y->parent()->nodeRight();
        }
        else if (y == y->parent()->nodeRight())
        {
            return y->parent()->nodeLeft();
        }
        else
        {
            FatalErrorInFunction
                << "wrong addressing of the initial node"
                << exit(FatalError);
            return nullptr;
        }
    }
    return nullptr;
}

template <class CompType, class ThermoType>
Foam::binaryNodeS<CompType, ThermoType> *
Foam::chemistryTabulationMethodSs::binaryTreeS<CompType, ThermoType>::nodeSibling(chP *x)
{
    if (size_ > 1)
    {
        if (x == x->node()->leafLeft())
        {
            // x is on the left, return right side
            return x->node()->nodeRight();
        }
        else if (x == x->node()->leafRight())
        {
            // x is on the right, return left side
            return x->node()->nodeLeft();
        }
        else
        {
            FatalErrorInFunction
                << "wrong addressing of the initial leaf"
                << exit(FatalError);
            return nullptr;
        }
    }
    return nullptr;
}

template <class CompType, class ThermoType>
void Foam::chemistryTabulationMethodSs::binaryTreeS<CompType, ThermoType>::deleteAllNode(bn *subTreeRoot)
{
    if (subTreeRoot != nullptr)
    {
        deleteAllNode(subTreeRoot->nodeLeft());
        deleteAllNode(subTreeRoot->nodeRight());
        deleteDemandDrivenData(subTreeRoot);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class CompType, class ThermoType>
Foam::chemistryTabulationMethodSs::binaryTreeS<CompType, ThermoType>::binaryTreeS(
    STDACChemistryModel<CompType, ThermoType> &chemistry,
    dictionary coeffsDict, ISAT_manager<CompType, ThermoType> &ISATmanager_in)
    : chemistry_(chemistry),
      ISATmanager(ISATmanager_in),
      root_(nullptr),
      maxNLeafs_(readLabel(coeffsDict.lookup("maxNLeafs"))),
      size_(0),
      n2ndSearch_(0),
      max2ndSearch_(coeffsDict.lookupOrDefault("max2ndSearch", 0)),
      nRetrieved_(0),
      nAdd_(0),
      cleaningRequired_(false),
      coeffsDict_(coeffsDict)

{
    tolerance_ = ISATmanager.tolerance();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class CompType, class ThermoType>
Foam::label Foam::chemistryTabulationMethodSs::binaryTreeS<CompType, ThermoType>::depth(bn *subTreeRoot)
{
    // when we reach the leaf, we return 0
    if (subTreeRoot == nullptr)
    {
        return 0;
    }
    else
    {
        return 1 + max(
                       depth(subTreeRoot->nodeLeft()),
                       depth(subTreeRoot->nodeRight()));
    }
}

template <class CompType, class ThermoType>
void Foam::chemistryTabulationMethodSs::binaryTreeS<CompType, ThermoType>::insertNewLeaf(
    const scalarField &phiq,
    const scalarField &Rphiq,
    const scalarSquareMatrix &A,
    const scalarField &scaleFactor,
    const scalar &epsTol,
    const label nCols,
    chP *&phi0)
{
    if (size_ == 0) // no points are stored
    {
        // create an empty binary node and point root_ to it
        root_ = new bn();
        // create the new chemPoint which holds the composition point
        // phiq and the data to initialize the EOA
        chP *newChemPoint =
            new chP(
                chemistry_,
                phiq,
                Rphiq,
                A,
                scaleFactor,
                epsTol,
                nCols,
                coeffsDict_,
                root_);
        root_->leafLeft() = newChemPoint;
    }
    else // at least one point stored
    {
        // no reference chemPoint, a BT search is required
        if (phi0 == nullptr)
        {
            binaryTreeSearch(phiq, root_, phi0);
        }
        // access to the parent node of the chemPoint
        bn *parentNode = phi0->node();

        // create the new chemPoint which holds the composition point
        // phiq and the data to initialize the EOA
        chP *newChemPoint =
            new chP(
                chemistry_,
                phiq,
                Rphiq,
                A,
                scaleFactor,
                epsTol,
                nCols,
                coeffsDict_);
        // insert new node on the parent node in the position of the
        // previously stored leaf (phi0)
        // the new node contains phi0 on the left and phiq on the right
        // the hyper plane is computed in the binaryNodeS constructor
        bn *newNode;
        if (size_ > 1)
        {
            newNode = new bn(phi0, newChemPoint, parentNode);
            // make the parent of phi0 point to the newly created node
            insertNode(phi0, newNode);
        }
        else // size_ == 1 (because not equal to 0)
        {
            // when size is 1, the binaryNodeS is without hyperplane
            deleteDemandDrivenData(root_);
            newNode = new bn(phi0, newChemPoint, nullptr);
            root_ = newNode;
        }

        phi0->node() = newNode;
        newChemPoint->node() = newNode;
    }
    size_++;
}

template <class CompType, class ThermoType>
void Foam::chemistryTabulationMethodSs::binaryTreeS<CompType, ThermoType>::binaryTreeSearch(
    const scalarField &phiq,
    bn *node,
    chP *&nearest)
{
    if (size_ > 1)
    {
        scalar vPhi = 0.0;
        const scalarField &v = node->v();
        const scalar &a = node->a();
        // compute v*phi
        for (label i = 0; i < phiq.size(); i++)
            vPhi += phiq[i] * v[i];

        if (vPhi > a) // on right side (side of the newly added point)
        {
            if (node->nodeRight() != nullptr)
            {
                binaryTreeSearch(phiq, node->nodeRight(), nearest);
            }
            else // the terminal node is reached, store leaf on the right
            {
                nearest = node->leafRight();
                return;
            }
        }
        else // on left side (side of the previously stored point)
        {
            if (node->nodeLeft() != nullptr)
            {
                binaryTreeSearch(phiq, node->nodeLeft(), nearest);
            }
            else // the terminal node is reached, return element on right
            {
                nearest = node->leafLeft();
                return;
            }
        }
    }
    // only one point stored (left element of the root)
    else if (size_ == 1)
    {
        nearest = root_->leafLeft();
    }
    else // no point stored
    {
        nearest = nullptr;
    }
}

template <class CompType, class ThermoType>
bool Foam::chemistryTabulationMethodSs::binaryTreeS<CompType, ThermoType>::secondaryBTSearch(
    const scalarField &phiq,
    chP *&x)
{
    // initialize n2ndSearch_
    n2ndSearch_ = 0;
    if ((n2ndSearch_ < max2ndSearch_) && (size_ > 1))
    {
        chP *xS = chemPSibling(x);
        if (xS != nullptr)
        {
            n2ndSearch_++;
            if (xS->inEOA(phiq))
            {
                x = xS;
                return true;
            }
        }
        else if (inSubTree(phiq, nodeSibling(x), x))
        {
            return true;
        }
        // if we reach this point, no leafs were found at this depth or lower
        // we move upward in the tree
        bn *y = x->node();
        while ((y->parent() != nullptr) && (n2ndSearch_ < max2ndSearch_))
        {
            xS = chemPSibling(y);
            if (xS != nullptr)
            {
                n2ndSearch_++;
                if (xS->inEOA(phiq))
                {
                    x = xS;
                    return true;
                }
            }
            else if (inSubTree(phiq, nodeSibling(y), x))
            {
                return true;
            }
            y = y->parent();
        }
        // if we reach this point it is either because
        // we did not find another covering EOA in the entire tree or
        // we reach the maximum number of secondary search
        return false;
    }
    else
    {
        return false;
    }
}

template <class CompType, class ThermoType>
void Foam::chemistryTabulationMethodSs::binaryTreeS<CompType, ThermoType>::deleteLeaf(chP *&phi0)
{
    if (size_ == 1) // only one point is stored
    {
        deleteDemandDrivenData(phi0);
        deleteDemandDrivenData(root_);
    }
    else if (size_ > 1)
    {
        bn *z = phi0->node();
        bn *x;
        chP *siblingPhi0 = chemPSibling(phi0);

        if (siblingPhi0 != nullptr) // the sibling of phi0 is a chemPoint
        {
            // z was root (only two chemPoints in the tree)
            if (z->parent() == nullptr)
            {
                root_ = new bn();
                root_->leafLeft() = siblingPhi0;
                siblingPhi0->node() = root_;
            }
            else if (z == z->parent()->nodeLeft())
            {
                z->parent()->leafLeft() = siblingPhi0;
                z->parent()->nodeLeft() = nullptr;
                siblingPhi0->node() = z->parent();
            }
            else if (z == z->parent()->nodeRight())
            {
                z->parent()->leafRight() = siblingPhi0;
                z->parent()->nodeRight() = nullptr;
                siblingPhi0->node() = z->parent();
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
            x = nodeSibling(phi0);
            if (x != nullptr)
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
        deleteDemandDrivenData(phi0);
        deleteDemandDrivenData(z);
    }
    size_--;
}

template <class CompType, class ThermoType>
void Foam::chemistryTabulationMethodSs::binaryTreeS<CompType, ThermoType>::balance()
{
    scalarField mean(chemistry_.nEqns(), 0.0);

    //1) walk through the entire tree by starting with the tree's most left
    // chemPoint
    chP *x = treeMin();
    List<chP *> chemPoints(size_);
    label chPi = 0;
    //2) compute the mean composition
    while (x != nullptr)
    {
        const scalarField &phij = x->phi();
        mean += phij;
        chemPoints[chPi++] = x;
        x = treeSuccessor(x);
    }
    mean /= size_;

    //3) compute the variance for each space direction
    List<scalar> variance(chemistry_.nEqns(), 0.0);
    forAll(chemPoints, j)
    {
        const scalarField &phij = chemPoints[j]->phi();
        forAll(variance, vi)
        {
            variance[vi] += sqr(phij[vi] - mean[vi]);
        }
    }

    //4) analyze what is the direction of the maximal variance
    scalar maxVariance(-1.0);
    label maxDir(-1);
    forAll(variance, vi)
    {
        if (maxVariance < variance[vi])
        {
            maxVariance = variance[vi];
            maxDir = vi;
        }
    }
    // maxDir indicates the direction of maximum variance
    // we create the new root node by taking the two extreme points
    // in this direction if these extreme points were not deleted in the
    // cleaning that come before the balance function they are still important
    // and the tree should therefore take them into account
    SortableList<scalar> phiMaxDir(chemPoints.size(), 0.0);
    forAll(chemPoints, j)
    {
        phiMaxDir[j] = chemPoints[j]->phi()[maxDir];
    }

    phiMaxDir.sort();
    // delete reference to all node since the tree is reshaped
    deleteAllNode();
    root_ = nullptr;

    // add the node for the two extremum
    bn *newNode = new bn(
        chemPoints[phiMaxDir.indices()[0]],
        chemPoints[phiMaxDir.indices()[phiMaxDir.size() - 1]],
        nullptr);
    root_ = newNode;

    chemPoints[phiMaxDir.indices()[0]]->node() = newNode;
    chemPoints[phiMaxDir.indices()[phiMaxDir.size() - 1]]->node() = newNode;

    for (label cpi = 1; cpi < chemPoints.size() - 1; cpi++)
    {
        chP *phi0;
        binaryTreeSearch(
            chemPoints[phiMaxDir.indices()[cpi]]->phi(),
            root_,
            phi0);
        // add the chemPoint
        bn *nodeToAdd =
            new bn(phi0, chemPoints[phiMaxDir.indices()[cpi]], phi0->node());
        // make the parent of phi0 point to the newly created node
        insertNode(phi0, nodeToAdd);
        phi0->node() = nodeToAdd;
        chemPoints[phiMaxDir.indices()[cpi]]->node() = nodeToAdd;
    }
}

template <class CompType, class ThermoType>
Foam::chemPointISATS<CompType, ThermoType> *
Foam::chemistryTabulationMethodSs::binaryTreeS<CompType, ThermoType>::treeMin(bn *subTreeRoot)
{
    if (subTreeRoot != nullptr)
    {
        while (subTreeRoot->nodeLeft() != nullptr)
        {
            subTreeRoot = subTreeRoot->nodeLeft();
        }
        return subTreeRoot->leafLeft();
    }
    else
    {
        return nullptr;
    }
}

template <class CompType, class ThermoType>
Foam::chemPointISATS<CompType, ThermoType> *
Foam::chemistryTabulationMethodSs::binaryTreeS<CompType, ThermoType>::treeSuccessor(chP *x)
{
    if (size_ > 1)
    {
        if (x == x->node()->leafLeft())
        {
            if (x->node()->nodeRight() == nullptr)
            {
                return x->node()->leafRight();
            }
            else
            {
                return treeMin(x->node()->nodeRight());
            }
        }
        else if (x == x->node()->leafRight())
        {
            bn *y = x->node();
            while ((y->parent() != nullptr))
            {
                if (y == y->parent()->nodeLeft())
                {
                    if (y->parent()->nodeRight() == nullptr)
                    {
                        return y->parent()->leafRight();
                    }
                    else
                    {
                        return treeMin(y->parent()->nodeRight());
                    }
                }
                y = y->parent();
            }
            // when we reach this point, y points to the root and
            // never entered in the if loop (coming from the right)
            // so we are at the tree maximum and there is no successor
            return nullptr;
        }
        else
        {
            FatalErrorInFunction
                << "inconsistent structure of the tree, no leaf and no node"
                << exit(FatalError);
            return nullptr;
        }
    }

    return nullptr;
}

template <class CompType, class ThermoType>
void Foam::chemistryTabulationMethodSs::binaryTreeS<CompType, ThermoType>::clear()
{
    // Recursively delete the element in the subTree
    deleteSubTree();

    // Reset root node (should already be nullptr)
    root_ = nullptr;

    // Reset size_
    size_ = 0;
}

template <class CompType, class ThermoType>
bool Foam::chemistryTabulationMethodSs::binaryTreeS<CompType, ThermoType>::isFull()
{
    return size_ >= maxNLeafs_;
}

template <class CompType, class ThermoType>
void Foam::chemistryTabulationMethodSs::binaryTreeS<CompType, ThermoType>::resetNumRetrieve()
{
    // Has to go along each chP of the tree
    if (size_ > 0)
    {
        // First finds the first leaf
        chP *chP0 = treeMin();
        chP0->resetNumRetrieve();

        // Then go to each successor
        chP *nextChP = treeSuccessor(chP0);
        while (nextChP != nullptr)
        {
            nextChP->resetNumRetrieve();
            nextChP = treeSuccessor(nextChP);
        }
    }
}

template <class CompType, class ThermoType>
bool Foam::chemistryTabulationMethodSs::binaryTreeS<CompType, ThermoType>::retrieve(
    const scalarField &phiq,
    scalarField &Rphiq)
{
    bool retrieved(false);
    chemPointISATS<CompType, ThermoType> *phi0;
    // If the tree is not empty
    if (size())
    {
        binaryTreeSearch(phiq, root_, phi0);

        // lastSearch keeps track of the chemPoint we obtain by the regular
        // binary tree search
        //lastSearch_ = phi0;
        if (phi0->inEOA(phiq))
        {
            retrieved = true;
        }
        // After a successful secondarySearch, phi0 store a pointer to the
        // found chemPoint
        else if (secondaryBTSearch(phiq, phi0))
        {
            retrieved = true;
        }
        //else if (MRURetrieve_)
        //{
        //    typename SLList<
        //        chemPointISATS<CompType, ThermoType> *>::iterator iter = MRUList_.begin();
        //    for (; iter != MRUList_.end(); ++iter)
        //    {
        //        phi0 = iter();
        //        if (phi0->inEOA(phiq))
        //        {
        //            retrieved = true;
        //            break;
        //        }
        //    }
        //}
    }
    // The tree is empty, retrieved is still false
    //else
    //{
    // There is no chempoints that we can try to grow
    //lastSearch_ = nullptr;
    //}

    if (retrieved)
    {
        phi0->increaseNumRetrieve();
        scalar elapsedTimeSteps =
            this->chemistry_.timeSteps() - phi0->timeTag();

        // Raise a flag when the chemPoint has been used more than the allowed
        // number of time steps
        if (elapsedTimeSteps > chPMaxLifeTime_ && !phi0->toRemove())
        {
            cleaningRequired_ = true;
            phi0->toRemove() = true;
        }
        //lastSearch_->lastTimeUsed() = this->chemistry_.timeSteps();
        //addToMRU(phi0);
        //calcNewC(phi0, phiq, Rphiq);
        phi0->retrieve(phiq, Rphiq);
        nRetrieved_++;
        return true;
    }
    else
    {
        // This point is reached when every retrieve trials have failed
        // or if the tree is empty
        return false;
    }
}

template <class CompType, class ThermoType>
Foam::label Foam::chemistryTabulationMethodSs::binaryTreeS<CompType, ThermoType>::add(
    const scalarField &phiq,
    const scalarField &Rphiq,
    const scalar rho,
    const scalar deltaT)
{
    label growthOrAddFlag = 1;
    // If lastSearch_ holds a valid pointer to a chemPoint AND the growPoints_
    // option is on, the code first tries to grow the point hold by lastSearch_
    //if (lastSearch_ && growPoints_)
    //{
    chemPointISATS<CompType, ThermoType> *phi0 = nullptr;
    binaryTreeSearch(phiq, root_, phi0);
    if (grow(phi0, phiq, Rphiq))
    {
        nGrowth_++;
        growthOrAddFlag = 0;
        //addToMRU(lastSearch_);
        //the structure of the tree is not modified, return false
        return growthOrAddFlag;
    }
    //}

    // If the code reach this point, it is either because lastSearch_ is not
    // valid, OR because growPoints_ is not on, OR because the grow operation
    // has failed. In the three cases, a new point is added to the tree.
    if (isFull())
    {
        if (!cleanAndBalance())
        {
            clear();
        }
    }
    /* if (chemisTree_.isFull())
                {
                    // If cleanAndBalance operation do not result in a reduction of the tree
                    // size, the last possibility is to delete completely the tree.
                    // It can be partially rebuild with the MRU list if this is used.
                    if (!cleanAndBalance())
                    {
                        DynamicList<chemPointISATS<CompType, ThermoType> *> tempList;
                        if (maxMRUSize_ > 0)
                        {
                            // Create a copy of each chemPointISATS of the MRUList_ before
                            // they are deleted
                            typename SLList<
                                chemPointISATS<CompType, ThermoType> *>::iterator iter = MRUList_.begin();
                            for (; iter != MRUList_.end(); ++iter)
                            {
                                tempList.append(
                                    new chemPointISATS<CompType, ThermoType>(*iter()));
                            }
                        }
                        chemisTree_.clear();

                        // Pointers to chemPoint are not valid anymore, clear the list
                        //MRUList_.clear();

                        // Construct the tree without giving a reference to attach to it
                        // since the structure has been completely discarded
                        chemPointISATS<CompType, ThermoType> *nulPhi = 0;
                        forAll(tempList, i)
                        {
                            chemisTree_.insertNewLeaf(
                                tempList[i]->phi(),
                                tempList[i]->Rphi(),
                                tempList[i]->A(),
                                scaleFactor_,
                                this->tolerance(),
                                scaleFactor_.size(),
                                nulPhi);
                            deleteDemandDrivenData(tempList[i]);
                        }
                    }

                    // The structure has been changed, it will force the binary tree to
                    // perform a new search and find the most appropriate point still stored
                    //lastSearch_ = nullptr;
                } */

    // Compute the A matrix needed to store the chemPoint.
    label ASize = this->chemistry_.nEqns() + ISATmanager.nAdditionalEqns_ - 2;
    scalarSquareMatrix A(ASize, Zero);
    ISATmanager.computeA(A, Rphiq, rho, deltaT);
    chemPointISATS<CompType, ThermoType> *NP = nullptr;
    insertNewLeaf(
        phiq,
        Rphiq,
        A,
        ISATmanager.scaleFactor_,
        tolerance_,
        ISATmanager.scaleFactor_.size(),
        NP // lastSearch_ may be nullptr (handled by binaryTreeS)
    );
    /*                 if (lastSearch_ != nullptr)
                {
                    addToMRU(lastSearch_);
                } */
    nAdd_++;

    return growthOrAddFlag;
}
// ************************************************************************* //
template <class CompType, class ThermoType>
bool Foam::chemistryTabulationMethodSs::binaryTreeS<CompType, ThermoType>::grow(
    chemPointISATS<CompType, ThermoType> *phi0,
    const scalarField &phiq,
    const scalarField &Rphiq)
{
    // If the pointer to the chemPoint is nullptr, the function stops
    if (!phi0)
    {
        return false;
    }

    // Raise a flag when the chemPoint used has been grown more than the
    // allowed number of time
    if (phi0->nGrowth() > maxGrowth_)
    {
        cleaningRequired_ = true;
        phi0->toRemove() = true;
        return false;
    }

    // If the solution RphiQ is still within the tolerance we try to grow it
    // in some cases this might result in a failure and the grow function of
    // the chemPoint returns false
    if (phi0->checkSolution(phiq, Rphiq))
    {
        return phi0->grow(phiq);
    }
    // The actual solution and the approximation given by ISAT are too different
    else
    {
        return false;
    }
}

template <class CompType, class ThermoType>
bool Foam::chemistryTabulationMethodSs::binaryTreeS<CompType, ThermoType>::cleanAndBalance()
{

    bool treeModified(false);

    {
        chemPointISATS<CompType, ThermoType> *x = treeMin();
        while (x != nullptr)
        {
            chemPointISATS<CompType, ThermoType> *xtmp =
                treeSuccessor(x);

            scalar elapsedTimeSteps = this->chemistry_.timeSteps() - x->timeTag();

            if ((elapsedTimeSteps > chPMaxLifeTime_) || (x->nGrowth() > maxGrowth_))
            {
                deleteLeaf(x);
                treeModified = true;
            }
            x = xtmp;
        }

        if (
            size() > minBalanceThreshold_ && depth() >
                                                 maxDepthFactor_ * log(scalar(size())) / log(2.0))
        {
            balance();
            treeModified = true;
        }
    }

    return (treeModified && isFull());
}