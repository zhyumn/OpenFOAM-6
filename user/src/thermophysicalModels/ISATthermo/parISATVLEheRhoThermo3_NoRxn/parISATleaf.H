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
    Foam::ISATleaf

Description
    Leaf of the binary tree.
    The chemPoint stores the composition 'phi', the mapping of this
    composition Rphi, the mapping gradient matrix A and the matrix describing
    the Ellipsoid Of Accuracy (EOA).

  1)When the chemPoint is created the region of accuracy is approximated by
    an ellipsoid E centered in 'phi' (obtained with the constant):
        E = {x| ||L^T.(x-phi)|| <= 1},
    with x a point in the composition space and L^T the transpose of an upper
    triangular matrix describing the EOA (see below: "Computation of L" ).

  2)To RETRIEVE the mapping from the chemPoint phi, the query point phiq has to
    be in the EOA of phi. It follows that, dphi=phiq-phi and to test if phiq
    is in the ellipsoid there are two methods. First, compare r=||dphi|| with
    rmin and rmax. If r < rmin, phiq is in the EOA. If r > rmax, phiq is out of
    the EOA. This operations is O(completeSpaceSize) and is performed first.
    If rmin < r < rmax, then the second method is used:
        ||L^T.dphi|| <= 1 to be in the EOA.

    If phiq is in the EOA, Rphiq is obtained by linear interpolation:
        Rphiq= Rphi + A.dphi.

  3)If phiq is not in the EOA, then the mapping is computed. But as the EOA
    is a conservative approximation of the region of accuracy surrounding the
    point phi, we could expand it by comparing the computed results with the
    one obtained by linear interpolation. The error epsGrow is calculated:
        epsGrow = ||B.(dR - dRl)||,
    with dR = Rphiq - Rphi, dRl = A.dphi and B the diagonal scale factor
    matrix.
    If epsGrow <= tolerance, the EOA is too conservative and a GROW is perforned
    otherwise, the newly computed mapping is associated to the initial
    composition and added to the tree.

  4)To GROW the EOA, we expand it to include the previous EOA and the query
    point phiq. The rank-one matrix method is used. The EOA is transformed
    to a hypersphere centered at the origin. Then it is expanded to include
    the transformed point phiq' on its boundary. Then the inverse transformation
    give the modified matrix L' (see below: "Grow the EOA").


  Computation of L :
    In [1], the EOA of the constant approximation is given by
        E = {x| ||B.A/tolerance.(x-phi)|| <= 1},
    with B a scale factor diagonal matrix, A the mapping gradient matrix and
    tolerance the absolute tolerance. If we take the QR decomposition of
    (B.A)/tolerance= Q.R, with Q an orthogonal matrix and R an upper triangular
    matrix such that the EOA is described by
    (phiq-phi0)^T.R^T.R.(phiq-phi0) <= 1
    L^T = R, both Cholesky decomposition of A^T.B^T.B.A/tolerance^2
    This representation of the ellipsoid is used in [2] and in order to avoid
    large value of semi-axe length in certain direction, a Singular Value
    Decomposition (SVD) is performed on the L matrix:
        L = UDV^T,
    with the orthogonal matrix U giving the directions of the principal axes
    and 1/di the inverse of the element of the diagonal matrix D giving the
    length of the principal semi-axes. To avoid very large value of those
    length,
    di' = max(di, 1/(alphaEOA*sqrt(tolerance))), with alphaEOA = 0.1 (see [2])
    di' = max(di, 1/2), see [1]. The latter will be used in this implementation.
    And L' = UD'V^T, with D' the diagonal matrix with the modified di'.

  Grow the EOA :
    More details about the minimum-volume ellipsoid covering an ellispoid E and
    a point p are found in [3]. Here is the main steps to obtain the modified
    matrix L' describind the new ellipsoid.
        1) calculate the point p' in the transformed space :
            p' = L^T.(p-phi)
        2) compute the rank-one decomposition:
            G = I + gamma.p'.p'^T,
           with gamma = (1/|p'|-1)*1/|p'|^2
        3) compute L':
            L' = L.G.

    References:
    \verbatim
        [1] Pope, S. B. (1997).
        Computationally efficient implementation of combustion chemistry using
        in situ adaptive tabulation.
        Combustion Theory and Modelling, 1, 41-63.

        [2] Lu, L., & Pope, S. B. (2009).
        An improved algorithm for in situ adaptive tabulation.
        Journal of Computational Physics, 228(2), 361-386.

        [3] Pope, S. B. (2008).
        Algorithms for ellipsoids.
        Cornell University Report No. FDA, 08-01.
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef PARISATleaf_H
#define PARISATleaf_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#include "scalarList.H"
#include "scalarMatrices.H"
#include "messageStream.H"
#include "List2D.H"
#include "slab.H"
#include "parData.H"
namespace Foam
{

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    class parISATNode;

    //template<class CompType, class ThermoType>
    //class TDACChemistryModel;

    /*---------------------------------------------------------------------------*\
                       Class ISATleaf Declaration
\*---------------------------------------------------------------------------*/

    class parISATleaf
    {

    public:
        SharedPointer<parISATNode> node_;
        SList<scalar> value_;
        SList<scalar> data_;
        SMatrix<scalar> A_;
        SMatrix<scalar> EOA_;

        label lastUsed;
        label numRetrieve_;
        label SleafN;
        //typename List2D<ISATleaf>::Node* pTimeTagList_;
        // Constructors

        //- Construct from components
        parISATleaf(int n_in, int n_out, const scalarList &v);
        parISATleaf(int n_in, int n_out, const scalarList &v, SharedPointer<parISATNode> &node);
        parISATleaf(int n_in, int n_out, const scalarList &v, SharedPointer<parISATNode> &node, const scalarList &data_in);
        void init(int n_in, int n_out);
        void reuse();

        void set(const parISATleaf &x)
        {
            lastUsed = 0;
            node_ = x.node_;
            forAll(value_, i)
            {
                value_[i] = x.value_[i];
            }
            forAll(data_, i)
            {
                data_[i] = x.data_[i];
            }
            for (int i = 0; i < EOA_.n_in_; i++)
                for (int j = 0; j < EOA_.n_out_; j++)
                {
                    EOA_(i, j) = x.EOA_(i, j);
                }
            for (int i = 0; i < A_.n_in_; i++)
                for (int j = 0; j < A_.n_out_; j++)
                {
                    A_(i, j) = x.A_(i, j);
                }
        }

        void set(const scalarList &v, SharedPointer<parISATNode> &node, const scalarList &data_in)
        {
            lastUsed = 0;
            node_ = node;
            forAll(value_, i)
            {
                value_[i] = v[i];
            }
            forAll(data_, i)
            {
                data_[i] = data_in[i];
            }
        }
        void Delete()
        {
            value_.Delete();
            data_.Delete();
            A_.Delete();
            EOA_.Delete();
        }

        void setNext(label offset)
        {
            node_.offset = offset;
        }
        label next()
        {
            return node_.offset;
        }

        /*         inline scalarList &value()
        {
            return value_;
        }
        inline const scalarList &value() const
        {
            return value_;
        }
        inline scalarList &data()
        {
            return data_;
        }
        inline const scalarList &data() const
        {
            return data_;
        }

        inline scalarRectangularMatrix &A()
        {
            return A_;
        }

        inline const scalarRectangularMatrix &A() const
        {
            return A_;
        }
        inline scalarRectangularMatrix &EOA()
        {
            return EOA_;
        }; */

        /*         inline typename List2D<ISATleaf>::Node *&pTimeTagList()
        {
            return pTimeTagList_;
        }; */

        bool inEOA(const scalarList &point, const scalarRectangularMatrix &scaleIn);

        void eval(const scalarList &value, scalarList &ret);
        void grow(const scalarList &point, const scalarRectangularMatrix &scaleIn);

        /*         inline ISATNode *&node()
        {
            return node_;
        } */

        //- Increases the number of retrieves the chempoint has generated
        void increaseNumRetrieve();

        //- Resets the number of retrieves at each time step
        void resetNumRetrieve();

        inline const label &numRetrieve()
        {
            return numRetrieve_;
        }

        template <typename T>
        friend class parISATmanager;

        friend class parISATbinaryTree;
    };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#ifdef NoRepository
//#include "ISATleaf.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
