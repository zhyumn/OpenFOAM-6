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
    Foam::ISATNode

Description
    Node of the binary tree

SourceFile
    ISATNode.C

\*---------------------------------------------------------------------------*/

#ifndef parData_H
#define parData_H

#include "slab.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace chemistryTabulationMethodSs
    {

        template <class Type>
        struct SList
        {
            int size_;
            SharedPointer<Type> v_;
            void init(int size_in)
            {
                size_ = size_in;
                v_ = pslab->alloc(sizeof(Type) * size_);
                if (v_.offset == sptr_NULL)
                {
                    FatalErrorInFunction << "Run out of memory\n"
                                         << exit(FatalError);
                }
            }
            void init(const SList<Type> &x)
            {
                init(x.size_);
                this->operator=(x);
            }
            inline Type &operator[](int x)
            {
                /*                 if (x >= size_)
                    FatalErrorInFunction
                        << abort(FatalError); */
                //assert(x < size_);
                return v_()[x];
            }
            inline const Type &operator[](int x) const
            {
                /*                 if (x >= size_)
                    FatalErrorInFunction
                        << "x=" << x
                        << ",size=" << size_
                        << abort(FatalError); */
                //assert(x < size_);
                return v_()[x];
            }

            inline const SList<Type> &operator=(const SList<Type> &x)
            {
                for (int i = 0; i < size_; i++)
                    v_()[i] = x[i];
                return *this;
            }

            inline const SList<Type> &operator=(const List<Type> &x)
            {
                for (int i = 0; i < size_; i++)
                    v_()[i] = x[i];
                return *this;
            }
            template <int Datasize>
            inline const SList<Type> &operator=(const FixedList<scalar, Datasize> &x)
            {
                if (Datasize >= size_)
                {
                    for (int i = 0; i < size_; i++)
                        v_()[i] = x[i];
                }
                return *this;
            }
            inline void setSize(int size_in)
            {
                SharedPointer<Type> tmp;
                /*if (UNLIKELY(size_in > size_))
                {
                    tmp = pslab->alloc(sizeof(Type) * size_in);
                    size_t ms = min(size_in, size_);
                    memcpy(tmp(), v_(), ms * sizeof(Type));
                    //for (size_t i = 0; i < ms; i++)
                    //    tmp()[i] = v_()[i];
                    Delete();

                    size_ = size_in;
                    v_ = tmp;
                }*/
                //init(size_in);
            }
            const int &size() const
            {
                return size_;
            }
            void Delete()
            {
                if (size_ != 0)
                    pslab->free(sizeof(Type) * size_, v_.offset);
            }
            ~SList()
            {
                pslab->free(sizeof(Type) * size_, v_.offset);
            }
        };

        template <class Type>
        struct SSquareMatrix
        {
            int size_;
            int total_size_;
            SharedPointer<Type> v_;

            void init()
            {
                size_ = 0;
            }
            void init(int size_in)
            {
                size_ = size_in;
                total_size_ = size_in * size_in;
                v_ = pslab->alloc(sizeof(Type) * total_size_);
                if (v_.offset == sptr_NULL)
                {
                    FatalErrorInFunction << "Run out of memory\n"
                                         << exit(FatalError);
                }
            }
            void init(const SSquareMatrix<Type> &x)
            {
                init(x.size_);
                this->operator=(x);
            }
            inline Type &operator()(int x, int y)
            {
                //assert(x < size_);
                //assert(y < size_);
                return v_()[x * size_ + y];
            }
            inline const Type &operator()(int x, int y) const
            {
                //assert(x < size_);
                //assert(y < size_);
                return v_()[x * size_ + y];
            }

            inline const SSquareMatrix<Type> &operator=(const SSquareMatrix<Type> &x)
            {
                for (int i = 0; i < total_size_; i++)
                    v_()[i] = x.v_()[i];
                return *this;
            }
            inline void setSize(int size_in)
            {
                Delete();
                init(size_in);
            }
            const int &size() const
            {
                return total_size_;
            }
            void Delete()
            {
                if (size_ != 0)
                    pslab->free(sizeof(Type) * total_size_, v_.offset);
            }
            ~SSquareMatrix()
            {
                pslab->free(sizeof(Type) * total_size_, v_.offset);
            }
        };

        template <class Type>
        struct FSSquareMatrix
        {
            Type *v_;
            int size_;
            int total_size_;

            SSquareMatrix<Type> *p;
            FSSquareMatrix(SSquareMatrix<Type> &a) : v_(a.v_()), size_(a.size_), total_size_(a.total_size_), p(&a) {}
            inline Type &operator()(int x, int y)
            {
                //assert(x < size_);
                //assert(y < size_);
                return v_[x * size_ + y];
            }
            inline const Type &operator()(int x, int y) const
            {
                //assert(x < size_);
                //assert(y < size_);
                return v_[x * size_ + y];
            }
            inline void setSize(int size_in)
            {
                p->setSize(size_in);
                size_ = p->size_;
                total_size_ = p->total_size_;
                v_ = p->v_();
            }
            const int &size() const
            {
                return total_size_;
            }
        };
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#ifdef NoRepository
//    #include "ISATNode.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
