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

Class
    Foam::parallelISAT_chem

Description

SourceFiles
    parallelISAT_chemI.H
    parallelISAT_chem.C
    parallelISAT_chemIO.C

\*---------------------------------------------------------------------------*/

#ifndef parallelISAT_chem_H
#define parallelISAT_chem_H

#include "parallelISAT.H"
#include <iostream>
//#include "FixedSquareMatrix.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#define Datasize 56

#define DateSize1 53
#define DateSize2 56
#define LIKELY(exp) __builtin_expect(exp, 1)
#define UNLIKELY(exp) __builtin_expect(exp, 0)
namespace Foam
{
    namespace chemistryTabulationMethodSs
    {

        template <class CompType, class ThermoType>
        class parallelISAT_chem;

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

            inline const SList<Type> &operator=(const FixedList<scalar, Datasize> &x)
            {
                for (int i = 0; i < size_; i++)
                    v_()[i] = x[i];
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

        /*---------------------------------------------------------------------------*\
                                 Class parallelISAT_chem Declaration
        \*---------------------------------------------------------------------------*/
        template <class CompType, class ThermoType>
        struct ISAT_chem
        {
            typedef SList<scalar> inputType;
            typedef SList<label> inputType_label;
            typedef SList<scalar> outputType;
            typedef SSquareMatrix<scalar> gradientType;

            static parallelISAT_chem<CompType, ThermoType> *pISAT;

            //static label out;
            class nodeData;
            struct leafData
            {
            private:
                inputType phi_;
                outputType Rphi_;
                gradientType A_; //gradient
                gradientType LT_;

                inputType_label simplifiedToCompleteIndex_;
                inputType_label completeToSimplifiedIndex_;

                label nActiveSpecies_;
                outputType scaleFactor_;
                label nGrowth_;
                scalar timeTag_;
                scalar EOA;

            public:
                void set(const scalarField &x, const scalarField &y, const scalar rhoi, const scalar dt);

                void init(const leafData &x)
                {
                    phi_.init(x.phi_);
                    Rphi_.init(x.Rphi_);
                    A_.init(x.A_);
                    LT_.init(x.LT_);

                    simplifiedToCompleteIndex_.init(x.simplifiedToCompleteIndex_);
                    completeToSimplifiedIndex_.init(x.completeToSimplifiedIndex_);

                    nActiveSpecies_ = x.nActiveSpecies_;
                    scaleFactor_.init(x.scaleFactor_);
                    nGrowth_ = x.nGrowth_;
                    timeTag_ = x.timeTag_;
                    EOA = x.EOA;
                }

                static void computeA(
                    const scalarField &phiq,
                    const scalarField &Rphiq,
                    gradientType &A_out,
                    scalarSquareMatrix &A,
                    const scalar rhoi,
                    const scalar dt);

                void qrDecompose(
                    const label nCols,
                    scalarSquareMatrix &R);
                template <class M>
                void qrUpdate(
                    M &R,
                    const label n,
                    const Foam::scalarField &u,
                    const Foam::scalarField &v);
                template <class M>
                void rotate(
                    M &R,
                    const label i,
                    const scalar a,
                    const scalar b,
                    label n);

                bool checkSolution(
                    const scalarField &phiq,
                    const scalarField &Rphiq);

                bool grow(const scalarField &x);
                static outputType func(const inputType &x);
                static gradientType gradFunc(const inputType &x);
                void retrieve(const scalarField &x, scalarField &y);
                void Delete()
                {
                    phi_.Delete();
                    Rphi_.Delete();
                    A_.Delete();
                    LT_.Delete();
                    simplifiedToCompleteIndex_.Delete();
                    completeToSimplifiedIndex_.Delete();
                    scaleFactor_.Delete();
                }
                void init(int size_in)
                {
                    phi_.init(size_in);
                    Rphi_.init(size_in);
                    A_.init();
                    LT_.init();
                    simplifiedToCompleteIndex_.init(size_in);
                    completeToSimplifiedIndex_.init(size_in);
                    scaleFactor_.init(size_in);
                }
                bool inEOA(const scalarField &x);

                void computeLT(scalarSquareMatrix &A);
                label simplifiedToCompleteIndex(const label i);
                inline const scalar &timeTag() { return timeTag_; }
                inline const label &nGrowth() { return nGrowth_; }
                inline const inputType &phi() { return phi_; };
                friend nodeData;
                template <class CompType_, class ThermoType_>
                friend std::ostream &operator<<(std::ostream &out, typename ISAT_chem<CompType_, ThermoType_>::leafData &A);
            };
            struct nodeData
            {
            private:
                inputType v_;
                scalar a_;

            public:
                void set() {}

                void set(leafData &x, leafData &y)
                {
                    calcV(x, y, v_);
                    a_ = calcA(x, y);
                }
                void calcV(
                    leafData &elementLeft,
                    leafData &elementRight,
                    inputType &v);
                scalar calcA(
                    const leafData &elementLeft,
                    const leafData &elementRight);

                bool goLeft(const inputType &x)
                {
                    scalar vPhi = 0.0;
                    for (label i = 0; i < Datasize; i++)
                        vPhi += x[i] * v_[i];
                    return vPhi < a_;
                }

                bool goLeft(const scalarField &x)
                {
                    scalar vPhi = 0.0;
                    for (label i = 0; i < x.size(); i++)
                        vPhi += x[i] * v_[i];
                    return vPhi < a_;
                }
                bool goLeft(leafData &x)
                {
                    scalar vPhi = 0.0;
                    for (label i = 0; i < Datasize; i++)
                        vPhi += x.phi_[i] * v_[i];
                    return vPhi < a_;
                    //return true;
                }
                void Delete()
                {
                    v_.Delete();
                }
                void init(int size_in)
                {
                    v_.init(size_in);
                }
                friend leafData;
                template <class CompType_, class ThermoType_>
                friend std::ostream &operator<<(std::ostream &out, typename ISAT_chem<CompType_, ThermoType_>::nodeData &A);
            };
        };
        template <class CompType, class ThermoType>
        std::ostream &operator<<(std::ostream &out, typename ISAT_chem<CompType, ThermoType>::leafData &A);
        template <class CompType, class ThermoType>
        std::ostream &operator<<(std::ostream &out, typename ISAT_chem<CompType, ThermoType>::nodeData &A);
        struct emptyClass
        {
        };
        template <class CompType, class ThermoType>
        class parallelISAT_chem
            : public chemistryTabulationMethodS<CompType, ThermoType>,
              public parallelISAT<ISAT_chem<CompType, ThermoType>, emptyClass>
        {
            // Private data

            //- Number of equations in addition to the species eqs.
            label nAdditionalEqns_;

            const Time &runTime_;

            //scalar tolerance_;

            label completeSpaceSize_;

            label iddeltaT_;

            scalarField scaleFactor_;

            bool printProportion_;

            label maxNumNewDim_;

            label idT_;
            label idp_;

            label chPMaxLifeTime_;
            label maxGrowth_;

            scalar minBalanceThreshold_;

            scalar maxDepthFactor_;

            autoPtr<OFstream> nRetrievedFile_;
            autoPtr<OFstream> nGrowthFile_;
            autoPtr<OFstream> nAddFile_;
            autoPtr<OFstream> sizeFile_;
            autoPtr<OFstream> totalGrowthFile_;

            //SUPstream::mpi_mutex test_lock;

            // Private Member Functions

            //- Disallow default bitwise copy construct

            //- Disallow default bitwise assignment

        public:
            TypeName("parallelISAT");

            // Static data members

            //- Static data staticData

            typedef ISAT_chem<CompType, ThermoType> DataType;
            // Constructors

            //- Construct null

            //- Construct from components
            parallelISAT_chem(
                const dictionary &dict,
                STDACChemistryModel<CompType, ThermoType> &chemistry);

            //- Construct from Istream
            //parallelISAT_chem(Istream&);

            //- Construct as copy
            //parallelISAT_chem(const parallelISAT_chem&);

            // Selectors

            //- Select null constructed
            //static autoPtr<parallelISAT_chem> New();

            //- Destructor
            ~parallelISAT_chem();

            // Member Functions

            // Access

            // Check

            // Edit
            //void computeA(scalarSquareMatrix& A, const calarField& Rphiq, const scalar rhoi, const scalar dt);
            // Write
            void balance();

            bool cleanAndBalance();

            bool update()
            {
                bool ret = false;
                this->sync_.sync();
                if (this->manager_.rank == 0)
                {
                    //parallelISAT<ISAT_chem<CompType, ThermoType>, emptyClass>::update();
                    if (this->chemistry_.timeSteps() % 20 == 0)
                        ret = cleanAndBalance();
                    if (this->isFull())
                    {
                        this->clear();
                        ret = true;
                    }
                    parallelISAT<ISAT_chem<CompType, ThermoType>, emptyClass>::update();
                }
                this->sync_.sync();
                return ret;
            }
            label size() { return this->size_leaf(); };
            void writePerformance()
            {

                if (this->log())
                {
                    nRetrievedFile_()
                        << runTime_.timeOutputValue() << "    " << this->nRetrieved_ << endl;
                    this->nRetrieved_ = 0;

                    nGrowthFile_()
                        << runTime_.timeOutputValue() << "    " << this->nGrowth_ << endl;
                    this->nGrowth_ = 0;

                    nAddFile_()
                        << runTime_.timeOutputValue() << "    " << this->nAdd_ << endl;
                    this->nAdd_ = 0;

                    sizeFile_()
                        << runTime_.timeOutputValue() << "    " << size() << endl;

                    totalGrowthFile_()
                        << runTime_.timeOutputValue() << "    " << this->totalGrowth_ << endl;
                    this->sync_.sync();
                    this->totalGrowth_ = 0;
                }
            }
            bool retrieve(
                const scalarField &phiQ,
                scalarField &RphiQ)
            {
                //bool flag = parallelISAT<ISAT_chem<CompType, ThermoType>, emptyClass>::retrieve(phiQ, RphiQ);
                return parallelISAT<ISAT_chem<CompType, ThermoType>, emptyClass>::retrieve(phiQ, RphiQ);
            };

            bool retrieve(
                const FixedList<scalar, DateSize2> &phiQ,
                FixedList<scalar, DateSize2> &RphiQ)
            {
                scalarField phiQ_2(phiQ.size()), RphiQ_2(phiQ.size());
                for (label i = 0; i < phiQ.size(); i++)
                {
                    phiQ_2[i] = phiQ[i];
                    RphiQ_2[i] = RphiQ[i];
                }
                bool ret = parallelISAT<ISAT_chem<CompType, ThermoType>, emptyClass>::retrieve(phiQ_2, RphiQ_2);

                for (label i = 0; i < phiQ.size(); i++)
                {
                    phiQ[i] = phiQ_2[i];
                    RphiQ[i] = RphiQ_2[i];
                }
                //bool flag = parallelISAT<ISAT_chem<CompType, ThermoType>, emptyClass>::retrieve(phiQ, RphiQ);
                return ret;
            };

            label add(
                const scalarField &phiQ,
                const scalarField &RphiQ,
                const scalar rho,
                const scalar deltaT)
            {
                //label ret = parallelISAT<ISAT_chem<CompType, ThermoType>, emptyClass>::add(phiQ, RphiQ, rho, deltaT);
                return parallelISAT<ISAT_chem<CompType, ThermoType>, emptyClass>::add(phiQ, RphiQ, rho, deltaT);
            };

            label add(
                const FixedList<scalar, DateSize2> &phiQ,
                const FixedList<scalar, DateSize2> &RphiQ,
                const scalar rho,
                const scalar deltaT)
            {
                //label ret = parallelISAT<ISAT_chem<CompType, ThermoType>, emptyClass>::add
                scalarField phiQ_2(phiQ.size()), RphiQ_2(phiQ.size());
                for (label i = 0; i < phiQ.size(); i++)
                {
                    phiQ_2[i] = phiQ[i];
                    RphiQ_2[i] = RphiQ[i];
                }
                bool ret = parallelISAT<ISAT_chem<CompType, ThermoType>, emptyClass>::add(phiQ_2, RphiQ_2, rho, deltaT);
                for (label i = 0; i < phiQ.size(); i++)
                {
                    phiQ[i] = phiQ_2[i];
                    RphiQ[i] = RphiQ_2[i];
                }

                return ret;
            };

            // Member Operators

            //void operator=(const parallelISAT_chem&);

            // Friend Functions

            // Friend Operators

            // IOstream Operators

            //friend Istream& operator>>(Istream&, parallelISAT_chem&);
            //friend Ostream& operator<<(Ostream&, const parallelISAT_chem&);
            friend ISAT_chem<CompType, ThermoType>;
        };

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    } // End namespace chemistryTabulationMethodSs
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#ifdef NoRepository
//#include "parallelISAT_chem.C"
//#endif

#include "parallelISAT_chemI.H"
//#include "parallelISAT_chemI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
