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
#include "FixedSquareMatrix.H"
#include "slab.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#define Datasize 56
namespace Foam
{
    namespace chemistryTabulationMethods
    {

        template<class CompType, class ThermoType>
        class parallelISAT_chem;
        /*---------------------------------------------------------------------------*\
                                 Class parallelISAT_chem Declaration
        \*---------------------------------------------------------------------------*/
        template<class CompType, class ThermoType>
        struct ISAT_chem
        {
            typedef FixedList<scalar, Datasize> inputType;
            typedef FixedList<scalar, Datasize> outputType;
            typedef FixedSquareMatrix<scalar, Datasize> gradientType;

            static parallelISAT_chem<CompType, ThermoType>* pISAT;

            //static label out;
            class nodeData;
            struct leafData
            {
            private:
                inputType phi_;
                outputType Rphi_;
                gradientType A_;//gradient
                gradientType LT_;


                inputType simplifiedToCompleteIndex_;
                inputType completeToSimplifiedIndex_;

                label nActiveSpecies_;
                outputType scaleFactor_;
                label nGrowth_;
                scalar timeTag_;
                scalar EOA;
            public:
                //void set(const inputType& x) { v = x; Rv = x * x; gradient = 2 * x; EOA = pISAT->tolerance_ / max(fabs(gradient), 1); }
                //void set(const inputType& x);

                void set(const inputType& x, const outputType& y, const scalar rhoi, const scalar dt);
                //template<typename ...Args>
                //static void computeA(const inputType& x, const outputType& y, gradientType& A, Args ...args);

                static void computeA
                (
                    const inputType& phiq,
                    const outputType& Rphiq,
                    gradientType& A_out,
                    scalarSquareMatrix& A,
                    const scalar rhoi,
                    const scalar dt
                );

                void qrDecompose
                (
                    const label nCols,
                    scalarSquareMatrix& R
                );
                template<class M>
                void qrUpdate
                (
                    M& R,
                    const label n,
                    const Foam::scalarField& u,
                    const Foam::scalarField& v
                );
                template<class M>
                void rotate
                (
                    M& R,
                    const label i,
                    const scalar a,
                    const scalar b,
                    label n
                );

                bool checkSolution
                (
                    const inputType& phiq,
                    const outputType& Rphiq
                );

                bool grow(const inputType& x);
                static outputType func(const inputType& x);
                static gradientType gradFunc(const inputType& x);
                //bool collide(const inputType& x) { return x == phi_; }
                void retrieve(const inputType& x, outputType& y);
                bool inEOA(const inputType& x);
                void computeLT(scalarSquareMatrix& A);
                label simplifiedToCompleteIndex(const label i);
                inline const scalar& timeTag() { return timeTag_; }
                inline const label& nGrowth() { return nGrowth_; }
                inline const inputType& phi() { return phi_; };
                friend nodeData;
                template<class CompType_, class ThermoType_>
                friend std::ostream& operator<<(std::ostream& out, typename ISAT_chem<CompType_, ThermoType_>::leafData& A);
            };
            struct nodeData
            {
            private:
                inputType v_;
                scalar a_;
            public:
                void set() { }
                //void set(const inputType& x, const inputType& y) {  }
                void set(leafData& x, leafData& y)
                {
                    calcV(x, y, v_);
                    a_ = calcA(x, y);
                }
                void calcV
                (
                    leafData& elementLeft,
                    leafData& elementRight,
                    inputType& v
                );
                scalar calcA
                (
                    const leafData& elementLeft,
                    const leafData& elementRight
                );

                bool goLeft(inputType& x) {
                    scalar vPhi = 0.0;
                    for (label i = 0; i < Datasize; i++)
                        vPhi += x[i] * v_[i];
                    return vPhi < a_;
                }
                bool goLeft(leafData& x)
                {
                    scalar vPhi = 0.0;
                    for (label i = 0; i < Datasize; i++)
                        vPhi += x.phi_[i] * v_[i];
                    return vPhi < a_;
                    //return true;
                }
                friend leafData;
                template<class CompType_, class ThermoType_>
                friend std::ostream& operator<<(std::ostream& out, typename ISAT_chem<CompType_, ThermoType_>::nodeData& A);
            };
        };
        template<class CompType, class ThermoType>
        std::ostream& operator<<(std::ostream& out, typename ISAT_chem<CompType, ThermoType>::leafData& A);
        template<class CompType, class ThermoType>
        std::ostream& operator<<(std::ostream& out, typename ISAT_chem<CompType, ThermoType>::nodeData& A);
        struct emptyClass
        {};
        template<class CompType, class ThermoType>
        class parallelISAT_chem
            :public chemistryTabulationMethod<CompType, ThermoType>,
            public parallelISAT<ISAT_chem<CompType, ThermoType>, emptyClass>
        {
            // Private data

                //- Description of data_
            //dataType data_;

            //- Number of equations in addition to the species eqs.
            label nAdditionalEqns_;

            const Time& runTime_;

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

            //label memsize_;
            //label maxmem_;
            //MPI_Shared::Slab slab;

            autoPtr<OFstream> nRetrievedFile_;
            autoPtr<OFstream> nGrowthFile_;
            autoPtr<OFstream> nAddFile_;
            autoPtr<OFstream> sizeFile_;
            autoPtr<OFstream> totalGrowthFile_;

            // Private Member Functions

                //- Disallow default bitwise copy construct


            //- Disallow default bitwise assignment



        public:

            TypeName("parallelISAT");

            // Static data members

            //- Static data staticData
        //static const dataType staticData;

            typedef ISAT_chem<CompType, ThermoType> DataType;
            // Constructors

                //- Construct null


            //- Construct from components
            parallelISAT_chem
            (
                const dictionary& dict,
                TDACChemistryModel<CompType, ThermoType>& chemistry
            );

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
                    parallelISAT<ISAT_chem<CompType, ThermoType>, emptyClass>::update();
                    //if (!this->treeCheck())
                    //{
                    //    FatalErrorInFunction
                    //        << "broken tree 1\n"
                    //        << "out = " << DataType::out
                    //        << abort(FatalError);
                    //}

                    ret = cleanAndBalance();
                    if (this->isFull())
                    {
                        this->clear();
                        ret = true;
                    }
                    parallelISAT<ISAT_chem<CompType, ThermoType>, emptyClass>::update();
                    //if (!this->treeCheck())
                    //{
                    //    FatalErrorInFunction
                    //        << "broken tree 2"
                    //        << abort(FatalError);
                    //}
                }
                return ret;
            }
            label size() { return this->size_leaf(); };
            void writePerformance() {

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

                //DataType::out++;

                /*if (DataType::out == 1)
                {

                    FatalErrorInFunction
                        << "out=" << DataType::out
                        << abort(FatalError);
                }*/
            }
            bool retrieve
            (
                const scalarField& phiQ,
                scalarField& RphiQ
            )
            {
                /*if (!this->treeCheck())
                {
                    FatalErrorInFunction
                        << "broken tree retr0\n"
                        << "out = " << DataType::out
                        << abort(FatalError);
                }*/
                //static label out0 = 0, out1 = 0;
                /*if (out0 == 3000)
                {
                    FatalErrorInFunction
                        << "out0=" << out0
                        << "\nout1=" << out1
                        << abort(FatalError);
                }*/
                //out0++;
                typename DataType::inputType phiQ_fixed, RphiQ_fixed;
                for (label i = 0;i < Datasize;i++)
                {
                    phiQ_fixed[i] = phiQ[i];
                }
                bool flag = parallelISAT<ISAT_chem<CompType, ThermoType>, emptyClass>::retrieve(phiQ_fixed, RphiQ_fixed);
                for (label i = 0;i < Datasize;i++)
                {
                    RphiQ[i] = RphiQ_fixed[i];
                }
                //if (flag == true)
                //{
                //    out1++;
                //}
                /*if (!this->treeCheck())
                {
                    FatalErrorInFunction
                        << "broken tree retr1\n"
                        << "out = " << DataType::out
                        << abort(FatalError);
                }*/
                return flag;
            };

            label add
            (
                const scalarField& phiQ,
                const scalarField& RphiQ,
                const scalar rho,
                const scalar deltaT
            )
            {

                /*if (!this->treeCheck())
                {
                    FatalErrorInFunction
                        << "broken tree add0\n"
                        << "out = " << DataType::out
                        << abort(FatalError);
                }*/

                typename DataType::inputType phiQ_fixed, RphiQ_fixed;
                for (label i = 0;i < Datasize;i++)
                {
                    phiQ_fixed[i] = phiQ[i];
                    RphiQ_fixed[i] = RphiQ[i];
                }

                label ret = parallelISAT<ISAT_chem<CompType, ThermoType>, emptyClass>::add(phiQ_fixed, RphiQ_fixed, rho, deltaT);
                /*if (!this->treeCheck())
                {
                    FatalErrorInFunction
                        << "broken tree add1\n"
                        << "out = " << DataType::out
                        << abort(FatalError);
                }*/
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
    } // End namespace chemistryTabulationMethods
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "parallelISAT_chem.C"
#endif
//#include "parallelISAT_chemI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
