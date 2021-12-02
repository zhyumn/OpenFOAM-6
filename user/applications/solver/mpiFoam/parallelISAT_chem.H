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
#include "SVD.H"
#include "chemistryTabulationMethod.H"
#include <iostream>
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    template<class CompType, class ThermoType>
    class parallelISAT_chem;
    /*---------------------------------------------------------------------------*\
                             Class parallelISAT_chem Declaration
    \*---------------------------------------------------------------------------*/
    template<class CompType, class ThermoType>
    struct ISAT_chem
    {
        typedef FixedList<scalar, 3> inputType;
        typedef FixedList<scalar, 3> outputType;
        typedef FixedList<scalar, 9> gradientType;
        static parallelISAT_chem<CompType, ThermoType>* pISAT;
        class nodeData;
        struct leafData
        {
        private:
            inputType phi_;
            outputType Rphi_;
            gradientType A_;//gradient
            FixedList<scalar, 9> LT_;
            FixedList<scalar, 3> scaleFactor_;
            scalar EOA;
        public:
            //void set(const inputType& x) { v = x; Rv = x * x; gradient = 2 * x; EOA = pISAT->tolerance_ / max(fabs(gradient), 1); }
            void set(const inputType& x) {};
            template<typename ...Args>
            void set(const inputType& x, const outputType& y, Args ...args) {};
            static void computeA(gradientType& A, const inputType& Rphiq, const scalar rhoi, const scalar dt) {}
            bool grow(const inputType& phiq);
            static outputType func(const inputType& x);
            static void gradFunc(const inputType& x, gradientType& dxdy);
            //bool collide(const inputType& x) { return x == v; }
            void retrieve(const inputType& x, outputType& y);
            bool inEOA(const inputType& x);
            friend nodeData;
            template<class T1, class T2>
            friend std::ostream& operator<<(std::ostream& out, typename ISAT_chem<T1, T2>::leafData& A);
        };
        struct nodeData
        {
        private:
            scalarField v_;
            scalar a_;
        public:
            void set() { v_[0] = 0;a_ = 0; }

            void calcV
            (
                const leafData& elementLeft,
                const leafData& elementRight,
                scalarField& v
            );

            scalar calcA
            (
                const leafData& elementLeft,
                const leafData& elementRight
            );

            void set(const inputType& x, const inputType& y) { v_[0] = (x[0] + y[0]) / 2; }
            void set(const leafData& x, const leafData& y) { v_[0] = (x.v_[0] + y.v_[0]) / 2; }
            bool goLeft(inputType& x) { return x[0] < a_; }
            bool goLeft(leafData& x) { return x.a_ < a_; }
            friend leafData;
            template<class T1, class T2>
            friend std::ostream& operator<<(std::ostream& out, typename ISAT_chem<T1, T2>::nodeData& A);
        };
    };
    template<class CompType, class ThermoType>
    std::ostream& operator<<(std::ostream& out, typename ISAT_chem<CompType, ThermoType>::leafData& A);
    template<class CompType, class ThermoType>
    std::ostream& operator<<(std::ostream& out, typename ISAT_chem<CompType, ThermoType>::nodeData& A);

    template<class CompType, class ThermoType>
    class parallelISAT_chem
        :
        public parallelISAT<ISAT_chem<CompType, ThermoType>, chemistryTabulationMethod<CompType, ThermoType>>
    {
        // Private data

            //- Description of data_
        //dataType data_;


        // Private Member Functions

            //- Disallow default bitwise copy construct


        //- Disallow default bitwise assignment



    public:

        // Static data members

            //- Static data staticData
        //static const dataType staticData;

        typedef ISAT_chem<CompType, ThermoType> DataType;
        // Constructors

            //- Construct null


        //- Construct from components
        parallelISAT_chem(SUPstream::mpi_manager& manager_in, label nmem, SUPstream::mpi_sync& sync_in, const dictionary& dict,
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
        label size() { return 0; }
        void writePerformance() { return; }

        // Check

        // Edit
        bool update() { return 0; }

        label add
        (
            const scalarField& phiq,
            const scalarField& Rphiq,
            const scalar rho,
            const scalar deltaT
        )
        {
            return parallelISAT<ISAT_chem<CompType, ThermoType>, chemistryTabulationMethod<CompType, ThermoType>>::add(phiq, Rphiq, rho, deltaT);
        }

        bool retrieve
        (
            const Foam::scalarField& phiq,
            scalarField& Rphiq
        );
        //void computeA(scalarSquareMatrix& A, const calarField& Rphiq, const scalar rhoi, const scalar dt);
        // Write



        // Member Operators

        //void operator=(const parallelISAT_chem&);


        // Friend Functions

        // Friend Operators

        // IOstream Operators

        //friend Istream& operator>>(Istream&, parallelISAT_chem&);
        //friend Ostream& operator<<(Ostream&, const parallelISAT_chem&);
    };


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "parallelISAT_chem.C"
#endif
//#include "parallelISAT_chemI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
