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
    Foam::functionTabulation

Description

SourceFiles
    functionTabulationI.H
    functionTabulation.C
    functionTabulationIO.C

\*---------------------------------------------------------------------------*/

#ifndef functionTabulation_H
#define functionTabulation_H

#include "parallelISAT_1D.H"
#include<fstream>
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    // Forward declaration of classes
    //class Istream;
    //class Ostream;

    // Forward declaration of friend functions and operators
    //class functionTabulation;
    //Istream& operator>>(Istream&, functionTabulation&);
    //Ostream& operator<<(Ostream&, const functionTabulation&);


    /*---------------------------------------------------------------------------*\
                             Class functionTabulation Declaration
    \*---------------------------------------------------------------------------*/

    class functionTabulation
    {
        // Private data

            //- Description of data_
            //dataType data_;


        // Private Member Functions

            //- Disallow default bitwise copy construct
            //functionTabulation(const functionTabulation&);

            //- Disallow default bitwise assignment
            //void operator=(const functionTabulation&);

        parallelISAT_1D tabulation_;


    public:

        // Static data members

            //- Static data staticData
        //static const dataType staticData;


        // Constructors

            //- Construct null
        functionTabulation();

        //- Construct from components
        //functionTabulation(const dataType& data);

        //- Construct from Istream
        //functionTabulation(Istream&);

        //- Construct as copy
        //functionTabulation(const functionTabulation&);


        // Selectors

            //- Select null constructed
        //static autoPtr<functionTabulation> New();


        //- Destructor
        ~functionTabulation();


        // Member Functions

            // Access
        scalar gradient(scalar x) { return cos(x); }
        scalar func(scalar x) { return sin(x); }
        scalar calculate(scalar x);
        void writeDot(string name);
        // Check

        // Edit

        // Write


    // Member Operators

    //void operator=(const functionTabulation&);


    // Friend Functions

    // Friend Operators

    // IOstream Operators

    //friend Istream& operator>>(Istream&, functionTabulation&);
    //friend Ostream& operator<<(Ostream&, const functionTabulation&);
    };


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#include "functionTabulationI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
