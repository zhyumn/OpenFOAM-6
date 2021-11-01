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
    Foam::ISAT

Description

SourceFiles
    ISATI.H
    ISAT.C
    ISATIO.C

\*---------------------------------------------------------------------------*/

#ifndef ISAT_H
#define ISAT_H

#include "SUPstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class ISAT;
Istream& operator>>(Istream&, ISAT&);
Ostream& operator<<(Ostream&, const ISAT&);


/*---------------------------------------------------------------------------*\
                         Class ISAT Declaration
\*---------------------------------------------------------------------------*/

class ISAT
{
    // Private data

        //- Description of data_
        //dataType data_;
        label nmem_leaf_;
        label nmem_node_;
        



    // Private Member Functions

        //- Disallow default bitwise copy construct
        ISAT(const ISAT&);

        //- Disallow default bitwise assignment
        void operator=(const ISAT&);
        struct Node;
        struct Leaf
        {
            //Node* node_;
            label v;

            //- Pointer to the chemistryModel object
            //TDACChemistryModel<CompType, ThermoType>& chemistry_;

            //- Vector storing the composition, temperature and pressure
            //  and deltaT if a variable time step is set on
            //scalarField phi_;

            //- Vector storing the mapping of the composition phi
            //scalarField Rphi_;

            //- LT the transpose of the L matrix describing the Ellipsoid Of
            //  Accuracy use List of Lists to be able to change size if DAC is used
            //scalarSquareMatrix LT_;

            //- A the mapping gradient matrix
            //scalarSquareMatrix A_;

            //- Vector storing the scale factor
            //scalarField scaleFactor_;

            //- Reference to the node in the binary tree holding this chemPoint
            Node* node_;

            //- The size of the composition space (size of the vector phi)
            //label completeSpaceSize_;

            //- Number of times the element has been grown
            //label nGrowth_;

            //- Tolerance for the Ellipsoid of accuracy
            //static scalar tolerance_;

            //- Number of active species stored in the chemPoint
            //label nActiveSpecies_;

            //- Vectors that store the index conversion between the simplified
            //  and the complete chemical mechanism
            //List<label> simplifiedToCompleteIndex_;

            //label timeTag_;
            //label lastTimeUsed_;

            //bool toRemove_;

            //label maxNumNewDim_;

            //Switch printProportion_;

            //- Variable to store the number of retrieves the chemPoint
            //  will generate at each time step
            //label numRetrieve_;

            //- Variable to store the number of time steps the chempoint is allowed
            //   to still live according to the maxChPLifeTime_ parameter
            //label nLifeTime_;

            //List<label> completeToSimplifiedIndex_;

            //- Number of equations in addition to the species eqs.
            //label nAdditionalEqns_;

            //label idT_;
            //label idp_;
            //label iddeltaT_;
            
        };
        struct Node
        {
            Leaf* leaf_;
            label v;

            Leaf* leafLeft_;

            //- Element on the right
            Leaf* leafRight_;

            //- Node which follows on the left
            Node* nodeLeft_;

            //- Node which follows on the right
            Node* nodeRight_;

            //- Parent node
            Node* parent_;

            //- Number of equations in addition to the species eqs.
            label nAdditionalEqns_;

            scalarField v_;

            scalar a_;
        };
        SUPstream::shared_data<Leaf> leafmem;
        SUPstream::shared_data<label> emptylist_leaf;
        SUPstream::shared_data<Node> nodemem;
        SUPstream::shared_data<label> emptylist_node;
        struct Private_variable
        {
            //         head      temp_tail       tail
            //         |         |               |
            //    | | |X|X|X|X|X|X|Y|Y|Y|Y|Y|Y|Y|Y| | | | 
            label size_leaf;
            label head_leaf;   
            label tail_leaf;
            label temp_tail_leaf;

            label size_node;
            label heaf_node;
            label tail_node;
            label temp_tail_node;
        };
        SUPstream::shared_data<Private_variable> pv_;
        label& size_leaf;
        label& head_leaf;
        label& tail_leaf;
        label& temp_tail_leaf;

        label& size_node;
        label& heaf_node;
        label& tail_node;
        label& temp_tail_node;

public:
        Node tnode;
        Leaf tleaf;
    // Static data members

        //- Static data staticData
        //static const dataType staticData;


    // Constructors

        //- Construct null
        //ISAT();

        //- Construct from components
        ISAT(SUPstream::mpi_manager& manager_in,label nmem);

        //- Construct from Istream
        //ISAT(Istream&);

        //- Construct as copy
        //ISAT(const ISAT&);


    // Selectors

        //- Select null constructed
        //static autoPtr<ISAT> New();


    //- Destructor
    ~ISAT();


    // Member Functions

        // Access

        // Check

        // Edit

        // Write


    // Member Operators

        //void operator=(const ISAT&);


    // Friend Functions

    // Friend Operators

    // IOstream Operators

        friend Istream& operator>>(Istream&, ISAT&);
        friend Ostream& operator<<(Ostream&, const ISAT&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#include "ISATI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
