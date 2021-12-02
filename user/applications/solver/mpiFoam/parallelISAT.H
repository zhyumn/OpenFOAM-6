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
    Foam::parallelISAT

Description

SourceFiles
    ISATI.H
    ISAT.C
    ISATIO.C

\*---------------------------------------------------------------------------*/

#ifndef PARALLELISAT_H
#define PARALLELISAT_H

#include "SUPstream.H"
#include<fstream>
#include "IOstreams.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    // Forward declaration of classes
    class Istream;
    class Ostream;

    // Forward declaration of friend functions and operators
    template<typename T, typename T2>
    class parallelISAT;
    //Istream& operator>>(Istream&, parallelISAT&);
    //Ostream& operator<<(Ostream&, const parallelISAT&);


    /*---------------------------------------------------------------------------*\
                             Class parallelISAT Declaration
    \*---------------------------------------------------------------------------*/
    template<typename T, typename T2>
    class parallelISAT :public T2
    {

    public:
        typedef typename T::inputType inputType;
        typedef typename T::outputType outputType;
    private:
        // Private data

            //- Description of data_
            //dataType data_;
        label nmem_leaf_;
        label nmem_node_;

        scalar tolerance_;




        // Private Member Functions

            //- Disallow default bitwise copy construct
        parallelISAT(const parallelISAT&);

        //- Disallow default bitwise assignment
        void operator=(const parallelISAT&);
        struct Node;
        struct Leaf :public T::leafData
        {
            //Node* node_;
            //label v;
            //scalar v;

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
            // Only for write
            label node_;

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
        struct Node :public T::nodeData
        {
            //Leaf* leaf_;
            //label v;
            //scalar v;

            label leafLeft_;

            //- Element on the right
            label leafRight_;

            //- Node which follows on the left

            //first check node
            label nodeLeft_;

            //- Node which follows on the right
            label nodeRight_;

            //- Parent node  Only for write
            label parent_;

            //- Number of equations in addition to the species eqs.
            //label nAdditionalEqns_;

            //scalarField v_;

            //scalar a_;
        };
        SUPstream::shared_data<Leaf> leafmem;
        SUPstream::shared_data<label> emptylist_leaf;
        SUPstream::shared_data<Node> nodemem;
        SUPstream::shared_data<label> emptylist_node;

        class Mem_manager
        {
            SUPstream::shared_data<label>& emptylist_;
            //label& size_;

            label& head_;
            label& temp_tail_;
            label& tail_;
            label& nmem_;
        public:
            Mem_manager(
                SUPstream::shared_data<label>& emptylist,
                label& head, label& tail, label& temp_tail, label& nmem
            ) :emptylist_(emptylist), head_(head),
                temp_tail_(temp_tail), tail_(tail), nmem_(nmem)
            {
                head_ = 0;
                temp_tail_ = nmem - 1;
                tail_ = nmem - 1;
                for (label i = 0;i < nmem;i++)
                {
                    emptylist_[i] = i;
                }
            }
            inline label New() { if (head_ == temp_tail_)return -1; label ret = head_;head_ = (head_ + 1) % nmem_;return ret; }
            inline void Delete(label index) { tail_ = (tail_ + 1) % nmem_;emptylist_[tail_] = index; }
            inline void update() { temp_tail_ = tail_; }
        };
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

            label root;
        };
        SUPstream::shared_data<Private_variable> pv_;
        label& size_leaf;
        label& head_leaf;
        label& tail_leaf;
        label& temp_tail_leaf;

        label& size_node;
        label& head_node;
        label& tail_node;
        label& temp_tail_node;

        label& root_;
        SUPstream::mpi_mutex write_lock;

        Mem_manager leaf_manager, node_manager;

    public:
        //Node tnode;
        //Leaf tleaf;
    // Static data members

        //- Static data staticData
        //static const dataType staticData;


    // Constructors

        //- Construct null
        //parallelISAT();

        //- Construct from components
        template<typename ...Args>
        parallelISAT(SUPstream::mpi_manager& manager_in, label nmem, SUPstream::mpi_sync& sync_in, Args ...args);

        //- Construct from Istream
        //ISAT(Istream&);

        //- Construct as copy
        //ISAT(const ISAT&);


    // Selectors

        //- Select null constructed
        //static autoPtr<ISAT> New();


    //- Destructor
        ~parallelISAT();


        // Member Functions

        // Access
        scalar tolerance() { return tolerance_; }
        label binaryTreeSearch(inputType value, label node_index);
        label binaryTreeSearch(inputType value);
        void writeDot(string name);
        void writeDot(label nodeindex, std::ofstream& fout);
        scalar leafvalue(label index) { if (index == -1)return -1;return leafmem[index].v; };

        bool retrieve(const inputType& x, outputType& y);
        // Check

        // Edit
        bool grow(const inputType& value);
        //void insert(const inputType& value);
        template<typename ...Args>
        label add(const inputType& x, const outputType& y, Args ...args);
        label add(const inputType& x);
        /* virtual label add
         (
             const scalarField& phiq,
             const scalarField& Rphiq,
             const scalar rho,
             const scalar deltaT
         );
         */
         // Write


     // Member Operators

         //void operator=(const ISAT&);


     // Friend Functions

     // Friend Operators

     // IOstream Operators

         //friend Istream& operator>>(Istream&, parallelISAT&);
         //friend Ostream& operator<<(Ostream&, const parallelISAT&);
    };


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "parallelISATI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
