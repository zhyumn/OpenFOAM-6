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
    Foam::chemistryReductionMethodS

Description
    An abstract class for methods of chemical mechanism reduction

SourceFiles
    chemistryReductionMethodS.C

\*---------------------------------------------------------------------------*/

#ifndef chemistryReductionMethodS_H
#define chemistryReductionMethodS_H

#include "IOdictionary.H"
#include "Switch.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class CompType, class ThermoType>
class STDACChemistryModel;

/*---------------------------------------------------------------------------*\
                           Class chemistrySolver Declaration
\*---------------------------------------------------------------------------*/

template<class CompType, class ThermoType>
class chemistryReductionMethodS
{

protected:

    const IOdictionary& dict_;

    //- Dictionary that store the algorithm data
    const dictionary coeffsDict_;

    //- Is mechanism reduction active?
    Switch active_;

    //- Switch to select performance logging
    Switch log_;

    STDACChemistryModel<CompType, ThermoType>& chemistry_;

    //- List of active species (active = true)
    List<bool> activeSpecies_;

    //- Number of active species
    label NsSimp_;

    //- Number of species
    const label nSpecie_;

    //- Tolerance for the mechanism reduction algorithm
    scalar tolerance_;


public:

    //- Runtime type information
    TypeName("chemistryReductionMethodS");


    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        chemistryReductionMethodS,
        dictionary,
        (
            const IOdictionary& dict,
            STDACChemistryModel<CompType, ThermoType>& chemistry
        ),
        (dict, chemistry)
    );


    // Constructors

        //- Construct from components
        chemistryReductionMethodS
        (
            const IOdictionary& dict,
            STDACChemistryModel<CompType, ThermoType>& chemistry
        );


    // Selector

        static autoPtr<chemistryReductionMethodS<CompType, ThermoType>> New
        (
            const IOdictionary& dict,
            STDACChemistryModel<CompType, ThermoType>& chemistry
        );


    //- Destructor
    virtual ~chemistryReductionMethodS();


    // Member Functions

        //- Is mechanism reduction active?
        inline bool active() const;

        //- Is performance data logging enabled?
        inline bool log() const;

        //- Return the active species
        inline const List<bool>& activeSpecies() const;

        //- Return the number of active species
        inline label NsSimp();

        //- Return the initial number of species
        inline label nSpecie();

        //- Return the tolerance
        inline scalar tolerance() const;

        //- Reduce the mechanism
        virtual void reduceMechanism
        (
            const scalarField &c,
            const scalar T,
            const scalar p
        ) = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "chemistryReductionMethodSI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "chemistryReductionMethodS.C"
    #include "chemistryReductionMethodSNew.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
