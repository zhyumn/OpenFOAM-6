/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::chemistryReaderS

Description
    Abstract class for reading chemistry

SourceFiles
    chemistryReaderS.C

\*---------------------------------------------------------------------------*/

#ifndef chemistryReaderS_H
#define chemistryReaderS_H

#include "typeInfo.H"
#include "specieElement.H"
#include "Reaction.H"
#include "ReactionList.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

typedef HashTable<List<specieElement>> speciesCompositionTable;


/*---------------------------------------------------------------------------*\
                      Class chemistryReaderS Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class chemistryReaderS
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        chemistryReaderS(const chemistryReaderS&);

        //- Disallow default bitwise assignment
        void operator=(const chemistryReaderS&);


public:

    //- Runtime type information
    TypeName("chemistryReaderS");

    //- The type of thermo package the reader was instantiated for
    typedef ThermoType thermoType;


    // Constructors

        //- Construct null
        chemistryReaderS()
        {}


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            chemistryReaderS,
            dictionary,
            (
                const dictionary& thermoDict,
                speciesTable& species
            ),
            (thermoDict, species)
        );


    // Selectors

        //- Select constructed from dictionary
        static autoPtr<chemistryReaderS> New
        (
            const dictionary& thermoDict,
            speciesTable& species
        );


    //- Destructor
    virtual ~chemistryReaderS()
    {}


    // Member Functions

        //- Return access to the list of species
        virtual const speciesTable& species() const = 0;

        //- Table of species composition
        virtual const speciesCompositionTable& specieComposition() const = 0;

        //- Return access to the thermo packages
        virtual const HashPtrTable<ThermoType>& speciesThermo() const = 0;

        //- Return access to the list of reactions
        virtual const ReactionList<ThermoType>& reactions() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "chemistryReaderS.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#define makeChemistryReaderS(Thermo)                                            \
    defineTemplateTypeNameAndDebug(chemistryReaderS<Thermo>, 0);                \
    defineTemplateRunTimeSelectionTable(chemistryReaderS<Thermo>, dictionary)


#define makeChemistryReaderSType(Reader, Thermo)                                \
    defineNamedTemplateTypeNameAndDebug(Reader<Thermo>, 0);                    \
    chemistryReaderS<Thermo>::adddictionaryConstructorToTable<Reader<Thermo>> \
        add##Reader##Thermo##ConstructorToTable_


// for non-templated chemistry readers
#define addChemistryReaderSType(Reader, Thermo)                                 \
    defineTypeNameAndDebug(Reader, 0);                                         \
    chemistryReaderS<Thermo>::adddictionaryConstructorToTable<Reader>           \
        add##Reader##Thermo##ConstructorToTable_


// for templated chemistry readers
#define addTemplateChemistryReaderSType(Reader, Thermo)                         \
    defineNamedTemplateTypeNameAndDebug(Reader, 0);                            \
    chemistryReaderS<Thermo>::adddictionaryConstructorToTable<Reader>           \
        add##Reader##Thermo##ConstructorToTable_


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
