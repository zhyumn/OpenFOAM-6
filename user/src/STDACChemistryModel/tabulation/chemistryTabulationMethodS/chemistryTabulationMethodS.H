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
    Foam::chemistryTabulationMethodS

Description
    An abstract class for chemistry tabulation.

SourceFiles
    chemistryTabulationMethodS.C

\*---------------------------------------------------------------------------*/

#ifndef chemistryTabulationMethodS_H
#define chemistryTabulationMethodS_H

#include "IOdictionary.H"
#include "scalarField.H"
#include "Switch.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
class STDACChemistryModel;

/*---------------------------------------------------------------------------*\
                           Class chemistryTabulationMethodS Declaration
\*---------------------------------------------------------------------------*/

template<class CompType, class ThermoType>
class chemistryTabulationMethodS
{

protected:

    const dictionary& dict_;

    const dictionary coeffsDict_;

    //- Is tabulation active?
    Switch active_;

    //- Switch to select performance logging
    Switch log_;

    STDACChemistryModel<CompType, ThermoType>& chemistry_;

    scalar tolerance_;


public:

    //- Runtime type information
    TypeName("chemistryTabulationMethodS");


    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        chemistryTabulationMethodS,
        dictionary,
        (
            const dictionary& dict,
            STDACChemistryModel<CompType, ThermoType>& chemistry
        ),
        (dict, chemistry)
    );


    // Constructors

        //- Construct from components
        chemistryTabulationMethodS
        (
            const dictionary& dict,
            STDACChemistryModel<CompType, ThermoType>& chemistry
        );


    // Selectors

        static autoPtr<chemistryTabulationMethodS> New
        (
            const IOdictionary& dict,
            STDACChemistryModel<CompType, ThermoType>& chemistry
        );


    //- Destructor
    virtual ~chemistryTabulationMethodS();


    // Member Functions

        inline bool active()
        {
            return active_;
        }

        inline bool log()
        {
            return active_ && log_;
        }

        inline bool variableTimeStep()
        {
            return chemistry_.variableTimeStep();
        }

        inline scalar tolerance() const
        {
            return tolerance_;
        }

        virtual label size() = 0;

        virtual void writePerformance() = 0;

        // Retrieve function: (only virtual here)
        // Try to retrieve a stored point close enough (according to tolerance)
        // to a stored point. If successful, it returns true and store the
        // results in RphiQ, i.e. the result of the integration of phiQ
        virtual bool retrieve
        (
             const scalarField& phiQ,
             scalarField& RphiQ
        ) = 0;

        // Add function: (only virtual here)
        // Add information to the tabulation algorithm. Give the reference for
        // future retrieve (phiQ) and the corresponding result (RphiQ).
        virtual label add
        (
            const scalarField& phiQ,
            const scalarField& RphiQ,
            const scalar rho,
            const scalar deltaT
        ) = 0;

        // Update function: (only virtual here)
        // The underlying structure of the tabulation is updated/cleaned
        // to increase the performance of the retrieve
        virtual bool update() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "chemistryTabulationMethodS.C"
    #include "chemistryTabulationMethodSNew.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
