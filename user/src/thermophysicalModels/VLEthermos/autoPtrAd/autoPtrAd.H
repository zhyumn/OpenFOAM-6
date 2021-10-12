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
    Foam::autoPtrAd

Description

SourceFiles
    autoPtrAdI.H
    autoPtrAd.C
    autoPtrAdIO.C

\*---------------------------------------------------------------------------*/

#ifndef autoPtrAd_H
#define autoPtrAd_H

#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes


// Forward declaration of friend functions and operators



/*---------------------------------------------------------------------------*\
                         Class autoPtrAd Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class autoPtrAd
:
    public autoPtr
{
    // Private data



    // Private Member Functions

        //- Disallow default bitwise copy construct


        //- Disallow default bitwise assignment


public:

    // Static data members

        //- Static data someStaticData


    // Constructors

        //- Construct null
        inline explicit autoPtrAd(T* p= nullptr):autoPtr<T>(p){};

        //- Construct as copy by transferring pointer to this autoPtr and
        //  setting the arguments pointer to nullptr
        inline autoPtrAd(const autoPtr<T>& ap):autoPtr<T>(ap){};

        //- Construct either by transferring pointer or cloning.
        //  Should only be called with type that supports cloning
        inline autoPtrAd(const autoPtr<T>& ap, const bool reuse):autoPtr<T>(ap,reuse){};


    // Selectors

        //- Select null constructed
        inline void reset(const autoPtr<T>&);


    //- Destructor
    ~autoPtrAd();


    // Member Functions

    // Member Operators



    // Friend Functions

    // Friend Operators

    // IOstream Operators

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "autoPtrAdI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
