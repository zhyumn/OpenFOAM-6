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
    Foam::chemistryReductionMethodPT

Description
    An abstract class for reducing chemical mechanism

SourceFiles
    chemistryReductionMethodPT.C

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
inline bool
Foam::chemistryReductionMethodPT<CompType, ThermoType>::active() const
{
    return active_;
}


template<class CompType, class ThermoType>
inline bool
Foam::chemistryReductionMethodPT<CompType, ThermoType>::log() const
{
    return active_ && log_;
}


template<class CompType, class ThermoType>
inline const Foam::List<bool>&
Foam::chemistryReductionMethodPT<CompType, ThermoType>::activeSpecies() const
{
    return activeSpecies_;
}

template<class CompType, class ThermoType>
inline Foam::label
Foam::chemistryReductionMethodPT<CompType, ThermoType>::NsSimp()
{
    return NsSimp_;
}


template<class CompType, class ThermoType>
inline Foam::label
Foam::chemistryReductionMethodPT<CompType, ThermoType>::nSpecie()
{
    return nSpecie_;
}


template<class CompType, class ThermoType>
inline Foam::scalar
Foam::chemistryReductionMethodPT<CompType, ThermoType>::tolerance() const
{
    return tolerance_;
}


// ************************************************************************* //
