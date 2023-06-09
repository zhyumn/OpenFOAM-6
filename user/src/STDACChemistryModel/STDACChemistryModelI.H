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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
inline bool
Foam::STDACChemistryModel<ReactionThermo, ThermoType>::variableTimeStep() const
{
    return variableTimeStep_;
}


template<class ReactionThermo, class ThermoType>
inline Foam::label
Foam::STDACChemistryModel<ReactionThermo, ThermoType>::timeSteps() const
{
    return timeSteps_;
}


template<class ReactionThermo, class ThermoType>
inline Foam::autoPtr<Foam::OFstream>
Foam::STDACChemistryModel<ReactionThermo, ThermoType>::
logFile(const word& name) const
{
    mkDir(this->mesh().time().path()/"STDAC"/this->group());
    return autoPtr<OFstream>
    (
        new OFstream
        (
            this->mesh().time().path()/"STDAC"/this->group()/name
        )
    );
}


template<class ReactionThermo, class ThermoType>
inline Foam::PtrList<Foam::volScalarField>&
Foam::STDACChemistryModel<ReactionThermo, ThermoType>::Y()
{
    return this->Y_;
}


template<class ReactionThermo, class ThermoType>
inline
Foam::autoPtr<Foam::chemistryReductionMethodS<ReactionThermo, ThermoType>>&
Foam::STDACChemistryModel<ReactionThermo, ThermoType>::mechRed()
{
    return mechRed_;
}


template<class ReactionThermo, class ThermoType>
inline void Foam::STDACChemistryModel<ReactionThermo, ThermoType>::setActive
(
    const label i
)
{
    this->thermo().composition().setActive(i);
}


template<class ReactionThermo, class ThermoType>
inline bool Foam::STDACChemistryModel<ReactionThermo, ThermoType>::active
(
    const label i
) const
{
    return this->thermo().composition().active(i);
}


template<class ReactionThermo, class ThermoType>
inline void Foam::STDACChemistryModel<ReactionThermo, ThermoType>::
setNsDAC(const label newNsDAC)
{
    NsDAC_ = newNsDAC;
}


template<class ReactionThermo, class ThermoType>
inline void Foam::STDACChemistryModel<ReactionThermo, ThermoType>::
setNSpecie(const label newNs)
{
    this->nSpecie_ = newNs;
}


template<class ReactionThermo, class ThermoType>
inline Foam::DynamicList<Foam::label>&
Foam::STDACChemistryModel<ReactionThermo, ThermoType>::
simplifiedToCompleteIndex()
{
    return simplifiedToCompleteIndex_;
}


template<class ReactionThermo, class ThermoType>
inline Foam::Field<Foam::label>&
Foam::STDACChemistryModel<ReactionThermo, ThermoType>::
completeToSimplifiedIndex()
{
    return completeToSimplifiedIndex_;
}


template<class ReactionThermo, class ThermoType>
inline const Foam::Field<Foam::label>&
Foam::STDACChemistryModel<ReactionThermo, ThermoType>::
completeToSimplifiedIndex() const
{
    return completeToSimplifiedIndex_;
}


template<class ReactionThermo, class ThermoType>
inline Foam::Field<bool>&
Foam::STDACChemistryModel<ReactionThermo, ThermoType>::reactionsDisabled()
{
    return reactionsDisabled_;
}


template<class ReactionThermo, class ThermoType>
inline Foam::scalarField&
Foam::STDACChemistryModel<ReactionThermo, ThermoType>::completeC()
{
    return completeC_;
}


template<class ReactionThermo, class ThermoType>
inline Foam::scalarField&
Foam::STDACChemistryModel<ReactionThermo, ThermoType>::simplifiedC()
{
    return simplifiedC_;
}


template<class ReactionThermo, class ThermoType>
inline Foam::List<Foam::List<Foam::specieElement>>&
Foam::STDACChemistryModel<ReactionThermo, ThermoType>::specieComp()
{
    return specieComp_;
}


template<class ReactionThermo, class ThermoType>
void Foam::STDACChemistryModel<ReactionThermo, ThermoType>::
resetTabulationResults()
{
    forAll(tabulationResults_, tabi)
    {
        tabulationResults_[tabi] = 2;
    }
}


// ************************************************************************* //
