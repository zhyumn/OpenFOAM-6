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

\*---------------------------------------------------------------------------*/

//#include "specie.H"
//#include "IOstreams.H"
#include "makeReactionThermo.H"
#include "multithermo.H"
#include "thermo.H"
#include "janafThermoAd.H"
#include "sensibleInternalEnergy.H"
#include "VLE.H"
#include "chungTransportMixture.H"
#include "chungTransport.H"
#include "multispecie.H"
#include "nLreactingMixture.H"
#include "VLEhePsiThermo.H"
#include "psiReactionThermo.H"
#include "foamChemistryReader.H"
#include "makeReactionUser.H"
//#include "Reaction.H"
//#include "makeReaction.H"
//#include "thermoPhysicsTypes.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{
    


    typedef 
    chungTransport
    <
        species::thermo
        <  
            janafThermoAd
            <
                PengRobinson<specie>>
            ,
            sensibleInternalEnergy
        >
    > VLEChungPREThermoPhysics;

    typedef Reaction<VLEChungPREThermoPhysics> VLEChungPREThermoPhysicsReaction;

    makeReactions(VLEChungPREThermoPhysics, VLEChungPREThermoPhysicsReaction)
    makeChemistryReader(VLEChungPREThermoPhysics);
    makeChemistryReaderType(foamChemistryReader, VLEChungPREThermoPhysics);
    template <class Thermo>
    using  Emultithermo = species::multithermo<Thermo, sensibleInternalEnergy>;
    template <class Thermo>
    using VLEthermoE = Emultithermo<VLE<chungTransportMixture<PengRobinsonMixture<multispecie<Thermo>>>>>;
    template <class Thermo>
    using  nLreactingMixtureEChungPR = nLreactingMixture<Thermo, VLEthermoE>;
    makeThermoPhysicsReactionThermos
    (
        psiThermo,
        psiReactionThermo,
        VLEhePsiThermo,
        nLreactingMixtureEChungPR,
        VLEChungPREThermoPhysics
    );

}


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //




// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


// ************************************************************************* //
