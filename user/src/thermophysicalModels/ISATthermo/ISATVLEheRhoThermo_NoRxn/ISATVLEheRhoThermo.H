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
    Foam::ISATVLEheRhoThermo

Description
    Energy for a mixture based on compressibility

SourceFiles
    ISATVLEheRhoThermo.C

\*---------------------------------------------------------------------------*/

#ifndef ISATVLEheRhoThermo_H
#define ISATVLEheRhoThermo_H

//#include "psiThermo.H"
#include <tuple>
#include "heThermo.H"
#include "clockTime.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    /*---------------------------------------------------------------------------*\
                             Class ISATVLEheRhoThermo Declaration
    \*---------------------------------------------------------------------------*/

    template <class BasicPsiThermo, class MixtureType>
    class ISATVLEheRhoThermo
        : public heThermo<BasicPsiThermo, MixtureType>
    {

        volScalarField vaporfrac_;
        volScalarField soundspeed_;
        volScalarField kappa_;
        //volScalarField rho_G_;
        //volScalarField Z_;

        //- Store list of mass diffusion coefficients
        PtrList<volScalarField> Dimix_; //

        PtrList<volScalarField> heList_;
        PtrList<volScalarField> WList_;
        //PtrList<volScalarField> Y_G_List_;

        bool inviscid_;
        bool DF_;
        label nloop;
        const clockTime clockTime_;

        autoPtr<OFstream> cpuISAT_VLE_;
        bool ISATlog_;

        // bool noVLE_;

        // Private Member Functions

        //- Calculate the thermo variables
        void calculate_init();

        void calculate();

        //- Construct as copy (not implemented)
        ISATVLEheRhoThermo(const ISATVLEheRhoThermo<BasicPsiThermo, MixtureType> &);

    public:
        //- Runtime type information
        TypeName("ISATVLEheRhoThermo");

        // Constructors

        //- Construct from mesh and phase name
        ISATVLEheRhoThermo(
            const fvMesh &,
            const word &phaseName);

        //- Destructor
        virtual ~ISATVLEheRhoThermo();

        inline autoPtr<OFstream> logFile(const word &name, const fvMesh &mesh) const;

        // Member functions

        //- Update properties
        virtual void correct();
    };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "ISATVLEheRhoThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
