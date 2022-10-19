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
        volScalarField rho_G_;

        //- Store list of mass diffusion coefficients
        PtrList<volScalarField> Dimix_; //

        PtrList<volScalarField> heList_;
        PtrList<volScalarField> WList_;
        PtrList<volScalarField> Y_G_List_;

        bool inviscid_;
        bool noVLE_;

        

        // Private Member Functions

        //- Calculate the thermo variables
        void calculate_init();

        void calculate();

        //- Construct as copy (not implemented)
        ISATVLEheRhoThermo(const ISATVLEheRhoThermo<BasicPsiThermo, MixtureType> &);

        const double R_r = 8.314 * 1E7; //CGS
        double T_r, P_r, rho_r;
        double c_r[6],Ygas_r[6];
        double reactRate_r[6];
        double W_r[6];

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

        // Member functions

        //- Update properties
        virtual void correct();

        void customChemistry(scalar temperature, scalar pressure, scalar density);
        void setThermo(scalar temperature, scalar pressure, scalar density);
        void setY(scalar Y, label speciesLabel);
        void setConc();
        void setRate();
        scalar RR(label speciesLabel);
        scalar Qdot();
        scalar gasPhaseRho(scalar pressure, scalar temperature);
        scalar gasPhaseY(label celli, label speciesLabel);
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
