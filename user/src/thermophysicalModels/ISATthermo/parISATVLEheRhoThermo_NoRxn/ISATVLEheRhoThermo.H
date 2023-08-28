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
#include "SUPstream.H"
#define DateSize1 4
#define DateSize2 4
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
        volScalarField rho_d;
        //volScalarField rho_G_;
        //volScalarField Z_;

        //- Store list of mass diffusion coefficients
        PtrList<volScalarField> Dimix_; //

        PtrList<volScalarField> heList_;
        PtrList<volScalarField> WList_;
        //PtrList<volScalarField> Y_G_List_;

        bool inviscid_;
        //bool DF_;
        word scheme_;
        label nloop;
        const clockTime clockTime_;

        autoPtr<OFstream> cpuISAT_VLE_;
        bool ISATlog_;
        bool loadBalance_;

        volScalarField *FCcell;

        struct jobInput
        {

            scalar p0;
            scalar T0;
            FixedList<scalar, DateSize1> in;
            FixedList<scalar, DateSize2> out;
        };

        struct jobArray
        {
            label rank;
            label next;
            label first_i;
            label size_i;
            jobInput jobs[1];
        };
        label n_block;
        int Batch_Size;
        SUPstream::shared_data2<jobArray> ja;

        bool spare;
        template <typename T>
        struct padded_atomic
        {
            std::atomic<T> var;
            char m_pad[CacheLineSize - sizeof(std::atomic<T>)];
        };
        SUPstream::shared_data<padded_atomic<int>> sspare_cpu;
        std::atomic<int> &spare_cpu;
        int nfinished_block;
        SUPstream::shared_data<std::atomic<int>> sfinished_head;
        std::atomic<int> *finished_head;

        SUPstream::mpi_mutex l_empty;
        SUPstream::mpi_mutex l_filled;
        SUPstream::mpi_mutex l_finished;

        SUPstream::mpi_mutex l_receiver;
        SUPstream::mpi_mutex l_sender;

        SUPstream::shared_data<std::atomic<int>> sempty_head, sempty_tail;
        std::atomic<int> &empty_head, &empty_tail;
        SUPstream::shared_data<std::atomic<int>> sfilled_head, sfilled_tail;
        std::atomic<int> &filled_head, &filled_tail;

        int &rank;
        label n_cpu;

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
