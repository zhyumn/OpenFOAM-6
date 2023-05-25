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

#include <stdlib.h>
#include <stdio.h>
#include "/scratch.nike/srini237/OpenFOAM/temp/OpenFOAM-6/user/tutorials/python/UVFlash/ANN/onnx_runtime/onnxruntime-linux-x64-1.14.1/include/onnxruntime_cxx_api.h"
#include <vector>
#include <iostream>

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

        volScalarField *FCcell; 

        // bool noVLE_;

        // Private Member Functions

        //- Calculate the thermo variables
        void calculate_init();

        void calculate();

        //- Construct as copy (not implemented)
        ISATVLEheRhoThermo(const ISATVLEheRhoThermo<BasicPsiThermo, MixtureType> &);
              
        //struct ANN_predict{
        // TF_Graph* Graph; 
        // TF_Status* Status; 
        // TF_SessionOptions* SessionOpts; 
        // TF_Buffer* RunOpts;
        // TF_Session* Session;

        // TF_Graph* Graph_BC; 
        // TF_Status* Status_BC; 
        // TF_SessionOptions* SessionOpts_BC; 
        // TF_Buffer* RunOpts_BC;
        // TF_Session* Session_BC;
       // };

        Ort::Env env; //(ORT_LOGGING_LEVEL_WARNING, "example");
        Ort::SessionOptions session_options;
        Ort::AllocatorWithDefaultOptions allocator;
        // const char* model_dir = "/scratch.nike/srini237/OpenFOAM/temp/OpenFOAM-6/user/tutorials/python/UVFlash/ANN/OnnxModel.onnx";
        //Ort::Session session{nullptr}; //(Ort::Env env, model_dir, session_options);
        Ort::Session session{env,"/scratch.nike/srini237/OpenFOAM/temp/OpenFOAM-6/user/tutorials/python/UVFlash/ANN/OnnxModel2.onnx",Ort::SessionOptions{nullptr}};
        // auto* inputname = "input";
        // auto* outputname0 = "T";//session.GetOutputName(0,env);
        // auto* outputname1 = "p";//session.GetOutputName(1);
        // auto* outputname2 = "phi";//session.GetOutputName(2);
        // auto* outputname3 = "c";//session.GetOutputName(3);
        // auto memoryinfo; //= Ort::MemoryInfo::CreateCpu(OrtDeviceAllocator,OrtMemTypeCPU);
        // auto inputOnnxTensor;
        // const char* inputNames[] = {"input"};
        // const char* outputNames[] = {"T","P","phi","c"};
        // auto outputValues;

        //ANN_predict ANN(saved_model_dir,tags,ntags);

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

        Ort::Session ANN_load();
        std::tuple<scalar,scalar,scalar,scalar> ANN_predict(float u, float v, float x1, float x2);

        // void ANN_load_BC(char* saved_model_dir_BC, char* tags_BC,int ntags_BC);
        // std::tuple<scalar,scalar,scalar,scalar> ANN_predict_BC(float T, float p, float x1, float x2);
        //namespace std{
    //}
    };
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void NoOpDeallocator(void* data, size_t a, void* b) {}
#include "ANN.C"

#ifdef NoRepository
#include "ISATVLEheRhoThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
