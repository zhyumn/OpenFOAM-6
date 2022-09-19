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

Application
    rhoCentralFoam

Description
    Density-based compressible flow solver based on central-upwind schemes of
    Kurganov and Tadmor.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "rhoReactionThermo.H"
#include "turbulentFluidThermoModel.H"
#include "fixedRhoFvPatchScalarField.H"
#include "directionInterpolate.H"
#include "localEulerDdtScheme.H"
#include "CombustionModel.H"
#include "fvcSmooth.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#define NO_CONTROL
#include "postProcess.H"

#include "setRootCaseLists.H"
#include "createTime.H"
#include "createMesh.H"
#include "createFields.H"
#include "createFieldRefs.H"
#include "createTimeControls.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "readFluxScheme.H"

    dimensionedScalar v_zero("v_zero", dimVolume / dimTime, 0.0);

    // Courant numbers used to adjust the time-step
    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;

    /* int procID;

     if (UPstream::parRun())
     {
         procID = UPstream::myProcNo();
     }
     */

    Info << "\nStarting time loop\n"
         << endl;

    /*if (UPstream::parRun())
    {
        Pout << procID << endl;
    }
    */

    while (runTime.run())
    {
        // --- Directed interpolation of primitive fields onto faces

        surfaceScalarField rho_pos(interpolate(rho, pos));
        surfaceScalarField rho_neg(interpolate(rho, neg));

        surfaceVectorField rhoU_pos(interpolate(rhoU, pos, U.name()));
        surfaceVectorField rhoU_neg(interpolate(rhoU, neg, U.name()));

        volScalarField rPsi("rPsi", 1.0 / psi);
        surfaceScalarField rPsi_pos(interpolate(rPsi, pos, T.name()));
        surfaceScalarField rPsi_neg(interpolate(rPsi, neg, T.name()));

        surfaceScalarField gammaStar_pos(interpolate(gammaStar, pos));
        surfaceScalarField gammaStar_neg(interpolate(gammaStar, neg));

        surfaceScalarField eStar_pos(interpolate(eStar, pos));
        surfaceScalarField eStar_neg(interpolate(eStar, neg));

        surfaceScalarField e_pos(interpolate(e, pos, T.name()));
        surfaceScalarField e_neg(interpolate(e, neg, T.name()));

        surfaceVectorField U_pos("U_pos", rhoU_pos / rho_pos);
        surfaceVectorField U_neg("U_neg", rhoU_neg / rho_neg);

        surfaceScalarField p_pos("p_pos", rho_pos * rPsi_pos);
        surfaceScalarField p_neg("p_neg", rho_neg * rPsi_neg);

        surfaceScalarField e_pos_pos("e_pos_pos", p_pos / (rho_pos * (gammaStar_pos - 1)) + eStar_pos); // = e_pos
        surfaceScalarField e_pos_neg("e_pos_neg", p_neg / (rho_neg * (gammaStar_pos - 1)) + eStar_pos);

        surfaceScalarField e_neg_pos("e_neg_pos", p_pos / (rho_pos * (gammaStar_neg - 1)) + eStar_neg);
        surfaceScalarField e_neg_neg("e_neg_neg", p_neg / (rho_neg * (gammaStar_neg - 1)) + eStar_neg); // = e_neg

        surfaceScalarField phiv_pos("phiv_pos", U_pos & mesh.Sf());
        surfaceScalarField phiv_neg("phiv_neg", U_neg & mesh.Sf());

        surfaceScalarField phiv_pos_pos("phiv_pos_pos", U_pos & mesh.Sf());
        surfaceScalarField phiv_pos_neg("phiv_pos_neg", U_neg & mesh.Sf());

        surfaceScalarField phiv_neg_pos("phiv_neg_pos", U_pos & mesh.Sf());
        surfaceScalarField phiv_neg_neg("phiv_neg_neg", U_neg & mesh.Sf());

        // volScalarField c("c", sqrt(thermo.Cp() / thermo.Cv() * rPsi)); //REalgas
        // volScalarField c("c", cc); //REalgas
        // Info << "!!!!!!!!!!!c=" << cc[0] << endl;
        // FatalErrorInFunction << "break " << exit(FatalError);
        surfaceScalarField cSf_pos(
            "cSf_pos",
            interpolate(c, pos, T.name()) * mesh.magSf());
        surfaceScalarField cSf_neg(
            "cSf_neg",
            interpolate(c, neg, T.name()) * mesh.magSf());

        surfaceScalarField cSf_pos_pos(
            "cSf_pos_pos",
            sqrt(gammaStar_pos * rPsi_pos) * mesh.magSf()); // = cSf_pos
        surfaceScalarField cSf_pos_neg(
            "cSf_pos_neg",
            sqrt(gammaStar_pos * rPsi_neg) * mesh.magSf());

        surfaceScalarField cSf_neg_pos(
            "cSf_neg_pos",
            sqrt(gammaStar_neg * rPsi_pos) * mesh.magSf());
        surfaceScalarField cSf_neg_neg(
            "cSf_neg_neg",
            sqrt(gammaStar_neg * rPsi_neg) * mesh.magSf()); // = cSf_neg

        surfaceScalarField ap(
            "ap",
            max(max(phiv_pos + cSf_pos, phiv_neg + cSf_neg), v_zero));
        surfaceScalarField am(
            "am",
            min(min(phiv_pos - cSf_pos, phiv_neg - cSf_neg), v_zero));

        surfaceScalarField ap_pos(
            "ap_pos",
            max(max(phiv_pos + cSf_pos_pos, phiv_neg + cSf_pos_neg), v_zero));

        surfaceScalarField ap_neg(
            "ap_neg",
            max(max(phiv_pos + cSf_neg_pos, phiv_neg + cSf_neg_neg), v_zero));
        surfaceScalarField am_pos(
            "am_pos",
            min(min(phiv_pos - cSf_pos_pos, phiv_neg - cSf_pos_neg), v_zero));
        surfaceScalarField am_neg(
            "am_neg",
            min(min(phiv_pos - cSf_neg_pos, phiv_neg - cSf_neg_neg), v_zero));

        surfaceScalarField a_pos("a_pos", ap / (ap - am));

        surfaceScalarField a_pos_pos("a_pos_pos", ap_pos / (ap_pos - am_pos));

        surfaceScalarField a_neg_pos("a_neg_pos", ap_neg / (ap_neg - am_neg));

        surfaceScalarField amaxSf("amaxSf", max(mag(am), mag(ap)));

        surfaceScalarField amaxSf_pos("amaxSf_pos", max(mag(am_pos), mag(ap_pos)));

        surfaceScalarField amaxSf_neg("amaxSf_neg", max(mag(am_neg), mag(ap_neg)));

        surfaceScalarField aSf("aSf", am * a_pos);

        surfaceScalarField aSf_pos("aSf_pos", am_pos * a_pos_pos);

        surfaceScalarField aSf_neg("aSf_neg", am_neg * a_neg_pos);

        if (fluxScheme == "Tadmor")
        {
            aSf = -0.5 * amaxSf;
            a_pos = 0.5;

            aSf_pos = -0.5 * amaxSf_pos;
            a_pos_pos = 0.5;

            aSf_neg = -0.5 * amaxSf_neg;
            a_neg_pos = 0.5;
        }

        surfaceScalarField a_neg("a_neg", 1.0 - a_pos);

        surfaceScalarField a_pos_neg("a_pos_neg", 1.0 - a_pos_pos);

        surfaceScalarField a_neg_neg("a_neg_neg", 1.0 - a_neg_pos);

        phiv_pos *= a_pos;
        phiv_neg *= a_neg;

        phiv_pos_pos *= a_pos_pos;
        phiv_pos_neg *= a_pos_neg;

        phiv_neg_pos *= a_neg_pos;
        phiv_neg_neg *= a_neg_neg;

        surfaceScalarField aphiv_pos("aphiv_pos", phiv_pos - aSf);
        surfaceScalarField aphiv_neg("aphiv_neg", phiv_neg + aSf);

        surfaceScalarField aphiv_pos_pos("aphiv_pos_pos", phiv_pos_pos - aSf_pos);
        surfaceScalarField aphiv_pos_neg("aphiv_pos_neg", phiv_pos_neg + aSf_pos);

        surfaceScalarField aphiv_neg_pos("aphiv_neg_pos", phiv_neg_pos - aSf_neg);
        surfaceScalarField aphiv_neg_neg("aphiv_neg_neg", phiv_neg_neg + aSf_neg);

        // Reuse amaxSf for the maximum positive and negative fluxes
        // estimated by the central scheme
        amaxSf = max(mag(aphiv_pos), mag(aphiv_neg));

        amaxSf_pos = max(mag(aphiv_pos_pos), mag(aphiv_pos_neg));

        amaxSf_neg = max(mag(aphiv_neg_pos), mag(aphiv_neg_neg));

#include "centralCourantNo.H"
#include "readTimeControls.H"

        if (LTS)
        {
#include "setRDeltaT.H"
        }
        else
        {
#include "setDeltaT.H"
        }

        runTime++;

        Info << "Time = " << runTime.timeName() << nl << endl;
        scalar t = runTime.value();

        phi = aphiv_pos * rho_pos + aphiv_neg * rho_neg;

        surfaceScalarField phi_pos("phi_pos", aphiv_pos_pos * rho_pos + aphiv_pos_neg * rho_neg);
        surfaceScalarField phi_neg("phi_neg", aphiv_neg_pos * rho_pos + aphiv_neg_neg * rho_neg);

        surfaceVectorField phiUp(
            (aphiv_pos * rhoU_pos + aphiv_neg * rhoU_neg) + (a_pos * p_pos + a_neg * p_neg) * mesh.Sf());

        surfaceVectorField phiUp_pos(
            (aphiv_pos_pos * rhoU_pos + aphiv_pos_neg * rhoU_neg) + (a_pos_pos * p_pos + a_pos_neg * p_neg) * mesh.Sf());

        surfaceVectorField phiUp_neg(
            (aphiv_neg_pos * rhoU_pos + aphiv_neg_neg * rhoU_neg) + (a_neg_pos * p_pos + a_neg_neg * p_neg) * mesh.Sf());

        surfaceScalarField phiEp(
            "phiEp",
            aphiv_pos * (rho_pos * (e_pos + 0.5 * magSqr(U_pos)) + p_pos) + aphiv_neg * (rho_neg * (e_neg + 0.5 * magSqr(U_neg)) + p_neg) + aSf * p_pos - aSf * p_neg);
        surfaceScalarField phiEp_pos(
            "phiEp_pos",
            aphiv_pos_pos * (rho_pos * (e_pos_pos + 0.5 * magSqr(U_pos)) + p_pos) + aphiv_pos_neg * (rho_neg * (e_pos_neg + 0.5 * magSqr(U_neg)) + p_neg) + aSf_pos * p_pos - aSf_pos * p_neg);

        surfaceScalarField phiEp_neg(
            "phiEp_pos",
            aphiv_neg_pos * (rho_pos * (e_neg_pos + 0.5 * magSqr(U_pos)) + p_pos) + aphiv_neg_neg * (rho_neg * (e_neg_neg + 0.5 * magSqr(U_neg)) + p_neg) + aSf_neg * p_pos - aSf_neg * p_neg);

        volScalarField muEff("muEff", turbulence->muEff());
        volTensorField tauMC("tauMC", muEff * dev2(Foam::T(fvc::grad(U))));

        if (divScheme == "Doubleflux")
        {
            // --- Solve density
            solve(fvm::ddt(rho) + fvc::div_doubleflux(phi_pos, phi_neg));

            // --- Solve momentum
            solve(fvm::ddt(rhoU) + fvc::div_doubleflux(phiUp_pos, phiUp_neg));
        }
        else if (divScheme == "Conservativeflux")
        {

            // --- Solve density
            solve(fvm::ddt(rho) + fvc::div(phi));

            // --- Solve momentum
            solve(fvm::ddt(rhoU) + fvc::div(phiUp));
        }

        U.ref() =
            rhoU() / rho();
        U.correctBoundaryConditions();

        if (!inviscid)
        {
            solve(
                fvm::ddt(rho, U) - fvc::ddt(rho, U) - fvm::laplacian(muEff, U) - fvc::div(tauMC));
            rhoU = rho * U;
        }
        if (!inviscid)
        {
            sumHeatDiffusion *= 0.;
            sumHeatDiffusion2 *= 0.;

            forAll(Y, k)
            {
                sumHeatDiffusion += fvc::laplacian(alpha * hei[k], Y[k]);
                sumHeatDiffusion2 += fvc::div(hei[k] * rho * YVi[k]);
            }
        }
#include "YEqn.H"

        // --- Solve energy
        surfaceScalarField sigmaDotU(
            "sigmaDotU",
            (
                fvc::interpolate(muEff) * mesh.magSf() * fvc::snGrad(U) + fvc::dotInterpolate(mesh.Sf(), tauMC)) &
                (a_pos * U_pos + a_neg * U_neg));

        fvScalarMatrix rhoEEqn(
            fvm::ddt(rhoE));

        if (divScheme == "Doubleflux")
        {
            rhoEEqn += fvc::div_doubleflux(phiEp_pos, phiEp_neg);
            /*solve
            (
                fvm::ddt(rhoE)
                //+ fvc::div(phiEp)
                + fvc::div_doubleflux(phiEp_pos, phiEp_neg)
                - fvc::div(sigmaDotU)
            );*/
        }
        else if (divScheme == "Conservativeflux")
        {
            rhoEEqn += fvc::div(phiEp);
            /*solve
            (
                fvm::ddt(rhoE)
                + fvc::div(phiEp)
                //+ fvc::div_doubleflux(phiEp_pos, phiEp_neg)
                - fvc::div(sigmaDotU)
            );*/
        }
        if (!inviscid)
        {
            rhoEEqn -= fvc::div(sigmaDotU);
        }
        solve(rhoEEqn);

        e = rhoE / rho - 0.5 * magSqr(U);
        e.correctBoundaryConditions();
        // rho_thermo=rho;
        // thermo.correct();

        if (!inviscid)
        {
            volScalarField &he = thermo.he();
            fvScalarMatrix EEqn(
                fvm::ddt(rho, e) - fvc::ddt(rho, e)
                    //- fvm::laplacian(turbulence->alphaEff(), e)
                    //- fvm::laplacian(kappa / Cp, he) //Todo kappa Cp
                    - fvm::laplacian(alpha, he) ==
                -sumHeatDiffusion - sumHeatDiffusion2
                //== fvOptions(rho, e)
            );
            /*forAll(Y, k)
            {
                //EEqn -= fvc::laplacian((kappa / Cp) * hei[k], Y[k]);//Todo hei
                EEqn -= fvc::laplacian((alpha)*hei[k], Y[k]);//Todo hei
                EEqn -= fvc::div(hei[k] * rho * YVi[k]);//Todo YVi
            }*/
            EEqn.relax();
            // fvOptions.constrain(EEqn);
            EEqn.solve();
            // fvOptions.correct(e);
            // thermo.correct();
        }

        if (divScheme == "Doubleflux")
        {
            p.ref() = (e() - eStar()) * rho() * (gammaStar() - 1);
        }
        else if (divScheme == "Conservativeflux")
        {
            p.ref() =
                rho() / psi();
        }
        else
        {
            FatalErrorInFunction << "divScheme must be Doubleflux or Conservativeflux." << exit(FatalError);
        }
        p.correctBoundaryConditions();
        thermo.correct();
        T.correctBoundaryConditions();
        thermo.correct();
        rho.boundaryFieldRef() == psi.boundaryField() * p.boundaryField();
        rhoU.boundaryFieldRef() == rho.boundaryField() * U.boundaryField();

        volScalarField Yt(0.0 * Y[0]);
        forAll(Y, i)
        {
            if (i != inertIndex && composition.active(i))
            {
                volScalarField &Yi = Y[i];
                volScalarField &rhoYi = rhoY[i];

                rhoYi.boundaryFieldRef() = rho.boundaryField() * Yi.boundaryField();

                Yi.max(0.0);
                Yi.min(1.0);

                Yt += Yi;
            }
        }
        Y[inertIndex] = scalar(1) - Yt;
        Y[inertIndex].max(0.0);

        //rhoE.boundaryFieldRef() ==
        //    rho.boundaryField() *
        //        (e.boundaryField() + 0.5 * magSqr(U.boundaryField()));

        rhoE = rho * (e + 0.5 * magSqr(U));

        turbulence->correct();
        gammaStar = rho * c * c / p;
        eStar = e - p / (rho * (gammaStar - 1));
        Wmix = thermo.W();
        forAll(Dimix, i)
        {
            Dimix[i] = Dimix_t[i];
            // Dimix[i] = thermo.Dimix(i);
            hei[i] = hei_t[i];
            // hei[i] = thermo.hei(i);
        }
        // U*=0.1;
        rho_write = rho;
        runTime.write();

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << endl;
    }

    Info << "End\n"
         << endl;

    return 0;
}

// ************************************************************************* //
