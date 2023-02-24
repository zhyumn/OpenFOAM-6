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
#include "minMod.H"
#include "rhoReactionThermo.H"
#include "turbulentFluidThermoModel.H"
#include "fixedRhoFvPatchScalarField.H"
#include "directionInterpolate.H"
#include "localEulerDdtScheme.H"
#include "CombustionModel.H"
#include "fvcSmooth.H"
#include "mathematicalConstants.H"
#include "thermodynamicConstants.H"
#include "clockTime.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#define NO_CONTROL
    Foam::argList::addBoolOption(
        "init",
        "initailize");

    Foam::argList::addBoolOption(
        "updateZero",
        "updateZero");
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

    Info << "\nStarting time loop\n"
         << endl;
    //Z.write();

    if (args.optionFound("init"))
    {
        double rho_bottom = rho[0];
        double rho_top = rho[rho.size() - 1];
        forAll(rho, i)
        {
            rho[i] = (Y[0][i] + Y[1][i]) / (Y[0][i] / rho_top + Y[1][i] / rho_bottom);
            //rho[i] = Y[0][i]*rho_top+Y[1][i]*rho_bottom;
        }
        thermo.correct();
        /*forAll (rho,i)
        {
            T[i]= p[i]*Wmix[i]/rho[i]/8.314;
        }*/
        T.write();
        rho.write();
        frac.write();
        return 0;
    }

    if (args.optionFound("updateZero"))
    {
        e.write();
        rho.write();
        frac.write();
        return 0;
    }
    int nloop = 0;

    double cputime;

    while (runTime.run())
    {
        cputime = 0;
        clockTime_.timeIncrement();
        // --- Directed interpolation of primitive fields onto faces
        //volScalarField rho_d("rho_d", rho);

        surfaceScalarField rho_pos(interpolate(rho, pos));
        surfaceScalarField rho_neg(interpolate(rho, neg));

        surfaceVectorField rhoU_pos(interpolate(rhoU, pos, U.name()));
        surfaceVectorField rhoU_neg(interpolate(rhoU, neg, U.name()));

        volScalarField rPsi("rPsi", 1.0 / psi);
        surfaceScalarField rPsi_pos(interpolate(rPsi, pos));
        surfaceScalarField rPsi_neg(interpolate(rPsi, neg));

        surfaceScalarField gammaStar_pos(interpolate(gammaStar, pos, DF_name));
        surfaceScalarField gammaStar_neg(interpolate(gammaStar, neg, DF_name));

        surfaceScalarField eStar_pos(interpolate(eStar, pos, DF_name));
        surfaceScalarField eStar_neg(interpolate(eStar, neg, DF_name));

        surfaceScalarField e_pos(interpolate(e, pos));
        surfaceScalarField e_neg(interpolate(e, neg));

        surfaceVectorField U_pos("U_pos", rhoU_pos / rho_pos);
        surfaceVectorField U_neg("U_neg", rhoU_neg / rho_neg);

        surfaceScalarField p_pos(interpolate(p, pos));
        surfaceScalarField p_neg(interpolate(p, neg));

        surfaceScalarField e_pos_pos("e_pos_pos", p_pos / (rho_pos * (gammaStar_pos - 1)) + eStar_pos); // = e_pos
        surfaceScalarField e_pos_neg("e_pos_neg", p_neg / (rho_neg * (gammaStar_pos - 1)) + eStar_pos);

        surfaceScalarField e_neg_pos("e_neg_pos", p_pos / (rho_pos * (gammaStar_neg - 1)) + eStar_neg);
        surfaceScalarField e_neg_neg("e_neg_neg", p_neg / (rho_neg * (gammaStar_neg - 1)) + eStar_neg); // = e_neg

        surfaceScalarField phiv_pos("phiv_pos", U_pos & mesh.Sf());
        surfaceScalarField phiv_neg("phiv_neg", U_neg & mesh.Sf());

        surfaceScalarField cSf_pos(
            "cSf_pos",
            interpolate(c, pos) * mesh.magSf());
        surfaceScalarField cSf_neg(
            "cSf_neg",
            interpolate(c, neg) * mesh.magSf());

        surfaceScalarField ap(
            "ap",
            max(max(phiv_pos + cSf_pos, phiv_neg + cSf_neg), v_zero));
        surfaceScalarField am(
            "am",
            min(min(phiv_pos - cSf_pos, phiv_neg - cSf_neg), v_zero));

        surfaceScalarField a_pos("a_pos", ap / (ap - am));

        surfaceScalarField amaxSf("amaxSf", max(mag(am), mag(ap)));

        surfaceScalarField aSf("aSf", am * a_pos);

        if (fluxScheme == "Tadmor")
        {
            aSf = -0.5 * amaxSf;
            a_pos = 0.5;
        }

        surfaceScalarField a_neg("a_neg", 1.0 - a_pos);

        phiv_pos *= a_pos;
        phiv_neg *= a_neg;

        surfaceScalarField aphiv_pos("aphiv_pos", phiv_pos - aSf);
        surfaceScalarField aphiv_neg("aphiv_neg", phiv_neg + aSf);

        // Reuse amaxSf for the maximum positive and negative fluxes
        // estimated by the central scheme
        amaxSf = max(mag(aphiv_pos), mag(aphiv_neg));

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
        //scalar t = runTime.value();

        volScalarField muEff("muEff", turbulence->muEff());
        volTensorField tauMC("tauMC", muEff * dev2(Foam::T(fvc::grad(U))));

        surfaceScalarField rhostar((ap * rho_neg - am * rho_pos - ((rhoU_neg & mesh.Sf()) - (rhoU_pos & mesh.Sf()))) / (ap - am));
        surfaceScalarField rhoQ = minMod(rho_neg - rhostar, rhostar - rho_pos) * ap * am / (ap - am);

        surfaceVectorField rhoUstar((ap * rhoU_neg - am * rhoU_pos - (rhoU_neg * (U_neg & mesh.Sf()) + p_neg * mesh.Sf() - rhoU_pos * (U_pos & mesh.Sf()) - p_pos * mesh.Sf())) / (ap - am));
        surfaceVectorField rhoUQ = minMod(rhoU_neg - rhoUstar, rhoUstar - rhoU_pos) * ap * am / (ap - am);

        phi = aphiv_pos * rho_pos + aphiv_neg * rho_neg - 0 * rhoQ;

        surfaceVectorField phiUp((aphiv_pos * rhoU_pos + aphiv_neg * rhoU_neg) + (a_pos * p_pos + a_neg * p_neg) * mesh.Sf() - 0 * rhoUQ);

        // --- Solve density
        solve(fvm::ddt(rho) + fvc::div(phi));

        // --- Solve momentum
        solve(fvm::ddt(rhoU) + fvc::div(phiUp));

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

        fvScalarMatrix rhoEEqn(fvm::ddt(rhoE));
        fvScalarMatrix rhoEEqn_tmp(fvm::ddt(rhoE_tmp));
        //fvScalarMatrix rhoeEqn(fvm::ddt(rhoe));

        /*         surfaceScalarField rhoEstar_pos((ap * rho_neg * (e_pos_neg + 0.5 * magSqr(U_neg)) - am * rho_pos * (e_pos_pos + 0.5 * magSqr(U_pos)) - ((U_neg & mesh.Sf()) * (rho_neg * (e_pos_neg + 0.5 * magSqr(U_neg)) + p_neg) - (U_pos & mesh.Sf()) * (rho_pos * (e_pos_pos + 0.5 * magSqr(U_pos)) + p_pos))) / (ap - am));
        surfaceScalarField rhoEQ_pos = minMod(rho_neg * (e_pos_neg + 0.5 * magSqr(U_neg)) - rhoEstar_pos, rhoEstar_pos - rho_pos * (e_pos_pos + 0.5 * magSqr(U_pos))) * ap * am / (ap - am);

        surfaceScalarField rhoEstar_neg((ap * rho_neg * (e_neg_neg + 0.5 * magSqr(U_neg)) - am * rho_pos * (e_neg_pos + 0.5 * magSqr(U_pos)) - ((U_neg & mesh.Sf()) * (rho_neg * (e_neg_neg + 0.5 * magSqr(U_neg)) + p_neg) - (U_pos & mesh.Sf()) * (rho_pos * (e_neg_pos + 0.5 * magSqr(U_pos)) + p_pos))) / (ap - am));
        surfaceScalarField rhoEQ_neg = minMod(rho_neg * (e_neg_neg + 0.5 * magSqr(U_neg)) - rhoEstar_neg, rhoEstar_neg - rho_pos * (e_neg_pos + 0.5 * magSqr(U_pos))) * ap * am / (ap - am);

        surfaceScalarField phiEp_pos(
            "phiEp_pos",
            aphiv_pos * (rho_pos * (e_pos_pos + 0.5 * magSqr(U_pos)) + p_pos) + aphiv_neg * (rho_neg * (e_pos_neg + 0.5 * magSqr(U_neg)) + p_neg) + aSf * p_pos - aSf * p_neg); //- 0 * rhoEQ_pos)

        surfaceScalarField phiEp_neg(
            "phiEp_pos",
            aphiv_pos * (rho_pos * (e_neg_pos + 0.5 * magSqr(U_pos)) + p_pos) + aphiv_neg * (rho_neg * (e_neg_neg + 0.5 * magSqr(U_neg)) + p_neg) + aSf * p_pos - aSf * p_neg); //- 0 * rhoEQ_neg)

        surfaceScalarField phiEp(
            "phiEp",
            aphiv_pos * (rho_pos * (e_pos + 0.5 * magSqr(U_pos)) + p_pos) + aphiv_neg * (rho_neg * (e_neg + 0.5 * magSqr(U_neg)) + p_neg) + aSf * p_pos - aSf * p_neg); */
        //surfaceScalarField phiep_pos("phiep_pos", aphiv_pos * (rho_pos * e_pos_pos) + aphiv_neg * (rho_neg * e_pos_neg));

        //surfaceScalarField phiep_neg("phiep_pos", aphiv_pos * (rho_pos * e_neg_pos) + aphiv_neg * (rho_neg * e_neg_neg));

        if (scheme == "doubleFlux")
        {
            surfaceScalarField phiEp_pos(
                "phiEp_pos",
                aphiv_pos * (rho_pos * (e_pos_pos + 0.5 * magSqr(U_pos)) + p_pos) + aphiv_neg * (rho_neg * (e_pos_neg + 0.5 * magSqr(U_neg)) + p_neg) + aSf * p_pos - aSf * p_neg); //- 0 * rhoEQ_pos)

            surfaceScalarField phiEp_neg(
                "phiEp_pos",
                aphiv_pos * (rho_pos * (e_neg_pos + 0.5 * magSqr(U_pos)) + p_pos) + aphiv_neg * (rho_neg * (e_neg_neg + 0.5 * magSqr(U_neg)) + p_neg) + aSf * p_pos - aSf * p_neg); //- 0 * rhoEQ_neg)

            rhoEEqn += fvc::div_doubleflux(phiEp_pos, phiEp_neg);
            //rhoeEqn += fvc::div_doubleflux(phiep_pos, phiep_neg);
        }
        else if (scheme == "doubleFlux++")
        {
            surfaceScalarField phiEp_pos(
                "phiEp_pos",
                aphiv_pos * (rho_pos * (e_pos_pos + 0.5 * magSqr(U_pos)) + p_pos) + aphiv_neg * (rho_neg * (e_pos_neg + 0.5 * magSqr(U_neg)) + p_neg) + aSf * p_pos - aSf * p_neg); //- 0 * rhoEQ_pos)

            surfaceScalarField phiEp_neg(
                "phiEp_pos",
                aphiv_pos * (rho_pos * (e_neg_pos + 0.5 * magSqr(U_pos)) + p_pos) + aphiv_neg * (rho_neg * (e_neg_neg + 0.5 * magSqr(U_neg)) + p_neg) + aSf * p_pos - aSf * p_neg); //- 0 * rhoEQ_neg)

            rhoEEqn_tmp += fvc::div_doubleflux(phiEp_pos, phiEp_neg);

            surfaceScalarField phiEp(
                "phiEp",
                aphiv_pos * (rho_pos * (e_pos + 0.5 * magSqr(U_pos)) + p_pos) + aphiv_neg * (rho_neg * (e_neg + 0.5 * magSqr(U_neg)) + p_neg) + aSf * p_pos - aSf * p_neg);

            rhoEEqn += fvc::div(phiEp);
        }
        else if (scheme == "conservativeFlux")
        {
            surfaceScalarField phiEp(
                "phiEp",
                aphiv_pos * (rho_pos * (e_pos + 0.5 * magSqr(U_pos)) + p_pos) + aphiv_neg * (rho_neg * (e_neg + 0.5 * magSqr(U_neg)) + p_neg) + aSf * p_pos - aSf * p_neg);

            rhoEEqn += fvc::div(phiEp);
        }
        else if (scheme == "mix")
        {
            surfaceScalarField phiEp_pos(
                "phiEp_pos",
                aphiv_pos * (rho_pos * (e_pos_pos + 0.5 * magSqr(U_pos)) + p_pos) + aphiv_neg * (rho_neg * (e_pos_neg + 0.5 * magSqr(U_neg)) + p_neg) + aSf * p_pos - aSf * p_neg); //- 0 * rhoEQ_pos)

            surfaceScalarField phiEp_neg(
                "phiEp_pos",
                aphiv_pos * (rho_pos * (e_neg_pos + 0.5 * magSqr(U_pos)) + p_pos) + aphiv_neg * (rho_neg * (e_neg_neg + 0.5 * magSqr(U_neg)) + p_neg) + aSf * p_pos - aSf * p_neg); //- 0 * rhoEQ_neg)
            surfaceScalarField phiEp(
                "phiEp",
                aphiv_pos * (rho_pos * (e_pos + 0.5 * magSqr(U_pos)) + p_pos) + aphiv_neg * (rho_neg * (e_neg + 0.5 * magSqr(U_neg)) + p_neg) + aSf * p_pos - aSf * p_neg);

            forAll(owner, facei)
            {
                if (FCcell[owner[facei]] == 1)
                {
                    phiEp_pos[facei] = phiEp[facei];
                    //phiEp_neg[facei] = phiEp[facei];
                }
                if (FCcell[neighbour[facei]] == 1)
                {
                    //phiEp_pos[facei] = phiEp[facei];
                    phiEp_neg[facei] = phiEp[facei];
                }
            }
            rhoEEqn += fvc::div_doubleflux(phiEp_pos, phiEp_neg);
        }
        if (!inviscid)
        {

            surfaceScalarField sigmaDotU(
                "sigmaDotU",
                (
                    fvc::interpolate(muEff) * mesh.magSf() * fvc::snGrad(U) + fvc::dotInterpolate(mesh.Sf(), tauMC)) &
                    (a_pos * U_pos + a_neg * U_neg));
            rhoEEqn -= fvc::div(sigmaDotU);
            rhoEEqn_tmp -= fvc::div(sigmaDotU);
        }

        solve(rhoEEqn);
        //solve(rhoeEqn);

        e = rhoE / rho - 0.5 * magSqr(U);
        //e = rhoe / rho;
        e.correctBoundaryConditions();

        if (scheme == "doubleFlux++")
        {
            solve(rhoEEqn_tmp);
            e_tmp = rhoE_tmp / rho - 0.5 * magSqr(U);
            e_tmp.correctBoundaryConditions();
        }

        if (!inviscid)
        {
            volScalarField &he = thermo.he();
            fvScalarMatrix EEqn(
                fvm::ddt(rho, e) - fvc::ddt(rho, e)
                //- fvm::laplacian(turbulence->alphaEff(), e)
                //- fvm::laplacian(kappa / Cp, he) //Todo kappa Cp
                //- fvm::laplacian(alpha, he) ==
                - fvc::laplacian(kappa, T) //==
                //-sumHeatDiffusion - sumHeatDiffusion2
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

            if (scheme == "doubleFlux++")
            {
                fvScalarMatrix EEqn_tmp(
                    fvm::ddt(rho, e) - fvc::ddt(rho, e) - fvc::laplacian(kappa, T));
                EEqn_tmp.relax();
                EEqn_tmp.solve();
            }
        }

        if (scheme == "doubleFlux" || scheme == "mix")
        {
            p.ref() = (e() - eStar()) * rho() * (gammaStar() - 1);
            p.max(1e3);
            //p.min(1e7);
        }
        if (scheme == "doubleFlux++")
        {
            p.ref() = (e_tmp() - eStar()) * rho() * (gammaStar() - 1);
        }

        //if (scheme == "doubleFlux++")
        //{
        //    rho_d = rho;
        //}

        thermo.correct();
        p.correctBoundaryConditions();
        T.correctBoundaryConditions();
        thermo.correct();

        if (scheme == "doubleFlux++")
        {
            rho_d -= rho;
        }

        rho.boundaryFieldRef() == psi.boundaryField() * p.boundaryField();
        U.correctBoundaryConditions();
        rho.correctBoundaryConditions();
        rhoU.boundaryFieldRef() == rho.boundaryField() * U.boundaryField();

        volScalarField Yt(0.0 * Y[0]);
        forAll(Y, i)
        {
            if (i != inertIndex && composition.active(i))
            {
                volScalarField &Yi = Y[i];
                volScalarField &rhoYi = rhoY[i];

                rhoYi.boundaryFieldRef() = rho.boundaryField() * Yi.boundaryField();

                Yi.max(1e-5);
                Yi.min(1.0);

                Yt += Yi;
            }
        }
        Y[inertIndex] = scalar(1) - Yt;
        Y[inertIndex].max(1e-5);

        //rhoE.boundaryFieldRef() ==
        //    rho.boundaryField() *
        //        (e.boundaryField() + 0.5 * magSqr(U.boundaryField()));

        rhoE = rho * (e + 0.5 * magSqr(U));
        //rhoe = rho * e;

        turbulence->correct();
        gammaStar = rho * c * c / p;
        eStar = e - p / (rho * (gammaStar - 1));
        Wmix = thermo.W();
	vol = Wmix/rho;
	e_mol = e*Wmix;
        //forAll(Dimix, i)
        //{
        //Dimix[i] = Dimix_t[i];
        // Dimix[i] = thermo.Dimix(i);
        //hei[i] = hei_t[i];
        // hei[i] = thermo.hei(i);
        //}
        surfalpha = frac * (1 - frac);

        Drho *= 0;
        forAll(owner, facei)
        {
            drho = mag(rho[owner[facei]] - rho[neighbour[facei]]);
            if (drho > rho[owner[facei]] * Drho[owner[facei]])
            {
                Drho[owner[facei]] = drho / rho[owner[facei]];
            }
            if (drho > rho[neighbour[facei]] * Drho[neighbour[facei]])
            {
                Drho[neighbour[facei]] = drho / rho[neighbour[facei]];
            }
        }

        forAll(FCcell, i)
        {
            if (Drho[i] > Drho_min)
            {
                FCcell[i] = 1;
            }
            else
            {
                FCcell[i] = 0;
            }
        }

        if (logflag == true)
        {
            cputime += clockTime_.timeIncrement();
            cputotal()
                << runTime.timeOutputValue()
                << ",    " << cputime << endl;
        }

        runTime.write();
        //rho_d.write();

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << endl;
    }

    Info << "End\n"
         << endl;

    return 0;
}

// ************************************************************************* //
