/*
tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi_h)")
    )
);
*/

PtrList<volScalarField> tempY(Y.size());

{
    forAll(Y, i)
    {
        tempY.set(
            i,
            new volScalarField(
                IOobject(
                    "tempY",
                    mesh.time().timeName(),
                    mesh),
                mesh,
                scalar(0.0)));
    }

    forAll(YVi, i)
    {
        YVi[i] *= 0.;
    }
    volScalarField Yt(0.0 * Y[0]);
    // add turbulence effect on D, assuming constant turbulence Lewis number
    // turbulence thermal diffusivity alphat
    volScalarField alphat = turbulence->alphat();
    const scalar Let(thermo.lookupOrDefault<scalar>("Let", 1.0));
    volScalarField Dt = 0 * alphat / (rho * Let);
    forAll(Y, i)
    {
        Dimix[i] = Dimix[i] + Dt;
    }
    //- Reset diffusion velocity to zero for each time step

    //-Calculate diffusion velocity of inert species
    volVectorField YVt(0.0 * YVi[0]);
    forAll(YVi, i)
    {
        //- Update diffusion velocity for all species
        YVi[i] = (-(Dimix[i] / Wmix) * Wmix * fvc::grad(Y[i]) - (Dimix[i] / Wmix) * Y[i] * fvc::grad(Wmix));
        if (i != inertIndex && composition.active(i))
        {
            YVt += YVi[i];
        }
    }

    YVi[inertIndex] = -YVt;

    //reaction->correct();
    //Qdot = reaction->Qdot();
    //volScalarField Yt(0.0 * Y[0]);
    tempY = Y;
    forAll(Y, i)
    {

        // calculate correction term
        sumDiffusionCorrections1 *= 0.;
        sumDiffusionCorrections2 *= 0.;
        forAll(Y, k)
        {
            sumDiffusionCorrections1 += fvc::laplacian(Y[i] * rho * Dimix[k], Y[k]);
            sumDiffusionCorrections2 += fvc::laplacian(Y[i] * rho * Dimix[k] * Y[k] / Wmix, Wmix);
        }

        if (i != inertIndex && composition.active(i))
        {
            volScalarField &Yi = Y[i];
            volScalarField &rhoYi = rhoY[i];

            // Central upwind scheme
            // interpolation of Yi to cell faces
            surfaceScalarField rhoYi_pos(interpolate(rhoYi, pos, rhoYi_name));
            surfaceScalarField rhoYi_neg(interpolate(rhoYi, neg, rhoYi_name));

            surfaceScalarField phiYi("phiYi", aphiv_pos * rhoYi_pos + aphiv_neg * rhoYi_neg);

            solve(fvm::ddt(rhoYi) + fvc::div(phiYi));

            Yi = rhoYi / rho;
            Yi.correctBoundaryConditions();
            Yi.max(0.0);

            rhoYi.boundaryFieldRef() = rho.boundaryField() * Yi.boundaryField();
            if (!inviscid)
            {
                fvScalarMatrix YiEqn(
                    fvm::ddt(rho, Yi) - fvc::ddt(rho, Yi) + fvm::div(phi, Yi) - fvc::div(phi, Yi)
                    //- fvm::laplacian(rho * D[i], Yi)
                    //+ fvc::div(rho*YVi[i])

                    // - fvm::laplacian(turbulence->muEff(), Yi)
                    ==
                    //reaction->R(Yi)
                    fvm::laplacian(rho * Dimix[i], Yi) + fvc::laplacian(rho * Dimix[i] * Yi / Wmix, Wmix) - sumDiffusionCorrections1 - sumDiffusionCorrections2
                    //+ fvOptions(rho, Yi)
                );

                YiEqn.relax();

                //fvOptions.constrain(YiEqn);

                YiEqn.solve(mesh.solver("Yi"));

                //fvOptions.correct(Yi);

                Yi.correctBoundaryConditions();
            }
            Yi.max(0.0);
            Yi.min(1.0);

            rhoYi = rho * Yi;
            Yt += Yi;
        }
    }
    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].max(0);
    // update oldY here for CoDAC
    //oldY = tempY;
}
