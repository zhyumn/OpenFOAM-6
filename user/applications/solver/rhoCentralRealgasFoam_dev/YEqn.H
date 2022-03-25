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



/*
PtrList<volVectorField> Dm =
dynamic_cast<const reactingMixture<gasHThermoPhysics>&>
(composition).Dm(p, T, Y, rho);

PtrList<volScalarField> hsi(Y.size());
PtrList<surfaceScalarField> J(Y.size());
PtrList<volScalarField> D(Y.size());
PtrList<volScalarField> tempY(Y.size());

forAll(Y, i)
{
    hsi.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "hsi",
                mesh.time().timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar
            (
                "hsi",
                dimEnergy / dimMass,
                Zero
            )
        )
    );

    tempY.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "tempY",
                mesh.time().timeName(),
                mesh
            ),
            mesh,
            scalar(0.0)
        )
    );

    D.set
    (
        i,
        mag(Dm[i]) / (dimensionedScalar("SMALL", dimensionSet(0, -1, 0, 0, 0), Foam::SMALL)
            + mag(fvc::grad(Y[i]).ref()))
    );

    J.set
    (
        i,
        fvc::interpolate(Dm[i] * rho) & mesh.Sf()
    );

}

forAll(Y, i)
{
    volScalarField& tHsi = hsi[i];
    forAll(tHsi, celli)
    {
        tHsi[celli] = composition.Hs(i, p[celli], T[celli]);
    }
    volScalarField::Boundary& Bf = tHsi.boundaryFieldRef();
    forAll(Bf, patchi)
    {
        forAll(Bf[patchi], facei)
        {
            Bf[patchi][facei] =
                composition.Hs
                (
                    i,
                    p.boundaryField()[patchi][facei],
                    T.boundaryField()[patchi][facei]
                );
        }
    }
}
*/


PtrList<volScalarField> tempY(Y.size());

{
    forAll(Y, i)
    {
        tempY.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "tempY",
                    mesh.time().timeName(),
                    mesh
                ),
                mesh,
                scalar(0.0)
            )
        );
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
    volScalarField Dt = alphat / (rho * Let);
    forAll(Y, i)
    {
        Dimix[i] = Dimix[i] + Dt;
    }
    //- Reset diffusion velocity to zero for each time step

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
            volScalarField& Yi = Y[i];
            volScalarField& rhoYi = rhoY[i];
            // Central upwind scheme
                // interpolation of Yi to cell faces
            surfaceScalarField Yi_pos = interpolate(Yi, pos, T.name());
            surfaceScalarField Yi_neg = interpolate(Yi, neg, T.name());
            surfaceScalarField phiYi
            (
                "phiYi",
                aphiv_pos * rho_pos * Yi_pos + aphiv_neg * rho_neg * Yi_neg
            );

            solve(fvm::ddt(rhoYi) + fvc::div(phiYi));

            Yi = rhoYi / rho;
            Yi.correctBoundaryConditions();
            Yi.max(0.0);

            rhoYi.boundaryFieldRef() = rho.boundaryField() * Yi.boundaryField();

            fvScalarMatrix YiEqn
            (
                fvm::ddt(rho, Yi) - fvc::ddt(rho, Yi)
                + fvm::div(phi, Yi) - fvc::div(phi, Yi)
                //- fvm::laplacian(rho * D[i], Yi)
                - fvm::laplacian(rho * Dimix[i], Yi)
                // - fvm::laplacian(turbulence->muEff(), Yi)
                ==
                //reaction->R(Yi)
                fvc::laplacian(rho * Dimix[i] * Yi / Wmix, Wmix)
                - sumDiffusionCorrections1
                - sumDiffusionCorrections2
                //+ fvOptions(rho, Yi)
            );

            YiEqn.relax();

            //fvOptions.constrain(YiEqn);

            YiEqn.solve(mesh.solver("Yi"));

            //fvOptions.correct(Yi);

            Yi.correctBoundaryConditions();
            Yi.max(0.0);
            Yi.min(1.0);

            rhoYi = rho * Yi;
            Yt += Yi;
        }
    }
    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].max(0.0);
    // update oldY here for CoDAC
    //oldY = tempY;

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
}
