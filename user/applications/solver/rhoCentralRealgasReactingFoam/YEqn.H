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
    /*Info << "\n93:" << YVi[0][93]
        << "\n94:" << YVi[0][94]
        << "\n95:" << YVi[0][95]
        << "\n96:" << YVi[0][96]
        << "\n97:" << YVi[0][97] << endl;*/
    YVi[inertIndex] = -YVt;

    //reaction->correct();
    //Qdot = reaction->Qdot();

    // thermo.customChemistry(T,p,rho);
    // thermo.setY(Y);
    // forAll(Y,i)
    // {
    //     reactRate[i] = thermo.RR(i);
    // }
    // Qdot = thermo.Qdot();

    forAll(T,celli){
        double T_i = T[celli];
        double P_i = p[celli];
        double rho_i = rho_G[celli]; //thermo.gasPhaseRho(T_i, P_i);
        
        thermo.customChemistry(T_i,P_i,rho_i);

        forAll(Y,i){ //(int i=0;i<5;i++){
            double Ygas_i = Y_G[i][celli]; //thermo.gasPhaseY(i,celli);
            thermo.setY(Ygas_i,i);
        }

        thermo.setRate();

        forAll(Y,i){
            reactRate[i][celli] = thermo.RR(i);
        }

        Qdot[celli] = thermo.Qdot();
    }

    forAll(Qdot_boundary,patchi){
        forAll(Qdot_boundary[patchi],facei){

            double T_i = T.boundaryField()[patchi][facei];
            double P_i = p.boundaryField()[patchi][facei];
            double rho_i = rho.boundaryField()[patchi][facei]; //thermo.gasPhaseRho(T_i, P_i);
            
            thermo.customChemistry(T_i,P_i,rho_i);

            for(int i=0;i<6;i++){
                double Ygas_i = Y[i].boundaryField()[patchi][facei]; //thermo.gasPhaseY(i,celli);
                thermo.setY(Ygas_i,i);
            }

            thermo.setRate();
            
            for(int i=0;i<6;i++){
                reactRate[i].boundaryFieldRef()[patchi][facei] = thermo.RR(i);
            }

            Qdot.boundaryFieldRef()[patchi][facei] = thermo.Qdot();
        }
    }
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
            surfaceScalarField Yi_pos(interpolate(Yi, pos, Yi_name));
            surfaceScalarField Yi_neg(interpolate(Yi, neg, Yi_name));
            surfaceScalarField phiYi(
                "phiYi",
                aphiv_pos * rho_pos * Yi_pos + aphiv_neg * rho_neg * Yi_neg);

            surfaceScalarField phiYi_pos(
                "phiYi_pos",
                aphiv_pos_pos * rho_pos * Yi_pos + aphiv_pos_neg * rho_neg * Yi_neg);

            surfaceScalarField phiYi_neg(
                "phiYi_neg",
                aphiv_neg_pos * rho_pos * Yi_pos + aphiv_neg_neg * rho_neg * Yi_neg);

            if (divScheme == "Doubleflux")
            {
                solve(fvm::ddt(rhoYi) + fvc::div_doubleflux(phiYi_pos, phiYi_neg) == reactRate[i]);
            }
            else if (divScheme == "Conservativeflux")
            {

                solve(fvm::ddt(rhoYi) + fvc::div(phiYi) == reactRate[i]);
            }

            Yi = rhoYi / rho;
            Yi.correctBoundaryConditions();
            Yi.max(0.000001);

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
                    reactRate[i] + fvm::laplacian(rho * Dimix[i], Yi) + fvc::laplacian(rho * Dimix[i] * Yi / Wmix, Wmix) - sumDiffusionCorrections1 - sumDiffusionCorrections2
                    //+ fvOptions(rho, Yi)
                );

                YiEqn.relax();

                //fvOptions.constrain(YiEqn);

                YiEqn.solve(mesh.solver("Yi"));

                //fvOptions.correct(Yi);

                Yi.correctBoundaryConditions();
            }
            Yi.max(0.000001);
            Yi.min(1.0);

            rhoYi = rho * Yi;
            Yt += Yi;
        }
    }
    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].max(0.000001);
    // update oldY here for CoDAC
    //oldY = tempY;
}
