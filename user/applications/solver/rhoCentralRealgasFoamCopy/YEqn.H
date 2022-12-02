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
    forAll(T,celli){
        double T_i = T[celli];
        double P_i = p[celli];
        scalar W_r[6];

        W_r[0] = 0.170338;
        W_r[1] = 0.032;
        W_r[2] = 0.04401;
        W_r[3] =  0.0280101;
        W_r[4] = 0.018;
        W_r[5] = 0.0280134;
        scalar Mavg = 0.0;
        scalar Mavg_G = 0.0;

        forAll(Y,i){
            Mavg += Y[i][celli]/W_r[i];
            Mavg_G += Y_G[i][celli]/W_r[i];
        }
        Mavg = 1.0/Mavg;
        Mavg_G = 1.0/Mavg_G;
        double mass_frac = Mavg_G*frac[celli]/Mavg;

        double rho_i = mass_frac*rho[celli]; //thermo.gasPhaseRho(T_i, P_i);
        rhoGas[celli] = rho_i;

        thermo.customChemistry(T_i,P_i,rho_i);

        forAll(Y,i){ //(int i=0;i<5;i++){
            double Ygas_i = Y_G[i][celli]; //thermo.gasPhaseY(i,celli);
            thermo.setY(Ygas_i,i);
        }

        //EqR = (Y_G[0][celli]/Y_G[1][celli])/((1*W_r[0])/(18.5*W_r[1]));

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

            scalar W_r[6];

            W_r[0] = 0.170338;
            W_r[1] = 0.032;
            W_r[2] = 0.04401;
            W_r[3] =  0.0280101;
            W_r[4] = 0.018;
            W_r[5] = 0.0280134;
            scalar Mavg = 0.0;
            scalar Mavg_G = 0.0;

            forAll(Y,i){
                Mavg += Y[i].boundaryField()[patchi][facei]/W_r[i];
                Mavg_G += Y_G[i].boundaryField()[patchi][facei]/W_r[i];
            }
            Mavg = 1.0/Mavg;
            Mavg_G = 1.0/Mavg_G;
            double mass_frac = Mavg_G*frac.boundaryField()[patchi][facei]/Mavg;

            double rho_i = mass_frac*rho.boundaryField()[patchi][facei]; //thermo.gasPhaseRho(T_i, P_i);

            thermo.customChemistry(T_i,P_i,rho_i);

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

            surfaceScalarField rhoYstar((ap * rhoYi_neg - am * rhoYi_pos - ((U_neg * rhoYi_neg & mesh.Sf()) - (U_pos * rhoYi_pos & mesh.Sf()))) / (ap - am));
            surfaceScalarField rhoYQ = minMod(rhoYi_neg - rhoYstar, rhoYstar - rhoYi_pos) * ap * am / (ap - am);

            //phiY.set(i, aphiv_pos * rhoYi_pos + aphiv_neg * rhoYi_neg - rhoYQ);

            surfaceScalarField phiYi("phiYi", aphiv_pos * rhoYi_pos + aphiv_neg * rhoYi_neg - 0*rhoYQ);

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
                    reactRate[i]+fvm::laplacian(rho * Dimix[i], Yi) + fvc::laplacian(rho * Dimix[i] * Yi / Wmix, Wmix) - sumDiffusionCorrections1 - sumDiffusionCorrections2
                    //+ fvOptions(rho, Yi)
                );

                YiEqn.relax();

                //fvOptions.constrain(YiEqn);

                YiEqn.solve(mesh.solver("Yi"));

                //fvOptions.correct(Yi);

                Yi.correctBoundaryConditions();
            }
            Yi.max(1e-5);
            Yi.min(1.0);

            rhoYi = rho * Yi;
            Yt += Yi;
        }
    }
    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].max(1e-5);
    // update oldY here for CoDAC
    //oldY = tempY;
}
