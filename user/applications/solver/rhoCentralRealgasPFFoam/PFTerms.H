volScalarField rho_diff = rho_l - rho_v;
rho_diff.max(1e-20);
volScalarField phi_l = 1-frac; //liq frac
volScalarField psi_l = frac;

// forAll(T,i){
//     Mw1[i] = ( composition.W(0)*composition.W(1) )/ (composition.W(0)*Y_v[1][i] + composition.W(1)*Y_v[0][i]);
//     phi_l[i] = 1 - (sqr(Mw1[i]*1E-3)*frac[i]/(rho_v[i]+1e-8))/(sqr(Wmix[i])/rho[i]);
//     //Info<<psi_l[i]<<endl;
// }
// phi_l.max(0);

dimensionedScalar minVol = min(mesh.V());
dimensionedScalar minLen = cbrt(minVol);

scalar interface_thickness = int_points*(0.002/512.0);

volScalarField eps(
IOobject
(
"delta",
runTime.timeName(),
runTime,
IOobject::NO_READ,
IOobject::NO_WRITE
),
mesh,
dimensionedScalar("cellDim", dimLength, interface_thickness) //int_points*(minLen) //dimensionedScalar("cellDim", dimLength, 0.5e-7)
);

volScalarField magU = mag(U);

dimensionedScalar lambda
(
"lambda",
dimensionSet(0,1,-1,0,0,0,0),
gMax(magU)
);

dimensionedScalar lenDim
(
"lenDim",
dimensionSet(0,1,0,0,0,0,0),
1
);

volVectorField gradphi_l = fvc::grad(phi_l);
volScalarField alpha = phi_l * (1.0 - phi_l) ;
volScalarField mag_gradphi_l(mag(gradphi_l));
mag_gradphi_l.max(1e-20);
Info<<"interface thickness = "<<interface_thickness<<endl;
volVectorField al = lambda *( eps - alpha /(mag_gradphi_l) ) * gradphi_l; //lambda*( eps - (alpha/(mag_gradphi_l))) * gradphi_l; lambda * ( eps
volScalarField mag_al = mag(al);
Info<<gMax(mag_al)<<endl;
//thermo calculations
volVectorField f = (rho_l-rho_v)*al;


// //Regulation terms

// mass
volScalarField massReg = fvc::div(f);

// momentum
surfaceVectorField f_pos(interpolate(f, pos, U.name()));
surfaceVectorField f_neg(interpolate(f, neg, U.name()));
surfaceVectorField fUp((aphiv_pos * f_pos + aphiv_neg * f_neg));
volVectorField momReg = fvc::div(fUp);

// species
PtrList<volScalarField> div_Rs(Y.size()); //div(rhoYi*al)
forAll(Y, i)
{
    div_Rs.set(
        i,
        new volScalarField(
            IOobject(
                "div_Rs"+species[i],
                mesh.time().timeName(),
                mesh),
            mesh,
            dimensionedScalar("div_Rs"+species[i], dimDensity/dimTime, 0.0))
            );
}

forAll(Y,i){
    div_Rs[i] = fvc::div( (rho_l*Y_l[i] - rho_v*Y_v[i])*al );
    div_Rs[i].rename("Yreg_" + species[i]);
    div_Rs[i].writeOpt() = IOobject::AUTO_WRITE;
}

// energy
volScalarField energyReg = fvc::div( (e + magSqr(U)/2)*f + (rho_l*h_l - rho_v*h_v)*al);  //magSqr(U)/2)*f )+

// surfaceScalarField energyRegp(
//                 "energyRegp",
//                 aphiv_pos * (rho_pos * (0.5 * magSqr(U_pos)) ) + aphiv_neg * (rho_neg * (0.5 * magSqr(U_neg)) ) );


massReg.rename("massReg");
phi_l.rename("phi_l");
momReg.rename("momReg");
al.rename("al_Reg");
energyReg.rename("energyReg");
massReg.writeOpt() = IOobject::AUTO_WRITE;
phi_l.writeOpt() = IOobject::AUTO_WRITE;
momReg.writeOpt() = IOobject::AUTO_WRITE;
al.writeOpt() = IOobject::AUTO_WRITE;
energyReg.writeOpt() = IOobject::AUTO_WRITE;


