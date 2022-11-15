{

    volScalarField& he = thermo.he();

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + mvConvection->fvmDiv(phi, he)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + (
            he.name() == "e"
          ? fvc::div
            (
                fvc::absolute(phi/fvc::interpolate(rho), U),
                p,
                "div(phiv,p)"
            )
          : -dpdt
        )
     // - fvm::laplacian(turbulence->alphaEff(), he)
     ==
        Qdot
      //+ fvOptions(rho, he)
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    volScalarField Gp=mag(fvc::grad(p)/(pu*2e3));
    Gp.min(1);
    thermo.correct();
    thermo.correct();
    psi= 1/((he-eStar)*(gammaStar-1))*(1-Gp) + thermo.psi()*Gp;
    rho=p*psi;

    Info<< "min/max(T) = "
        << min(T).value() << ", " << max(T).value() << endl;
}
