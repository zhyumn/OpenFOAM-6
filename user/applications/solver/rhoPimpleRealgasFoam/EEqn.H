{
    Info<<"1"<<endl;
    volScalarField& he = thermo.he();
Info<<"2"<<endl;
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
      - fvm::laplacian(turbulence->alphaEff(), he)
     ==
        Qdot
      + fvOptions(rho, he)
    );
Info<<"3"<<endl;
    EEqn.relax();

    fvOptions.constrain(EEqn);
Info<<"4"<<endl;
    EEqn.solve();
Info<<"5"<<endl;
    fvOptions.correct(he);
Info<<"6"<<endl;
    thermo.correct();
Info<<"7"<<endl;
    Info<< "min/max(T) = "
        << min(T).value() << ", " << max(T).value() << endl;
}
