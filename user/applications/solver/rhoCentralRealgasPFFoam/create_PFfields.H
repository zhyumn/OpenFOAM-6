// PF term fields

scalar int_points = readScalar(thermo.lookup("interfacePoints"));

volScalarField &rho_l = mesh.objectRegistry::lookupObjectRef<volScalarField>("thermo:rho_l");
volScalarField &rho_v = mesh.objectRegistry::lookupObjectRef<volScalarField>("thermo:rho_v");
volScalarField &h_l = mesh.objectRegistry::lookupObjectRef<volScalarField>("thermo:h_l");
volScalarField &h_v = mesh.objectRegistry::lookupObjectRef<volScalarField>("thermo:h_v");

UPtrList<volScalarField> Y_l(Y.size());
UPtrList<volScalarField> Y_v(Y.size());
forAll(Y_l, i)
{
    Y_l.set(i, &mesh.objectRegistry::lookupObjectRef<volScalarField>("Y_l:" + species[i]));
    Y_v.set(i, &mesh.objectRegistry::lookupObjectRef<volScalarField>("Y_v:" + species[i]));
}

volScalarField Mw1(
    IOobject(
        "Mw1",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("Mw1", dimMass/dimMoles, 0.0));

volScalarField Mw2(
    IOobject(
        "Mw2",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("Mw2", dimMass/dimMoles, 0.0));

// volScalarField rho_l
// (
//     IOobject
//     (
//         "rho_l",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::AUTO_WRITE
//     ),
//     mesh,
//     dimensionedScalar("rho_l", dimDensity, 0.0)
// );

// volScalarField rho_v
// (
//     IOobject
//     (
//         "rho_v",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::AUTO_WRITE
//     ),
//     mesh,
//     dimensionedScalar("rho_v", dimDensity, 0.0)
// );

// volScalarField h_l
// (
//     IOobject
//     (
//         "h_l",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::AUTO_WRITE
//     ),
//     mesh,
//     dimensionedScalar("h_l", dimEnergy/dimMass, 0.0)
// );

// volScalarField h_v
// (
//     IOobject
//     (
//         "h_v",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::AUTO_WRITE
//     ),
//     mesh,
//     dimensionedScalar("h_v", dimDensity/dimMass, 0.0)
// );

// PtrList<volScalarField> Y_l(Y.size());
// forAll(Y, i)
// {
//     Y_l.set(
//         i,
//         new volScalarField(
//             IOobject(
//                 "Yi_l"+ species[i],
//                 runTime.timeName(),
//                 mesh,
//                 IOobject::NO_READ,
//                 IOobject::AUTO_WRITE),
//             Y[i]));
// }

// PtrList<volScalarField> Y_v(Y.size());
// forAll(Y, i)
// {
//     Y_v.set(
//         i,
//         new volScalarField(
//             IOobject(
//                 "Yi_v"+ species[i],
//                 runTime.timeName(),
//                 mesh,
//                 IOobject::NO_READ,
//                 IOobject::AUTO_WRITE),
//             Y[i]));
// }

// volScalarField eps
// (
//     IOobject
//     (
//     "delta",
//     runTime.timeName(),
//     runTime,
//     IOobject::NO_READ,
//     IOobject::NO_WRITE
//     ),
//     mesh,
//     dimensionedScalar("cellDim", dimLength, 0)
// );
// scalarField cbrtvol = cbrt(mesh.V());
// eps.primitiveFieldRef() = cbrtvol;

