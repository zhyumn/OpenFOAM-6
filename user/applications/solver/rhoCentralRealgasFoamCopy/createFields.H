#include "createRDeltaT.H"

const word DF_name("DF");
const word Yi_name("Yi");
const word rhoYi_name("rhoYi");
Info << "Reading thermophysical properties\n"
     << endl;

autoPtr<rhoReactionThermo> pThermo(
    rhoReactionThermo::New(mesh));
rhoReactionThermo &thermo = pThermo();

volScalarField &c = mesh.objectRegistry::lookupObjectRef<volScalarField>("soundspeed");
volScalarField &p = thermo.p();
basicSpecieMixture &composition = thermo.composition();
PtrList<volScalarField> &Y = composition.Y();
volScalarField &T = thermo.T();
// initialize oldY for CoDAC
///PtrList<volScalarField>& oldY = composition.oldY();
//forAll(Y, i)
//{
//    oldY[i] = 1000.0;
//}

volScalarField &e = thermo.he();

const word inertSpecie(thermo.lookup("inertSpecie"));
// correction terms for diffusion velocity
volScalarField sumDiffusionCorrections1(
    IOobject(
        "gas_sumCorrections1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("zero", dimensionSet(1, -3, -1, 0, 0), 0.));

volScalarField sumDiffusionCorrections2(
    IOobject(
        "gas_sumCorrections2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("zero", dimensionSet(1, -3, -1, 0, 0), 0.));
//- Creates List of species name
speciesTable species = composition.species();
//const List<word>& speciesName = thermo.ListOfSpeciesName();

PtrList<volVectorField> YVi(Y.size());
forAll(YVi, i)
{
    YVi.set(
        i,
        new volVectorField(
            IOobject(
                "YVi_" + species[i],
                runTime.timeName(),
                mesh,
                IOobject::NO_READ, //IOobject::MUST_READ,
                IOobject::NO_WRITE),
            mesh,
            dimensionedVector("YVi_" + species[i], dimensionSet(0, 1, -1, 0, 0), vector(0., 0., 0.))));
}

//- Creates List of mass diffusivity field Dimix_i[m^2/s] of species ith in the mixture
UPtrList<volScalarField> Dimix_t(species.size());
//List<volScalarField*> Dimix_t(species.size());
forAll(species, i)
{
    Dimix_t.set(i, &mesh.objectRegistry::lookupObjectRef<volScalarField>("thermo:Dimix:" + species[i]));
    //Dimix_t[i] = &mesh.objectRegistry::lookupObjectRef<volScalarField>("thermo:Dimix:" + species[i]);
}
Info << "\nCreating field mass diffusivity Dimix_i" << endl;
PtrList<volScalarField> Dimix(Y.size());
forAll(Dimix, i)
{
    Dimix.set(
        i,
        new volScalarField(
            IOobject(
                "Dimix_" + species[i],
                runTime.timeName(),
                mesh,
                IOobject::NO_READ, //IOobject::MUST_READ,
                IOobject::NO_WRITE
                //IOobject::AUTO_WRITE
                ),
            Dimix_t[i]));
}

volScalarField &kappa = mesh.objectRegistry::lookupObjectRef<volScalarField>("thermo:kappa");
//- Create Thermal Conductivity field"kappa" [W/mK]
/*
Info << "\nCreating field thermal conductivity - kappa" << endl;
volScalarField kappa(
    IOobject(
        "kappa",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ, //IOobject::MUST_READ,
        IOobject::NO_WRITE),
    kappa_t);
*/

//- Creates Heat Capacity at constant pressure "Cp" [J/kgK]
/*Info << "\nCreating field heat capacity Cp" << endl;
volScalarField Cp
(
    IOobject
    (
        "Cp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,//IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.Cp()
);*/

//- Creates List of Enthalpy/Internal Energy field he_i[J/kg] of individual species
UPtrList<volScalarField> hei_t(species.size());
forAll(species, i)
{
    hei_t.set(i, &mesh.objectRegistry::lookupObjectRef<volScalarField>("hei:" + species[i]));
}
Info << "\nCreating field individual enthalpy he_i" << endl;
PtrList<volScalarField> hei(Y.size());
forAll(hei, i)
{
    hei.set(
        i,
        new volScalarField(
            IOobject(
                "hei_" + species[i],
                runTime.timeName(),
                mesh,
                IOobject::NO_READ, //IOobject::MUST_READ,
                IOobject::NO_WRITE),
            hei_t[i]));
}

Info
    << "Reading field U" << endl;
volVectorField U(
    IOobject(
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE),
    mesh);
volScalarField &rho = mesh.objectRegistry::lookupObjectRef<volScalarField>("thermo:rho");
rho.rename("rho");
rho.writeOpt() = IOobject::AUTO_WRITE;
/*volScalarField rho_write
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    rho
);*/
//rho=thermo.rho();
volVectorField rhoU(
    IOobject(
        "rhoU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    rho *U);

volScalarField rhoE(
    IOobject(
        "rhoE",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    rho *(e + 0.5 * magSqr(U)));
/* volScalarField rhoe(
    IOobject(
        "rhoe",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    rho *e ); */

PtrList<volScalarField> rhoY(Y.size());
forAll(Y, i)
{
    rhoY.set(
        i,
        new volScalarField(
            IOobject(
                "rhoYi",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE),
            rho * Y[i]));
}

surfaceScalarField pos(
    IOobject(
        "pos",
        runTime.timeName(),
        mesh),
    mesh,
    dimensionedScalar("pos", dimless, 1.0));

surfaceScalarField neg(
    IOobject(
        "neg",
        runTime.timeName(),
        mesh),
    mesh,
    dimensionedScalar("neg", dimless, -1.0));

surfaceScalarField phi("phi", fvc::flux(rhoU));

Info << "Creating turbulence model\n"
     << endl;
autoPtr<compressible::turbulenceModel> turbulence(
    compressible::turbulenceModel::New(
        rho,
        U,
        phi,
        thermo));

/*
Info << "\nCreating reaction model\n" << endl;
autoPtr<CombustionModel<rhoReactionThermo>> reaction
(
    CombustionModel<rhoReactionThermo>::New(thermo, turbulence())
);
*/
volScalarField Qdot(
    IOobject(
        "Qdot",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("Qdot", dimEnergy / dimVolume / dimTime, 0.0));

// volScalarField EqR(
//     IOobject(
//         "EqR",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::AUTO_WRITE),
//     mesh,
//     dimensionedScalar("EqR", dimless, 0.0));

volScalarField::Boundary& Qdot_boundary = Qdot.boundaryFieldRef();
//Nav added
PtrList<volScalarField> reactRate(Y.size());
forAll(Y, i)
{
    reactRate.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "reactRate_" + species[i],
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("reactRate" + species[i], dimMass / dimVolume / dimTime , 0.0)
        )
    );
}

volScalarField rhoGas
(
    IOobject
    (
        "rhoGas",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("rhoGas", dimMass / dimVolume, 0.0)
);

//- Creates Molecular weight field [kg/kmol]
Info << "Creating field molecular weight - W" << endl;
volScalarField Wmix(
    IOobject(
        "Wmix",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ, //IOobject::MUST_READ,
        //IOobject::AUTO_WRITE
        IOobject::NO_WRITE),
    thermo.W());
dimensionedScalar unitW("unit1", dimMass / dimMoles, 1);
PtrList<volScalarField> X(Y.size());
forAll(Y, i)
{
    X.set(
        i,
        new volScalarField(
            IOobject(
                "X." + species[i],
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE),
            Wmix * Y[i] * 1e3 / unitW / composition.W(i)));
}

// - for Energy Equation
volScalarField sumHeatDiffusion(
    IOobject(
        "sumHeatDiffusion",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("zero", dimensionSet(1, -1, -3, 0, 0), 0.));

volScalarField sumHeatDiffusion2(
    IOobject(
        "sumHeatDiffusion2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("zero", dimensionSet(1, -1, -3, 0, 0), 0.));

volScalarField gammaStar(
    IOobject(
        "gammaStar",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    rho *c *c / p);

volScalarField eStar(
    IOobject(
        "eStar",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    e - p / (rho * (gammaStar - 1)));

volScalarField rRhoStar(
    IOobject(
        "rStar",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    p / T);

volScalarField S(
    IOobject(
        "S",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    mesh, dimensionedScalar("zero", dimless, 0));
volScalarField Smin(
    IOobject(
        "Smin",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    mesh, dimensionedScalar("zero", dimless, 0));

surfaceScalarField troubled(
    IOobject(
        "troubled",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    mesh, dimensionedScalar("zero", dimless, 0));

volScalarField troubledVol(
    IOobject(
        "troubledVol",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE),
    mesh, dimensionedScalar("zero", dimless, 0));
dimensionedScalar unit1("unit1", dimPressure, 1);
dimensionedScalar unit2("unit2", dimensionSet(0, -3, 0, 0, 1, 0, 0), 1);
S = log(p / unit1) - gammaStar * log(constant::thermodynamic::RR * 1e3 / Wmix * rho / unit2);
forAll(Y, i)
{
    S += (gammaStar - 1) * X[i] * log(X[i]);
}

const labelUList &owner = mesh.owner();
const labelUList &neighbour = mesh.neighbour();
Smin = 1e50;

forAll(owner, facei)
{
    if (Smin[owner[facei]] > S[neighbour[facei]])
        Smin[owner[facei]] = S[neighbour[facei]];
    if (Smin[neighbour[facei]] > S[owner[facei]])
        Smin[neighbour[facei]] = S[owner[facei]];
}

volScalarField surfbeta(
    IOobject(
        "surfbeta",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE),
    mesh, dimensionedScalar("zero", dimless, 0));


//#include "createFvOptions.H"
