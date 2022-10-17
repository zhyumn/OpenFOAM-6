#include "createRDeltaT.H"

const word DF_name("DF");
const word Yi_name("Yi");
Info << "Reading thermophysical properties\n"
     << endl;

autoPtr<rhoReactionThermo> pThermo(
    rhoReactionThermo::New(mesh));
rhoReactionThermo &thermo = pThermo();

volScalarField &c = mesh.objectRegistry::lookupObjectRef<volScalarField>("soundspeed");
volScalarField &p = thermo.p();
basicSpecieMixture &composition = thermo.composition();
PtrList<volScalarField> &Y = composition.Y();

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

volScalarField &kappa_t = mesh.objectRegistry::lookupObjectRef<volScalarField>("thermo:kappa");
//- Create Thermal Conductivity field"kappa" [W/mK]
Info << "\nCreating field thermal conductivity - kappa" << endl;
volScalarField kappa(
    IOobject(
        "kappa",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ, //IOobject::MUST_READ,
        IOobject::NO_WRITE),
    kappa_t);

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
        IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("Qdot", dimEnergy / dimVolume / dimTime, 0.0));

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

//#include "createFvOptions.H"