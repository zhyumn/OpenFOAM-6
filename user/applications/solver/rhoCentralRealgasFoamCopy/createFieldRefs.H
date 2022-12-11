//volScalarField& p = thermo.p();

const volScalarField &psi = thermo.psi();
const volScalarField &mu = thermo.mu();
const label inertIndex(composition.species()[inertSpecie]);

bool inviscid(true);

if (max(mu.primitiveField()) > 0.0)
{
    inviscid = false;
}
inviscid = false;//true;

volScalarField& frac = mesh.objectRegistry::lookupObjectRef<volScalarField>("vaporfrac");

volScalarField &alpha = mesh.objectRegistry::lookupObjectRef<volScalarField>("thermo:alpha");
volScalarField &rho_G = mesh.objectRegistry::lookupObjectRef<volScalarField>("thermo:rho_G");
rho_G.rename("rho_G");
rho_G.writeOpt() = IOobject::AUTO_WRITE;

//volScalarField &Z = mesh.objectRegistry::lookupObjectRef<volScalarField>("thermo:Z");

UPtrList<volScalarField> Y_G(Y.size());
forAll(Y_G, i)
{
    Y_G.set(i, &mesh.objectRegistry::lookupObjectRef<volScalarField>("Y_G:" + species[i]));
    Y_G[i].rename("Y_G_" + species[i]);
    Y_G[i].writeOpt()=IOobject::AUTO_WRITE;
}

//volScalarField& rho_thermo = mesh.objectRegistry::lookupObjectRef<volScalarField>("rho_thermo");
//Info << frac[0] << endl;

//FatalErrorInFunction << "break " << exit(FatalError);
