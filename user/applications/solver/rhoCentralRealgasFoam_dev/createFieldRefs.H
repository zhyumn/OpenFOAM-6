volScalarField& p = thermo.p();
const volScalarField& T = thermo.T();
const volScalarField& psi = thermo.psi();
const volScalarField& mu = thermo.mu();
const label inertIndex(composition.species()[inertSpecie]);

bool inviscid(true);

if (max(mu.primitiveField()) > 0.0)
{
    inviscid = false;
}



volScalarField& frac = mesh.objectRegistry::lookupObjectRef<volScalarField>("vaporfrac");
volScalarField& c = mesh.objectRegistry::lookupObjectRef<volScalarField>("soundspeed");
volScalarField& alpha = mesh.objectRegistry::lookupObjectRef<volScalarField>("thermo:alpha");
//volScalarField& rho_thermo = mesh.objectRegistry::lookupObjectRef<volScalarField>("rho_thermo");
//Info << frac[0] << endl;

//FatalErrorInFunction << "break " << exit(FatalError);