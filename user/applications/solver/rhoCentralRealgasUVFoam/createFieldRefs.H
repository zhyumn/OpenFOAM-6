//volScalarField& p = thermo.p();

const volScalarField &psi = thermo.psi();
const volScalarField &mu = thermo.mu();
const label inertIndex(composition.species()[inertSpecie]);
/*
bool inviscid(true);


inviscid true;

if (max(mu.primitiveField()) > 0.0)
{
    inviscid = false;
}
*/
bool inviscid(thermo.lookup("inviscid"));
//inviscid = true;

volScalarField &frac = mesh.objectRegistry::lookupObjectRef<volScalarField>("vaporfrac");

volScalarField &alpha = mesh.objectRegistry::lookupObjectRef<volScalarField>("thermo:alpha");
//volScalarField &rho_G = mesh.objectRegistry::lookupObjectRef<volScalarField>("thermo:rho_G");

//volScalarField &Z = mesh.objectRegistry::lookupObjectRef<volScalarField>("thermo:Z");
/*
UPtrList<volScalarField> Y_G(Y.size());
forAll(Y_G, i)
{
    Y_G.set(i, &mesh.objectRegistry::lookupObjectRef<volScalarField>("Y_G:" + species[i]));
}
*/
//volScalarField& rho_thermo = mesh.objectRegistry::lookupObjectRef<volScalarField>("rho_thermo");
//Info << frac[0] << endl;
bool logflag = thermo.lookupOrDefault<bool>("ISATlog", false);
word scheme(thermo.lookup("scheme"));
Info << "scheme is " << scheme << endl;
autoPtr<OFstream> cputotal;
if (logflag == true)
{
    cputotal.set(new OFstream(mesh.time().path() / "ISAT_VLE" / "cputotal"));
}
//autoPtr<OFstream> cputotal(new OFstream(mesh.time().path() / "ISAT_VLE" / "cputotal"));
const clockTime clockTime_;

//FatalErrorInFunction << "break " << exit(FatalError);
