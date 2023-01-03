word fluxScheme("Kurganov");
if (mesh.schemesDict().readIfPresent("fluxScheme", fluxScheme))
{
    if ((fluxScheme == "Tadmor") || (fluxScheme == "Kurganov"))
    {
        Info<< "fluxScheme: " << fluxScheme << endl;
    }
    else
    {
        FatalErrorInFunction
            << "fluxScheme: " << fluxScheme
            << " is not a valid choice. "
            << "Options are: Tadmor, Kurganov"
            << abort(FatalError);
    }
}
/*
word divScheme("Doubleflux");
if (mesh.schemesDict().readIfPresent("divScheme", divScheme))
{
    if ((divScheme == "Conservativeflux") || (divScheme == "Doubleflux"))
    {
        Info<< "divScheme: " << divScheme << endl;
    }
    else
    {
        FatalErrorInFunction
            << "divScheme: " << fluxScheme
            << " is not a valid choice. "
            << "Options are: Conservativeflux, Doubleflux"
            << abort(FatalError);
    }
}
*/
