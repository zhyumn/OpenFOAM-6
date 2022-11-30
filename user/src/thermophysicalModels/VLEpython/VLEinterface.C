#include "VLEinterface.h"
template <>
bool VLE<chungTransportMixture<PengRobinsonMixture<multispecie<Stype>>>>::noVLE = false;
/*
void solver::update()
{

    dict = new Foam::dictionary(Foam::IFstream(path.c_str())());
    int n = specie.size();

    char** slist = new char* [n];
    for (int i = 0;i < n;i++)
    {
        slist[i] = new char[20];
        strcpy(slist[i], specie[i].c_str());
        //std::cout << slist[i] << std::endl;
    }
    s = new Foam::speciesTable(n, (const char**)slist);
    PR = new Foam::PengRobinsonM<Foam::specie>(*s, *dict);
    comp_liq.resize(n);
    comp_gas.resize(n);
    equalconstant.resize(n);
    delete slist;

}
solver::solver(std::string path_i) :path(path_i), dict(nullptr), PR(nullptr) {}
solver::~solver() {
    if (dict)
        delete dict;
    if (s)
        delete s;
    if (PR)
        delete PR;
}
void solver::solve(bool flag)
{
    vaporfra = 0.9;
    Foam::scalarList equalconstant_of(X.size(), Foam::Zero);
    Foam::scalarList comp_of(X.size(), Foam::Zero);
    Foam::scalarList comp_liq_of(X.size(), Foam::Zero);
    Foam::scalarList comp_gas_of(X.size(), Foam::Zero);
    for (int i = 0;i < X.size();i++)
    {
        comp_of[i] = X[i];
        equalconstant_of[i] = equalconstant[i];
    }
    PR->TPn_flash(P_, T_, comp_of, comp_liq_of, comp_gas_of, vaporfra, equalconstant_of, flag);
    for (int i = 0;i < X.size();i++)
    {
        comp_liq[i] = comp_liq_of[i];
        comp_gas[i] = comp_gas_of[i];
        equalconstant[i] = equalconstant_of[i];
    }

}
double solver::Z()
{
    Foam::scalarList comp_liq_of(X.size(), Foam::Zero);
    Foam::scalarList comp_gas_of(X.size(), Foam::Zero);
    Foam::scalarList comp_of(X.size(), Foam::Zero);
    for (int i = 0;i < X.size();i++)
    {
        comp_liq_of[i] = comp_liq[i];
        comp_gas_of[i] = comp_gas[i];
        comp_of[i] = X[i];
    }
    double alphagas = PR->Evaluate_alpha(P_, T_, vaporfra, comp_liq_of, comp_gas_of, comp_of);
    double VspecificGas = PR->volmmix_phase(0, P_, T_, comp_gas_of); //m3/mol
    double VspecificLiq = PR->volmmix_phase(1, P_, T_, comp_liq_of);
    double ZGas = P * VspecificGas / (RR * 1.0e-03 * T);
    double ZLiq = P * VspecificLiq / (RR * 1.0e-03 * T);
    double ZMixture = ZGas * alphagas + ZLiq * (1.0 - alphagas);
    return  ZMixture;

}
double solver::density()
{
    Foam::scalarList comp_liq_of(X.size(), Foam::Zero);
    Foam::scalarList comp_gas_of(X.size(), Foam::Zero);
    Foam::scalarList comp_of(X.size(), Foam::Zero);
    for (int i = 0;i < X.size();i++)
    {
        comp_liq_of[i] = comp_liq[i];
        comp_gas_of[i] = comp_gas[i];
        comp_of[i] = X[i];
    }

    double mw_gas = PR->mwmix(comp_gas_of); //kg/mol
    double mw_liq = PR->mwmix(comp_liq_of); //kg/mol
    double mw_mixture = mw_gas * vaporfra + mw_liq * (1.0 - vaporfra);
    double rhoMixture = P * mw_mixture / (Z() * RR * 1.0e-03 * T);
    return rhoMixture;

}
bool solver::twophase(double rho, double lt, double rt)
{
    double tm = (lt + rt) / 2;

    Foam::scalarList comp_liq_of(X.size(), Foam::Zero);
    Foam::scalarList comp_gas_of(X.size(), Foam::Zero);
    Foam::scalarList comp_of(X.size(), Foam::Zero);
    Foam::scalarList equalconstant_of(X.size(), Foam::Zero);
    for (int i = 0;i < X.size();i++)
    {
        comp_of[i] = X[i];
    }
    while (rt - lt > 1e-4)
    {
        tm = (lt + rt) / 2;
        PR->TPn_flash(P, tm, comp_of, comp_liq_of, comp_gas_of, vaporfra, equalconstant_of, true);

        double rho_gas = PR->rhomix_phase(0, P, tm, comp_gas_of); //kg/m3
        double rho_liq = PR->rhomix_phase(1, P, tm, comp_liq_of);
        if (rho_gas < rho && rho < rho_liq && rho_liq - rho_gas>1e-3)
            return true;
        if (rho_gas < rho && rho > rho_liq)
            rt = tm;
        else
        {
            lt = tm;
        }

    }
    T = tm;
    return false;

}
*/

solver_new::solver_new(std::string path_i) : path(path_i), dict(IFstream(path + "system/thermotableDict")()), thermoDict(IFstream(path + "system/thermo")()), thermoDictM(IFstream(path + "system/thermoMixture")()), thermo(nullptr)
{
    TPn_flag = 0;
    wordList s(dict.lookup("species"));
    species.transfer(s);
    //Info << species << endl;
    HashPtrTable<Stype> speciesThermo(thermoDict);
    speciesData.resize(species.size());
    updated = false;
    n_species = species.size();
    forAll(species, i)
    {
        speciesData.set(
            i,
            new Stype(*speciesThermo[species[i]]));
    }
    thermo = new Mtype("test", speciesData, species, thermoDictM);
}

void solver_new::reset_specie(std::vector<std::string> m_specie)
{
    if (thermo)
        delete thermo;

    wordList s(dict.lookup("species"));
    s.resize(m_specie.size());
    for (int i = 0; i < m_specie.size(); i++)
    {
        s[i] = m_specie[i];
    }
    species.transfer(s);
    //Info << species << endl;
    HashPtrTable<Stype> speciesThermo(thermoDict);
    speciesData.resize(species.size());
    forAll(species, i)
    {
        speciesData.set(
            i,
            new Stype(*speciesThermo[species[i]]));
    }
    thermo = new Mtype("test", speciesData, species, thermoDictM);
    updated = false;
    n_species = species.size();
}

const std::vector<std::string> &solver_new::specie()
{
    const speciesTable &table = thermo->species();
    specie_.resize(n_species);
    for (unsigned int i = 0; i < n_species; i++)
    {
        specie_[i] = table[i];
    }
    return specie_;
}
/*
void solver_new::update()
{

}
*/

solver_new::~solver_new()
{
    if (thermo)
        delete thermo;
}

double solver_new::Z(int flag)
{
    //Foam::scalarList comp_liq_of(X.size(), Foam::Zero);
    //Foam::scalarList comp_gas_of(X.size(), Foam::Zero);
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (unsigned int i = 0; i < comp.size(); i++)
    {
        if (flag == 1)
            comp_of[i] = comp_gas[i];
        else if (flag == 0)
            comp_of[i] = comp_liq[i];
    }

    return ((PengRobinsonMixture<multispecie<Stype>> *)thermo)->Z(P_, T_, comp_of, flag);
}

double solver_new::dZdT(int flag)
{
    //Foam::scalarList comp_liq_of(X.size(), Foam::Zero);
    //Foam::scalarList comp_gas_of(X.size(), Foam::Zero);
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }

    return ((PengRobinsonMixture<multispecie<Stype>> *)thermo)->dZdT(P_, T_, comp_of, flag);
}
void solver_new::setT(double Tout)
{
    T_ = Tout;
    updated = false;
}
void solver_new::setP(double Pout)
{
    P_ = Pout;
    updated = false;
}
void solver_new::setX(std::vector<double> Xout)
{
    Foam::scalarList comp_X_of(n_species, Foam::Zero);
    for (unsigned int i = 0; i < n_species; i++)
    {
        comp_X_of[i] = Xout[i];
    }
    thermo->setX(comp_X_of);
    updated = false;
}
void solver_new::setY(std::vector<double> Yout)
{
    Foam::scalarList comp_Y_of(n_species, Foam::Zero);
    for (unsigned int i = 0; i < n_species; i++)
    {
        comp_Y_of[i] = Yout[i];
    }
    thermo->setY(comp_Y_of);
    updated = false;
}
void solver_new::setTPn_flag(int flagout)
{
    TPn_flag = flagout;
    updated = false;
}
void solver_new::TPn_flash_old()
{
    //Foam::scalarList comp_liq_of(comp.size(), Foam::Zero);
    //Foam::scalarList comp_gas_of(comp.size(), Foam::Zero);
    //Foam::scalarList comp_of(comp.size(), Foam::Zero);
    //for (unsigned int i = 0;i < comp.size();i++)
    //{
    //    comp_of[i] = comp[i];
    //}
    //thermo->setX(comp_of);
    Mtype::solution ret(thermo->Mtype::TPn_flash(P_, T_)());
    comp_liq.resize(comp.size());
    comp_gas.resize(comp.size());
    equalconstant.resize(comp.size());
    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_liq[i] = ret.X_liq()[i];
        comp_gas[i] = ret.X_gas()[i];
        equalconstant[i] = ret.equalconstant()[i];
    }
    //vaporfra = ret.vaporfra;
}

void solver_new::TPn_flash_update()
{
    if (!updated)
        TPn_flash(TPn_flag);
}
void solver_new::TPn_flash(int flag)
{
    if (flag == TPN_old)
        //Mtype::solution ret(thermo->Mtype::TPn_flash(P_, T_)());
        sol = thermo->Mtype::TPn_flash(P_, T_)();
    else if (flag == TPN_v2)
        sol = thermo->Mtype::TPn_flash_New(P_, T_)();
    else if (flag == TPN_TPD)
        sol = thermo->Mtype::TPn_flash_New_TPD(P_, T_)();
    else if (flag == TPN_TPD_Tud)
        sol = thermo->Mtype::TPn_flash_New_TPD_Tudisco(P_, T_)();
    if (flag == TPn_flag)
        updated = true;
    else
        updated = false;
}

void solver_new::TPn_flash_New_TPD()
{
    Foam::scalarList comp_liq_of(comp.size(), Foam::Zero);
    Foam::scalarList comp_gas_of(comp.size(), Foam::Zero);
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution ret(thermo->Mtype::TPn_flash_New_TPD(P_, T_)());
    comp_liq.resize(comp.size());
    comp_gas.resize(comp.size());
    equalconstant.resize(comp.size());
    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_liq[i] = ret.X_liq()[i];
        comp_gas[i] = ret.X_gas()[i];
        equalconstant[i] = ret.equalconstant()[i];
    }
    //vaporfra = ret.vaporfra;
}

void solver_new::TPn_flash_New()
{
    Foam::scalarList comp_liq_of(comp.size(), Foam::Zero);
    Foam::scalarList comp_gas_of(comp.size(), Foam::Zero);
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution ret(thermo->Mtype::TPn_flash_New(P_, T_)());
    comp_liq.resize(comp.size());
    comp_gas.resize(comp.size());
    equalconstant.resize(comp.size());
    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_liq[i] = ret.X_liq()[i];
        comp_gas[i] = ret.X_gas()[i];
        equalconstant[i] = ret.equalconstant()[i];
    }
    //vaporfra = ret.vaporfra;
}
std::vector<double> solver_new::Z_coe()
{
    std::vector<double> A(3);
    std::tie(A[0], A[1], A[2]) = thermo->Mtype::Z_coe(P_, T_, thermo->X_);
    return A;
}
std::vector<double> solver_new::Z_coe(std::vector<double> Xout)
{
    std::vector<double> A(3);
    Foam::scalarList comp_of(Xout.size(), Foam::Zero);
    for (unsigned int i = 0; i < Xout.size(); i++)
    {
        comp_of[i] = Xout[i];
    }
    std::tie(A[0], A[1], A[2]) = thermo->Mtype::Z_coe(P_, T_, comp_of);
    return A;
}

double solver_new::A()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }

    return ((PengRobinsonMixture<multispecie<Stype>> *)thermo)->A(P_, T_, comp_of);
}
double solver_new::dAdT()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }

    return ((PengRobinsonMixture<multispecie<Stype>> *)thermo)->dAdT(P_, T_, comp_of);
}

double solver_new::B()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }

    return ((PengRobinsonMixture<multispecie<Stype>> *)thermo)->B(P_, T_, comp_of);
}
double solver_new::dBdT()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }

    return ((PengRobinsonMixture<multispecie<Stype>> *)thermo)->dBdT(P_, T_, comp_of);
}

double solver_new::W(std::vector<double> &in)
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = in[i];
    }
    return thermo->W(comp_of);
}

void solver_new::fugacityCoefficient(int flag, std::vector<double> &in)
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = in[i];
    }
    Foam::autoPtr<scalarList> pret(thermo->PengRobinsonMixture<multispecie<Stype>>::fugacityCoefficient(P_, T_, comp_of, flag));
    ret.resize(comp.size());
    for (unsigned int i = 0; i < comp.size(); i++)
    {
        ret[i] = pret()[i];
    }
}
double solver_new::Ha()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);

    return thermo->Ha(P_, T_);
}

double solver_new::Es()
{
    return Hs() - P_ / rho();
}
double solver_new::Ha_singlePhase(int flag, std::vector<double> &in)
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = in[i];
    }
    return thermo->PengRobinsonMixture<multispecie<Stype>>::Ha(P_, T_, comp_of, flag);
}
double solver_new::dHadT_singlePhase(int flag, std::vector<double> &in)
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = in[i];
    }
    return thermo->PengRobinsonMixture<multispecie<Stype>>::dHadT(P_, T_, comp_of, flag);
}

double solver_new::Hideal(std::vector<double> &in)
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = in[i];
    }
    return thermo->PengRobinsonMixture<multispecie<Stype>>::Hideal(P_, T_, comp_of);
}
double solver_new::dHidealdT(std::vector<double> &in)
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = in[i];
    }
    return thermo->PengRobinsonMixture<multispecie<Stype>>::dHidealdT(P_, T_, comp_of);
}

double solver_new::Cp()
{
    TPn_flash_update();
    return thermo->Cp(P_, T_, sol);
}

double solver_new::c()
{
    TPn_flash_update();
    return thermo->c_opt(P_, T_, rho(), sol);
}

double solver_new::kappaS()
{
    TPn_flash_update();
    return thermo->kappaS(P_, T_, sol);
}
double solver_new::kappaT()
{
    TPn_flash_update();
    return thermo->kappaT(P_, T_, sol);
}
double solver_new::alphaP()
{
    TPn_flash_update();
    return thermo->alphaP(P_, T_, sol);
}
void solver_new::setKinit(const std::vector<double> &Kinit)
{
    if (Kinit.size() != n_species)
    {
        WarningInFunction
            << "Kinit size = "
            << Kinit.size()
            << ", which is not equal to number of species: "
            << n_species
            << ". \'inputK\' is reset to false."
            << endl;
    }
    else
    {
        thermo->inputK = true;
        thermo->Kinit.resize(n_species);
        for (unsigned int i = 0; i < n_species; i++)
        {
            thermo->Kinit[i] = Kinit[i];
        }
    }
}
const std::vector<double> &solver_new::K()
{
    TPn_flash_update();
    equalconstant.resize(n_species);
    for (unsigned int i = 0; i < n_species; i++)
    {
        equalconstant[i] = sol.equalconstant()[i];
    }
    return equalconstant;
}
double solver_new::dHadP()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution ret(thermo->Mtype::TPn_flash(P_, T_)());
    return thermo->dHadP(P_, T_, ret);
}

double solver_new::dHadXi(int di)
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution ret(thermo->Mtype::TPn_flash(P_, T_)());
    return thermo->dHadXi(P_, T_, di, ret);
}

double solver_new::dHsdXi(int di)
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution ret(thermo->Mtype::TPn_flash(P_, T_)());
    return thermo->dHsdXi(P_, T_, di, ret);
}

std::vector<double> solver_new::eps(std::vector<double> K)
{
    //autoPtr<scalarList> ret, Kof;
    Foam::scalarList ret(K.size(), Foam::Zero), Kof(K.size(), Foam::Zero);
    for (unsigned int i = 0; i < ret.size(); i++)
    {
        Kof[i] = K[i];
    }
    bool suc = thermo->TPn_flash_Matheis_test(P_, T_, thermo->X_, Kof, ret);
    if (!suc)
        for (unsigned int i = 0; i < ret.size(); i++)
        {
            ret[i] = 0;
        }
    std::vector<double> ret2(K.size());
    for (unsigned int i = 0; i < ret.size(); i++)
    {
        ret2[i] = ret[i];
    }
    return ret2;
}
void solver_new::Ln_fugacityCoefficient()
{
    autoPtr<scalarList> fugcoef;
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    double z = thermo->Z_gibbs(P_, T_, comp_of);
    fugcoef.reset((thermo->ln_fugacityCoefficient(P_, T_, z, comp_of)).ptr());
    ret.resize(comp.size());
    for (unsigned int i = 0; i < comp.size(); i++)
    {
        ret[i] = fugcoef()[i];
    }
}
void solver_new::Ln_fugacityCoefficient(int flag)
{
    autoPtr<scalarList> fugcoef;
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    fugcoef.reset((thermo->fugacityCoefficient(P_, T_, comp_of, flag)).ptr());
    ret.resize(comp.size());
    for (unsigned int i = 0; i < comp.size(); i++)
    {
        ret[i] = ::log(fugcoef()[i]);
    }
}

void solver_new::ddT_Ln_fugacityCoefficient(int flag)
{
    autoPtr<scalarList> ddT_Ln_fugcoef;
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    ddT_Ln_fugcoef.reset((thermo->ddT_Ln_fugacityCoefficient(P_, T_, comp_of, flag)).ptr());
    ret.resize(comp.size());
    for (unsigned int i = 0; i < comp.size(); i++)
    {
        ret[i] = ddT_Ln_fugcoef()[i];
    }
}

void solver_new::ddxi_Ln_fugacityCoefficient(int di, int flag)
{
    autoPtr<scalarList> ddT_Ln_fugcoef;
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    ddT_Ln_fugcoef.reset((thermo->ddxi_Ln_fugacityCoefficient(P_, T_, di, comp_of, flag)).ptr());
    ret.resize(comp.size());
    for (unsigned int i = 0; i < comp.size(); i++)
    {
        ret[i] = ddT_Ln_fugcoef()[i];
    }
}
void solver_new::setY()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setY(comp_of);
    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp[i] = thermo->X_[i];
    }
}
void solver_new::dvidT()
{
    Foam::scalarList comp_liq_of(comp.size(), Foam::Zero);
    Foam::scalarList comp_gas_of(comp.size(), Foam::Zero);
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }

    //std::cout << "!!comp";
    //for (unsigned int i = 0;i < comp.size();i++)
    //{
    //    std::cout << comp[i] << " ";
    //}
    //std::cout << std::endl;
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P_, T_)());
    //for (unsigned int i = 0;i < comp.size();i++)
    //{
    //    std::cout << "liq=" << so.X_liq()[i] << ",gas=" << so.X_gas()[i] << std::endl;
    //}
    solu.reset((thermo->dvidT(P_, T_, so)).ptr());

    ret.resize(comp.size());
    for (unsigned int i = 0; i < comp.size(); i++)
    {
        ret[i] = solu()[i];
    }
}

double solver_new::muideal_Mole(int di)
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    return thermo->muideal_Mole(P_, T_, di, &comp_of);
}

bool solver_new::solveTPD_BFGS()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    bool rb, isVapor;
    Foam::autoPtr<Foam::scalarList> rl;
    std::tie(rb, rl, isVapor) = thermo->solveTPD_BFGS(P_, T_);
    return rb;
}

double solver_new::Gideal()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    return thermo->Gideal(P_, T_, comp_of);
}

double solver_new::Gideal_Mole()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    return thermo->Gideal_Mole(P_, T_, comp_of);
}

double solver_new::G_Mole()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    return thermo->G_Mole(P_, T_, comp_of);
}

double solver_new::G()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    return thermo->G(P_, T_);
}
/*
double solver_new::S()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    return thermo->S(P_, T_);
}
*/
double solver_new::G_departure_Mole()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    return thermo->G_departure_Mole(P_, T_, comp_of);
}

double solver_new::Gibbs_single()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    return thermo->G_TPD(P_, T_, comp_of);
}

double solver_new::A_single()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    return thermo->A_TPD(P_, T_, comp_of);
}
double solver_new::z_single()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    return thermo->Z_gibbs(P_, T_, comp_of);
}
void solver_new::dvidP()
{
    Foam::scalarList comp_liq_of(comp.size(), Foam::Zero);
    Foam::scalarList comp_gas_of(comp.size(), Foam::Zero);
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }

    //std::cout << "!!comp";
    //for (unsigned int i = 0;i < comp.size();i++)
    //{
    //    std::cout << comp[i] << " ";
    //}
    //std::cout << std::endl;
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P_, T_)());
    //for (unsigned int i = 0;i < comp.size();i++)
    //{
    //    std::cout << "liq=" << so.X_liq()[i] << ",gas=" << so.X_gas()[i] << std::endl;
    //}
    solu.reset((thermo->dvidP(P_, T_, so)).ptr());

    ret.resize(comp.size());
    for (unsigned int i = 0; i < comp.size(); i++)
    {
        ret[i] = solu()[i];
    }
}

void solver_new::dvidXi(int di)
{
    Foam::scalarList comp_liq_of(comp.size(), Foam::Zero);
    Foam::scalarList comp_gas_of(comp.size(), Foam::Zero);
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P_, T_)());

    solu.reset((thermo->dvidXi(P_, T_, di, so)).ptr());

    ret.resize(comp.size());
    for (unsigned int i = 0; i < comp.size(); i++)
    {
        ret[i] = solu()[i];
    }
}

double solver_new::dTdP_HP()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P_, T_)());
    return thermo->dTdP_HP(P_, T_, so);
}
double solver_new::dEdT()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P_, T_)());
    return thermo->dEdT(P_, T_, so);
}

double solver_new::dEdP()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P_, T_)());
    return thermo->dEdP(P_, T_, so);
}

double solver_new::dEdXi(int di)
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P_, T_)());
    return thermo->dEdXi(P_, T_, di, so);
}

double solver_new::dTdE_rhoX()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P_, T_)());
    return thermo->dTdE_rhoX(P_, T_, so);
}

double solver_new::dTdrho_EX()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P_, T_)());
    return thermo->dTdrho_EX(P_, T_, so);
}

double solver_new::dPdE_rhoX()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P_, T_)());
    return thermo->dPdE_rhoX(P_, T_, so);
}

double solver_new::dPdrho_EX()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P_, T_)());
    return thermo->dPdrho_EX(P_, T_, so);
}

double solver_new::dvfdE_rhoX()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P_, T_)());
    return thermo->dvfdE_rhoX(P_, T_, so);
}

double solver_new::dvfdrho_EX()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P_, T_)());
    return thermo->dvfdrho_EX(P_, T_, so);
}

double solver_new::dTdH_HP()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P_, T_)());
    return thermo->dTdH_HP(P_, T_, so);
}
double solver_new::dTdXi_HP(int di)
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P_, T_)());
    return thermo->dTdXi_HP(P_, T_, di, so);
}

double solver_new::dTdXi_Erho(int di)
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P_, T_)());
    return thermo->dTdXi_Erho(P_, T_, di, so);
}

double solver_new::dPdXi_Erho(int di)
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P_, T_)());
    return thermo->dPdXi_Erho(P_, T_, di, so);
}

double solver_new::dvfdXi_Erho(int di)
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P_, T_)());
    return thermo->dvfdXi_Erho(P_, T_, di, so);
}

double solver_new::vaporfra()
{
    TPn_flash_update();
    return sol.vaporfra;
}

std::vector<double> solver_new::X_gas()
{
    TPn_flash_update();
    std::vector<double> ret(sol.X_gas().size());
    for (unsigned int i = 0; i < sol.X_gas().size(); i++)
    {
        ret[i] = sol.X_gas()[i];
    }
    return ret;
}
std::vector<double> solver_new::X_liq()
{
    TPn_flash_update();
    std::vector<double> ret(sol.X_gas().size());
    for (unsigned int i = 0; i < sol.X_gas().size(); i++)
    {
        ret[i] = sol.X_liq()[i];
    }
    return ret;
}

double solver_new::rho()
{
    TPn_flash_update();
    return thermo->rho(P_, T_, sol);
    /*
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (unsigned int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P_, T_)());
    return thermo->rho(P_, T_);
    */
}

double solver_new::S()
{
    return thermo->S(P_, T_);
}

double solver_new::W()
{

    return thermo->W(thermo->X_);
}

double solver_new::mu()
{
    TPn_flash_update();
    return thermo->mu_rho(P_, T_, thermo->rho(P_, T_, sol));
}

std::vector<double> solver_new::Ln_fugacityCoefficient(std::vector<double> Xout, int flag)
{
    autoPtr<scalarList> fugcoef;
    Foam::scalarList comp_of(Xout.size(), Foam::Zero);
    for (unsigned int i = 0; i < Xout.size(); i++)
    {
        comp_of[i] = Xout[i];
    }
    fugcoef.reset((thermo->Ln_fugacityCoefficient(P_, T_, comp_of, flag)).ptr());
    std::vector<double> ret(comp_of.size());

    for (unsigned int i = 0; i < comp_of.size(); i++)
    {
        ret[i] = fugcoef()[i];
    }
    return ret;
}

std::vector<double> solver_new::ddT_Ln_fugacityCoefficient(std::vector<double> Xout, int flag)
{
    autoPtr<scalarList> fugcoef;
    Foam::scalarList comp_of(Xout.size(), Foam::Zero);
    for (unsigned int i = 0; i < Xout.size(); i++)
    {
        comp_of[i] = Xout[i];
    }
    fugcoef.reset((thermo->ddT_Ln_fugacityCoefficient(P_, T_, comp_of, flag)).ptr());
    std::vector<double> ret(comp_of.size());

    for (unsigned int i = 0; i < comp_of.size(); i++)
    {
        ret[i] = fugcoef()[i];
    }
    return ret;
}

double solver_new::alphah_dev()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    //Mtype::solution so(thermo->Mtype::TPn_flash(P_, T_)());
    return thermo->alphah_dev(P_, T_);
}

double solver_new::mu_dev()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    //Mtype::solution so(thermo->Mtype::TPn_flash(P_, T_)());
    return thermo->mu_dev(P_, T_);
}

double solver_new::Dimix(int di)
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    //Mtype::solution so(thermo->Mtype::TPn_flash(P_, T_)());
    return thermo->Dimix(P_, T_, di);
}

double solver_new::dZdT()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P_, T_)());
    return thermo->dZdT(P_, T_, so);
}
double solver_new::dSdT()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P_, T_)());
    return thermo->dSdT(P_, T_, so);
}
double solver_new::dSdP()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P_, T_)());
    return thermo->dSdP(P_, T_, so);
}
/*
double solver_new::Z()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P_, T_)());
    return thermo->Z(P_, T_, so);
}*/
double solver_new::drhodT()
{
    TPn_flash_update();
    return thermo->drhodT(P_, T_, sol);
}

double solver_new::drhodP()
{
    TPn_flash_update();
    return thermo->drhodP(P_, T_, sol);
}
double solver_new::dZdP()
{
    TPn_flash_update();
    return thermo->dZdP(P_, T_, sol);
}
double solver_new::Z()
{
    TPn_flash_update();
    return thermo->Z(P_, T_, sol);
}

std::vector<std::vector<double>> solver_new::dTHvfc_G_rhoY_dXrhoP()
{
    TPn_flash_update();
    std::vector<std::vector<double>> Grad;
    Grad.resize(n_species + 2, std::vector<double>(n_species + 5));
    autoPtr<scalarRectangularMatrix> grad(thermo->dTHvfc_G_rhoY_dXrhoP(P_, T_, sol));
    for (int i = 0; i < n_species + 2; i++)
    {
        for (int j = 0; j < n_species + 5; j++)
        {
            Grad[i][j] = grad()[i][j];
        }
    }
    return Grad;
}
std::vector<std::vector<double>> solver_new::dErhovfc_G_rhoY_dXTP()
{
    TPn_flash_update();
    std::vector<std::vector<double>> Grad;
    Grad.resize(n_species + 2, std::vector<double>(n_species + 5));
    autoPtr<scalarRectangularMatrix> grad(thermo->dErhovfc_G_rhoY_dXTP(P_, T_, sol));
    for (int i = 0; i < n_species + 2; i++)
    {
        for (int j = 0; j < n_species + 5; j++)
        {
            Grad[i][j] = grad()[i][j];
        }
    }
    return Grad;
}

double solver_new::Hs()
{
    TPn_flash_update();
    return thermo->Hs(P_, T_, sol);
}

double solver_new::rho_G()
{
    TPn_flash_update();
    return P_ * thermo->W(sol.X_gas()) / (thermo->singlePhaseMixtureThermoType::Z(P_, T_, sol.X_gas(), 1) * RR * 1.0e-03 * T_);
}

double solver_new::W_G()
{
    TPn_flash_update();
    return thermo->W(sol.X_gas());
}

std::vector<double> solver_new::Y_G()
{
    TPn_flash_update();
    double W_G = thermo->W(sol.X_gas());
    std::vector<double> ret(n_species);
    for (int i = 0; i < n_species; i++)
    {
        ret[i] = (*thermo)[i].W() * 1e-3 * sol.X_gas()[i] / W_G;
    }
    return ret;
}

double solver_new::drhodXi(int di)
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P_, T_)());
    return thermo->drhodXi(P_, T_, di, so);
}

double solver_new::dZdXi(int di)
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P_, T_)());
    return thermo->dZdXi(P_, T_, di, so);
}

double solver_new::drhoPdH_HP()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P_, T_)());
    return thermo->drhoPdH_HP(P_, T_, so);
}
double solver_new::drhoPdP_HP()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P_, T_)());
    return thermo->drhoPdP_HP(P_, T_, so);
}

double solver_new::drhoPdH_HsP()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P_, T_)());
    return thermo->drhoPdH_HsP(P_, T_, so);
}

double solver_new::dTdXi_HsP(int di)
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P_, T_)());
    return thermo->dTdXi_HsP(P_, T_, di, so);
}
double solver_new::drhoPdP_HsP()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P_, T_)());
    return thermo->drhoPdP_HsP(P_, T_, so);
}

double solver_new::drhodP_HP()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P_, T_)());
    return thermo->drhodP_HP(P_, T_, so);
}
double solver_new::drhoPdXi_HP(int di)
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P_, T_)());
    return thermo->drhoPdXi_HP(P_, T_, di, so);
}
double solver_new::drhoPdXi_HsP(int di)
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P_, T_)());
    return thermo->drhoPdXi_HsP(P_, T_, di, so);
}
double solver_new::drhodXi_HP(int di)
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P_, T_)());
    return thermo->drhodXi_HP(P_, T_, di, so);
}

double solver_new::dvfdP_HP()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P_, T_)());
    return thermo->dvfdP_HP(P_, T_, so);
}
double solver_new::dvfdH_HP()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P_, T_)());
    return thermo->dvfdH_HP(P_, T_, so);
}

double solver_new::dvfdXi_HP(int di)
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P_, T_)());
    return thermo->dvfdXi_HP(P_, T_, di, so);
}

double solver_new::dvfdP_HsP()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P_, T_)());
    return thermo->dvfdP_HsP(P_, T_, so);
}
double solver_new::dvfdH_HsP()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P_, T_)());
    return thermo->dvfdH_HsP(P_, T_, so);
}

double solver_new::dvfdXi_HsP(int di)
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P_, T_)());
    return thermo->dvfdXi_HsP(P_, T_, di, so);
}

double solver_new::T_HsP(double h, double P, double T0)
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);

    double Ttemp = T0;
    double htemp = thermo->Hs(P_, Ttemp);
    while (fabs(htemp - h) > 1e-5)
    {
        //cout << Ttemp << "," << htemp << std::endl;
        Mtype::solution so(thermo->Mtype::TPn_flash(P_, Ttemp)());
        cout << thermo->Cp(P_, Ttemp) << "," << thermo->Cp_Hs(P_, Ttemp, so) << std::endl;
        Ttemp -= (htemp - h) / thermo->Cp_Hs(P_, Ttemp, so);
        htemp = thermo->Hs(P_, Ttemp);
    }
    return Ttemp;
}
double solver_new::P()
{
    return P_;
}
double solver_new::T()
{
    return T_;
}
const std::vector<double> &solver_new::X()
{
    X_.resize(n_species);
    for (unsigned int i = 0; i < n_species; i++)
    {
        X_[i] = thermo->X_[i];
    }
    return X_;
}

void solver_new::drhoPdXHP_HsP()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (unsigned int i = 0; i < comp.size(); i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P_, T_)());
    autoPtr<scalarList> grad = thermo->drhoPdXHP_HsP(P_, T_, so);
    ret.resize(grad().size());
    for (unsigned int i = 0; i < ret.size(); i++)
    {
        ret[i] = grad()[i];
    }
}

double solver_new::S(int flag)
{
    //Foam::scalarList comp_liq_of(X.size(), Foam::Zero);
    //Foam::scalarList comp_gas_of(X.size(), Foam::Zero);
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (unsigned int i = 0; i < comp.size(); i++)
    {
        if (flag == 1)
            comp_of[i] = comp_gas[i];
        else if (flag == 0)
            comp_of[i] = comp_liq[i];
    }

    return ((PengRobinsonMixture<multispecie<Stype>> *)thermo)->S(P_, T_, comp_of, flag);
}

double solver_new::dSdT(int flag)
{
    //Foam::scalarList comp_liq_of(X.size(), Foam::Zero);
    //Foam::scalarList comp_gas_of(X.size(), Foam::Zero);
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (unsigned int i = 0; i < comp.size(); i++)
    {
        if (flag == 1)
            comp_of[i] = comp_gas[i];
        else if (flag == 0)
            comp_of[i] = comp_liq[i];
    }

    return ((PengRobinsonMixture<multispecie<Stype>> *)thermo)->dSdT(P_, T_, comp_of, flag);
}

double solver_new::dSdP(int flag)
{
    //Foam::scalarList comp_liq_of(X.size(), Foam::Zero);
    //Foam::scalarList comp_gas_of(X.size(), Foam::Zero);
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (unsigned int i = 0; i < comp.size(); i++)
    {
        if (flag == 1)
            comp_of[i] = comp_gas[i];
        else if (flag == 0)
            comp_of[i] = comp_liq[i];
    }

    return ((PengRobinsonMixture<multispecie<Stype>> *)thermo)->dSdP(P_, T_, comp_of, flag);
}
