#include "dictionary.H"
#include "IFstream.H"
#include "specie.H"
#include "perfectGas.H"
#include "hConstThermo.H"
#include "speciesTable.H"
#include "sensibleEnthalpy.H"
#include "thermo.H"
#include "constTransport.H"
//#include "PengRobinsonS.H"
//#include "PengRobinsonM.H"


#include <vector>
#include <string>
#include <iostream>
#include <cmath>


#include "dictionary.H"
#include "IFstream.H"
#include "specie.H"
#include "multispecie.H"
#include "perfectGas.H"
#include "hConstThermo.H"
#include "speciesTable.H"
#include "sensibleEnthalpy.H"
#include "thermo.H"
#include "constTransport.H"
#include "PengRobinson.H"
#include "PengRobinsonMixture.H"
#include "VLE.H"
#include "janafThermoAd.H"
#include "chungTransport.H"
#include "chungTransportMixture.H"
#include "HashPtrTable.H"
#include "multithermo.H"
#include "IOmanip.H"

/*struct solver
{
    std::string path;
    double P;
    double T;
    std::vector<double> X;
    std::vector<std::string> specie;
    Foam::speciesTable* s;
    Foam::dictionary* dict;
    void update();
    Foam::PengRobinsonM<Foam::specie>* PR;
    solver(std::string);
    ~solver();
    void solve(bool flag);
    double vaporfra;
    double density();
    double Z();
    bool twophase(double rho, double lt, double rt);
    std::vector<double> equalconstant;
    std::vector<double> comp_liq;
    std::vector<double> comp_gas;
};
*/
typedef chungTransport<species::thermo<janafThermoAd<PengRobinson<specie>>, sensibleEnthalpy>> Stype;

const int TPN_old = 0;
const int TPN_v2 = 1;
const int TPN_TPD = 2;
const int TPN_TPD_Tud = 3;
struct solver_new
{
    //typedef chungTransport<species::thermo<janafThermoAd<PengRobinson<specie>>, sensibleEnthalpy>> Stype;
    private:
    typedef species::multithermo<VLE<chungTransportMixture<PengRobinsonMixture<multispecie<Stype>>>>, sensibleEnthalpy> Mtype;
    std::string path;
    Mtype::solution sol;
    bool updated;
    dictionary dict;
    speciesTable species;
    int n_species;
    int TPn_flag;
    dictionary thermoDict;
    dictionary thermoDictM;
    Mtype* thermo;
    PtrList<Stype> speciesData;
    double P_;
    double T_;

    
    void setY();
    //std::vector<double> X;
    //std::vector<std::string> m_specie;
    //Foam::speciesTable* s;
    //Foam::dictionary* dict;
    //void update();
    //Foam::PengRobinsonM<Foam::specie>* PR;
    
    public:
    solver_new(std::string file);
    void reset_specie(std::vector<std::string>);

    ~solver_new();

    // settings 
    void setT(double);
    void setP(double);
    void setX(std::vector<double>);
    void setTPn_flag(int);

    // output
    double P();
    double T();
    double vaporfra();
    double rho();
    double Cp();

    double drhodP();
    double drhodT();
    const std::vector<double>& X();
    const std::vector<double>& K();
    void setKinit(const std::vector<double>&);
    const std::vector<std::string>& specie();

    private:
    //void solve(bool flag);
    
    //double density();

    
    double dSdT();
    double dSdP();
    double drhodXi(int);
    double dZdT();
    double dZdXi(int);
    double drhodP_HP();
    double drhoPdH_HP();
    double drhoPdP_HP();
    double drhoPdXi_HP(int);
    double drhoPdH_HsP();
    double drhoPdP_HsP();
    double drhoPdXi_HsP(int);
    void drhoPdXHP_HsP();
    double dTdXi_HsP(int);
    double dvfdH_HP();
    double dvfdP_HP();
    double dvfdXi_HP(int);
    double c();
    double dvfdH_HsP();
    double dvfdP_HsP();
    double dvfdXi_HsP(int);
    double drhodXi_HP(int);
    double Z();
    double Gibbs_single();
    double G_departure_Mole();
    double Gideal_Mole();
    double G_Mole();
    double G();
    double S();
    bool solveTPD_BFGS();
    double A_single();
    double z_single();
    double muideal_Mole(int);
    double Gideal();
    double Z(int);
    double S(int);
    double dSdP(int);
    double dSdT(int);
    void TPn_flash(int);
    void TPn_flash_old();
    void TPn_flash_New();
    void TPn_flash_New_TPD();
    void TPn_flash_update();
    double dZdT(int);
    double A();
    double dAdT();
    double B();
    double dBdT();
    double Ha();
    double Es();
    double Hs();
    double Ha_singlePhase(int, std::vector<double>&);
    double dHadT_singlePhase(int, std::vector<double>&);
    double Hideal(std::vector<double>&);
    double dHidealdT(std::vector<double>&);

    double dHadP();
    double dHadXi(int);
    double dHsdXi(int);
    double W();
    double alphah_dev();
    double mu_dev();
    double Dimix(int);
    double W(std::vector<double>&);
    double T_HsP(double h, double p, double T0);
    void Ln_fugacityCoefficient();
    void Ln_fugacityCoefficient(int);
    void ddT_Ln_fugacityCoefficient(int);
    void ddxi_Ln_fugacityCoefficient(int, int);
    void dvidT();
    void dvidP();
    void dvidXi(int);
    void  fugacityCoefficient(int, std::vector<double>&);

    double dTdP_HP();
    double dTdH_HP();
    double dTdXi_HP(int);

    double dTdE_rhoX();
    double dTdrho_EX();
    double dPdE_rhoX();
    double dPdrho_EX();
    double dvfdE_rhoX();
    double dvfdrho_EX();
    double dTdXi_Erho(int);
    double dPdXi_Erho(int);
    double dvfdXi_Erho(int);

    double dEdT();
    double dEdP();
    double dEdXi(int);

    //double drhodT();
    //double drhodP();
    //double drhodXi(int);
    //bool twophase(double rho, double lt, double rt);
    std::vector<double> equalconstant;
    std::vector<double> comp_liq;
    std::vector<double> comp_gas;
    std::vector<double> comp;
    std::vector<double> ret;
    std::vector<double> X_;
    std::vector<std::string> specie_;
};