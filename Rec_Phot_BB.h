//====================================================================================================================
// Author: Jens Chluba 
// last modification: Aug 2020
// 
// Purpose: compute the recombination & photoionzation rates of a hydrogenic atom in a
// blackbody ambient radiation field of temperature Tg. The electrons have temperature Te.
//====================================================================================================================

//====================================================================================================================
// class based on Storey & Hummer, 1991 
// these routines have been checked for n<=300 but should work until n~1000
//====================================================================================================================
// 28.08.2020: changed spline setup to do it in two steps. This is important for parallel computations [JC]

#ifndef REC_PHOT_BB_H
#define REC_PHOT_BB_H

#include "Photoionization_cross_section.h"
#include "physical_consts.h"
#include <gsl/gsl_spline.h>
#include <vector>

using namespace std;

//====================================================================================================================
//
// class with basic information
//
//====================================================================================================================
class Rec_Phot_BB
{
private:
    
protected:

    double g_l;                                                    // statistical weight
    double xi_small;

public:
    //==================================================================================
    // Konstructors and Destructors
    //==================================================================================
    Rec_Phot_BB(){}
    ~Rec_Phot_BB(){}

    double xx2nu(double x, double T_g){ return x*T_g/const_h_kb; }
    double nu2xx(double nu, double T_g){ return const_h_kb*nu/T_g; }
};


//====================================================================================================================
//
// class based on Storey & Hummer, 1991 
// the recomination rates for given (n, l) can be calculated.
//
//====================================================================================================================
class Rec_Phot_BB_SH: public Rec_Phot_BB, public Photoionization_cross_section_SH
{
private:
    //==================================================================================
    // Photoinzation rate for bb photons at temperature T_g
    //==================================================================================
    double R_nl_c_JC(double T_g);                                         // in 1/sec 
    //==================================================================================
    // Recombination rate for bb photons at temperature T_g and rho=T_g/T_M
    //==================================================================================
    double R_c_nl_JC(double T_g, double rho);
    //==================================================================================
    // derivative with respect to Te
    //==================================================================================
    double dR_c_nl_dTe_JC(double T_g, double rho);                        // in cm^3/sec
    
public:
    //==================================================================================
    // Konstructors and Destructors
    //==================================================================================
    Rec_Phot_BB_SH();
    Rec_Phot_BB_SH(int nv, int lv, int Z, double Np, int mflag=1);          
    ~Rec_Phot_BB_SH(){}

    void init(int nv, int lv, int Z, double Np, int mflag=1);
    void clear(){ Photoionization_cross_section_SH::clear(); return; }
    
    double Get_xi_small(){ return xi_small; }
    double Get_nu_small(){ return Get_nu_ionization()*xi_small; }

    //==================================================================================
    // Photoinzation rate for bb photons at temperature T_g
    //==================================================================================
    double R_nl_c_Int(double T_g);                                        // in 1/sec
    //==================================================================================
    // Recombination rate for bb photons at temperature T_g and rho=T_g/T_M
    //==================================================================================
    double R_c_nl_Int(double T_g, double rho=1.0);                        // in cm^3/sec
    //==================================================================================
    // derivative of Recombination rate for bb photons at temperature T_g and rho=T_g/T_M
    //==================================================================================
    double dR_c_nl_dTe_Int(double T_g, double rho=1.0);                   // in cm^3/sec

    double sig_phot_ion_lim(double nu);
    double g_phot_ion_lim(double nu);
};

//====================================================================================================================
//
// using qubic spline interpolation for the Gaunt-factors
//
//====================================================================================================================
class Rec_Phot_BB_SH_QSP: public Rec_Phot_BB, public Photoionization_cross_section_SH
{
private:

    int nxi, lmin_spline;
    vector<double> xi, lgxi, gaunt;

    int spline_is_allocated, spline_is_setup;
    gsl_interp_accel *acc;
    gsl_spline *spline;
    
    int calc_data_for_spline();
    double calc_gaunt_fac_spline(double xiv);

    //==================================================================================
    // Photoinzation rate for bb photons at temperature T_g
    //==================================================================================
    double R_nl_c_JC(double T_g);                                         // in 1/sec 
    //==================================================================================
    // Recombination rate for bb photons at temperature T_g and rho=T_g/T_M
    //==================================================================================
    double R_c_nl_JC(double T_g, double rho);                             // in cm^3/sec
    //==================================================================================
    // derivative with respect to Te
    //==================================================================================
    double dR_c_nl_dTe_JC(double T_g, double rho);                        // in cm^3/sec
    
public:
    //==================================================================================
    // Konstructors and Destructors
    //==================================================================================
    Rec_Phot_BB_SH_QSP();
    Rec_Phot_BB_SH_QSP(int nv, int lv, int Z, double Np, int mflag=1);
    ~Rec_Phot_BB_SH_QSP();

    void init(int nv, int lv, int Z, double Np, int mflag=1);
    void clear();

    //==================================================================================
    // for parallel setup splines have to be done in serial. Hence first compute
    // data and then setup splines
    //==================================================================================
    void init_parallel(int nv, int lv, int Z, double Np, int mflag=1);
    void arm_spline_parallel();

    double Get_xi_small(){ return xi_small; }
    double Get_nu_small(){ return Get_nu_ionization()*xi_small; }
    //==================================================================================
    // Photoinzation rate for bb photons at temperature T_g
    //==================================================================================
    double R_nl_c_Int(double T_g);                                             // in 1/sec
    //==================================================================================
    // Recombination rate for bb photons at temperature T_g and rho=T_g/T_M
    //==================================================================================
    double R_c_nl_Int(double T_g, double rho=1.0);                             // in cm^3/sec
    //==================================================================================
    // derivative of Recombination rate for bb photons at temperature T_g and rho=T_g/T_M
    //==================================================================================
    double dR_c_nl_dTe_Int(double T_g, double rho=1.0);                        // in cm^3/sec
    
    double sig_phot_ion_lim(double nu);
    double g_phot_ion_lim(double nu);
};

#endif
//====================================================================================================================
//====================================================================================================================
