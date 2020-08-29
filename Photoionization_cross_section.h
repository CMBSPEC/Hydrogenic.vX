//======================================================================================
// Author: Jens Chluba 
// last modification: Oct 2010
// purpose: compute photoionization cross sections for hydrogenic atoms
// 
// The routines are based on Storey & Hummer, 1991. Notation and 
// recursion relations can be found there. However, we modified these 
// slightly to achieve better stability.
//======================================================================================
// 06.03.2017: added rescaling with alpha and me for SH routine [JC]
// 04.08.2014: added gaunt factor version that returns values for l=0..n-1

#ifndef PHOTO_IONIZATION_CROSS_SECTION_H
#define PHOTO_IONIZATION_CROSS_SECTION_H

#include <vector>
#include <iostream>

#include "physical_consts.h"

using namespace std;

//======================================================================================
//
// class with basic properties of photonionization cross sections
//
//======================================================================================
class Photoionization_cross_section
{
private:
    
protected:
    //==================================================================================
    // important data that is set when initializing
    //==================================================================================
    int nn, ll;
    int nc;                       // above this value for n simplified versions for the cross-sections are used
                                  // e.g. Kramers-formula with g=1 

    int Z;                        // atomic charge
    double Np;                    // effective number of protons in the nucleus (Np=M_nucleus/mp)
    double mu_red;                // reduced mass in units of me
    double Ry;                    // Rydberg of the atom: R_inf/(1+me/M) in ergs
    double Ry_eV;                 // Rydberg of the atom: R_inf/(1+me/M) in eV
    double Ry_Hz;                 // Rydberg of the atom: R_inf/(1+me/M) in Hz
    
    double nu_ionization, E_ionization, E_ionization_eV;
    
    //==================================================================================
    // auxilliary variables
    //==================================================================================
    int mess_flag;
    
    double E_n(int n){ return Z*Z*Ry/n/n; }                                   // in ergs
    double E_n_eV(int n){ return Z*Z*Ry_eV/n/n; }                             // in eV
    double nu_n(int n){ return Z*Z*Ry_Hz/n/n; }                               // in Hz
    double nu2E_e(double nu){ return (nu/nu_ionization-1.0)/nn/nn; }          // free electron energy 
                                                                              // in Rydbergs of the atom
                                                                              // with charge Z, and Np protons 
                                                                              // i.e. R_inf*Z^2/(1+me/(Np*mp))
    
    double E_e2nu(double E){ return (E*nn*nn+1.0)*nu_ionization; }            // emitted photon energy in Hz
    double nu_threshold();
    
    void Set_Ry(double Np);
    
public:
    
    //==================================================================================
    // Konstructors and Destructors
    //==================================================================================
    Photoionization_cross_section(){}
    Photoionization_cross_section(int nv, int lv, int Z, double Np, int mflag=1);
    ~Photoionization_cross_section(){}
    void init(int nv, int lv, int Z, double Np, int mflag=1);
    
    //==================================================================================
    int Get_n(){return nn;}
    int Get_nc(){return nc;}                                     // for n>nc the Kramers-formula will be applied
    int Get_l(){return ll;}
    int Get_Z(){return Z;}
    double Get_mu_red(){return mu_red;}
    double Get_Np(){return Np;}
    double Get_nu_ionization(){ return nu_ionization;}           // in Hz
    double Get_E_ionization(){ return E_ionization;}             // in ergs
    double Get_E_ionization_eV(){ return E_ionization_eV;}       // in eV
    
    //==================================================================================
    double sig_phot_ion_Kramers(double nu);
    double sig_phot_ion_Kramers_E(double E){ return sig_phot_ion_Kramers(E_e2nu(E));}
};



//======================================================================================
//
// class to calculate the photoionization cross sections for given n & l
//
//======================================================================================
class Photoionization_cross_section_SH: public Photoionization_cross_section
{
private:
    
    vector<double> Rp1;
    vector<double> Rm1;
    
    //==================================================================================
    // to norm of the emission profile (int sig_nl_c(nu) d nu )
    //==================================================================================
    double sig_phot_ion(const vector<double> &Rp1, const vector<double> &Rm1);
    
    double sigma_nucval;
    double gaunt_nucval;

    //==================================================================================
    // added by JC to rescale cross section with alp and me
    //==================================================================================
    double FSC_scale;          // Scaling of the fine structure constant alpha/alpha_ref
    double ME_scale;           // Scaling of the electron mass me/me_ref
    double energy_scale;       // Scaling of the energy based on FSC and ME scaling
    double sig_scale;          // Scaling of the photoionisation cross section

public:
    
    //==================================================================================
    // Konstructors and Destructors
    //==================================================================================
    Photoionization_cross_section_SH(){}
    Photoionization_cross_section_SH(int nv, int lv, int Z, double Np, int mflag=1);
    ~Photoionization_cross_section_SH(){}
    void init(int nv, int lv, int Z, double Np, int mflag=1);
    void clear(){ Rp1.clear(); Rm1.clear(); return; }
    
    //==================================================================================
    double sig_phot_ion(double nu);
    double sig_phot_ion_nuc(){ return sigma_nucval*sig_scale; };
    double gaunt_nuc(){ return gaunt_nucval; };
    double nu_sig_phot_ion_small(double eps);
    
    double g_phot_ion(double nu);

    //===================================================================================
    // return all photoionization cross section gaunt factor values at nu for l=0..n-1
    //===================================================================================
    void g_phot_ion(double nu, vector<double> &g_phot_l);

    //==================================================================================
    // same as above but as a function of the free electron energy
    //==================================================================================
    double sig_phot_ion_E(double E){ return sig_phot_ion(E_e2nu(E));}
    double g_phot_ion_E(double E){ return g_phot_ion(E_e2nu(E));}
    
    //===================================================================================
    // to include variations of alpha & me [March, 2017, JC]
    //===================================================================================
    void rescale_phot(double alp_scale, double me_scale);
    void reset_phot(){ this->rescale_phot(1.,1.); }
    
    // Access functions for the fine structure constant scaling and the mass scaling
    double get_FSC_scale() const { return this->FSC_scale; }
    double get_ME_scale() const { return this->ME_scale; }

    // override relevant base functions
    double Get_nu_ionization();           // in Hz
    double Get_E_ionization();            // in ergs
    double Get_E_ionization_eV();         // in eV
    double sig_phot_ion_Kramers(double nu);
};



//======================================================================================
//
// class to calculate the photoionization cross sections for ground state of hydrogen
//
//======================================================================================
class Photoionization_Lyc
{
private:
    
protected:
    double nu_ionization, E_ionization, E_ionization_eV;
    
    //==================================================================================
    // auxilliary variables
    //==================================================================================
    int mess_flag;
    
    double E_n(){ return const_EH_inf_ergs/(1.0+const_me_mp); }  // in ergs
    double E_n_eV(){ return const_EH_inf/(1.0+const_me_mp); }    // in eV
    double nu_n(){ return const_EH_inf_Hz/(1.0+const_me_mp); }   // in Hz
    double nu2E_e(double nu){ return (nu/nu_ionization-1.0); }   // free electron energy in Rydbergs of the atom
    double E_e2nu(double E){ return (E+1.0)*nu_ionization; }     // emitted photon energy in Hz
    
    double nu_threshold();
    
public:
    
    //==================================================================================
    // Konstructors and Destructors
    //==================================================================================
    Photoionization_Lyc(){}
    Photoionization_Lyc(int mflag=1);
    ~Photoionization_Lyc(){}
    void init(int mflag=1);
    
    //==================================================================================
    double Get_nu_ionization(){ return nu_ionization;}           // in Hz
    double Get_E_ionization(){ return E_ionization;}             // in ergs
    double Get_E_ionization_eV(){ return E_ionization_eV;}       // in eV
    
    double sig_phot_ion_Kramers(double nu);
    double sig_phot_ion_Kramers_E(double E){ return sig_phot_ion_Kramers(E_e2nu(E));}
    
    double sig_phot_ion(double nu);
    double sig_phot_ion_nuc(){ return sig_phot_ion(Get_nu_ionization()); };
    double g_phot_ion(double nu);
    
    //==================================================================================
    // same as above but as a function of the free electron energy
    //==================================================================================
    double sig_phot_ion_E(double E){ return sig_phot_ion(E_e2nu(E));}
    double g_phot_ion_E(double E){ return g_phot_ion(E_e2nu(E));}
};

#endif
//======================================================================================
//======================================================================================
