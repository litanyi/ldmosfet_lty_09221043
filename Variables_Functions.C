/*
 * Variables_Functions.C
 *
 *  Created on: May 10, 2017
 *      Author: tong
 */

#include "Variables_Functions.h"
#include "math.h"
#include "DoubleVector.h"
#include "VectorVariable.h"
#include "VectorData.h"
#include "Vector.h"

#include "material_database.h"

#include <fstream>

using namespace std;


/**
 * @brief Variables_Functions::Variables_Functions
 */
Variables_Functions::Variables_Functions(){

}

/**
 * @brief Variables_Functions::~Variables_Functions
 */
Variables_Functions::~Variables_Functions(){

}

/**
 * @brief Variables_Functions::FermiInteg: fermi integration
 * @param ef
 * @param fermi_flag
 * @param fermi_order
 * @return
 */
double Variables_Functions::FermiInteg(double ef, int fermi_flag,
                            int fermi_order){
    double fermii;
    if(fermi_order == 1){
        if(fermi_flag == 1){
            double exp_fac = exp(-0.17*pow((ef+1.0),2.0));
            double nu = pow(ef,4.0)+50.0+33.6*ef*(1.0-0.68*exp_fac);
            double zeta = 3.0*sqrt(M_PI)/(4.0*pow(nu,0.375));
            fermii = exp(ef)/(1.0+zeta*exp(ef));
        } else if(fermi_flag == 0){
            fermii = exp(ef);
        }
    } else if(fermi_order == 0){
        if(fermi_flag == 1){
            fermii = log(1.0+exp(ef));
        } else if(fermi_flag == 0){
            fermii = exp(ef);
        }
    } else if(fermi_order == -1){
        if(fermi_flag == 1){
            double exp_fac = exp(-0.17*pow((ef+1.0),2.0));
            double nu = pow(ef,4.0)+50.0+33.6*ef*(1.0-0.68*exp_fac);
            double zeta = 3.0*sqrt(M_PI)/(4.0*pow(nu,0.375));
            double nu_prime = 4.0*pow(ef,3.0)+33.6-22.848*exp_fac*(1-0.34*(ef+ef*ef));
            double zeta_prime = -(9.0*sqrt(M_PI)/32.0)*pow(nu,-11.0/8.0)*nu_prime;
            fermii = (exp(-ef)-zeta_prime)/pow((exp(-ef)+zeta),2.0);
        } else if(fermi_flag == 0){
            fermii = exp(ef);
        }
    }
    return fermii;
}

/**
 * @brief Variables_Functions::antiFermi: anti fermi integration
 * @param fermi
 * @param dummy_flag
 * @param fermi_flag
 * @param ef
 */
void Variables_Functions::antiFermi(double fermi, int dummy_flag,
                          int fermi_flag, double &ef){
    if(dummy_flag == 0){
        if(fermi_flag == 0){
            ef = log(fermi);
        } else if(fermi_flag == 1){
            ef = log(exp(fermi)-1);
        }
    } else if(dummy_flag == 1){
        if(fermi_flag == 0){
            ef = log(fermi);
        } else if(fermi_flag == 1){
            if(fermi < 8.463){
                double a = 0.35355339059327379;
                double b = -4.9500897298752622e-3;
                double c = 1.4838577128872821e-4;
                double d = -4.4256301190009895e-6;
                ef = log(fermi) + a*fermi + b*pow(fermi,2.0) +
                        c*pow(fermi,3.0) + d*pow(fermi,4.0);
            }else{
                double tmp_ef = pow(0.75*sqrt(M_PI)*fermi, 0.75) - pow(M_PI,2.0)/6;
                ef = sqrt(tmp_ef);
            }
        }
    }
}


/**
 * @brief Variables_Functions::dummyElectron: calculate electron density
 * @param volt
 * @param fermin
 * @param deriv_flag
 * @param affinity
 * @param T
 * @param nc
 * @param n_dens
 */
void Variables_Functions::dummyElectron(double volt, double fermin, bool deriv_flag,
                     double affinity, double T, double kb, double e, double nc, double &n_dens){
    double thermo_volt = kb*T/e;
    double Ec = -(volt+affinity);
    double zeta = (fermin-Ec)/thermo_volt;
    double fermii;
    if(deriv_flag){
        fermii = FermiInteg(zeta, 0, -1);
        n_dens = (nc/thermo_volt)*fermii;
    }else{
        fermii = FermiInteg(zeta, 0, 1);
        n_dens = nc*fermii;
    }
}

/**
 * @brief Variables_Functions::dummyElectron: calculate electron density
 * @param volt
 * @param fermin
 * @param deriv_flag
 * @param affinity
 * @param T
 * @param nc
 * @param n_dens
 */
adtl::AutoDScalar Variables_Functions::dummyElectron(adtl::AutoDScalar volt, double fermin, bool deriv_flag,
                     double affinity, double T, double kb, double e, double nc){
    adtl::AutoDScalar n_dens;
    double thermo_volt = kb*T/e;
    adtl::AutoDScalar Ec = -(volt+affinity);
    adtl::AutoDScalar zeta = (fermin-Ec)/thermo_volt;
    adtl::AutoDScalar fermii;
    if(deriv_flag){
        fermii = exp(zeta);
        n_dens = (nc/thermo_volt)*fermii;
    }else{
        fermii = exp(zeta);
        n_dens = nc*fermii;
    }
    return n_dens;
}

/**
 * @brief Variables_Functions::dummyHole: calculate hole density
 * @param volt
 * @param fermip
 * @param eg
 * @param deriv_flag
 * @param affinity
 * @param T
 * @param nv
 * @param p_dens
 */
adtl::AutoDScalar Variables_Functions::dummyHole(adtl::AutoDScalar volt, double fermip, double eg, bool deriv_flag,
                       double affinity, double T, double kb, double e, double nv){
    // thermal voltage
    adtl::AutoDScalar p_dens;
    double thermo_volt = kb*T/e;
    adtl::AutoDScalar  Ev= -(volt+affinity + eg);
    adtl::AutoDScalar  zeta = (Ev - fermip)/thermo_volt;
    adtl::AutoDScalar  fermii;
    if(deriv_flag){
        fermii = exp(zeta);
        p_dens = -nv/thermo_volt*fermii;
    }else {
        fermii = exp(zeta);
        p_dens = nv*fermii;
    }
    return p_dens;
}

/**
 * @brief Variables_Functions::dummyElectron: calculate electron density
 * @param volt
 * @param fermin
 * @param deriv_flag
 * @param affinity
 * @param T
 * @param nc
 * @param n_dens
 */
double Variables_Functions::dummyElectron(double volt, double fermin, bool deriv_flag,
                     double affinity, double T, double kb, double e, double nc){
    double n_dens;
    dummyElectron(volt, fermin, deriv_flag, affinity, T, kb, e, nc, n_dens);
    return n_dens;
}

/**
 * @brief Variables_Functions::dummyHole: calculate hole density
 * @param volt
 * @param fermip
 * @param eg
 * @param deriv_flag
 * @param affinity
 * @param T
 * @param nv
 * @param p_dens
 */
void Variables_Functions::dummyHole(double volt, double fermip, double eg, bool deriv_flag,
                       double affinity, double T, double kb, double e, double nv, double &p_dens){
    // thermal voltage
    double thermo_volt = kb*T/e;
    double Ev= -(volt+affinity + eg);
    double zeta = (Ev - fermip)/thermo_volt;
    double fermii;
    if(deriv_flag){
        fermii = FermiInteg(zeta,0,-1);
        p_dens = -nv/thermo_volt*fermii;
    }else {
        fermii = FermiInteg(zeta,0,1);
        p_dens = nv*fermii;
    }
}

/**
 * @brief Variables_Functions::dummyHole: calculate hole density
 * @param volt
 * @param fermip
 * @param eg
 * @param deriv_flag
 * @param affinity
 * @param T
 * @param nv
 * @param p_dens
 */
double Variables_Functions::dummyHole(double volt, double fermip, double eg, bool deriv_flag,
                       double affinity, double T, double kb, double e, double nv){
    // thermal voltage
    double p_dens;
    dummyHole(volt, fermip, eg, deriv_flag, affinity, T, kb, e, nv, p_dens);
    return p_dens;
}

/**
 * @brief Variables_Functions::antiDummyElectron: quasi fermi level for electron
 * @param n_dens
 * @param volt
 * @param affinity
 * @param T
 * @param nc
 * @param fermin
 */
void Variables_Functions::antiDummyElectron(double n_dens, double volt,
                       double affinity, double T, double kb, double e,  double nc, double &fermin){
    // thermal voltage
    double thermo_volt = kb*T/e;
    double zeta;
    antiFermi(n_dens/nc,1,0,zeta);  //zeta=log(n_dens/nc)
    fermin = -zeta*thermo_volt + volt + affinity;
}

/**
 * @brief Variables_Functions::antiDummyHole: quasi fermi level for hole
 * @param p_dens
 * @param volt
 * @param eg
 * @param affinity
 * @param T
 * @param nv
 * @param fermip
 */
void Variables_Functions::antiDummyHole(double p_dens, double volt, double eg,
                        double affinity, double T, double kb, double e, double nv, double &fermip){
    // thermal voltage
    double thermo_volt = kb*T/e;
    double zeta;
    antiFermi(p_dens/nv,1,0,zeta);
    fermip = zeta*thermo_volt + volt + affinity + eg;
}

/**
 * @brief Variables_Functions::intrinsicFermilevel: intrinsic fermi level
 * @param volt
 * @param affinity
 * @param eg
 * @param T
 * @param nc
 * @param nv
 * @return
 */
double Variables_Functions::intrinsicFermilevel(double volt, double affinity, double eg,
                                      double T, double nc, double nv){
    double thermo_volt = kB*T/q_e;
    double fi = volt + affinity + eg/2 + (thermo_volt/2)*log(nc/nv);
    return fi;
}

/**
 * @brief Variables_Functions::gama_n_p:  gama for non-degeneration semi-conductor
 * @param eta
 * @return
 */
double Variables_Functions::gama_n_p(double eta){
    double gama;
    gama = FermiInteg(eta,0,1)/FermiInteg(eta,0,-1);
    return gama;
}

/**
 * @brief Variables_Functions::equilibrium_potential
 * @param dop
 * @param ni
 * @param nc
 * @param nv
 * @param affin
 * @param eg
 * @param T_in
 * @return
 */
double Variables_Functions::equilibrium_potential( double dop, double ni, double nc, double nv, double affin, double eg, double T_in, double kb, double e){
    double phi_eq;
    double alfa, carn, carp;
    double thermo_v = kb*T_in/e;
    if (ni == 0){
        alfa = 0;
    } else {
        alfa = dop/ni;
    }

    if (alfa > 0){
        carn = 0.5*alfa + sqrt(pow(alfa, 0.5) + 1.0);
        carp = 1/carn;
    }else {
        carp =-0.5*alfa + sqrt(pow(alfa, 0.5) + 1.0);
        carn = 1/carp;
    }

    carn = carn*ni;
    carp = carp*ni;

    phi_eq = thermo_v*asinh(dop/(2*ni)) - affin -eg/2 - (thermo_v/2)*log(nc/nv);
    return phi_eq;
}

/** *********************** Temperature dependent semi variables ******************* **/
/**
 * @brief Variables_Functions::carrierMobility
 * @param mu_o
 * @param mun_min
 * @param T_in
 * @param dop
 * @param e_field
 * @param v_sat
 * @return
 */
double Variables_Functions::carrierMobility(double mu_o, double mun_min, double T_in, double dop, double e_field, double v_sat){
    double e_ref = 2.2e7;    // electrical field limit
    double eps = 0.5;
    double beta = 1.0;
    double c_ref = 1e16;     // guessed // low field mobility
    double mu_low = carrierMobility_low_e(mu_o, mun_min, T_in, c_ref, dop);
    double mu_T;
    if (e_field < e_ref){
        mu_T = mu_low;
    }else{
        double mu_low_t = 1*mu_low;
        double mu_low_u = 0.1*mu_low;
        double v_f_t = 3e5;
        double v_f_u = 1e5;
        double mu_t = carrierMobility_high_e(mu_low_t, e_field, v_f_t, eps, beta);
        double mu_u = carrierMobility_high_e(mu_low_u, e_field, v_f_u, eps, beta);
        double PDD = carrierMobility_PDD(T_in, e_field, e_ref);
        mu_T = (mu_t + mu_u*PDD)/(1 + PDD);
    }

    return mu_T;
}

/**
 * @brief Variables_Functions::carrierMobility_low_e
 * @param mu_o
 * @param mu_min
 * @param T_in
 * @param c_ref
 * @param dop
 * @return
 */
double Variables_Functions::carrierMobility_low_e(double mu_o, double mu_min, double T_in, double c_ref, double dop){
    double T_0 = 300.0;
    double gama_0 = 1.3;
    double gama_1 =-1.5;
    double gama_2 =-0.2;
    double alfa = 0.7;

    double T_ratio = T_in/T_0;
    double c_ref_T = c_ref*pow(T_ratio, gama_0);
    double mu_o_T = mu_o*pow(T_ratio, gama_1);
    double mu_min_T = mu_min*pow(T_ratio, gama_2);
    double c_ratio = dop/c_ref_T;

    double mu_low_e = mu_min_T + (mu_o_T - mu_min_T)/(1 + pow(c_ratio, alfa));

    return mu_low_e;
}

/**
 * @brief Variables_Functions::carrierMobility_high_e
 * @param mu_low
 * @param e_field
 * @param v_sat
 * @param eps
 * @param beta
 * @return
 */
double Variables_Functions::carrierMobility_high_e(double mu_low, double e_field, double v_sat, double eps, double beta){
    double v_ratio = mu_low*e_field/v_sat;
    double mid = pow(1-eps,beta) + pow(v_ratio, beta);
    double mu_high_e = mu_low/(eps + pow(mid, 1/beta));

    return mu_high_e;
}

/**
 * @brief Variables_Functions::carrierMobility_PDD
 * @param T_in
 * @param e_field
 * @param e_ref
 * @return
 */
double Variables_Functions::carrierMobility_PDD(double T_in, double e_field, double e_ref){
    double Mu_Mt = 6.0; /** M number of equivalent vallyes **/
    double mu_mt = 1.5; /** m effective masses for low and upper valleys **/
    double delt_ec = 1.4; /** difference in conduction bands (eV) **/
    double thermo_v = kB*T_in/q_e;
    double e_ratio = e_field/e_ref;

    double PDD = Mu_Mt*pow(mu_mt, 1.5)*exp(-delt_ec/(thermo_v*(1 + e_ratio)));

    return PDD;
}

/**
 * @brief Variables_Functions::bandGap
 * @param eg_o
 * @param T_in
 * @param eg_alfa
 * @param eg_beta
 * @return
 */
double Variables_Functions::bandGap(double eg_o, double T_in,
                       double eg_alfa, double eg_beta){
    double eg;
    eg = eg_o - eg_alfa*pow(T_in,2.0)/(T_in + eg_beta);
    return eg;
}

/**
 * @brief Variables_Functions::effectivebandGap
 * @param eg_o
 * @param T_in
 * @param eg_alfa
 * @param eg_beta
 * @return
 */
double Variables_Functions::effectivebandGap(double eg_o, double T_in,
                       double eg_alfa, double eg_beta){
    double eg;
    eg = eg_o - eg_alfa*pow(T_in,2.0)/(T_in + eg_beta);
    return eg;
}

/**
 * @brief Variables_Functions::effectiveNc
 * @param T
 * @return
 */
double Variables_Functions::effectiveNc(double T){
    double me = 1.18*m_o;
    double tmp = (2*M_PI*me*kB*T/pow(2*M_PI*hbar,2.0));
    double nc = 2*pow(tmp,1.5);
    return nc;
}

/**
 * @brief Variables_Functions::effectiveNv
 * @param T
 * @return
 */
double Variables_Functions::effectiveNv(double T){
    double nv;
    double me = 0.81*m_o;
    double tmp = (2*M_PI*me*kB*T/pow(2*M_PI*hbar,2.0));
    nv = 2*pow(tmp,1.5);
    return nv;
}

/**
 * @brief Variables_Functions::effective_state_dens: effective state density
 * @param s_o
 * @param T_in
 * @return
 */
double Variables_Functions::effective_state_dens(double s_o, double T_in){
    double index = 1.5;
    double T_o = 300.0;
    double tmp_ratio = T_in/T_o;
    double dens = s_o*pow(tmp_ratio, index);
    return dens;
}

/**
 * @brief Variables_Functions::effective_mass: effective carrier mass
 * @param npo
 * @return
 */
double Variables_Functions::effective_mass(double npo){
    double me = pow(npo/2.54e19,2/3);
    return me;
}

/**
 * @brief Variables_Functions::comp_mat_param: linear composed material parameters
 * @param x
 * @param var_1
 * @param var_2
 * @return
 */
double Variables_Functions::comp_mat_param(double x, double var_1, double var_2){
    double var = var_1*x + var_2*(1-x);
    return var;
}

/**
 * @brief Variables_Functions::intrinsicCarDens:intrinsic carrier density
 *        (effective intrinsice carrier density: input effective bandgap)
 * @param nc
 * @param nv
 * @param eg
 * @param T_in
 * @return
 */
double Variables_Functions::intrinsicCarDens(double nc,
                          double nv, double eg, double T_in){
    double thermo_volt = kB*T_in/q_e;
    double ni = sqrt(nc*nv)*exp(-eg/(2*thermo_volt));
    return ni;
}

/**
 * @brief Variables_Functions::electron_affinity: electron affinity on temperature
 * @param affin_o
 * @param T_in
 * @param delt_eg
 * @param alfa
 * @param beta
 * @param gama
 * @return
 */
double Variables_Functions::electron_affinity(double affin_o, double T_in, double delt_eg,
                          double alfa, double beta, double gama){
    double affin = affin_o + gama*(alfa*pow(T_in,2.0)/(T_in + beta) + delt_eg);
    return affin;
}

/**
 * @brief Variables_Functions::thermo_emis_velo: thermo emission velocity
 * @param mx
 * @param T_in
 * @param Nx
 * @return
 */
double Variables_Functions::thermo_emis_velo(double mx, double T_in, double Nx){
    double Ax = q_e*4*M_PI*mx*m_o*pow(kB,2.0)/pow(2*M_PI*hbar,3.0);
    double vx = Ax*pow(T_in,2.0)/(q_e*Nx);
    return vx;
}

/**
 * @brief Variables_Functions::ElemVariable_T_homo
 * @param coe
 * @param T_elem
 * @param var
 */
void Variables_Functions::ElemVariable_T_homo(tbox::Array< double > coe,
                                             double T_elem,
                                             double &var){
    // size of the coe array
    int num_coe = coe.getSize();
    double tmp;
        double tmp_d;
    var = 0.0;
    for(int i=0; i<num_coe; ++i){
        tmp = i;
        tmp_d = coe[i];
        var += (tmp_d*pow(T_elem,tmp));
    }
}

/**
 * @brief Variables_Functions::ElemVariable_T_inhomo: inhomogeneous variable (T)
 * @param coe
 * @param T_elem
 * @param e_var
 */
void Variables_Functions::ElemVariable_T_inhomo(tbox::Array< double > coe,
                                               double T_elem,
                                               tbox::Vector<double> e_var){
    double tmp;
    int num_coe = coe.getSize();
    for(int i=0; i<NDIM; ++i){
        e_var[i] = 0.0;
        for(int l=0; l<num_coe; ++l){
            tmp = l;
            e_var[i] +=coe[l]*pow(T_elem,tmp);
        }
    }
}

/**
 * @brief Variables_Functions::recombGenerationSRH
 * @param carn
 * @param carp
 * @param zetan
 * @param zetap
 * @param T_in
 * @param ni
 * @param e_trap
 * @param fermi
 * @return
 */
double Variables_Functions::recombGenerationSRH(double carn, double carp, double zetan, double zetap, double T_in, double kb, double e, double ni, double e_trap, bool fermi){

    double thermov = kb*T_in/e;
    double n_o = ni*exp(e_trap/thermov);
    double p_o = ni*exp(-e_trap/thermov);

    double tao_n = 1e-7*1e12;
    double tao_p = 1e-7*1e12;

    double gaman, gamap;
    if (fermi){
        gaman = FermiInteg(zetan, 1, 1)/FermiInteg(zetan, 0, 1);
        gamap = FermiInteg(zetap, 1, 1)/FermiInteg(zetap, 0, 1);
    } else{
        gaman = 1;
        gamap = 1;
    }
    double recomG = (carn*carp-gaman*gamap*pow(ni, 2.0))/(tao_p*(carn+gaman*n_o) + tao_n*(carp+gamap*p_o));
    //陷阱辅助复合率
    return recomG;
}

/**
 * @brief Variables_Functions::recombGenerationSRH
 * @param carn
 * @param carp
 * @param ni
 * @param e_trap
 * @return
 */
adtl::AutoDScalar Variables_Functions::recombGenerationSRH(adtl::AutoDScalar carn, adtl::AutoDScalar carp,
                                                           double Vt, double Ni, double e_trap){

    double thermov = Vt;
    double n_o = Ni*exp(e_trap/thermov);
    double p_o = Ni*exp(-e_trap/thermov);

    double tao_n = 1e-7*1e12;
    double tao_p = 1e-7*1e12;

    double gaman, gamap;
    gaman = 1;
    gamap = 1;

    adtl::AutoDScalar recomG = (carn*carp-gaman*gamap*pow(Ni, 2.0))/(tao_p*(carn+gaman*n_o) + tao_n*(carp+gamap*p_o));
    adtl::AutoDScalar value = (carn*carp);
    //陷阱辅助复合率
    return recomG;
}

adtl::AutoDScalar Variables_Functions::pre_Jn(adtl::AutoDScalar Vi, adtl::AutoDScalar Vj,
                                              adtl::AutoDScalar ni, adtl::AutoDScalar nj,
                                              double mun, double Vt, double e){
    adtl::AutoDScalar pre_Jn;
    adtl::AutoDScalar deta = (Vj - Vi)/(Vt);
    pre_Jn = mun * Vt * (nj*Bern(deta)-ni*Bern(-deta));
    return pre_Jn;
}

adtl::AutoDScalar Variables_Functions::pre_Jp(adtl::AutoDScalar Vi, adtl::AutoDScalar Vj,
                                              adtl::AutoDScalar pi, adtl::AutoDScalar pj,
                                              double mup, double Vt, double e){
    adtl::AutoDScalar pre_Jp;
    adtl::AutoDScalar deta = (Vj - Vi)/(Vt);
    pre_Jp = mup * Vt * (-pj*Bern(-deta)+pi*Bern(deta));
    return pre_Jp;
}

adtl::AutoDScalar Variables_Functions::Bern(adtl::AutoDScalar delt){
    adtl::AutoDScalar data_out;
    if (fabs(delt)>1e-5){
        data_out = delt/(exp(delt)-1);
    }else{
        data_out = -0.5*delt+1.0;
        //data_out = 1.0;
    }
    return data_out;
}

adtl::AutoDScalar Variables_Functions::get_fn(adtl::AutoDScalar p_index, adtl::AutoDScalar n_index,
                                              int index,
                                              adtl::AutoDScalar V0, adtl::AutoDScalar V1,
                                              adtl::AutoDScalar V2, adtl::AutoDScalar V3,
                                              adtl::AutoDScalar V4, adtl::AutoDScalar V5,
                                              adtl::AutoDScalar V6, adtl::AutoDScalar V7,
                                              adtl::AutoDScalar n0, adtl::AutoDScalar n1,
                                              adtl::AutoDScalar n2, adtl::AutoDScalar n3,
                                              adtl::AutoDScalar n4, adtl::AutoDScalar n5,
                                              adtl::AutoDScalar n6, adtl::AutoDScalar n7,
                                              double mun, double Vt, double e,
                                              double Ni, double e_trap,
                                              tbox::Pointer<tbox::Array<hier::DoubleVector<8> > > integ,
                                              tbox::Pointer<tbox::Vector<double> > ele_volume){
    adtl::AutoDScalar fn;
    //adtl::AutoDScalar pfn = pre_Jn(V4, V5, n4, n5, mun, Vt, e);
    fn = pre_Jn(V4, V5, n4, n5, mun, Vt, e) * (*integ)[0][index]
            + pre_Jn(V7, V6, n7, n6, mun, Vt, e) * (*integ)[1][index]
            + pre_Jn(V0, V1, n0, n1, mun, Vt, e) * (*integ)[2][index]
            + pre_Jn(V3, V2, n3, n2, mun, Vt, e) * (*integ)[3][index]
            + pre_Jn(V4, V7, n4, n7, mun, Vt, e) * (*integ)[4][index]
            + pre_Jn(V0, V3, n0, n3, mun, Vt, e) * (*integ)[5][index]
            + pre_Jn(V5, V6, n5, n6, mun, Vt, e) * (*integ)[6][index]
            + pre_Jn(V1, V2, n1, n2, mun, Vt, e) * (*integ)[7][index]
            + pre_Jn(V4, V0, n4, n0, mun, Vt, e) * (*integ)[8][index]
            + pre_Jn(V5, V1, n5, n1, mun, Vt, e) * (*integ)[9][index]
            + pre_Jn(V7, V3, n7, n3, mun, Vt, e) * (*integ)[10][index]
            + pre_Jn(V6, V2, n6, n2, mun, Vt, e) * (*integ)[11][index];
            //- recombGenerationSRH(n_index, p_index, Vt, Ni, e_trap) * (*ele_volume)[index];
    #if 0
    cout<<"integ:"<<endl;
    for(int i_e=0;i_e<12;i_e++){
        cout<<(*integ)[i_e][index]<<" ";
    }
    cout<<"\nvolume:"<<(*ele_volume)[index]<<endl;
    #endif
    adtl::AutoDScalar fnl = pre_Jn(V4, V5, n4, n5, mun, Vt, e) * (*integ)[0][index]
            + pre_Jn(V7, V6, n7, n6, mun, Vt, e) * (*integ)[1][index]
            + pre_Jn(V0, V1, n0, n1, mun, Vt, e) * (*integ)[2][index]
            + pre_Jn(V3, V2, n3, n2, mun, Vt, e) * (*integ)[3][index]
            + pre_Jn(V4, V7, n4, n7, mun, Vt, e) * (*integ)[4][index]
            + pre_Jn(V0, V3, n0, n3, mun, Vt, e) * (*integ)[5][index]
            + pre_Jn(V5, V6, n5, n6, mun, Vt, e) * (*integ)[6][index]
            + pre_Jn(V1, V2, n1, n2, mun, Vt, e) * (*integ)[7][index]
            + pre_Jn(V4, V0, n4, n0, mun, Vt, e) * (*integ)[8][index]
            + pre_Jn(V5, V1, n5, n1, mun, Vt, e) * (*integ)[9][index]
            + pre_Jn(V7, V3, n7, n3, mun, Vt, e) * (*integ)[10][index]
            + pre_Jn(V6, V2, n6, n2, mun, Vt, e) * (*integ)[11][index];
    adtl::AutoDScalar fnr = recombGenerationSRH(n_index, p_index, Vt, Ni, e_trap) * (*ele_volume)[index];
    //cout<<"SRH:"<<recombGenerationSRH(n_index, p_index, Vt, Ni, e_trap)<<endl;
    return fn;
}

adtl::AutoDScalar Variables_Functions::get_fp(adtl::AutoDScalar p_index, adtl::AutoDScalar n_index,
                                              int index,
                                              adtl::AutoDScalar V0, adtl::AutoDScalar V1,
                                              adtl::AutoDScalar V2, adtl::AutoDScalar V3,
                                              adtl::AutoDScalar V4, adtl::AutoDScalar V5,
                                              adtl::AutoDScalar V6, adtl::AutoDScalar V7,
                                              adtl::AutoDScalar p0, adtl::AutoDScalar p1,
                                              adtl::AutoDScalar p2, adtl::AutoDScalar p3,
                                              adtl::AutoDScalar p4, adtl::AutoDScalar p5,
                                              adtl::AutoDScalar p6, adtl::AutoDScalar p7,
                                              double mup, double Vt, double e,
                                              double Ni, double e_trap,
                                              tbox::Pointer<tbox::Array<hier::DoubleVector<8> > > integ,
                                              tbox::Pointer<tbox::Vector<double> > ele_volume){
    adtl::AutoDScalar fp;
    fp = pre_Jp(V4, V5, p4, p5, mup, Vt, e) * (*integ)[0][index]
            + pre_Jp(V7, V6, p7, p6, mup, Vt, e) * (*integ)[1][index]
            + pre_Jp(V0, V1, p0, p1, mup, Vt, e) * (*integ)[2][index]
            + pre_Jp(V3, V2, p3, p2, mup, Vt, e) * (*integ)[3][index]
            + pre_Jp(V4, V7, p4, p7, mup, Vt, e) * (*integ)[4][index]
            + pre_Jp(V0, V3, p0, p3, mup, Vt, e) * (*integ)[5][index]
            + pre_Jp(V5, V6, p5, p6, mup, Vt, e) * (*integ)[6][index]
            + pre_Jp(V1, V2, p1, p2, mup, Vt, e) * (*integ)[7][index]
            + pre_Jp(V4, V0, p4, p0, mup, Vt, e) * (*integ)[8][index]
            + pre_Jp(V5, V1, p5, p1, mup, Vt, e) * (*integ)[9][index]
            + pre_Jp(V7, V3, p7, p3, mup, Vt, e) * (*integ)[10][index];
            + pre_Jp(V6, V2, p6, p2, mup, Vt, e) * (*integ)[11][index];
            //+ recombGenerationSRH(n_index, p_index, Vt, Ni, e_trap) * (*ele_volume)[index];
    adtl::AutoDScalar fpl = pre_Jp(V4, V5, p4, p5, mup, Vt, e) * (*integ)[0][index]
            + pre_Jp(V7, V6, p7, p6, mup, Vt, e) * (*integ)[1][index]
            + pre_Jp(V0, V1, p0, p1, mup, Vt, e) * (*integ)[2][index]
            + pre_Jp(V3, V2, p3, p2, mup, Vt, e) * (*integ)[3][index]
            + pre_Jp(V4, V7, p4, p7, mup, Vt, e) * (*integ)[4][index]
            + pre_Jp(V0, V3, p0, p3, mup, Vt, e) * (*integ)[5][index]
            + pre_Jp(V5, V6, p5, p6, mup, Vt, e) * (*integ)[6][index]
            + pre_Jp(V1, V2, p1, p2, mup, Vt, e) * (*integ)[7][index]
            + pre_Jp(V4, V0, p4, p0, mup, Vt, e) * (*integ)[8][index]
            + pre_Jp(V5, V1, p5, p1, mup, Vt, e) * (*integ)[9][index]
            + pre_Jp(V7, V3, p7, p3, mup, Vt, e) * (*integ)[10][index]
            + pre_Jp(V6, V2, p6, p2, mup, Vt, e) * (*integ)[11][index];
    adtl::AutoDScalar fpr = recombGenerationSRH(n_index, p_index, Vt, Ni, e_trap) * (*ele_volume)[index];
    return fp;
}

adtl::AutoDScalar Variables_Functions::get_fV(adtl::AutoDScalar p_index, adtl::AutoDScalar n_index,
                                              int index,
                                              adtl::AutoDScalar V0, adtl::AutoDScalar V1,
                                              adtl::AutoDScalar V2, adtl::AutoDScalar V3,
                                              adtl::AutoDScalar V4, adtl::AutoDScalar V5,
                                              adtl::AutoDScalar V6, adtl::AutoDScalar V7,
                                              double e, double eps0, double cell_epsr,
                                              tbox::Array<double> dop_node,
                                              tbox::Pointer<tbox::Matrix<double> > ele_mat_E,
                                              tbox::Pointer<tbox::Vector<double> > ele_volume){
    adtl::AutoDScalar fV;
    fV = cell_epsr*(*ele_mat_E)(index,0) * V0
         + cell_epsr*(*ele_mat_E)(index,1) * V1
         + cell_epsr*(*ele_mat_E)(index,2) * V2
         + cell_epsr*(*ele_mat_E)(index,3) * V3
         + cell_epsr*(*ele_mat_E)(index,4) * V4
         + cell_epsr*(*ele_mat_E)(index,5) * V5
         + cell_epsr*(*ele_mat_E)(index,6) * V6
         + cell_epsr*(*ele_mat_E)(index,7) * V7
        - (e/eps0)*(n_index - p_index - dop_node[index]) * (*ele_volume)[index];
    #if 0
    for(int i=0; i<8; i++){
        cout<<(*ele_volume)[index]<<" ";
    }
    cout<<"\n";
    for(int i=0;i<8;i++){
        cout<<(*ele_mat_E)(index,i)<<" ";
    }
    #endif
    adtl::AutoDScalar fr = -(e/eps0)*(n_index - p_index - dop_node[index]) * (*ele_volume)[index];
    adtl::AutoDScalar fl = cell_epsr*(*ele_mat_E)(index,0) * V0
            + cell_epsr*(*ele_mat_E)(index,1) * V1
            + cell_epsr*(*ele_mat_E)(index,2) * V2
            + cell_epsr*(*ele_mat_E)(index,3) * V3
            + cell_epsr*(*ele_mat_E)(index,4) * V4
            + cell_epsr*(*ele_mat_E)(index,5) * V5
            + cell_epsr*(*ele_mat_E)(index,6) * V6
            + cell_epsr*(*ele_mat_E)(index,7) * V7;
    return fV;
}

adtl::AutoDScalar Variables_Functions::get_fV(double a, double b, double u_t, double dt,
                                              tbox::Array<adtl::AutoDScalar> V_node,
                                              adtl::AutoDScalar h,
                                              double ele_volume,
                                              int index,
                                              tbox::Pointer<tbox::Matrix<double> > ele_mat_E){
    adtl::AutoDScalar fV=0;
    fV = (a*(V_node[index]-u_t)/dt+h) * (ele_volume);
    int num_ver = V_node.getSize();
    for(int j=0;j<num_ver;j++){
        fV = fV + b*(*ele_mat_E)(index,j)*V_node[j];
    }
    #if 0
    for(int j=0;j<num_ver;j++){
        cout<<"get_fV,V_node:"<<V_node[j].getValue()<<endl;
    }
    #endif
    return fV;
}

adtl::AutoDScalar Variables_Functions::get_fV(double a, double b, double u_t, double dt,
                                              tbox::Array<double> V_node,
                                              adtl::AutoDScalar h,
                                              double ele_volume,
                                              int index,
                                              tbox::Pointer<tbox::Matrix<double> > ele_mat_E){
    adtl::AutoDScalar fV=0;
    adtl::AutoDScalar V_node_index = V_node[index]; V_node_index.setADValue(index*3,1.0);
    fV = (a*(V_node_index-u_t)/dt+h) * (ele_volume);
    int num_ver = V_node.getSize();
    for(int j=0;j<num_ver;j++){
        adtl::AutoDScalar V_node_j = V_node[j]; V_node_j.setADValue(j*3,1.0);
        fV = fV + b*(*ele_mat_E)(index,j)*V_node_j;
    }
    #if 0
    for(int j=0;j<num_ver;j++){
        cout<<"get_fV,V_node:"<<V_node[j].getValue()<<endl;
    }
    #endif
    return fV;
}

adtl::AutoDScalar Variables_Functions::get_fn(double a, double u_t, double dt,
                                              tbox::Array<adtl::AutoDScalar> n_node,
                                              tbox::Array<adtl::AutoDScalar> V_node,
                                              adtl::AutoDScalar h,
                                              double ele_volume,
                                              int index,
                                              tbox::Array<tbox::Array<int> > r_node,
                                              tbox::Array<tbox::Array<double> > r_integ,
                                              double mun, double Vt, double e){
    adtl::AutoDScalar fn =0;
    fn = (a*(n_node[index]-u_t)/dt+h) * (ele_volume);
    double num_edge = r_node.getSize();
    for(int i_e = 0; i_e < num_edge; ++i_e){
        int n_1 = r_node[i_e][0];
        int n_2 = r_node[i_e][1];
        adtl::AutoDScalar sub_Jn=
                pre_Jn(V_node[n_1], V_node[n_2], n_node[n_1], n_node[n_2], mun, Vt, e) * r_integ[i_e][index];
        fn = fn + sub_Jn;
    }
#if 0
for(int j=0;j<8;j++){
    cout<<"get_fn,n_node:"<<n_node[j].getValue()<<endl;
}
#endif
    return fn;
}

adtl::AutoDScalar Variables_Functions::get_fn(double a, double u_t, double dt,
                                              tbox::Array<double> n_node,
                                              tbox::Array<double> V_node,
                                              adtl::AutoDScalar h,
                                              double ele_volume,
                                              int index,
                                              tbox::Array<tbox::Array<int> > r_node,
                                              tbox::Array<tbox::Array<double> > r_integ,
                                              double mun, double Vt, double e){
    adtl::AutoDScalar fn =0;
    adtl::AutoDScalar n_node_index = n_node[index]; n_node_index.setADValue(index*3+1,1.0);
    fn = (a*(n_node_index-u_t)/dt+h) * (ele_volume);
    double num_edge = r_node.getSize();
    for(int i_e = 0; i_e < num_edge; ++i_e){
        int n_1 = r_node[i_e][0];
        int n_2 = r_node[i_e][1];
        adtl::AutoDScalar n_node_n_1 = n_node[n_1]; n_node_n_1.setADValue(n_1*3+1,1.0);
        adtl::AutoDScalar n_node_n_2 = n_node[n_2]; n_node_n_2.setADValue(n_2*3+1,1.0);
        adtl::AutoDScalar V_node_n_1 = V_node[n_1]; V_node_n_1.setADValue(n_1*3,1.0);
        adtl::AutoDScalar V_node_n_2 = V_node[n_2]; V_node_n_2.setADValue(n_2*3,1.0);
        adtl::AutoDScalar sub_Jn=
                pre_Jn(V_node_n_1, V_node_n_2, n_node_n_1, n_node_n_2, mun, Vt, e) * r_integ[i_e][index];
        fn = fn + sub_Jn;
    }
#if 0
for(int j=0;j<8;j++){
    cout<<"get_fn,n_node:"<<n_node[j].getValue()<<endl;
}
#endif
    return fn;
}

adtl::AutoDScalar Variables_Functions::get_fp(double a, double u_t, double dt,
                                              tbox::Array<adtl::AutoDScalar> p_node,
                                              tbox::Array<adtl::AutoDScalar> V_node,
                                              adtl::AutoDScalar h,
                                              double ele_volume,
                                              int index,
                                              tbox::Array<tbox::Array<int> > r_node,
                                              tbox::Array<tbox::Array<double> > r_integ,
                                              double mup, double Vt, double e){
    adtl::AutoDScalar fp =0;
    fp = (a*(p_node[index]-u_t)/dt+h) * (ele_volume);
    double num_edge = r_node.getSize();
    for(int i_e = 0; i_e < num_edge; ++i_e){
        int n_1 = r_node[i_e][0];
        int n_2 = r_node[i_e][1];
        adtl::AutoDScalar sub_Jp=
                pre_Jp(V_node[n_1], V_node[n_2], p_node[n_1], p_node[n_2], mup, Vt, e) * r_integ[i_e][index];
        fp = fp + sub_Jp;
    }
    return fp;
}

adtl::AutoDScalar Variables_Functions::get_fp(double a, double u_t, double dt,
                                              tbox::Array<double> p_node,
                                              tbox::Array<double> V_node,
                                              adtl::AutoDScalar h,
                                              double ele_volume,
                                              int index,
                                              tbox::Array<tbox::Array<int> > r_node,
                                              tbox::Array<tbox::Array<double> > r_integ,
                                              double mup, double Vt, double e){
    adtl::AutoDScalar fp =0;
    adtl::AutoDScalar p_node_index = p_node[index]; p_node_index.setADValue(index*3+2,1.0);
    fp = (a*(p_node_index-u_t)/dt+h) * (ele_volume);
    double num_edge = r_node.getSize();
    for(int i_e = 0; i_e < num_edge; ++i_e){
        int n_1 = r_node[i_e][0];
        int n_2 = r_node[i_e][1];
        adtl::AutoDScalar p_node_n_1 = p_node[n_1]; p_node_n_1.setADValue(n_1*3+2,1.0);
        adtl::AutoDScalar p_node_n_2 = p_node[n_2]; p_node_n_2.setADValue(n_2*3+2,1.0);
        adtl::AutoDScalar V_node_n_1 = V_node[n_1]; V_node_n_1.setADValue(n_1*3,1.0);
        adtl::AutoDScalar V_node_n_2 = V_node[n_2]; V_node_n_2.setADValue(n_2*3,1.0);
        adtl::AutoDScalar sub_Jp=
                pre_Jp(V_node_n_1, V_node_n_2, p_node_n_1, p_node_n_2, mup, Vt, e) * r_integ[i_e][index];
        fp = fp + sub_Jp;
    }
    return fp;
}

adtl::AutoDScalar Variables_Functions::get_fVogamma(double a, double u_t, double dt,
                                                    double b,
                                                    adtl::AutoDScalar Vogamma_node,
                                                    adtl::AutoDScalar h,
                                                    double ele_volume){
    adtl::AutoDScalar fVogamma=0;
    fVogamma = (a*(Vogamma_node-u_t)/dt+b*Vogamma_node+h) * (ele_volume);

    #if 0
    for(int i=0; i<8; i++){
        cout<<(*ele_volume)[index]<<" ";
    }
    #endif
    return fVogamma;
}

/** *************************************************************************************************** **/
/**
 * @brief Variables_Functions::PecletNum
 * @param e_field
 * @param mup
 * @param Dp
 * @param elem_size
 * @return
 */
double Variables_Functions::PecletNum(double e_field, double mup,
                          double Dp, double elem_size){
    double Pe = mup*e_field*elem_size/(2*Dp);
    return Pe;

}

/**
 * @brief Variables_Functions::Bern
 * @param delt
 * @return
 */
double Variables_Functions::Bern(double delt){
    double data_out;
    if (fabs(delt)>1e-5){
        data_out = delt/(exp(delt)-1);
    }else{
        data_out = 1.0;
        //data_out = -0.5*delt+1.0;
    }
    return data_out;
}

/**
 * @brief Variables_Functions::face_CV_integ
 * @param vertex
 * @return
 */
tbox::Array<double> Variables_Functions::face_CV_integ(tbox::Array<hier::DoubleVector<NDIM> > vertex){
    tbox::Array<double> integ;
    int n_vertex = vertex.getSize();
    int n_edge = 4;
    tbox::Array<hier::DoubleVector<2> > r_node;
    tbox::Array<hier::DoubleVector<2> > r_edge;
    r_node.resizeArray(n_vertex);
    r_edge.resizeArray(n_edge);
    integ.resizeArray(n_vertex);

    r_node[0][0] = 4; r_node[0][1] = 5;
    r_node[1][0] = 7; r_node[1][1] = 6;
    r_node[2][0] = 0; r_node[2][1] = 1;
    r_node[3][0] = 3; r_node[3][1] = 2;
    r_node[4][0] = 4; r_node[4][1] = 7;
    r_node[5][0] = 0; r_node[5][1] = 3;
    r_node[6][0] = 5; r_node[6][1] = 6;
    r_node[7][0] = 1; r_node[7][1] = 2;
    r_node[8][0] = 4; r_node[8][1] = 0;
    r_node[9][0] = 5; r_node[9][1] = 1;
    r_node[10][0] = 7; r_node[10][1] = 3;
    r_node[11][0] = 6; r_node[11][1] = 2;

    r_edge[0][0] = 2; r_edge[0][1] = 0;
    r_edge[1][0] = 0; r_edge[1][1] = 1;
    r_edge[2][0] = 1; r_edge[2][1] = 2;
    r_edge[3][0] = 2; r_edge[3][1] = 0;
    r_edge[4][0] = 0; r_edge[4][1] = 1;
    r_edge[5][0] = 1; r_edge[5][1] = 2;
    r_edge[6][0] = 2; r_edge[6][1] = 0;
    r_edge[7][0] = 0; r_edge[7][1] = 1;
    r_edge[8][0] = 1; r_edge[8][1] = 2;
    r_edge[9][0] = 2; r_edge[9][1] = 0;
    r_edge[10][0] = 0; r_edge[10][1] = 1;
    r_edge[11][0] = 1; r_edge[11][1] = 2;

    tbox::Array<hier::DoubleVector<NDIM> > e_w_crd(n_vertex);
    tbox::Array<hier::DoubleVector<NDIM> > f_w_crd(1);

    for (int i_d = 0; i_d < NDIM; ++i_d){
        f_w_crd[0][i_d] = 0.0;
        for (int i_edge = 0; i_edge < 4; ++i_edge){
            e_w_crd[i_edge][i_d] = 0.0;
        }
    }

    for (int i_d = 0; i_d < NDIM; ++i_d){
        for (int i_n = 0; i_n < n_vertex; ++i_n){
            f_w_crd[0][i_d] += vertex[i_n][i_d]/4;
        }
        for (int i_edge = 0; i_edge < 4; ++i_edge){
            e_w_crd[i_edge][i_d] = (vertex[r_node[i_edge][0]][i_d] + vertex[r_node[i_edge][1]][i_d])/2;
        }
    }

    for (int i_v = 0; i_v < n_vertex; ++i_v){
        tbox::Array<hier::DoubleVector<NDIM> > tmp_crd(4);
        tmp_crd[0] = vertex[i_v];
        tmp_crd[1] = e_w_crd[r_edge[i_v][1]];
        tmp_crd[2] = f_w_crd[0];
        tmp_crd[3] = e_w_crd[r_edge[i_v][0]];
        integ[i_v] = quadrilateral_area(tmp_crd);
    }
    return integ;
}

/**
 * @brief Variables_Functions::quadrilateral_area
 * @param crd
 * @return
 */
double Variables_Functions::quadrilateral_area(tbox::Array<hier::DoubleVector<NDIM> > crd){
    hier::DoubleVector<NDIM> vec_1, vec_2, vec_norm;
    for (int i_d = 0; i_d < NDIM; ++i_d){
        vec_1[i_d] = crd[2][i_d] - crd[0][i_d];
        vec_2[i_d] = crd[3][i_d] - crd[1][i_d];
    }
    vec_norm = vector_product(vec_1, vec_2);
    double area = sqrt(pow(vec_norm[0],2.0) + pow(vec_norm[1],2.0) + pow(vec_norm[2],2.0));
    return area;
}
/** *************************************************************************************************** **/
/**
 * @brief Variables_Functions::vecNorm
 * @param vec
 * @return
 */
double Variables_Functions::vecNorm(tbox::Array<double> vec){
    int num_e = vec.getSize();
    double sum = 0;
    for (int i=0; i<num_e; ++i){
        sum += vec[i]*vec[i];
    }
    double vec_norm = sqrt(sum);
    return vec_norm;
}

double Variables_Functions::vecNorm(tbox::Pointer<tbox::Array<double> > vec){
    int num_e = vec->getSize();
    double sum = 0.0;
    for (int i=0; i<num_e; ++i){
        sum += (*vec)[i]*(*vec)[i];
    }
    double vec_norm = sqrt(sum);
    return vec_norm;
}

/**
 * @brief Variables_Functions::vector_product
 * @param vec_1
 * @param vec_2
 * @return
 */
hier::DoubleVector<NDIM> Variables_Functions::vector_product(hier::DoubleVector<NDIM> vec_1,
                                                             hier::DoubleVector<NDIM> vec_2){
    hier::DoubleVector<NDIM> vec_norm;
    vec_norm[0] = vec_1[1]*vec_2[2] - vec_1[2]*vec_2[1];
    vec_norm[1] = vec_1[2]*vec_2[0] - vec_1[0]*vec_2[2];
    vec_norm[2] = vec_1[0]*vec_2[1] - vec_1[1]*vec_2[0];
    return vec_norm;
}

/**
 * @brief Variables_Functions::node2ElemVariable: interpolate node variables into element
 * @param real_vertex
 * @param var
 * @param dt
 * @param time
 * @param e_var
 */
double Variables_Functions::node2ElemVariable(tbox::Array< hier::DoubleVector<NDIM> > real_vertex,
                                      tbox::Array< double > var,
                                      const double dt,
                                      const double time){

    int num_var=var.getSize();

    tbox::Pointer<IntegratorManager<NDIM> > integrator_manager
      = IntegratorManager<NDIM>::getManager();
    tbox::Pointer<BaseIntegrator<NDIM> > integrator;

    tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager
      = ShapeFunctionManager<NDIM>::getManager();
    tbox::Pointer<BaseShapeFunction<NDIM> > shape_func;

    if(num_var==8){
        integrator   = integrator_manager->getIntegrator("Hexahedron");
        shape_func = shape_manager->getShapeFunction("Hexahedron");
    }else if(num_var==4){
        integrator   = integrator_manager->getIntegrator("Tetrahedron");
        shape_func = shape_manager->getShapeFunction("Tetrahedron");
    }else{
        TBOX_ERROR( num_var << ": "
                   << " No matched number of node value in node2ElemVariable." << endl);
    }

    int n_dof = shape_func->getNumberOfDof();
    int num_quad_pnts = integrator->getNumberOfQuadraturePoints();
    tbox::Array<hier::DoubleVector<NDIM> >
      quad_pnt = integrator->getQuadraturePoints(real_vertex);

    tbox::Array<double> weight = integrator->getQuadratureWeights();

    tbox::Array<tbox::Array<double> >
      bas_val  = shape_func->value(real_vertex,quad_pnt);

    double e_var = 0.0;
    for (int i = 0; i < n_dof; ++i) {
        for (int l = 0; l < num_quad_pnts; ++l) {
            e_var += var[i]*weight[l]*bas_val[l][i];
            #if 0
            if(var[i]!=300){
                cout<<"温度不是300K"<<endl;
            }else{
                //cout<<"温度是300K"<<endl;
            }
            #endif
        }
    }
    return e_var;
}

/**
 * @brief Variables_Functions::debyeLength
 * @param L
 * @param eps
 * @param n_o
 * @param debye_len
 */
void Variables_Functions::debyeLength(double L, double eps,
                                   double n_o, double &debye_len){
    double thermal_volt = T_amb*kB/q_e; //thermal voltage
    debye_len = (thermal_volt*eps*eps_o)/(q_e*L*L*n_o);
}

/**
 * @brief Variables_Functions::DataOutput
 * @param filename
 * @param data
 */
void Variables_Functions::DataOutput(string filename,
                                     tbox::Pointer<pdat::VectorData<NDIM, double> > data){
    int num_e = data->getVectorLength();
    ofstream fout;
    fout.open(filename.data());
    for (int i=0;i<num_e;++i){
        fout<<(*data)(i)<<endl;
    }
    fout.close();
}

/**
 * @brief Variables_Functions::Constraint
 * @param data
 */
double Variables_Functions::Constraint(double data, double constraint_value){
    if(data<=0){
        data = constraint_value;
        //cout<<"nagetive value for concentration!"<<endl;
    }
    //cout<<"concentration:"<<data<<endl;
    return data;
}

void Variables_Functions::calculateDop(tbox::Pointer<pdat::NodeData<NDIM, double> > node_coord,
                                       int num_nodes,
                                       tbox::Pointer<pdat::NodeData<NDIM, double> > ND_data,
                                       tbox::Pointer<pdat::NodeData<NDIM, double> > NA_data,
                                       tbox::Array<doplist> Dopassignment){
    for(int i=0; i<num_nodes; i++){
        double x = (*node_coord)(0, i);
        double y = (*node_coord)(1, i);
        double z = (*node_coord)(2, i);
        double NDprev = 0;
        double NAprev = 0;
        // p doping
        double NA0 = Dopassignment[0].NA0;
        double NA = NAprev + NA0;
        double ND = NDprev;
        // n doping 1
        ND = ND + getND(Dopassignment,1, NA0, x, y, z);
        // n doping 2
        ND = ND + getND(Dopassignment,2, NA0, x, y, z);
        // assignment the ND and NA data
        ND_data->getPointer()[i] = ND;
        NA_data->getPointer()[i] = NA;
    }
}

double Variables_Functions::getND(tbox::Array<doplist> Dopassignment,int index,
                                  double NA0, double x, double y, double z){
    double ND01 = Dopassignment[index].ND0;
    double r0x1,r0y1,r0z1; double W1,D1,H1;
    double rxminus1,rxplus1; double ryminus1,ryplus1; double rzminus1,rzplus1;
    double djx1,djy1,djz1;
    double k = 1;
    r0x1 = Dopassignment[index].r0[0];
    r0y1 = Dopassignment[index].r0[1];
    r0z1 = Dopassignment[index].r0[2];
    W1 = Dopassignment[index].W;
    D1 = Dopassignment[index].D;
    H1 = Dopassignment[index].H;
    rxminus1 = r0x1;
    rxplus1 = r0x1+W1;
    ryminus1 = r0y1;
    ryplus1 = r0y1+D1;
    rzminus1 = r0z1;
    rzplus1 = r0z1+H1;
    djx1 = Dopassignment[index].dj[0];
    djy1 = Dopassignment[index].dj[1];
    djz1 = Dopassignment[index].dj[2];
    double Nb1 = NA0;
    double lx1,ly1,lz1;
    lx1 = djx1/sqrt(log(ND01/Nb1));
    ly1 = djy1/sqrt(log(ND01/Nb1));
    lz1 = djz1/sqrt(log(ND01/Nb1));
    double rxminusreal,ryminusreal,rzminusreal; double rxplusreal,ryplusreal,rzplusreal;
    if(x>=rxminus1) rxminusreal = 0;
    else rxminusreal = k*(rxminus1 - x);
    if(x<=rxplus1) rxplusreal = 0;
    else rxplusreal = k*(x - rxplus1);
    if(y>=ryminus1) ryminusreal = 0;
    else ryminusreal = k*(ryminus1 - y);
    if(y<=ryplus1) ryplusreal = 0;
    else ryplusreal = k*(y - ryplus1);
    if(z>=rzminus1) rzminusreal = 0;
    else rzminusreal = k*(rzminus1 - z);
    if(z<=rzplus1) rzplusreal = 0;
    else rzplusreal = k*(z - rzplus1);
    double ND1 = ND01*exp(-(pow(rxplusreal/lx1,2)+pow(rxminusreal/lx1,2)+
                         pow(ryplusreal/ly1,2)+pow(ryminusreal/ly1,2)+
                         pow(rzplusreal/lz1,2)+pow(rzminusreal/lz1,2)));
    return ND1;
}

