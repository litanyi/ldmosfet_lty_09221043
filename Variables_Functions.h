/*
 * Variables_Functions.h
 *
 *  Created on: May 10, 2017
 *      Author: tong
 */

#ifndef SOURCE_DIRECTORY__VARIABLES_FUNCTIONS_H_
#define SOURCE_DIRECTORY__VARIABLES_FUNCTIONS_H_

#include <string>
#include "math.h"
#include "../interface/BaseElement.h"
#include "../interface/BaseIntegrator.h"
#include "../interface/BaseShapeFunction.h"
#include "DoubleVector.h"
//#include "JVector.h"
#include "VectorVariable.h"
#include "VectorData.h"
#include "Vector.h"

#include "IntegratorManager.h"
#include "ShapeFunctionManager.h"
#include "ElementManager.h"
#include "Bound_condition.h"

#include "material_database.h"
#include "Global_variable.h"
#include "adolc.h"

using namespace std;
using namespace JAUMIN;

const double T_amb = 300.0;

class Variables_Functions{
public:
    /**
     * @brief construction function
     */
    Variables_Functions();

    /**
     * @brief anti construction function
     */
    ~Variables_Functions();

    /**
    *  @param coord  input, Pointer, point to vertex coord
    *  @param var    input, Pointer, point to vertex variables
    *  @param dt     input, double, time step
    *  @param time   input, double, current time
    */
    double node2ElemVariable(tbox::Array< hier::DoubleVector<NDIM> > coord,
                                   tbox::Array<double> var,
                                   const double dt,
                                   const double time);

    /**
    *  @param coord  input, Pointer, point to vertex coord
    *  @param e_var    input, Pointer, point to vertex variables
    *  @param dt     input, double, time step
    *  @param time   input, double, current time
    *  @param var  output, interpolating for effective variable on the whole elem
    */
    void elem2NodeVariable(tbox::Array< hier::DoubleVector<NDIM> > coord,
                                   tbox::Array<double> e_var,
                                   const double dt,
                                   const double time,
                                   double & var);

    /**
     * @param N    input, Pointer, node carrier density
     * @param n_o  output, calculated carrier density scalor
     */
    double effectiveCarScalor(tbox::Array<material_pro> mat_list);

    /**
     * @param T    input, double, element temperature
     * @param L    input, length
     * @param eps  input, relative permittivity
     * @param n_o  input, carrier density scalor
     * @param debye_len output, calculated debye length
     */
    void debyeLength(double L, double eps, double n_o, double &debye_len);

    /**
     * @brief carrier density <-> quasi-fermi level
     * @param ef           input, double, the normalized relative fermi level position
     * @param fermi_flag   input, int, the degeneracy flag: 1 for degenerate calculation;
     *                     0 for non-degenerate calculation
     * @param fermi_order  input, int, the integral order: 1, 0, -1
     * @param fermii       output, double, the calculated result
     */
    double FermiInteg(double ef, int fermi_flag, int fermi_order);

    /**
     * @brief anti-fermi integral
     * @param ef           output, double
     * @param fermi_flag   input, int, the degeneracy flag: 1 for degenerate calculation;
     *                     0 for non-degenerate calculation
     * @param fermi_order  input, int, the integral order
     * @param fermi        input, double, exponential part of carrier density
     */
    void antiFermi(double fermi, int fermi_flag, int fermi_order, double &ef);


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
    adtl::AutoDScalar dummyElectron(adtl::AutoDScalar volt, double fermin, bool deriv_flag,
                         double affinity, double T, double kb, double e, double nc);

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
    adtl::AutoDScalar dummyHole(adtl::AutoDScalar volt, double fermip, double eg, bool deriv_flag,
                           double affinity, double T, double kb, double e, double nv);

    /**
     * @brief Maxwell-Boltzmann (MB) statistic for electron and hole
     * @param volt         input, the potential
     * @param fermin       input, the electron fermi level
     * @param fermip       input, the hole quasi-fermi level
     * @param eg           input, the band gap
     * @param afffinity    input, the electron affinity
     * @param T_node       input, the temperature for node
     * @param n_dens       output, the calculated electron density
     * @param p_dens       output, the calculated hole density
     */
    void dummyElectron(double volt, double fermin, bool deriv_flag,
                       double affinity, double T, double kb, double e, double nc, double &n_dens);

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
    double dummyElectron(double volt, double fermin, bool deriv_flag,
                         double affinity, double T, double kb, double e, double nc);

    void dummyHole(double volt, double fermip, double eg, bool deriv_flag,
                       double affinity, double T, double kb, double e, double nv, double &p_dens);

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
    double dummyHole(double volt, double fermip, double eg, bool deriv_flag,
                           double affinity, double T, double kb, double e, double nv);


    /**
     * @brief anti-dummy process, used to convert the Ne to quasi fermi energy
     * @param n_dens      input, double, the electron density
     * @param p_dens      input, double, the hole density
     * @param volt        input, double, electrical potential
     * @param eg          input, double, bandgap (temperature dependent parameter)
     * @param affinity    input, double, electron affinity (temperature dependent parameter)
     * @param T_node      input, double, temperature at nodes
     * @param fermin/p    output, double, the quasi-fermi level
     */
    void antiDummyElectron(double n_dens, double volt,
                       double affinity, double T, double kb, double e, double nc, double &fermin);

    void antiDummyHole(double p_dens, double volt, double eg, double T, double kb, double e,
                       double affinity, double nv, double &fermip);

    /**
     * @brief intrinsicFermilevel
     * @param volt
     * @param affinity
     * @param eg
     * @param T
     * @param nc
     * @param nv
     * @return
     */
    double intrinsicFermilevel(double volt, double affinity, double eg,
                               double T, double nc, double nv);

    /**
     * @brief equilibrium_potential
     * @param dop
     * @param ni
     * @param nc
     * @param nv
     * @param affin
     * @param eg
     * @param T_in
     * @return
     */
    double equilibrium_potential( double dop, double ni, double nc, double nv,
                                  double affin, double eg, double T_in, double kb, double e);

    /**
     * @brief gama_n_p
     * @param eta
     * @return
     */
    double gama_n_p(double eta);

    /** **********************semi-conduction variable models********************* **/
    /**
     * @brief carrierMobility
     * @param mu_o
     * @param mun_min
     * @param T_in
     * @param dop
     * @param e_field
     * @param v_sat
     * @return
     */
    double carrierMobility(double mu_o, double mun_min, double T_in,
                           double dop, double e_field, double v_sat);

    /**
     * @brief carrierMobility_low_e
     * @param mu_o
     * @param mun_min
     * @param T_in
     * @param c_ref
     * @param dop
     * @return
     */
    double carrierMobility_low_e(double mu_o, double mun_min, double T_in,
                                 double c_ref, double dop);

    /**
     * @brief carrierMobility_high_e
     * @param mu_low
     * @param e_field
     * @param v_sat
     * @param eps
     * @param beta
     * @return
     */
    double carrierMobility_high_e(double mu_low, double e_field, double v_sat,
                                  double eps, double beta);

    /**
     * @brief carrierMobility_PDD
     * @param T_in
     * @param e_field
     * @param e_ref
     * @return
     */
    double carrierMobility_PDD(double T_in, double e_field, double e_ref);

    /**
     * @brief recombGenerationSRH : Recombination-Generation model
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
    double recombGenerationSRH(double carn, double carp, double zetan, double zetap,
                               double T_in, double kb, double e, double ni, double e_trap, bool fermi);

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
    adtl::AutoDScalar recombGenerationSRH(adtl::AutoDScalar carn, adtl::AutoDScalar carp,
                                          double Vt, double ni, double e_trap);
    /**
     * @brief Variables_Functions::recombGenerationSRH
     * @param carn
     * @param carp
     * @param ni
     * @param e_trap
     * @return
     */
    double recombGenerationSRH(double carn, double carp,
                               double Vt, double Ni,
                               double e_trap);

    /**
     * @brief bandGap
     * @param eg_o
     * @param T_in
     * @param eg_alfa
     * @param eg_beta
     * @return
     */
    double bandGap(double eg_o, double T_in, double eg_alfa, double eg_beta);

    /**
     * @brief effectivebandGap
     * @param eg_o
     * @param T_in
     * @param eg_alfa
     * @param eg_beta
     * @return
     */
    double effectivebandGap(double eg_o, double T_in, double eg_alfa, double eg_beta);

    /**
     * @brief effectiveNc
     * @param T
     * @return
     */
    double effectiveNc(double T);

    /**
     * @brief effectiveNv
     * @param T
     * @return
     */
    double effectiveNv(double T);

    /**
     * @brief effective_state_dens
     * @param s_o
     * @param T_in
     * @return
     */
    double effective_state_dens(double s_o, double T_in);

    /**
     * @brief effective_mass
     * @param npo
     * @return
     */
    double effective_mass(double npo);

    /**
     * @brief comp_mat_param
     * @param x
     * @param var_1
     * @param var_2
     * @return
     */
    double comp_mat_param(double x, double var_1, double var_2);

    /**
     * @brief thermo_emis_velo
     * @param mx
     * @param T_in
     * @param Nx
     * @return
     */
    double thermo_emis_velo(double mx, double T_in, double Nx);

    /**
     * @brief intrinsicCarDens
     * @param nc
     * @param nv
     * @param eg
     * @param T_in
     * @return
     */
    double intrinsicCarDens(double nc, double nv, double eg, double T_in);

    /**
     * @brief electron_affinity
     * @param affin_o
     * @param T_in
     * @param delt_eg
     * @param alfa
     * @param beta
     * @param gama
     * @return
     */
    double electron_affinity(double affin_o, double T_in, double delt_eg,
                             double alfa, double beta, double gama);

    /** ***********************thermal dependent material variables******************* **/
    /**
     * @brief linear polynomial interpolation for homogeneous variables
     * @param coe     input, Pointer, point to thermal interpolating coefficient
     * @param T_elem  input, double, temperature
     * @param var     output, double, variable at T
     */
    void ElemVariable_T_homo(tbox::Array< double > coe,
                            double T_elem, double &e_var);

    /**
     * @brief linear polynomial interpolation for in-homogeneous variables
     * @param coe       input, Pointer, point to thermal interpolating coefficient
     * @param T         input, double, temperature
     * @param e_var     output, double, variable at T
     */
    void ElemVariable_T_inhomo(tbox::Array< double >  coe,
                               double T_elem,
                               tbox::Vector<double> e_var);

    /**
     * @brief vecNorm
     * @param e_field
     * @return
     */
    double vecNorm(tbox::Array<double> e_field);

    double vecNorm(tbox::Pointer<tbox::Array<double> > e_field);

    /**
     * @brief vector_product
     * @param vec_1
     * @param vec_2
     * @return
     */
    hier::DoubleVector<NDIM> vector_product(hier::DoubleVector<NDIM> vec_1,
                                            hier::DoubleVector<NDIM> vec_2);

    /**
     * @brief PecletNum
     * @param e_field
     * @param mup
     * @param Dp
     * @param elem_size
     * @return
     */
    double PecletNum(double e_field, double mup, double Dp, double elem_size);

    /**
     * @brief Bern
     * @param delt
     * @return
     */
    double Bern(double delt);

    adtl::AutoDScalar pre_Jn(adtl::AutoDScalar Vi, adtl::AutoDScalar Vj,
                             adtl::AutoDScalar ni, adtl::AutoDScalar nj,
                             double mun, double Vt, double e);

    adtl::AutoDScalar pre_Jp(adtl::AutoDScalar Vi, adtl::AutoDScalar Vj,
                             adtl::AutoDScalar pi, adtl::AutoDScalar pj,
                             double mup, double Vt, double e);

    adtl::AutoDScalar Bern(adtl::AutoDScalar delt);

    adtl::AutoDScalar get_fn(adtl::AutoDScalar p_index, adtl::AutoDScalar n_index,
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
                             double ni, double e_trap,
                             tbox::Pointer<tbox::Array<hier::DoubleVector<8> > > integ,
                             tbox::Pointer<tbox::Vector<double> > ele_volume);

    adtl::AutoDScalar get_fp(adtl::AutoDScalar p_index, adtl::AutoDScalar n_index,
                             int index,
                             adtl::AutoDScalar V0, adtl::AutoDScalar V1,
                             adtl::AutoDScalar V2, adtl::AutoDScalar V3,
                             adtl::AutoDScalar V4, adtl::AutoDScalar V5,
                             adtl::AutoDScalar V6, adtl::AutoDScalar V7,
                             adtl::AutoDScalar n0, adtl::AutoDScalar n1,
                             adtl::AutoDScalar n2, adtl::AutoDScalar n3,
                             adtl::AutoDScalar n4, adtl::AutoDScalar n5,
                             adtl::AutoDScalar n6, adtl::AutoDScalar n7,
                             double mup, double Vt, double e,
                             double ni, double e_trap,
                             tbox::Pointer<tbox::Array<hier::DoubleVector<8> > > integ,
                             tbox::Pointer<tbox::Vector<double> > ele_volume);

    adtl::AutoDScalar get_fn(double a, double u_t, double dt,
                                                  tbox::Array<adtl::AutoDScalar> n_node,
                                                  tbox::Array<adtl::AutoDScalar> V_node,
                                                  adtl::AutoDScalar h,
                                                  double ele_volume,
                                                  int index,
                                                  tbox::Array<tbox::Array<int> > r_node,
                                                  tbox::Array<tbox::Array<double> > r_integ,
                                                  double mun, double Vt, double e);

    adtl::AutoDScalar get_fp(double a, double u_t, double dt,
                                                  tbox::Array<adtl::AutoDScalar> p_node,
                                                  tbox::Array<adtl::AutoDScalar> V_node,
                                                  adtl::AutoDScalar h,
                                                  double ele_volume,
                                                  int index,
                                                  tbox::Array<tbox::Array<int> > r_node,
                                                  tbox::Array<tbox::Array<double> > r_integ,
                                                  double mup, double Vt, double e);

    adtl::AutoDScalar get_fV(adtl::AutoDScalar p_index, adtl::AutoDScalar n_index,
                             int index,
                             adtl::AutoDScalar V0, adtl::AutoDScalar V1,
                             adtl::AutoDScalar V2, adtl::AutoDScalar V3,
                             adtl::AutoDScalar V4, adtl::AutoDScalar V5,
                             adtl::AutoDScalar V6, adtl::AutoDScalar V7,
                             double e, double eps0, double cell_epsr,
                             tbox::Array<double> dop_node,
                             tbox::Pointer<tbox::Matrix<double> > ele_mat_E,
                             tbox::Pointer<tbox::Vector<double> > ele_volume);

    adtl::AutoDScalar get_fV(double a, double b, double u_t, double dt,
                                                  tbox::Array<adtl::AutoDScalar> V_node,
                                                  adtl::AutoDScalar h,
                                                  double ele_volume,
                                                  int index,
                                                  tbox::Pointer<tbox::Matrix<double> > ele_mat_E);

    adtl::AutoDScalar get_fn(double a, double u_t, double dt,
                                                  tbox::Array<double> n_node,
                                                  tbox::Array<double> V_node,
                                                  adtl::AutoDScalar h,
                                                  double ele_volume,
                                                  int index,
                                                  tbox::Array<tbox::Array<int> > r_node,
                                                  tbox::Array<tbox::Array<double> > r_integ,
                                                  double mun, double Vt, double e);

    adtl::AutoDScalar get_fp(double a, double u_t, double dt,
                                                  tbox::Array<double> p_node,
                                                  tbox::Array<double> V_node,
                                                  adtl::AutoDScalar h,
                                                  double ele_volume,
                                                  int index,
                                                  tbox::Array<tbox::Array<int> > r_node,
                                                  tbox::Array<tbox::Array<double> > r_integ,
                                                  double mup, double Vt, double e);

    adtl::AutoDScalar get_fV(double a, double b, double u_t, double dt,
                                                  tbox::Array<double> V_node,
                                                  adtl::AutoDScalar h,
                                                  double ele_volume,
                                                  int index,
                                                  tbox::Pointer<tbox::Matrix<double> > ele_mat_E);

    adtl::AutoDScalar get_fVogamma(double a, double u_t, double dt,
                                                        double b,
                                                        adtl::AutoDScalar Vogamma_node,
                                                        adtl::AutoDScalar h,
                                                        double ele_volume);


    /**
     * @brief DataOutput
     * @param filename
     * @param data
     */
    void DataOutput(string filename, tbox::Pointer<pdat::VectorData<NDIM, double> > data);

    /**
     * @brief Variables_Functions::DataOutput
     * @param filename
     * @param data
     */
    void DataOutput(string outputname,
                    int i_bc,
                    double data,
                    const double time, const double dt){
        stringstream file;
        file << outputname << "_" << i_bc;
        string file_dir = file.str();
        // open the file stream
        ofstream streamout;
        if(int(time/dt)==0){
            streamout.open(file_dir,ios::out);
        }else{
            streamout.open(file_dir,ios::app);
        }
        //输出电流信息及其对应的电压、时间
        streamout<<data<<"\n";
    };

    /**
     * @brief Variables_Functions::Errorutput
     * @param filename
     * @param iter
     * @param error
     */
    void ErrorOutput(string outputname,
                    int iter,
                    double error,
                     const double time, const double dt){
        stringstream file;
        file << outputname;
        string file_dir = file.str();
        // open the file stream
        ofstream streamout;
        if(int(time/dt)==0){
            streamout.open(file_dir,ios::out);
        }else{
            streamout.open(file_dir,ios::app);
        }
        //输出电流信息及其对应的电压、时间
        streamout<<iter<<" "<<error<<"\n";
    };

    /**
     * @brief Variables_Functions::Constraint
     * @param data
     */
    double Constraint(double data, double constraint_value);

    void calculateDop(tbox::Pointer<pdat::NodeData<NDIM, double> > node_coord,
                      int num_nodes,
                      tbox::Pointer<pdat::NodeData<NDIM, double> > ND_data,
                      tbox::Pointer<pdat::NodeData<NDIM, double> > NA_data,
                      tbox::Array<doplist> Dopassignment);
    double getND(tbox::Array<doplist> Dopassignment,int index,
                 double Nb, double x, double y, double z);
    double getNA(tbox::Array<doplist> Dopassignment,int index,
                 double Nb, double x, double y, double z);

    /**
     * @brief face_CV_integ
     * @param vertex
     * @return
     */
    tbox::Array<double> face_CV_integ(tbox::Array<hier::DoubleVector<NDIM> > vertex);

    /**
     * @brief quadrilateral_area
     * @param crd
     * @return
     */
    double quadrilateral_area(tbox::Array<hier::DoubleVector<NDIM> > crd);

    void set_DD_type(string DD_type){
        my_DD_type = DD_type;
    }

private:
    string my_DD_type;
};

#endif /* SOURCE_DIRECTORY__VARIABLES_FUNCTIONS_H_ */
