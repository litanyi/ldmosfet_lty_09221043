/*
 * source_input.h
 *
 *  Created on: May 12, 2017
 *      Author: tong
 */

#ifndef SOURCE_DIRECTORY__BOUND_CONDITION_H_
#define SOURCE_DIRECTORY__BOUND_CONDITION_H_

#include <string>
#include "math.h"
#include "DoubleVector.h"

#include "Global_variable.h"

using namespace std;
using namespace JAUMIN;

struct electrobc{
    string bc_type;
    double bc_value;
    int bc_no;
    tbox::Array<int> bc_face;
    string sc_type;
    tbox::Array<double> sc_coe;
};

struct doplist{
    string dop_type;
    tbox::Array<double> r0;
    double ND0;
    double NA0;
    double W;
    double D;
    double H;
    tbox::Array<double> dj;
    double Nb;
};

struct carrierbc{
    string bc_type;
    double bc_value;
    int bc_no;
    tbox::Array<int> bc_face;
    string sc_type;
    tbox::Array<double> sc_coe;
};

struct thermobc{
    string bc_type;
    double bc_value;
    int bc_no;
    tbox::Array<int> bc_face;
    tbox::Array<double> bc_coe;
};

// tbox::Array<pulsetype> BCassignment;

class bound_condition{
public:

    // construction functions
    bound_condition();

    // anti construction functions
    ~bound_condition();

    // source assignment
    // "coe" contains the variable require by the source function
    double source_assign(string source_name,
            const double time,
            const double dt,
            tbox::Array<double> coe);

    // "coe" contains the variables required
    double e_bound_value(string contact_type,
             const double time,
             const double dt,
             double kb,
             double e,
             tbox::Array<double> coe);

    // Gaussian pulse
    /*
     * @param coe array double, the coefficients used to define the gaussian pulse
     */
    void gaussian_pulse(tbox::Array<double> coe,
            const double time,
            const double dt,
            double &amp);

    // Rectangular pulse
    /*
     * @param coe array double, the coefficients used to define the rectangular pulse
     */
    void rectangular_pulse(tbox::Array<double> coe,
            const double time,
            const double dt,
            double &amp);

    // saw pulse
    /*
     * @param coe array double, the coefficients used to define the saw pulse
     */
    void saw_pulse(tbox::Array<double> coe,
            const double time,
            const double dt,
            double &amp);

    // cosine pulse
    /*
     * @param coe array double, the coefficients used to define the cosine pulse
     */
    void cosine_pulse(tbox::Array<double> coe,
                  const double time,
              const double dt,
              double &amp);

    // step pulse
    /*
     * @param coe array double, the coefficients used to define the steady pulse
     */
    void step_pulse(tbox::Array<double> coe,
            const double time,
            const double dt,
            double &amp);

    // steady voltage source
    void steady_volt(tbox::Array<double> coe,
            const double time,
            const double dt,
            double &volt);

    // step pulse
    /*
     * @param coe array double, the coefficients used to define the steady pulse
     */
    void mixstep_pulse(tbox::Array<double> coe,
            const double time,
            const double dt,
            double &amp);

    // steady voltage source
    void mixsteady_volt(tbox::Array<double> coe,
            const double time,
            const double dt,
            double &volt);

    // GND
    void grouded(double &amp);

    /** carrier density boundary assignment **/
    void c_bound_value(string contact_type,
                   const double time,
               const double dt,
               tbox::Array<double> coe,
               double &fixed_e,
               double &fixed_h);

    /** electron density boundary assignment **/
    double carn_bound_value(string contact_type, bool Merger_semi,
                            const double time,
                            const double dt,
                            double kb,
                            double e,
                            tbox::Array<double> coe);

    /** hole density boundary assignment **/
    double carh_bound_value(string contact_type, bool Merger_semi,
                            const double time,
                            const double dt,
                            double kb,
                            double e,
                            tbox::Array<double> coe);


    /** read electrical boundary condition assignment **/
    tbox::Array<electrobc> getEBCassignment(tbox::Pointer<tbox::Database> db);

    /** read dopant condition assignment **/
    tbox::Array<doplist> getDopassignment(tbox::Pointer<tbox::Database> db);

    /** read carrier continuity boundary information **/
    tbox::Array<carrierbc> getCBCassignment(tbox::Pointer<tbox::Database> db);

    /** read thermal boundary condition assignment **/
    tbox::Array<thermobc> getTBCassignment(tbox::Pointer<tbox::Database> db);

private:
    tbox::Array<electrobc> EBCassignment;
    tbox::Array<doplist> Dopassignment;
    tbox::Array<carrierbc> CBCassignment;
    tbox::Array<thermobc> TBCassignment;

    vector<double> Vg_array;
    vector<double> Vd_array;
};


#endif /* SOURCE_DIRECTORY__BOUND_CONDITION_H_ */
