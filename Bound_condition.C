/*
 * source_input.C
 *
 *  Created on: May 12, 2017
 *      Author: tong
 */

#include "Bound_condition.h"
#include "Variables_Functions.h"

#include "InputManager.h"

using namespace std;
using namespace JAUMIN;

bound_condition::bound_condition(){
    string voltage_dir = "../input/Vgd.txt";

    ////////////////////////////////////////////////////////
    std::fstream input(voltage_dir);

    if(!input) {
        TBOX_ERROR("fail to open file: "<< voltage_dir << endl);
    }

    string line;

    while( getline(input, line) )
    {
        stringstream buf(line);
        double Vg;
        double Vd;

        //read and convert data
        buf >> Vg;
        buf >> Vd;

        Vg_array.push_back(Vg);
        Vd_array.push_back(Vd);
    }
    input.close();
}

bound_condition::~bound_condition(){

}

/** pulse source models **/
double bound_condition::source_assign(string source_name,
        const double time,
        const double dt,
        tbox::Array<double> coe){
    double amp;

    if(source_name == "GaussianPulse"){
        gaussian_pulse(coe, time, dt, amp);
    }else if(source_name == "RectangularPulse"){
        rectangular_pulse(coe, time, dt, amp);
    }else if(source_name == "SawPulse"){
        saw_pulse(coe, time, dt, amp);
    }else if(source_name == "CosinePulse"){
        cosine_pulse(coe, time, dt, amp);
    }else if(source_name == "StepPulse"){
        step_pulse(coe, time, dt, amp);
    }else if(source_name == "SteadyVolt"){
        steady_volt(coe, time, dt, amp);
    }else if(source_name == "GND"){
        grouded(amp);
    }else if(source_name == "mixSteadyVolt"){
        mixsteady_volt(coe, time, dt, amp);
    }else if(source_name == "mixStepPulse"){
        mixstep_pulse(coe, time, dt, amp);
    }else if(source_name == "Vg"){
        mixsteady_volt(coe, time, dt, amp);
    }else if(source_name == "Vd"){
        mixstep_pulse(coe, time, dt, amp);
    }
    return amp;
}

/** Poisson Equation Direchlet Boundary Value **/
double bound_condition::e_bound_value(string contact_type,
                                      const double time,
                                      const double dt,
                                      double kb,
                                      double e,
                                      tbox::Array<double> coe){
    double phi_contact;
    if (contact_type == "Ohmic_contact"){
        double Vapp = coe[0];
        double T_node = coe[1];
        double affinity = coe[2];
        double eg = coe[3];
        double nc = coe[4];
        double nv = coe[5];
        double ND = coe[6];
        double NA = coe[7];
        double ni = coe[8];

        double thermo_volt = kb*T_node/e;
        double l = asinh((ND-NA)/ni/2);
        double l1 = Vapp + thermo_volt*asinh((ND-NA)/ni/2) - affinity - eg/2 - (thermo_volt/2)*log(nc/nv);
        phi_contact = Vapp + thermo_volt*asinh((ND-NA)/ni/2) - affinity - eg/2 - (thermo_volt/2)*log(nc/nv);

    }else if (contact_type == "Schottky_contact"){
        double Vapp = coe[0];
        double phi_B = coe[2];
        phi_contact = Vapp - phi_B;
    }else if (contact_type == "Gate_contact"){
        double Vapp = coe[0];
        double phi_B = coe[2];
        phi_contact = Vapp - phi_B;
    }else{
        cout<<"wrong contact type!"<<endl;
    }
    return phi_contact;
}

/** Current Continuity Equation Boundary Value **/
void bound_condition::c_bound_value(string contact_type,
                                    const double time,
                                    const double dt,
                                    tbox::Array<double> coe,
                                    double &fixed_e,
                                    double &fixed_h) {
    if (contact_type == "Ohmic_contact"){
        double Vapp = coe[0];
        double T_node = coe[1];
        double affinity = coe[2]; // V
        double eg = coe[3]; //V
        double nc = coe[4];
        double nv = coe[5];
        double phi = coe[6];

        double thermo_volt = kB*T_node/q_e;
        double eta_n = (phi - affinity - Vapp)/thermo_volt;
        double eta_h = ((Vapp + eg + affinity) - phi)/thermo_volt;
        fixed_e = nc*exp(eta_n);
        fixed_h = nv*exp(eta_h);
    }else if (contact_type == "Schottky_contact"){

    }else{
        cout<<"wrong contact type!"<<endl;
    }
}

double bound_condition::carn_bound_value(string contact_type,
                                         bool Merger_semi,
                                         const double time,
                                         const double dt,
                                         double kb,
                                         double e,
                                         tbox::Array<double> coe){
    double fixed_carn;
    tbox::Pointer<Variables_Functions> Var_func = new Variables_Functions();

    if (contact_type == "Ohmic_contact"){
        double Vapp = coe[0];
        double T_node = coe[1];
        double affinity = coe[2];
        double eg = coe[3];
        double phin = coe[4];
        double phih = coe[5];
        double phi = coe[6];
        double dop = coe[7];
        double ni = coe[8];

        if (Merger_semi){
            double thermo_volt = kb*T_node/e;
            double eta_n = (phi + affinity - phin)/thermo_volt;
            double eta_h = -((phi + affinity + eg) - phih)/thermo_volt;
            double gama_n = Var_func->gama_n_p(eta_n);
            double gama_p = Var_func->gama_n_p(eta_h);
            fixed_carn = 0.5*dop + 0.5*sqrt(pow(dop,2.0)+4*gama_n*gama_p*pow(ni,2.0));
        }else{
            double thermo_volt = kb*T_node/e;
            double eta_n = (phi + affinity - phin)/thermo_volt;
            double eta_h = -((phi + affinity + eg) - phih)/thermo_volt;
            double gama_n = Var_func->gama_n_p(eta_n);
            double gama_p = Var_func->gama_n_p(eta_h);
            if (dop>0){
                fixed_carn = 0.5*dop + 0.5*sqrt(pow(dop,2.0)+4*gama_n*gama_p*pow(ni,2.0));
            }else{
                double tmp_carp = -0.5*dop + 0.5*sqrt(pow(dop,2.0)+4*gama_n*gama_p*pow(ni,2.0));
                fixed_carn = pow(ni,2.0)/tmp_carp;
            }
        }
    }else if (contact_type == "Schottky_contact"){

    }else{
    cout<<"wrong contact type!"<<endl;
    }
    return fixed_carn;
}

double bound_condition::carh_bound_value(string contact_type,
                                         bool Merger_semi,
                                         const double time,
                                         const double dt,
                                         double kb,
                                         double e,
                                         tbox::Array<double> coe){
    double fixed_carh;
    tbox::Pointer<Variables_Functions> Var_func = new Variables_Functions();
    if (contact_type == "Ohmic_contact"){
    double Vapp = coe[0];
    double T_node = coe[1];
    double affinity = coe[2]; // V
    double eg = coe[3]; //V
        double phin = coe[4];
        double phih = coe[5];
    double phi = coe[6];
        double dop = coe[7];
        double ni = coe[8];
        if (Merger_semi){
            double thermo_volt = kb*T_node/e;
            double eta_n = ((phi + affinity) - phin)/thermo_volt;
            double eta_h = -((phi + affinity + eg) - phih)/thermo_volt;
            double gama_n = Var_func->gama_n_p(eta_n);  //1
            double gama_p = Var_func->gama_n_p(eta_h);  //1
            fixed_carh = -0.5*dop + 0.5*sqrt(pow(dop,2.0)+4*gama_n*gama_p*pow(ni,2.0));
        }else{
            double thermo_volt = kb*T_node/e;
            double eta_n = ((phi + affinity) - phin)/thermo_volt;
            double eta_h = -((phi + affinity + eg) - phih)/thermo_volt;
            double gama_n = Var_func->gama_n_p(eta_n);
            double gama_p = Var_func->gama_n_p(eta_h);
            if(dop>0){
                double tmp_carn = 0.5*dop + 0.5*sqrt(pow(dop,2.0)+4*gama_n*gama_p*pow(ni,2.0));
                fixed_carh = pow(ni,2.0)/tmp_carn;
            }else{
                fixed_carh = -0.5*dop + 0.5*sqrt(pow(dop,2.0)+4*gama_n*gama_p*pow(ni,2.0));
            }
        }
    }else if (contact_type == "Schottky_contact"){

    }else{
    cout<<"wrong contact type!"<<endl;
    }
    return fixed_carh;
}

/** ************** pulse functions*************** **/
/** gaussian pulse **/
void bound_condition::gaussian_pulse(tbox::Array<double> coe,
        const double time,
        const double dt,
        double &amp){

}

/** rectangular pulse **/
void bound_condition::rectangular_pulse(tbox::Array<double> coe,
        const double time,
        const double dt,
        double &amp){
    double width = coe[0];
    double hight = coe[1];
    double first = coe[2];
    int tmp_n = (time/width);
    if(tmp_n%2 < dt){
        amp = (1-first)*hight;
    }else{
        amp = first*hight;
    }
}

/** saw wave pulse **/
void bound_condition::saw_pulse(tbox::Array<double> coe,
        const double time,
        const double dt,
        double &amp){

}

/** cosine pulse **/
void bound_condition::cosine_pulse(tbox::Array<double> coe,
        const double time,
        const double dt,
        double &amp){
    double freq = coe[0];
    double amplitude = coe[1];
    double phase = coe[2];
    amp = cos(2*M_PI*freq*time+phase);
}

/** stair pulse **/
void bound_condition::step_pulse(tbox::Array<double> coe,  /// change to structure
        const double time,
        const double dt,
        double &amp){
    double volt = coe[0];
    double ramp = coe[1];
    double ramp_time = coe[2];

    if (time < ramp_time){
        amp = ramp*time;
    }else{
        amp = volt;
    }
}

/** stair pulse **/
void bound_condition::mixstep_pulse(tbox::Array<double> coe,  /// change to structure
        const double time,
        const double dt,
        double &amp){
//    double volt = coe[0];
//    double ramp = coe[1];
//    double ramp_time = coe[2];

    int step = int(time/dt);
    amp = Vd_array[step];

    #if 0
    if (time <= ramp_time){
        amp = ramp*time;
    }else if (time > ramp_time){
        int step = (time-ramp_time)/0.05;
        if(step==0)
            amp = 2.0;
        else if(step==1)
            amp = 3.0;
        else if(step>=2)
            amp = 4.0;
    }
    #endif
}

/** steady voltage **/
void bound_condition::mixsteady_volt(tbox::Array<double> coe,
        const double time,
        const double dt,
        double &volt){
    int step = int(time/dt);
    volt = Vg_array[step];
//    if (time <= 2.15){
//        volt = 0.01;
//    }else if (time > 2.15 & time <=12.999){
//        volt = (time - 2.1)*0.5 + 0.01;
//    }else if (time > 12.999){
//		int step = (time - 13)/0.05+1;
//		volt = step * 10;
//	}
}

/** steady voltage **/
void bound_condition::steady_volt(tbox::Array<double> coe,
        const double time,
        const double dt,
        double &volt){
    volt = coe[0];
}

/** GND **/
void bound_condition::grouded(double &amp){
    amp = 0.0;
}

/** ***************load boundary condtion assign from input files*************** **/

/** load electrical boundary assignment **/
tbox::Array<electrobc> bound_condition::getEBCassignment(tbox::Pointer<tbox::Database> db){

    /** get the route for file containing boundary conditions **/
    string bc_dir = db->getDatabase("EBClist")->getString("ebc_list_route");

    /** setup database for the ebc list **/
    tbox::Pointer<tbox::Database> bclist_db = new tbox::InputDatabase("bclist_db");
    tbox::InputManager::getManager()->parseInputFile(bc_dir, bclist_db);

    tbox::Pointer<tbox::Database> bc_list_db = bclist_db->getDatabase("bc_list");
    int num_bc = bc_list_db->getInteger("num_bc");
    EBCassignment.resizeArray(num_bc);

    for(int i=0;i<num_bc;++i){
        stringstream bc_file;
        bc_file << "electrobc_" << i;
        string ebc_file_dir = bc_list_db->getString(bc_file.str());
        /** setup database for electro bc **/
        tbox::Pointer<tbox::Database> ebc_assign_db = new tbox::InputDatabase("ebc_assign_db");
        tbox::InputManager::getManager()->parseInputFile(ebc_file_dir, ebc_assign_db);

        //get database mat_domain
        tbox::Pointer<tbox::Database> sub_ebc_assign_db = ebc_assign_db->getDatabase("ebc_assign");

        EBCassignment[i].bc_type = sub_ebc_assign_db->getString("bc_type");
        EBCassignment[i].bc_value = sub_ebc_assign_db->getDouble("bc_value");
        EBCassignment[i].bc_no = sub_ebc_assign_db->getInteger("bc_num");
        EBCassignment[i].bc_face = sub_ebc_assign_db->getIntegerArray("bc_face");
        EBCassignment[i].sc_type = sub_ebc_assign_db->getString("sc_type");
        EBCassignment[i].sc_coe = sub_ebc_assign_db->getDoubleArray("sc_coe");
    }
    return EBCassignment;
}

/** load dopant condition assignment **/
tbox::Array<doplist> bound_condition::getDopassignment(tbox::Pointer<tbox::Database> db){

    /** get the route for file containing boundary conditions **/
    string dop_dir = db->getDatabase("Doplist")->getString("dop_list_route");

    /** setup database for the ebc list **/
    tbox::Pointer<tbox::Database> doplist_db = new tbox::InputDatabase("doplist_db");
    tbox::InputManager::getManager()->parseInputFile(dop_dir, doplist_db);

    tbox::Pointer<tbox::Database> dop_list_db = doplist_db->getDatabase("dop_list");
    int num_dop = dop_list_db->getInteger("num_dop");
    Dopassignment.resizeArray(num_dop);

    for(int i=0;i<num_dop;++i){
        stringstream dop_file;
        dop_file << "dop_" << i;
        string dop_file_dir = dop_list_db->getString(dop_file.str());
        /** setup database for electro bc **/
        tbox::Pointer<tbox::Database> dop_assign_db = new tbox::InputDatabase("dop_assign_db");
        tbox::InputManager::getManager()->parseInputFile(dop_file_dir, dop_assign_db);

        //get database mat_domain
        tbox::Pointer<tbox::Database> sub_dop_assign_db = dop_assign_db->getDatabase("dop_assign");

        Dopassignment[i].dop_type = sub_dop_assign_db->getString("dop_type");
        if(i==0){
            if(strcmp((Dopassignment[i].dop_type).c_str(),"n")==0){
                Dopassignment[i].ND0 = sub_dop_assign_db->getDouble("ND0");
                Dopassignment[i].Nb = Dopassignment[i].ND0;
            }else if(strcmp((Dopassignment[i].dop_type).c_str(),"p")==0){
                Dopassignment[i].NA0 = sub_dop_assign_db->getDouble("NA0");
                Dopassignment[i].Nb = Dopassignment[i].NA0;
            }
        }else{
            if(strcmp((Dopassignment[i].dop_type).c_str(),"n")==0){
                Dopassignment[i].r0 = sub_dop_assign_db->getDoubleArray("r0");
                Dopassignment[i].ND0 = sub_dop_assign_db->getDouble("ND0");
                Dopassignment[i].W = sub_dop_assign_db->getDouble("W");
                Dopassignment[i].D = sub_dop_assign_db->getDouble("D");
                Dopassignment[i].H = sub_dop_assign_db->getDouble("H");
                Dopassignment[i].dj = sub_dop_assign_db->getDoubleArray("dj");
                Dopassignment[i].Nb = Dopassignment[0].Nb;
            }else if(strcmp((Dopassignment[i].dop_type).c_str(),"p")==0){
                Dopassignment[i].r0 = sub_dop_assign_db->getDoubleArray("r0");
                Dopassignment[i].NA0 = sub_dop_assign_db->getDouble("NA0");
                Dopassignment[i].W = sub_dop_assign_db->getDouble("W");
                Dopassignment[i].D = sub_dop_assign_db->getDouble("D");
                Dopassignment[i].H = sub_dop_assign_db->getDouble("H");
                Dopassignment[i].dj = sub_dop_assign_db->getDoubleArray("dj");
                Dopassignment[i].Nb = Dopassignment[0].Nb;
            }else{
                cout<<"error dopant type!"<<endl;
            }
        }
    }
    return Dopassignment;
}

/** load carrier density boundary assignment **/
tbox::Array<carrierbc> bound_condition::getCBCassignment(tbox::Pointer<tbox::Database> db){

    /** get the route for file containing boundary conditions **/
    string bc_dir = db->getDatabase("CBClist")->getString("cbc_list_route");

    /** setup database for the ebc list **/
    tbox::Pointer<tbox::Database> bclist_db = new tbox::InputDatabase("bclist_db");
    tbox::InputManager::getManager()->parseInputFile(bc_dir, bclist_db);

    tbox::Pointer<tbox::Database> bc_list_db = bclist_db->getDatabase("bc_list");
    int num_bc = bc_list_db->getInteger("num_bc");
    CBCassignment.resizeArray(num_bc);

    for(int i=0;i<num_bc;++i){
        stringstream bc_file;
        bc_file << "carrierbc_" << i;
        string cbc_file_dir = bc_list_db->getString(bc_file.str());
        /** setup database for electro bc **/
        tbox::Pointer<tbox::Database> cbc_assign_db = new tbox::InputDatabase("cbc_assign_db");
        tbox::InputManager::getManager()->parseInputFile(cbc_file_dir, cbc_assign_db);

        tbox::Pointer<tbox::Database> sub_cbc_assign_db = cbc_assign_db->getDatabase("cbc_assign");
        CBCassignment[i].bc_type = sub_cbc_assign_db->getString("bc_type");
        CBCassignment[i].bc_value = sub_cbc_assign_db->getDouble("bc_value");
        CBCassignment[i].bc_no = sub_cbc_assign_db->getInteger("bc_num");
        CBCassignment[i].bc_face = sub_cbc_assign_db->getIntegerArray("bc_face");
        CBCassignment[i].sc_type = sub_cbc_assign_db->getString("sc_type");
        CBCassignment[i].sc_coe = sub_cbc_assign_db->getDoubleArray("sc_coe");
    }
    return CBCassignment;
}

/** load thermal boundary assignment **/
tbox::Array<thermobc> bound_condition::getTBCassignment(tbox::Pointer<tbox::Database> db){

    /** get the route for file containing boundary conditions **/
    string bc_dir = db->getDatabase("TBClist")->getString("tbc_list_route");

    /** setup database for the thermo bc **/
    tbox::Pointer<tbox::Database> bclist_db = new tbox::InputDatabase("bclist_db");
    tbox::InputManager::getManager()->parseInputFile(bc_dir, bclist_db);

    tbox::Pointer<tbox::Database> bc_list_db = bclist_db->getDatabase("bc_list");
    int num_bc = bc_list_db->getInteger("num_bc");
    TBCassignment.resizeArray(num_bc);

    for(int i=0;i<num_bc;++i){
        stringstream bc_file;
        bc_file << "thermobc_" << i;
        string tbc_file_dir = bc_list_db->getString(bc_file.str());

        /** setup database for thermo bc **/
        tbox::Pointer<tbox::Database> tbc_assign_db = new tbox::InputDatabase("tbc_assign_db");
        tbox::InputManager::getManager()->parseInputFile(tbc_file_dir, tbc_assign_db);

        tbox::Pointer<tbox::Database> sub_tbc_assign_db = tbc_assign_db->getDatabase("tbc_assign");
        TBCassignment[i].bc_type = sub_tbc_assign_db->getString("bc_type");
        TBCassignment[i].bc_value = sub_tbc_assign_db->getDouble("bc_value");
        TBCassignment[i].bc_no = sub_tbc_assign_db->getInteger("bc_num");
        TBCassignment[i].bc_face = sub_tbc_assign_db->getIntegerArray("bc_face");
        TBCassignment[i].bc_coe = sub_tbc_assign_db->getDoubleArray("bc_coe");
    }
    return TBCassignment;
}



