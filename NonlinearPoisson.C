//
// 文件名:     NonlinearPoisson.C
// 软件包:     JAUMIN
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: Tue Apr 14 16:39:38 2015 $
// 描述　:     非线性Poisson方程网格片策略类实现
// 类别　:     %Internal File% ( Don't delete this line )
//

#include "VectorVariable.h"
#include "CSRMatrixVariable.h"
#include "VectorData.h"
#include "CSRMatrixData.h"
#include "ElementManager.h"
#include "ShapeFunctionManager.h"
#include "IntegratorManager.h"
#include "BaseElement.h"
#include "PatchTopology.h"
#include "PatchGeometry.h"
#include "Patch.h"
#include "NonlinearPoisson.h"
#include "Matrix.h"
#include "Vector.h"
#include "NodeVariable.h"
#include "CellVariable.h"
#include "FaceVariable.h"
#include "NodeData.h"
#include "CellData.h"
#include "FaceData.h"
#include "LinearQuadrangle.h"
#include "LinearTriangle.h"
#include "LinearHexahedron.h"
#include "LinearTetrahedron.h"
#include "IntegratorManager.h"
#include "ShapeFunctionManager.h"
#include "ElementManager.h"
#include "ShapeFunctionManager.h"
#include "IntegratorManager.h"
#include "QuadrangleIntegrator.h"
#include "QuadrangleShapeFunction.h"
#include "TriangleIntegrator.h"
#include "TriangleShapeFunction.h"
#include "HexahedronShapeFunction.h"
#include "HexahedronIntegrator.h"
#include "TetrahedronShapeFunction.h"
#include "TetrahedronIntegrator.h"
#include "math.h"
#include "VectorVariable.h"

#include "MaterialManager.h"
#include "matRegister.h"

#include "Bound_condition.h"
#include "Variables_Functions.h"
#include "SerialRead.h"
#include "Global_variable.h"

#include <cmath>
#include <fstream>

/** ***********************************************************************
 *
 * 构造函数.
 *
 *********************************************************************** **/

NonlinearPoisson::NonlinearPoisson(const string& object_name,
                                   bool is_from_restart,
                                   tbox::Pointer<tbox::Database> input_db)
{
    /** (形函数, 积分器, 单元)管理器 **/
    tbox::Pointer<ElementManager<NDIM> >
            ele_manager = ElementManager<NDIM>::getManager();
    tbox::Pointer<ShapeFunctionManager<NDIM> >
            func_manager = ShapeFunctionManager<NDIM>::getManager();
    tbox::Pointer<IntegratorManager<NDIM> >
            integrator_manager = IntegratorManager<NDIM>::getManager();

    /** material manager **/
    tbox::Pointer<matRegister> mat_regist = new matRegister();
    tbox::Pointer<MaterialManager<NDIM> > mat_manager = mat_regist->registerMatLib();

#if (NDIM == 2)                                         /** 三角形 **/
    tbox::Pointer<LinearTriangle> linear_ele =
        new LinearTriangle("LinearTriangle");
    tbox::Pointer<TriangleIntegrator> linear_int =
        new TriangleIntegrator(1, "Triangle");
    tbox::Pointer<Triangle_1_ShapeFunction> linear_sha =
        new Triangle_1_ShapeFunction("Triangle");
#else
    /** 四边形 **/
    tbox::Pointer<LinearQuadrangle> linear_ele_quadrangle =
        new LinearQuadrangle("LinearQuadrangle");
    tbox::Pointer<QuadrangleIntegrator> linear_int_quadrangle =
        new QuadrangleIntegrator(2, "Quadrangle");
    tbox::Pointer<QuadrangleShapeFunction> linear_sha_quadrangle =
        new QuadrangleShapeFunction("Quadrangle");
    /** 六面体 **/
    tbox::Pointer<LinearHexahedron> linear_ele =
        new LinearHexahedron("LinearHexahedron");
    tbox::Pointer<HexahedronIntegrator> linear_int =
        new HexahedronIntegrator(2, "Hexahedron");
    tbox::Pointer<HexahedronShapeFunction> linear_sha =
        new HexahedronShapeFunction("Hexahedron");
    /** 四面体 **/
    tbox::Pointer<LinearTetrahedron> linear_ele_tetrahedron =
        new LinearTetrahedron("LinearTetrahedron");
    tbox::Pointer<TetrahedronIntegrator> linear_int_tetrahedron =
        new TetrahedronIntegrator(2, "Tetrahedron");
    tbox::Pointer<TetrahedronShapeFunction> linear_sha_tetrahedron =
        new TetrahedronShapeFunction("Tetrahedron");
#endif
    /** 将六面体单元添加到单元管理器 **/
    ele_manager->addElement(linear_ele);
    integrator_manager->addIntegrator(linear_int);
    func_manager->addShapeFunction(linear_sha);
    d_object_name = object_name;

    /** 将四面体单元添加到单元管理器 **/
    ele_manager->addElement(linear_ele_tetrahedron);
    integrator_manager->addIntegrator(linear_int_tetrahedron);
    func_manager->addShapeFunction(linear_sha_tetrahedron);

    /** 将四边形单元添加到单元管理器 **/
    ele_manager->addElement(linear_ele_quadrangle);
    integrator_manager->addIntegrator(linear_int_quadrangle);
    func_manager->addShapeFunction(linear_sha_quadrangle);

    /** read material assignment **/
    tbox::Pointer<mat_management> mat_list = new mat_management();
    material_list = mat_list->getmateriallist(input_db);
    n_mat = material_list.getSize();

    /** entity mat list **/
    int num_en = 0;
    for (int i_m = 0; i_m < n_mat; ++i_m) num_en = num_en + material_list[i_m].num_entity;

    entity_mat_list.resizeArray(num_en);
    for (int i_m = 0; i_m < n_mat; ++i_m){
        for (int i_en = 0; i_en < material_list[i_m].num_entity; ++i_en){
            entity_mat_list[material_list[i_m].entity_no[i_en]-1] = i_m+1;
        }
    }
    n_entity = num_en;

    /** load data from input files && restart files **/
    if (is_from_restart) {
      getFromRestart(input_db);
    } else {
      getFromInput(input_db);
    }

    /** register model variables **/
    registerModelVariable();

    /** 将当前类对象注册为重启动对象 **/
    tbox::RestartManager::getManager()->registerRestartItem(d_object_name, this);

    /** obtain the pointor of the database **/
    GetPointor_db = input_db.getPointer();

    /** read boundary condition information **/
    tbox::Pointer<bound_condition> bc_manager = new bound_condition();
    EBoundassignment = bc_manager->getEBCassignment(input_db);
    Dopassignment = bc_manager->getDopassignment(input_db);
    //CBoundassignment = bc_manager->getCBCassignment(input_db);

    /** read polarized charge density from 'input_db' **/
    polarized_p = input_db->getDatabase("Polarized_charge")->getDouble("polarized_p");

    Dt=input_db->getDatabase("Timecontroller")->getDouble("time_step");
    dopant_file_name = input_db->getDatabase("Dopant")->getString("dopant_file_name");
    sol_file_name = input_db->getDatabase("Sol")->getString("sol_file_name");
    num_sol_file_name = input_db->getDatabase("Sol")->getDouble("num_sol_file_name");
    transient_state = input_db->getDatabase("Timecontroller")->getBool("transient_state");
    constraint_value = input_db->getDatabase("Constraint")->getDouble("constraint_value");
    enable_unit_scalling = input_db->getDatabase("Variable_scalling")->getBool("enable_unit_scalling");
    //Qt = input_db->getDatabase("Radiation")->getDouble("Qt");
    Not = input_db->getDatabase("Radiation")->getDouble("Not");
    car_bound_mode = input_db->getBool("car_bound_mode");
    save_bas = input_db->getBool("save_bas");
    save_carrier = input_db->getBool("save_carrier");
    save_sol = input_db->getBool("save_sol");
    enable_tid = input_db->getBool("enable_tid");
    get_InitialSol = input_db->getBool("get_InitialSol");
    get_PoissonSolcur = input_db->getBool("get_PoissonSolcur");
    enable_selfDop = input_db->getDatabase("selfDop")->getBool("enable_selfDop");
    num_PoissonSolcur_filename = input_db->getDatabase("PoissonSolcur")->getDouble("num_PoissonSolcur_filename");
    PoissonSolcur_filename = input_db->getDatabase("PoissonSolcur")->getString("PoissonSolcur_filename");

}

/** ***********************************************************************
 * anti construction function
 ********************************************************************** **/
NonlinearPoisson::~NonlinearPoisson() {
  tbox::RestartManager::getManager()->unregisterRestartItem(d_object_name);
}

/** ***********************************************************************
 * 注册变量和数据片.
 ********************************************************************** **/
void NonlinearPoisson::registerModelVariable() {
  /** 获取有限元变量的数据库 **/
  hier::VariableDatabase<NDIM> *
    variable_db = hier::VariableDatabase<NDIM>::getDatabase();

  /** 创建自由度信息，参数信息（四个逻辑型表示点，边，面，体上自由度是否非零 **/
  d_dof_info = new solv::DOFInfo<NDIM>(true,false,false,false);
  d_dof_semi_info = new solv::DOFInfo<NDIM>(true,false,false,false);
  d_dof_Vogamma_info = new solv::DOFInfo<NDIM>(true,false,false,false);

  /** 变量上下文 **/
  tbox::Pointer<hier::VariableContext> nnew = variable_db->getContext("NEW");
  tbox::Pointer<hier::VariableContext> current = variable_db->getContext("CURRENT");
  tbox::Pointer<hier::VariableContext> scratch = variable_db->getContext("SCRATCH");

  /**
   * @brief current continuity equation && poisson equation
   **/
  /** current conservation equation **/
  tbox::Pointer<pdat::CSRMatrixVariable<NDIM,double> >
    CC_matrix   = new pdat::CSRMatrixVariable<NDIM,double>("cc_matrix",d_dof_info);
  tbox::Pointer<pdat::VectorVariable<NDIM,double> >
    CC_rhs = new pdat::VectorVariable<NDIM,double>("cc_rhs",d_dof_info);
  tbox::Pointer<pdat::VectorVariable<NDIM,double> >
    cc_sol_vec = new pdat::VectorVariable<NDIM,double>("cc_sol_vec",d_dof_info);

  /** Poisson Equation -- initialization **/
  tbox::Pointer<pdat::CSRMatrixVariable<NDIM,double> >
    P_matrix   = new pdat::CSRMatrixVariable<NDIM,double>("p_matrix",d_dof_info);
  tbox::Pointer<pdat::VectorVariable<NDIM,double> >
    p_rhs = new pdat::VectorVariable<NDIM,double>("p_rhs",d_dof_info);
  tbox::Pointer<pdat::VectorVariable<NDIM,double> >
    p_sol = new pdat::VectorVariable<NDIM,double>("p_sol",d_dof_info);


  /** Poisson Equation -- Newton **/
  tbox::Pointer<pdat::CSRMatrixVariable<NDIM,double> >
    P_newton_matrix = new pdat::CSRMatrixVariable<NDIM,double>("p_newton_matrix",d_dof_semi_info);
  tbox::Pointer<pdat::VectorVariable<NDIM,double> >
    p_newton_rhs = new pdat::VectorVariable<NDIM,double>("p_newton_rhs",d_dof_semi_info);
  tbox::Pointer<pdat::VectorVariable<NDIM,double> >
    p_newton_sol_vec = new pdat::VectorVariable<NDIM,double>("p_newton_sol_vec",d_dof_semi_info);
  tbox::Pointer<pdat::VectorVariable<NDIM,double> >
    p_stetch = new pdat::VectorVariable<NDIM,double>("p_stetch",d_dof_semi_info);
  tbox::Pointer<pdat::VectorVariable<NDIM,double> >
    p_sol_vec = new pdat::VectorVariable<NDIM,double>("p_sol_vec",d_dof_semi_info);

  /** Vogamma Equation -- Newton **/
  tbox::Pointer<pdat::CSRMatrixVariable<NDIM,double> >
    Vogamma_newton_matrix = new pdat::CSRMatrixVariable<NDIM,double>("Vogamma_newton_matrix",d_dof_Vogamma_info);
  tbox::Pointer<pdat::VectorVariable<NDIM,double> >
    Vogamma_newton_rhs = new pdat::VectorVariable<NDIM,double>("Vogamma_newton_rhs",d_dof_Vogamma_info);
  tbox::Pointer<pdat::VectorVariable<NDIM,double> >
    Vogamma_newton_sol_vec = new pdat::VectorVariable<NDIM,double>("Vogamma_newton_sol_vec",d_dof_Vogamma_info);
  tbox::Pointer<pdat::VectorVariable<NDIM,double> >
    Vogamma_sol_vec = new pdat::VectorVariable<NDIM,double>("Vogamma_sol_vec",d_dof_Vogamma_info);

  tbox::Pointer<pdat::NodeVariable<NDIM,double> >
    Vogamma_plot = new pdat::NodeVariable<NDIM,double>("Vogamma_plot",1);
  tbox::Pointer<pdat::NodeVariable<NDIM,double> >
    Vogammaplus_plot = new pdat::NodeVariable<NDIM,double>("Vogammaplus_plot",1);
  tbox::Pointer<pdat::NodeVariable<NDIM,double> >
    VogammaHplus_plot = new pdat::NodeVariable<NDIM,double>("VogammaHplus_plot",1);
  tbox::Pointer<pdat::NodeVariable<NDIM,double> >
    VogammaH2plus_plot = new pdat::NodeVariable<NDIM,double>("VogammaH2plus_plot",1);
  tbox::Pointer<pdat::NodeVariable<NDIM,double> >
    VogammaH_plot = new pdat::NodeVariable<NDIM,double>("VogammaH_plot",1);
  tbox::Pointer<pdat::NodeVariable<NDIM,double> >
    VogammaH2_plot = new pdat::NodeVariable<NDIM,double>("VogammaH2_plot",1);

  /** plot variable (node) **/
  tbox::Pointer<pdat::NodeVariable<NDIM,double> >
    p_plot = new pdat::NodeVariable<NDIM,double>("p_plot",1);

  /**
   * @brief Carrier transport equation: electron
   **/
  /** matrix_electron **/
  tbox::Pointer<pdat::CSRMatrixVariable<NDIM,double> >
    N_matrix   = new pdat::CSRMatrixVariable<NDIM,double>("n_matrix",d_dof_info);
  tbox::Pointer<pdat::VectorVariable<NDIM,double> >
    n_rhs = new pdat::VectorVariable<NDIM,double>("n_rhs",d_dof_info);
  /** vector variable-electron (solution) **/
  tbox::Pointer<pdat::VectorVariable<NDIM,double> >
    n_stetch = new pdat::VectorVariable<NDIM,double>("n_stetch",d_dof_info);
  tbox::Pointer<pdat::VectorVariable<NDIM,double> >
    n_sol_vec = new pdat::VectorVariable<NDIM,double>("n_sol_vec",d_dof_info);
  tbox::Pointer<pdat::VectorVariable<NDIM,double> >
    logn_sol_vec = new pdat::VectorVariable<NDIM,double>("logn_sol_vec",d_dof_info);
  tbox::Pointer<pdat::VectorVariable<NDIM,double> >
    fermin_vec = new pdat::VectorVariable<NDIM,double>("fermin_vec",d_dof_info);
  tbox::Pointer<pdat::VectorVariable<NDIM,double> >
    fermin_stetch = new pdat::VectorVariable<NDIM,double>("fermin_stetch",d_dof_info);

  /** plot variable-electron (node) **/
  tbox::Pointer<pdat::NodeVariable<NDIM,double> >
    n_plot = new pdat::NodeVariable<NDIM,double>("n_plot",1);
  tbox::Pointer<pdat::NodeVariable<NDIM,double> >
    logn_plot = new pdat::NodeVariable<NDIM,double>("logn_plot",1);
  tbox::Pointer<pdat::NodeVariable<NDIM,double> >
    fermin_plot = new pdat::NodeVariable<NDIM,double>("fermin_plot",1);

  /**
   * @brief Carrier transport equation: hole
   **/
  /** stiff matrix **/
  tbox::Pointer<pdat::CSRMatrixVariable<NDIM,double> >
    H_matrix   = new pdat::CSRMatrixVariable<NDIM,double>("h_matrix",d_dof_info);
  tbox::Pointer<pdat::VectorVariable<NDIM,double> >
    h_rhs = new pdat::VectorVariable<NDIM,double>("h_rhs",d_dof_info);
  /** vector variable-hole (solution) **/
  tbox::Pointer<pdat::VectorVariable<NDIM,double> >
    h_stetch = new pdat::VectorVariable<NDIM,double>("h_stetch",d_dof_info);
  tbox::Pointer<pdat::VectorVariable<NDIM,double> >
    h_sol_vec = new pdat::VectorVariable<NDIM,double>("h_sol_vec",d_dof_info);
  tbox::Pointer<pdat::VectorVariable<NDIM,double> >
    logh_sol_vec = new pdat::VectorVariable<NDIM,double>("logh_sol_vec",d_dof_info);
  tbox::Pointer<pdat::VectorVariable<NDIM,double> >
    fermih_vec = new pdat::VectorVariable<NDIM,double>("fermih_vec",d_dof_info);
  tbox::Pointer<pdat::VectorVariable<NDIM,double> >
    fermih_stetch = new pdat::VectorVariable<NDIM,double>("fermih_stetch",d_dof_info);

  /** plot variable-hole (node) **/
  tbox::Pointer<pdat::NodeVariable<NDIM,double> >
    h_plot = new pdat::NodeVariable<NDIM,double>("h_plot",1);
  tbox::Pointer<pdat::NodeVariable<NDIM,double> >
    logh_plot = new pdat::NodeVariable<NDIM,double>("logh_plot",1);
  tbox::Pointer<pdat::NodeVariable<NDIM,double> >
    fermih_plot = new pdat::NodeVariable<NDIM,double>("fermih_plot",1);


  /**
   * @brief electric field (edge)
   **/

  tbox::Pointer<pdat::CellVariable<NDIM,double> >
    E_vec = new pdat::CellVariable<NDIM,double>("Ex_vec",cell_edge);
  tbox::Pointer<pdat::CellVariable<NDIM,double> >
    E_cell = new pdat::CellVariable<NDIM,double>("E_cell",3);

  /**
   * @brief current density
   **/
  tbox::Pointer<pdat::CellVariable<NDIM,double> >
    Jn_vec = new pdat::CellVariable<NDIM,double>("Jn_vec",3);
  tbox::Pointer<pdat::CellVariable<NDIM,double> >
    Jp_vec = new pdat::CellVariable<NDIM,double>("Jp_vec",3);
  tbox::Pointer<pdat::CellVariable<NDIM,double> >
    J_vec = new pdat::CellVariable<NDIM,double>("J_vec",3);

  tbox::Pointer<pdat::CellVariable<NDIM,double> >
    J_edge = new pdat::CellVariable<NDIM,double>("J_edge",cell_edge);

  /**
   * @brief cell/node to entity no. and material no. relationship
   **/
  tbox::Pointer<pdat::CellVariable<NDIM,int> >
    cell_entity = new pdat::CellVariable<NDIM,int>("cell_entity",1);
  tbox::Pointer<pdat::CellVariable<NDIM,int> >
    cell_mat = new pdat::CellVariable<NDIM,int>("cell_mat",1);
  tbox::Pointer<pdat::NodeVariable<NDIM,int> >
    node_entity = new pdat::NodeVariable<NDIM,int>("node_entity",2);
  tbox::Pointer<pdat::NodeVariable<NDIM,int> >
    node_mat = new pdat::NodeVariable<NDIM,int>("node_mat",n_mat);

  /**
   * @brief semi-region node mark
   **/
  tbox::Pointer<pdat::NodeVariable<NDIM,int> >
    semi_node = new pdat::NodeVariable<NDIM,int>("semi_node",1);

  /**
   * @brief metal eletrode region node mark
   **/
  tbox::Pointer<pdat::NodeVariable<NDIM,int> >
    metal_node = new pdat::NodeVariable<NDIM,int>("metal_node",1);

  /**
   * @brief sio2 region node mark
   **/
  tbox::Pointer<pdat::NodeVariable<NDIM,int> >
    sio2_node = new pdat::NodeVariable<NDIM,int>("sio2_node",1);

  /**
   * @brief siclass region node mark
   **/
  tbox::Pointer<pdat::NodeVariable<NDIM,int> >
    siclass_node = new pdat::NodeVariable<NDIM,int>("siclass_node",1);

  /**
   * @brief semi_insulator boundary node mark
   **/
  tbox::Pointer<pdat::NodeVariable<NDIM,int> >
    semiinsbc_node = new pdat::NodeVariable<NDIM,int>("semiinsbc_node",1);

  /**
   * @brief contact boundary node mark
   **/
  tbox::Pointer<pdat::NodeVariable<NDIM,int> >
    contact_node = new pdat::NodeVariable<NDIM,int>("contact_node",1);

  /**
   * @brief heterojunction node mark
   **/
  tbox::Pointer<pdat::NodeVariable<NDIM,int> >
    interface_node = new pdat::NodeVariable<NDIM,int>("interface_node",1);
  tbox::Pointer<pdat::NodeVariable<NDIM,int> >
    heter_node = new pdat::NodeVariable<NDIM,int>("heter_node",1);
  tbox::Pointer<pdat::NodeVariable<NDIM,int> >
    hnode_entity = new pdat::NodeVariable<NDIM,int>("hnode_entity",2);
  tbox::Pointer<pdat::NodeVariable<NDIM,int> >
    hnode_mat = new pdat::NodeVariable<NDIM,int>("hnode_mat",2);
  tbox::Pointer<pdat::CellVariable<NDIM,int> >
    cell_hnode = new pdat::CellVariable<NDIM,int>("cell_hnode",1);

  /**
   * @brief boundary face - adjacent entity
   **/
  tbox::Pointer<pdat::FaceVariable<NDIM,int> >
    face_entity = new pdat::FaceVariable<NDIM,int>("face_entity", 2);
  tbox::Pointer<pdat::FaceVariable<NDIM,int> >
    face_mat = new pdat::FaceVariable<NDIM,int>("face_mat", 2);

  /**
   * @brief dopant information
   **/
  tbox::Pointer<pdat::NodeVariable<NDIM,double> >
    dopant_node = new pdat::NodeVariable<NDIM,double>("dopant_node",1);
  tbox::Pointer<pdat::NodeVariable<NDIM,double> >
    maxdoping_node = new pdat::NodeVariable<NDIM,double>("maxdoping_node",1);
  tbox::Pointer<pdat::NodeVariable<NDIM,double> >
    ND_node = new pdat::NodeVariable<NDIM,double>("ND_node",1);
  tbox::Pointer<pdat::NodeVariable<NDIM,double> >
    NA_node = new pdat::NodeVariable<NDIM,double>("NA_node",1);

  tbox::Pointer<pdat::NodeVariable<NDIM,double> >
    sol_node = new pdat::NodeVariable<NDIM,double>("sol_node",4);

  /**
   * @brief polarized charge for heterojunction
   **/
  tbox::Pointer<pdat::NodeVariable<NDIM,double> >
    polarizedp_node = new pdat::NodeVariable<NDIM,double>("polarizedp_node",1);

  /**
   * @brief node variable for debug
   **/
  tbox::Pointer<pdat::NodeVariable<NDIM,double> >
    node_variable = new pdat::NodeVariable<NDIM,double>("node_variable",1);

  /**
   * @brief anti-scaled varibles
   **/
  tbox::Pointer<pdat::VectorVariable<NDIM,double> >
    p_sol_antiscal = new pdat::VectorVariable<NDIM,double>("p_sol_antiscal",d_dof_info);
  tbox::Pointer<pdat::VectorVariable<NDIM,double> >
    n_sol_antiscal = new pdat::VectorVariable<NDIM,double>("n_sol_antiscal",d_dof_info);
  tbox::Pointer<pdat::VectorVariable<NDIM,double> >
    h_sol_antiscal = new pdat::VectorVariable<NDIM,double>("h_sol_antiscal",d_dof_info);
  tbox::Pointer<pdat::VectorVariable<NDIM,double> >
    fermin_antiscal = new pdat::VectorVariable<NDIM,double>("fermin_antiscal",d_dof_info);
  tbox::Pointer<pdat::VectorVariable<NDIM,double> >
    fermih_antiscal = new pdat::VectorVariable<NDIM,double>("fermih_antiscal",d_dof_info);

  /**
   *@brief solution storage
   **/
  tbox::Pointer<pdat::NodeVariable<NDIM,double> >
    carn_sol_1 = new pdat::NodeVariable<NDIM,double>("carn_sol_1",2);
  tbox::Pointer<pdat::NodeVariable<NDIM,double> >
    fn_sol_1 = new pdat::NodeVariable<NDIM,double>("fn_sol_1",2);
  tbox::Pointer<pdat::NodeVariable<NDIM,double> >
    fn_sol_0 = new pdat::NodeVariable<NDIM,double>("fn_sol_0",2);
  tbox::Pointer<pdat::NodeVariable<NDIM,double> >
    carh_sol_1 = new pdat::NodeVariable<NDIM,double>("carh_sol_1",2);
  tbox::Pointer<pdat::NodeVariable<NDIM,double> >
    fh_sol_1 = new pdat::NodeVariable<NDIM,double>("fh_sol_1",2);
  tbox::Pointer<pdat::NodeVariable<NDIM,double> >
    fh_sol_0 = new pdat::NodeVariable<NDIM,double>("fh_sol_0",2);
  tbox::Pointer<pdat::NodeVariable<NDIM,double> >
    phi_sol_0 = new pdat::NodeVariable<NDIM,double>("phi_sol_0",1);
  tbox::Pointer<pdat::NodeVariable<NDIM,double> >
    phi_sol_1 = new pdat::NodeVariable<NDIM,double>("phi_sol_1",2);

  //bas vector
  tbox::Pointer<pdat::CellVariable<NDIM, double> > bas_integ = new pdat::CellVariable<NDIM, double>("bas_integ",cell_edge*cell_point);
  tbox::Pointer<pdat::CellVariable<NDIM, double> > bas_rhs = new pdat::CellVariable<NDIM, double>("bas_rhs",cell_point);

  /**
   * @brief register variables and context to the variable database
   **/

  /** current conservation **/
  d_curconserve_mat_id    = variable_db->registerVariableAndContext(CC_matrix, current,1);
  d_curconserve_fval_id   = variable_db->registerVariableAndContext(CC_rhs,current,1);
  d_curconserve_sol_id    = variable_db->registerVariableAndContext(cc_sol_vec, current,1);

  /** poison equation **/
  d_poisson_mat_id        = variable_db->registerVariableAndContext(P_matrix, current,1);
  d_poisson_fval_id       = variable_db->registerVariableAndContext(p_rhs, current,1);
  d_poisson_sol_id              = variable_db->registerVariableAndContext(p_sol, current,1);
  d_poisson_sol_cur_id    = variable_db->registerVariableAndContext(p_sol_vec, current,1);
  d_poisson_stetch_id     = variable_db->registerVariableAndContext(p_stetch, current,1);

  /** newton -- poison **/
  d_poisson_newton_mat_id = variable_db->registerVariableAndContext(P_newton_matrix, current, 1);
  d_poisson_newton_fval_id= variable_db->registerVariableAndContext(p_newton_rhs, current, 1);
  d_poisson_newton_sol_id = variable_db->registerVariableAndContext(p_newton_sol_vec, current, 1);

  /** newton -- Vogamma **/
  d_Vogamma_newton_mat_id = variable_db->registerVariableAndContext(Vogamma_newton_matrix, current, 1);
  d_Vogamma_newton_fval_id= variable_db->registerVariableAndContext(Vogamma_newton_rhs, current, 1);
  d_Vogamma_newton_sol_id = variable_db->registerVariableAndContext(Vogamma_newton_sol_vec, current, 1);
  d_Vogamma_sol_cur_id    = variable_db->registerVariableAndContext(Vogamma_sol_vec, current,1);

  d_Vogamma_plot_id       = variable_db->registerVariableAndContext(Vogamma_plot, current,1);
  d_Vogammaplus_plot_id   = variable_db->registerVariableAndContext(Vogammaplus_plot, current,1);
  d_VogammaHplus_plot_id  = variable_db->registerVariableAndContext(VogammaHplus_plot, current,1);
  d_VogammaH2plus_plot_id = variable_db->registerVariableAndContext(VogammaH2plus_plot, current,1);
  d_VogammaH_plot_id      = variable_db->registerVariableAndContext(VogammaH_plot, current,1);
  d_VogammaH2_plot_id     = variable_db->registerVariableAndContext(VogammaH2_plot, current,1);

  d_poisson_plot_id       = variable_db->registerVariableAndContext(p_plot, current,1);

  /** carrier transprot equation -- electron **/
  d_carn_mat_id        = variable_db->registerVariableAndContext(N_matrix, current,1);
  d_carn_fval_id       = variable_db->registerVariableAndContext(n_rhs,current,1);
  d_carn_sol_cur_id    = variable_db->registerVariableAndContext(n_sol_vec, current,1);
  d_carn_stetch_id     = variable_db->registerVariableAndContext(n_stetch, current,1);
  d_logn_sol_id        = variable_db->registerVariableAndContext(logn_sol_vec, current,1);

  d_carn_plot_id       = variable_db->registerVariableAndContext(n_plot, current,1);
  d_logn_plot_id       = variable_db->registerVariableAndContext(logn_plot, current,1);

  d_fermin_sol_id      = variable_db->registerVariableAndContext(fermin_vec, current,1);
  d_fermin_stetch_id   = variable_db->registerVariableAndContext(fermin_stetch, current,1);
  d_fermin_plot_id     = variable_db->registerVariableAndContext(fermin_plot, current,1);

  /** carrier transport equation -- hole **/
  d_carh_mat_id        = variable_db->registerVariableAndContext(H_matrix, current,1);
  d_carh_fval_id       = variable_db->registerVariableAndContext(h_rhs,current,1);
  d_carh_sol_cur_id    = variable_db->registerVariableAndContext(h_sol_vec, current,1);
  d_carh_stetch_id     = variable_db->registerVariableAndContext(h_stetch, current,1);
  d_logh_sol_id        = variable_db->registerVariableAndContext(logh_sol_vec, current,1);
  d_carh_plot_id       = variable_db->registerVariableAndContext(h_plot, current,1);
  d_logh_plot_id       = variable_db->registerVariableAndContext(logh_plot, current,1);

  d_fermih_sol_id      = variable_db->registerVariableAndContext(fermih_vec, current,1);
  d_fermih_stetch_id   = variable_db->registerVariableAndContext(fermih_stetch, current,1);
  d_fermih_plot_id     = variable_db->registerVariableAndContext(fermih_plot, current,1);

  /** unscaled solutions **/
  d_poisson_unscal_id    = variable_db->registerVariableAndContext(p_sol_antiscal,current,1);
  d_carn_unscal_id       = variable_db->registerVariableAndContext(n_sol_antiscal,current,1);
  d_carh_unscal_id       = variable_db->registerVariableAndContext(h_sol_antiscal,current,1);
  d_fermin_unscal_id     = variable_db->registerVariableAndContext(fermin_antiscal,current,1);
  d_fermih_unscal_id     = variable_db->registerVariableAndContext(fermih_antiscal,current,1);

  /** fields **/
  d_E_vec_id             = variable_db->registerVariableAndContext(E_vec, current,1);
  d_E_cell_id            = variable_db->registerVariableAndContext(E_cell, current,1);
  d_Jn_vec_id            = variable_db->registerVariableAndContext(Jn_vec, current,1);
  d_Jp_vec_id            = variable_db->registerVariableAndContext(Jp_vec, current,1);
  d_J_vec_id             = variable_db->registerVariableAndContext(J_vec, current,1);
  d_J_edge_id            = variable_db->registerVariableAndContext(J_edge, current,1);

  /** structure information **/
  d_cell_entity_id       = variable_db->registerVariableAndContext(cell_entity, current,1);
  d_cell_mat_id          = variable_db->registerVariableAndContext(cell_mat, current,1);
  d_node_entity_id       = variable_db->registerVariableAndContext(node_entity, current,1);
  d_node_mat_id          = variable_db->registerVariableAndContext(node_mat, current,1);
  d_semi_node_id         = variable_db->registerVariableAndContext(semi_node, current,1);
  d_metal_node_id        = variable_db->registerVariableAndContext(metal_node, current,1);
  d_sio2_node_id         = variable_db->registerVariableAndContext(sio2_node, current,1);
  d_siclass_node_id           = variable_db->registerVariableAndContext(siclass_node, current,1);
  d_contact_node_id      = variable_db->registerVariableAndContext(contact_node, current,1);

  d_interface_node_id    = variable_db->registerVariableAndContext(interface_node, current,1);
  d_heter_node_id        = variable_db->registerVariableAndContext(heter_node, current,1);
  d_hnode_entity_id      = variable_db->registerVariableAndContext(hnode_entity, current,1);
  d_hnode_mat_id         = variable_db->registerVariableAndContext(hnode_mat, current,1);
  d_cell_hnode_id        = variable_db->registerVariableAndContext(cell_hnode, current,1);
  d_semiinsbc_node_id    = variable_db->registerVariableAndContext(semiinsbc_node, current,1);
  d_face_entity_id       = variable_db->registerVariableAndContext(face_entity, current,1);
  d_face_mat_id          = variable_db->registerVariableAndContext(face_mat, current,1);

  /** dopant information **/
  d_dopant_id            = variable_db->registerVariableAndContext(dopant_node, current,1);
  d_maxdoping_id         = variable_db->registerVariableAndContext(maxdoping_node, current,1);
  d_polarizedp_id        = variable_db->registerVariableAndContext(polarizedp_node, current,1);

  d_ND_id                = variable_db->registerVariableAndContext(ND_node, current,1);
  d_NA_id                = variable_db->registerVariableAndContext(NA_node, current,1);

  d_sol_id            = variable_db->registerVariableAndContext(sol_node, current,1);

  /** debug **/
  d_node_var_plot_id     = variable_db->registerVariableAndContext(node_variable, current,1);

  /** solution storage **/
  d_carn_sol_1_id        = variable_db->registerVariableAndContext(carn_sol_1, current,1);
  d_fn_sol_1_id          = variable_db->registerVariableAndContext(fn_sol_1, current,1);
  d_fn_sol_0_id          = variable_db->registerVariableAndContext(fn_sol_0, current,1);
  d_carh_sol_1_id        = variable_db->registerVariableAndContext(carh_sol_1, current,1);
  d_fh_sol_1_id          = variable_db->registerVariableAndContext(fh_sol_1, current,1);
  d_fh_sol_0_id          = variable_db->registerVariableAndContext(fh_sol_0, current,1);
  d_phi_sol_0_id         = variable_db->registerVariableAndContext(phi_sol_0, current,1);
  d_phi_sol_1_id         = variable_db->registerVariableAndContext(phi_sol_1, current,1);

  d_bas_integ_id = variable_db->registerVariableAndContext(bas_integ,current,1);
  d_bas_rhs_id = variable_db->registerVariableAndContext(bas_rhs,current,1);

  tbox::pout << "***** At ending of registerModelVariable! *****"<< std::endl;
}

/** ***********************************************************************
 *  初始化指定的积分构件.
 *
 *  注册待填充的数据片或待调度内存空间的数据片到积分构件.
 *********************************************************************** **/
/**8.11 613 622**/
void NonlinearPoisson::initializeComponent(
        algs::IntegratorComponent<NDIM> * component) const {

  const string &component_name = component->getName();
  /** initialization component **/
  if ( component_name=="INIT") {
    component->registerInitPatchData(d_dopant_id);
    component->registerInitPatchData(d_maxdoping_id);
  }
  /** 8.11 **/
  else if ( component_name=="INIT_SOL"){
      component->registerInitPatchData(d_sol_id);
  }
  else if (component_name == "DOP_INIT"){
      component->registerInitPatchData(d_ND_id);
      component->registerInitPatchData(d_NA_id);
      if(enable_selfDop){
          component->registerInitPatchData(d_dopant_id);
          component->registerInitPatchData(d_maxdoping_id);
      }
  }
  else if(component_name == "TID_INIT"){
      component->registerInitPatchData(d_Vogamma_plot_id);
      component->registerInitPatchData(d_Vogammaplus_plot_id);
      component->registerInitPatchData(d_VogammaHplus_plot_id);
      component->registerInitPatchData(d_VogammaH2plus_plot_id);
      component->registerInitPatchData(d_VogammaH_plot_id);
      component->registerInitPatchData(d_VogammaH2_plot_id);
      d_dof_Vogamma_info->registerToInitComponent(component);
  }
  /** node or cell data **/
  else if (component_name == "DOF_INIT"){
      component->registerInitPatchData(d_poisson_plot_id);

      component->registerInitPatchData(d_carn_plot_id);
      component->registerInitPatchData(d_logn_plot_id);
      component->registerInitPatchData(d_fermin_plot_id);

      component->registerInitPatchData(d_carh_plot_id);
      component->registerInitPatchData(d_logh_plot_id);
      component->registerInitPatchData(d_fermih_plot_id);

      component->registerInitPatchData(d_node_var_plot_id);

      component->registerInitPatchData(d_E_vec_id);
      component->registerInitPatchData(d_E_cell_id);
      component->registerInitPatchData(d_Jn_vec_id);
      component->registerInitPatchData(d_Jp_vec_id);
      component->registerInitPatchData(d_J_vec_id);
      component->registerInitPatchData(d_J_edge_id);

      component->registerInitPatchData(d_polarizedp_id);

      component->registerInitPatchData(d_cell_entity_id);
      component->registerInitPatchData(d_cell_mat_id);
      component->registerInitPatchData(d_node_entity_id);
      component->registerInitPatchData(d_node_mat_id);
      component->registerInitPatchData(d_semi_node_id);
      component->registerInitPatchData(d_metal_node_id);
      component->registerInitPatchData(d_sio2_node_id);
      component->registerInitPatchData(d_siclass_node_id);
      component->registerInitPatchData(d_contact_node_id);

      component->registerInitPatchData(d_semiinsbc_node_id);

      component->registerInitPatchData(d_heter_node_id);
      component->registerInitPatchData(d_hnode_entity_id);
      component->registerInitPatchData(d_hnode_mat_id);
      component->registerInitPatchData(d_cell_hnode_id);

      component->registerInitPatchData(d_face_entity_id);
      component->registerInitPatchData(d_face_mat_id);

      component->registerInitPatchData(d_carn_sol_1_id);
      component->registerInitPatchData(d_fn_sol_1_id);
      component->registerInitPatchData(d_fn_sol_0_id);
      component->registerInitPatchData(d_carh_sol_1_id);
      component->registerInitPatchData(d_fh_sol_1_id);
      component->registerInitPatchData(d_fh_sol_0_id);
      component->registerInitPatchData(d_phi_sol_0_id);
      component->registerInitPatchData(d_phi_sol_1_id);

      if(save_bas){
          component->registerInitPatchData(d_bas_integ_id);
          component->registerInitPatchData(d_bas_rhs_id);
      }

      /** 将自由度信息中的若干数据片注册到初始化构件 **/
      d_dof_info->registerToInitComponent(component);
      d_dof_semi_info->registerToInitComponent(component);
  }

  /** allocate variable vectors **/
  else if ( component_name=="AllocVec") {
      component->registerPatchData(d_curconserve_sol_id);
      component->registerPatchData(d_poisson_sol_cur_id);
      component->registerPatchData(d_poisson_sol_id);
      component->registerPatchData(d_poisson_stetch_id);
      component->registerPatchData(d_poisson_newton_sol_id);

      component->registerPatchData(d_carn_sol_cur_id);
      component->registerPatchData(d_carn_stetch_id);
      component->registerPatchData(d_fermin_sol_id);
      component->registerPatchData(d_fermin_stetch_id);
      component->registerPatchData(d_logn_sol_id);

      component->registerPatchData(d_carh_sol_cur_id);
      component->registerPatchData(d_carh_stetch_id);
      component->registerPatchData(d_fermih_sol_id);
      component->registerPatchData(d_fermih_stetch_id);
      component->registerPatchData(d_logh_sol_id);

      component->registerPatchData(d_poisson_unscal_id);
      component->registerPatchData(d_carn_unscal_id);
      component->registerPatchData(d_carh_unscal_id);
      component->registerPatchData(d_fermin_unscal_id);
      component->registerPatchData(d_fermih_unscal_id);

      if(enable_tid){
          component->registerPatchData(d_Vogamma_sol_cur_id);
      }
  }

  /** initialization **/
  else if (component_name == "INIT_REBALANCE") {
      if(!enable_selfDop){
          component->registerPatchData(d_dopant_id);
          component->registerPatchData(d_maxdoping_id);
      }
      component->registerPatchData(d_sol_id);
  }

  /** sulotion vector initialization **/
  else if ( component_name=="SOLINIT") {
  }

  else if ( component_name=="CHECKORDER") {

  }

  /** get neutral status for initial **/
  else if (component_name=="EQUILIBRIUM"){
      component->registerCommunicationPatchData(d_poisson_sol_cur_id, d_poisson_sol_cur_id);
  }

  /** current conservation equation **/
  else if ( component_name=="AllocCCMat"){    /** allocate matrix **/
      component->registerPatchData(d_curconserve_mat_id);
  }
  else if ( component_name=="AllocCCRhs"){    /** allocate rhs **/
      component->registerPatchData(d_curconserve_fval_id);
  }
  else if ( component_name=="CCFX"){          /** rhs construction  **/

  }
  else if ( component_name=="CCMAT"){         /** matrix construction **/

  }
  else if ( component_name=="CCPHYBC"){       /** matrix **/

  }
  else if ( component_name=="CCAPPBC"){       /** rhs **/

  }
  /** Poisson Equation **/
  else if ( component_name=="AllocPMat") {    /** allocate matrix **/
      component->registerPatchData(d_poisson_mat_id);
  }
  else if (component_name =="AllocPRhs") {    /** allocate rhs **/
      component->registerPatchData(d_poisson_fval_id);
  }
  else if ( component_name=="PMAT") {         /** matrix construction **/
      component->registerCommunicationPatchData(d_poisson_sol_cur_id, d_poisson_sol_cur_id);
  }
  else if ( component_name=="PFX") {          /** rhs vector construction **/
      component->registerCommunicationPatchData(d_poisson_sol_cur_id, d_poisson_sol_cur_id);
  }
  else if ( component_name=="PPHYBC") {       /** matrix **/
      component->registerCommunicationPatchData(d_poisson_sol_cur_id, d_poisson_sol_cur_id);
  }
  else if ( component_name=="PAppBC") {       /** rhs **/
      component->registerCommunicationPatchData(d_poisson_sol_cur_id, d_poisson_sol_cur_id);
      component->registerCommunicationPatchData(d_curconserve_sol_id, d_curconserve_sol_id);
  }

  /** Newton-Raphson for Poison **/
  else if ( component_name=="AllocNewtonPMat") { /** allocate matrix **/
      component->registerPatchData(d_poisson_newton_mat_id);
  }
  else if (component_name =="AllocNewtonPRhs") { /** allocate rhs **/
      component->registerPatchData(d_poisson_newton_fval_id);
  }
  /** Newton-Raphson for Vogamma **/
  else if ( component_name=="AllocNewtonVogammaMat") { /** allocate matrix **/
      component->registerPatchData(d_Vogamma_newton_mat_id);
  }
  else if (component_name =="AllocNewtonVogammaRhs") { /** allocate rhs **/
      component->registerPatchData(d_Vogamma_newton_fval_id);
  }
  else if ( component_name=="ADTEST"){

  }
  else if ( component_name=="NewtonPMAT") {      /** matrix construction **/
      component->registerCommunicationPatchData(d_poisson_sol_cur_id, d_poisson_sol_cur_id);
      component->registerCommunicationPatchData(d_fermin_sol_id, d_fermin_sol_id);
      component->registerCommunicationPatchData(d_fermih_sol_id, d_fermih_sol_id);
  }
  else if ( component_name=="NewtonPFX"){        /** rhs construction **/
      component->registerCommunicationPatchData(d_poisson_sol_cur_id, d_poisson_sol_cur_id);
  }
  else if ( component_name=="NewtonPPHYBC") {    /** matrix **/
      component->registerCommunicationPatchData(d_poisson_sol_cur_id, d_poisson_sol_cur_id);
  }
  else if ( component_name=="NewtonPAppBC") {    /** rhs **/
      component->registerCommunicationPatchData(d_curconserve_sol_id, d_curconserve_sol_id);
      component->registerCommunicationPatchData(d_poisson_sol_cur_id, d_poisson_sol_cur_id);
  }

  /** carrier transport equation **/
  else if ( component_name=="AllocCMat") {     /** allocate matrix **/
      component->registerPatchData(d_carn_mat_id);
      component->registerPatchData(d_carh_mat_id);
  }
  else if ( component_name=="AllocCRhs") {     /** allocate rhs **/
      component->registerPatchData(d_carn_fval_id);
      component->registerPatchData(d_carh_fval_id);
  }
  else if (component_name == "CVMAT"){         /** matrix construction for both electron & hole **/
      component->registerCommunicationPatchData(d_poisson_sol_cur_id, d_poisson_sol_cur_id);
      component->registerCommunicationPatchData(d_fermin_sol_id, d_fermin_sol_id);
      component->registerCommunicationPatchData(d_fermih_sol_id, d_fermih_sol_id);
  }
  /** electron **/
  else if ( component_name=="NFX") {           /** rhs construction **/
      component->registerCommunicationPatchData(d_poisson_sol_cur_id, d_poisson_sol_cur_id);
      component->registerCommunicationPatchData(d_fermin_sol_id, d_fermin_sol_id);
      component->registerCommunicationPatchData(d_fermih_sol_id, d_fermih_sol_id);
  }
  else if ( component_name=="NPHYBC") {        /** matrix **/

  }
  else if ( component_name=="NAppBC") {        /** rhs **/

  }
  /** hole **/
  else if ( component_name=="HFX") {           /** rhs construction **/
      component->registerCommunicationPatchData(d_poisson_sol_cur_id, d_poisson_sol_cur_id);
      component->registerCommunicationPatchData(d_fermin_sol_id, d_fermin_sol_id);
      component->registerCommunicationPatchData(d_fermih_sol_id, d_fermih_sol_id);
  }
  else if ( component_name=="HPHYBC") {        /** matrix **/

  }
  else if ( component_name=="HAppBC") {        /** rhs **/

  }

  /** electric field & current density & Joule heat generation **/
  else if ( component_name=="JHEAT") {         /** calculate joule heating **/
      component->registerCommunicationPatchData(d_poisson_sol_cur_id, d_poisson_sol_cur_id);
      component->registerCommunicationPatchData(d_carn_sol_cur_id, d_carn_sol_cur_id);
      component->registerCommunicationPatchData(d_carh_sol_cur_id, d_carh_sol_cur_id);
  }
  else if ( component_name=="CURRDENS") {      /** calculate current density **/
      component->registerCommunicationPatchData(d_poisson_sol_cur_id, d_poisson_sol_cur_id);
      component->registerCommunicationPatchData(d_carn_sol_cur_id, d_carn_sol_cur_id);
      component->registerCommunicationPatchData(d_carh_sol_cur_id, d_carh_sol_cur_id);
  }
  else if ( component_name=="TFX") {           /** rhs constrcution **/
      component->registerCommunicationPatchData(d_poisson_sol_cur_id, d_poisson_sol_cur_id);
      component->registerCommunicationPatchData(d_carn_sol_cur_id, d_carn_sol_cur_id);
      component->registerCommunicationPatchData(d_carh_sol_cur_id, d_carh_sol_cur_id);
      component->registerCommunicationPatchData(d_fermin_sol_id, d_fermin_sol_id);
      component->registerCommunicationPatchData(d_fermih_sol_id, d_fermih_sol_id);
  }
  else if ( component_name=="TPHYBC") {        /** matrix **/

  }

  /** semiconductor related component **/
  else if ( component_name=="FINTRIN") {       /** calculate intrinsic fermi level **/

  }
  else if ( component_name=="FLEVEL") {        /** calculate fermi level **/
      component->registerCommunicationPatchData(d_carn_sol_cur_id, d_carn_sol_cur_id);
      component->registerCommunicationPatchData(d_carh_sol_cur_id, d_carh_sol_cur_id);
      component->registerCommunicationPatchData(d_poisson_sol_cur_id, d_poisson_sol_cur_id);
  }
  else if ( component_name=="NPCONS") {
      component->registerCommunicationPatchData(d_carn_sol_cur_id, d_carn_sol_cur_id);
      component->registerCommunicationPatchData(d_carh_sol_cur_id, d_carh_sol_cur_id);
      component->registerCommunicationPatchData(d_poisson_sol_cur_id, d_poisson_sol_cur_id);
  }

  /** transport newly calculated solution for storage **/
  else if ( component_name=="PSOLTRANS") {      /** potential **/
      component->registerCommunicationPatchData(d_poisson_newton_sol_id, d_poisson_newton_sol_id);
      component->registerCommunicationPatchData(d_poisson_sol_cur_id, d_poisson_sol_cur_id);
  }
  else if ( component_name=="GummelPSOLTRANS") {      /** potential **/
      component->registerCommunicationPatchData(d_poisson_sol_id, d_poisson_sol_id);
      component->registerCommunicationPatchData(d_poisson_sol_cur_id, d_poisson_sol_cur_id);
  }
  else if ( component_name=="SEMISOLTRANS") {   /** fermi level **/
      component->registerCommunicationPatchData(d_fermin_sol_id, d_fermin_sol_id);
      component->registerCommunicationPatchData(d_fermih_sol_id, d_fermih_sol_id);
      component->registerCommunicationPatchData(d_carn_sol_cur_id, d_carn_sol_cur_id);
      component->registerCommunicationPatchData(d_carh_sol_cur_id, d_carh_sol_cur_id);
      component->registerCommunicationPatchData(d_poisson_sol_cur_id, d_poisson_sol_cur_id);
  }
  else if ( component_name=="GummelSEMISOLTRANS") {   /** fermi level **/
      component->registerCommunicationPatchData(d_fermin_sol_id, d_fermin_sol_id);
      component->registerCommunicationPatchData(d_fermih_sol_id, d_fermih_sol_id);
      component->registerCommunicationPatchData(d_carn_sol_cur_id, d_carn_sol_cur_id);
      component->registerCommunicationPatchData(d_carh_sol_cur_id, d_carh_sol_cur_id);
      component->registerCommunicationPatchData(d_poisson_sol_cur_id, d_poisson_sol_cur_id);
  }
  else if ( component_name=="SSOLTRANS") {     /** system -- potential, carrier (fermi) & temperature **/

  }

  /** maximum error searching **/
  else if ( component_name=="POISSON") {       /** potential **/

  }
  else if ( component_name=="GummelPOISSON") {       /** potential **/

  }
  else if ( component_name=="THERMAL") {       /** temperature **/

  }
  else if ( component_name=="SYSTEM") {        /** system--temperature **/

  }
  else if ( component_name=="SEMI") {          /** DDM--fermi level **/

  }
  else if ( component_name=="GummelSEMI") {          /** DDM--fermi level **/

  }
  else if ( component_name=="get_I") {

  }
  else if ( component_name=="saveSol") {

  }

  /** solution transportation to plot data **/
  else if ( component_name=="PPOST") {          /** potential **/
      component->registerCommunicationPatchData(d_poisson_sol_cur_id, d_poisson_sol_cur_id);
  }
  else if ( component_name=="NPOST") {          /** electron density, fermi level **/
      component->registerCommunicationPatchData(d_carn_sol_cur_id, d_carn_sol_cur_id);
      component->registerCommunicationPatchData(d_poisson_sol_cur_id, d_poisson_sol_cur_id);
  }
  else if ( component_name=="HPOST") {          /** hole density, fermi level **/
      component->registerCommunicationPatchData(d_carh_sol_cur_id, d_carh_sol_cur_id);
      component->registerCommunicationPatchData(d_poisson_sol_cur_id, d_poisson_sol_cur_id);
  }
  else if ( component_name=="calculate_E") {
      component->registerCommunicationPatchData(d_poisson_sol_cur_id,d_poisson_sol_cur_id);
  }
  else if ( component_name=="calculate_J") {
      component->registerCommunicationPatchData(d_poisson_sol_cur_id,d_poisson_sol_cur_id);
      component->registerCommunicationPatchData(d_carn_sol_cur_id,d_carn_sol_cur_id);
      component->registerCommunicationPatchData(d_carh_sol_cur_id,d_carh_sol_cur_id);
  }
  else if ( component_name=="get_IV"){
      component->registerCommunicationPatchData(d_poisson_sol_cur_id,d_poisson_sol_cur_id);
      component->registerCommunicationPatchData(d_carn_sol_cur_id,d_carn_sol_cur_id);
      component->registerCommunicationPatchData(d_carh_sol_cur_id,d_carh_sol_cur_id);
  }
  else if ( component_name=="MAX_DOP"){

  }
  /** get time step for transient **/
  else if ( component_name=="Dt") {

  }

  /** step solution for transient **/
  else if ( component_name=="SOLPROCESS") {     /** transport solution for next time step **/

  }

  /** scaling && anti scaling process **/
  else if ( component_name=="SCALOR") {         /** scaling variables **/

  }
  else if ( component_name=="ANTISCALOR") {     /** anti scaling variables **/

  }
  else if ( component_name=="ANTICARRIER") {    /** anti scaling carrier density **/

  }
  /** other matrix **/
  else if ( component_name=="MV") {

  }
  else {
    TBOX_ERROR("\n::initializeComponent() : component "
               << component_name <<" is not matched. "<<endl);
  }

}

/** ***********************************************************************
 *  初始化数据片（支持初值构件）
 *********************************************************************** **/
void NonlinearPoisson::initializePatchData(hier::Patch<NDIM> & patch,
                                           const double  time,
                                           const bool    initial_time,
                                           const string& component_name) {

    NULL_USE(time);             /**< 初始化中没有用到time */

    /** 8.11 **/
    if(component_name == "INIT_SOL"){
        tbox::Pointer< hier::PatchTopology<NDIM> > patch_top = patch.getPatchTopology();

        tbox::Pointer< hier::PatchGeometry<NDIM> > patch_geo = patch.getPatchGeometry();

        int num_nodes = patch.getNumberOfNodes(0);

        tbox::Pointer<pdat::NodeData<NDIM, double> >
                sol_data = patch.getPatchData(d_sol_id);

        for(int index_node=0;index_node<num_nodes;index_node++){
            (*sol_data)(3,index_node) = index_node;
        }

        if(get_InitialSol){
            for (int i = 0; i < num_sol_file_name; ++i) {
                stringstream index_sol_file_name;
                index_sol_file_name<<sol_file_name<<"_"<<i;
                SolRead(patch,index_sol_file_name.str());
            }
        }
        if(get_PoissonSolcur){
            for (int i = 0; i < num_PoissonSolcur_filename; ++i) {
                stringstream index_PoissonSolcur_filename;
                if(i<10)
                    index_PoissonSolcur_filename<<PoissonSolcur_filename<<"_"<<"0000000"<<i;
                else if(i<100&&i>=10)
                    index_PoissonSolcur_filename<<PoissonSolcur_filename<<"_"<<"000000"<<i;
                PoissonSolcurRead(patch,index_PoissonSolcur_filename.str());
            }
        }
    }

    else if(component_name == "INIT"){

        if (initial_time) {
            tbox::Pointer<pdat::NodeData<NDIM, double> >
                    dopant_data = patch.getPatchData(d_dopant_id);
            tbox::Pointer<pdat::NodeData<NDIM, double> >
                    maxdoping_data = patch.getPatchData(d_maxdoping_id);

            double *dopant = dopant_data->getPointer();
            double *max_doping = maxdoping_data->getPointer();

            tbox::pout<<"dopant_file_name: "<<dopant_file_name<<endl;
            std::string dop_dir(dopant_file_name);

            SerialRead(dop_dir,patch,dopant);   /** read dopant information **/

            /** get the max doping information**/
            Get_MaxDop(dop_dir,patch,max_doping);
        }
    }else if(component_name == "DOP_INIT"){
        if(initial_time){
            /** Get Local PatchTopology **/
            tbox::Pointer< hier::PatchTopology<NDIM> > patch_top = patch.getPatchTopology();

            tbox::Pointer< hier::PatchGeometry<NDIM> > patch_geo = patch.getPatchGeometry();
            tbox::Pointer<pdat::NodeData<NDIM, double> >
                    dopant_data = patch.getPatchData(d_dopant_id);
            tbox::Pointer<pdat::NodeData<NDIM, double> >
                    maxdoping_data = patch.getPatchData(d_maxdoping_id);

            tbox::Pointer<pdat::NodeData<NDIM, double> >
                    ND_data = patch.getPatchData(d_ND_id);
            tbox::Pointer<pdat::NodeData<NDIM, double> >
                    NA_data = patch.getPatchData(d_NA_id);

            /** get Patch node coord **/
            tbox::Pointer<pdat::NodeData<NDIM, double> >
                    node_coord = patch_geo->getNodeCoordinates();
            int num_nodes = patch.getNumberOfNodes(1);
            if(enable_selfDop){
                /** calculate dopant information **/
                tbox::Pointer<Variables_Functions> Var_func = new Variables_Functions();
                Var_func->calculateDop(node_coord,num_nodes,ND_data,NA_data,Dopassignment);
            }
            for (int i = 0; i < num_nodes; ++i){
                if(enable_selfDop)
                    dopant_data->getPointer()[i] = ((*ND_data)(0,i)-(*NA_data)(0,i));
            }
        }
    }
    else if (component_name == "TID_INIT"){
        if(initial_time){
          /** Get Local PatchTopology **/
          tbox::Pointer< hier::PatchTopology<NDIM> > patch_top = patch.getPatchTopology();

          tbox::Pointer< hier::PatchGeometry<NDIM> > patch_geo = patch.getPatchGeometry();

          /** get Patch node coord **/
          tbox::Pointer<pdat::NodeData<NDIM, double> >
                  node_coord = patch_geo->getNodeCoordinates();

          tbox::Pointer<pdat::NodeData<NDIM,double> >
                  Vogamma_plot_data = patch.getPatchData(d_Vogamma_plot_id);
          tbox::Pointer<pdat::NodeData<NDIM,double> >
                  Vogammaplus_plot_data = patch.getPatchData(d_Vogammaplus_plot_id);
          tbox::Pointer<pdat::NodeData<NDIM,double> >
                  VogammaHplus_plot_data = patch.getPatchData(d_VogammaHplus_plot_id);
          tbox::Pointer<pdat::NodeData<NDIM,double> >
                  VogammaH2plus_plot_data = patch.getPatchData(d_VogammaH2plus_plot_id);
          tbox::Pointer<pdat::NodeData<NDIM,double> >
                  VogammaH_plot_data = patch.getPatchData(d_VogammaH_plot_id);
          tbox::Pointer<pdat::NodeData<NDIM,double> >
                  VogammaH2_plot_data = patch.getPatchData(d_VogammaH2_plot_id);

          int num_nodes = patch.getNumberOfNodes(1);
          int num_cells = patch.getNumberOfCells(1);

          tbox::Pointer<pdat::CellData<NDIM, int> >
                  cell_mat_data = patch.getPatchData(d_cell_mat_id);

          /** material list **/
          tbox::Pointer<MaterialManager<NDIM> >
                      mat_manager = MaterialManager<NDIM>::getManager();

          /** dof information for Vogamma**/
          int * dis_Vogamma_ptr = d_dof_Vogamma_info->getDOFDistribution(patch, hier::EntityUtilities::NODE);
          for (int i = 0; i < num_nodes; ++i){
              dis_Vogamma_ptr[i] = 0;
          }

          /** dof information for semi**/
          tbox::Array<int> can_extent, can_indices;
          patch_top->getCellAdjacencyNodes(can_extent, can_indices);

          for(int i_e = 0; i_e < num_cells; ++i_e){
              if(material_list[(*cell_mat_data)(0,i_e)].mat_type == "dielectric"){
                  for(int i_n = 0; i_n < cell_point; ++i_n){
                      int tmp_num = can_indices[can_extent[i_e] + i_n];
                      dis_Vogamma_ptr[tmp_num] = 6;
                  }
              }
          }

          /** build dof mapping **/
          d_dof_Vogamma_info->buildPatchDOFMapping(patch);

          for(int i = 0; i < num_nodes; ++i){
              (*Vogamma_plot_data)(0,i) = 1e20;
              (*Vogammaplus_plot_data)(0,i) = 0;
              (*VogammaHplus_plot_data)(0,i) = 0;
              (*VogammaH2plus_plot_data)(0,i) = 0;
              (*VogammaH_plot_data)(0,i) = 1e20;
              (*VogammaH2_plot_data)(0,i) = 1e20;
          }
      }
    }
    else if (component_name == "DOF_INIT"){
        if(initial_time){
          /** Get Local PatchTopology **/
          tbox::Pointer< hier::PatchTopology<NDIM> > patch_top = patch.getPatchTopology();

          tbox::Pointer< hier::PatchGeometry<NDIM> > patch_geo = patch.getPatchGeometry();

          /** get Patch node coord **/
          tbox::Pointer<pdat::NodeData<NDIM, double> >
                  node_coord = patch_geo->getNodeCoordinates();

          tbox::Pointer<pdat::NodeData<NDIM,double> >
                  poisson_plot_data = patch.getPatchData(d_poisson_plot_id);

          tbox::Pointer<pdat::NodeData<NDIM,double> >
                  carn_plot_data = patch.getPatchData(d_carn_plot_id);
          tbox::Pointer<pdat::NodeData<NDIM,double> >
                  logn_plot_data = patch.getPatchData(d_logn_plot_id);
          tbox::Pointer<pdat::NodeData<NDIM,double> >
                  fermin_plot_data = patch.getPatchData(d_fermin_plot_id);

          tbox::Pointer<pdat::NodeData<NDIM,double> >
                  carh_plot_data = patch.getPatchData(d_carh_plot_id);
          tbox::Pointer<pdat::NodeData<NDIM,double> >
                  logh_plot_data = patch.getPatchData(d_logh_plot_id);
          tbox::Pointer<pdat::NodeData<NDIM,double> >
                  fermih_plot_data = patch.getPatchData(d_fermih_plot_id);

          tbox::Pointer<pdat::NodeData<NDIM,double> >
                  node_plot_data = patch.getPatchData(d_node_var_plot_id);

          tbox::Pointer<pdat::CellData<NDIM, double> >
                  e_cell_data =patch.getPatchData(d_E_cell_id);
          tbox::Pointer<pdat::CellData<NDIM, double> >
                  e_field_data = patch.getPatchData(d_E_vec_id);
          tbox::Pointer<pdat::CellData<NDIM, double> >
                  jn_field_data = patch.getPatchData(d_Jn_vec_id);
          tbox::Pointer<pdat::CellData<NDIM, double> >
                  jp_field_data = patch.getPatchData(d_Jp_vec_id);
          tbox::Pointer<pdat::CellData<NDIM, double> >
                  j_field_data = patch.getPatchData(d_J_vec_id);
          tbox::Pointer<pdat::CellData<NDIM, double> >
                  j_edge_data = patch.getPatchData(d_J_edge_id);

          tbox::Pointer<pdat::NodeData<NDIM, double> >
                  pp_data = patch.getPatchData(d_polarizedp_id);
          tbox::Pointer<pdat::NodeData<NDIM, double> >
                  dopant_data = patch.getPatchData(d_dopant_id);

          tbox::Pointer<pdat::CellData<NDIM, int> >
                  cell_entity_data = patch.getPatchData(d_cell_entity_id);
          tbox::Pointer<pdat::CellData<NDIM, int> >
                  cell_mat_data = patch.getPatchData(d_cell_mat_id);
          tbox::Pointer<pdat::NodeData<NDIM, int> >
                  node_entity_data = patch.getPatchData(d_node_entity_id);
          tbox::Pointer<pdat::NodeData<NDIM, int> >
                  node_mat_data = patch.getPatchData(d_node_mat_id);
          tbox::Pointer<pdat::NodeData<NDIM, int> >
                  semi_node_data = patch.getPatchData(d_semi_node_id);
          tbox::Pointer<pdat::NodeData<NDIM, int> >
                  metal_node_data = patch.getPatchData(d_metal_node_id);
          tbox::Pointer<pdat::NodeData<NDIM, int> >
                  sio2_node_data = patch.getPatchData(d_sio2_node_id);
          tbox::Pointer<pdat::NodeData<NDIM, int> >
                  siclass_node_data = patch.getPatchData(d_siclass_node_id);
          tbox::Pointer<pdat::NodeData<NDIM, int> >
                  contact_node_data = patch.getPatchData(d_contact_node_id);
          tbox::Pointer<pdat::NodeData<NDIM, int> >
                  semiinsbc_node_data = patch.getPatchData(d_semiinsbc_node_id);
          tbox::Pointer<pdat::NodeData<NDIM, int> >
                  heter_node_data = patch.getPatchData(d_heter_node_id);
          tbox::Pointer<pdat::NodeData<NDIM, int> >
                  hnode_entity_data = patch.getPatchData(d_hnode_entity_id);
          tbox::Pointer<pdat::NodeData<NDIM, int> >
                  hnode_mat_data = patch.getPatchData(d_hnode_mat_id);
          tbox::Pointer<pdat::CellData<NDIM, int> >
                  cell_hnode_data = patch.getPatchData(d_cell_hnode_id);
          tbox::Pointer<pdat::FaceData<NDIM, int> >
                  face_entity_data = patch.getPatchData(d_face_entity_id);
          tbox::Pointer<pdat::FaceData<NDIM, int> >
                  face_mat_data = patch.getPatchData(d_face_mat_id);

          tbox::Pointer<pdat::NodeData<NDIM, double> >
                  carn_sol_1_data = patch.getPatchData(d_carn_sol_1_id);
          tbox::Pointer<pdat::NodeData<NDIM, double> >
                  fn_sol_1_data = patch.getPatchData(d_fn_sol_1_id);
          tbox::Pointer<pdat::NodeData<NDIM, double> >
                  fn_sol_0_data = patch.getPatchData(d_fn_sol_0_id);
          tbox::Pointer<pdat::NodeData<NDIM, double> >
                  carh_sol_1_data = patch.getPatchData(d_carh_sol_1_id);
          tbox::Pointer<pdat::NodeData<NDIM, double> >
                  fh_sol_1_data = patch.getPatchData(d_fh_sol_1_id);
          tbox::Pointer<pdat::NodeData<NDIM, double> >
                  fh_sol_0_data = patch.getPatchData(d_fh_sol_0_id);
          tbox::Pointer<pdat::NodeData<NDIM, double> >
                  phi_sol_0_data = patch.getPatchData(d_phi_sol_0_id);
          tbox::Pointer<pdat::NodeData<NDIM, double> >
                  phi_sol_1_data = patch.getPatchData(d_phi_sol_1_id);

          int num_nodes = patch.getNumberOfNodes(1);
          int num_cells = patch.getNumberOfCells(1);
          int num_faces = patch.getNumberOfFaces(1);

          /** material list **/
          tbox::Pointer<MaterialManager<NDIM> >
                      mat_manager = MaterialManager<NDIM>::getManager();

          /** initial for node or cell variables **/
          for (int i = 0; i < num_nodes; ++i){
              (*node_entity_data)(0,i) = 0;
              (*node_entity_data)(1,i) = 0;
              (*semi_node_data)(0,i) = 0;
              (*metal_node_data)(0,i) = 0;
              (*sio2_node_data)(0,i) = 0;
              (*siclass_node_data)(0,i) = 0;
              (*contact_node_data)(0,i) = 0;
              (*semiinsbc_node_data)(0,i) = 0;
              (*heter_node_data)(0,i) = 0;
              (*hnode_entity_data)(0,i) = 0;
              (*hnode_entity_data)(1,i) = 0;
              (*hnode_mat_data)(0,i) = 0;
              (*hnode_mat_data)(1,i) = 0;
              (*pp_data)(0,i) = 0;
              for(int i_m = 0; i_m < n_mat; ++i_m){
                  (*node_mat_data)(i_m, i) = 0;
              }
          }
          for (int i = 0; i < num_cells; ++i){
              (*cell_entity_data)(0,i) = 0;
              (*cell_mat_data)(0,i) = 0; //0
              (*cell_hnode_data)(0,i) = 0;
          }

          for (int i = 0; i < num_faces; ++i){
              (*face_entity_data)(0,i) = -1;
              (*face_entity_data)(1,i) = -1;
              (*face_mat_data)(0,i) = -1;
              (*face_mat_data)(1,i) = -1;
          }

          /** cell - entity // node - entity relationship **/
          getCellEntityRelation(patch, component_name);
          getCellMatRelation(patch, component_name);
          getNodeEntityRelation(patch, component_name);
          getNodeMatRelation(patch, component_name);

          /** dof information **/
          int * dis_ptr = d_dof_info->getDOFDistribution(patch, hier::EntityUtilities::NODE);
          for (int i = 0; i < num_nodes; ++i){
              dis_ptr[i] = 0;
          }

          /** dof information for semi**/
          tbox::Array<int> can_extent, can_indices;
          patch_top->getCellAdjacencyNodes(can_extent, can_indices);
          int * dis_semi_ptr = d_dof_semi_info->getDOFDistribution(patch, hier::EntityUtilities::NODE);
          for(int i = 0; i < num_nodes; ++i){
              dis_semi_ptr[i] = 0;
          }

          for(int i_e = 0; i_e < num_cells; ++i_e){
              if(material_list[(*cell_mat_data)(0,i_e)].mat_type == "semiconductor"){
                  for(int i_n = 0; i_n < cell_point; ++i_n){
                      int tmp_num = can_indices[can_extent[i_e] + i_n];
                      dis_semi_ptr[tmp_num] = 3;
                      dis_ptr[tmp_num] = 1;
                      (*semi_node_data)(0,tmp_num) = 1;
                      (*siclass_node_data)(0,tmp_num) = 1;
                  }
              }else if(material_list[(*cell_mat_data)(0,i_e)].mat_type == "metal"){
                  for(int i_n = 0; i_n < cell_point; ++i_n){
                      int tmp_num = can_indices[can_extent[i_e] + i_n];
                      (*metal_node_data)(0,tmp_num) = 1;
                  }
              }else if(material_list[(*cell_mat_data)(0,i_e)].mat_type == "dielectric"){
                  for(int i_n = 0; i_n < cell_point; ++i_n){
                      int tmp_num = can_indices[can_extent[i_e] + i_n];
                      dis_semi_ptr[tmp_num] = 3;
                      dis_ptr[tmp_num] = 1;
                      (*sio2_node_data)(0,tmp_num) = 1;
                      (*siclass_node_data)(0,tmp_num) = 1;
                  }
              }
          }

#if 0
          for(int i_e = 0; i_e < num_cells; ++i_e){
              for(int i_n = 0; i_n < cell_point; ++i_n){
                  int tmp_num = can_indices[can_extent[i_e] + i_n];
                  if(((*semi_node_data)(0,tmp_num) == 1)&&((*sio2_node_data)(0,tmp_num) == 1)){
                      dis_semi_ptr[tmp_num] = 6;
                      dis_ptr[tmp_num] = 2;
                      (*semiinsbc_node_data)(0,tmp_num) = 1;
                  }
              }
          }

#endif
#if 0
          int num_node = patch.getNumberOfNodes(1);
          for (int i = 0; i < num_node; ++i){
              if(((*semi_node_data)(0,i) == 1)&&((*sio2_node_data)(0,i) == 1)){
                  dis_semi_ptr[i] = 6;
                  dis_ptr[i] = 2;
                  (*semiinsbc_node_data)(0,i) = 1;
              }
          }
#endif

          /** assigned dof for inner boundary **/
          int num_bc = EBoundassignment.getSize();
          for(int i_bc=0; i_bc<num_bc; ++i_bc){
              if (EBoundassignment[i_bc].bc_type == "heterojunction"){
                  int bound_no = EBoundassignment[i_bc].bc_no;
                  if (patch_geo->hasEntitySet(bound_no, hier::EntityUtilities::NODE)){
                      const tbox::Array<int> &entity_idx = patch_geo->getEntityIndicesInSet(bound_no,
                                                           hier::EntityUtilities::NODE,
                                                           num_nodes);
                      int size = entity_idx.getSize();
                      for (int i_n = 0; i_n < size; ++i_n){
                          dis_semi_ptr[entity_idx[i_n]] = 2;
                          dis_ptr[entity_idx[i_n]] = 2;
                      }
                  }
              }
          }

          /** get the contact id boundary node list and marked **/
          /** contact id boundary mark **/

          for(int i_bc=0; i_bc<num_bc; ++i_bc){
              int bound_no = EBoundassignment[i_bc].bc_no;
              if (patch_geo->hasEntitySet(bound_no, hier::EntityUtilities::NODE)){
                  const tbox::Array<int> &entity_idx = patch_geo->getEntityIndicesInSet(bound_no,
                                                       hier::EntityUtilities::NODE,
                                                       num_nodes);
                  int size = entity_idx.getSize();
                  for(int i_n = 0; i_n < size; ++i_n){
                      (*contact_node_data)(0,entity_idx[i_n]) = bound_no;
                  }
              }
          }

          /** get semiconductor-insulator boundary node list and marked **/
          /** heter boundary node mark**/

          for(int i_bc=0; i_bc<num_bc; ++i_bc){
              if (EBoundassignment[i_bc].bc_type == "semiinsbc"){
                  int bound_no = EBoundassignment[i_bc].bc_no;
                  if (patch_geo->hasEntitySet(bound_no, hier::EntityUtilities::NODE)){
                      const tbox::Array<int> &entity_idx = patch_geo->getEntityIndicesInSet(bound_no,
                                                           hier::EntityUtilities::NODE,
                                                           num_nodes);
                      int size = entity_idx.getSize();
                      for (int i_n = 0; i_n < size; ++i_n) {
                          (*semiinsbc_node_data)(0,entity_idx[i_n]) = 1;
                          dis_ptr[entity_idx[i_n]] = 2;
                          dis_semi_ptr[entity_idx[i_n]] = 6;
                      }
                  }
              }
              else if (EBoundassignment[i_bc].bc_type == "heterojunction"){
                  int bound_no = EBoundassignment[i_bc].bc_no;
                  if (patch_geo->hasEntitySet(bound_no, hier::EntityUtilities::NODE)){
                      const tbox::Array<int> &entity_idx = patch_geo->getEntityIndicesInSet(bound_no,
                                                           hier::EntityUtilities::NODE,
                                                           num_nodes);
                      int size = entity_idx.getSize();
                      for (int i_n = 0; i_n < size; ++i_n) {
                          (*pp_data)(0,entity_idx[i_n]) = EBoundassignment[i_bc].bc_value;

                          (*dopant_data)(0,entity_idx[i_n]) = EBoundassignment[i_bc].bc_value;
                          (*heter_node_data)(0,entity_idx[i_n]) = 1;
                      }
                  }
              }
          }

          /** build dof mapping **/
          d_dof_info->buildPatchDOFMapping(patch);
          d_dof_semi_info->buildPatchDOFMapping(patch);

          /** face - entity relationship: larger affinity write first **/
          tbox::Array<int> face_can_ext, face_can_idx;
          patch_top->getFaceAdjacencyCells(face_can_ext, face_can_idx);
          for (int i_face = 0; i_face < num_faces; ++i_face){
              int tmp_n = face_can_ext[i_face + 1] - face_can_ext[i_face];
              if (tmp_n == 1){
                  int tmp_cell = face_can_idx[face_can_ext[i_face]];
                  if (material_list[(*cell_mat_data)(0,tmp_cell)].mat_type == "semiconductor"){
                      (*face_entity_data)(0,i_face) = (*cell_entity_data)(0,tmp_cell);
                      (*face_mat_data)(0,i_face) = (*cell_mat_data)(0,tmp_cell);
                  }else{
                      (*face_entity_data)(1,i_face) = (*cell_entity_data)(0,tmp_cell);
                      (*face_mat_data)(1,i_face) = (*cell_mat_data)(0,tmp_cell);
                  }
              }else if (tmp_n == 2){
                  int tmp_cell_0 = face_can_idx[face_can_ext[i_face]];
                  int tmp_cell_1 = face_can_idx[face_can_ext[i_face]+1];
                  if (tmp_cell_0 == tmp_cell_1){
                      (*face_entity_data)(0,i_face) = (*cell_entity_data)(0,tmp_cell_0);
                      (*face_entity_data)(1,i_face) = (*cell_entity_data)(0,tmp_cell_1);
                      (*face_mat_data)(0,i_face) = (*cell_mat_data)(0,tmp_cell_0);
                      (*face_mat_data)(1,i_face) = (*cell_mat_data)(0,tmp_cell_1);
                  }else{
                      /** get material for cell **/
                      tbox::Pointer<BaseMaterial<NDIM> > mat_cell_0 =
                              mat_manager->getMaterial(material_list[(*cell_mat_data)(0,tmp_cell_0)].mat_name);
                      tbox::Pointer<BaseMaterial<NDIM> > mat_cell_1 =
                              mat_manager->getMaterial(material_list[(*cell_mat_data)(0,tmp_cell_1)].mat_name);

                      if (mat_cell_0->effectiveAffinity(300.0) > mat_cell_1->effectiveAffinity(300.0)){
                          (*face_entity_data)(0,i_face) = (*cell_entity_data)(0,tmp_cell_0);
                          (*face_entity_data)(1,i_face) = (*cell_entity_data)(0,tmp_cell_1);

                          (*face_mat_data)(0,i_face) = (*cell_mat_data)(0,tmp_cell_0);
                          (*face_mat_data)(1,i_face) = (*cell_mat_data)(0,tmp_cell_1);
                      }else {
                          (*face_entity_data)(0,i_face) = (*cell_entity_data)(0,tmp_cell_1);
                          (*face_entity_data)(1,i_face) = (*cell_entity_data)(0,tmp_cell_0);

                          (*face_mat_data)(0,i_face) = (*cell_mat_data)(0,tmp_cell_1);
                          (*face_mat_data)(1,i_face) = (*cell_mat_data)(0,tmp_cell_0);
                      }
                  }
              }
          }

          /** heterojunction node-entity: larger affinity written first **/
          tbox::Array<int> face_node_ext, face_node_idx;
          patch_top->getFaceAdjacencyNodes(face_node_ext, face_node_idx);
          for(int i_bc=0; i_bc<num_bc; ++i_bc){
              if (EBoundassignment[i_bc].bc_type == "heterojunction"){
                  int num_face = EBoundassignment[i_bc].bc_face.getSize();
                  for (int i_f = 0; i_f < num_face; ++i_f){
                      if(patch_geo->hasEntitySet(EBoundassignment[i_bc].bc_face[i_f], hier::EntityUtilities::FACE)){
                          const tbox::Array<int> & face_id = patch_geo->getEntityIndicesInSet(EBoundassignment[i_bc].bc_face[i_f],
                                                             hier::EntityUtilities::FACE,
                                                             num_faces);
                          int size_f = face_id.size();
                          for (int i_e = 0; i_e < size_f; ++i_e){
                              int n_fvertex = face_node_ext[face_id[i_e] + 1] - face_node_ext[face_id[i_e]];
                              for (int i1 = 0, j = face_node_ext[face_id[i_e]]; i1 < n_fvertex; ++i1, ++j) {
                                  int node_n = face_node_idx[j];
                                  (*hnode_entity_data)(0,node_n) = (*face_entity_data)(0,face_id[i_e]);
                                  (*hnode_entity_data)(1,node_n) = (*face_entity_data)(1,face_id[i_e]);

                                  (*hnode_entity_data)(0,node_n) = entity_mat_list[(*face_entity_data)(0,face_id[i_e])-1];
                                  (*hnode_entity_data)(1,node_n) = entity_mat_list[(*face_entity_data)(1,face_id[i_e])-1];
                              }
                          }
                      }
                  }
              }
          }

          /** cell - heterojunction node relationship **/
          for (int i_e = 0; i_e < num_cells; ++i_e){
              int n_vertex = can_extent[i_e+1] - can_extent[i_e];
              for (int i1 = 0, j = can_extent[i_e]; i1 < n_vertex; ++i1, ++j) {
                  if ((*heter_node_data)(0,can_indices[j]) == 1){
                      if((*cell_entity_data)(0,i_e) == (*hnode_entity_data)(0, can_indices[j])){
                          (*cell_hnode_data)(0,i_e) = 1;  /** GaN **/
                      }else{
                          (*cell_hnode_data)(0,i_e) = 2;  /** AlGaN **/
                      }
                  }
              }
          }

          for(int i = 0; i < num_nodes; ++i){
              tbox::Pointer<Variables_Functions> Var_func = new Variables_Functions();
              (*poisson_plot_data)(0, i) = 0.0;
              (*carn_plot_data)(0, i) = 0.0;
              (*logn_plot_data)(0, i) = 0.0;
              (*carh_plot_data)(0, i) = 0.0;
              (*logh_plot_data)(0, i) = 0.0;
              (*fermin_plot_data)(0, i) = 0.0;
              (*fermih_plot_data)(0, i) = 0.0;
              (*node_plot_data)(0,i) = 0.0;

              (*carn_sol_1_data)(0,i) = 0.0;
              (*carn_sol_1_data)(1,i) = 0.0;
              (*fn_sol_1_data)(0,i) = 0.0;
              (*fn_sol_1_data)(1,i) = 0.0;
              (*fn_sol_0_data)(0,i) = 0.0;
              (*fn_sol_0_data)(1,i) = 0.0;

              (*carh_sol_1_data)(0,i) = 0.0;
              (*carh_sol_1_data)(1,i) = 0.0;
              (*fh_sol_1_data)(0,i) = 0.0;
              (*fh_sol_1_data)(1,i) = 0.0;
              (*fh_sol_0_data)(0,i) = 0.0;
              (*fh_sol_0_data)(1,i) = 0.0;

              (*phi_sol_0_data)(0,i) = 0.0;
              (*phi_sol_1_data)(0,i) = 0.0;
          }
          for(int l=0;l<num_cells;++l){
              for(int i_d = 0; i_d < 3; ++i_d){
                  (*e_field_data)(i_d,l) = 0.0;
                  (*e_cell_data)(i_d,l) = 0.0;
                  (*jn_field_data)(i_d,l) = 0.0;
                  (*jp_field_data)(i_d,l) = 0.0;
                  (*j_field_data)(i_d,l) = 0.0;
              }
              for(int ii = 0; ii < cell_edge; ++ii){
                  (*j_edge_data)(ii,l) = 0.0;
              }
          }
      }
    }
}

/** ***********************************************************************
 * 完成单个网格片上的数值计算（支持数值构件）.
 ********************************************************************** **/
void NonlinearPoisson::computeOnPatch(hier::Patch<NDIM>& patch,
                                      const double  time,
                                      const double  dt,
                                      const bool    initial_time,
                                      const string& component_name)
{
    /** current conservation equation **/
    if(component_name=="CCMAT"){
        buildCurConserveMatOnPatch(patch, time, dt, component_name);
    }
    else if(component_name=="CCPHYBC"){
        setCurConservePhysicalBC(patch, time, dt, component_name);
    }
    else if(component_name=="CCAPPBC"){
        applyCurConserveBCOnPatch(patch, time, dt, component_name);
    }

    /** Poisson Equations: initialization **/
    else if (component_name=="PMAT"){
        buildPoissonMatOnPatch(patch, time, dt, component_name);
    }
    else if (component_name=="PFX"){
        buildPoissonFXOnPatch(patch, time, dt, component_name);
    }
    else if (component_name=="PPHYBC"){
        setPoissonPhysicalBC(patch, time, dt, component_name);
    }
    else if (component_name=="PAppBC"){
        applyPoissonBCOnPatch(patch, time, dt, component_name);
    }

    /** Poisson Equation: Newton-Raphson **/
    else if (component_name=="ADTEST"){
        adTest(patch, time, dt, component_name);
    }
    else if (component_name=="NewtonPMAT"){
        buildNewtonPMatOnPatch(patch, time, dt, component_name);
    }
    else if (component_name=="NewtonPFX"){
        buildNewtonPFXOnPatch(patch, time, dt, component_name);
    }
    else if (component_name=="NewtonPPHYBC"){
        setNewtonPPhysicalBC(patch, time, dt, component_name);
    }
    else if (component_name=="NewtonPAppBC"){
        applyNewtonPBCOnPatch(patch, time, dt, component_name);
    }
    else if (component_name=="PPOST"){
        transferPoissonSoln(patch, time, dt, component_name);
    }

    /** Carrier Continuity Equations: CV-FEM-SG **/
    else if (component_name == "CVMAT"){
        buildCarMatOnPatch(patch, time, dt, component_name);
    }

    /** Electron Continuity Equations: CV-FEM-SG **/
    else if (component_name=="NFX"){
        buildCarnFXOnPatch(patch, time, dt, component_name);
    }
    else if (component_name=="NPHYBC"){
        setCarnPhysicalBC(patch, time, dt, component_name);
    }
    else if (component_name=="NAppBC"){
        applyCarnBCOnPatch(patch, time, dt, component_name);
    }
    else if (component_name=="NPOST"){
        transferCarnSoln(patch, time, dt, component_name);
    }

    /** Hole Continuity Equations: CV-FEM-SG **/
    else if (component_name=="HFX"){
        buildCarhFXOnPatch(patch, time, dt, component_name);
    }
    else if (component_name=="HPHYBC"){
        setCarhPhysicalBC(patch, time, dt, component_name);
    }
    else if (component_name=="HAppBC"){
        applyCarhBCOnPatch(patch, time, dt, component_name);
    }
    else if (component_name=="HPOST"){
        transferCarhSoln(patch, time, dt, component_name);
    }

    /** Qusi fermi level **/
    else if (component_name=="FLEVEL"){
        calculateQuasiFermilevel(patch, time, dt, component_name);
    }
    else if (component_name=="NPCONS"){
        applyConstraint2NP(patch, time, dt, component_name);
    }
    else if (component_name=="FINTRIN"){
        initialQuasiFermilevel(patch, time, dt, component_name);
    }

    else if (component_name=="CHECKORDER"){
        checkOrder(patch, time, dt, component_name);
    }

    /** initialization **/
    else if (component_name=="SOLINIT"){
        solVecInitialization(patch, time, dt, component_name);
    }
    else if (component_name=="EQUILIBRIUM"){
        neutral_status(patch, time, dt, component_name);
    }

    /** solution transport: current to stetch **/
    else if (component_name=="PSOLTRANS"){
        transportPoissonSol(patch, time, dt, component_name);
    }
    else if (component_name=="GummelPSOLTRANS"){
        transportGummelPoissonSol(patch, time, dt, component_name);
    }
    else if (component_name=="SEMISOLTRANS"){
        transportSemiSol(patch, time, dt, component_name);
    }
    else if (component_name=="GummelSEMISOLTRANS"){
        transportGummelSemiSol(patch, time, dt, component_name);
    }
    else if ( component_name=="SOLPROCESS") {
        stepSolPostProcess(patch, time, dt, component_name);
    }
    else if (component_name=="calculate_E"){
        calculateEx(patch, time, dt, component_name);
    }
    else if (component_name=="calculate_J"){
        calculateElemJ(patch, time, dt, component_name);
    }
    else if (component_name=="get_IV"){
        getIV(patch, time, dt, component_name);
    }

    /** variable scaling and anti-scaling **/
    else if (component_name=="ANTISCALOR"){
        solutionAntiScal(patch, time, dt, component_name);
    }
    else if (component_name=="SCALOR"){
        solutionScal(patch, time, dt, component_name);
    }
    else if (component_name=="ANTICARRIER"){
        carrierAntiScal(patch, time, dt, component_name);
    }

    /** *********888888888888888******* **/
    else {
        TBOX_ERROR(" NonlinearPoisson :: component name is error! " );
    }
}

/**
 * @brief solVecInitialization: initialize solution vectors: data from plot vectors
 * @param patch
 * @param time
 * @param dt
 * @param component_name
 */
void NonlinearPoisson::solVecInitialization(hier::Patch<NDIM> & patch,
                                            const double  time,
                                            const double  dt,
                                            const string& component_name){
    /** get the solution vectors **/
    //cout<<"get the solution vectors"<<endl;

    tbox::Pointer<pdat::NodeData<NDIM, double> >
          potential_sol_0_data = patch.getPatchData(d_phi_sol_0_id);
    tbox::Pointer<pdat::NodeData<NDIM, double> >
          potential_sol_1_data = patch.getPatchData(d_phi_sol_1_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
          potential_sol_cur_data = patch.getPatchData(d_poisson_sol_cur_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
          potential_stetch_data = patch.getPatchData(d_poisson_stetch_id);

    /** tbox::Pointer<pdat::NodeData<NDIM, double> >
          carn_sol_1_data = patch.getPatchData(d_carn_sol_1_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
          carn_sol_cur_data = patch.getPatchData(d_carn_sol_cur_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
          carn_stetch_data = patch.getPatchData(d_carn_stetch_id); **/

    tbox::Pointer<pdat::NodeData<NDIM, double> >
          fermin_sol_1_data = patch.getPatchData(d_fn_sol_1_id);
    tbox::Pointer<pdat::NodeData<NDIM, double> >
          fermin_sol_0_data = patch.getPatchData(d_fn_sol_0_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
          fermin_sol_cur_data = patch.getPatchData(d_fermin_sol_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
          fermin_stetch_data = patch.getPatchData(d_fermin_stetch_id);

    tbox::Pointer<pdat::VectorData<NDIM, double> >
          carn_stetch_data = patch.getPatchData(d_carn_stetch_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
          carh_stetch_data = patch.getPatchData(d_carh_stetch_id);

    tbox::Pointer<pdat::NodeData<NDIM, double> >
          carn_sol_1_data = patch.getPatchData(d_carn_sol_1_id);
    tbox::Pointer<pdat::NodeData<NDIM, double> >
          carh_sol_1_data = patch.getPatchData(d_carh_sol_1_id);

    /** tbox::Pointer<pdat::NodeData<NDIM, double> >
          carh_sol_1_data = patch.getPatchData(d_carh_sol_1_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
          carh_sol_cur_data = patch.getPatchData(d_carh_sol_cur_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
          carh_stetch_data = patch.getPatchData(d_carh_stetch_id); **/

    tbox::Pointer<pdat::NodeData<NDIM, double> >
          fermih_sol_1_data = patch.getPatchData(d_fh_sol_1_id);
    tbox::Pointer<pdat::NodeData<NDIM, double> >
          fermih_sol_0_data = patch.getPatchData(d_fh_sol_0_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
          fermih_sol_cur_data = patch.getPatchData(d_fermih_sol_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
          fermih_stetch_data = patch.getPatchData(d_fermih_stetch_id);

    tbox::Pointer<pdat::NodeData<NDIM, int> >
          heter_node = patch.getPatchData(d_heter_node_id);
    tbox::Pointer<pdat::NodeData<NDIM, int> >
          semiinsbc_node = patch.getPatchData(d_semiinsbc_node_id);
    tbox::Pointer<pdat::NodeData<NDIM, int> >
          semi_node = patch.getPatchData(d_semi_node_id);
    tbox::Pointer<pdat::NodeData<NDIM, int> >
          siclass_node = patch.getPatchData(d_siclass_node_id);

    tbox::Pointer<pdat::NodeData<NDIM, double> >
            dopant_data = patch.getPatchData(d_dopant_id);

    int * dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);
    int * dof_map_semi = d_dof_semi_info->getDOFMapping(patch, hier::EntityUtilities::NODE);

    int num_n = patch.getNumberOfNodes(1);

    /** setup the source target **/
    tbox::Pointer<bound_condition> bound_source = new bound_condition();

    int num_bound = EBoundassignment.getSize();
    for (int i_b = 0; i_b < num_bound; ++i_b){
            double phi = bound_source->source_assign(EBoundassignment[i_b].sc_type,
                                                       time, dt,
                                                       EBoundassignment[i_b].sc_coe);

            tbox::pout<<"time:"<<time<<" electrode:"<<i_b<<" type:"<<EBoundassignment[i_b].sc_type<<" voltage:"<<phi<<endl;
    }

    if (time == 0){
        if(enable_unit_scalling){
            //cout<<"scalling unit"<<endl;
            //max_dop = Get_MaxDop(dopant_file_name)*1e-6;  //*1e-6是为了将单位从m-3转化成cm-3
            //double max_dop = max_doping;
            double max_dop = parameter.getmaxdop()*1e-6;  //*1e-6是为了将单位从m-3转化成cm-3
            double length = pow(max_dop,1.0/3.0);  //pow(max_dop,1/3)会出现类型强制转化问题，变成计算pow(max_dop,0)，ni=6621424696833250e-6
            PhysicalUnit::set_unit(length);

            /** material list **/
            tbox::Pointer<MaterialManager<NDIM> >
                        mat_manager = MaterialManager<NDIM>::getManager();
            /** get material for cell **/
            tbox::Pointer<BaseMaterial<NDIM> > mat_cell =
                    mat_manager->getMaterial("matSilicon");
            T = 300 * K;
        }
        if(save_bas){
            /** Get Local PatchTopology **/
            tbox::Pointer< hier::PatchTopology<NDIM> > patch_top = patch.getPatchTopology();

            tbox::Pointer< hier::PatchGeometry<NDIM> > patch_geo = patch.getPatchGeometry();
            /** get Patch node coord **/
            tbox::Pointer<pdat::NodeData<NDIM, double> >
                    node_coord = patch_geo->getNodeCoordinates();
            //bas vector
            /// (形函数, 积分器, 单元)管理器
            tbox::Pointer<ElementManager<NDIM> > ele_manager =
                ElementManager<NDIM>::getManager();
            tbox::Pointer<BaseElement<NDIM> > ele =
                    ele_manager->getElement(d_element_type);

            tbox::Pointer<pdat::CellData<NDIM, double> > bas_integ = patch.getPatchData(d_bas_integ_id);
            tbox::Pointer<pdat::CellData<NDIM, double> > bas_rhs = patch.getPatchData(d_bas_rhs_id);

            int num_cells = patch.getNumberOfCells(1);
            tbox::Array<int> can_extent, can_indices;
            patch_top->getCellAdjacencyNodes(can_extent, can_indices);
            //extent 数组首尾地址。 can_indices是节点编号
            for (int i = 0; i < num_cells; ++i) {
                int cell=i;
                int n_vertex = can_extent[i + 1] - can_extent[i];
                tbox::Array<hier::DoubleVector<NDIM> > vertex(n_vertex);
                /// 下面的循环做一件事情：1.取出结点坐标。
                for (int i1 = 0, j = can_extent[i]; i1 < n_vertex; ++i1, ++j) {
                    for (int k = 0; k < NDIM; ++k) {
                        vertex[i1][k] = (*node_coord)(k, can_indices[j])*m;
                    }
                }
                /// 声明ele_bas_value、ele_bas_grad和ele_bas_vec，用来存储单元基函数。
                tbox::Pointer<tbox::Vector<double> > ele_integ = new tbox::Vector<double>();
                tbox::Pointer<tbox::Vector<double> > ele_rhs = new tbox::Vector<double>();
                ele_integ->resize(cell_edge*cell_point);
                ele_rhs->resize(cell_point);
                /// 计算单元基函数，并组装。
                ele->buildBas(vertex,ele_integ,ele_rhs);
                //将计算得到的值赋值给BAS_VALUE变量
                for(int j1=0;j1<cell_edge*cell_point;j1++){
                    //cout<<"cell:"<<cell<<endl;
                    (*bas_integ)(j1,cell)=(*ele_integ)[j1];
                }
                //将计算得到的值赋值给BAS_GRAD变量
                for(int j1=0;j1<cell_point;j1++){
                    //cout<<"cell:"<<cell<<endl;
                    (*bas_rhs)(j1,cell)=(*ele_rhs)[j1];
                }
            }
            //cout<<"save bas already."<<endl;
        }

        for (int i_n = 0; i_n < num_n; ++i_n){
            dopant_data->getPointer()[i_n] = dopant_data->getPointer()[i_n] * pow(m,-3);
            //int mapping = dof_map[i_n];
            /** if((*semi_node)(0,i_n) == 1){

                int semi_mapping = dof_map_semi[i_n];

                (*carn_stetch_data)(semi_mapping) = (*carn_sol_cur_data)(semi_mapping);
                (*carh_stetch_data)(semi_mapping) = (*carh_sol_cur_data)(semi_mapping);
            } **/

            if ((*heter_node)(0,i_n) == 1){
                /** int semi_mapping = dof_map_semi[i_n];
                if((*semi_node)(0,i_n) == 1){
                    (*carn_stetch_data)(semi_mapping+1) = (*carn_sol_cur_data)(semi_mapping+1);
                    (*carh_stetch_data)(semi_mapping+1) = (*carh_sol_cur_data)(semi_mapping+1);
                } **/
            }
        }

        /** initialize potential **/
        neutral_status(patch, time, dt, component_name);
        //MPI_Barrier(MPI_COMM_WORLD);
        //cout<<"initialize potential"<<endl;

        for (int i_n = 0; i_n < num_n; ++i_n) {
            if((*siclass_node)(0,i_n)==0) continue;

            int mapping_semi = dof_map_semi[i_n];
            int mapping = dof_map[i_n];
            (*potential_stetch_data)(mapping_semi) = (*potential_sol_cur_data)(mapping_semi);
            (*potential_stetch_data)(mapping_semi+1) = (*potential_sol_cur_data)(mapping_semi+1);
            (*potential_stetch_data)(mapping_semi+2) = (*potential_sol_cur_data)(mapping_semi+2);

            (*carn_stetch_data)(mapping) = (*potential_sol_cur_data)(mapping_semi+1);
            (*carh_stetch_data)(mapping) = (*potential_sol_cur_data)(mapping_semi+2);

            if((*carn_stetch_data)(mapping)<=-1){
                cout<<"carn_stetch_data:"<<(*carn_stetch_data)(mapping)<<endl;
            }

            if((*semiinsbc_node)(0,i_n)==1){
                (*potential_stetch_data)(mapping_semi+3) = (*potential_sol_cur_data)(mapping_semi+3);
                (*potential_stetch_data)(mapping_semi+4) = (*potential_sol_cur_data)(mapping_semi+4);
                (*potential_stetch_data)(mapping_semi+5) = (*potential_sol_cur_data)(mapping_semi+5);

                (*carn_stetch_data)(mapping+1) = (*potential_sol_cur_data)(mapping_semi+4);
                (*carh_stetch_data)(mapping+1) = (*potential_sol_cur_data)(mapping_semi+5);
            }
            //cout<<i_n<<" "<<(*fermih_stetch_data)(mapping)<<endl;
        }

    }else if(time == dt){

        int num_bound = EBoundassignment.getSize();
        double t_p  = time - dt;
        double delt_new;

        for (int i_b = 0; i_b < num_bound; ++i_b){
            double phi_s = bound_source->source_assign(EBoundassignment[i_b].sc_type,
                                                       time, dt,
                                                       EBoundassignment[i_b].sc_coe);

            double phi_p = bound_source->source_assign(EBoundassignment[i_b].sc_type,
                                                       t_p, dt,
                                                       EBoundassignment[i_b].sc_coe);

            delt_new = phi_s - phi_p;

            if (abs(delt_new)>0) break;
        }

        for (int i_n = 0; i_n < num_n; ++i_n){
            if((*siclass_node)(0,i_n)==0) continue;

            int mapping = dof_map[i_n];
            int mapping_semi = dof_map_semi[i_n];

            (*potential_stetch_data)(mapping_semi) = (*potential_sol_1_data)(0,i_n);
            (*potential_sol_cur_data)(mapping_semi) = (*potential_sol_1_data)(0,i_n);
            (*potential_sol_cur_data)(mapping_semi+1) = (*carn_sol_1_data)(0,i_n);
            (*potential_sol_cur_data)(mapping_semi+2) = (*carh_sol_1_data)(0,i_n);

            (*fermin_stetch_data)(mapping) = (*fermin_sol_1_data)(0,i_n);
            (*fermin_sol_cur_data)(mapping) = (*fermin_sol_1_data)(0,i_n);

            (*fermih_stetch_data)(mapping) = (*fermih_sol_1_data)(0,i_n);
            (*fermih_sol_cur_data)(mapping) = (*fermih_sol_1_data)(0,i_n);

            (*carn_stetch_data)(mapping) = (*carn_sol_1_data)(0,i_n);
            (*carh_stetch_data)(mapping) = (*carh_sol_1_data)(0,i_n);

            if((*semiinsbc_node)(0,i_n)==1){
                mapping = dof_map[i_n]+1;
                mapping_semi = dof_map_semi[i_n]+3;
                (*potential_stetch_data)(mapping_semi) = (*potential_sol_1_data)(1,i_n);
                (*potential_sol_cur_data)(mapping_semi) = (*potential_sol_1_data)(1,i_n);
                (*potential_sol_cur_data)(mapping_semi+1) = (*carn_sol_1_data)(1,i_n);
                (*potential_sol_cur_data)(mapping_semi+2) = (*carh_sol_1_data)(1,i_n);

                (*fermin_stetch_data)(mapping) = (*fermin_sol_1_data)(1,i_n);
                (*fermin_sol_cur_data)(mapping) = (*fermin_sol_1_data)(1,i_n);

                (*fermih_stetch_data)(mapping) = (*fermih_sol_1_data)(1,i_n);
                (*fermih_sol_cur_data)(mapping) = (*fermih_sol_1_data)(1,i_n);

                (*carn_stetch_data)(mapping) = (*carn_sol_1_data)(1,i_n);
                (*carh_stetch_data)(mapping) = (*carh_sol_1_data)(1,i_n);
            }
        }

        /** transfer solution for initialization of next step **/
        for (int i_n = 0; i_n < num_n; ++ i_n){
            (*potential_sol_0_data)(0,i_n) = (*potential_sol_1_data)(0,i_n);
            (*fermin_sol_0_data)(0,i_n) = (*fermin_sol_1_data)(0,i_n);
            (*fermih_sol_0_data)(0,i_n) = (*fermih_sol_1_data)(0,i_n);
        }

        /** update potential boundaries **/
        potential_bound(patch, time, dt, component_name);

    }else{
        /** time **/
        double t_p  = time - dt;
        double t_pp = time - 2*dt;

        int num_bound = EBoundassignment.size();
        double delt_pre, delt_new;

        for (int i_b = 0; i_b < num_bound; ++i_b){
            double phi_s = bound_source->source_assign(EBoundassignment[i_b].sc_type,
                                                       time, dt,
                                                       EBoundassignment[i_b].sc_coe);
            double phi_p = bound_source->source_assign(EBoundassignment[i_b].sc_type,
                                                       t_p, dt,
                                                       EBoundassignment[i_b].sc_coe);
            double phi_pp = bound_source->source_assign(EBoundassignment[i_b].sc_type,
                                                        t_pp, dt,
                                                        EBoundassignment[i_b].sc_coe);

            delt_new = phi_s - phi_p;
            delt_pre = phi_p - phi_pp;

            if(abs(delt_new)>0) break;
        }

        if(delt_pre == 0){
            for (int i_n = 0; i_n < num_n; ++i_n){
                if((*siclass_node)(0,i_n)==0) continue;

                int mapping = dof_map[i_n];
                int mapping_semi = dof_map_semi[i_n];

                (*potential_stetch_data)(mapping_semi) = (*potential_sol_1_data)(0,i_n);
                (*potential_sol_cur_data)(mapping_semi) = (*potential_sol_1_data)(0,i_n);
                (*potential_sol_cur_data)(mapping_semi+1) = (*carn_sol_1_data)(0,i_n);
                (*potential_sol_cur_data)(mapping_semi+2) = (*carh_sol_1_data)(0,i_n);

                (*fermin_stetch_data)(mapping) = (*fermin_sol_1_data)(0,i_n);
                (*fermin_sol_cur_data)(mapping) = (*fermin_sol_1_data)(0,i_n);

                (*fermih_stetch_data)(mapping) = (*fermih_sol_1_data)(0,i_n);
                (*fermih_sol_cur_data)(mapping) = (*fermih_sol_1_data)(0,i_n);

                (*carn_stetch_data)(mapping) = (*carn_sol_1_data)(0,i_n);
                (*carh_stetch_data)(mapping) = (*carh_sol_1_data)(0,i_n);

                if((*semiinsbc_node)(0,i_n)==1){
                    mapping = dof_map[i_n]+1;
                    mapping_semi = dof_map_semi[i_n]+3;
                    (*potential_stetch_data)(mapping_semi) = (*potential_sol_1_data)(1,i_n);
                    (*potential_sol_cur_data)(mapping_semi) = (*potential_sol_1_data)(1,i_n);
                    (*potential_sol_cur_data)(mapping_semi+1) = (*carn_sol_1_data)(1,i_n);
                    (*potential_sol_cur_data)(mapping_semi+2) = (*carh_sol_1_data)(1,i_n);

                    (*fermin_stetch_data)(mapping) = (*fermin_sol_1_data)(1,i_n);
                    (*fermin_sol_cur_data)(mapping) = (*fermin_sol_1_data)(1,i_n);

                    (*fermih_stetch_data)(mapping) = (*fermih_sol_1_data)(1,i_n);
                    (*fermih_sol_cur_data)(mapping) = (*fermih_sol_1_data)(1,i_n);

                    (*carn_stetch_data)(mapping) = (*carn_sol_1_data)(1,i_n);
                    (*carh_stetch_data)(mapping) = (*carh_sol_1_data)(1,i_n);
                }
            }

            /** transfer solution for initialization of next step **/
            for (int i_n = 0; i_n < num_n; ++ i_n){
                if((*siclass_node)(0,i_n)==0) continue;
                (*potential_sol_0_data)(0,i_n) = (*potential_sol_1_data)(0,i_n);
                (*fermin_sol_0_data)(0,i_n) = (*fermin_sol_1_data)(0,i_n);
                (*fermih_sol_0_data)(0,i_n) = (*fermih_sol_1_data)(0,i_n);
                if((*heter_node)(0,i_n) == 1){
                    (*fermin_sol_0_data)(1,i_n) = (*fermin_sol_1_data)(1,i_n);
                    (*fermih_sol_0_data)(1,i_n) = (*fermih_sol_1_data)(1,i_n);
                }
            }

        }else{
            for (int i_n = 0; i_n < num_n; ++i_n){
                if((*siclass_node)(0,i_n)==0) continue;

                double coe = delt_new/delt_pre;
                coe = 0;
                int mapping = dof_map[i_n];
                int mapping_semi = dof_map_semi[i_n];

                (*potential_stetch_data)(mapping_semi) = (*potential_sol_1_data)(0,i_n) + coe*((*potential_sol_1_data)(0,i_n) - (*potential_sol_0_data)(0,i_n));
                (*potential_sol_cur_data)(mapping_semi) = (*potential_sol_1_data)(0,i_n) + coe*((*potential_sol_1_data)(0,i_n) - (*potential_sol_0_data)(0,i_n));

                (*potential_sol_cur_data)(mapping_semi+1) = (*carn_sol_1_data)(0,i_n);
                (*potential_sol_cur_data)(mapping_semi+2) = (*carh_sol_1_data)(0,i_n);

                (*fermin_stetch_data)(mapping) = (*fermin_sol_1_data)(0,i_n) + coe*((*fermin_sol_1_data)(0,i_n) - (*fermin_sol_0_data)(0,i_n));
                (*fermin_sol_cur_data)(mapping) = (*fermin_sol_1_data)(0,i_n) + coe*((*fermin_sol_1_data)(0,i_n) - (*fermin_sol_0_data)(0,i_n));

                (*fermih_stetch_data)(mapping) = (*fermih_sol_1_data)(0,i_n) + coe*((*fermih_sol_1_data)(0,i_n) - (*fermih_sol_0_data)(0,i_n));
                (*fermih_sol_cur_data)(mapping) = (*fermih_sol_1_data)(0,i_n) + coe*((*fermih_sol_1_data)(0,i_n) - (*fermih_sol_0_data)(0,i_n));

                (*carn_stetch_data)(mapping) = (*carn_sol_1_data)(0,i_n);
                (*carh_stetch_data)(mapping) = (*carh_sol_1_data)(0,i_n);

                if((*semiinsbc_node)(0,i_n)==1){
                    mapping = dof_map[i_n]+1;
                    mapping_semi = dof_map_semi[i_n]+3;
                    (*potential_stetch_data)(mapping_semi) = (*potential_sol_1_data)(1,i_n);
                    (*potential_sol_cur_data)(mapping_semi) = (*potential_sol_1_data)(1,i_n);
                    (*potential_sol_cur_data)(mapping_semi+1) = (*carn_sol_1_data)(1,i_n);
                    (*potential_sol_cur_data)(mapping_semi+2) = (*carh_sol_1_data)(1,i_n);

                    (*fermin_stetch_data)(mapping) = (*fermin_sol_1_data)(1,i_n);
                    (*fermin_sol_cur_data)(mapping) = (*fermin_sol_1_data)(1,i_n);

                    (*fermih_stetch_data)(mapping) = (*fermih_sol_1_data)(1,i_n);
                    (*fermih_sol_cur_data)(mapping) = (*fermih_sol_1_data)(1,i_n);

                    (*carn_stetch_data)(mapping) = (*carn_sol_1_data)(1,i_n);
                    (*carh_stetch_data)(mapping) = (*carh_sol_1_data)(1,i_n);
                }
            }

            /** transfer solution for initialization of next step **/
            for (int i_n = 0; i_n < num_n; ++ i_n){
                (*potential_sol_0_data)(0,i_n) = (*potential_sol_1_data)(0,i_n);
                (*fermin_sol_0_data)(0,i_n) = (*fermin_sol_1_data)(0,i_n);
                (*fermih_sol_0_data)(0,i_n) = (*fermih_sol_1_data)(0,i_n);
                if((*heter_node)(0,i_n) == 1){
                    (*fermin_sol_0_data)(1,i_n) = (*fermin_sol_1_data)(1,i_n);
                    (*fermih_sol_0_data)(1,i_n) = (*fermih_sol_1_data)(1,i_n);
                }
            }
        }

        /** update potential boundaries **/
        potential_bound(patch, time, dt, component_name);
    }
}

/**
 * @brief NonlinearPoisson::neutral_status: neutral status for semiconductor domains
 * @param patch
 * @param time
 * @param dt
 * @param component_name
 */
void NonlinearPoisson::neutral_status(hier::Patch<NDIM> &patch,
                                      const double  time,
                                      const double dt,
                                      const string& component_name){
    /** gPatchGeometry **/
    tbox::Pointer< hier::PatchGeometry<NDIM> >
            patch_geo = patch.getPatchGeometry();
    /** get PatchTopology **/
    tbox::Pointer< hier::PatchTopology<NDIM> >
            patch_top = patch.getPatchTopology();

    tbox::Pointer<pdat::VectorData<NDIM, double> >
      V_soln = patch.getPatchData(d_poisson_sol_cur_id);

    tbox::Pointer<pdat::VectorData<NDIM, double> >
      Fn_soln = patch.getPatchData(d_fermin_sol_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
      Fp_soln = patch.getPatchData(d_fermih_sol_id);

    tbox::Pointer<pdat::NodeData<NDIM, double> >
      dop = patch.getPatchData(d_dopant_id);

    tbox::Pointer<pdat::NodeData<NDIM, int> >
      node_mat = patch.getPatchData(d_node_mat_id);

    tbox::Pointer<pdat::CellData<NDIM, int> >
      cell_mat = patch.getPatchData(d_cell_mat_id);

    tbox::Pointer<pdat::CellData<NDIM, int> >
      cell_hnode = patch.getPatchData(d_cell_hnode_id);

    tbox::Pointer<pdat::NodeData<NDIM, int> >
      heter_node = patch.getPatchData(d_heter_node_id);
    tbox::Pointer<pdat::NodeData<NDIM, int> >
      semiinsbc_node = patch.getPatchData(d_semiinsbc_node_id);
    tbox::Pointer<pdat::NodeData<NDIM, int> >
      siclass_node = patch.getPatchData(d_siclass_node_id);

    tbox::Pointer<pdat::NodeData<NDIM, double> > sol;
    if(get_InitialSol||get_PoissonSolcur)    sol = patch.getPatchData(d_sol_id);

    /** 获取本地单元(边, 面)周围结点的索引关系 **/
    tbox::Array<int> can_extent, can_indices;
    patch_top->getCellAdjacencyNodes(can_extent, can_indices);

    /** material list **/
    tbox::Pointer<MaterialManager<NDIM> >
            mat_manager = MaterialManager<NDIM>::getManager();

    tbox::Pointer<Variables_Functions>
            Var_func = new Variables_Functions();

    int * dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);
    int * dof_map_semi = d_dof_semi_info->getDOFMapping(patch, hier::EntityUtilities::NODE);

    int num_n = patch.getNumberOfNodes(1);
    int num_e = patch.getNumberOfCells(1);


    for (int i_e = 0; i_e < num_e; ++i_e){

        if(!((material_list[(*cell_mat)(0, i_e)].mat_type == "semiconductor")||(material_list[(*cell_mat)(0, i_e)].mat_type == "dielectric"))) continue;

        /** get node information **/
        int n_vertex = can_extent[i_e+1] - can_extent[i_e];
        tbox::Array<int> e_node(n_vertex);
        tbox::Array<int> mapping(n_vertex);
        tbox::Array<int> semi_mapping(n_vertex);

        for (int i1 = 0, j = can_extent[i_e]; i1 < n_vertex; ++i1, ++j) {
            mapping[i1] = dof_map[can_indices[j]];
            semi_mapping[i1] = dof_map_semi[can_indices[j]];
            e_node[i1] = can_indices[j];
        }

        /** setup cell material **/
        tbox::Pointer<BaseMaterial<NDIM> > mat_cell =
                mat_manager->getMaterial(material_list[(*cell_mat)(0, i_e)].mat_name);

        /** calculate initial potential at nodes **/
        double tmp_ni = mat_cell->intrinsicDensityNi(T_amb)*pow(m,-3);
        double tmp_nc = mat_cell->effectiveDensityNc(T_amb)*pow(m,-3);
        double tmp_nv = mat_cell->effectiveDensityNv(T_amb)*pow(m,-3);
        double tmp_affin = mat_cell->effectiveAffinity(T_amb);
        double tmp_eg = mat_cell->bandGap(T_amb);

        if(material_list[(*cell_mat)(0, i_e)].mat_type == "semiconductor"){
            for(int l=0; l<n_vertex; ++l){
                double tmp_dop = (*dop)(0,e_node[l]);
                (*V_soln)(semi_mapping[l]) = Var_func->equilibrium_potential(tmp_dop,
                               tmp_ni, tmp_nc, tmp_nv, tmp_affin, tmp_eg, T, kb, e);
                //(*V_soln)(semi_mapping[l]+1) = 0;
                //(*V_soln)(semi_mapping[l]+2) = 0;
                #if 1
                (*V_soln)(semi_mapping[l]+1) = Var_func->dummyElectron((*V_soln)(semi_mapping[l]),
                               (*Fn_soln)(mapping[l]), false, tmp_affin, T, kb, e, tmp_nc);
                (*V_soln)(semi_mapping[l]+2) = Var_func->dummyHole((*V_soln)(semi_mapping[l]),
                               (*Fp_soln)(mapping[l]), tmp_eg, false, tmp_affin, T, kb, e, tmp_nv);
                #endif
                //cout<<"V_soln:"<<(*V_soln)(semi_mapping[l])<<" "<<(*V_soln)(semi_mapping[l]+1)<<" "<<(*V_soln)(semi_mapping[l]+2)<<endl;
                //cout<<"Fn_soln:"<<(*Fn_soln)(mapping[l])<<endl;

                //cout<<"equilibrium_potential"<<Var_func->equilibrium_potential(tmp_dop,tmp_ni, tmp_nc, tmp_nv, tmp_affin, tmp_eg, T_amb)<<endl;
            }
        }else if(material_list[(*cell_mat)(0, i_e)].mat_type == "dielectric"){
            for(int l=0; l<n_vertex; ++l){
                //double tmp_dop = (*dop)(0,e_node[l]);
                double tmp_dop = 0;
                #if 0
                (*V_soln)(semi_mapping[l]) = Var_func->equilibrium_potential(tmp_dop,
                               tmp_ni, tmp_nc, tmp_nv, tmp_affin, tmp_eg, T, kb, e)+tmp_eg/2;
                #else
                if((*semiinsbc_node)(0,e_node[l])==0){
                    (*V_soln)(semi_mapping[l]) = -4.31549;
                    (*V_soln)(semi_mapping[l]+1) = 1e-10*pow(m,-3);
                    (*V_soln)(semi_mapping[l]+2) = 1e-10*pow(m,-3);
                }
                #endif
                #if 0
                (*V_soln)(semi_mapping[l]+1) = Var_func->dummyElectron((*V_soln)(semi_mapping[l]),
                               (*Fn_soln)(mapping[l]), false, tmp_affin, T, kb, e, tmp_nc);
                (*V_soln)(semi_mapping[l]+2) = Var_func->dummyHole((*V_soln)(semi_mapping[l]),
                               (*Fp_soln)(mapping[l]), tmp_eg, false, tmp_affin, T, kb, e, tmp_nv);
                double V=(*V_soln)(semi_mapping[l]); double n = (*V_soln)(semi_mapping[l]+1); double p= (*V_soln)(semi_mapping[l]+2);
                printf("V n p:%f,%f,%f\n",(*V_soln)(semi_mapping[l]),(*V_soln)(semi_mapping[l]+1),(*V_soln)(semi_mapping[l]+2));
                #else

                #endif
            }
        }

        if(get_InitialSol||get_PoissonSolcur){
            for(int l=0; l<n_vertex; ++l){
                (*V_soln)(semi_mapping[l]) = (*sol)(0,e_node[l]);
                (*V_soln)(semi_mapping[l]+1) = (*sol)(1,e_node[l])*pow(m,-3);
                (*V_soln)(semi_mapping[l]+2) = (*sol)(2,e_node[l])*pow(m,-3);
            }
        }
    }

    for(int i=0;i<num_n;i++){
        if((*semiinsbc_node)(0,i)==1){
            int index = dof_map_semi[i];
            (*V_soln)(index+3) = (*V_soln)(index);
            (*V_soln)(index+4) = 1e-10*pow(m,-3);
            (*V_soln)(index+5) = 1e-10*pow(m,-3);
        }
    }

    /** setup the source target **/
    tbox::Pointer<bound_condition> bound_source = new bound_condition();

    /** boundary condition source **/
    int num_bc = EBoundassignment.getSize();

    tbox::Array<double> bound_coe;
    tbox::Array<double> bound_n;
    tbox::Array<double> bound_p;
    bound_coe.resizeArray(9);
    bound_n.resizeArray(9);
    bound_p.resizeArray(9);

    for(int i_bc = 0; i_bc < num_bc; ++i_bc){
        if (patch_geo->hasEntitySet(EBoundassignment[i_bc].bc_no, hier::EntityUtilities::NODE)){
            const tbox::Array<int> &entity_idx = patch_geo->getEntityIndicesInSet(EBoundassignment[i_bc].bc_no,
                                                 hier::EntityUtilities::NODE,
                                                 num_n);
            /** load fixed volt value **/
            double app_volt = bound_source->source_assign(EBoundassignment[i_bc].sc_type,
                                                          time, dt,
                                                          EBoundassignment[i_bc].sc_coe);

            tbox::pout<<i_bc<<" bc type is "<<EBoundassignment[i_bc].bc_type<<" sc type is "<<EBoundassignment[i_bc].sc_type<<endl;
            tbox::pout<<"applied voltage is "<<app_volt<<" time is "<<time<<endl;

            int size = entity_idx.getSize();
            string contact_type;

            /** assign boundary value to solution vector **/
            if(EBoundassignment[i_bc].bc_type == "ohmic"){

                contact_type = "Ohmic_contact";
                for (int i = 0; i < size; ++i) {

                    /** setup material target **/
                    int dom_no;
                    dom_no = 0;
                    for(int i_m = 0; i_m < material_list.size(); ++i_m){
                        if(material_list[i_m].mat_type == "semiconductor"){
                            if((*node_mat)(i_m,entity_idx[i])==1){
                                dom_no = i_m;
                            }
                        }
                    }

                    tbox::Pointer<BaseMaterial<NDIM> >
                            mat_dom = mat_manager->getMaterial(material_list[dom_no].mat_name);

                    int index = dof_map[entity_idx[i]];
                    int index_s = dof_map_semi[entity_idx[i]];

                    bound_coe[0] = app_volt;
                    bound_coe[1] = T;
                    bound_coe[2] = mat_dom->effectiveAffinity(T_amb);
                    bound_coe[3] = mat_dom->bandGap(T_amb);
                    bound_coe[4] = mat_dom->effectiveDensityNc(T_amb)*pow(m,-3);
                    bound_coe[5] = mat_dom->effectiveDensityNv(T_amb)*pow(m,-3);
                    bound_coe[6] = (*dop)(0,entity_idx[i]);
                    bound_coe[7] = 0.0;
                    bound_coe[8] = mat_dom->intrinsicDensityNi(T_amb)*pow(m,-3);

                    bound_n[0] = app_volt;
                    bound_n[1] = T;
                    bound_n[2] = mat_dom->effectiveAffinity(T_amb);
                    bound_n[3] = mat_dom->bandGap(T_amb);
                    bound_n[4] = (*Fn_soln)(index);
                    bound_n[5] = (*Fp_soln)(index);
                    bound_n[6] = (*V_soln)(index_s);
                    bound_n[7] = (*dop)(0,entity_idx[i]);
                    bound_n[8] = mat_dom->intrinsicDensityNi(T_amb)*pow(m,-3);

                    bound_p[0] = app_volt;
                    bound_p[1] = T;
                    bound_p[2] = mat_dom->effectiveAffinity(T_amb);
                    bound_p[3] = mat_dom->bandGap(T_amb);
                    bound_p[4] = (*Fn_soln)(index);
                    bound_p[5] = (*Fp_soln)(index);
                    bound_p[6] = (*V_soln)(index_s);
                    bound_p[7] = (*dop)(0,entity_idx[i]);
                    bound_p[8] = mat_dom->intrinsicDensityNi(T_amb)*pow(m,-3);

                    (*V_soln)(index_s) = bound_source->e_bound_value(contact_type, time, dt, kb, e, bound_coe);
                    (*V_soln)(index_s+1) = bound_source->carn_bound_value(contact_type, car_bound_mode, time, dt, kb, e, bound_n);
                    (*V_soln)(index_s+2) = bound_source->carh_bound_value(contact_type, car_bound_mode, time, dt, kb, e, bound_p);
                    //cout<<"V_soln:"<<(*V_soln)(index_s)<<" "<<(*V_soln)(index_s+1)<<" "<<(*V_soln)(index_s+2)<<endl;
                }
            }else if(EBoundassignment[i_bc].bc_type == "gate"){
                for (int i = 0; i < size; ++i) {
                    contact_type = "Gate_contact";
                    int index_s = dof_map_semi[entity_idx[i]];
                    bound_coe[0] = app_volt;
                    bound_coe[1] = T;
                    bound_coe[2] = EBoundassignment[i_bc].sc_coe[3];
                    (*V_soln)(index_s) = bound_source->e_bound_value(contact_type, time, dt, kb, e, bound_coe);
                    (*V_soln)(index_s+1) = 1e-10*pow(m,-3);
                    (*V_soln)(index_s+2) = 1e-10*pow(m,-3);
                    double V=(*V_soln)(index_s); double n = (*V_soln)(index_s+1); double p= (*V_soln)(index_s+2);
                }
            }
        }
    }

    #if 0
    for (int i_e = 0; i_e < num_e; ++i_e){

        if(!((material_list[(*cell_mat)(0, i_e)].mat_type == "semiconductor")||(material_list[(*cell_mat)(0, i_e)].mat_type == "dielectric"))) continue;

        /** get node information **/
        int n_vertex = can_extent[i_e+1] - can_extent[i_e];
        tbox::Array<int> e_node(n_vertex);
        tbox::Array<int> mapping(n_vertex);
        tbox::Array<int> semi_mapping(n_vertex);

        for (int i1 = 0, j = can_extent[i_e]; i1 < n_vertex; ++i1, ++j) {
            mapping[i1] = dof_map[can_indices[j]];
            semi_mapping[i1] = dof_map_semi[can_indices[j]];
            e_node[i1] = can_indices[j];
        }

        if(get_InitialSol){
            for(int l=0; l<n_vertex; ++l){
                (*V_soln)(semi_mapping[l]) = (*sol)(0,e_node[l]);
                (*V_soln)(semi_mapping[l]+1) = (*sol)(1,e_node[l])*pow(m,-3);
                (*V_soln)(semi_mapping[l]+2) = (*sol)(2,e_node[l])*pow(m,-3);
                cout<<"sol_data:"<<(*sol_data)(0,node_id)<<" "<<(*sol_data)(1,node_id)<<" "<<(*sol_data)(2,node_id)<<endl;
            }
        }
    }
    #else
    if(get_InitialSol||get_PoissonSolcur){
        for(int i=0;i<num_n;i++){
            int index = dof_map_semi[i];
            (*V_soln)(index) = (*sol)(0,i);
            (*V_soln)(index+1) = (*sol)(1,i)*pow(m,-3);
            (*V_soln)(index+2) = (*sol)(2,i)*pow(m,-3);
            //cout<<"sol_data:"<<(*sol)(0,i)<<" "<<(*sol)(1,i)<<" "<<(*sol)(2,i)<<endl;
        }
    }
    #endif
}

/**
 * @brief NonlinearPoisson::neutral_status: neutral status for semiconductor domains
 * @param patch
 * @param time
 * @param dt
 * @param component_name
 */
void NonlinearPoisson::potential_bound(hier::Patch<NDIM> &patch,
                                       const double  time,
                                       const double dt,
                                       const string& component_name){

    /** gPatchGeometry **/
    tbox::Pointer< hier::PatchGeometry<NDIM> >
            patch_geo = patch.getPatchGeometry();
    /** get PatchTopology **/
    tbox::Pointer< hier::PatchTopology<NDIM> >
            patch_top = patch.getPatchTopology();

    tbox::Pointer<pdat::VectorData<NDIM, double> >
      V_soln = patch.getPatchData(d_poisson_sol_cur_id);

    tbox::Pointer<pdat::VectorData<NDIM, double> >
      Fn_soln = patch.getPatchData(d_fermin_sol_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
      Fp_soln = patch.getPatchData(d_fermih_sol_id);

    tbox::Pointer<pdat::NodeData<NDIM, double> >
      dop = patch.getPatchData(d_dopant_id);

    tbox::Pointer<pdat::NodeData<NDIM, int> >
      node_mat = patch.getPatchData(d_node_mat_id);

    /** 获取本地单元(边, 面)周围结点的索引关系 **/
    tbox::Array<int> can_extent, can_indices;
    patch_top->getCellAdjacencyNodes(can_extent, can_indices);

    /** material list **/
    tbox::Pointer<MaterialManager<NDIM> >
            mat_manager = MaterialManager<NDIM>::getManager();

    int * dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);
    int * dof_map_semi = d_dof_semi_info->getDOFMapping(patch, hier::EntityUtilities::NODE);

    int num_n = patch.getNumberOfNodes(1);

    /** setup the source target **/
    tbox::Pointer<bound_condition> bound_source = new bound_condition();

    /** boundary condition source **/
    int num_bc = EBoundassignment.getSize();

    tbox::Array<double> bound_coe;
    tbox::Array<double> bound_n;
    tbox::Array<double> bound_p;
    bound_coe.resizeArray(9);
    bound_n.resizeArray(9);
    bound_p.resizeArray(9);

    for(int i_bc = 0; i_bc < num_bc; ++i_bc){
        if (patch_geo->hasEntitySet(EBoundassignment[i_bc].bc_no, hier::EntityUtilities::NODE)){
            const tbox::Array<int> &entity_idx = patch_geo->getEntityIndicesInSet(EBoundassignment[i_bc].bc_no,
                                                 hier::EntityUtilities::NODE,
                                                 num_n);
            /** load fixed volt value **/
            double app_volt = bound_source->source_assign(EBoundassignment[i_bc].sc_type,
                                                          time, dt,
                                                          EBoundassignment[i_bc].sc_coe);

            tbox::pout<<i_bc<<" bc type is "<<EBoundassignment[i_bc].bc_type<<" sc type is "<<EBoundassignment[i_bc].sc_type<<endl;
            tbox::pout<<"applied voltage is "<<app_volt<<" time is "<<time<<endl;

            int size = entity_idx.getSize();
            string contact_type;

            /** assign boundary value to solution vector **/
            if(EBoundassignment[i_bc].bc_type == "ohmic"){

                contact_type = "Ohmic_contact";
                for (int i = 0; i < size; ++i) {

                    /** setup material target **/
                    int dom_no;
                    dom_no = 0;
                    for(int i_m = 0; i_m < material_list.size(); ++i_m){
                        if(material_list[i_m].mat_type == "semiconductor"){
                            if((*node_mat)(i_m,entity_idx[i])==1){
                                dom_no = i_m;
                            }
                        }
                    }

                    tbox::Pointer<BaseMaterial<NDIM> >
                            mat_dom = mat_manager->getMaterial(material_list[dom_no].mat_name);

                    int index = dof_map[entity_idx[i]];
                    int index_s = dof_map_semi[entity_idx[i]];

                    bound_coe[0] = app_volt; //app_volt=time
                    bound_coe[1] = T;
                    bound_coe[2] = mat_dom->effectiveAffinity(T_amb);
                    bound_coe[3] = mat_dom->bandGap(T_amb);
                    bound_coe[4] = mat_dom->effectiveDensityNc(T_amb)*pow(m,-3);
                    bound_coe[5] = mat_dom->effectiveDensityNv(T_amb)*pow(m,-3);
                    bound_coe[6] = (*dop)(0,entity_idx[i]);
                    bound_coe[7] = 0.0;
                    bound_coe[8] = mat_dom->intrinsicDensityNi(T_amb)*pow(m,-3);

                    bound_n[0] = app_volt;
                    bound_n[1] = T;
                    bound_n[2] = mat_dom->effectiveAffinity(T_amb);
                    bound_n[3] = mat_dom->bandGap(T_amb);
                    bound_n[4] = (*Fn_soln)(index);
                    bound_n[5] = (*Fp_soln)(index);
                    bound_n[6] = (*V_soln)(index_s);
                    bound_n[7] = (*dop)(0,entity_idx[i]);
                    bound_n[8] = mat_dom->intrinsicDensityNi(T_amb)*pow(m,-3);

                    bound_p[0] = app_volt;
                    bound_p[1] = T;
                    bound_p[2] = mat_dom->effectiveAffinity(T_amb);
                    bound_p[3] = mat_dom->bandGap(T_amb);
                    bound_p[4] = (*Fn_soln)(index);
                    bound_p[5] = (*Fp_soln)(index);
                    bound_p[6] = (*V_soln)(index_s);
                    bound_p[7] = (*dop)(0,entity_idx[i]);
                    bound_p[8] = mat_dom->intrinsicDensityNi(T_amb)*pow(m,-3);

                    (*V_soln)(index_s) = bound_source->e_bound_value(contact_type, time, dt, kb, e, bound_coe);
                    (*V_soln)(index_s+1) = bound_source->carn_bound_value(contact_type, car_bound_mode, time, dt, kb, e, bound_n);
                    (*V_soln)(index_s+2) = bound_source->carh_bound_value(contact_type, car_bound_mode, time, dt, kb, e, bound_p);
                }
            }else if(EBoundassignment[i_bc].bc_type == "gate"){
                for (int i = 0; i < size; ++i) {
                    contact_type = "Gate_contact";
                    int index_s = dof_map_semi[entity_idx[i]];
                    bound_coe[0] = app_volt;
                    bound_coe[1] = T;
                    bound_coe[2] = EBoundassignment[i_bc].sc_coe[3];
                    (*V_soln)(index_s) = bound_source->e_bound_value(contact_type, time, dt, kb, e, bound_coe);
                    (*V_soln)(index_s+1) = 1e-10*pow(m,-3);
                    (*V_soln)(index_s+2) = 1e-10*pow(m,-3);
                    double V=(*V_soln)(index_s); double n = (*V_soln)(index_s+1); double p= (*V_soln)(index_s+2);
                }
            }
        }
    }
}

/**
 * @brief NonlinearPoisson::reduceOnPatch ::  get system level varibales
 * @param vector
 * @param len
 * @param patch
 * @param time
 * @param dt
 * @param component_name
 */
void NonlinearPoisson::reduceOnPatch(double* vector, int len,
                                    hier::Patch<NDIM>& patch,
                                    const double time, const double dt,
                                    const string& component_name){
    if (component_name=="POISSON"){
        PotentialMaxError(vector, patch, time, dt, component_name);
    } else if(component_name=="GummelPOISSON"){
        PoissonMaxError(vector, patch, time, dt, component_name);
    } else if(component_name == "SEMI"){
        SemiMaxError(vector, len, patch, time, dt, component_name);
    } else if(component_name == "GummelSEMI"){
        GummelSemiMaxError(vector, len, patch, time, dt, component_name);
    } else if(component_name == "get_I"){
        getI(vector, len, patch, time, dt, component_name);
    } else if(component_name == "MAX_DOP"){
        MaxDop(vector, len, patch, time, dt, component_name);
    } else if(component_name == "saveSol"){
        saveSol(vector, len, patch, time, dt, component_name);
    }
}

/**
 * @brief NonlinearPoisson::PotentialMaxError: to find the maximum change in potential
 * @param vector
 * @param patch
 * @param time
 * @param dt
 * @param component_name
 */
void NonlinearPoisson::PotentialMaxError(double* vector,
                                    hier::Patch<NDIM>& patch,
                                    const double time, const double dt,
                                    const string& component_name){

    tbox::Pointer<pdat::VectorData<NDIM, double> >
        poisson_sol = patch.getPatchData(d_poisson_newton_sol_id);
    tbox::Pointer<pdat::NodeData<NDIM, int> >
            semiinsbc_node = patch.getPatchData(d_semiinsbc_node_id);
    tbox::Pointer<pdat::NodeData<NDIM, int> >
          siclass_node = patch.getPatchData(d_siclass_node_id);
    tbox::Pointer<pdat::NodeData<NDIM, int> >
          sio2_node = patch.getPatchData(d_sio2_node_id);

    int * dof_map_semi = d_dof_semi_info->getDOFMapping(patch, hier::EntityUtilities::NODE);

    int num_n = patch.getNumberOfNodes(1);
    double error;
    error = 0;
    for (int i_n = 0; i_n<num_n; ++i_n){
        if((*siclass_node)(0,i_n)==0) continue;
        int semi_mapping = dof_map_semi[i_n];
        double tmp_error = fabs((*poisson_sol)(semi_mapping));
        if(error<tmp_error){
            error = tmp_error;
        }
        if(((*sio2_node)(0,i_n)==1)&&((*semiinsbc_node)(0,i_n)==0)){
            int semi_mapping = dof_map_semi[i_n];
            double tmp_error = fabs((*poisson_sol)(semi_mapping));
            double tmp_error1 = fabs((*poisson_sol)(semi_mapping+1));
            double tmp_error2 = fabs((*poisson_sol)(semi_mapping+2));
            double tol = 1e-2;
            if((tmp_error1>tol)||(tmp_error2>tol)){
                (*poisson_sol)(semi_mapping+1) = 0;
                (*poisson_sol)(semi_mapping+2) = 0;
            }
        }
        if((*semiinsbc_node)(0,i_n)==1){
            int semi_mapping = dof_map_semi[i_n]+3;
            double tmp_error = fabs((*poisson_sol)(semi_mapping));
            double tmp_error1 = fabs((*poisson_sol)(semi_mapping+1));
            double tmp_error2 = fabs((*poisson_sol)(semi_mapping+2));
            double tol = 1e-2;
            if((tmp_error1>tol)||(tmp_error2>tol)){
                (*poisson_sol)(semi_mapping+1) = 0;
                (*poisson_sol)(semi_mapping+2) = 0;
            }
            if(error<tmp_error){
                error = tmp_error;
            }
        }
    }
    vector[0] = error;
}

/**
 * @brief NonlinearPoisson::SemiMaxError: to find the maximum change in quasi fermi level
 * @param vector
 * @param len
 * @param patch
 * @param time
 * @param dt
 * @param component_name
 */
void NonlinearPoisson::SemiMaxError(double* vector, int len,
                                    hier::Patch<NDIM>& patch,
                                    const double time, const double dt,
                                    const string& component_name){

    tbox::Pointer<pdat::VectorData<NDIM, double> >
        poisson_sol = patch.getPatchData(d_poisson_newton_sol_id);
    tbox::Pointer<pdat::NodeData<NDIM, int> >
            semiinsbc_node = patch.getPatchData(d_semiinsbc_node_id);
    tbox::Pointer<pdat::NodeData<NDIM, int> >
          siclass_node = patch.getPatchData(d_siclass_node_id);

    int * dof_map_semi = d_dof_semi_info->getDOFMapping(patch, hier::EntityUtilities::NODE);

    int num_n = patch.getNumberOfNodes(1);
    double error;
    error = 0;
    for (int i_n = 0; i_n<num_n; ++i_n){
        if((*siclass_node)(0,i_n)==0) continue;
        int semi_mapping = dof_map_semi[i_n];
        double tmp_error_n = fabs((*poisson_sol)(semi_mapping+1));
        double tmp_error_p = fabs((*poisson_sol)(semi_mapping+2));
        if(error<tmp_error_n){
            error = tmp_error_n;
        }else if(error<tmp_error_p){
            error = tmp_error_p;
        }
        if(error>5){
            //cout<<"error:"<<error<<" "<<"i_n:"<<i_n<<" "<<semi_mapping<<endl;
        }
        if((*semiinsbc_node)(0,i_n)==1){
            int semi_mapping = dof_map_semi[i_n]+3;
            double tmp_error_n = fabs((*poisson_sol)(semi_mapping+1));
            double tmp_error_p = fabs((*poisson_sol)(semi_mapping+2));
            if(error<tmp_error_n){
                error = tmp_error_n;
            }else if(error<tmp_error_p){
                error = tmp_error_p;
            }
            if(error>5){
                cout<<"error:"<<error<<" "<<"i_n:"<<i_n<<" "<<semi_mapping<<endl;
            }
        }
    }
    vector[0] = error;
}

/**
 * @brief NonlinearPoisson::MaxDop: to find the maximum doping
 * @param vector
 * @param len
 * @param patch
 * @param time
 * @param dt
 * @param component_name
 */
void NonlinearPoisson::MaxDop(double* vector, int len,
                                    hier::Patch<NDIM>& patch,
                                    const double time, const double dt,
                                    const string& component_name){

    tbox::Pointer<pdat::NodeData<NDIM, double> >
        dopant_data = patch.getPatchData(d_dopant_id);

    int num_n = patch.getNumberOfNodes(1);
    double MAX_DOP = 0;
    for (int i_n = 0; i_n<num_n; ++i_n){
        if(MAX_DOP<abs((*dopant_data)(0,i_n))){
            MAX_DOP=abs((*dopant_data)(0,i_n));
        }
    }
    vector[0] = MAX_DOP;
}

/**
 * @brief NonlinearPoisson::saveSol
 * @param vector
 * @param len
 * @param patch
 * @param time
 * @param dt
 * @param component_name
 */
void NonlinearPoisson::saveSol(double* vector, int len,
                                    hier::Patch<NDIM>& patch,
                                    const double time, const double dt,
                                    const string& component_name){

    tbox::Pointer<pdat::NodeData<NDIM, double> >
        dopant_data = patch.getPatchData(d_dopant_id);

    int num_n = patch.getNumberOfNodes(1);
    double MAX_DOP = 0;
    tbox::Pointer<pdat::VectorData<NDIM, double> >
      solution = patch.getPatchData(d_poisson_sol_cur_id);
    tbox::Pointer<pdat::NodeData<NDIM, double> >
      plotn_data = patch.getPatchData(d_carn_plot_id);
    tbox::Pointer<pdat::NodeData<NDIM, double> >
      logn_plot_data = patch.getPatchData(d_logn_plot_id);
    tbox::Pointer<pdat::NodeData<NDIM, double> >
      fermin_plot_data = patch.getPatchData(d_fermin_plot_id);

    tbox::Pointer<pdat::NodeData<NDIM, double> >
      carn_sol_1_data = patch.getPatchData(d_carn_sol_1_id);
    tbox::Pointer<pdat::NodeData<NDIM, double> >
      fn_sol_1_data = patch.getPatchData(d_fn_sol_1_id);

    tbox::Pointer<pdat::VectorData<NDIM, double> >
      solution_n = patch.getPatchData(d_carn_sol_cur_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
      logn_solution = patch.getPatchData(d_logn_sol_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
      fermin_solution = patch.getPatchData(d_fermin_sol_id);

    tbox::Pointer<pdat::NodeData<NDIM, double> >
      ploth_data = patch.getPatchData(d_carh_plot_id);
    tbox::Pointer<pdat::NodeData<NDIM, double> >
      logh_plot_data = patch.getPatchData(d_logh_plot_id);
    tbox::Pointer<pdat::NodeData<NDIM, double> >
      fermih_plot_data = patch.getPatchData(d_fermih_plot_id);

    tbox::Pointer<pdat::NodeData<NDIM, double> >
      carh_sol_1_data = patch.getPatchData(d_carh_sol_1_id);
    tbox::Pointer<pdat::NodeData<NDIM, double> >
      fh_sol_1_data = patch.getPatchData(d_fh_sol_1_id);

    tbox::Pointer<pdat::VectorData<NDIM, double> >
      solution_h = patch.getPatchData(d_carh_sol_cur_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
      logh_solution = patch.getPatchData(d_logh_sol_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
      fermih_solution = patch.getPatchData(d_fermih_sol_id);

    tbox::Pointer<pdat::NodeData<NDIM, int> >
      semi_node = patch.getPatchData(d_semi_node_id);
    tbox::Pointer<pdat::NodeData<NDIM, int> >
      sio2_node = patch.getPatchData(d_sio2_node_id);
    tbox::Pointer<pdat::NodeData<NDIM, int> >
      semiinsbc_node = patch.getPatchData(d_semiinsbc_node_id);

    tbox::Pointer<pdat::NodeData<NDIM, int> >
          siclass_node = patch.getPatchData(d_siclass_node_id);

    tbox::Pointer<pdat::NodeData<NDIM, double> >
            sol_data = patch.getPatchData(d_sol_id);

    int *dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);
    int *dof_map_semi = d_dof_semi_info->getDOFMapping(patch, hier::EntityUtilities::NODE);

    /** 取出本地PatchGeometry **/
    tbox::Pointer< hier::PatchGeometry<NDIM> > patch_geo = patch.getPatchGeometry();
    /** get PatchTopology **/
    tbox::Pointer< hier::PatchTopology<NDIM> > patch_top = patch.getPatchTopology();
    tbox::Pointer<pdat::NodeData<NDIM, double> >
      poi_sol_1_data = patch.getPatchData(d_phi_sol_1_id);

    int num = patch.getNumberOfNodes(0);
    if(save_carrier){
        ofstream streamout;
        streamout.open("carrier",ios::out);

        /** 取出本地Patch的结点坐标数组 **/
        tbox::Pointer<pdat::NodeData<NDIM, double> > node_coord = patch_geo->getNodeCoordinates();
        for(int i = 0; i < num; ++i){
            if((*siclass_node)(0,i)==0) continue;
            streamout<<(*node_coord)(0,i)<<" "<<(*node_coord)(1,i)<<" "<<(*node_coord)(2,i)<<" ";
            streamout<<(*carn_sol_1_data)(0,i)<<" "<<(*carh_sol_1_data)(0,i)<<" "<<(*poi_sol_1_data)(0,i)<<"\n";
            if((*semiinsbc_node)(0,i)==1){
                streamout<<(*node_coord)(0,i)<<" "<<(*node_coord)(1,i)<<" "<<(*node_coord)(2,i)<<" ";
                streamout<<(*carn_sol_1_data)(1,i)<<" "<<(*carh_sol_1_data)(1,i)<<" "<<(*poi_sol_1_data)(1,i)<<"\n";
            }
        }
        streamout.close();
    }

    if(save_sol){
        stringstream sol_file;
        sol_file << "sol_" << int(time/dt)<<"_"<<patch.getIndex();//patch.getPatchId()<<
        ofstream streamout_sol;
        streamout_sol.open(sol_file.str(),ios::out);
        if(!streamout_sol) cout<<"can't open file "<<sol_file.str()<<endl;
        streamout_sol<<"nd_id volt n p\n";
        for(int i = 0; i < num; ++i){
            if((*siclass_node)(0,i)==0) continue;
            int index_s = dof_map_semi[i];
            streamout_sol<<(*sol_data)(3,i)<<" "<<solution->getPointer()[index_s]
                       <<" "<<solution->getPointer()[index_s+1]/pow(m,-3)
                      <<" "<<solution->getPointer()[index_s+2]/pow(m,-3)<<"\n";
            if((*semiinsbc_node)(0,i)==1){
                index_s = dof_map_semi[i]+3;
                streamout_sol<<(*sol_data)(3,i)<<" "<<solution->getPointer()[index_s]
                           <<" "<<solution->getPointer()[index_s+1]/pow(m,-3)
                          <<" "<<solution->getPointer()[index_s+2]/pow(m,-3)<<"\n";
            }
        }
        streamout_sol.close();
    }
    //cout<<"patch.getIndex():"<<patch.getIndex()<<endl;
    vector[0] = MAX_DOP;
}

/**
 * @brief NonlinearPoisson::GummelSemiMaxError: to find the maximum change in quasi fermi level
 * @param vector
 * @param len
 * @param patch
 * @param time
 * @param dt
 * @param component_name
 */
void NonlinearPoisson::GummelSemiMaxError(double* vector, int len,
                                    hier::Patch<NDIM>& patch,
                                    const double time, const double dt,
                                    const string& component_name){

    #if 1
    tbox::Pointer<pdat::VectorData<NDIM, double> >
            carn_sol_cur_data = patch.getPatchData(d_carn_sol_cur_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
            carn_stetch_data = patch.getPatchData(d_carn_stetch_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
            carh_sol_cur_data = patch.getPatchData(d_carh_sol_cur_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
            carh_stetch_data = patch.getPatchData(d_carh_stetch_id);

    /** tbox::Pointer<pdat::VectorData<NDIM, double> >
            carh_sol_cur = patch.getPatchData(d_carh_sol_cur_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
            carh_sol_scr = patch.getPatchData(d_carh_stetch_id); **/

    tbox::Pointer<pdat::NodeData<NDIM, int> >
            semi_node = patch.getPatchData(d_semi_node_id);
    tbox::Pointer<pdat::NodeData<NDIM, int> >
            sio2_node = patch.getPatchData(d_sio2_node_id);

    tbox::Pointer<pdat::NodeData<NDIM, int> >
          siclass_node = patch.getPatchData(d_siclass_node_id);

    /** tbox::Pointer<pdat::VectorData<NDIM, double> >
            poisson_sol_cur = patch.getPatchData(d_poisson_sol_cur_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
            poisson_sol_scr = patch.getPatchData(d_poisson_stetch_id); **/

    tbox::Pointer<pdat::NodeData<NDIM,double> >
            node_plot_data = patch.getPatchData(d_node_var_plot_id);

    int * dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);

    tbox::Pointer<pdat::NodeData<NDIM, int> >
            semiinsbc_node = patch.getPatchData(d_semiinsbc_node_id);

    /** int * dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE); **/

    int num_n = patch.getNumberOfNodes(1);

    double error;
    error = 0;

    for (int i_n = 0; i_n < num_n; ++i_n){
        if((*siclass_node)(0,i_n)==0) continue;

        int mapping = dof_map[i_n];
        double tmp_error_n = fabs((*carn_sol_cur_data)(mapping) - (*carn_stetch_data)(mapping));
        double tmp_error_p = fabs((*carh_sol_cur_data)(mapping) - (*carh_stetch_data)(mapping));
        if(error<tmp_error_n){
            error = tmp_error_n;
        }else if(error<tmp_error_p){
            error = tmp_error_p;
        }
        if((*semiinsbc_node)(0,i_n)==1){
            int mapping = dof_map[i_n]+1;
            double tmp_error_n = fabs((*carn_sol_cur_data)(mapping) - (*carn_stetch_data)(mapping));
            double tmp_error_p = fabs((*carh_sol_cur_data)(mapping) - (*carh_stetch_data)(mapping));
            if(error<tmp_error_n){
                error = tmp_error_n;
            }else if(error<tmp_error_p){
                error = tmp_error_p;
            }
        }
    }
    vector[0] = error;
    #endif
    #if 0
    tbox::Pointer<pdat::VectorData<NDIM, double> >
            fermin_sol_cur = patch.getPatchData(d_fermin_sol_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
            fermin_sol_scr = patch.getPatchData(d_fermin_stetch_id);

    /** tbox::Pointer<pdat::VectorData<NDIM, double> >
            carh_sol_cur = patch.getPatchData(d_carh_sol_cur_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
            carh_sol_scr = patch.getPatchData(d_carh_stetch_id); **/

    tbox::Pointer<pdat::NodeData<NDIM, int> >
            semi_node = patch.getPatchData(d_semi_node_id);

    /** tbox::Pointer<pdat::VectorData<NDIM, double> >
            poisson_sol_cur = patch.getPatchData(d_poisson_sol_cur_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
            poisson_sol_scr = patch.getPatchData(d_poisson_stetch_id); **/

    tbox::Pointer<pdat::NodeData<NDIM,double> >
            node_plot_data = patch.getPatchData(d_node_var_plot_id);

    int * dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);

    /** int * dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE); **/

    int num_n = patch.getNumberOfNodes(1);

    double error;
    error = 0;

    for (int i_n = 0; i_n < num_n; ++i_n){
        if((*semi_node)(0,i_n) == 1){
            int mapping = dof_map[i_n];
            double tmp_elem = (*fermin_sol_cur)(mapping) - (*fermin_sol_scr)(mapping);
            double tmp_error = fabs(tmp_elem);
            (*node_plot_data)(0,i_n) = tmp_error;
            if(error<tmp_error){
                error = tmp_error;
                //cout<<"error:"<<error<<" "<<(*fermin_sol_cur)(mapping)<<" "<<(*fermin_sol_scr)(mapping);
            }
        }
    }
    vector[0] = error;
    #endif
}

/**
 * @brief NonlinearPoisson::getI: to find the current
 * @param vector
 * @param len
 * @param patch
 * @param time
 * @param dt
 * @param component_name
 */
void NonlinearPoisson::getI(double* vector, int len,
                            hier::Patch<NDIM>& patch,
                            const double time, const double dt,
                            const string& component_name){
    /** (形函数, 积分器, 单元)管理器 **/
    tbox::Pointer<ElementManager<NDIM> >
      ele_manager = ElementManager<NDIM>::getManager();
    tbox::Pointer<BaseElement<NDIM> >
      ele = ele_manager->getElement(d_element_type);

    /** 取出本地PatchGeometry **/
    tbox::Pointer< hier::PatchGeometry<NDIM> > patch_geo = patch.getPatchGeometry();
    /** 取出本地PatchTopology **/
    tbox::Pointer< hier::PatchTopology<NDIM> > patch_top = patch.getPatchTopology();
    /** 取出本地Patch的结点坐标数组 **/
    tbox::Pointer<pdat::NodeData<NDIM, double> >
      node_coord = patch_geo->getNodeCoordinates();

    /** get node values **/
    tbox::Pointer<pdat::VectorData<NDIM, double> >
      V_soln = patch.getPatchData(d_poisson_sol_cur_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
      N_soln = patch.getPatchData(d_carn_sol_cur_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
      H_soln = patch.getPatchData(d_carh_sol_cur_id);

    /** 取出本地cotact id值 **/
    tbox::Pointer<pdat::NodeData<NDIM, int> >
      contact_node_data = patch.getPatchData(d_contact_node_id);

    tbox::Pointer<pdat::CellData<NDIM, double> > bas_integ;
    tbox::Pointer<pdat::CellData<NDIM, double> > bas_rhs;
    if(save_bas){
        bas_integ = patch.getPatchData(d_bas_integ_id);
        bas_rhs = patch.getPatchData(d_bas_rhs_id);
    }

    tbox::Pointer<Variables_Functions>
            Var_func = new Variables_Functions();

    /** 获取本地单元(边, 面)周围结点的索引关系 **/
    tbox::Array<int> can_extent, can_indices;
    patch_top->getCellAdjacencyNodes(can_extent, can_indices);

    /** setup the source target **/
    tbox::Pointer<bound_condition> bound_source = new bound_condition();

    //int num_faces = patch.getNumberOfFaces(1);
    int num_cells = patch.getNumberOfCells(0);  //不能考虑影像区！否则电流会计算出错。

    /** 取出自由度映射信息 **/
    int * dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);
    int * dof_map_semi = d_dof_semi_info->getDOFMapping(patch, hier::EntityUtilities::NODE);

    /** boundary condition source **/
    int num_bc = EBoundassignment.getSize();

    /** cell - entity relationship **/
    tbox::Pointer<pdat::CellData<NDIM, int> >
            cell_mat = patch.getPatchData(d_cell_mat_id);

    /** material list **/
    tbox::Pointer<MaterialManager<NDIM> >
            mat_manager = MaterialManager<NDIM>::getManager();

    tbox::Array<double> current_coe;
    current_coe.resizeArray(num_bc);
    for(int i_bc=0; i_bc < num_bc; ++i_bc){
        current_coe[i_bc] = 0;
    }

    for(int i=0; i<num_cells; i++){
        if(!((material_list[(*cell_mat)(0, i)].mat_type == "semiconductor")||(material_list[(*cell_mat)(0, i)].mat_type == "dielectric"))) continue;

        int n_vertex = can_extent[i+1] - can_extent[i];
        tbox::Array<hier::DoubleVector<NDIM> > vertex(n_vertex);
        tbox::Array<int> mapping(n_vertex);
        tbox::Array<int> mapping_semi(n_vertex);
        tbox::Array<int> node_n(n_vertex);
        tbox::Array<int> contact_node(n_vertex);
        tbox::Array<double> T_node(n_vertex);
        tbox::Array<double> V_node(n_vertex);
        tbox::Array<double> N_node(n_vertex);
        tbox::Array<double> H_node(n_vertex);

        /** get cell node coord **/
        for (int i1 = 0, j = can_extent[i]; i1 < n_vertex; ++i1, ++j) {
            mapping[i1] = dof_map[can_indices[j]];
            mapping_semi[i1] = dof_map_semi[can_indices[j]];
            for(int k = 0; k< NDIM; ++k){
                vertex[i1][k] = (*node_coord)(k, can_indices[j])*m;
            }
            node_n[i1] = can_indices[j];
            contact_node[i1] = (*contact_node_data)(0, node_n[i1]);
        }

        for (int i_n = 0; i_n < n_vertex; ++i_n){
            T_node[i_n]=T;
            V_node[i_n]=(*V_soln)(mapping_semi[i_n]);
            if(material_list[(*cell_mat)(0, i)].mat_type == "semiconductor"){
                N_node[i_n]=(*V_soln)(mapping_semi[i_n]+1);
                H_node[i_n]=(*V_soln)(mapping_semi[i_n]+2);
            }
        }

        /// 初始化contact面单元的面积和外法向量
        tbox::Pointer<tbox::Array<double> > norm_vec = new tbox::Array<double>(3);
        tbox::Pointer<tbox::Array<double> > face_area = new tbox::Array<double>(1);
        (*face_area)[0]=0.0;
        for (int i = 0; i < 3; ++i) {
          (*norm_vec)[i] = 0.0;
        }

        for(int i_bc=0; i_bc<num_bc; i_bc++){
            int bc_no = EBoundassignment[i_bc].bc_no;
            for(int i_v=0; i_v<n_vertex; i_v++){
                if(contact_node[i_v]==bc_no){
                    tbox::Array<tbox::Array<int> > r_node = ele->get_r_node();
                    tbox::Pointer<tbox::Array<double> > edge_field = new tbox::Array<double >(cell_edge);
                    ele->calculateEdgeE(vertex, V_node, dt, time, edge_field);

                    /** get material for cell **/
                    tbox::Pointer<BaseMaterial<NDIM> > mat_cell =
                            mat_manager->getMaterial(material_list[(*cell_mat)(0,i)].mat_name);

                    tbox::Array<tbox::Array<double> > r_integ = getInteg(n_vertex,save_bas,i,ele,bas_integ,vertex,time,dt);

                    for(int i_e = 0; i_e < cell_edge; ++i_e){

                        int n_1 = r_node[i_e][0];
                        int n_2 = r_node[i_e][1];

                        double N1 = N_node[n_1];
                        double N2 = N_node[n_2];
                        double H1 = H_node[n_1];
                        double H2 = H_node[n_2];

                        double T_mid = (T_node[n_2] + T_node[n_1])/2;
                        double mun = mat_cell->electronMobility()*pow(m,2)/V/s;
                        double mup = mat_cell->holeMobility()*pow(m,2)/V/s;
                        //double affin = mat_cell->effectiveAffinity(T_mid);
                        //double eg = mat_cell->bandGap(T_mid);

                        double Vt = kb*T_mid/e;
                        double delt = -(V_node[n_2] - V_node[n_1])/Vt;

                        double cn_1;
                        double cn_2;
                        double cp_1;
                        double cp_2;
                        if(DD_type == "SUPG"){
                            /* SUPG */
                            double delt_v = -(V_node[n_2] - V_node[n_1]);

                            double art_dn = mun*Vt;
                            double art_dp = mup*Vt;
                            if(delt == 0.0){
                                art_dn = mun*Vt;
                                art_dp = mup*Vt;
                            } else{
                                art_dn = mun*Vt*(-delt/2)/tanh(-delt/2);
                                art_dp = mup*Vt*(-delt/2)/tanh(-delt/2);
                            }

                            cn_1 = -1*(0.5*mun*delt_v - art_dn)*(r_integ[i_e][i_v])*e;
                            cn_2 = (0.5*mun*delt_v + art_dn)*(r_integ[i_e][i_v])*e;
                            cp_1 = (0.5*mup*delt_v + art_dp)*(r_integ[i_e][i_v])*e;
                            cp_2 = -1*(0.5*mup*delt_v - art_dp)*(r_integ[i_e][i_v])*e;
                        }else if(DD_type=="SG"){
                            /* SG */
                            double tmp_bern_n1 = Var_func->Bern( delt);
                            double tmp_bern_n2 = Var_func->Bern(-delt);
                            double tmp_bern_p1 = Var_func->Bern(-delt);
                            double tmp_bern_p2 = Var_func->Bern( delt);

                            cn_1 = (mun*(kb*T_node[n_1]/e))*tmp_bern_n1*(r_integ[i_e][i_v])*e;
                            cn_2 = (mun*(kb*T_node[n_2]/e))*tmp_bern_n2*(r_integ[i_e][i_v])*e;
                            cp_1 = (mup*(kb*T_node[n_1]/e))*tmp_bern_p1*(r_integ[i_e][i_v])*e;
                            cp_2 = (mup*(kb*T_node[n_2]/e))*tmp_bern_p2*(r_integ[i_e][i_v])*e;
                        }

                        double Jn = N2 * cn_2 - N1 * cn_1;
                        double Jp = H1 * cp_1 - H2 * cp_2;

                        double J = Jn + Jp;
                        J = J / A;

                        current_coe[i_bc] = current_coe[i_bc] + J;
                    }
                }
            }
        }
    }

    for(int i_bc = 0; i_bc < num_bc; ++i_bc){
        /** output the current information **/
        vector[i_bc] = current_coe[i_bc];
    }
}

/**
 * @brief NonlinearPoisson::transportPoissonSol: potential-from current to stetch
 * @param patch
 * @param time
 * @param dt
 * @param component_name
 */
void NonlinearPoisson::transportPoissonSol(hier::Patch<NDIM>& patch,
                                        const double time,
                                        const double dt,
                                        const string& component_name){

    tbox::Pointer<pdat::VectorData<NDIM, double> >
      poisson_sol_cur = patch.getPatchData(d_poisson_sol_cur_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
      poisson_sol_delta = patch.getPatchData(d_poisson_newton_sol_id);
    tbox::Pointer<pdat::NodeData<NDIM, int> >
      semiinsbc_node = patch.getPatchData(d_semiinsbc_node_id);

    tbox::Pointer<pdat::NodeData<NDIM, int> >
          siclass_node = patch.getPatchData(d_siclass_node_id);

    int * dof_map_semi = d_dof_semi_info->getDOFMapping(patch, hier::EntityUtilities::NODE);

    int num_n = patch.getNumberOfNodes(1);

    for (int i_n = 0; i_n < num_n; ++i_n){
        if((*siclass_node)(0,i_n)==0) continue;
        int semi_mapping = dof_map_semi[i_n];
        double deta_phi = (*poisson_sol_delta)(semi_mapping);
        double deta_n = (*poisson_sol_delta)(semi_mapping+1);
        double deta_p = (*poisson_sol_delta)(semi_mapping+2);

        #if 1
        if(((*poisson_sol_cur)(semi_mapping+1) + deta_n)<=0){
            deta_n = 0;
        }
        if(((*poisson_sol_cur)(semi_mapping+2) + deta_p)<=0){
            deta_p = 0;
        }
        #endif

        (*poisson_sol_cur)(semi_mapping) = (*poisson_sol_cur)(semi_mapping) + deta_phi;
        (*poisson_sol_cur)(semi_mapping+1) = (*poisson_sol_cur)(semi_mapping+1) + deta_n;
        (*poisson_sol_cur)(semi_mapping+2) = (*poisson_sol_cur)(semi_mapping+2) + deta_p;

        if((*semiinsbc_node)(0,i_n)==1){
            int semi_mapping = dof_map_semi[i_n]+3;
            double deta_phi = (*poisson_sol_delta)(semi_mapping);
            double deta_n = (*poisson_sol_delta)(semi_mapping+1);
            double deta_p = (*poisson_sol_delta)(semi_mapping+2);

            #if 1
            if(((*poisson_sol_cur)(semi_mapping+1) + deta_n)<=0){
                deta_n = 0;
            }
            if(((*poisson_sol_cur)(semi_mapping+2) + deta_p)<=0){
                deta_p = 0;
            }
            #endif

            (*poisson_sol_cur)(semi_mapping) = (*poisson_sol_cur)(semi_mapping) + deta_phi;
            (*poisson_sol_cur)(semi_mapping+1) = (*poisson_sol_cur)(semi_mapping+1) + deta_n;
            (*poisson_sol_cur)(semi_mapping+2) = (*poisson_sol_cur)(semi_mapping+2) + deta_p;
        }
    }
}

/**
 * @brief NonlinearPoisson::transportSemiSol: fermi level-from current to stetch
 * @param patch
 * @param time
 * @param dt
 * @param component_name
 */
void NonlinearPoisson::transportSemiSol(hier::Patch<NDIM>& patch,
                                        const double time,
                                        const double dt,
                                        const string& component_name){
    tbox::Pointer<pdat::VectorData<NDIM, double> >
            fermin_sol_cur = patch.getPatchData(d_fermin_sol_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
            fermin_sol_scr = patch.getPatchData(d_fermin_stetch_id);

    tbox::Pointer<pdat::VectorData<NDIM, double> >
            carn_sol_cur_data = patch.getPatchData(d_carn_sol_cur_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
            carn_stetch_data = patch.getPatchData(d_carn_stetch_id);

    tbox::Pointer<pdat::VectorData<NDIM, double> >
            fermih_sol_cur = patch.getPatchData(d_fermih_sol_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
            fermih_sol_scr = patch.getPatchData(d_fermih_stetch_id);

    tbox::Pointer<pdat::VectorData<NDIM, double> >
            carh_sol_cur_data = patch.getPatchData(d_carh_sol_cur_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
            carh_stetch_data = patch.getPatchData(d_carh_stetch_id);

    tbox::Pointer<pdat::VectorData<NDIM, double> >
            poisson_sol_cur = patch.getPatchData(d_poisson_sol_cur_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
            poisson_sol_scr = patch.getPatchData(d_poisson_stetch_id);

    tbox::Pointer<pdat::NodeData<NDIM, int> >
            semi_node = patch.getPatchData(d_semi_node_id);
    tbox::Pointer<pdat::NodeData<NDIM, int> >
      semiinsbc_node = patch.getPatchData(d_semiinsbc_node_id);

    tbox::Pointer<pdat::NodeData<NDIM, int> >
          siclass_node = patch.getPatchData(d_siclass_node_id);

    int * dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);
    int * dof_map_semi = d_dof_semi_info->getDOFMapping(patch, hier::EntityUtilities::NODE);

    int num_n = patch.getNumberOfNodes(1);

    for (int i_n = 0; i_n < num_n; ++i_n){
        if((*siclass_node)(0,i_n)==0) continue;
        int mapping = dof_map[i_n];
        int semi_mapping = dof_map_semi[i_n];
        (*poisson_sol_scr)(semi_mapping) = (*poisson_sol_cur)(semi_mapping);
        (*poisson_sol_scr)(semi_mapping+1) = (*poisson_sol_cur)(semi_mapping+1);
        (*poisson_sol_scr)(semi_mapping+2) = (*poisson_sol_cur)(semi_mapping+2);

        (*fermin_sol_scr)(mapping) = (*fermin_sol_cur)(mapping);
        (*fermih_sol_scr)(mapping) = (*fermih_sol_cur)(mapping);

        (*carn_stetch_data)(mapping) = (*carn_sol_cur_data)(mapping);
        (*carh_stetch_data)(mapping) = (*carh_sol_cur_data)(mapping);
        if((*semiinsbc_node)(0,i_n)==1){
            int mapping = dof_map[i_n]+1;
            int semi_mapping = dof_map_semi[i_n]+3;
            (*poisson_sol_scr)(semi_mapping) = (*poisson_sol_cur)(semi_mapping);
            (*poisson_sol_scr)(semi_mapping+1) = (*poisson_sol_cur)(semi_mapping+1);
            (*poisson_sol_scr)(semi_mapping+2) = (*poisson_sol_cur)(semi_mapping+2);

            (*fermin_sol_scr)(mapping) = (*fermin_sol_cur)(mapping);
            (*fermih_sol_scr)(mapping) = (*fermih_sol_cur)(mapping);

            (*carn_stetch_data)(mapping) = (*carn_sol_cur_data)(mapping);
            (*carh_stetch_data)(mapping) = (*carh_sol_cur_data)(mapping);
        }
    }
    //cout<<"(*poisson_sol_cur)(4):"<<(*poisson_sol_cur)(4)<<endl;
}

/**
 * @brief NonlinearPoisson::transportGummelSemiSol: fermi level-from current to stetch
 * @param patch
 * @param time
 * @param dt
 * @param component_name
 */
void NonlinearPoisson::transportGummelSemiSol(hier::Patch<NDIM>& patch,
                                        const double time,
                                        const double dt,
                                        const string& component_name){
    tbox::Pointer<pdat::VectorData<NDIM, double> >
            carn_sol_cur_data = patch.getPatchData(d_carn_sol_cur_id);

    tbox::Pointer<pdat::VectorData<NDIM, double> >
            carh_sol_cur_data = patch.getPatchData(d_carh_sol_cur_id);

    tbox::Pointer<pdat::VectorData<NDIM, double> >
            poisson_sol_cur = patch.getPatchData(d_poisson_sol_cur_id);

    tbox::Pointer<pdat::NodeData<NDIM, int> >
            semi_node = patch.getPatchData(d_semi_node_id);
    tbox::Pointer<pdat::NodeData<NDIM, int> >
      semiinsbc_node = patch.getPatchData(d_semiinsbc_node_id);

    tbox::Pointer<pdat::NodeData<NDIM, int> >
          siclass_node = patch.getPatchData(d_siclass_node_id);

    int * dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);
    int * dof_map_semi = d_dof_semi_info->getDOFMapping(patch, hier::EntityUtilities::NODE);

    int num_n = patch.getNumberOfNodes(1);

    for (int i_n = 0; i_n < num_n; ++i_n){
        if((*siclass_node)(0,i_n)==0) continue;
        int mapping = dof_map[i_n];
        int semi_mapping = dof_map_semi[i_n];
        (*poisson_sol_cur)(semi_mapping+1) = (*carn_sol_cur_data)(mapping);
        (*poisson_sol_cur)(semi_mapping+2) = (*carh_sol_cur_data)(mapping);
        if((*semiinsbc_node)(0,i_n)==1){
            int mapping = dof_map[i_n]+1;
            int semi_mapping = dof_map_semi[i_n]+3;
            (*poisson_sol_cur)(semi_mapping+1) = (*carn_sol_cur_data)(mapping);
            (*poisson_sol_cur)(semi_mapping+2) = (*carh_sol_cur_data)(mapping);
        }
    }
}

/**
 * @brief NonlinearPoisson::stepSolPostProcess: log of carrier density
 * @param patch
 * @param time
 * @param dt
 * @param component_name
 */
void NonlinearPoisson::stepSolPostProcess(hier::Patch<NDIM>& patch,
                                          const double time,
                                          const double dt,
                                          const string& component_name){
    tbox::Pointer<pdat::VectorData<NDIM, double> >
      carn_sol_data = patch.getPatchData(d_carn_sol_cur_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
      carh_sol_data = patch.getPatchData(d_carh_sol_cur_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
      logn_sol_data = patch.getPatchData(d_logn_sol_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
      logh_sol_data = patch.getPatchData(d_logh_sol_id);

    tbox::Pointer<pdat::NodeData<NDIM, int> >
      semi_node = patch.getPatchData(d_semi_node_id);
    tbox::Pointer<pdat::NodeData<NDIM, int> >
      sio2_node = patch.getPatchData(d_sio2_node_id);
    tbox::Pointer<pdat::NodeData<NDIM, int> >
      semiinsbc_node = patch.getPatchData(d_semiinsbc_node_id);

    tbox::Pointer<pdat::NodeData<NDIM, int> >
          siclass_node = patch.getPatchData(d_siclass_node_id);

    tbox::Pointer<pdat::NodeData<NDIM, int> >
      heter_node = patch.getPatchData(d_heter_node_id);

    int * dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);

    int num_n = patch.getNumberOfNodes(1);

    for (int i = 0; i < num_n; ++i){
        if((*siclass_node)(0,i)==0) continue;
        int mapping = dof_map[i];
        (*logn_sol_data)(mapping) = log10(fabs((*carn_sol_data)(mapping))/pow(m,-3));
        (*logh_sol_data)(mapping) = log10(fabs((*carh_sol_data)(mapping))/pow(m,-3));
        if((*semiinsbc_node)(0,i)==1){
            int mapping = dof_map[i]+1;
            (*logn_sol_data)(mapping) = log10(fabs((*carn_sol_data)(mapping))/pow(m,-3));
            (*logh_sol_data)(mapping) = log10(fabs((*carh_sol_data)(mapping))/pow(m,-3));
        }
    }
}



/**
 * @brief NonlinearPoisson::getCellEntityRelation: cell to entity relationship
 * @param patch
 * @param component_name
 */
void NonlinearPoisson::getCellEntityRelation(hier::Patch<NDIM> & patch,
                                       const string& component_name){

    tbox::Pointer< pdat::CellData<NDIM, int> >
      cell_entity = patch.getPatchData(d_cell_entity_id);

    int num_mat = material_list.size();
    for (int im=0; im<num_mat; ++im){
        int num_entity = material_list[im].num_entity;
        for (int i_e = 0; i_e < num_entity; ++i_e){
            int id = material_list[im].entity_no[i_e];
            if (!patch.getPatchGeometry()->hasEntitySet
                (id, hier::EntityUtilities::CELL, patch.getNumberOfCells(1)))
                continue;
            tbox::Array<int> cells =  patch.getPatchGeometry()->getEntityIndicesInSet
                (id, hier::EntityUtilities::CELL, patch.getNumberOfCells(1));
            int num_cell = cells.getSize();
            for (int i=0; i<num_cell; ++i){
                (*cell_entity)(0,cells[i]) = id-1;
            }
        }
    }
}

/**
 * @brief NonlinearPoisson::getCellMatRelation: cell to material no. relationship
 * @param patch
 * @param component_name
 */
void NonlinearPoisson::getCellMatRelation(hier::Patch<NDIM> &patch,
                                          const string& component_name){
    tbox::Pointer< pdat::CellData<NDIM, int> >
      cell_mat = patch.getPatchData(d_cell_mat_id);

    int num_mat = material_list.size();
    for (int im=0; im<num_mat; ++im){
        int num_entity = material_list[im].num_entity;
        for (int i_e = 0; i_e < num_entity; ++i_e){
            int id = material_list[im].entity_no[i_e];
            if (!patch.getPatchGeometry()->hasEntitySet
                (id, hier::EntityUtilities::CELL, patch.getNumberOfCells(1)))
                continue;
            tbox::Array<int> cells =  patch.getPatchGeometry()->getEntityIndicesInSet
                (id, hier::EntityUtilities::CELL, patch.getNumberOfCells(1));
            int num_cell = cells.getSize();
            for (int i=0; i<num_cell; ++i){
                (*cell_mat)(0,cells[i]) = im;  /** material number starts from 0 **/
            }
        }
    }
}

/**
 * @brief node to entity relationship
 **/
/** this variable is now not that usefule due to that it does not contain
 *  boundary and vertex (shared nodes) informations. It should be motified
 *  into data structures such as 'chain table'[state--vertex/face/inner][num_of_entity][enetity no.]
 **/

void NonlinearPoisson::getNodeEntityRelation(hier::Patch<NDIM> &patch,
                                   const string& component_name){
    /** get patch information **/
    tbox::Pointer< hier::PatchTopology<NDIM> >
            patch_top = patch.getPatchTopology();

    tbox::Pointer< pdat::NodeData<NDIM, int> >
            node_entity = patch.getPatchData(d_node_entity_id);

    /** get indices information **/
    tbox::Array<int> can_extent, can_indices;
    patch_top->getCellAdjacencyNodes(can_extent, can_indices);

    int n_vertex = cell_point;

    int num_entity = entity_mat_list.size();
    for (int id=1; id<num_entity+1; ++id){
        if (!patch.getPatchGeometry()->hasEntitySet
            (id, hier::EntityUtilities::CELL, patch.getNumberOfCells(1)))
            continue;
        tbox::Array<int> cells =  patch.getPatchGeometry()->getEntityIndicesInSet
            (id, hier::EntityUtilities::CELL, patch.getNumberOfCells(1));

        int num_cell = cells.getSize();

        /** **/
        for (int i_e = 0; i_e < num_cell; ++i_e){
            int tmp_ie = cells[i_e];
            tbox::Array<int> node_n(n_vertex);
            for (int i1 = 0, j = can_extent[tmp_ie]; i1 < n_vertex; ++i1, ++j) {
                node_n[i1]=can_indices[j];
            }
            for (int i_n = 0; i_n < n_vertex; ++i_n){
                if((*node_entity)(0,node_n[i_n]) == 0){
                    (*node_entity)(0,node_n[i_n]) = id-1;
                }else if ((*node_entity)(0,node_n[i_n]) > 0 && (*node_entity)(0,node_n[i_n]) < id){
                    (*node_entity)(1,node_n[i_n]) = id-1;
                }
            }
        }
    }
}

/**
 * @brief NonlinearPoisson::getNodeMatRelation: node to material no. relationship
 * @param patch
 * @param component_name
 */
void NonlinearPoisson::getNodeMatRelation(hier::Patch<NDIM> &patch,
                                   const string& component_name){
    /** get patch information **/
    tbox::Pointer< hier::PatchTopology<NDIM> >
            patch_top = patch.getPatchTopology();

    tbox::Pointer< pdat::NodeData<NDIM, int> >
            node_mat = patch.getPatchData(d_node_mat_id);

    /** get indices information **/
    tbox::Array<int> can_extent, can_indices;
    patch_top->getCellAdjacencyNodes(can_extent, can_indices);

    int num_n = patch.getNumberOfNodes(1);
    int n_vertex = cell_point;

    /** initialization to be 0 **/
    for(int i_n = 0; i_n < num_n; ++i_n){
        for(int i_m = 0; i_m < n_mat; ++i_m){
            (*node_mat)(i_m,i_n) = 0;
        }
    }

    for (int id=0; id<n_mat; ++id){
        int num_d = material_list[id].num_entity;
        for (int ii_md = 0; ii_md < num_d; ++ii_md){
            if (!patch.getPatchGeometry()->hasEntitySet
                    (material_list[id].entity_no[ii_md], hier::EntityUtilities::CELL, patch.getNumberOfCells(1)))
                continue;

            tbox::Array<int> cells =  patch.getPatchGeometry()->getEntityIndicesInSet
                    (material_list[id].entity_no[ii_md], hier::EntityUtilities::CELL, patch.getNumberOfCells(1));

            int num_cell = cells.getSize();

            for (int i_e = 0; i_e < num_cell; ++i_e){
                int tmp_ie = cells[i_e];
                tbox::Array<int> node_n(n_vertex);
                for (int i1 = 0, j = can_extent[tmp_ie]; i1 < n_vertex; ++i1, ++j) {
                    node_n[i1]=can_indices[j];
                }

                for (int i_n = 0; i_n < n_vertex; ++i_n){
                    if((*node_mat)(id,node_n[i_n]) == 0){
                        (*node_mat)(id,node_n[i_n]) = 1;
                    }
                }
            }
        }
    }
}

void NonlinearPoisson::adTest(hier::Patch<NDIM> & patch,
                              const double  time,
                              const double  dt,
                              const string& component_name)
{
    adtl::AutoDScalar::numdir = 3*2;
    adtl::AutoDScalar V0 = 1; V0.setADValue(0,1.0);
    adtl::AutoDScalar n0 = 1; n0.setADValue(1,1.0);
    adtl::AutoDScalar p0 = 1; p0.setADValue(2,1.0);
    adtl::AutoDScalar V1 = 1; V1.setADValue(3,1.0);
    adtl::AutoDScalar n1 = 1; n1.setADValue(4,1.0);
    adtl::AutoDScalar p1 = 1; p1.setADValue(5,1.0);
    adtl::AutoDScalar F = 2*exp(V0)+n0+p0+V1+n1*n1;
    //cout<<F.getADValue(0)<<" "<<F.getADValue(1)<<" "<<F.getADValue(2)<<" "<<F.getADValue(3)<<" "<<F.getADValue(4)<<endl;
    double F_v = F.getValue();
    int i=0;
}

#if 0
/** **********************************************************************
 * @brief Newton poisson method
 ********************************************************************* **/

/**
 * @brief NonlinearPoisson::buildNewtonPMatOnPatch: Newton-Raphson
 * @param patch
 * @param time
 * @param dt
 * @param component_name
 */
void NonlinearPoisson::buildNewtonPMatOnPatch(hier::Patch<NDIM> & patch,
                                          const double  time,
                                          const double  dt,
                                          const string& component_name)
{
    /** manager for shapefuction, integrator & element **/
    tbox::Pointer<ElementManager<NDIM> >
            ele_manager = ElementManager<NDIM>::getManager();
    tbox::Pointer<BaseElement<NDIM> >
            ele = ele_manager->getElement(d_element_type);

    /** gPatchGeometry **/
    tbox::Pointer< hier::PatchGeometry<NDIM> >
            patch_geo = patch.getPatchGeometry();
    /** get PatchTopology **/
    tbox::Pointer< hier::PatchTopology<NDIM> >
            patch_top = patch.getPatchTopology();

    /** get Patch node coord **/
    tbox::Pointer<pdat::NodeData<NDIM, double> >
            node_coord = patch_geo->getNodeCoordinates();
    /** stiff matrix on patch **/
    tbox::Pointer<pdat::CSRMatrixData<NDIM,double> >
            mat_data = patch.getPatchData(d_poisson_newton_mat_id);
    /** right side vector **/
    tbox::Pointer<pdat::VectorData<NDIM,double> >
            vec_data = patch.getPatchData(d_poisson_newton_fval_id);

    tbox::Pointer<pdat::VectorData<NDIM, double> >
            V_soln = patch.getPatchData(d_poisson_sol_cur_id);

    /** cell - entity relationship **/
    tbox::Pointer<pdat::CellData<NDIM, int> >
            cell_mat = patch.getPatchData(d_cell_mat_id);
    /** cell - heter node relationship **/
    tbox::Pointer<pdat::CellData<NDIM, int> >
            cell_hnode = patch.getPatchData(d_cell_hnode_id);
    /** node - entity relationship **/
    tbox::Pointer<pdat::NodeData<NDIM, int> >
            heter_node = patch.getPatchData(d_heter_node_id);
    tbox::Pointer<pdat::NodeData<NDIM, double> >
            dop = patch.getPatchData(d_dopant_id);
    tbox::Pointer<pdat::NodeData<NDIM, int> >
            semiinsbc_node_data = patch.getPatchData(d_semiinsbc_node_id);

    /** quasi fermi level **/
    tbox::Pointer<pdat::VectorData<NDIM,double> >
      Fn_soln = patch.getPatchData(d_fermin_sol_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
      Fh_soln = patch.getPatchData(d_fermih_sol_id);

    tbox::Pointer<pdat::NodeData<NDIM, double> >
          carn_sol_1_data = patch.getPatchData(d_carn_sol_1_id);
    tbox::Pointer<pdat::NodeData<NDIM, double> >
          carh_sol_1_data = patch.getPatchData(d_carh_sol_1_id);
    tbox::Pointer<pdat::NodeData<NDIM, double> >
          V_plot_data = patch.getPatchData(d_poisson_plot_id);

    /** 获取本地单元(边, 面)周围结点的索引关系 **/
    tbox::Array<int> can_extent, can_indices;
    patch_top->getCellAdjacencyNodes(can_extent, can_indices);

    int num_cells = patch.getNumberOfCells(1);
    int num_nodes = patch.getNumberOfNodes(1);

    int * dof_map_semi = d_dof_semi_info->getDOFMapping(patch, hier::EntityUtilities::NODE);
    int * dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);
    double * fval_val = vec_data->getPointer();

    /** material list **/
    tbox::Pointer<MaterialManager<NDIM> >
            mat_manager = MaterialManager<NDIM>::getManager();

    /** variable function **/
    tbox::Pointer<Variables_Functions> Var_func = new Variables_Functions();
    Var_func->set_DD_type(DD_type);

    #if 0
    cout<<"m:"<<m<<endl;
    TBOX_ERROR("break");
    #endif

    adtl::AutoDScalar::numdir = 3*8;

    /** Traversing Unit **/
    for(int i = 0; i < num_cells; ++i){
        /** get material for cell **/
        tbox::Pointer<BaseMaterial<NDIM> > mat_cell =
                mat_manager->getMaterial(material_list[(*cell_mat)(0,i)].mat_name);
        //cout<<"mat_name:"<<material_list[(*cell_mat)(0,i)].mat_name<<endl;

        /** get node information **/
        int n_vertex = can_extent[i+1] - can_extent[i];
        tbox::Array<hier::DoubleVector<NDIM> > vertex(n_vertex);
        tbox::Array<int> node_n(n_vertex);
        tbox::Array<int> mapping_semi(n_vertex);
        tbox::Array<int> mapping(n_vertex);
        tbox::Array<double> value(n_vertex);
        tbox::Array<double> dop_node(n_vertex);
        tbox::Array<double> fermin_node(n_vertex);
        tbox::Array<double> fermip_node(n_vertex);

        for (int i1 = 0, j = can_extent[i]; i1 < n_vertex; ++i1, ++j) {
            mapping_semi[i1] = dof_map_semi[can_indices[j]];
            mapping[i1] = dof_map[can_indices[j]];
            for(int k = 0; k< NDIM; ++k){
                vertex[i1][k] = (*node_coord)(k, can_indices[j])*m;
            }
            if(material_list[(*cell_mat)(0,i)].mat_type == "semiconductor"){
                dop_node[i1]=(*dop)(0,can_indices[j]);
            }else {
                dop_node[i1]=0;
            }
            node_n[i1] = can_indices[j];
            fermin_node[i1] = (*Fn_soln)(mapping[i1]);
            fermip_node[i1] = (*Fh_soln)(mapping[i1]);
        }

        #if 0
        cout<<"cell:"<<i<<endl;
        for(int i1=0;i1<8;i1++){
            cout<<vertex[i1][0]<<" "<<vertex[i1][1]<<" "<<vertex[i1][2]<<"\n";
        }
        #endif

        double T_cell = T;

        double cell_epsr = mat_cell->electroPermittivity(T_cell);

        double tmp_ni = mat_cell->intrinsicDensityNi(T_amb)*pow(m,-3);
        double e_trap = mat_cell->electronTrap(T_cell);

        /*including
         * migration_CaugheyThomas_up(double Tl,double ds,double Fn,double Fp)
         * migration_Lombardi_up(double Tl,double Na,double Nd,double En,double Ep)
         * migration_Fletcher_up(double Tl,double n,double p)
         * migration_Arora_up(double Tl,double Na,double Nd)
         * migration_milv_up(double Tl)
         * electronMobility()/holeMobility()
        */
        double mun = mat_cell->electronMobility()*pow(m,2)/V/s;
        double mup = mat_cell->holeMobility()*pow(m,2)/V/s;


        tbox::Array<double> mun_node(n_vertex);
        tbox::Array<double> mup_node(n_vertex);
        for(int i1=0;i1<n_vertex;i1++){
            double dop_data = dop_node[i1]/pow(m,-3);
            double T_data = T/K;
            mun_node[i1]  = mat_cell->migration_milv_un(T)*pow(m,2)/V/s;
        }

        double Vt = kb*T/e;

        tbox::Array<tbox::Array<double> > r_integ = ele->buildCVFEMInteg(vertex, dt, time);

        tbox::Array<tbox::Array<int> > r_node = ele->get_r_node();

        /** initialize element matrix **/
        tbox::Pointer<tbox::Matrix<double> > ele_mat_E = new tbox::Matrix<double>();

        ele_mat_E->resize(n_vertex, n_vertex);

        for (int l = 0; l < n_vertex; ++l){
            for (int j = 0; j < n_vertex; ++j){
                (*ele_mat_E)(l,j) = 0.0;
            }
        }

        for(int i_v = 0; i_v < cell_point; ++i_v){
            for(int i_e = 0; i_e < cell_edge; ++i_e){
                int n_1 = r_node[i_e][0];
                int n_2 = r_node[i_e][1];

                double cn_1 = (r_integ[i_e][i_v]);
                double cn_2 = (r_integ[i_e][i_v]);

                (*ele_mat_E)(i_v,n_1) = (*ele_mat_E)(i_v,n_1) - cn_1;
                (*ele_mat_E)(i_v,n_2) = (*ele_mat_E)(i_v,n_2) + cn_2;
            }
        }

        /** element volume **/
        tbox::Pointer<tbox::Vector<double> > ele_volume = new tbox::Vector<double>();
        ele_volume->resize(n_vertex);

        for (int i_n = 0; i_n < n_vertex; ++i_n){
            value[i_n] = 1;
        }

        /** building element volume **/
        ele->buildCVFEMRhs(vertex, value, dt, time, ele_volume);

        #if 0
        for(int i_v = 0; i_v < cell_point; ++i_v){
            for(int i_n = 0; i_n < cell_point; ++i_n){
                cout<<"(*ele_mat_E)(i_v,i_n):"<<i_v<<":"<<i_n<<":"<<(*ele_mat_E)(i_v,i_n)<<endl;
            }
        }
        #endif

        tbox::Array<adtl::AutoDScalar> V_node(n_vertex);
        tbox::Array<adtl::AutoDScalar> n_node(n_vertex);
        tbox::Array<adtl::AutoDScalar> p_node(n_vertex);
        tbox::Array<adtl::AutoDScalar> fermi_n_node(n_vertex);
        tbox::Array<adtl::AutoDScalar> fermi_p_node(n_vertex);
        adtl::AutoDScalar fV_node[n_vertex];
        adtl::AutoDScalar fn_node[n_vertex];
        adtl::AutoDScalar fp_node[n_vertex];

        if(material_list[(*cell_mat)(0,i)].mat_type == "semiconductor"){
            for (int i1 = 0; i1 < n_vertex; ++i1) {
                V_node[i1] = (*V_soln)(mapping_semi[i1]); V_node[i1].setADValue(i1*3,1.0);
                n_node[i1] = (*V_soln)(mapping_semi[i1]+1); n_node[i1].setADValue(i1*3+1,1.0);
                p_node[i1] = (*V_soln)(mapping_semi[i1]+2); p_node[i1].setADValue(i1*3+2,1.0);

                //get quasi-fermi
                double nc_tmp = mat_cell->effectiveDensityNc(T_amb)*pow(m,-3);
                double nv_tmp = mat_cell->effectiveDensityNv(T_amb)*pow(m,-3);
                double affin_tmp = mat_cell->effectiveAffinity(T);
                double eg_tmp = mat_cell->bandGap(T);

                adtl::AutoDScalar n_tmp = Var_func->dummyElectron(V_node[i1], fermin_node[i1], false, affin_tmp, T, kb, e, nc_tmp);
                adtl::AutoDScalar p_tmp = Var_func->dummyHole(V_node[i1], fermip_node[i1], eg_tmp, false, affin_tmp, T, kb, e, nv_tmp);

                fermi_n_node[i1] = n_tmp; fermi_n_node[i1].setADValue(i1*3+1,1.0);
                fermi_p_node[i1] = p_tmp; fermi_p_node[i1].setADValue(i1*3+2,1.0);
            }

            for (int i1 = 0; i1 < n_vertex; ++i1) {
                /* get a b h */
                double a1[3];
                if(!transient_state||int(time/dt)==0){
                     for(int index=0;index<3;index++){
                         a1[index]=0;
                     }
                }else{
                    a1[0]=0;
                    a1[1]=-1;
                    a1[2]=1;
                }
                adtl::AutoDScalar R = 0;
                double U_radiation = 0;
                //R = Var_func->recombGenerationSRH(n_node[i1], p_node[i1], Vt, tmp_ni, e_trap);
                double b1[]={cell_epsr,1,1};
                adtl::AutoDScalar h1[]= {-(e/eps0)*(n_node[i1] - p_node[i1] - dop_node[i1]),
                                       U_radiation-R,
                                       R-U_radiation};

                /* get ut */
                double V_plot = (*V_plot_data)(0,node_n[i1]);
                double n_plot = (*carn_sol_1_data)(0,node_n[i1]);
                double p_plot = (*carh_sol_1_data)(0,node_n[i1]);

                //get fV
                fV_node[i1] = Var_func->get_fV(a1[0], b1[0], V_plot, dt,
                                               V_node,h1[0],
                                               (*ele_volume)[i1],
                                               i1,ele_mat_E);
                //get fn
                fn_node[i1] = Var_func->get_fn(a1[1], n_plot, dt,
                                               n_node, V_node, h1[1],
                                               (*ele_volume)[i1],i1,
                                               r_node,r_integ,
                                               mun, Vt, e);
                //get fp
                fp_node[i1] = Var_func->get_fp(a1[2], p_plot, dt,
                                               p_node, V_node, h1[2],
                                               (*ele_volume)[i1],i1,
                                               r_node, r_integ,
                                               mup, Vt, e);
            }
            for(int l = 0; l < n_vertex; ++l){
                int index_l = mapping_semi[l];
                for(int j = 0; j < n_vertex; ++j){
                    int index_j = mapping_semi[j];
                    //get dfv/dv,dfv/dn,dfv/dp
                    mat_data->addMatrixValue(index_l, index_j, fV_node[l].getADValue(3*j));
                    mat_data->addMatrixValue(index_l, index_j+1, fV_node[l].getADValue(3*j+1));
                    mat_data->addMatrixValue(index_l, index_j+2, fV_node[l].getADValue(3*j+2));
                    //get dfn/dv,dfn/dn,dfn/dp
                    mat_data->addMatrixValue(index_l+1, index_j, fn_node[l].getADValue(3*j));
                    mat_data->addMatrixValue(index_l+1, index_j+1, fn_node[l].getADValue(3*j+1));
                    mat_data->addMatrixValue(index_l+1, index_j+2, fn_node[l].getADValue(3*j+2));
                    //get dfp/dv,dfp/dn,dfp/dp
                    mat_data->addMatrixValue(index_l+2, index_j, fp_node[l].getADValue(3*j));
                    mat_data->addMatrixValue(index_l+2, index_j+1, fp_node[l].getADValue(3*j+1));
                    mat_data->addMatrixValue(index_l+2, index_j+2, fp_node[l].getADValue(3*j+2));
                }
                vec_data->addVectorValue(index_l, -fV_node[l].getValue());
                vec_data->addVectorValue(index_l+1, -fn_node[l].getValue());
                vec_data->addVectorValue(index_l+2, -fp_node[l].getValue());
            }
        }else if(material_list[(*cell_mat)(0,i)].mat_type == "dielectric"){
            for (int i1 = 0; i1 < n_vertex; ++i1) {
                //cout<<mapping_semi[i1]<<endl;
                if((*semiinsbc_node_data)(0,node_n[i1])==1){
                    //int index = mapping_semi[i1];
                    V_node[i1] = (*V_soln)(mapping_semi[i1]+3); V_node[i1].setADValue(i1*3,1.0);
                    n_node[i1] = (*V_soln)(mapping_semi[i1]+4); n_node[i1].setADValue(i1*3+1,1.0);
                    p_node[i1] = (*V_soln)(mapping_semi[i1]+5); p_node[i1].setADValue(i1*3+2,1.0);
                    //cout<<"semiinsbc_node_data,v n p:"<<V_node[i1]<<" "<<n_node[i1]<<" "<<p_node[i1]<<endl;
                }else{
                    V_node[i1] = (*V_soln)(mapping_semi[i1]); V_node[i1].setADValue(i1*3,1.0);
                    n_node[i1] = (*V_soln)(mapping_semi[i1]+1); n_node[i1].setADValue(i1*3+1,1.0);
                    p_node[i1] = (*V_soln)(mapping_semi[i1]+2); p_node[i1].setADValue(i1*3+2,1.0);
                }
            }
            for (int i1 = 0; i1 < n_vertex; ++i1) {
                /* get a b h */
                double a1[3];
                if(!transient_state||int(time/dt)==0){
                     for(int index=0;index<3;index++){
                         a1[index]=0;
                     }
                }else{
                    a1[0]=0;
                    a1[1]=-1;
                    a1[2]=1;
                }
                adtl::AutoDScalar R = 0;
                double U_radiation = 0;
                double Rd = 10;
                U_radiation = 0.01*8.1e12*Rd*pow(cm,-3)/s; //U=Y*g0*Rd;
                double b1[]={cell_epsr,1,1};
                adtl::AutoDScalar h1[]= {-(e/eps0)*(n_node[i1] - p_node[i1] - dop_node[i1]),
                                       U_radiation-R,
                                       R-U_radiation};

                /* get ut */
                double V_plot = (*V_plot_data)(0,node_n[i1]);
                double n_plot = (*carn_sol_1_data)(0,node_n[i1]);
                double p_plot = (*carh_sol_1_data)(0,node_n[i1]);
                if((*semiinsbc_node_data)(0,node_n[i1])==1){
                    n_plot = (*carn_sol_1_data)(1,node_n[i1]);
                    p_plot = (*carh_sol_1_data)(1,node_n[i1]);
                }

                //get fV
                fV_node[i1] = Var_func->get_fV(a1[0], b1[0], V_plot, dt,
                                               V_node,h1[0],
                                               (*ele_volume)[i1],
                                               i1,ele_mat_E);
                //get fn
                fn_node[i1] = Var_func->get_fn(a1[1], n_plot, dt,
                                               n_node, V_node, h1[1],
                                               (*ele_volume)[i1],i1,
                                               r_node,r_integ,
                                               mun, Vt, e);
                //get fp
                fp_node[i1] = Var_func->get_fp(a1[2], p_plot, dt,
                                               p_node, V_node, h1[2],
                                               (*ele_volume)[i1],i1,
                                               r_node, r_integ,
                                               mup, Vt, e);
            }
            for(int l = 0; l < n_vertex; ++l){
                int index_l = mapping_semi[l];
                if((*semiinsbc_node_data)(0,node_n[l])==1){
                    index_l = mapping_semi[l]+3;
                }
                for(int j = 0; j < n_vertex; ++j){
                    int index_j = mapping_semi[j];
                    if((*semiinsbc_node_data)(0,node_n[j])==1){
                        index_j = mapping_semi[j]+3;
                    }
                    //get dfv/dv,dfv/dn,dfv/dp
                    mat_data->addMatrixValue(mapping_semi[l], index_j, fV_node[l].getADValue(3*j));
                    mat_data->addMatrixValue(mapping_semi[l], index_j+1, fV_node[l].getADValue(3*j+1));
                    mat_data->addMatrixValue(mapping_semi[l], index_j+2, fV_node[l].getADValue(3*j+2));
                    //get dfn/dv,dfn/dn,dfn/dp
                    mat_data->addMatrixValue(index_l+1, index_j, fn_node[l].getADValue(3*j));
                    mat_data->addMatrixValue(index_l+1, index_j+1, fn_node[l].getADValue(3*j+1));
                    mat_data->addMatrixValue(index_l+1, index_j+2, fn_node[l].getADValue(3*j+2));
                    //get dfp/dv,dfp/dn,dfp/dp
                    mat_data->addMatrixValue(index_l+2, index_j, fp_node[l].getADValue(3*j));
                    mat_data->addMatrixValue(index_l+2, index_j+1, fp_node[l].getADValue(3*j+1));
                    mat_data->addMatrixValue(index_l+2, index_j+2, fp_node[l].getADValue(3*j+2));
                }
                vec_data->addVectorValue(mapping_semi[l], -fV_node[l].getValue());
                vec_data->addVectorValue(index_l+1, -fn_node[l].getValue());
                vec_data->addVectorValue(index_l+2, -fp_node[l].getValue());
                #if 0
                if(fV_node[l].getValue()>=1){
                    cout<<"sio2,fV_node[l].getValue(),"<<fV_node[l].getValue()<<endl;
                }
                #endif
            }
        }
    }
    /** assign boundary to si-sio2 nodes **/
    for (int i_n = 0; i_n < num_nodes; ++i_n){
        int index = dof_map_semi[i_n];
        if ((*semiinsbc_node_data)(0,i_n) == 1){
            mat_data->addMatrixValue(index+3,index+3,1.0);
            mat_data->addMatrixValue(index+3,index,-1.0);
            fval_val[index+3] = 0.0;
        }
    }
    mat_data->assemble();
}
#else
/**
 * @brief NonlinearPoisson::buildNewtonPMatOnPatch: Newton-Raphson
 * @param patch
 * @param time
 * @param dt
 * @param component_name
 */
void NonlinearPoisson::buildNewtonPMatOnPatch(hier::Patch<NDIM> & patch,
                                          const double  time,
                                          const double  dt,
                                          const string& component_name)
{
    /** manager for shapefuction, integrator & element **/
    tbox::Pointer<ElementManager<NDIM> >
            ele_manager = ElementManager<NDIM>::getManager();
    tbox::Pointer<BaseElement<NDIM> >
            ele = ele_manager->getElement(d_element_type);

    /** gPatchGeometry **/
    tbox::Pointer< hier::PatchGeometry<NDIM> >
            patch_geo = patch.getPatchGeometry();
    /** get PatchTopology **/
    tbox::Pointer< hier::PatchTopology<NDIM> >
            patch_top = patch.getPatchTopology();

    /** get Patch node coord **/
    tbox::Pointer<pdat::NodeData<NDIM, double> >
            node_coord = patch_geo->getNodeCoordinates();
    /** stiff matrix on patch **/
    tbox::Pointer<pdat::CSRMatrixData<NDIM,double> >
            mat_data = patch.getPatchData(d_poisson_newton_mat_id);
    /** right side vector **/
    tbox::Pointer<pdat::VectorData<NDIM,double> >
            vec_data = patch.getPatchData(d_poisson_newton_fval_id);

    tbox::Pointer<pdat::VectorData<NDIM, double> >
            V_soln = patch.getPatchData(d_poisson_sol_cur_id);

    /** cell - entity relationship **/
    tbox::Pointer<pdat::CellData<NDIM, int> >
            cell_mat = patch.getPatchData(d_cell_mat_id);
    /** cell - heter node relationship **/
    tbox::Pointer<pdat::CellData<NDIM, int> >
            cell_hnode = patch.getPatchData(d_cell_hnode_id);
    /** node - entity relationship **/
    tbox::Pointer<pdat::NodeData<NDIM, int> >
            heter_node = patch.getPatchData(d_heter_node_id);
    tbox::Pointer<pdat::NodeData<NDIM, double> >
            dop = patch.getPatchData(d_dopant_id);
    tbox::Pointer<pdat::NodeData<NDIM, int> >
            semiinsbc_node_data = patch.getPatchData(d_semiinsbc_node_id);

    /** quasi fermi level **/
    tbox::Pointer<pdat::VectorData<NDIM,double> >
      Fn_soln = patch.getPatchData(d_fermin_sol_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
      Fh_soln = patch.getPatchData(d_fermih_sol_id);

    tbox::Pointer<pdat::NodeData<NDIM, double> >
          carn_sol_1_data = patch.getPatchData(d_carn_sol_1_id);
    tbox::Pointer<pdat::NodeData<NDIM, double> >
          carh_sol_1_data = patch.getPatchData(d_carh_sol_1_id);
    tbox::Pointer<pdat::NodeData<NDIM, double> >
          V_plot_data = patch.getPatchData(d_poisson_plot_id);
    tbox::Pointer<pdat::CellData<NDIM, double> > Ex_data =
        patch.getPatchData(d_E_cell_id);

    /** 获取本地单元(边, 面)周围结点的索引关系 **/
    tbox::Array<int> can_extent, can_indices;
    patch_top->getCellAdjacencyNodes(can_extent, can_indices);

    int num_cells = patch.getNumberOfCells(1);
    int num_nodes = patch.getNumberOfNodes(1);

    int * dof_map_semi = d_dof_semi_info->getDOFMapping(patch, hier::EntityUtilities::NODE);
    int * dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);
    double * fval_val = vec_data->getPointer();

    /** material list **/
    tbox::Pointer<MaterialManager<NDIM> >
            mat_manager = MaterialManager<NDIM>::getManager();

    /** variable function **/
    tbox::Pointer<Variables_Functions> Var_func = new Variables_Functions();

    #if 0
    cout<<"m:"<<m<<endl;
    TBOX_ERROR("break");
    #endif

    adtl::AutoDScalar::numdir = 3*8;

    /** Traversing Unit **/
    for(int i = 0; i < num_cells; ++i){
        /** get material for cell **/
        tbox::Pointer<BaseMaterial<NDIM> > mat_cell =
                mat_manager->getMaterial(material_list[(*cell_mat)(0,i)].mat_name);
        //cout<<"mat_name:"<<material_list[(*cell_mat)(0,i)].mat_name<<endl;

        /** get node information **/
        int n_vertex = can_extent[i+1] - can_extent[i];
        tbox::Array<hier::DoubleVector<NDIM> > vertex(n_vertex);
        tbox::Array<int> node_n(n_vertex);
        tbox::Array<int> mapping_semi(n_vertex);
        tbox::Array<int> mapping(n_vertex);
        tbox::Array<double> value(n_vertex);
        tbox::Array<double> dop_node(n_vertex);
        tbox::Array<double> fermin_node(n_vertex);
        tbox::Array<double> fermip_node(n_vertex);

        for (int i1 = 0, j = can_extent[i]; i1 < n_vertex; ++i1, ++j) {
            mapping_semi[i1] = dof_map_semi[can_indices[j]];
            mapping[i1] = dof_map[can_indices[j]];
            for(int k = 0; k< NDIM; ++k){
                vertex[i1][k] = (*node_coord)(k, can_indices[j])*m;
            }
            if(material_list[(*cell_mat)(0,i)].mat_type == "semiconductor"){
                dop_node[i1]=(*dop)(0,can_indices[j]);
            }else {
                dop_node[i1]=0;
            }
            node_n[i1] = can_indices[j];
            fermin_node[i1] = (*Fn_soln)(mapping[i1]);
            fermip_node[i1] = (*Fh_soln)(mapping[i1]);
        }

        #if 0
        cout<<"cell:"<<i<<endl;
        for(int i1=0;i1<8;i1++){
            cout<<vertex[i1][0]<<" "<<vertex[i1][1]<<" "<<vertex[i1][2]<<"\n";
        }
        #endif

        double T_cell = T;

        double cell_epsr = mat_cell->electroPermittivity(T_cell);

        double tmp_ni = mat_cell->intrinsicDensityNi(T_amb)*pow(m,-3);
        double e_trap = mat_cell->electronTrap(T_cell);

        /*including
         * migration_CaugheyThomas_up(double Tl,double ds,double Fn,double Fp)
         * migration_Lombardi_up(double Tl,double Na,double Nd,double En,double Ep)
         * migration_Fletcher_up(double Tl,double n,double p)
         * migration_Arora_up(double Tl,double Na,double Nd)
         * migration_milv_up(double Tl)
         * electronMobility()/holeMobility()
        */
        double mun = mat_cell->electronMobility()*pow(m,2)/V/s;
        double mup = mat_cell->holeMobility()*pow(m,2)/V/s;


        tbox::Array<double> mun_node(n_vertex);
        tbox::Array<double> mup_node(n_vertex);
        for(int i1=0;i1<n_vertex;i1++){
            double dop_data = dop_node[i1]/pow(m,-3);
            double T_data = T/K;
            mun_node[i1]  = mat_cell->migration_milv_un(T)*pow(m,2)/V/s;
        }

        double Vt = kb*T/e;

        tbox::Array<tbox::Array<double> > r_integ = ele->buildCVFEMInteg(vertex, dt, time);

        tbox::Array<tbox::Array<int> > r_node = ele->get_r_node();

        /** initialize element matrix **/
        tbox::Pointer<tbox::Matrix<double> > ele_mat_E = new tbox::Matrix<double>();

        ele_mat_E->resize(n_vertex, n_vertex);

        for (int l = 0; l < n_vertex; ++l){
            for (int j = 0; j < n_vertex; ++j){
                (*ele_mat_E)(l,j) = 0.0;
            }
        }

        for(int i_v = 0; i_v < cell_point; ++i_v){
            for(int i_e = 0; i_e < cell_edge; ++i_e){
                int n_1 = r_node[i_e][0];
                int n_2 = r_node[i_e][1];

                double cn_1 = (r_integ[i_e][i_v]);
                double cn_2 = (r_integ[i_e][i_v]);

                (*ele_mat_E)(i_v,n_1) = (*ele_mat_E)(i_v,n_1) - cn_1;
                (*ele_mat_E)(i_v,n_2) = (*ele_mat_E)(i_v,n_2) + cn_2;
            }
        }

        /** element volume **/
        tbox::Pointer<tbox::Vector<double> > ele_volume = new tbox::Vector<double>();
        ele_volume->resize(n_vertex);

        for (int i_n = 0; i_n < n_vertex; ++i_n){
            value[i_n] = 1;
        }

        /** building element volume **/
        ele->buildCVFEMRhs(vertex, value, dt, time, ele_volume);

        #if 0
        for(int i_v = 0; i_v < cell_point; ++i_v){
            for(int i_n = 0; i_n < cell_point; ++i_n){
                cout<<"(*ele_mat_E)(i_v,i_n):"<<i_v<<":"<<i_n<<":"<<(*ele_mat_E)(i_v,i_n)<<endl;
            }
        }
        #endif

        tbox::Array<double> V_node(n_vertex);
        tbox::Array<double> n_node(n_vertex);
        tbox::Array<double> p_node(n_vertex);
        tbox::Array<adtl::AutoDScalar> fermi_n_node(n_vertex);
        tbox::Array<adtl::AutoDScalar> fermi_p_node(n_vertex);

        if(material_list[(*cell_mat)(0,i)].mat_type == "semiconductor"){
            for (int i1 = 0; i1 < n_vertex; ++i1) {
                V_node[i1] = (*V_soln)(mapping_semi[i1]);
                n_node[i1] = (*V_soln)(mapping_semi[i1]+1);
                p_node[i1] = (*V_soln)(mapping_semi[i1]+2);

                double V = (*V_soln)(mapping_semi[i1]); double n = (*V_soln)(mapping_semi[i1]+1); double p = (*V_soln)(mapping_semi[i1]+2);

                //get quasi-fermi
                double nc_tmp = mat_cell->effectiveDensityNc(T_amb)*pow(m,-3);
                double nv_tmp = mat_cell->effectiveDensityNv(T_amb)*pow(m,-3);
                double affin_tmp = mat_cell->effectiveAffinity(T);
                double eg_tmp = mat_cell->bandGap(T);

                adtl::AutoDScalar n_tmp = Var_func->dummyElectron(V_node[i1], fermin_node[i1], false, affin_tmp, T, kb, e, nc_tmp);
                adtl::AutoDScalar p_tmp = Var_func->dummyHole(V_node[i1], fermip_node[i1], eg_tmp, false, affin_tmp, T, kb, e, nv_tmp);

                fermi_n_node[i1] = n_tmp; fermi_n_node[i1].setADValue(i1*3+1,1.0);
                fermi_p_node[i1] = p_tmp; fermi_p_node[i1].setADValue(i1*3+2,1.0);
            }

            for (int i1 = 0; i1 < n_vertex; ++i1) {
                /* get a b h */
                double a1[3];
                if(!transient_state||int(time/dt)==0){
                     for(int index=0;index<3;index++){
                         a1[index]=0;
                     }
                }else{
                    a1[0]=0;
                    a1[1]=-1;
                    a1[2]=1;
                }
                adtl::AutoDScalar R = 0;
                double U_radiation = 0;
                //R = Var_func->recombGenerationSRH(n_node[i1], p_node[i1], Vt, tmp_ni, e_trap);

                //zxy
                double T_in = T/K;
                double e_field = sqrt(pow((*Ex_data)(0,i),2)+pow((*Ex_data)(1,i),2)+pow((*Ex_data)(2,i),2));

                double n_des = n_node[i1]/pow(m,-3);
                double p_des = p_node[i1]/pow(m,-3);
                double mun_node = mun / (pow(m,2)/V/s);
                double mup_node = mup / (pow(m,2)/V/s);
                double G = 0;//mat_cell->AvalancheGeneration( alpha_n,  alpha_p, n_des, p_des, e_field, mun_node, mup_node)*pow(m,-3)*pow(s,-1);
                //zxy

                double b1[]={cell_epsr,1,1};
                adtl::AutoDScalar n_node_i1 = n_node[i1]; n_node_i1.setADValue(i1*3+1,1.0);
                adtl::AutoDScalar p_node_i1 = p_node[i1]; p_node_i1.setADValue(i1*3+2,1.0);
                adtl::AutoDScalar h1[]= {-(e/eps0)*(n_node_i1 - p_node_i1 - dop_node[i1]),
                                       G+U_radiation-R,
                                       R-G-U_radiation};

                /* get ut */
                double V_plot = (*V_plot_data)(0,node_n[i1]);
                double n_plot = (*carn_sol_1_data)(0,node_n[i1]);
                double p_plot = (*carh_sol_1_data)(0,node_n[i1]);

                //get fV
                adtl::AutoDScalar fV_node = Var_func->get_fV(a1[0], b1[0], V_plot, dt,
                                               V_node,h1[0],
                                               (*ele_volume)[i1],
                                               i1,ele_mat_E);
                //get fn
                adtl::AutoDScalar fn_node = Var_func->get_fn(a1[1], n_plot, dt,
                                               n_node, V_node, h1[1],
                                               (*ele_volume)[i1],i1,
                                               r_node,r_integ,
                                               mun, Vt, e);
                //get fp
                adtl::AutoDScalar fp_node = Var_func->get_fp(a1[2], p_plot, dt,
                                               p_node, V_node, h1[2],
                                               (*ele_volume)[i1],i1,
                                               r_node, r_integ,
                                               mup, Vt, e);
                int l=i1;
                int index_l = mapping_semi[l];
                for(int j = 0; j < n_vertex; ++j){
                    int index_j = mapping_semi[j];
                    //get dfv/dv,dfv/dn,dfv/dp
                    mat_data->addMatrixValue(index_l, index_j, fV_node.getADValue(3*j));
                    mat_data->addMatrixValue(index_l, index_j+1, fV_node.getADValue(3*j+1));
                    mat_data->addMatrixValue(index_l, index_j+2, fV_node.getADValue(3*j+2));
                    //get dfn/dv,dfn/dn,dfn/dp
                    mat_data->addMatrixValue(index_l+1, index_j, fn_node.getADValue(3*j));
                    mat_data->addMatrixValue(index_l+1, index_j+1, fn_node.getADValue(3*j+1));
                    mat_data->addMatrixValue(index_l+1, index_j+2, fn_node.getADValue(3*j+2));
                    //get dfp/dv,dfp/dn,dfp/dp
                    mat_data->addMatrixValue(index_l+2, index_j, fp_node.getADValue(3*j));
                    mat_data->addMatrixValue(index_l+2, index_j+1, fp_node.getADValue(3*j+1));
                    mat_data->addMatrixValue(index_l+2, index_j+2, fp_node.getADValue(3*j+2));
                }
                vec_data->addVectorValue(index_l, -fV_node.getValue());
                vec_data->addVectorValue(index_l+1, -fn_node.getValue());
                vec_data->addVectorValue(index_l+2, -fp_node.getValue());
            }
        }else if(material_list[(*cell_mat)(0,i)].mat_type == "dielectric"){
            for (int i1 = 0; i1 < n_vertex; ++i1) {
                //cout<<mapping_semi[i1]<<endl;
                if((*semiinsbc_node_data)(0,node_n[i1])==1){
                    //int index = mapping_semi[i1];
                    V_node[i1] = (*V_soln)(mapping_semi[i1]+3);
                    n_node[i1] = (*V_soln)(mapping_semi[i1]+4);
                    p_node[i1] = (*V_soln)(mapping_semi[i1]+5);
                    //cout<<"semiinsbc_node_data,v n p:"<<V_node[i1]<<" "<<n_node[i1]<<" "<<p_node[i1]<<endl;
                }else{
                    V_node[i1] = (*V_soln)(mapping_semi[i1]);
                    n_node[i1] = (*V_soln)(mapping_semi[i1]+1);
                    p_node[i1] = (*V_soln)(mapping_semi[i1]+2);
                }
            }
            for (int i1 = 0; i1 < n_vertex; ++i1) {
                /* get a b h */
                double a1[3];
                if(!transient_state||int(time/dt)==0){
                     for(int index=0;index<3;index++){
                         a1[index]=0;
                     }
                }else{
                    a1[0]=0;
                    a1[1]=-1;
                    a1[2]=1;
                }
                adtl::AutoDScalar R = 0;
                double U_radiation = 0;
                double Rd = 10;
                U_radiation = 0.01*8.1e12*Rd*pow(cm,-3)/s; //U=Y*g0*Rd;
                double b1[]={cell_epsr,1,1};
                adtl::AutoDScalar h1[]= {-(e/eps0)*(n_node[i1] - p_node[i1] - dop_node[i1]),
                                       U_radiation-R,
                                       R-U_radiation};

                /* get ut */
                double V_plot = (*V_plot_data)(0,node_n[i1]);
                double n_plot = (*carn_sol_1_data)(0,node_n[i1]);
                double p_plot = (*carh_sol_1_data)(0,node_n[i1]);
                if((*semiinsbc_node_data)(0,node_n[i1])==1){
                    n_plot = (*carn_sol_1_data)(1,node_n[i1]);
                    p_plot = (*carh_sol_1_data)(1,node_n[i1]);
                }

                //get fV
                adtl::AutoDScalar fV_node = Var_func->get_fV(a1[0], b1[0], V_plot, dt,
                                               V_node,h1[0],
                                               (*ele_volume)[i1],
                                               i1,ele_mat_E);
                //get fn
                adtl::AutoDScalar fn_node = Var_func->get_fn(a1[1], n_plot, dt,
                                               n_node, V_node, h1[1],
                                               (*ele_volume)[i1],i1,
                                               r_node,r_integ,
                                               mun, Vt, e);
                //get fp
                adtl::AutoDScalar fp_node = Var_func->get_fp(a1[2], p_plot, dt,
                                               p_node, V_node, h1[2],
                                               (*ele_volume)[i1],i1,
                                               r_node, r_integ,
                                               mup, Vt, e);

                int l = i;
                int index_l = mapping_semi[l];
                if((*semiinsbc_node_data)(0,node_n[l])==1){
                    index_l = mapping_semi[l]+3;
                }
                for(int j = 0; j < n_vertex; ++j){
                    int index_j = mapping_semi[j];
                    if((*semiinsbc_node_data)(0,node_n[j])==1){
                        index_j = mapping_semi[j]+3;
                    }
                    //get dfv/dv,dfv/dn,dfv/dp
                    mat_data->addMatrixValue(mapping_semi[l], index_j, fV_node.getADValue(3*j));
                    mat_data->addMatrixValue(mapping_semi[l], index_j+1, fV_node.getADValue(3*j+1));
                    mat_data->addMatrixValue(mapping_semi[l], index_j+2, fV_node.getADValue(3*j+2));
                    //get dfn/dv,dfn/dn,dfn/dp
                    mat_data->addMatrixValue(index_l+1, index_j, fn_node.getADValue(3*j));
                    mat_data->addMatrixValue(index_l+1, index_j+1, fn_node.getADValue(3*j+1));
                    mat_data->addMatrixValue(index_l+1, index_j+2, fn_node.getADValue(3*j+2));
                    //get dfp/dv,dfp/dn,dfp/dp
                    mat_data->addMatrixValue(index_l+2, index_j, fp_node.getADValue(3*j));
                    mat_data->addMatrixValue(index_l+2, index_j+1, fp_node.getADValue(3*j+1));
                    mat_data->addMatrixValue(index_l+2, index_j+2, fp_node.getADValue(3*j+2));
                }
                vec_data->addVectorValue(mapping_semi[l], -fV_node.getValue());
                vec_data->addVectorValue(index_l+1, -fn_node.getValue());
                vec_data->addVectorValue(index_l+2, -fp_node.getValue());
            }
            for(int l = 0; l < n_vertex; ++l){
                #if 0
                if(fV_node[l].getValue()>=1){
                    cout<<"sio2,fV_node[l].getValue(),"<<fV_node[l].getValue()<<endl;
                }
                #endif
            }
        }
    }
    /** assign boundary to si-sio2 nodes **/
    for (int i_n = 0; i_n < num_nodes; ++i_n){
        int index = dof_map_semi[i_n];
        if ((*semiinsbc_node_data)(0,i_n) == 1){
            mat_data->addMatrixValue(index+3,index+3,1.0);
            mat_data->addMatrixValue(index+3,index,-1.0);
            fval_val[index+3] = 0.0;
        }
    }
    mat_data->assemble();
}
#endif


/**
 * @brief NonlinearPoisson::setNewtonPPhysicalBC: matrix
 * @param patch
 * @param time
 * @param dt
 * @param component_name
 */
void NonlinearPoisson::setNewtonPPhysicalBC(hier::Patch<NDIM> & patch,
                          const double  time,
                          const double  dt,
                          const string& component_name){
    /** (形函数, 积分器, 单元)管理器 **/
    tbox::Pointer<ElementManager<NDIM> >
            ele_manager = ElementManager<NDIM>::getManager();

    /** 取出本地PatchGeometry **/
    tbox::Pointer< hier::PatchGeometry<NDIM> > patch_geo = patch.getPatchGeometry();
    /** get PatchTopology **/
    tbox::Pointer< hier::PatchTopology<NDIM> > patch_top = patch.getPatchTopology();
    /** 获取网格片矩阵对象 **/
    tbox::Pointer<pdat::CSRMatrixData<NDIM,double> >
            mat_data = patch.getPatchData(d_poisson_newton_mat_id);

    /** node - entity relationship **/
    tbox::Pointer<pdat::NodeData<NDIM, int> >
            heter_node = patch.getPatchData(d_heter_node_id);

    /** get node values **/
    tbox::Pointer<pdat::VectorData<NDIM, double> >
            V_soln = patch.getPatchData(d_poisson_sol_cur_id);

    tbox::Array<int> face_node_ext, face_node_idx;
    patch_top->getFaceAdjacencyNodes(face_node_ext, face_node_idx);

    int num_nodes = patch.getNumberOfNodes(1);

    /** obtain the right side vector **/
    tbox::Pointer<pdat::VectorData<NDIM,double> >
            fval_data = patch.getPatchData(d_poisson_newton_fval_id);
    double * fval_val = fval_data->getPointer();

    int * dof_map_semi = d_dof_semi_info->getDOFMapping(patch, hier::EntityUtilities::NODE);
    /** 取出矩阵向量的核心数据结构的指针 **/
    int * row_start = mat_data->getRowStartPointer();
    int * col_idx   = mat_data->getColumnIndicesPointer();
    double * mat_val    = mat_data->getValuePointer();

    /** 取出本地Patch的结点坐标数组 **/
    tbox::Pointer<pdat::NodeData<NDIM, double> > node_coord = patch_geo->getNodeCoordinates();

    /** setup the source target **/
    tbox::Pointer<bound_condition> bound_source = new bound_condition();

    /** assign boundary condition to system matrix **/
    int num_bc = EBoundassignment.getSize();
    for(int i_bc=0; i_bc<num_bc; ++i_bc){
        if (patch_geo->hasEntitySet(EBoundassignment[i_bc].bc_no, hier::EntityUtilities::NODE)){
            const tbox::Array<int> &entity_idx = patch_geo->getEntityIndicesInSet(EBoundassignment[i_bc].bc_no,
                                                 hier::EntityUtilities::NODE,
                                                 num_nodes);
            if ((EBoundassignment[i_bc].bc_type == "ohmic")||(EBoundassignment[i_bc].bc_type == "gate")){
                int size = entity_idx.getSize();
                for (int i = 0; i < size; ++i) {
                    int index = dof_map_semi[entity_idx[i]];
                    fval_val[index] = 0.0;
                    fval_val[index+1] = 0.0;
                    fval_val[index+2] = 0.0;
                    /** <mat_value at boundary> **/
                    for (int j = row_start[index]; j < row_start[index+1]; ++j){
                        //cout<<mat_val[j]<<" ";
                        if (col_idx[j] == index){
                            mat_val[j] = 1.0;   /** 对角线元素 **/
                        }
                        else {
                            mat_val[j] = 0.0;   /**< 非对角线元素 */
                        }
                    }
                    //cout<<endl;
                    int index_1 = dof_map_semi[entity_idx[i]]+1;
                    /** <mat_value at boundary> **/
                    for (int j = row_start[index_1]; j < row_start[index_1+1]; ++j){
                        //cout<<mat_val[j]<<" ";
                        if (col_idx[j] == index_1){
                            mat_val[j] = 1.0;   /** 对角线元素 **/
                        }
                        else {
                            mat_val[j] = 0.0;   /**< 非对角线元素 */
                        }
                    }
                    int index_2 = dof_map_semi[entity_idx[i]]+2;
                    /** <mat_value at boundary> **/
                    for (int j = row_start[index_2]; j < row_start[index_2+1]; ++j){
                        //cout<<mat_val[j]<<" ";
                        if (col_idx[j] == index_2){
                            mat_val[j] = 1.0;   /** 对角线元素 **/
                        }
                        else {
                            mat_val[j] = 0.0;   /**< 非对角线元素 */
                        }
                    }
                }
            }else if (EBoundassignment[i_bc].bc_type == "gate-simplified"){
                /** load fixed volt value **/
                double dins = EBoundassignment[i_bc].sc_coe[4]*m;
                double epsilon_ins = EBoundassignment[i_bc].sc_coe[5];
                double Vapp = bound_source->source_assign(EBoundassignment[i_bc].sc_type,time, dt,EBoundassignment[i_bc].sc_coe);
                double WORKFUNC = EBoundassignment[i_bc].sc_coe[3];
                double local_Qt = EBoundassignment[i_bc].sc_coe[6];
                double Qins = local_Qt*pow(m,-2);
                if(int(time/dt)==0){
                    Qins = 0;
                    //Vapp = 0;
                }
                double value_rhs = (e*Qins/eps0 + epsilon_ins * (Vapp - WORKFUNC)/dins);
                //double Qins = Qt*pow(m,-2);
                tbox::Pointer<BaseElement<NDIM> > ele_quadrangle = ele_manager->getElement("LinearQuadrangle");
                int num_faces = patch.getNumberOfEntities(hier::EntityUtilities::FACE, 1);
                int num_face = EBoundassignment[i_bc].bc_face.getSize();
                for (int i_f = 0; i_f < num_face; ++i_f){
                    if(patch_geo->hasEntitySet(EBoundassignment[i_bc].bc_face[i_f], hier::EntityUtilities::FACE)){
                        //cout<<"get the gate-simplified face:"<<EBoundassignment[i_bc].bc_face[i_f]<<endl;
                        const tbox::Array<int> & face_id = patch_geo->getEntityIndicesInSet(EBoundassignment[i_bc].bc_face[i_f],
                                                           hier::EntityUtilities::FACE,
                                                           num_faces);
                        int size_f = face_id.size();
                        for (int i_e = 0; i_e < size_f; ++i_e){
                            int n_fvertex = face_node_ext[face_id[i_e] + 1] - face_node_ext[face_id[i_e]];
                            tbox::Array<hier::DoubleVector<NDIM> > face_vertex(n_fvertex);
                            tbox::Array<int> node_n(n_fvertex);
                            tbox::Array<int> mapping(n_fvertex);
                            tbox::Array<int> mapping_semi(n_fvertex);
                            tbox::Array<double> value(n_fvertex);
                            tbox::Array<double> V_node(n_fvertex);

                            /** inform for node **/
                            for (int i1 = 0, j = face_node_ext[face_id[i_e]]; i1 < n_fvertex; ++i1, ++j) {
                                node_n[i1] = face_node_idx[j];
                                mapping_semi[i1] = dof_map_semi[face_node_idx[j]];
                                V_node[i1] = (*V_soln)(mapping_semi[i1]);
                                for(int k = 0; k< NDIM; ++k){
                                    face_vertex[i1][k] = (*node_coord)(k, face_node_idx[j])*m;
                                }
                            }

                            tbox::Pointer<tbox::Vector<double> > ele_vec = new tbox::Vector<double>();
                            ele_vec->resize(n_fvertex);

                            for (int i_n = 0; i_n < n_fvertex; ++i_n){
                                value[i_n] = -epsilon_ins/dins;
                            }
                            /** building element vector **/
                            ele_quadrangle->buildCVFEMRhs(face_vertex, value, dt, time, ele_vec);

                            /** add up the element matrix for system matrix **/
                            int index = 0;
                            for (int l = 0; l < n_fvertex; ++l){
                                index = mapping_semi[l];
                                (*mat_data)(index,index) = (*mat_data)(index,index) + (*ele_vec)[l];
                            }

                            /** element vector RHS **/

                            for (int i_n = 0; i_n < n_fvertex; ++i_n){
                                value[i_n] = - value_rhs + epsilon_ins*V_node[i_n]/dins;
                            }
                            /** building element vector **/
                            ele_quadrangle->buildCVFEMRhs(face_vertex, value, dt, time, ele_vec);

                            /** construct fval vector **/
                            for(int i2= 0; i2 < n_fvertex; ++i2)
                                fval_val[mapping_semi[i2]] += (*ele_vec)[i2];
                        }
                    }else{
                        //cout<<"can't get the gate-simplified face:"<<EBoundassignment[i_bc].bc_face[i_f]<<endl;
                    }
                }
            }else if (EBoundassignment[i_bc].bc_type == "diode"){
                double voltage = 0;
                int index_semi = dof_map_semi[entity_idx[0]];
                voltage = (*V_soln)(index_semi);
                stringstream file;
                file << "diode_voltage";
                string file_dir = file.str();
                // open the file stream
                ofstream streamout;
                streamout.open(file_dir);
                //output the electric potential of diode
                streamout<<voltage;
            }else if (EBoundassignment[i_bc].bc_type == "gate-diode"){
                /** load fixed volt value **/
                //double Vapp = bound_source->source_assign(EBoundassignment[i_bc].sc_type,time, dt, EBoundassignment[i_bc].sc_coe);
                //double WORKFUNC = EBoundassignment[i_bc].sc_coe[3];
                double dins = EBoundassignment[i_bc].sc_coe[4]*m;
                double epsilon_ins = EBoundassignment[i_bc].sc_coe[5];
                //double Qins = Qt*pow(m,-2);
                tbox::Pointer<BaseElement<NDIM> > ele_quadrangle = ele_manager->getElement("LinearQuadrangle");
                int num_faces = patch.getNumberOfEntities(hier::EntityUtilities::FACE, 1);
                int num_face = EBoundassignment[i_bc].bc_face.getSize();
                for (int i_f = 0; i_f < num_face; ++i_f){
                    if(patch_geo->hasEntitySet(EBoundassignment[i_bc].bc_face[i_f], hier::EntityUtilities::FACE)){
                        //cout<<"get the gate-simplified face:"<<EBoundassignment[i_bc].bc_face[i_f]<<endl;
                        const tbox::Array<int> & face_id = patch_geo->getEntityIndicesInSet(EBoundassignment[i_bc].bc_face[i_f],
                                                           hier::EntityUtilities::FACE,
                                                           num_faces);
                        int size_f = face_id.size();
                        for (int i_e = 0; i_e < size_f; ++i_e){
                            int n_fvertex = face_node_ext[face_id[i_e] + 1] - face_node_ext[face_id[i_e]];
                            tbox::Array<hier::DoubleVector<NDIM> > face_vertex(n_fvertex);
                            tbox::Array<int> node_n(n_fvertex);
                            tbox::Array<int> mapping(n_fvertex);
                            tbox::Array<int> mapping_semi(n_fvertex);
                            tbox::Array<double> value(n_fvertex);

                            /** inform for node **/
                            for (int i1 = 0, j = face_node_ext[face_id[i_e]]; i1 < n_fvertex; ++i1, ++j) {
                                node_n[i1] = face_node_idx[j];
                                mapping_semi[i1] = dof_map_semi[face_node_idx[j]];
                                for(int k = 0; k< NDIM; ++k){
                                    face_vertex[i1][k] = (*node_coord)(k, face_node_idx[j])*m;
                                }
                            }

                            tbox::Pointer<tbox::Vector<double> > ele_vec = new tbox::Vector<double>();
                            ele_vec->resize(n_fvertex);

                            for (int i_n = 0; i_n < n_fvertex; ++i_n){
                                value[i_n] = -epsilon_ins/dins;
                            }
                            /** building element vector **/
                            ele_quadrangle->buildCVFEMRhs(face_vertex, value, dt, time, ele_vec);

                            /** add up the element matrix for system matrix **/
                            int index = 0;
                            for (int l = 0; l < n_fvertex; ++l){
                                index = mapping_semi[l];
                                (*mat_data)(index,index) = (*mat_data)(index,index) + (*ele_vec)[l];
                            }
                        }
                    }else{
                        cout<<"can't get the gate-diode face:"<<EBoundassignment[i_bc].bc_face[i_f]<<endl;
                    }
                }
            }
        }
    }
}

/**
 * @brief NonlinearPoisson::transferPoissonSoln: transfer from current to plot
 * @param patch
 * @param time
 * @param dt
 * @param component_name
 */
void NonlinearPoisson::transferPoissonSoln(hier::Patch<NDIM> & patch,
                                    const double  time,
                                    const double  dt,
                                    const string& component_name)
{
  tbox::Pointer<pdat::NodeData<NDIM, double> >
    plot_data = patch.getPatchData(d_poisson_plot_id);
  tbox::Pointer<pdat::NodeData<NDIM, double> >
    sol_1_data = patch.getPatchData(d_phi_sol_1_id);
  tbox::Pointer<pdat::VectorData<NDIM, double> >
    solution = patch.getPatchData(d_poisson_sol_cur_id);
  tbox::Pointer<pdat::NodeData<NDIM, int> >
    semiinsbc_node = patch.getPatchData(d_semiinsbc_node_id);
  tbox::Pointer<pdat::NodeData<NDIM, int> >
        siclass_node = patch.getPatchData(d_siclass_node_id);

  int num_nodes = patch.getNumberOfNodes(1);
  int *dof_map_semi = d_dof_semi_info->getDOFMapping(patch, hier::EntityUtilities::NODE);

  for(int i = 0; i < num_nodes; ++i){
      if((*siclass_node)(0,i)==0) continue;
      int index_s = dof_map_semi[i];
      plot_data->getPointer()[i] = solution->getPointer()[index_s];
      (*sol_1_data)(0,i) = solution->getPointer()[index_s];
      if((*semiinsbc_node)(0,i)==1){
          index_s = dof_map_semi[i]+3;
          (*sol_1_data)(1,i) = solution->getPointer()[index_s];
      }
  }
}

/**
 * @brief NonlinearPoisson::transferCarnSoln: transfer current solution to plot
 * @param patch
 * @param time
 * @param dt
 * @param component_name
 */
void NonlinearPoisson::transferCarnSoln(hier::Patch<NDIM> & patch,
                                      const double  time,
                                      const double  dt,
                                      const string& component_name)
{
    tbox::Pointer<pdat::VectorData<NDIM, double> >
      solution = patch.getPatchData(d_poisson_sol_cur_id);
    tbox::Pointer<pdat::NodeData<NDIM, double> >
      plotn_data = patch.getPatchData(d_carn_plot_id);
    tbox::Pointer<pdat::NodeData<NDIM, double> >
      logn_plot_data = patch.getPatchData(d_logn_plot_id);
    tbox::Pointer<pdat::NodeData<NDIM, double> >
      fermin_plot_data = patch.getPatchData(d_fermin_plot_id);

    tbox::Pointer<pdat::NodeData<NDIM, double> >
      carn_sol_1_data = patch.getPatchData(d_carn_sol_1_id);
    tbox::Pointer<pdat::NodeData<NDIM, double> >
      fn_sol_1_data = patch.getPatchData(d_fn_sol_1_id);

    tbox::Pointer<pdat::VectorData<NDIM, double> >
      solution_n = patch.getPatchData(d_carn_sol_cur_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
      logn_solution = patch.getPatchData(d_logn_sol_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
      fermin_solution = patch.getPatchData(d_fermin_sol_id);

    tbox::Pointer<pdat::NodeData<NDIM, double> >
      ploth_data = patch.getPatchData(d_carh_plot_id);
    tbox::Pointer<pdat::NodeData<NDIM, double> >
      logh_plot_data = patch.getPatchData(d_logh_plot_id);
    tbox::Pointer<pdat::NodeData<NDIM, double> >
      fermih_plot_data = patch.getPatchData(d_fermih_plot_id);

    tbox::Pointer<pdat::NodeData<NDIM, double> >
      carh_sol_1_data = patch.getPatchData(d_carh_sol_1_id);
    tbox::Pointer<pdat::NodeData<NDIM, double> >
      fh_sol_1_data = patch.getPatchData(d_fh_sol_1_id);

    tbox::Pointer<pdat::VectorData<NDIM, double> >
      solution_h = patch.getPatchData(d_carh_sol_cur_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
      logh_solution = patch.getPatchData(d_logh_sol_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
      fermih_solution = patch.getPatchData(d_fermih_sol_id);

    tbox::Pointer<pdat::NodeData<NDIM, int> >
      semi_node = patch.getPatchData(d_semi_node_id);
    tbox::Pointer<pdat::NodeData<NDIM, int> >
      sio2_node = patch.getPatchData(d_sio2_node_id);
    tbox::Pointer<pdat::NodeData<NDIM, int> >
      semiinsbc_node = patch.getPatchData(d_semiinsbc_node_id);

    tbox::Pointer<pdat::NodeData<NDIM, int> >
          siclass_node = patch.getPatchData(d_siclass_node_id);

    int *dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);
    int *dof_map_semi = d_dof_semi_info->getDOFMapping(patch, hier::EntityUtilities::NODE);

    /** 获取本地单元(边, 面)周围结点的索引关系 **/
    int num_nodes = patch.getNumberOfNodes(1);

    for(int i = 0; i < num_nodes; ++i){
        if((*siclass_node)(0,i)==0) continue;
        int index = dof_map[i];
        int index_s = dof_map_semi[i];
        plotn_data->getPointer()[i] = solution->getPointer()[index_s+1]/pow(m,-3);//solution_n->getPointer()[index]/pow(m,-3);
        logn_plot_data->getPointer()[i] = log10(abs(plotn_data->getPointer()[i]));//logn_solution->getPointer()[index];
        fermin_plot_data->getPointer()[i] = fermin_solution->getPointer()[index];
        (*carn_sol_1_data)(0,i) = solution->getPointer()[index_s+1];
        (*fn_sol_1_data)(0,i) =fermin_solution->getPointer()[index];

        ploth_data->getPointer()[i] = solution->getPointer()[index_s+2]/pow(m,-3);//solution_h->getPointer()[index]/pow(m,-3);
        logh_plot_data->getPointer()[i] = log10(abs(ploth_data->getPointer()[i]));//logh_solution->getPointer()[index];
        fermih_plot_data->getPointer()[i] = fermih_solution->getPointer()[index];
        (*carh_sol_1_data)(0,i) = solution->getPointer()[index_s+2];
        (*fh_sol_1_data)(0,i) = (*fermih_solution)(index);

        if((*semiinsbc_node)(0,i)==1){
            int index = dof_map[i]+1;
            (*carn_sol_1_data)(1,i) = solution->getPointer()[index_s+4];
            (*fn_sol_1_data)(1,i) =fermin_solution->getPointer()[index];

            (*carh_sol_1_data)(1,i) = solution->getPointer()[index_s+5];
            (*fh_sol_1_data)(1,i) = (*fermih_solution)(index);
        }
    }

#if 1
    /** 取出本地PatchGeometry **/
    tbox::Pointer< hier::PatchGeometry<NDIM> > patch_geo = patch.getPatchGeometry();
    /** get PatchTopology **/
    tbox::Pointer< hier::PatchTopology<NDIM> > patch_top = patch.getPatchTopology();
    tbox::Pointer<pdat::NodeData<NDIM, double> >
      poi_sol_1_data = patch.getPatchData(d_phi_sol_1_id);
    tbox::Pointer<pdat::NodeData<NDIM, double> >
            sol_data = patch.getPatchData(d_sol_id);

    int num = patch.getNumberOfNodes(0);
    if(save_carrier){
        ofstream streamout;
        streamout.open("carrier",ios::out);

        /** 取出本地Patch的结点坐标数组 **/
        tbox::Pointer<pdat::NodeData<NDIM, double> > node_coord = patch_geo->getNodeCoordinates();
        for(int i = 0; i < num; ++i){
            if((*siclass_node)(0,i)==0) continue;
            streamout<<(*node_coord)(0,i)<<" "<<(*node_coord)(1,i)<<" "<<(*node_coord)(2,i)<<" ";
            streamout<<(*carn_sol_1_data)(0,i)<<" "<<(*carh_sol_1_data)(0,i)<<" "<<(*poi_sol_1_data)(0,i)<<"\n";
            if((*semiinsbc_node)(0,i)==1){
                streamout<<(*node_coord)(0,i)<<" "<<(*node_coord)(1,i)<<" "<<(*node_coord)(2,i)<<" ";
                streamout<<(*carn_sol_1_data)(1,i)<<" "<<(*carh_sol_1_data)(1,i)<<" "<<(*poi_sol_1_data)(1,i)<<"\n";
            }
        }
        streamout.close();
    }

    if(save_sol){
        stringstream sol_file;
        sol_file << "sol_" << int(time/dt)<<"_"<<patch.getIndex();//patch.getPatchId()<<
        ofstream streamout_sol;
        streamout_sol.open(sol_file.str(),ios::out);
        if(!streamout_sol) cout<<"can't open file "<<sol_file.str()<<endl;
        streamout_sol<<"nd_id volt n p\n";
        for(int i = 0; i < num; ++i){
            if((*siclass_node)(0,i)==0) continue;
            int index_s = dof_map_semi[i];
            streamout_sol<<(*sol_data)(3,i)<<" "<<solution->getPointer()[index_s]
                       <<" "<<solution->getPointer()[index_s+1]/pow(m,-3)
                      <<" "<<solution->getPointer()[index_s+2]/pow(m,-3)<<"\n";
            if((*semiinsbc_node)(0,i)==1){
                index_s = dof_map_semi[i]+3;
                streamout_sol<<(*sol_data)(3,i)<<" "<<solution->getPointer()[index_s]
                           <<" "<<solution->getPointer()[index_s+1]/pow(m,-3)
                          <<" "<<solution->getPointer()[index_s+2]/pow(m,-3)<<"\n";
            }
        }
        streamout_sol.close();
    }
#endif
}

/**
* @brief NonlinearPoisson::calculateQuasiFermilevel: calculate quasi fermi level
* @param patch
* @param time
* @param dt
* @param component_name
*/
void NonlinearPoisson::calculateQuasiFermilevel(hier::Patch<NDIM>& patch,
                                        const double time,
                                        const double dt,
                                        const string& component_name){
   /** 取出本地PatchTopology **/
   tbox::Pointer< hier::PatchTopology<NDIM> > patch_top = patch.getPatchTopology();
   tbox::Pointer<pdat::NodeData<NDIM, int> >
     semi_node = patch.getPatchData(d_semi_node_id);
   tbox::Pointer<pdat::NodeData<NDIM, int> >
     sio2_node = patch.getPatchData(d_sio2_node_id);
   tbox::Pointer<pdat::NodeData<NDIM, int> >
     semiinsbc_node = patch.getPatchData(d_semiinsbc_node_id);

   /** get node values **/
   tbox::Pointer<pdat::VectorData<NDIM, double> >
     V_soln = patch.getPatchData(d_poisson_sol_cur_id);
   tbox::Pointer<pdat::VectorData<NDIM, double> >
     N_soln = patch.getPatchData(d_carn_sol_cur_id);
   tbox::Pointer<pdat::VectorData<NDIM, double> >
     H_soln = patch.getPatchData(d_carh_sol_cur_id);

   tbox::Pointer<pdat::NodeData<NDIM, int> >
     heter_node = patch.getPatchData(d_heter_node_id);
   tbox::Pointer<pdat::NodeData<NDIM, int> >
     hnode_mat = patch.getPatchData(d_hnode_mat_id);

   tbox::Pointer< pdat::NodeData<NDIM, int> >
           node_mat = patch.getPatchData(d_node_mat_id);

   /** quasi fermi level **/
   tbox::Pointer<pdat::VectorData<NDIM,double> >
     Fn_soln = patch.getPatchData(d_fermin_sol_id);
   tbox::Pointer<pdat::VectorData<NDIM, double> >
     Fh_soln = patch.getPatchData(d_fermih_sol_id);

   tbox::Pointer<pdat::NodeData<NDIM, int> >
         siclass_node = patch.getPatchData(d_siclass_node_id);

   int * dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);
   int * dof_map_semi = d_dof_semi_info->getDOFMapping(patch, hier::EntityUtilities::NODE);

   /** 获取本地单元(边, 面)周围结点的索引关系 **/
   tbox::Array<int> can_extent, can_indices;
   patch_top->getCellAdjacencyNodes(can_extent, can_indices);
   int num_nodes = patch.getNumberOfNodes(1);

   /** material list **/
   tbox::Pointer<MaterialManager<NDIM> >
           mat_manager = MaterialManager<NDIM>::getManager();

   /** setup target for class 'Variables_Functions' **/
   tbox::Pointer<Variables_Functions> Var_func = new Variables_Functions();

   /** traversing nodes **/
   double tmp_fermin, tmp_fermih;
   for(int i = 0; i < num_nodes; ++i){
       if((*siclass_node)(0,i)==0) continue;
       if((*semi_node)(0,i) == 1){
           int index = dof_map[i];
           int index_s = dof_map_semi[i];

           (*N_soln)(index) = (*V_soln)(index_s+1);
           (*H_soln)(index) = (*V_soln)(index_s+2);

           //cout<<index<<" "<<(*N_soln)(index)<<" "<<(*H_soln)(index)<<endl;

           tbox::Pointer<BaseMaterial<NDIM> > mat_node =
                   mat_manager->getMaterial("matSilicon");

           double carn_tmp = fabs((*N_soln)(index));
           double carh_tmp = fabs((*H_soln)(index));
           double v_tmp = (*V_soln)(index_s);
           double T_tmp = T;

           double affin_tmp = mat_node->effectiveAffinity(T_tmp);
           double nc_tmp = mat_node->effectiveDensityNc(T_amb)*pow(m,-3);;
           double nv_tmp = mat_node->effectiveDensityNv(T_amb)*pow(m,-3);;
           double eg_tmp = mat_node->bandGap(T_tmp);

           Var_func->antiDummyElectron(carn_tmp, v_tmp, affin_tmp, T_tmp, kb, e, nc_tmp, tmp_fermin);
           Var_func->antiDummyHole(carh_tmp, v_tmp, eg_tmp, affin_tmp, T_tmp, kb, e, nv_tmp, tmp_fermih);

           (*Fn_soln)(index) = (3.0/3)*tmp_fermin+(0.0/3)*(*Fn_soln)(index);
           (*Fh_soln)(index) = (3.0/3)*tmp_fermih+(0.0/3)*(*Fh_soln)(index);

           #if 0
           if(abs(tmp_fermih)>=1){
               cout<<(*N_soln)(index)<<" "<<(*H_soln)(index)<<endl;
           }
           #endif

           if((*semiinsbc_node)(0,i) == 1){
               int index = dof_map[i]+1;
               int index_s = dof_map_semi[i]+3;

               (*N_soln)(index) = (*V_soln)(index_s+1);
               (*H_soln)(index) = (*V_soln)(index_s+2);

               //cout<<index<<" "<<(*N_soln)(index)<<" "<<(*H_soln)(index)<<endl;

               tbox::Pointer<BaseMaterial<NDIM> > mat_node =
                       mat_manager->getMaterial("matSiO2");

               double carn_tmp = 1e-10*pow(m,-3);
               double carh_tmp = 1e-10*pow(m,-3);
               double v_tmp = (*V_soln)(index_s);
               double T_tmp = T;

               double affin_tmp = mat_node->effectiveAffinity(T_tmp);
               double nc_tmp = mat_node->effectiveDensityNc(T_amb)*pow(m,-3);;
               double nv_tmp = mat_node->effectiveDensityNv(T_amb)*pow(m,-3);;
               double eg_tmp = mat_node->bandGap(T_tmp);

               Var_func->antiDummyElectron(carn_tmp, v_tmp, affin_tmp, T_tmp, kb, e, nc_tmp, tmp_fermin);
               Var_func->antiDummyHole(carh_tmp, v_tmp, eg_tmp, affin_tmp, T_tmp, kb, e, nv_tmp, tmp_fermih);

               (*Fn_soln)(index) = (3.0/3)*tmp_fermin+(0.0/3)*(*Fn_soln)(index);
               (*Fh_soln)(index) = (3.0/3)*tmp_fermih+(0.0/3)*(*Fh_soln)(index);
           }
       }else if((((*sio2_node)(0,i) == 1)&&((*semiinsbc_node)(0,i) == 0))){
           int index = dof_map[i];
           int index_s = dof_map_semi[i];

           (*N_soln)(index) = (*V_soln)(index_s+1);
           (*H_soln)(index) = (*V_soln)(index_s+2);

           //cout<<index<<" "<<(*N_soln)(index)<<" "<<(*H_soln)(index)<<endl;

           tbox::Pointer<BaseMaterial<NDIM> > mat_node =
                   mat_manager->getMaterial("matSiO2");

           double carn_tmp = 1e-10*pow(m,-3);
           double carh_tmp = 1e-10*pow(m,-3);
           double v_tmp = (*V_soln)(index_s);
           double T_tmp = T;

           double affin_tmp = mat_node->effectiveAffinity(T_tmp);
           double nc_tmp = mat_node->effectiveDensityNc(T_amb)*pow(m,-3);;
           double nv_tmp = mat_node->effectiveDensityNv(T_amb)*pow(m,-3);;
           double eg_tmp = mat_node->bandGap(T_tmp);

           Var_func->antiDummyElectron(carn_tmp, v_tmp, affin_tmp, T_tmp, kb, e, nc_tmp, tmp_fermin);
           Var_func->antiDummyHole(carh_tmp, v_tmp, eg_tmp, affin_tmp, T_tmp, kb, e, nv_tmp, tmp_fermih);

           (*Fn_soln)(index) = (3.0/3)*tmp_fermin+(0.0/3)*(*Fn_soln)(index);
           (*Fh_soln)(index) = (3.0/3)*tmp_fermih+(0.0/3)*(*Fh_soln)(index);
       }
   }
}

/**
 * @brief NonlinearPoisson::initialQuasiFermilevel: initialize quasi fermi level
 * @param patch
 * @param time
 * @param initial_time
 * @param component_name
 */
void NonlinearPoisson::initialQuasiFermilevel(hier::Patch<NDIM> & patch,
                                       const double  time,
                                       const double    initial_time,
                                       const string& component_name){
    /** 取出本地PatchTopology **/
    tbox::Pointer< hier::PatchTopology<NDIM> > patch_top = patch.getPatchTopology();

    /** quasi fermi level **/
    tbox::Pointer<pdat::VectorData<NDIM,double> >
      Fn_soln = patch.getPatchData(d_fermin_sol_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
      Fh_soln = patch.getPatchData(d_fermih_sol_id);

    tbox::Pointer<pdat::NodeData<NDIM, int> >
      heter_node = patch.getPatchData(d_heter_node_id);
    tbox::Pointer<pdat::NodeData<NDIM, int> >
      semiinsbc_node = patch.getPatchData(d_semiinsbc_node_id);
    tbox::Pointer<pdat::NodeData<NDIM, int> >
      sio2_node = patch.getPatchData(d_sio2_node_id);

    tbox::Pointer<pdat::NodeData<NDIM, int> >
          siclass_node = patch.getPatchData(d_siclass_node_id);

    int * dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);

    /** 获取本地单元(边, 面)周围结点的索引关系. **/
    tbox::Array<int> can_extent, can_indices;
    patch_top->getCellAdjacencyNodes(can_extent, can_indices);

    int num_nodes = patch.getNumberOfNodes(1);

    /** traversing nodes **/
    for(int i = 0; i < num_nodes; ++i){
        if((*siclass_node)(0,i)==0) continue;
        int index_s = dof_map[i];
        (*Fn_soln)(index_s) = 0.0;
        (*Fh_soln)(index_s) = 0.0;
        if(((*sio2_node)(0,i)==1)&&((*semiinsbc_node)(0,i)==0)){
            (*Fn_soln)(index_s) = -5.0;
            (*Fh_soln)(index_s) = -5.0;
        }
        if((*semiinsbc_node)(0,i)==1){
            int index_s = dof_map[i]+1;
            (*Fn_soln)(index_s) = -5.0;
            (*Fh_soln)(index_s) = -5.0;
        }
    }
}

/*************************************************************************
 * calculate the electric field
 ************************************************************************/
void NonlinearPoisson::calculateEx(hier::Patch<NDIM>& patch,
                                               const double time,
                                               const double dt,
                                               const string& component_name) {

  /// (形函数, 积分器, 单元)管理器
  tbox::Pointer<ElementManager<NDIM> > ele_manager =
      ElementManager<NDIM>::getManager();
  tbox::Pointer<BaseElement<NDIM> > ele =
      ele_manager->getElement(d_element_type);

  /// 取出本地PatchGeometry.
  tbox::Pointer<hier::PatchGeometry<NDIM> > patch_geo =
      patch.getPatchGeometry();

  /// 取出本地PatchTopology.
  tbox::Pointer<hier::PatchTopology<NDIM> > patch_top =
      patch.getPatchTopology();
  /// 取出本地Patch的结点坐标数组.
  tbox::Pointer<pdat::NodeData<NDIM, double> > node_coord =
      patch_geo->getNodeCoordinates();
  //取得数据片
  tbox::Pointer<pdat::CellData<NDIM, double> > Ex_data =
      patch.getPatchData(d_E_cell_id);
  tbox::Pointer<pdat::VectorData<NDIM, double> >
      Phi_data = patch.getPatchData(d_poisson_sol_cur_id);

  tbox::Pointer<pdat::CellData<NDIM, int> >
    cell_mat = patch.getPatchData(d_cell_mat_id);

  /// 获取本地单元(边, 面)周围结点的索引关系.
  tbox::Array<int> can_extent, can_indices;
  patch_top->getCellAdjacencyNodes(can_extent, can_indices);
  int num_cells = patch.getNumberOfCells(1);//1表示影像区的宽度
  //cout << num_cells << endl;
  /// 取出自由度映射信息
  int* dof_map_semi = d_dof_semi_info->getDOFMapping(patch, hier::EntityUtilities::NODE);

  /// 遍历单元

  for (int i = 0; i < num_cells; ++i) {
    if(!((material_list[(*cell_mat)(0, i)].mat_type == "semiconductor")||(material_list[(*cell_mat)(0, i)].mat_type == "dielectric"))) continue;
    int cell=i;
    int n_vertex = can_extent[i + 1] - can_extent[i];
    tbox::Array<hier::DoubleVector<NDIM> > vertex(n_vertex);
    tbox::Array<double> phi(n_vertex);//节点电势
    tbox::Array<int> mapping_semi(n_vertex);

    /// 取出单元结点坐标，以及填写映射值。
    for (int i1 = 0, j = can_extent[i]; i1 < n_vertex; ++i1, ++j) {
      mapping_semi[i1] = dof_map_semi[can_indices[j]];
      for (int k = 0; k < NDIM; ++k) {
        vertex[i1][k] = (*node_coord)(k, can_indices[j])*m;
      }
      phi[i1]=(*Phi_data)(mapping_semi[i1]);//取得单元八个节点的电势值
    }

    tbox::Array<tbox::Array<int> > r_node = ele->get_r_node();
    tbox::Array<double> edge_E(cell_edge);//边电场
    for(int i_e=0;i_e<cell_edge;i_e++){
        int n_1 = r_node[i_e][0];
        int n_2 = r_node[i_e][1];
        double Vi=phi[n_1];
        double Vj=phi[n_2];
        double Lij=sqrt(pow((vertex[n_1][0]-vertex[n_2][0]),2)+pow((vertex[n_1][1]-vertex[n_2][1]),2)+pow((vertex[n_1][2]-vertex[n_2][2]),2));
        edge_E[i_e]=(Vj-Vi)/Lij;

    }

    /// 计算每个单元的Ex
    tbox::Pointer<tbox::Array<double> > ele_Ex = new tbox::Array<double >(3);
    for (int i = 0; i < 3; ++i) {
      (*ele_Ex)[i] = 0.0;
    }

    //ele->calculateElemE(vertex,phi,dt,time,ele_Ex);
    ele->calculateValueByEdge(vertex,edge_E,dt,time,ele_Ex);
    //赋值给相应的数据片单元
    (*Ex_data)(0,cell)=(*ele_Ex)[0]/(V/m);
    (*Ex_data)(1,cell)=(*ele_Ex)[1]/(V/m);
    (*Ex_data)(2,cell)=(*ele_Ex)[2]/(V/m);
    }
}

/**
 * @brief NonlinearPoisson::calculateElemJ: current density in each cell
 * @param patch
 * @param time
 * @param dt
 * @param component_name
 */
void NonlinearPoisson::calculateElemJ(hier::Patch<NDIM> &patch,
                                      const double time,
                                      const double dt,
                                      const string & component_name){
    /** (形函数, 积分器, 单元)管理器 **/
    tbox::Pointer<ElementManager<NDIM> >
      ele_manager = ElementManager<NDIM>::getManager();
    tbox::Pointer<BaseElement<NDIM> >
      ele = ele_manager->getElement(d_element_type);

    /** 取出本地PatchGeometry **/
    tbox::Pointer< hier::PatchGeometry<NDIM> > patch_geo = patch.getPatchGeometry();
    /** 取出本地PatchTopology **/
    tbox::Pointer< hier::PatchTopology<NDIM> > patch_top = patch.getPatchTopology();
    /** 取出本地Patch的结点坐标数组 **/
    tbox::Pointer<pdat::NodeData<NDIM, double> >
      node_coord = patch_geo->getNodeCoordinates();

    tbox::Pointer<pdat::NodeData<NDIM, int> >
            heter_node = patch.getPatchData(d_heter_node_id);
    tbox::Pointer<pdat::CellData<NDIM, int> >
            cell_hnode = patch.getPatchData(d_cell_hnode_id);

    /** get node values **/
    tbox::Pointer<pdat::VectorData<NDIM, double> >
      V_soln = patch.getPatchData(d_poisson_sol_cur_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
      N_soln = patch.getPatchData(d_carn_sol_cur_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
      H_soln = patch.getPatchData(d_carh_sol_cur_id);
    tbox::Pointer<pdat::VectorData<NDIM,double> >
      Fn_soln = patch.getPatchData(d_fermin_sol_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
      Fh_soln = patch.getPatchData(d_fermih_sol_id);
    tbox::Pointer<pdat::CellData<NDIM, int> >
      cell_mat = patch.getPatchData(d_cell_mat_id);

    /** current density vector **/
    tbox::Pointer<pdat::CellData<NDIM,double> >
      Jn_data = patch.getPatchData(d_Jn_vec_id);
    tbox::Pointer<pdat::CellData<NDIM,double> >
      Jp_data = patch.getPatchData(d_Jp_vec_id);
    tbox::Pointer<pdat::CellData<NDIM,double> >
      J_data = patch.getPatchData(d_J_vec_id);

    /** 获取本地单元(边, 面)周围结点的索引关系 **/
    tbox::Array<int> can_extent, can_indices;
    patch_top->getCellAdjacencyNodes(can_extent, can_indices);

    int num_cells = patch.getNumberOfCells(1);

    /** 取出自由度映射信息 **/
    int * dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);
    int * dof_map_semi = d_dof_semi_info->getDOFMapping(patch, hier::EntityUtilities::NODE);

    tbox::Pointer<MaterialManager<NDIM> >
            mat_manager = MaterialManager<NDIM>::getManager();

    /** setup target for class 'Variables_Functions' **/
    tbox::Pointer<Variables_Functions> Var_func = new Variables_Functions();

    tbox::Pointer<tbox::Array<double> > elem_E = new tbox::Array<double>(NDIM);
    tbox::Pointer<tbox::Array<double> > grad_N = new tbox::Array<double>(NDIM);
    tbox::Pointer<tbox::Array<double> > grad_H = new tbox::Array<double>(NDIM);
    tbox::Pointer<tbox::Array<double> > grad_Fn = new tbox::Array<double>(NDIM);
    tbox::Pointer<tbox::Array<double> > grad_Fh = new tbox::Array<double>(NDIM);

    /** 遍历单元 **/
    for(int i = 0; i < num_cells; ++i){
        int cell=i;
        if(!((material_list[(*cell_mat)(0, i)].mat_type == "semiconductor")||(material_list[(*cell_mat)(0, i)].mat_type == "dielectric"))) continue;

        int n_vertex = can_extent[i+1] - can_extent[i];
        tbox::Array<hier::DoubleVector<NDIM> > vertex(n_vertex);
        tbox::Array<int> mapping(n_vertex);
        tbox::Array<int> mapping_semi(n_vertex);
        tbox::Array<int> node_n(n_vertex);
        tbox::Array<double> T_node(n_vertex);
        tbox::Array<double> V_node(n_vertex);
        tbox::Array<double> Fn_node(n_vertex);
        tbox::Array<double> Fh_node(n_vertex);
        tbox::Array<double> N_node(n_vertex);
        tbox::Array<double> H_node(n_vertex);

        /** setup material target **/
        tbox::Pointer<BaseMaterial<NDIM> > mat_cell =
                mat_manager->getMaterial(material_list[(*cell_mat)(0,i)].mat_name);

        for(int l=0; l<NDIM; ++l){
            (*elem_E)[l] = 0.0;
            (*grad_N)[l] = 0.0;
            (*grad_H)[l] = 0.0;
            (*grad_Fn)[l] = 0.0;
            (*grad_Fh)[l] = 0.0;
        }

        if(material_list[(*cell_mat)(0,i)].mat_type == "semiconductor"){

            /** get cell node coord **/
            for (int i1 = 0, j = can_extent[i]; i1 < n_vertex; ++i1, ++j) {
                mapping[i1] = dof_map[can_indices[j]];
                mapping_semi[i1] = dof_map_semi[can_indices[j]];
                for(int k = 0; k< NDIM; ++k){
                    vertex[i1][k] = (*node_coord)(k, can_indices[j])*m;
                }
                node_n[i1] = can_indices[j];
            }

            for (int i_n = 0; i_n < n_vertex; ++i_n){
                T_node[i_n]=T;
                V_node[i_n]=(*V_soln)(mapping_semi[i_n]);
                N_node[i_n]=(*V_soln)(mapping_semi[i_n]+1);
                H_node[i_n]=(*V_soln)(mapping_semi[i_n]+2);
                Fn_node[i_n]=(*Fn_soln)(mapping[i_n]);
                Fh_node[i_n]=(*Fh_soln)(mapping[i_n]);
            }
#if 0
            ele->calculateElemE(vertex,V_node,dt,time,elem_E);
            ele->calculateElemGradCar(vertex,N_node,dt,time,grad_N);
            ele->calculateElemGradCar(vertex,H_node,dt,time,grad_H);
            ele->calculateElemGradCar(vertex,Fn_node,dt,time,grad_Fn);
            ele->calculateElemGradCar(vertex,Fh_node,dt,time,grad_Fh);

            double T_cell = Var_func->node2ElemVariable(vertex, T_node, dt, time);
            double elem_n = Var_func->node2ElemVariable(vertex, N_node, dt, time);
            double elem_h = Var_func->node2ElemVariable(vertex, H_node, dt, time);

            double e_field = Var_func->vecNorm(elem_E);

            double mun = mat_cell->electronMobility(T_cell, e_field)*pow(m,2)/V/s;
            double mup = mat_cell->holeMobility(T_cell, e_field)*pow(m,2)/V/s;
            double dn = mat_cell->electronDiffusion(T_amb, e_field)*pow(m,2)/V/s;
            double dp = mat_cell->holeDiffusion(T_amb, e_field)*pow(m,2)/V/s;///tmp/VMwareDnD/YV9N4w/code.txt
            //NDIM:number of dimension
            for(int j=0; j<NDIM; ++j){
                /*
                (*Jn_data)(j, i) = e*mun*elem_n*(*elem_E)[j]+e*dn*(*grad_N)[j];
                (*Jp_data)(j, i) = e*mup*elem_h*(*elem_E)[j]+e*dp*(*grad_H)[j];
                (*J_data)(j, i) = (*Jn_data)(j, i) + (*Jp_data)(j, i);
                */
                (*Jn_data)(j, i) = (*elem_E)[j];
                (*Jp_data)(j, i) = (*elem_E)[j];
                (*J_data)(j, i) = (*Jn_data)(j, i) + (*Jp_data)(j, i);

                (*Jn_data)(j, i) = (*Jn_data)(j, i) ;
                (*Jp_data)(j, i) = (*Jp_data)(j, i) ;
                (*J_data)(j, i) = (*J_data)(j, i) ;
            }
#endif
#if 1



            tbox::Array<tbox::Array<int> > r_node = ele->get_r_node();
            tbox::Array<double> edge_Jn(cell_edge);//电子电流密度
            tbox::Array<double> edge_Jp(cell_edge);//空穴电流密度

            for(int i_e=0;i_e<cell_edge;i_e++){
                int n_1 = r_node[i_e][0];
                int n_2 = r_node[i_e][1];
                double Vi=V_node[n_1];
                double Vj=V_node[n_2];

                double N1 = N_node[n_1];
                double N2 = N_node[n_2];
                double H1 = H_node[n_1];
                double H2 = H_node[n_2];

                double Tmid=(T_node[n_1]+T_node[n_2])/2;
                double Vt=kb*Tmid/e;
                double mun = mat_cell->electronMobility()*pow(m,2)/V/s; //电子迁移率
                double mup = mat_cell->holeMobility()*pow(m,2)/V/s; //空穴迁移率
                double alpha_ij=(Vj-Vi)/Vt;
                double B1=Var_func->Bern(alpha_ij), B2=Var_func->Bern(-alpha_ij);
                double Lij=sqrt(pow((vertex[n_1][0]-vertex[n_2][0]),2)+pow((vertex[n_1][1]-vertex[n_2][1]),2)+pow((vertex[n_1][2]-vertex[n_2][2]),2));
                //edge_Jn[i_e]=e*mun*Vt/Lij*(N2*B1-N1*B2) +e*mun*(Vj-Vi)/Lij*(N2+N1)/2;
                //edge_Jp[i_e]=e*mup*Vt/Lij*(-H2*B2+H1*B1)+e*mup*(Vj-Vi)/Lij*(H2+H1)/2;
                edge_Jn[i_e]=(Vj-Vi)/Lij/2;
                edge_Jp[i_e]=(Vj-Vi)/Lij/2;
            }

            /// 计算每个单元的Ex
            tbox::Pointer<tbox::Array<double> > ele_Jn = new tbox::Array<double >(3);
            tbox::Pointer<tbox::Array<double> > ele_Jp = new tbox::Array<double >(3);
            for (int i = 0; i < 3; ++i) {
              (*ele_Jp)[i] = 0.0;
              (*ele_Jn)[i] = 0.0;
            }

            //ele->calculateElemE(vertex,phi,dt,time,ele_Ex);
            ele->calculateValueByEdge(vertex,edge_Jp,dt,time,ele_Jp);
            ele->calculateValueByEdge(vertex,edge_Jn,dt,time,ele_Jn);
            //赋值给相应的数据片单元
            (*Jp_data)(0,cell)=(*ele_Jp)[0]/(V/m);
            (*Jp_data)(1,cell)=(*ele_Jp)[1]/(V/m);
            (*Jp_data)(2,cell)=(*ele_Jp)[2]/(V/m);

            (*Jn_data)(0,cell)=(*ele_Jn)[0]/(V/m);
            (*Jn_data)(1,cell)=(*ele_Jn)[1]/(V/m);
            (*Jn_data)(2,cell)=(*ele_Jn)[2]/(V/m);

            (*J_data)(0,cell)=(*Jp_data)(0,cell)+(*Jn_data)(0,cell);
            (*J_data)(1,cell)=(*Jp_data)(1,cell)+(*Jn_data)(1,cell);
            (*J_data)(2,cell)=(*Jp_data)(2,cell)+(*Jn_data)(2,cell);

#endif
#if 0
            ele->calculateElemE(vertex,V_node,dt,time,elem_E);
            ele->calculateElemGradCar(vertex,N_node,dt,time,grad_N);
            ele->calculateElemGradCar(vertex,H_node,dt,time,grad_H);
            ele->calculateElemGradCar(vertex,Fn_node,dt,time,grad_Fn);
            ele->calculateElemGradCar(vertex,Fh_node,dt,time,grad_Fh);

            double T_cell = Var_func->node2ElemVariable(vertex, T_node, dt, time);
            double elem_n = Var_func->node2ElemVariable(vertex, N_node, dt, time);
            double elem_h = Var_func->node2ElemVariable(vertex, H_node, dt, time);

            double e_field = Var_func->vecNorm(elem_E);

            double mun = mat_cell->electronMobility(T_cell, e_field)*pow(m,2)/V/s;
            double mup = mat_cell->holeMobility(T_cell, e_field)*pow(m,2)/V/s;
            double dn = mat_cell->electronDiffusion(T_amb, e_field)*pow(m,2)/V/s;
            double dp = mat_cell->holeDiffusion(T_amb, e_field)*pow(m,2)/V/s;

            for(int j=0; j<NDIM; ++j){
                (*Jn_data)(j, i) = e*mun*elem_n*(*elem_E)[j] + e*dn*(*grad_N)[j];
                (*Jp_data)(j, i) = e*mup*elem_h*(*elem_E)[j] - e*dp*(*grad_H)[j];
                (*J_data)(j, i) = (*Jn_data)(j, i) + (*Jp_data)(j, i);

                (*Jn_data)(j, i) = (*Jn_data)(j, i) / (A/m/m);
                (*Jp_data)(j, i) = (*Jp_data)(j, i) / (A/m/m);
                (*J_data)(j, i) = (*J_data)(j, i) / (A/m/m);
            }
#endif
        }else{
            /** get cell node coord **/
            for (int i1 = 0, j = can_extent[i]; i1 < n_vertex; ++i1, ++j) {
                mapping_semi[i1] = dof_map_semi[can_indices[j]];
                for(int k = 0; k< NDIM; ++k){
                    vertex[i1][k] = (*node_coord)(k, can_indices[j])*m;
                }
            }

            for (int i_n = 0; i_n < n_vertex; ++i_n){
                T_node[i_n]=T;
                V_node[i_n]=(*V_soln)(mapping_semi[i_n]);
            }

            ele->calculateElemE(vertex,V_node,time,dt,elem_E);

            double T_cell = Var_func->node2ElemVariable(vertex, T_node, dt, time);

            double sigma = mat_cell->electroSigma(T_cell);

            for (int j = 0; j < NDIM; ++j){
                (*Jn_data)(j,i) = sigma*(*elem_E)[j];
                (*Jp_data)(j,i) = sigma*(*elem_E)[j];
                (*J_data)(j,i) = sigma*(*elem_E)[j];

                (*Jn_data)(j, i) = (*Jn_data)(j, i) / (A/m/m);
                (*Jp_data)(j, i) = (*Jp_data)(j, i) / (A/m/m);
                (*J_data)(j, i) = (*J_data)(j, i) / (A/m/m);
            }
        }
    }
}

/** ****************************************************************************** **/
void NonlinearPoisson::registerPlotData(tbox::Pointer<appu::JaVisDataWriter<NDIM> > javis_data_writer)
{
  javis_data_writer->registerPlotQuantity("plot_volt", "SCALAR", d_poisson_plot_id);
  javis_data_writer->registerPlotQuantity("plot_electron", "SCALAR", d_carn_plot_id);
  javis_data_writer->registerPlotQuantity("plot_hole", "SCALAR", d_carh_plot_id);
  javis_data_writer->registerPlotQuantity("plot_logn", "SCALAR", d_logn_plot_id);
  javis_data_writer->registerPlotQuantity("plot_logh", "SCALAR", d_logh_plot_id);
  javis_data_writer->registerPlotQuantity("plot_fermin", "SCALAR", d_fermin_plot_id);
  javis_data_writer->registerPlotQuantity("plot_fermih", "SCALAR", d_fermih_plot_id);
  javis_data_writer->registerPlotQuantity("carh_sol_1", "COMPOSITE", d_carh_sol_1_id);
  javis_data_writer->registerPlotQuantity("carn_sol_1", "COMPOSITE", d_carn_sol_1_id);
//  javis_data_writer->registerPlotQuantity("plot_Vogamma", "SCALAR", d_Vogamma_plot_id);
//  javis_data_writer->registerPlotQuantity("plot_Vogammaplus", "SCALAR", d_Vogammaplus_plot_id);
//  javis_data_writer->registerPlotQuantity("plot_VogammaHplus", "SCALAR", d_VogammaHplus_plot_id);
//  javis_data_writer->registerPlotQuantity("plot_VogammaH2plus", "SCALAR", d_VogammaH2plus_plot_id);
//  javis_data_writer->registerPlotQuantity("plot_VogammaH", "SCALAR", d_VogammaH_plot_id);
//  javis_data_writer->registerPlotQuantity("plot_VogammaH2", "SCALAR", d_VogammaH2_plot_id);
//  javis_data_writer->registerPlotQuantity("plot_ND", "SCALAR", d_ND_id);
//  javis_data_writer->registerPlotQuantity("plot_NA", "SCALAR", d_NA_id);
  javis_data_writer->registerPlotQuantity("contact_node", "SCALAR", d_contact_node_id);
  javis_data_writer->registerPlotQuantity("SiO2", "SCALAR", d_sio2_node_id);
  javis_data_writer->registerPlotQuantity("Silicon", "SCALAR", d_semi_node_id);
  javis_data_writer->registerPlotQuantity("semiinsbc","SCALAR", d_semiinsbc_node_id);
  javis_data_writer->registerPlotQuantity("plot_dopant", "SCALAR", d_dopant_id);
  javis_data_writer->registerPlotQuantity("node_variable", "SCALAR", d_node_var_plot_id);
  javis_data_writer->registerPlotQuantity("plot_E", "VECTOR", d_E_cell_id);
  javis_data_writer->registerPlotQuantity("plot_E_value", "COMPOSITE", d_E_cell_id);
  javis_data_writer->registerPlotQuantity("plot_J", "VECTOR", d_J_vec_id);
  javis_data_writer->registerPlotQuantity("plot_Jn", "VECTOR", d_Jn_vec_id);
  javis_data_writer->registerPlotQuantity("plot_Jp", "VECTOR", d_Jp_vec_id);
}

void NonlinearPoisson::setParameter(tbox::Pointer<tbox::Database> input_db)
{
}

/** get time step **/
double NonlinearPoisson::getPatchDt(hier::Patch<NDIM>& patch, const double time,
                                 const bool initial_time,
                                 const int flag_last_dt, const double last_dt,
                                 const string& component_name) {
    return Dt;
}

/** ***********************************************************************
 *  从输入数据库读入数据.
 *********************************************************************** **/
void NonlinearPoisson::getFromInput(tbox::Pointer<tbox::Database> db) {
#ifdef DEBUG_CHECK_ASSERTIONS
  TBOX_ASSERT(!db.isNull());
#endif
  /// 从输入数据库中读入有限元计算的若干参数.
  if (db->keyExists("element_type")) {
    d_element_type = db->getString("element_type");
  } else {
    TBOX_ERROR(d_object_name << ": "
               << " No key `element_type' found in data." << endl);
  }

  /** cell_point and cell_edge information for semi**/
  tbox::pout<<"d_element_type: "<<d_element_type<<endl;
  if(strcmp(d_element_type.c_str(),"LinearHexahedron")==0){
        cell_point=8;
        cell_edge=12;
        tbox::pout<<"The number of cell_point and cell_edge is 8 and 12;"<<endl;
  }else if(strcmp(d_element_type.c_str(),"LinearTetrahedron")==0){
        cell_point=4;
        cell_edge=6;
        tbox::pout<<"The number of cell_point and cell_edge is 4 and 6;"<<endl;
  }else{
      TBOX_ERROR("\n::registerModelVariable() : d_element_type "
                 << d_element_type <<" is not matched. cell_point and cell_edge can't be initiallized!"<<endl);
  }

  if (db->keyExists("DD_type")) {
    DD_type = db->getString("DD_type");
  } else {
    TBOX_ERROR(d_object_name << ": "
               << " No key `DD_type' found in data." << endl);
  }
  if(DD_type=="SG"||DD_type=="SUPG"){
      tbox::pout<<"DD_type:"<<DD_type<<endl;
  }else{
      TBOX_ERROR(d_object_name << ": "
                 << " No key `DD_type' found in data." << endl);
  }
}

/** ***********************************************************************
 *  输出数据成员到重启动数据库.
 *********************************************************************** **/
void NonlinearPoisson::putToDatabase(tbox::Pointer<tbox::Database> db) {

#ifdef DEBUG_CHECK_ASSERTIONS
  TBOX_ASSERT(!db.isNull());
#endif
  //d_dof_info->putToDatabase(db);
}


/** ***********************************************************************
 *  从重启动数据库读取数据.
 *********************************************************************** **/
void NonlinearPoisson::getFromRestart(tbox::Pointer<tbox::Database> db) {
  getFromInput(db);
  tbox::Pointer<tbox::Database> root_db =
    tbox::RestartManager::getManager()->getRootDatabase();

  tbox::Pointer<tbox::Database> sub_db
    = root_db->getDatabase(d_object_name);
  //d_dof_info->getFromDatabase(sub_db);
}

/** ***********************************************************************
 *  建立网格片上的矩阵和右端项.
 *********************************************************************** **/
/**
 * @brief NonlinearPoisson::buildCurConserveMatOnPatch: FEM-current continuity equation
 * @param patch
 * @param time
 * @param dt
 * @param component_name
 */
void NonlinearPoisson::buildCurConserveMatOnPatch(hier::Patch<NDIM> & patch,
                                               const double  time,
                                               const double  dt,
                                               const string& component_name){
    /** (形函数, 积分器, 单元)管理器**/
    tbox::Pointer<ElementManager<NDIM> >
      ele_manager = ElementManager<NDIM>::getManager();
    tbox::Pointer<BaseElement<NDIM> >
      ele = ele_manager->getElement(d_element_type);
    /** 取出本地PatchGeometry **/
    tbox::Pointer< hier::PatchGeometry<NDIM> > patch_geo = patch.getPatchGeometry();
    /** 取出本地PatchTopology **/
    tbox::Pointer< hier::PatchTopology<NDIM> > patch_top = patch.getPatchTopology();

    /** 取出本地Patch的结点坐标数组 **/
    tbox::Pointer<pdat::NodeData<NDIM, double> >
      node_coord = patch_geo->getNodeCoordinates();

    /** stiff matrix on patch **/
    tbox::Pointer<pdat::CSRMatrixData<NDIM,double> >
      mat_data = patch.getPatchData(d_curconserve_mat_id);
    /** cell - entity relationship **/
    tbox::Pointer<pdat::CellData<NDIM, int> >
      cell_mat = patch.getPatchData(d_cell_mat_id);
    /** cell - heter node relationship **/
    tbox::Pointer<pdat::CellData<NDIM, int> >
      cell_hnode = patch.getPatchData(d_cell_hnode_id);
    /** heter node relationship **/
    tbox::Pointer<pdat::NodeData<NDIM, int> >
      heter_node = patch.getPatchData(d_heter_node_id);
    /** metal relationship **/
    tbox::Pointer<pdat::NodeData<NDIM, int> >
      metal_node = patch.getPatchData(d_metal_node_id);
    /** 获取本地单元(边, 面)周围结点的索引关系 **/
    tbox::Array<int> can_extent, can_indices;
    patch_top->getCellAdjacencyNodes(can_extent, can_indices);
    int num_cells = patch.getNumberOfCells(1);

    /** 取出自由度映射信息 **/
    int * dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);

    /** material list **/
    tbox::Pointer<MaterialManager<NDIM> >
            mat_manager = MaterialManager<NDIM>::getManager();

    /** variable function **/
    tbox::Pointer<Variables_Functions>
            Var_func = new Variables_Functions();

    /** Traversing Unit **/
    for(int i = 0; i < num_cells; ++i){

        int n_vertex = can_extent[i+1] - can_extent[i];
        tbox::Array<hier::DoubleVector<NDIM> > vertex(n_vertex);
        tbox::Array<double> T_node(n_vertex);
        tbox::Array<int> mapping(n_vertex);
        tbox::Array<int> node_n(n_vertex);

        /** 取出单元结点坐标，以及填写映射值 **/
        for (int i1 = 0, j = can_extent[i]; i1 < n_vertex; ++i1, ++j) {
            mapping[i1] = dof_map[can_indices[j]];
            for(int k = 0; k< NDIM; ++k){
              vertex[i1][k] = (*node_coord)(k, can_indices[j])*m;
            }
            T_node[i1]=T_amb;                    /** 取得单元节点的温度值 **/
            node_n[i1]=can_indices[j];
        }

        if(material_list[(*cell_mat)(0,i)].mat_type == "metal"){

            tbox::Pointer<tbox::Matrix<double> > ele_matG = new tbox::Matrix<double>();
            ele_matG->resize(n_vertex, n_vertex);
            for (int l = 0; l < n_vertex; ++l){
                for (int j = 0; j < n_vertex; ++j){
                    (*ele_matG)(l,j) = 0.0;
                }
            }
            ele->buildElemMatGG(vertex, dt, time, ele_matG);

            tbox::Pointer<BaseMaterial<NDIM> > mat_cell =
                    mat_manager->getMaterial(material_list[(*cell_mat)(0,i)].mat_name);

            double T_cell = Var_func->node2ElemVariable(vertex, T_node, dt, time);

            double sigma= mat_cell->electroSigma(T_cell);

            int row = 0, col = 0;
            for (int l = 0; l < n_vertex; ++l){
                row = mapping[l];
                for(int j = 0; j < n_vertex; ++j){
                    col = mapping[j];
                    double mat_tmp = sigma*(*ele_matG)(l,j);
                    mat_data->addMatrixValue(row, col, mat_tmp);
                }
            }
        }else{
            double diag_value = 1.0;
            int row;
            if((*cell_hnode)(0,i) == 2){
                for (int l = 0; l < n_vertex; ++l){
                    if((*metal_node)(0,node_n[l]) == 0){
                        if((*heter_node)(0,node_n[l]) == 1){
                            row = mapping[l] + 1;
                        }else{
                            row = mapping[l];
                        }
                        mat_data->addMatrixValue(row, row, diag_value);
                    }
                }
            }else{
                for (int l = 0; l < n_vertex; ++l){
                    row = mapping[l];
                    if((*metal_node)(0,node_n[l]) == 0){
                        mat_data->addMatrixValue(row, row, diag_value);
                    }
                }
            }
        }
    }
    mat_data->assemble();
}

/**
 * @brief NonlinearPoisson::setCurConservePhysicalBC: Matrix
 * @param patch
 * @param time
 * @param dt
 * @param component_name
 */
void NonlinearPoisson::setCurConservePhysicalBC(hier::Patch<NDIM> & patch,
                          const double  time,
                          const double  dt,
                          const string& component_name){
    /** 取出本地PatchGeometry **/
    tbox::Pointer< hier::PatchGeometry<NDIM> > patch_geo = patch.getPatchGeometry();
    /** 获取网格片矩阵对象 **/
    tbox::Pointer<pdat::CSRMatrixData<NDIM,double> >
      mat_data = patch.getPatchData(d_curconserve_mat_id);

    tbox::Pointer<pdat::NodeData<NDIM, int> >
      metal_node = patch.getPatchData(d_metal_node_id);

    tbox::Pointer<pdat::NodeData<NDIM, int> >
      heter_node = patch.getPatchData(d_heter_node_id);

    int num_nodes = patch.getNumberOfNodes(1);

    int * dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);

    /** 取出矩阵向量的核心数据结构的指针 **/
    int * row_start = mat_data->getRowStartPointer();
    int * col_idx   = mat_data->getColumnIndicesPointer();
    double * mat_val    = mat_data->getValuePointer();

    /** boundary condition source **/
    int num_bc = EBoundassignment.getSize();

    for(int i_bc=0; i_bc<num_bc; ++i_bc){
        if (patch_geo->hasEntitySet(EBoundassignment[i_bc].bc_no, hier::EntityUtilities::NODE)){
            const tbox::Array<int> &entity_idx = patch_geo->getEntityIndicesInSet(EBoundassignment[i_bc].bc_no,
                                                 hier::EntityUtilities::NODE,
                                                 num_nodes);
            if (EBoundassignment[i_bc].bc_type == "ohmic"){
                int size = entity_idx.getSize();
                for (int i = 0; i < size; ++i) {
                    int index = dof_map[entity_idx[i]];
                    /** <mat_value at boundary> **/
                    for (int j = row_start[index]; j < row_start[index+1]; ++j){
                        if (col_idx[j] == index){
                            mat_val[j] = 1.0;
                        }
                        else {
                            mat_val[j] = 0.0;
                        }
                    }
                }
            }else if (EBoundassignment[i_bc].bc_type == "schottky"){
                int size = entity_idx.getSize();
                for (int i = 0; i < size; ++i) {
                    int index = dof_map[entity_idx[i]];
                    /** <mat_value at boundary> **/
                    for (int j = row_start[index]; j < row_start[index+1]; ++j){
                        if (col_idx[j] == index){
                            mat_val[j] = 1.0;
                        }
                        else {
                            mat_val[j] = 0.0;
                        }
                    }
                }
            }
        }
    }

    /** non-metal entity **/
    for(int i_n = 0; i_n < num_nodes; ++i_n){
        if((*metal_node)(0,i_n) == 0){
            int index = dof_map[i_n];
            for(int j = row_start[index]; j < row_start[index+1]; ++j){
                if (col_idx[j] == index){
                    mat_val[j] = 1.0;
                }else {
                    mat_val[j] = 0.0;
                }
            }
            if((*heter_node)(0,i_n) == 1){
                int index = dof_map[i_n] + 1;
                for(int j = row_start[index]; j < row_start[index+1]; ++j){
                    if (col_idx[j] == index){
                        mat_val[j] = 1.0;
                    }else{
                        mat_val[j] = 0.0;
                    }
                }
            }
        }
    }
}

/**
 * @brief NonlinearPoisson::applyCurConserveBCOnPatch: Vector
 * @param patch
 * @param time
 * @param dt
 * @param component_name
 */
void NonlinearPoisson::applyCurConserveBCOnPatch(hier::Patch<NDIM> & patch,
                    const double  time,
                    const double  dt,
                    const string& component_name){
    /** 取出本地PatchGeometry **/
    tbox::Pointer< hier::PatchGeometry<NDIM> > patch_geo = patch.getPatchGeometry();

    tbox::Pointer<pdat::VectorData<NDIM,double> >
      fval_data = patch.getPatchData(d_curconserve_fval_id);

    tbox::Pointer<pdat::NodeData<NDIM, double> >
      dop = patch.getPatchData(d_dopant_id);

    tbox::Pointer<pdat::NodeData<NDIM, int> >
      node_mat = patch.getPatchData(d_node_mat_id);

    int num_nodes = patch.getNumberOfNodes(1);

    int * dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);
    /** 取出矩阵向量的核心数据结构的指针 **/
    double * fval_val = fval_data->getPointer();

    /** setup the source target **/
    tbox::Pointer<bound_condition> bound_source = new bound_condition();

    /** material list **/
    tbox::Pointer<MaterialManager<NDIM> >
            mat_manager = MaterialManager<NDIM>::getManager();

    /** boundary condition source **/
    int num_bc = EBoundassignment.getSize();
    tbox::Array<double> bound_coe;
    bound_coe.resizeArray(9);                               /** vector contains params for source **/

    for(int i_bc=0; i_bc<num_bc; ++i_bc){
        if (patch_geo->hasEntitySet(EBoundassignment[i_bc].bc_no, hier::EntityUtilities::NODE)){
            const tbox::Array<int> &entity_idx = patch_geo->getEntityIndicesInSet(EBoundassignment[i_bc].bc_no,
                                                 hier::EntityUtilities::NODE,
                                                 num_nodes);
            /** load fixed volt value **/
            double app_volt = bound_source->source_assign(EBoundassignment[i_bc].sc_type,
                                                          time, dt,
                                                          EBoundassignment[i_bc].sc_coe);
            int size = entity_idx.getSize();
            string contact_type;

            /** assign boundary value to solution vector **/
            if(EBoundassignment[i_bc].bc_type == "ohmic"){

                contact_type = "Ohmic_contact";

                for (int i = 0; i < size; ++i) {
                    /** setup material target **/
                    int dom_no;
                    dom_no = 0;
                    for(int i_m = 0; i_m < material_list.size(); ++i_m){
                        if(material_list[i_m].mat_type == "semiconductor"){
                            if((*node_mat)(i_m,entity_idx[i])==1){
                                dom_no = i_m;
                            }
                        }
                    }
                    tbox::Pointer<BaseMaterial<NDIM> > mat_dom =
                            mat_manager->getMaterial(material_list[dom_no].mat_name);

                    int index = dof_map[entity_idx[i]];

                    bound_coe[0] = app_volt;
                    bound_coe[1] = T_amb;
                    bound_coe[2] = mat_dom->effectiveAffinity(T_amb);
                    bound_coe[3] = mat_dom->bandGap(T_amb);
                    bound_coe[4] = mat_dom->effectiveDensityNc(T_amb);
                    bound_coe[5] = mat_dom->effectiveDensityNv(T_amb);
                    bound_coe[6] = (*dop)(0,entity_idx[i]);
                    bound_coe[7] = 0.0;
                    bound_coe[8] = mat_dom->intrinsicDensityNi(T_amb);

                    fval_val[index] = bound_source->e_bound_value(contact_type, time, dt, kb, e, bound_coe);
                }
            }else if (EBoundassignment[i_bc].bc_type == "schottky"){

                contact_type = "Schottky_contact";

                for (int i = 0; i < size; ++i) {
                    /** setup material target **/
                    int dom_no;
                    dom_no = 0;
                    for(int i_m = 0; i_m < material_list.size(); ++i_m){
                        if(material_list[i_m].mat_type == "semiconductor"){
                            if((*node_mat)(i_m,entity_idx[i])==1){
                                dom_no = i_m;
                            }
                        }
                    }
                    tbox::Pointer<BaseMaterial<NDIM> > mat_dom =
                            mat_manager->getMaterial(material_list[dom_no].mat_name);

                    int index = dof_map[entity_idx[i]];
                    bound_coe[0] = app_volt;
                    bound_coe[1] = T_amb;
                    bound_coe[2] = EBoundassignment[i_bc].bc_value;
                    bound_coe[3] = mat_dom->bandGap(T_amb);
                    bound_coe[4] = mat_dom->effectiveDensityNc(T_amb);
                    bound_coe[5] = mat_dom->effectiveDensityNv(T_amb);
                    bound_coe[6] = (*dop)(0,entity_idx[i]);
                    bound_coe[7] = 0.0;
                    bound_coe[8] = mat_dom->intrinsicDensityNi(T_amb);

                    fval_val[index] = bound_source->e_bound_value(contact_type, time, dt, kb, e, bound_coe);
                }
            }
        }
    }
}


/**
 * @brief NonlinearPoisson::buildPoissonMatOnPatch: for initialization
 * @param patch
 * @param time
 * @param dt
 * @param component_name
 */
void NonlinearPoisson::buildPoissonMatOnPatch(hier::Patch<NDIM> & patch,
                                          const double  time,
                                          const double  dt,
                                          const string& component_name)
{
    //MPI_Barrier(MPI_COMM_WORLD);
    //cout<<"build matrix"<<endl;

    /** manager for shapefuction, integrator & element **/
    tbox::Pointer<ElementManager<NDIM> >
            ele_manager = ElementManager<NDIM>::getManager();
    tbox::Pointer<BaseElement<NDIM> >
            ele = ele_manager->getElement(d_element_type);

    /** gPatchGeometry **/
    tbox::Pointer< hier::PatchGeometry<NDIM> >
            patch_geo = patch.getPatchGeometry();
    /** get PatchTopology **/
    tbox::Pointer< hier::PatchTopology<NDIM> >
            patch_top = patch.getPatchTopology();

    /** get Patch node coord **/
    tbox::Pointer<pdat::NodeData<NDIM, double> >
            node_coord = patch_geo->getNodeCoordinates();
    /** stiff matrix on patch **/
    tbox::Pointer<pdat::CSRMatrixData<NDIM,double> >
            mat_data = patch.getPatchData(d_poisson_mat_id);

    tbox::Pointer<pdat::VectorData<NDIM, double> >
            V_soln = patch.getPatchData(d_poisson_sol_cur_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
            Fn_soln = patch.getPatchData(d_fermin_sol_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
            Fp_soln = patch.getPatchData(d_fermih_sol_id);

    /** cell - entity relationship **/
    tbox::Pointer<pdat::CellData<NDIM, int> >
            cell_mat = patch.getPatchData(d_cell_mat_id);
    /** cell - heter node relationship **/
    tbox::Pointer<pdat::CellData<NDIM, int> >
            cell_hnode = patch.getPatchData(d_cell_hnode_id);
    /** node - entity relationship **/
    tbox::Pointer<pdat::NodeData<NDIM, int> >
            heter_node = patch.getPatchData(d_heter_node_id);
    tbox::Pointer<pdat::NodeData<NDIM, int> >
            semiinsbc_node_data = patch.getPatchData(d_semiinsbc_node_id);

    /** 获取本地单元(边, 面)周围结点的索引关系 **/
    tbox::Array<int> can_extent, can_indices;
    patch_top->getCellAdjacencyNodes(can_extent, can_indices);

    int num_cells = patch.getNumberOfCells(1);
    int num_nodes = patch.getNumberOfNodes(1);

    int * dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);
    int * dof_map_semi = d_dof_semi_info->getDOFMapping(patch, hier::EntityUtilities::NODE);

    tbox::Pointer<pdat::CellData<NDIM, double> > bas_integ;
    tbox::Pointer<pdat::CellData<NDIM, double> > bas_rhs;
    if(save_bas){
        bas_integ = patch.getPatchData(d_bas_integ_id);
        bas_rhs = patch.getPatchData(d_bas_rhs_id);
    }

    /** material list **/
    tbox::Pointer<MaterialManager<NDIM> >
            mat_manager = MaterialManager<NDIM>::getManager();

    /** variable function **/
    tbox::Pointer<Variables_Functions> Var_func = new Variables_Functions();

    #if 0
    cout<<"m:"<<m<<endl;
    TBOX_ERROR("break");
    #endif

    /** Traversing Unit **/
    for(int i = 0; i < num_cells; ++i){
        if(!((material_list[(*cell_mat)(0, i)].mat_type == "semiconductor")||(material_list[(*cell_mat)(0, i)].mat_type == "dielectric"))) continue;
        /** get node information **/
        int n_vertex = can_extent[i+1] - can_extent[i];
        tbox::Array<hier::DoubleVector<NDIM> > vertex(n_vertex);
        tbox::Array<int> node_n(n_vertex);
        tbox::Array<int> mapping(n_vertex);
        tbox::Array<int> mapping_semi(n_vertex);
        tbox::Array<double> T_node(n_vertex);

        for (int i1 = 0, j = can_extent[i]; i1 < n_vertex; ++i1, ++j) {
            mapping[i1] = dof_map[can_indices[j]];
            mapping_semi[i1] = dof_map_semi[can_indices[j]];
            for(int k = 0; k< NDIM; ++k){
                vertex[i1][k] = (*node_coord)(k, can_indices[j])*m;
            }
            node_n[i1] = can_indices[j];
            T_node[i1] = T;
        }

        /** initialize element matrix **/
        tbox::Pointer<tbox::Matrix<double> > ele_mat_E = new tbox::Matrix<double>();

        ele_mat_E->resize(n_vertex, n_vertex);

        for (int l = 0; l < n_vertex; ++l){
            for (int j = 0; j < n_vertex; ++j){
                (*ele_mat_E)(l,j) = 0.0;
            }
        }

        /** caclate element matrix **/

        /** get material for cell **/
        tbox::Pointer<BaseMaterial<NDIM> > mat_cell =
                mat_manager->getMaterial(material_list[(*cell_mat)(0,i)].mat_name);
        //cout<<"mat_name:"<<material_list[(*cell_mat)(0,i)].mat_name<<endl;

        double T_cell = Var_func->node2ElemVariable(vertex, T_node, dt, time);

        double cell_epsr = mat_cell->electroPermittivity(T_cell);

        double n_derive_tmp, p_derive_tmp;
        tbox::Array<double> n_derive_dens(n_vertex);
        tbox::Array<double> p_derive_dens(n_vertex);
        for(int i_v = 0; i_v < n_vertex; ++i_v){
            n_derive_dens[i_v] = 0.0;
            p_derive_dens[i_v] = 0.0;
        }

        double T_tmp, V_tmp, fn_tmp, fp_tmp, nc_tmp, nv_tmp, affin_tmp, eg_tmp;
        tbox::Array<tbox::Array<int> > r_node = ele->get_r_node();

        tbox::Array<tbox::Array<double> > r_integ = getInteg(n_vertex,save_bas,i,ele,bas_integ,vertex,time,dt);

        for(int i_v = 0; i_v < cell_point; ++i_v){
            for(int i_e = 0; i_e < cell_edge; ++i_e){
                int n_1 = r_node[i_e][0];
                int n_2 = r_node[i_e][1];

                double cn_1 = cell_epsr*(r_integ[i_e][i_v]);
                double cn_2 = cell_epsr*(r_integ[i_e][i_v]);

                (*ele_mat_E)(i_v,n_1) = (*ele_mat_E)(i_v,n_1) - cn_1;
                (*ele_mat_E)(i_v,n_2) = (*ele_mat_E)(i_v,n_2) + cn_2;
            }
        }

        for(int l=0; l<n_vertex; ++l){
            T_tmp = T;
            V_tmp = (*V_soln)(mapping_semi[l]);
            fn_tmp = (*Fn_soln)(mapping[l]);
            fp_tmp = (*Fp_soln)(mapping[l]);

            if((material_list[(*cell_mat)(0,i)].mat_type == "dielectric")&&((*semiinsbc_node_data)(0,node_n[l]) == 1)){
                V_tmp = (*V_soln)(mapping_semi[l]+3);
                fn_tmp = (*Fn_soln)(mapping[l]+1);
                fp_tmp = (*Fp_soln)(mapping[l]+1);
            }

            nc_tmp = mat_cell->effectiveDensityNc(T_amb)*pow(m,-3);
            nv_tmp = mat_cell->effectiveDensityNv(T_amb)*pow(m,-3);
            affin_tmp = mat_cell->effectiveAffinity(T_amb);
            eg_tmp = mat_cell->bandGap(T_amb);

            Var_func->dummyElectron(V_tmp, fn_tmp, true, affin_tmp, T_tmp, kb, e, nc_tmp, n_derive_tmp);
            Var_func->dummyHole(V_tmp, fp_tmp, eg_tmp, true, affin_tmp, T_tmp, kb, e, nv_tmp, p_derive_tmp);

            n_derive_dens[l] = n_derive_tmp;
            p_derive_dens[l] = p_derive_tmp;
        }

        tbox::Array<double> value(n_vertex);
        for (int i_n = 0; i_n < n_vertex; ++i_n){
            value[i_n] = (e/eps0)*(n_derive_dens[i_n] - p_derive_dens[i_n]);
        }

        tbox::Pointer<tbox::Vector<double> > cv_vol = getRhs(n_vertex,save_bas,i,ele,bas_rhs,value,vertex,time,dt);

        for(int i_n=0; i_n<n_vertex; i_n++){
            (*ele_mat_E)(i_n,i_n) = (*ele_mat_E)(i_n,i_n) - (*cv_vol)[i_n];
        }

        #if 0
        for(int i_n=0; i_n<n_vertex; i_n++){
            cout<<(*ele_mat_E)(i_n,i_n)<<endl;
        }
        #endif

        int row = 0, col = 0;
        for (int l = 0; l < n_vertex; ++l){
            row = mapping[l];
            for(int j = 0; j < n_vertex; ++j){
                col = mapping[j];
                if((material_list[(*cell_mat)(0,i)].mat_type == "dielectric")&&((*semiinsbc_node_data)(0,node_n[j]) == 1))
                    col = mapping[j]+1;
                double v_mat_tmp = (*ele_mat_E)(l,j);
                mat_data->addMatrixValue(row, col, v_mat_tmp);
            }
        }
    }

    /** assign boundary to heterjunction nodes **/
    for (int i_n = 0; i_n < num_nodes; ++i_n){
        int index = dof_map[i_n];
        if ((*heter_node)(0,i_n) == 1){
            mat_data->addMatrixValue(index+1,index+1,1.0);
            mat_data->addMatrixValue(index+1,index,-1.0);
            //cout<<"组装边界条件成功"<<endl;
        }
        if ((*semiinsbc_node_data)(0,i_n) == 1){
            mat_data->addMatrixValue(index+1,index+1,1.0);
            mat_data->addMatrixValue(index+1,index,-1.0);
            //cout<<"组装边界条件成功"<<endl;
        }
    }
    mat_data->assemble();
}

/**
 * @brief NonlinearPoisson::buildPoissonFXOnPatch
 * @param patch
 * @param time
 * @param dt
 * @param component_name
 */
void NonlinearPoisson::buildPoissonFXOnPatch(hier::Patch<NDIM> & patch,
                                      const double  time,
                                      const double  dt,
                                      const string& component_name)
{
    /** (形函数, 积分器, 单元)管理器 **/
    tbox::Pointer<ElementManager<NDIM> >
            ele_manager = ElementManager<NDIM>::getManager();
    tbox::Pointer<BaseElement<NDIM> > ele = ele_manager->getElement(d_element_type);

    /** 取出本地PatchGeometry **/
    tbox::Pointer< hier::PatchGeometry<NDIM> > patch_geo = patch.getPatchGeometry();
    /** 取出本地PatchTopology **/
    tbox::Pointer< hier::PatchTopology<NDIM> > patch_top = patch.getPatchTopology();
    /** 取出本地Patch的结点坐标数组 **/
    tbox::Pointer<pdat::NodeData<NDIM, double> >
            node_coord = patch_geo->getNodeCoordinates();
    /** right side vector **/
    tbox::Pointer<pdat::VectorData<NDIM,double> >
      vec_data = patch.getPatchData(d_poisson_fval_id);
    double *fval_val = vec_data->getPointer();

    /** get node values **/
    tbox::Pointer<pdat::VectorData<NDIM, double> >
            V_soln = patch.getPatchData(d_poisson_sol_cur_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
            Fn_soln = patch.getPatchData(d_fermin_sol_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
            Fp_soln = patch.getPatchData(d_fermih_sol_id);
    tbox::Pointer<pdat::NodeData<NDIM, double> >
            dop = patch.getPatchData(d_dopant_id);
    tbox::Pointer<pdat::NodeData<NDIM, double> >
            pp = patch.getPatchData(d_polarizedp_id);
    /** cell - entity relationship **/
    tbox::Pointer<pdat::CellData<NDIM, int> >
            cell_mat = patch.getPatchData(d_cell_mat_id);
    /** cell - entity relationship **/
    tbox::Pointer<pdat::CellData<NDIM, int> >
            cell_hnode = patch.getPatchData(d_cell_hnode_id);

    /** node - entity relationship **/
    tbox::Pointer<pdat::NodeData<NDIM, int> >
            heter_node = patch.getPatchData(d_heter_node_id);
    tbox::Pointer<pdat::NodeData<NDIM, int> >
            semiinsbc_node_data = patch.getPatchData(d_semiinsbc_node_id);

    tbox::Pointer<pdat::CellData<NDIM, double> > bas_integ;
    tbox::Pointer<pdat::CellData<NDIM, double> > bas_rhs;
    if(save_bas){
        bas_integ = patch.getPatchData(d_bas_integ_id);
        bas_rhs = patch.getPatchData(d_bas_rhs_id);
    }

    /** 获取本地单元(边, 面)周围结点的索引关系 **/
    tbox::Array<int> can_extent, can_indices;
    patch_top->getCellAdjacencyNodes(can_extent, can_indices);

    int num_cells = patch.getNumberOfCells(1);
    int num_nodes = patch.getNumberOfNodes(1);

    int * dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);
    int * dof_map_semi = d_dof_semi_info->getDOFMapping(patch, hier::EntityUtilities::NODE);

    /** material list **/
    tbox::Pointer<MaterialManager<NDIM> >
            mat_manager = MaterialManager<NDIM>::getManager();

    /** setup target for class 'Variables_Functions' **/
    tbox::Pointer<Variables_Functions> Var_func = new Variables_Functions();

    vec_data->fillAll(0.0);   // reset the vector data

    /** Trasversing Unit **/
    for(int i = 0; i < num_cells; ++i){
        if(!((material_list[(*cell_mat)(0, i)].mat_type == "semiconductor")||(material_list[(*cell_mat)(0, i)].mat_type == "dielectric"))) continue;

        int n_vertex = can_extent[i+1] - can_extent[i];
        tbox::Array<hier::DoubleVector<NDIM> > vertex(n_vertex);
        tbox::Array<double> dop_node(n_vertex);
        tbox::Array<double> pp_node(n_vertex);
        tbox::Array<int> e_node(n_vertex);
        tbox::Array<int> mapping(n_vertex);
        tbox::Array<int> mapping_semi(n_vertex);
        tbox::Array<double> T_node(n_vertex);
        tbox::Array<double> V_node(n_vertex);
        tbox::Array<double> fn_node(n_vertex);
        tbox::Array<double> fp_node(n_vertex);

        for (int i1 = 0, j = can_extent[i]; i1 < n_vertex; ++i1, ++j) {
            mapping[i1] = dof_map[can_indices[j]];
            mapping_semi[i1] = dof_map_semi[can_indices[j]];
            for(int k = 0; k< NDIM; ++k){
                vertex[i1][k] = (*node_coord)(k, can_indices[j])*m;
            }
            e_node[i1]=can_indices[j];
            dop_node[i1]=(*dop)(0,can_indices[j]);
            pp_node[i1]=(*pp)(0,can_indices[j]);
        }

        /** get material for cell **/
        tbox::Pointer<BaseMaterial<NDIM> > mat_cell =
                mat_manager->getMaterial(material_list[(*cell_mat)(0,i)].mat_name);

        double T_cell = 0.0;
        for (int i_v = 0; i_v < n_vertex; ++i_v){
            T_cell += T/n_vertex;
        }
        double cell_epsr = mat_cell->electroPermittivity(T_cell);

        /** carrier density **/
        double n_tmp, p_tmp;
        tbox::Array<double> n_dens(n_vertex);
        tbox::Array<double> p_dens(n_vertex);
        tbox::Array<double> d_dens(n_vertex);
        tbox::Array<double> pp_dens(n_vertex);
        tbox::Array<double> value(n_vertex);

        for (int i_v = 0; i_v < n_vertex; ++i_v){
            n_dens[i_v] = 0.0;
            p_dens[i_v] = 0.0;
            d_dens[i_v] = 0.0;
            pp_dens[i_v] = 0.0;
            value[i_v] = 0.0;
        }

        double nc_tmp, nv_tmp, affin_tmp, eg_tmp;
        for(int l=0; l<n_vertex; ++l){
            T_node[l] = T;
            V_node[l] = (*V_soln)(mapping_semi[l]);
            fn_node[l] = (*Fn_soln)(mapping[l]);
            fp_node[l] = (*Fp_soln)(mapping[l]);
            if((material_list[(*cell_mat)(0,i)].mat_type == "dielectric")&&((*semiinsbc_node_data)(0,e_node[l]) == 1)){
                V_node[l] = (*V_soln)(mapping_semi[l]+3);
                fn_node[l] = (*Fn_soln)(mapping[l]+1);
                fp_node[l] = (*Fp_soln)(mapping[l]+1);
            }
        }

        for (int i_n = 0; i_n < n_vertex; ++i_n){
            nc_tmp = mat_cell->effectiveDensityNc(T_amb)*pow(m,-3);
            nv_tmp = mat_cell->effectiveDensityNv(T_amb)*pow(m,-3);
            affin_tmp = mat_cell->effectiveAffinity(T_node[i_n]);
            eg_tmp = mat_cell->bandGap(T_node[i_n]);

            Var_func->dummyElectron(V_node[i_n], fn_node[i_n], false, affin_tmp, T_node[i_n], kb, e, nc_tmp, n_tmp);
            Var_func->dummyHole(V_node[i_n], fp_node[i_n], eg_tmp, false, affin_tmp, T_node[i_n], kb, e, nv_tmp, p_tmp);

            n_dens[i_n] = n_tmp;
            p_dens[i_n] = p_tmp;
            pp_dens[i_n] = pp_node[i_n];
            if(material_list[(*cell_mat)(0,i)].mat_type == "semiconductor") {
                d_dens[i_n] = dop_node[i_n];
            }else if (material_list[(*cell_mat)(0,i)].mat_type == "dielectric") {
                d_dens[i_n] = Not;
            }
        }


        for (int i_n = 0; i_n < n_vertex; ++i_n){
            value[i_n] = (e/eps0)*(n_dens[i_n] - p_dens[i_n] - d_dens[i_n]);
        }

        /** building element vector **/
        /** element vector RHS **/
        tbox::Pointer<tbox::Vector<double> > ele_vec = getRhs(n_vertex,save_bas,i,ele,bas_rhs,value,vertex,time,dt);

        #if 0
        for(int l=0;l<n_vertex;l++){
            cout<<"(*ele_vec)[l]:"<<(*ele_vec)[l]<<endl;
        }
        #endif

        tbox::Array<tbox::Array<double> > r_integ = getInteg(n_vertex,save_bas,i,ele,bas_integ,vertex,time,dt);

        tbox::Array<tbox::Array<int> > r_node = ele->get_r_node();

        /** initialize element matrix **/
        tbox::Pointer<tbox::Matrix<double> > ele_mat_E = new tbox::Matrix<double>();

        ele_mat_E->resize(n_vertex, n_vertex);

        for (int l = 0; l < n_vertex; ++l){
            for (int j = 0; j < n_vertex; ++j){
                (*ele_mat_E)(l,j) = 0.0;
            }
        }

        for(int i_v = 0; i_v < cell_point; ++i_v){
            for(int i_e = 0; i_e < cell_edge; ++i_e){
                int n_1 = r_node[i_e][0];
                int n_2 = r_node[i_e][1];

                double cn_1 = (r_integ[i_e][i_v]);
                double cn_2 = (r_integ[i_e][i_v]);

                (*ele_mat_E)(i_v,n_1) = (*ele_mat_E)(i_v,n_1) - cn_1;
                (*ele_mat_E)(i_v,n_2) = (*ele_mat_E)(i_v,n_2) + cn_2;
            }
        }

        #if 0
        for(int i_v = 0; i_v < cell_point; ++i_v){
            for(int i_n = 0; i_n < cell_point; ++i_n){
                cout<<"(*ele_mat_E)(i_v,i_n):"<<i_v<<":"<<i_n<<":"<<(*ele_mat_E)(i_v,i_n)<<endl;
            }
        }
        #endif

        for (int l = 0; l < n_vertex; ++l){
            double vec_tmp = 0.0;
            //double vec_test = 0.0;
            for(int j = 0; j < n_vertex; ++j){
                vec_tmp += - cell_epsr*(*ele_mat_E)(l,j)*V_node[j];
            }
            (*ele_vec)[l] += vec_tmp;
        }

        /** construct fval vector **/
        for(int i2= 0; i2 < n_vertex; ++i2)
            vec_data->addVectorValue(mapping[i2], (*ele_vec)[i2]);
    }
    /** assign boundary to heterjunction nodes **/
    for (int i_n = 0; i_n < num_nodes; ++i_n){
        if ((*semiinsbc_node_data)(0,i_n) == 1){
            int index = dof_map[i_n];
            fval_val[index+1] = 0;
            //cout<<"组装边界条件成功"<<endl;
        }
    }
}

/**
 * @brief NonlinearPoisson::setPoissonPhysicalBC: matrix
 * @param patch
 * @param time
 * @param dt
 * @param component_name
 */
void NonlinearPoisson::setPoissonPhysicalBC(hier::Patch<NDIM> & patch,
                          const double  time,
                          const double  dt,
                          const string& component_name){
    /** (形函数, 积分器, 单元)管理器 **/
    tbox::Pointer<ElementManager<NDIM> >
            ele_manager = ElementManager<NDIM>::getManager();

    /** 取出本地PatchGeometry **/
    tbox::Pointer< hier::PatchGeometry<NDIM> > patch_geo = patch.getPatchGeometry();
    /** get PatchTopology **/
    tbox::Pointer< hier::PatchTopology<NDIM> > patch_top = patch.getPatchTopology();
    /** 获取网格片矩阵对象 **/
    tbox::Pointer<pdat::CSRMatrixData<NDIM,double> >
            mat_data = patch.getPatchData(d_poisson_mat_id);

    /** obtain the right side vector **/
    tbox::Pointer<pdat::VectorData<NDIM,double> >
            fval_data = patch.getPatchData(d_poisson_fval_id);

    /** get node values **/
    tbox::Pointer<pdat::VectorData<NDIM, double> >
            V_soln = patch.getPatchData(d_poisson_sol_cur_id);

    /** assign boundary condition to right side vector **/
    double * fval_val = fval_data->getPointer();

    /** node - entity relationship **/
    tbox::Pointer<pdat::NodeData<NDIM, int> >
            heter_node = patch.getPatchData(d_heter_node_id);

    tbox::Array<int> face_node_ext, face_node_idx;
    patch_top->getFaceAdjacencyNodes(face_node_ext, face_node_idx);

    int num_nodes = patch.getNumberOfNodes(1);

    int * dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);
    int * dof_map_semi = d_dof_semi_info->getDOFMapping(patch, hier::EntityUtilities::NODE);
    /** 取出矩阵向量的核心数据结构的指针 **/
    int * row_start = mat_data->getRowStartPointer();
    int * col_idx   = mat_data->getColumnIndicesPointer();
    double * mat_val    = mat_data->getValuePointer();

    /** 取出本地Patch的结点坐标数组 **/
    tbox::Pointer<pdat::NodeData<NDIM, double> > node_coord = patch_geo->getNodeCoordinates();

    /** setup the source target **/
    tbox::Pointer<bound_condition> bound_source = new bound_condition();

    /** assign boundary condition to system matrix **/
    int num_bc = EBoundassignment.getSize();
    for(int i_bc=0; i_bc<num_bc; ++i_bc){
        if (patch_geo->hasEntitySet(EBoundassignment[i_bc].bc_no, hier::EntityUtilities::NODE)){
            const tbox::Array<int> &entity_idx = patch_geo->getEntityIndicesInSet(EBoundassignment[i_bc].bc_no,
                                                 hier::EntityUtilities::NODE,
                                                 num_nodes);
            if ((EBoundassignment[i_bc].bc_type == "ohmic")||(EBoundassignment[i_bc].bc_type == "gate")){
                int size = entity_idx.getSize();
                for (int i = 0; i < size; ++i) {
                    int index = dof_map[entity_idx[i]];
                    fval_val[index] = 0.0;
                    /** <mat_value at boundary> **/
                    for (int j = row_start[index]; j < row_start[index+1]; ++j){
                        if (col_idx[j] == index){
                            mat_val[j] = 1.0;   /** 对角线元素 **/
                        }
                        else {
                            mat_val[j] = 0.0;   /**< 非对角线元素 */
                        }
                    }
                }
            }else if (EBoundassignment[i_bc].bc_type == "schottky"){
                int size = entity_idx.getSize();
                for (int i = 0; i < size; ++i) {
                    int index = dof_map[entity_idx[i]];
                    fval_val[index] = 0.0;
                    /** <mat_value at boundary> **/
                    for (int j = row_start[index]; j < row_start[index+1]; ++j){
                        if (col_idx[j] == index){
                            mat_val[j] = 1.0;   /**< 对角线元素 */
                        }
                        else {
                            mat_val[j] = 0.0;   /**< 非对角线元素 */
                        }
                    }
                }
            }else if (EBoundassignment[i_bc].bc_type == "gate-simplified"){
                /** load fixed volt value **/
                double dins = EBoundassignment[i_bc].sc_coe[4]*m;
                double epsilon_ins = EBoundassignment[i_bc].sc_coe[5];
                double Vapp = bound_source->source_assign(EBoundassignment[i_bc].sc_type,time, dt,EBoundassignment[i_bc].sc_coe);
                double WORKFUNC = EBoundassignment[i_bc].sc_coe[3];
                double local_Qt = EBoundassignment[i_bc].sc_coe[6];
                double Qins = local_Qt*pow(m,-2);
                if(int(time/dt)==0){
                    Qins = 0;
                    //Vapp = 0;
                }
                double value_rhs = (e*Qins/eps0 + epsilon_ins * (Vapp - WORKFUNC)/dins);
                tbox::Pointer<BaseElement<NDIM> > ele_quadrangle = ele_manager->getElement("LinearQuadrangle");
                int num_faces = patch.getNumberOfEntities(hier::EntityUtilities::FACE, 1);
                int num_face = EBoundassignment[i_bc].bc_face.getSize();
                for (int i_f = 0; i_f < num_face; ++i_f){
                    if(patch_geo->hasEntitySet(EBoundassignment[i_bc].bc_face[i_f], hier::EntityUtilities::FACE)){
                        //cout<<"get the gate-simplified face:"<<EBoundassignment[i_bc].bc_face[i_f]<<endl;
                        const tbox::Array<int> & face_id = patch_geo->getEntityIndicesInSet(EBoundassignment[i_bc].bc_face[i_f],
                                                           hier::EntityUtilities::FACE,
                                                           num_faces);
                        int size_f = face_id.size();
                        for (int i_e = 0; i_e < size_f; ++i_e){
                            int n_fvertex = face_node_ext[face_id[i_e] + 1] - face_node_ext[face_id[i_e]];
                            tbox::Array<hier::DoubleVector<NDIM> > face_vertex(n_fvertex);
                            tbox::Array<int> node_n(n_fvertex);
                            tbox::Array<int> mapping(n_fvertex);
                            tbox::Array<int> mapping_semi(n_fvertex);
                            tbox::Array<double> value(n_fvertex);
                            tbox::Array<double> V_node(n_fvertex);

                            /** inform for node **/
                            for (int i1 = 0, j = face_node_ext[face_id[i_e]]; i1 < n_fvertex; ++i1, ++j) {
                                node_n[i1] = face_node_idx[j];
                                mapping[i1] = dof_map[face_node_idx[j]];
                                mapping_semi[i1] = dof_map_semi[face_node_idx[j]];
                                V_node[i1] = (*V_soln)(mapping_semi[i1]);
                                for(int k = 0; k< NDIM; ++k){
                                    face_vertex[i1][k] = (*node_coord)(k, face_node_idx[j])*m;
                                }
                            }

                            tbox::Pointer<tbox::Vector<double> > ele_vec = new tbox::Vector<double>();
                            ele_vec->resize(n_fvertex);

                            for (int i_n = 0; i_n < n_fvertex; ++i_n){
                                value[i_n] = -epsilon_ins/dins;
                            }
                            /** building element vector **/
                            ele_quadrangle->buildCVFEMRhs(face_vertex, value, dt, time, ele_vec);

                            /** add up the element matrix for system matrix **/
                            int index = 0;
                            for (int l = 0; l < n_fvertex; ++l){
                                index = mapping[l];
                                (*mat_data)(index,index) = (*mat_data)(index,index) + (*ele_vec)[l];
                            }

                            for (int i_n = 0; i_n < n_fvertex; ++i_n){
                                value[i_n] = - value_rhs + epsilon_ins*V_node[i_n]/dins;
                            }
                            /** building element vector **/
                            ele_quadrangle->buildCVFEMRhs(face_vertex, value, dt, time, ele_vec);

                            /** construct fval vector **/
                            for(int i2= 0; i2 < n_fvertex; ++i2)
                                fval_val[mapping[i2]] += (*ele_vec)[i2];
                        }
                    }else{
                        cout<<"can't get the gate-simplified face:"<<EBoundassignment[i_bc].bc_face[i_f]<<endl;
                    }
                }
            }
        }
    }
}

/**
 * @brief NonlinearPoisson::PoissonMaxError: to find the maximum change in potential
 * @param vector
 * @param patch
 * @param time
 * @param dt
 * @param component_name
 */
void NonlinearPoisson::PoissonMaxError(double* vector,
                                    hier::Patch<NDIM>& patch,
                                    const double time, const double dt,
                                    const string& component_name){

    tbox::Pointer<pdat::VectorData<NDIM, double> >
        poisson_sol = patch.getPatchData(d_poisson_sol_id);
    tbox::Pointer<pdat::NodeData<NDIM, int> >
          siclass_node = patch.getPatchData(d_siclass_node_id);
    tbox::Pointer<pdat::NodeData<NDIM, int> >
          semiinsbc = patch.getPatchData(d_semiinsbc_node_id);

    int * dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);

    int num_n = patch.getNumberOfNodes(1);
    double error;
    error = 0;
    for (int i_n = 0; i_n<num_n; ++i_n){
        if((*siclass_node)(0,i_n)==0) continue;
        int mapping = dof_map[i_n];
        double tmp_error = fabs((*poisson_sol)(mapping));
        if(error<tmp_error){
            error = tmp_error;
        }
        if((*semiinsbc)(0,i_n)==1){
            int mapping = dof_map[i_n]+1;
            double tmp_error = fabs((*poisson_sol)(mapping));
            if(error<tmp_error){
                error = tmp_error;
            }
        }
    }
    vector[0] = error;
}

/**
 * @brief NonlinearPoisson::transportGummelPoissonSol: potential-from current to stetch
 * @param patch
 * @param time
 * @param dt
 * @param component_name
 */
void NonlinearPoisson::transportGummelPoissonSol(hier::Patch<NDIM>& patch,
                                        const double time,
                                        const double dt,
                                        const string& component_name){

    tbox::Pointer<pdat::VectorData<NDIM, double> >
      poisson_sol_cur = patch.getPatchData(d_poisson_sol_cur_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
      poisson_sol_delta = patch.getPatchData(d_poisson_sol_id);
    tbox::Pointer<pdat::NodeData<NDIM, int> >
          siclass_node = patch.getPatchData(d_siclass_node_id);
    tbox::Pointer<pdat::NodeData<NDIM, int> >
          semiinsbc = patch.getPatchData(d_semiinsbc_node_id);

    int * dof_map_semi = d_dof_semi_info->getDOFMapping(patch, hier::EntityUtilities::NODE);
    int * dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);

    int num_n = patch.getNumberOfNodes(1);

    for (int i_n = 0; i_n < num_n; ++i_n){
        if((*siclass_node)(0,i_n)==0) continue;
        int semi_mapping = dof_map_semi[i_n];
        int mapping = dof_map[i_n];
        double deta_phi = (*poisson_sol_delta)(mapping);

        (*poisson_sol_cur)(semi_mapping) = (*poisson_sol_cur)(semi_mapping) + deta_phi;

        if((*semiinsbc)(0,i_n)==1){
            int semi_mapping = dof_map_semi[i_n]+3;
            int mapping = dof_map[i_n]+1;
            double deta_phi = (*poisson_sol_delta)(mapping);

            (*poisson_sol_cur)(semi_mapping) = (*poisson_sol_cur)(semi_mapping) + deta_phi;
        }
    }
}

/** ********************************************************************
 * leftside matrix for carrier transport equations
 ******************************************************************** **/

/**
 * @brief NonlinearPoisson::buildCarMatOnPatch: CV-FEM-SG for carrier continuity equation
 * @param patch
 * @param time
 * @param dt
 * @param component_name
 */
void NonlinearPoisson::buildCarMatOnPatch(hier::Patch<NDIM> & patch,
                                          const double  time,
                                          const double  dt,
                                          const string& component_name){
    /** (形函数, 积分器, 单元)管理器 **/
    tbox::Pointer<ElementManager<NDIM> >
      ele_manager = ElementManager<NDIM>::getManager();
    tbox::Pointer<BaseElement<NDIM> > ele = ele_manager->getElement(d_element_type);

    /** 取出本地PatchGeometry **/
    tbox::Pointer< hier::PatchGeometry<NDIM> > patch_geo = patch.getPatchGeometry();
    /** 取出本地PatchTopology **/
    tbox::Pointer< hier::PatchTopology<NDIM> > patch_top = patch.getPatchTopology();
    /** 取出本地Patch的结点坐标数组 **/
    tbox::Pointer<pdat::NodeData<NDIM, double> >
            node_coord = patch_geo->getNodeCoordinates();

    /** stiff matrix on patch **/
    tbox::Pointer<pdat::CSRMatrixData<NDIM,double> >
            n_mat_data = patch.getPatchData(d_carn_mat_id);
    tbox::Pointer<pdat::CSRMatrixData<NDIM,double> >
            h_mat_data = patch.getPatchData(d_carh_mat_id);

    /** get node values **/
    tbox::Pointer<pdat::VectorData<NDIM, double> >
            V_soln = patch.getPatchData(d_poisson_sol_cur_id);
    tbox::Pointer<pdat::VectorData<NDIM,double> >
            Fn_soln = patch.getPatchData(d_fermin_sol_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
            Fh_soln = patch.getPatchData(d_fermih_sol_id);

    tbox::Pointer<pdat::NodeData<NDIM, double> >
            dop_data = patch.getPatchData(d_dopant_id);
    /** heterojunction node inform **/
    tbox::Pointer<pdat::NodeData<NDIM, int> >
            heter_n = patch.getPatchData(d_heter_node_id);
    tbox::Pointer<pdat::NodeData<NDIM, int> >
            semiinsbc_n = patch.getPatchData(d_semiinsbc_node_id);
    /** cell - heterojunction relation **/
    tbox::Pointer<pdat::CellData<NDIM, int> >
            cell_hnode = patch.getPatchData(d_cell_hnode_id);
    /** cell - entity relationship **/
    tbox::Pointer<pdat::CellData<NDIM, int> >
            cell_mat = patch.getPatchData(d_cell_mat_id);
    /** face - entity relationship **/
    tbox::Pointer<pdat::FaceData<NDIM, int> >
            face_mat = patch.getPatchData(d_face_mat_id);

    tbox::Pointer<pdat::CellData<NDIM, double> > bas_integ;
    tbox::Pointer<pdat::CellData<NDIM, double> > bas_rhs;
    if(save_bas){
        bas_integ = patch.getPatchData(d_bas_integ_id);
        bas_rhs = patch.getPatchData(d_bas_rhs_id);
    }

    /** 获取本地单元(边, 面)周围结点的索引关系 **/
    tbox::Array<int> can_extent, can_indices;
    patch_top->getCellAdjacencyNodes(can_extent, can_indices);
    /** 面——单元邻接关系数据 **/
    tbox::Array<int>face_cell_ext,face_cell_idx;
    patch_top->getFaceAdjacencyCells(face_cell_ext,face_cell_idx);
    /** 单元节点邻接关系数据 **/
    tbox::Array<int>face_node_ext,face_node_idx;
    patch_top->getFaceAdjacencyNodes(face_node_ext,face_node_idx);

    //int num_faces = patch.getNumberOfFaces(1);
    int num_cells = patch.getNumberOfCells(1);

    int * dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);
    int * dof_map_semi = d_dof_semi_info->getDOFMapping(patch, hier::EntityUtilities::NODE);

    /** material list **/
    tbox::Pointer<MaterialManager<NDIM> >
            mat_manager = MaterialManager<NDIM>::getManager();

    tbox::Pointer<Variables_Functions> Var_func = new Variables_Functions();

    int num_edge = cell_edge;
    int num_ver = cell_point;
    /** Traversing Unit **/
    for(int i = 0; i < num_cells; ++i){
        if(!((material_list[(*cell_mat)(0, i)].mat_type == "semiconductor")||(material_list[(*cell_mat)(0, i)].mat_type == "dielectric"))) continue;

        int n_vertex = can_extent[i+1] - can_extent[i];
        tbox::Array<hier::DoubleVector<NDIM> > vertex(n_vertex);

        tbox::Array<int> mapping(n_vertex);
        tbox::Array<int> mapping_semi(n_vertex);
        tbox::Array<int> node_n(n_vertex);
        tbox::Array<double> T_node(n_vertex);
        tbox::Array<double> V_node(n_vertex);
        tbox::Array<double> Fn_node(n_vertex);
        tbox::Array<double> Fh_node(n_vertex);

        /** get material for cell **/
        tbox::Pointer<BaseMaterial<NDIM> > mat_cell =
                mat_manager->getMaterial(material_list[(*cell_mat)(0,i)].mat_name);

        /** 取出单元结点坐标，以及填写映射值 **/
        for (int i1 = 0, j = can_extent[i]; i1 < n_vertex; ++i1, ++j) {
            mapping[i1] = dof_map[can_indices[j]];
            mapping_semi[i1] = dof_map_semi[can_indices[j]];
            for(int k = 0; k< NDIM; ++k){
                vertex[i1][k] = (*node_coord)(k, can_indices[j])*m;
            }
            node_n[i1] = can_indices[j];
        }

        for (int i_n = 0; i_n < n_vertex; ++i_n){
            T_node[i_n]=T;
            V_node[i_n]=(*V_soln)(mapping_semi[i_n]);
            Fn_node[i_n] = (*Fn_soln)(mapping[i_n]);
            Fh_node[i_n] = (*Fh_soln)(mapping[i_n]);
            if((material_list[(*cell_mat)(0,i)].mat_type == "dielectric")&&((*semiinsbc_n)(0,node_n[i_n]) == 1)){
                V_node[i_n]=(*V_soln)(mapping_semi[i_n]+3);
                Fn_node[i_n] = (*Fn_soln)(mapping[i_n]+1);
                Fh_node[i_n] = (*Fh_soln)(mapping[i_n]+1);
            }
        }

        tbox::Pointer<tbox::Array<double> > edge_field = new tbox::Array<double >(cell_edge);
        ele->calculateEdgeE(vertex, V_node, dt, time, edge_field);

        /** initialize element matrix **/
        tbox::Pointer<tbox::Matrix<double> > ele_mat_N = new tbox::Matrix<double>();
        tbox::Pointer<tbox::Matrix<double> > ele_mat_H = new tbox::Matrix<double>();

        ele_mat_N->resize(n_vertex, n_vertex);
        ele_mat_H->resize(n_vertex, n_vertex);

        for (int l = 0; l < n_vertex; ++l){
            for (int j = 0; j < n_vertex; ++j){
                (*ele_mat_N)(l,j) = 0.0;
                (*ele_mat_H)(l,j) = 0.0;
            }
        }

        /** CV-FEM-SG **/
        tbox::Array<tbox::Array<int> > r_node = ele->get_r_node();
        tbox::Array<tbox::Array<double> > r_integ = getInteg(n_vertex,save_bas,i,ele,bas_integ,vertex,time,dt);

        for(int i_v = 0; i_v < num_ver; ++i_v){
            for(int i_e = 0; i_e < num_edge; ++i_e){
                int n_1 = r_node[i_e][0];
                int n_2 = r_node[i_e][1];

                double T_mid = (T_node[n_2] + T_node[n_1])/2;
                double mun = mat_cell->electronMobility()*pow(m,2)/V/s;
                double mup = mat_cell->holeMobility()*pow(m,2)/V/s;

                double Vt = kb*T_mid/e;
                double delt = -(V_node[n_2] - V_node[n_1])/Vt;

                double cn_1;
                double cn_2;
                double cp_1;
                double cp_2;
                if(DD_type == "SUPG"){
                    /* SUPG */
                    double delt_v = -(V_node[n_2] - V_node[n_1]);

                    double art_dn = mun*Vt;
                    double art_dp = mup*Vt;
                    if(delt == 0.0){
                        art_dn = mun*Vt;
                        art_dp = mup*Vt;
                    } else{
                        art_dn = mun*Vt*(-delt/2)/tanh(-delt/2);
                        art_dp = mup*Vt*(-delt/2)/tanh(-delt/2);
                    }

                    cn_1 = -1*(0.5*mun*delt_v - art_dn)*(r_integ[i_e][i_v]);
                    cn_2 = (0.5*mun*delt_v + art_dn)*(r_integ[i_e][i_v]);
                    cp_1 = (0.5*mup*delt_v + art_dp)*(r_integ[i_e][i_v]);
                    cp_2 = -1*(0.5*mup*delt_v - art_dp)*(r_integ[i_e][i_v]);
                }else if(DD_type=="SG"){
                    /* SG */
                    double tmp_bern_n1 = Var_func->Bern( delt);
                    double tmp_bern_n2 = Var_func->Bern(-delt);
                    double tmp_bern_p1 = Var_func->Bern(-delt);
                    double tmp_bern_p2 = Var_func->Bern( delt);

                    cn_1 = (mun*(kb*T_node[n_1]/e))*tmp_bern_n1*(r_integ[i_e][i_v]);
                    cn_2 = (mun*(kb*T_node[n_2]/e))*tmp_bern_n2*(r_integ[i_e][i_v]);
                    cp_1 = (mup*(kb*T_node[n_1]/e))*tmp_bern_p1*(r_integ[i_e][i_v]);
                    cp_2 = (mup*(kb*T_node[n_2]/e))*tmp_bern_p2*(r_integ[i_e][i_v]);
                }

                (*ele_mat_N)(i_v,n_1) = (*ele_mat_N)(i_v,n_1) - cn_1;
                (*ele_mat_N)(i_v,n_2) = (*ele_mat_N)(i_v,n_2) + cn_2;

                (*ele_mat_H)(i_v,n_1) = (*ele_mat_H)(i_v,n_1) + cp_1;
                (*ele_mat_H)(i_v,n_2) = (*ele_mat_H)(i_v,n_2) - cp_2;
            }
        }


        if(transient_state){
            //if transient state is true, the below will run, else the below will not run
            /** construct stiffness matrix **/
            tbox::Array<double> value(n_vertex);
            for (int i_n = 0; i_n < n_vertex; ++i_n){
                value[i_n] = 1;
            }

            /** building element vector **/
            tbox::Pointer<tbox::Vector<double> > cv_vol = getRhs(n_vertex,save_bas,i,ele,bas_rhs,value,vertex,time,dt);
            for(int i_n=0; i_n<n_vertex; i_n++){
                for(int j=0; j<n_vertex; j++){
                    (*ele_mat_N)(i_n,j) = (*ele_mat_N)(i_n,j)*dt*s;
                    (*ele_mat_H)(i_n,j) = (*ele_mat_H)(i_n,j)*dt*s;
                }
                if(time>0){
                    (*ele_mat_N)(i_n,i_n) = (*ele_mat_N)(i_n,i_n) - (*cv_vol)[i_n];
                    (*ele_mat_H)(i_n,i_n) = (*ele_mat_H)(i_n,i_n) + (*cv_vol)[i_n];
                }
            }
        }

        /** construct system matrix **/
        int row = 0, col = 0;
        for (int l = 0; l < n_vertex; ++l){
            row = mapping[l];
            if((material_list[(*cell_mat)(0,i)].mat_type == "dielectric")&&((*semiinsbc_n)(0,node_n[l]) == 1))
                row = mapping[l]+1;
            for(int j = 0; j < n_vertex; ++j){
                col = mapping[j];
                if((material_list[(*cell_mat)(0,i)].mat_type == "dielectric")&&((*semiinsbc_n)(0,node_n[j]) == 1))
                    col = mapping[j]+1;
                double n_mat_tmp = (*ele_mat_N)(l,j);
                double p_mat_tmp = (*ele_mat_H)(l,j);
                //cout<<l<<" "<<n_mat_tmp<<" "<<p_mat_tmp<<endl;
                n_mat_data->addMatrixValue(row, col, n_mat_tmp);
                h_mat_data->addMatrixValue(row, col, p_mat_tmp);
                #if 0
                if((row==4)&&(col==0)){
                    cout<<material_list[(*cell_mat)(0,i)].mat_type<<endl;
                    cout<<n_mat_tmp<<endl;
                    cout<<""<<endl;
                }
                #endif
            }
        }
    }
    /** assemble matrix **/
    n_mat_data->assemble();
    h_mat_data->assemble();
}

/**
 * @brief NonlinearPoisson::buildCarnFXOnPatch: CV-FEM-SG for electron continuity equation
 * @param patch
 * @param time
 * @param dt
 * @param component_name
 */
void NonlinearPoisson::buildCarnFXOnPatch(hier::Patch<NDIM> & patch,
                            const double  time,
                            const double  dt,
                            const string& component_name){
    /** (形函数, 积分器, 单元)管理器 **/
    tbox::Pointer<ElementManager<NDIM> >
      ele_manager = ElementManager<NDIM>::getManager();
    tbox::Pointer<BaseElement<NDIM> >
      ele = ele_manager->getElement(d_element_type);

    /** Get Local PatchGeometry **/
    tbox::Pointer< hier::PatchGeometry<NDIM> > patch_geo = patch.getPatchGeometry();
    /** Get Local PatchTopology **/
    tbox::Pointer< hier::PatchTopology<NDIM> > patch_top = patch.getPatchTopology();
    /** 取出本地Patch的结点坐标数组 **/
    tbox::Pointer<pdat::NodeData<NDIM, double> >
      node_coord = patch_geo->getNodeCoordinates();
    /** get data of RHS **/
    tbox::Pointer<pdat::VectorData<NDIM,double> >
      n_vec_data = patch.getPatchData(d_carn_fval_id);
    tbox::Pointer<pdat::VectorData<NDIM,double> >
      h_vec_data = patch.getPatchData(d_carh_fval_id);

    /** get node values **/
    tbox::Pointer<pdat::VectorData<NDIM,double> >
      N_soln = patch.getPatchData(d_carn_sol_cur_id);
    tbox::Pointer<pdat::VectorData<NDIM,double> >
      H_soln = patch.getPatchData(d_carh_sol_cur_id);

    /** cell - entity relationship **/
    tbox::Pointer<pdat::CellData<NDIM, int> >
      cell_mat = patch.getPatchData(d_cell_mat_id);
    /** cell - heter node relationship **/
    tbox::Pointer<pdat::CellData<NDIM, int> >
      cell_hnode = patch.getPatchData(d_cell_hnode_id);
    /** heterojunction - node relationship **/
    tbox::Pointer<pdat::NodeData<NDIM, int> >
      heter_node = patch.getPatchData(d_heter_node_id);
    tbox::Pointer<pdat::NodeData<NDIM, int> >
      semiinsbc_node = patch.getPatchData(d_semiinsbc_node_id);

    tbox::Pointer<pdat::NodeData<NDIM, double> >
          carn_sol_1_data = patch.getPatchData(d_carn_sol_1_id);
    tbox::Pointer<pdat::NodeData<NDIM, double> >
          carh_sol_1_data = patch.getPatchData(d_carh_sol_1_id);

    tbox::Pointer<pdat::CellData<NDIM, double> > bas_integ;
    tbox::Pointer<pdat::CellData<NDIM, double> > bas_rhs;
    if(save_bas){
        bas_integ = patch.getPatchData(d_bas_integ_id);
        bas_rhs = patch.getPatchData(d_bas_rhs_id);
    }

    /** get carn plot values **/
    tbox::Pointer<pdat::NodeData<NDIM, double> >
      carn_plot_data = patch.getPatchData(d_carn_plot_id);
    /** get carh plot values **/
    tbox::Pointer<pdat::NodeData<NDIM, double> >
      carh_plot_data = patch.getPatchData(d_carh_plot_id);

    /** 获取本地单元(边, 面)周围结点的索引关系 **/
    tbox::Array<int> can_extent, can_indices;
    patch_top->getCellAdjacencyNodes(can_extent, can_indices);

    int num_cells = patch.getNumberOfCells(1);

    int * dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);
    int * dof_map_semi = d_dof_semi_info->getDOFMapping(patch, hier::EntityUtilities::NODE);

    /** material list **/
    tbox::Pointer<MaterialManager<NDIM> >
            mat_manager = MaterialManager<NDIM>::getManager();

    tbox::Pointer<Variables_Functions> Var_func = new Variables_Functions();

    n_vec_data->fillAll(0.0);
    h_vec_data->fillAll(0.0);
    /** Traversing Unit **/
    for(int i = 0; i < num_cells; ++i){
        if(!((material_list[(*cell_mat)(0, i)].mat_type == "semiconductor")||(material_list[(*cell_mat)(0, i)].mat_type == "dielectric"))) continue;

        int n_vertex = can_extent[i+1] - can_extent[i];
        tbox::Array<hier::DoubleVector<NDIM> > vertex(n_vertex);
        tbox::Array<int> node_n(n_vertex);
        tbox::Array<int> mapping(n_vertex);
        tbox::Array<int> mapping_semi(n_vertex);
        tbox::Array<double> T_node(n_vertex);
        tbox::Array<double> N_node(n_vertex);
        tbox::Array<double> H_node(n_vertex);
        tbox::Array<double> n_RGrate(n_vertex);
        tbox::Array<double> h_RGrate(n_vertex);
        tbox::Array<double> carn_node(n_vertex);
        tbox::Array<double> carh_node(n_vertex);

        /** inform for node **/
        for (int i1 = 0, j = can_extent[i]; i1 < n_vertex; ++i1, ++j) {
            mapping[i1] = dof_map[can_indices[j]];
            mapping_semi[i1] = dof_map_semi[can_indices[j]];
            for(int k = 0; k< NDIM; ++k){
              vertex[i1][k] = (*node_coord)(k, can_indices[j])*m;
            }
            node_n[i1] = can_indices[j];
            carn_node[i1] = (*carn_plot_data)(0,can_indices[j])*pow(m,-3);
            carh_node[i1] = (*carh_plot_data)(0,can_indices[j])*pow(m,-3);
            if((material_list[(*cell_mat)(0, i)].mat_type == "dielectric")&&((*semiinsbc_node)(0,node_n[i1])==1)){
                carn_node[i1] = (*carn_sol_1_data)(1,node_n[i1]);
                carh_node[i1] = (*carh_sol_1_data)(1,node_n[i1]);
            }
        }

        /** get material for cell **/
        tbox::Pointer<BaseMaterial<NDIM> > mat_cell =
                mat_manager->getMaterial(material_list[(*cell_mat)(0,i)].mat_name);

        for (int i_n = 0; i_n < n_vertex; ++i_n){
            T_node[i_n]=T;
            N_node[i_n]=(*N_soln)(mapping[i_n]);
            H_node[i_n]=(*H_soln)(mapping[i_n]);
            if((material_list[(*cell_mat)(0, i)].mat_type == "dielectric")&&((*semiinsbc_node)(0,node_n[i_n])==1)){
                N_node[i_n]=(*N_soln)(mapping[i_n]+1);
                H_node[i_n]=(*H_soln)(mapping[i_n]+1);
            }
        }

        double T_cell = T;

        /** semi-variables for RG rate **/
        double e_trap = mat_cell->electronTrap(T_cell);
        bool fermi_st = false;

        double U_radiation = 0;
        if(material_list[(*cell_mat)(0, i)].mat_type == "dielectric")  {
            double Rd = 0;
            U_radiation = 0.01*8.1e12*Rd*pow(cm,-3)/s; //U=Y*g0*Rd;
        }

        /** two ways to calculate the recombination & generation rate
         1. get the cell carrier density the calculate the cell rate
         2. get the node rate and then calculate the element rate **/
        for (int i_n = 0; i_n < n_vertex; ++i_n){
            double ni_tmp = mat_cell->intrinsicDensityNi(T_amb)*pow(m,-3);
            double R = 0;
            if(material_list[(*cell_mat)(0, i)].mat_type == "semiconductor"){
                //R = Var_func->recombGenerationSRH(N_node[i_n], H_node[i_n], 0.0, 0.0, T_node[i_n], kb, e, ni_tmp, e_trap, fermi_st);
            }
            n_RGrate[i_n] = R - U_radiation;
            h_RGrate[i_n] = U_radiation - R;
        }

        /** building element vector **/
        tbox::Pointer<tbox::Vector<double> > n_ele_vec = getRhs(n_vertex,save_bas,i,ele,bas_rhs,n_RGrate,vertex,time,dt);
        tbox::Pointer<tbox::Vector<double> > h_ele_vec = getRhs(n_vertex,save_bas,i,ele,bas_rhs,h_RGrate,vertex,time,dt);

        if(transient_state){
            //if transient state is true, the below will run, else the below will not run
            tbox::Array<double> value(n_vertex);
            for (int i_n = 0; i_n < n_vertex; ++i_n){
                value[i_n] = 1;
            }
            tbox::Pointer<tbox::Vector<double> > cv_vol = getRhs(n_vertex,save_bas,i,ele,bas_rhs,value,vertex,time,dt);

            for(int i_n = 0; i_n < n_vertex; ++i_n){
                (*n_ele_vec)[i_n] = dt * s * (*n_ele_vec)[i_n] - (*cv_vol)[i_n] * carn_node[i_n];
                (*h_ele_vec)[i_n] = dt * s * (*h_ele_vec)[i_n] - (*cv_vol)[i_n] * carh_node[i_n];
            }
        }

        for(int i2= 0; i2 < n_vertex; ++i2){
            //double n1 = (*n_ele_vec)[i2]; double h1 = (*h_ele_vec)[i2];
            int index = mapping[i2];
            if((material_list[(*cell_mat)(0, i)].mat_type == "dielectric")&&((*semiinsbc_node)(0,node_n[i2])==1))
                index = mapping[i2]+1;
            n_vec_data->addVectorValue(index, (*n_ele_vec)[i2]);
            h_vec_data->addVectorValue(index, (*h_ele_vec)[i2]);
        }
    }
}

/**
 * @brief NonlinearPoisson::setCarnPhysicalBC: Matrix
 * @param patch
 * @param time
 * @param dt
 * @param component_name
 */
void NonlinearPoisson::setCarnPhysicalBC(hier::Patch<NDIM> & patch,
                            const double  time,
                            const double  dt,
                            const string& component_name){
    /** 取出本地PatchGeometry **/
    tbox::Pointer< hier::PatchGeometry<NDIM> > patch_geo = patch.getPatchGeometry();
    /** Get Local PatchTopology **/
    tbox::Pointer< hier::PatchTopology<NDIM> > patch_top = patch.getPatchTopology();
    /** 获取网格片矩阵对象 **/
    tbox::Pointer<pdat::CSRMatrixData<NDIM,double> >
      n_mat_data = patch.getPatchData(d_carn_mat_id);
    /** 获取网格片矩阵对象 **/
    tbox::Pointer<pdat::CSRMatrixData<NDIM,double> >
      h_mat_data = patch.getPatchData(d_carh_mat_id);

    /** face - entity relationship **/
    tbox::Pointer<pdat::FaceData<NDIM, int> >
      face_mat = patch.getPatchData(d_face_mat_id);

    /** heterojunction - node relationship **/
    tbox::Pointer<pdat::NodeData<NDIM, int> >
      heter_node = patch.getPatchData(d_heter_node_id);

    /** 取出本地Patch的结点坐标数组 **/
    tbox::Pointer<pdat::NodeData<NDIM, double> >
      node_coord = patch_geo->getNodeCoordinates();
    /** 面——单元邻接关系数据 **/
    tbox::Array<int>face_cell_ext,face_cell_idx;
    patch_top->getFaceAdjacencyCells(face_cell_ext,face_cell_idx);
    /** 单元节点邻接关系数据 **/
    tbox::Array<int>face_node_ext,face_node_idx;
    patch_top->getFaceAdjacencyNodes(face_node_ext,face_node_idx);

    /** obtain the current solution vector **/
    tbox::Pointer<pdat::VectorData<NDIM,double> >
      n_fval_data = patch.getPatchData(d_carn_fval_id);
    tbox::Pointer<pdat::VectorData<NDIM,double> >
      h_fval_data = patch.getPatchData(d_carh_fval_id);

    tbox::Pointer<pdat::VectorData<NDIM,double> >
      V_soln = patch.getPatchData(d_poisson_sol_cur_id);
    tbox::Pointer<pdat::NodeData<NDIM, double> >
      dop = patch.getPatchData(d_dopant_id);

    /** quasi fermi level **/
    tbox::Pointer<pdat::VectorData<NDIM,double> >
      Fn_soln = patch.getPatchData(d_fermin_sol_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
      Fh_soln = patch.getPatchData(d_fermih_sol_id);

    tbox::Pointer<pdat::NodeData<NDIM, int> >
      node_mat = patch.getPatchData(d_node_mat_id);

    int num_nodes = patch.getNumberOfNodes(1);
    //int num_faces = patch.getNumberOfEntities(hier::EntityUtilities::FACE,1);

    /** material list **/
    tbox::Pointer<MaterialManager<NDIM> >
            mat_manager = MaterialManager<NDIM>::getManager();

    tbox::Pointer<Variables_Functions> Var_func = new Variables_Functions();

    /** 自由度信息中的映射信息 **/
    int * dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);
    int * dof_map_semi = d_dof_semi_info->getDOFMapping(patch, hier::EntityUtilities::NODE);

    /** 取出矩阵向量的核心数据结构的指针 **/
    int * n_row_start = n_mat_data->getRowStartPointer();
    int * n_col_idx   = n_mat_data->getColumnIndicesPointer();
    double * n_mat_val    = n_mat_data->getValuePointer();
    int * h_row_start = h_mat_data->getRowStartPointer();
    int * h_col_idx   = h_mat_data->getColumnIndicesPointer();
    double * h_mat_val    = h_mat_data->getValuePointer();

    /** 取出矩阵向量的核心数据结构的指针 **/
    double *n_fval_val = n_fval_data->getPointer();
    double *h_fval_val = h_fval_data->getPointer();

    /** setup the source target **/
    tbox::Pointer<bound_condition> bound_source = new bound_condition();

    double n_tmp_val, h_tmp_val;

    /** boundary condition source **/
    int num_bc = EBoundassignment.getSize();
    tbox::Array<double> bound_coe;
    bound_coe.resizeArray(9);

    for(int i_bc=0; i_bc<num_bc; ++i_bc){
        if (EBoundassignment[i_bc].bc_type == "ohmic"){               // string.Equal(str1, str2);
            if (patch_geo->hasEntitySet(EBoundassignment[i_bc].bc_no, hier::EntityUtilities::NODE)){

                const tbox::Array<int> &entity_idx = patch_geo->getEntityIndicesInSet(EBoundassignment[i_bc].bc_no,
                                                     hier::EntityUtilities::NODE,
                                                     num_nodes);

                double app_volt = bound_source->source_assign(EBoundassignment[i_bc].sc_type,
                                                              time,
                                                              dt,
                                                              EBoundassignment[i_bc].sc_coe);
                string contact_type = "Ohmic_contact";
                int size = entity_idx.getSize();
                for (int i = 0; i < size; ++i) {
                    int index = dof_map[entity_idx[i]];
                    int index_s = dof_map_semi[entity_idx[i]];
                    /** setup material target **/
                    int dom_no;
                    dom_no = 0;
                    for(int i_m = 0; i_m < material_list.size(); ++i_m){
                        if(material_list[i_m].mat_type == "semiconductor"){
                            if((*node_mat)(i_m,entity_idx[i])==1){
                                dom_no = i_m;
                            }
                        }
                    }
                    tbox::Pointer<BaseMaterial<NDIM> > mat_dom =
                            mat_manager->getMaterial(material_list[dom_no].mat_name);
                    bound_coe[0] = app_volt;
                    bound_coe[1] = T;
                    bound_coe[2] = mat_dom->effectiveAffinity(T_amb);
                    bound_coe[3] = mat_dom->bandGap(T_amb);
                    bound_coe[4] = (*Fn_soln)(index);
                    bound_coe[5] = (*Fh_soln)(index);
                    bound_coe[6] = (*V_soln)(index_s);
                    bound_coe[7] = (*dop)(0,entity_idx[i]);
                    bound_coe[8] = mat_dom->intrinsicDensityNi(T_amb)*pow(m,-3);
//                    bound_coe[8] = 3.887199240115326E15*pow(m,-3); //这里采用COMSOL的有效本征载流子浓度
//                    double Ni = mat_dom->intrinsicDensityNi(T_amb);
//                    double Ni_eff = mat_dom->intrinsicDensityNi(T_amb)*pow(m,-3);
                    n_tmp_val = bound_source->carn_bound_value(contact_type, car_bound_mode, time, dt, kb, e, bound_coe);
                    n_fval_val[index] = n_tmp_val;
                    h_tmp_val = bound_source->carh_bound_value(contact_type, car_bound_mode, time, dt, kb, e, bound_coe);
                    h_fval_val[index] = h_tmp_val;
                    for (int j = n_row_start[index]; j < n_row_start[index+1]; ++j){
                        if (n_col_idx[j] == index){
                            n_mat_val[j] = 1.0;   /**< 对角线元素 **/
                        }
                        else {
                            n_mat_val[j] = 0.0;   /**< 非对角线元素 **/
                        }
                    }
                    for (int j = h_row_start[index]; j < h_row_start[index+1]; ++j){
                        if (h_col_idx[j] == index){
                            h_mat_val[j] = 1.0;   /**< 对角线元素 */
                        }
                        else {
                            h_mat_val[j] = 0.0;   /**< 非对角线元素 */
                        }
                    }
                }
            }
        } else if (EBoundassignment[i_bc].bc_type == "gate") {
            if (patch_geo->hasEntitySet(EBoundassignment[i_bc].bc_no, hier::EntityUtilities::NODE)){

                const tbox::Array<int> &entity_idx = patch_geo->getEntityIndicesInSet(EBoundassignment[i_bc].bc_no,
                                                     hier::EntityUtilities::NODE,
                                                     num_nodes);
                string contact_type = "Gate_contact";
                int size = entity_idx.getSize();
                for (int i = 0; i < size; ++i) {
                    int index = dof_map[entity_idx[i]];
                    n_fval_val[index] = 1e-10*pow(m,-3);
                    h_fval_val[index] = 1e-10*pow(m,-3);
                    for (int j = n_row_start[index]; j < n_row_start[index+1]; ++j){
                        if (n_col_idx[j] == index){
                            n_mat_val[j] = 1.0;   /**< 对角线元素 **/
                        }
                        else {
                            n_mat_val[j] = 0.0;   /**< 非对角线元素 **/
                        }
                    }
                    for (int j = h_row_start[index]; j < h_row_start[index+1]; ++j){
                        if (h_col_idx[j] == index){
                            h_mat_val[j] = 1.0;   /**< 对角线元素 **/
                        }
                        else {
                            h_mat_val[j] = 0.0;   /**< 非对角线元素 **/
                        }
                    }
                }
            }
        }
    }
}

/**
 * @brief NonlinearPoisson::solutionScal
 * @param patch
 * @param time
 * @param dt
 * @param component_name
 */
void NonlinearPoisson::solutionScal(hier::Patch<NDIM> & patch,
                                      const double  time,
                                      const double  dt,
                                      const string& component_name){

    /** scaled potential variables **/
    tbox::Pointer<pdat::VectorData<NDIM, double> >
      V_sol_scaled = patch.getPatchData(d_poisson_sol_cur_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
      Fn_sol_scaled = patch.getPatchData(d_fermin_sol_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
      Fh_sol_scaled = patch.getPatchData(d_fermih_sol_id);
    /** anti-scaled potential variables **/
    tbox::Pointer<pdat::VectorData<NDIM, double> >
      V_sol_unscaled = patch.getPatchData(d_poisson_unscal_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
      Fn_sol_unscaled = patch.getPatchData(d_fermin_unscal_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
      Fh_sol_unscaled = patch.getPatchData(d_fermih_unscal_id);

    tbox::Pointer<pdat::NodeData<NDIM, int> >
      heter_node = patch.getPatchData(d_heter_node_id);
    tbox::Pointer<pdat::NodeData<NDIM, int> >
      semiinsbc_node = patch.getPatchData(d_semiinsbc_node_id);

    int * dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);
    int * dof_map_semi = d_dof_semi_info->getDOFMapping(patch, hier::EntityUtilities::NODE);

    int num_node = patch.getNumberOfNodes(1);

    for (int i_n = 0; i_n < num_node; ++i_n){
        int mapping = dof_map[i_n];
        int mapping_semi = dof_map_semi[i_n];
        (*V_sol_scaled)(mapping_semi) = (*V_sol_unscaled)(mapping);
        (*Fn_sol_scaled)(mapping) = (*Fn_sol_unscaled)(mapping);
        (*Fh_sol_scaled)(mapping) = (*Fh_sol_unscaled)(mapping);
        if((*semiinsbc_node)(0,i_n)==1){
            (*V_sol_scaled)(mapping_semi+1) = (*V_sol_unscaled)(mapping+1);
            (*Fn_sol_scaled)(mapping+1) = (*Fn_sol_unscaled)(mapping+1);
            (*Fh_sol_scaled)(mapping+1) = (*Fh_sol_unscaled)(mapping+1);
        }
    }
}

/**
 * @brief NonlinearPoisson::solutionAntiScal
 * @param patch
 * @param time
 * @param dt
 * @param component_name
 */
void NonlinearPoisson::solutionAntiScal(hier::Patch<NDIM> & patch,
                                      const double  time,
                                      const double  dt,
                                      const string& component_name){

    /** scaled potential variables **/
    tbox::Pointer<pdat::VectorData<NDIM, double> >
      V_sol_scaled = patch.getPatchData(d_poisson_sol_cur_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
      Fn_sol_scaled = patch.getPatchData(d_fermin_sol_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
      Fh_sol_scaled = patch.getPatchData(d_fermih_sol_id);
    /** anti-scaled potential variables **/
    tbox::Pointer<pdat::VectorData<NDIM, double> >
      V_sol_unscaled = patch.getPatchData(d_poisson_unscal_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
      Fn_sol_unscaled = patch.getPatchData(d_fermin_unscal_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
      Fh_sol_unscaled = patch.getPatchData(d_fermih_unscal_id);

    tbox::Pointer<pdat::NodeData<NDIM, int> >
      heter_node = patch.getPatchData(d_heter_node_id);
    tbox::Pointer<pdat::NodeData<NDIM, int> >
      semiinsbc_node = patch.getPatchData(d_semiinsbc_node_id);

    int * dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);
    int * dof_map_semi = d_dof_semi_info->getDOFMapping(patch, hier::EntityUtilities::NODE);

    int num_node = patch.getNumberOfNodes(1);

    for (int i_n = 0; i_n < num_node; ++i_n){
        int mapping = dof_map[i_n];
        int mapping_semi = dof_map_semi[i_n];
        (*V_sol_unscaled)(mapping) = (*V_sol_scaled)(mapping_semi);
        (*Fn_sol_unscaled)(mapping) = (*Fn_sol_scaled)(mapping);
        (*Fh_sol_unscaled)(mapping) = (*Fh_sol_scaled)(mapping);
        if((*semiinsbc_node)(0,i_n)==1){
            (*V_sol_unscaled)(mapping+1) = (*V_sol_scaled)(mapping_semi+1);
            (*Fn_sol_unscaled)(mapping+1) = (*Fn_sol_scaled)(mapping+1);
            (*Fh_sol_unscaled)(mapping+1) = (*Fh_sol_scaled)(mapping+1);
        }
    }
}

/**
 * @brief NonlinearPoisson::carrierAntiScal
 * @param patch
 * @param time
 * @param dt
 * @param component_name
 */
void NonlinearPoisson::carrierAntiScal(hier::Patch<NDIM> & patch,
                                      const double  time,
                                      const double  dt,
                                      const string& component_name){

    /** scaled solution of carrier density **/
    tbox::Pointer<pdat::VectorData<NDIM, double> >
      N_sol_scaled = patch.getPatchData(d_carn_sol_cur_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
      H_sol_scaled = patch.getPatchData(d_carh_sol_cur_id);
    /** anti-scaled solution of carrier density **/
    tbox::Pointer<pdat::VectorData<NDIM, double> >
      N_sol_unscaled = patch.getPatchData(d_carn_unscal_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
      H_sol_unscaled = patch.getPatchData(d_carh_unscal_id);

    tbox::Pointer<pdat::NodeData<NDIM, int> >
      heter_node = patch.getPatchData(d_heter_node_id);
    tbox::Pointer<pdat::NodeData<NDIM, int> >
      semiinsbc_node = patch.getPatchData(d_semiinsbc_node_id);

    int * dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);

    int num_node = patch.getNumberOfNodes(1);
    for (int i_n = 0; i_n < num_node; ++i_n){
        int mapping = dof_map[i_n];
        (*N_sol_unscaled)(mapping) = (*N_sol_scaled)(mapping);
        (*H_sol_unscaled)(mapping) = (*H_sol_scaled)(mapping);
        if((*heter_node)(0,i_n) == 1){
            (*N_sol_unscaled)(mapping + 1) = (*N_sol_scaled)(mapping + 1);
            (*H_sol_unscaled)(mapping + 1) = (*H_sol_scaled)(mapping + 1);
        }
        if((*semiinsbc_node)(0,i_n) == 1){
            (*N_sol_unscaled)(mapping + 1) = (*N_sol_scaled)(mapping + 1);
            (*H_sol_unscaled)(mapping + 1) = (*H_sol_scaled)(mapping + 1);
        }
    }
}


/**
 * @brief NonlinearPoisson::buildNewtonPFXOnPatch: Newton-Raphson
 * @param patch
 * @param time
 * @param dt
 * @param component_name
 */
void NonlinearPoisson::buildNewtonPFXOnPatch(hier::Patch<NDIM> & patch,
                                      const double  time,
                                      const double  dt,
                                      const string& component_name)
{
    /** (形函数, 积分器, 单元)管理器 **/
    tbox::Pointer<ElementManager<NDIM> >
            ele_manager = ElementManager<NDIM>::getManager();
    tbox::Pointer<BaseElement<NDIM> > ele = ele_manager->getElement(d_element_type);

    /** 取出本地PatchGeometry **/
    tbox::Pointer< hier::PatchGeometry<NDIM> > patch_geo = patch.getPatchGeometry();
    /** 取出本地PatchTopology **/
    tbox::Pointer< hier::PatchTopology<NDIM> > patch_top = patch.getPatchTopology();
    /** 取出本地Patch的结点坐标数组 **/
    tbox::Pointer<pdat::NodeData<NDIM, double> >
            node_coord = patch_geo->getNodeCoordinates();
    /** right side vector **/
    tbox::Pointer<pdat::VectorData<NDIM,double> >
      vec_data = patch.getPatchData(d_poisson_newton_fval_id);

    /** get node values **/
    tbox::Pointer<pdat::VectorData<NDIM, double> >
            V_soln = patch.getPatchData(d_poisson_sol_cur_id);
    tbox::Pointer<pdat::NodeData<NDIM, double> >
            dop = patch.getPatchData(d_dopant_id);
    /** cell - entity relationship **/
    tbox::Pointer<pdat::CellData<NDIM, int> >
            cell_mat = patch.getPatchData(d_cell_mat_id);
    /** cell - entity relationship **/
    tbox::Pointer<pdat::CellData<NDIM, int> >
            cell_hnode = patch.getPatchData(d_cell_hnode_id);

    /** node - entity relationship **/
    tbox::Pointer<pdat::NodeData<NDIM, int> >
            heter_node = patch.getPatchData(d_heter_node_id);

    /** 获取本地单元(边, 面)周围结点的索引关系 **/
    tbox::Array<int> can_extent, can_indices;
    patch_top->getCellAdjacencyNodes(can_extent, can_indices);

    int num_cells = patch.getNumberOfCells(1);
    //int num_nodes = patch.getNumberOfNodes(1);

    int * dof_map_semi = d_dof_semi_info->getDOFMapping(patch, hier::EntityUtilities::NODE);
    /** material list **/
    tbox::Pointer<MaterialManager<NDIM> >
            mat_manager = MaterialManager<NDIM>::getManager();

    /** setup target for class 'Variables_Functions' **/
    tbox::Pointer<Variables_Functions> Var_func = new Variables_Functions();

    vec_data->fillAll(0.0);   // reset the vector data

    /** Trasversing Unit **/
    for(int i = 0; i < num_cells; ++i){
        int n_vertex = can_extent[i+1] - can_extent[i];
        tbox::Array<hier::DoubleVector<NDIM> > vertex(n_vertex);
        tbox::Array<double> dop_node(n_vertex);
        tbox::Array<int> e_node(n_vertex);
        tbox::Array<int> mapping_semi(n_vertex);
        tbox::Array<double> T_node(n_vertex);
        tbox::Array<double> V_node(n_vertex);
        /** carrier density **/
        tbox::Array<double> n_node(n_vertex);
        tbox::Array<double> p_node(n_vertex);
        tbox::Array<double> value(n_vertex);
        /** value of fV,fn,fp **/
        tbox::Array<double> fV_node(n_vertex);
        tbox::Array<double> fn_node(n_vertex);
        tbox::Array<double> fp_node(n_vertex);

        for (int i1 = 0, j = can_extent[i]; i1 < n_vertex; ++i1, ++j) {
            mapping_semi[i1] = dof_map_semi[can_indices[j]];
            for(int k = 0; k< NDIM; ++k){
                vertex[i1][k] = (*node_coord)(k, can_indices[j])*m;
            }
            dop_node[i1]=(*dop)(0,can_indices[j]);
            e_node[i1]=can_indices[j];
        }

        /** get material for cell **/
        tbox::Pointer<BaseMaterial<NDIM> > mat_cell =
                mat_manager->getMaterial(material_list[(*cell_mat)(0,i)].mat_name);

//        double T_cell = T;
//        double cell_epsr = mat_cell->electroPermittivity(T_cell);
//        double e_trap = mat_cell->electronTrap(T_cell);
//        double mun = mat_cell->electronMobility()*pow(m,2)/V/s;
//        double mup = mat_cell->holeMobility()*pow(m,2)/V/s;
//        double Vt = kb*T/e;

        tbox::Array<tbox::Array<double> > r_integ = ele->buildCVFEMInteg(vertex, dt, time);

        tbox::Array<tbox::Array<int> > r_node = ele->get_r_node();

        /** initialize element matrix **/
        tbox::Pointer<tbox::Matrix<double> > ele_mat_E = new tbox::Matrix<double>();

        ele_mat_E->resize(n_vertex, n_vertex);

        for (int l = 0; l < n_vertex; ++l){
            for (int j = 0; j < n_vertex; ++j){
                (*ele_mat_E)(l,j) = 0.0;
            }
        }

        for(int i_v = 0; i_v < cell_point; ++i_v){
            for(int i_e = 0; i_e < cell_edge; ++i_e){
                int n_1 = r_node[i_e][0];
                int n_2 = r_node[i_e][1];

                double cn_1 = (r_integ[i_e][i_v]);
                double cn_2 = (r_integ[i_e][i_v]);

                (*ele_mat_E)(i_v,n_1) = (*ele_mat_E)(i_v,n_1) - cn_1;
                (*ele_mat_E)(i_v,n_2) = (*ele_mat_E)(i_v,n_2) + cn_2;
            }
        }

        /** element volume **/
        tbox::Pointer<tbox::Vector<double> > ele_volume = new tbox::Vector<double>();
        ele_volume->resize(n_vertex);

        for (int i_n = 0; i_n < n_vertex; ++i_n){
            value[i_n] = 1;
        }

        /** building element volume **/
        ele->buildCVFEMRhs(vertex, value, dt, time, ele_volume);

        for(int l=0; l<n_vertex; ++l){
            V_node[l] = (*V_soln)(mapping_semi[l]);
            n_node[l] = (*V_soln)(mapping_semi[l]+1);
            p_node[l] = (*V_soln)(mapping_semi[l]+2);
        }

        for(int l=0; l<n_vertex; ++l){
        }


        #if 0
        for(int l=0;l<n_vertex;l++){
            cout<<"(*ele_vec)[l]:"<<(*ele_vec)[l]<<endl;
        }
        #endif

        #if 0
        for(int i_v = 0; i_v < cell_point; ++i_v){
            for(int i_n = 0; i_n < cell_point; ++i_n){
                cout<<"(*ele_mat_E)(i_v,i_n):"<<i_v<<":"<<i_n<<":"<<(*ele_mat_E)(i_v,i_n)<<endl;
            }
        }
        #endif

        /** construct fval vector **/
        for(int i2= 0; i2 < n_vertex; ++i2){
            vec_data->addVectorValue(mapping_semi[i2], -fV_node[i2]);
            vec_data->addVectorValue(mapping_semi[i2]+1, -fn_node[i2]);
            vec_data->addVectorValue(mapping_semi[i2]+2, -fp_node[i2]);
        }
    }
}

/** **********************************************************************
 * @brief Newton Vogamma method
 ********************************************************************* **/

/**
 * @brief NonlinearPoisson::buildNewtonVogammaMatOnPatch: Newton-Raphson
 * @param patch
 * @param time
 * @param dt
 * @param component_name
 */
void NonlinearPoisson::buildNewtonVogammaMatOnPatch(hier::Patch<NDIM> & patch,
                                          const double  time,
                                          const double  dt,
                                          const string& component_name)
{
    /** manager for shapefuction, integrator & element **/
    tbox::Pointer<ElementManager<NDIM> >
            ele_manager = ElementManager<NDIM>::getManager();
    tbox::Pointer<BaseElement<NDIM> >
            ele = ele_manager->getElement(d_element_type);

    /** gPatchGeometry **/
    tbox::Pointer< hier::PatchGeometry<NDIM> >
            patch_geo = patch.getPatchGeometry();
    /** get PatchTopology **/
    tbox::Pointer< hier::PatchTopology<NDIM> >
            patch_top = patch.getPatchTopology();

    /** get Patch node coord **/
    tbox::Pointer<pdat::NodeData<NDIM, double> >
            node_coord = patch_geo->getNodeCoordinates();
    /** stiff matrix on patch **/
    tbox::Pointer<pdat::CSRMatrixData<NDIM,double> >
            mat_data = patch.getPatchData(d_Vogamma_newton_mat_id);
    /** right side vector **/
    tbox::Pointer<pdat::VectorData<NDIM,double> >
            vec_data = patch.getPatchData(d_Vogamma_newton_fval_id);

    tbox::Pointer<pdat::VectorData<NDIM, double> >
            Vogamma_soln = patch.getPatchData(d_Vogamma_sol_cur_id);

    /** cell - entity relationship **/
    tbox::Pointer<pdat::CellData<NDIM, int> >
            cell_mat = patch.getPatchData(d_cell_mat_id);
    /** cell - heter node relationship **/
    tbox::Pointer<pdat::CellData<NDIM, int> >
            cell_hnode = patch.getPatchData(d_cell_hnode_id);
    /** node - entity relationship **/
    tbox::Pointer<pdat::NodeData<NDIM, int> >
            heter_node = patch.getPatchData(d_heter_node_id);

    /** get carn plot values **/
    tbox::Pointer<pdat::NodeData<NDIM, double> >
      Vogamma_plot_data = patch.getPatchData(d_Vogamma_plot_id);
    tbox::Pointer<pdat::NodeData<NDIM, double> >
      Vogammaplus_plot_data = patch.getPatchData(d_Vogammaplus_plot_id);
    tbox::Pointer<pdat::NodeData<NDIM, double> >
      VogammaHplus_plot_data = patch.getPatchData(d_VogammaHplus_plot_id);
    tbox::Pointer<pdat::NodeData<NDIM, double> >
      VogammaH2plus_plot_data = patch.getPatchData(d_VogammaH2plus_plot_id);
    tbox::Pointer<pdat::NodeData<NDIM, double> >
      VogammaH_plot_data = patch.getPatchData(d_VogammaH_plot_id);
    tbox::Pointer<pdat::NodeData<NDIM, double> >
      VogammaH2_plot_data = patch.getPatchData(d_VogammaH2_plot_id);

    /** 获取本地单元(边, 面)周围结点的索引关系 **/
    tbox::Array<int> can_extent, can_indices;
    patch_top->getCellAdjacencyNodes(can_extent, can_indices);

    int num_cells = patch.getNumberOfCells(1);
    int num_nodes = patch.getNumberOfNodes(1);

    int * dof_map_Vogamma = d_dof_Vogamma_info->getDOFMapping(patch, hier::EntityUtilities::NODE);

    /** material list **/
    tbox::Pointer<MaterialManager<NDIM> >
            mat_manager = MaterialManager<NDIM>::getManager();

    /** variable function **/
    tbox::Pointer<Variables_Functions> Var_func = new Variables_Functions();

    adtl::AutoDScalar::numdir = 6*8;

    /** Traversing Unit **/
    for(int i = 0; i < num_cells; ++i){

        /** get node information **/
        int n_vertex = can_extent[i+1] - can_extent[i];
        tbox::Array<hier::DoubleVector<NDIM> > vertex(n_vertex);
        tbox::Array<int> node_n(n_vertex);
        tbox::Array<int> mapping_Vogamma(n_vertex);
        tbox::Array<double> value(n_vertex);
        tbox::Array<double> Vogamma_plot(n_vertex);
        tbox::Array<double> Vogammaplus_plot(n_vertex);
        tbox::Array<double> VogammaHplus_plot(n_vertex);
        tbox::Array<double> VogammaH2plus_plot(n_vertex);
        tbox::Array<double> VogammaH_plot(n_vertex);
        tbox::Array<double> VogammaH2_plot(n_vertex);
        adtl::AutoDScalar Vogamma_node[n_vertex];
        adtl::AutoDScalar Vogammaplus_node[n_vertex];
        adtl::AutoDScalar VogammaHplus_node[n_vertex];
        adtl::AutoDScalar VogammaH2plus_node[n_vertex];
        adtl::AutoDScalar VogammaH_node[n_vertex];
        adtl::AutoDScalar VogammaH2_node[n_vertex];
        adtl::AutoDScalar fVogamma_node[n_vertex];
        adtl::AutoDScalar fVogammaplus_node[n_vertex];
        adtl::AutoDScalar fVogammaHplus_node[n_vertex];
        adtl::AutoDScalar fVogammaH2plus_node[n_vertex];
        adtl::AutoDScalar fVogammaH_node[n_vertex];
        adtl::AutoDScalar fVogammaH2_node[n_vertex];

        for (int i1 = 0, j = can_extent[i]; i1 < n_vertex; ++i1, ++j) {
            mapping_Vogamma[i1] = dof_map_Vogamma[can_indices[j]];
            for(int k = 0; k< NDIM; ++k){
                vertex[i1][k] = (*node_coord)(k, can_indices[j])*m;
            }
            node_n[i1] = can_indices[j];
            Vogamma_plot[i1] = (*Vogamma_plot_data)(0,can_indices[j])*pow(m,-3);
            Vogammaplus_plot[i1] = (*Vogammaplus_plot_data)(0,can_indices[j])*pow(m,-3);
            VogammaHplus_plot[i1] = (*VogammaHplus_plot_data)(0,can_indices[j])*pow(m,-3);
            VogammaH2plus_plot[i1] = (*VogammaH2plus_plot_data)(0,can_indices[j])*pow(m,-3);
            VogammaH_plot[i1] = (*VogammaH_plot_data)(0,can_indices[j])*pow(m,-3);
            VogammaH2_plot[i1] = (*VogammaH2_plot_data)(0,can_indices[j])*pow(m,-3);
        }

        #if 0
        cout<<"cell:"<<i<<endl;
        for(int i1=0;i1<8;i1++){
            cout<<vertex[i1][0]<<" "<<vertex[i1][1]<<" "<<vertex[i1][2]<<"\n";
        }
        #endif

        /** get material for cell **/
        tbox::Pointer<BaseMaterial<NDIM> > mat_cell =
                mat_manager->getMaterial(material_list[(*cell_mat)(0,i)].mat_name);
        //cout<<"mat_name:"<<material_list[(*cell_mat)(0,i)].mat_name<<endl;

        /** element volume **/
        tbox::Pointer<tbox::Vector<double> > ele_volume = new tbox::Vector<double>();
        ele_volume->resize(n_vertex);

        for (int i_n = 0; i_n < n_vertex; ++i_n){
            value[i_n] = 1;
        }

        /** building element volume **/
        ele->buildCVFEMRhs(vertex, value, dt, time, ele_volume);

        #if 0
        for(int i_v = 0; i_v < cell_point; ++i_v){
            for(int i_n = 0; i_n < cell_point; ++i_n){
                cout<<"(*ele_mat_E)(i_v,i_n):"<<i_v<<":"<<i_n<<":"<<(*ele_mat_E)(i_v,i_n)<<endl;
            }
        }
        #endif

        for (int i1 = 0; i1 < n_vertex; ++i1) {
            Vogamma_node[i1] = (*Vogamma_soln)(mapping_Vogamma[i1]); Vogamma_node[i1].setADValue(i1*6,1.0);
            Vogammaplus_node[i1] = (*Vogamma_soln)(mapping_Vogamma[i1]+1); Vogammaplus_node[i1].setADValue(i1*6+1,1.0);
            VogammaHplus_node[i1] = (*Vogamma_soln)(mapping_Vogamma[i1]+2); VogammaHplus_node[i1].setADValue(i1*6+2,1.0);
            VogammaH2plus_node[i1] = (*Vogamma_soln)(mapping_Vogamma[i1]+2); VogammaH2plus_node[i1].setADValue(i1*6+3,1.0);
            VogammaH_node[i1] = (*Vogamma_soln)(mapping_Vogamma[i1]+2); VogammaH_node[i1].setADValue(i1*6+4,1.0);
            VogammaH2_node[i1] = (*Vogamma_soln)(mapping_Vogamma[i1]+2); VogammaH2_node[i1].setADValue(i1*6+5,1.0);
        }

        for (int i1 = 0; i1 < n_vertex; ++i1) {
            /* get node data*/
            adtl::AutoDScalar Vogamma=Vogamma_node[i1];
            adtl::AutoDScalar Vogammaplus=Vogammaplus_node[i1];
            adtl::AutoDScalar VogammaHplus=VogammaHplus_node[i1];
            adtl::AutoDScalar VogammaH2plus=VogammaH2plus_node[i1];
            adtl::AutoDScalar VogammaH=VogammaH_node[i1];
            adtl::AutoDScalar VogammaH2=VogammaH2_node[i1];
            double hplus;
            double eminus;
            double Hplus;
            double H2;

            /* get a b h */
            double a2[6];
            if(!transient_state||int(time/dt)==0){
                 for(int index=0;index<6;index++){
                     a2[index]=0;
                 }
            }else{
                for(int index=0;index<6;index++){
                    a2[index]=1;
                }
            }
            double u1 = pow(cm,3)/s;  // unit cm^3/s
            double u2 = 1/s;  // unit 1/s
            double b2[]={1.03e-13*u1*hplus+3.21e-138*u2,
                1.96e-19*u1*H2+1.97e-14*u1*eminus+4.16e3*u2,
               5.04e-22*u2+2.06e-7*u1*eminus,
                2.06e-7*u1*eminus+5.75e5*u2,
               2.06e-19*u1*Hplus+1.03e-13*u1*hplus+5.07e-113*u2,
               1.03e-13*u1*hplus+3.21e-138*u2};

            adtl::AutoDScalar h2[]=
            {-1*(1.97e-14*u1*Vogammaplus*eminus  +4.16e3*u2*Vogammaplus   -8.21e-37*u1*Vogamma*Hplus  +5.04e-22*u2*VogammaHplus),
             -1*(1.03e-13*u1*Vogamma*hplus       +1.90e5*u2*VogammaH2plus +1.03e-19*u1*VogammaH*Hplus+3.21e-138*u2*Vogamma),
             -1*(1.03e-13*u1*VogammaH*hplus    +8.21e-37*u2*Vogamma*Hplus+5.07e-113*u2*VogammaH),
             -1*(1.92e-19*u1*Vogammaplus*H2      +3.81e5*u2*VogammaH2plus +1.26e-62*u2*VogammaHplus    +2.06e-7*u1*VogammaHplus*eminus),
             -1*( 2.06e-7*u1*VogammaH2plus*eminus+4.16e3*u2*VogammaH2plus)
            };

            //get fVogamma
            fVogamma_node[i1] = Var_func->get_fVogamma(a2[0], Vogamma_plot[i1], dt,
                                                       b2[0],
                                                       Vogamma,
                                                       h2[0],
                                                       (*ele_volume)[i1]);
            fVogammaplus_node[i1] = Var_func->get_fVogamma(a2[1], Vogammaplus_plot[i1], dt,
                                                       b2[1],
                                                       Vogammaplus,
                                                       h2[1],
                                                       (*ele_volume)[i1]);
            fVogammaHplus_node[i1] = Var_func->get_fVogamma(a2[2], VogammaHplus_plot[i1], dt,
                                                       b2[2],
                                                       VogammaHplus,
                                                       h2[2],
                                                       (*ele_volume)[i1]);
            fVogammaH2plus_node[i1] = Var_func->get_fVogamma(a2[3], VogammaH2plus_plot[i1], dt,
                                                       b2[3],
                                                       VogammaH2plus,
                                                       h2[3],
                                                       (*ele_volume)[i1]);
            fVogammaH_node[i1] = Var_func->get_fVogamma(a2[4], VogammaH_plot[i1], dt,
                                                       b2[4],
                                                       VogammaH,
                                                       h2[4],
                                                       (*ele_volume)[i1]);
            fVogammaH2_node[i1] = Var_func->get_fVogamma(a2[5], VogammaH2_plot[i1], dt,
                                                       b2[5],
                                                       VogammaH2,
                                                       h2[5],
                                                       (*ele_volume)[i1]);
        }

        for(int l = 0; l < n_vertex; ++l){
            for(int j = 0; j < n_vertex; ++j){
                //get dfVogamma/dx
                mat_data->addMatrixValue(mapping_Vogamma[l], mapping_Vogamma[j], fVogamma_node[l].getADValue(6*j));
                mat_data->addMatrixValue(mapping_Vogamma[l], mapping_Vogamma[j]+1, fVogamma_node[l].getADValue(6*j+1));
                mat_data->addMatrixValue(mapping_Vogamma[l], mapping_Vogamma[j]+2, fVogamma_node[l].getADValue(6*j+2));
                mat_data->addMatrixValue(mapping_Vogamma[l], mapping_Vogamma[j]+3, fVogamma_node[l].getADValue(6*j+3));
                mat_data->addMatrixValue(mapping_Vogamma[l], mapping_Vogamma[j]+4, fVogamma_node[l].getADValue(6*j+4));
                mat_data->addMatrixValue(mapping_Vogamma[l], mapping_Vogamma[j]+5, fVogamma_node[l].getADValue(6*j+5));
                //get dfVogammaplus/dx
                mat_data->addMatrixValue(mapping_Vogamma[l]+1, mapping_Vogamma[j], fVogammaplus_node[l].getADValue(6*j));
                mat_data->addMatrixValue(mapping_Vogamma[l]+1, mapping_Vogamma[j]+1, fVogammaplus_node[l].getADValue(6*j+1));
                mat_data->addMatrixValue(mapping_Vogamma[l]+1, mapping_Vogamma[j]+2, fVogammaplus_node[l].getADValue(6*j+2));
                mat_data->addMatrixValue(mapping_Vogamma[l]+1, mapping_Vogamma[j]+3, fVogammaplus_node[l].getADValue(6*j+3));
                mat_data->addMatrixValue(mapping_Vogamma[l]+1, mapping_Vogamma[j]+4, fVogammaplus_node[l].getADValue(6*j+4));
                mat_data->addMatrixValue(mapping_Vogamma[l]+1, mapping_Vogamma[j]+5, fVogammaplus_node[l].getADValue(6*j+5));
                //get dfVogammaHplus/dx
                mat_data->addMatrixValue(mapping_Vogamma[l]+2, mapping_Vogamma[j], fVogammaHplus_node[l].getADValue(6*j));
                mat_data->addMatrixValue(mapping_Vogamma[l]+2, mapping_Vogamma[j]+1, fVogammaHplus_node[l].getADValue(6*j+1));
                mat_data->addMatrixValue(mapping_Vogamma[l]+2, mapping_Vogamma[j]+2, fVogammaHplus_node[l].getADValue(6*j+2));
                mat_data->addMatrixValue(mapping_Vogamma[l]+2, mapping_Vogamma[j]+3, fVogammaHplus_node[l].getADValue(6*j+3));
                mat_data->addMatrixValue(mapping_Vogamma[l]+2, mapping_Vogamma[j]+4, fVogammaHplus_node[l].getADValue(6*j+4));
                mat_data->addMatrixValue(mapping_Vogamma[l]+2, mapping_Vogamma[j]+5, fVogammaHplus_node[l].getADValue(6*j+5));
                //get dfVogammaH2plus/dx
                mat_data->addMatrixValue(mapping_Vogamma[l]+3, mapping_Vogamma[j], fVogammaH2plus_node[l].getADValue(6*j));
                mat_data->addMatrixValue(mapping_Vogamma[l]+3, mapping_Vogamma[j]+1, fVogammaH2plus_node[l].getADValue(6*j+1));
                mat_data->addMatrixValue(mapping_Vogamma[l]+3, mapping_Vogamma[j]+2, fVogammaH2plus_node[l].getADValue(6*j+2));
                mat_data->addMatrixValue(mapping_Vogamma[l]+3, mapping_Vogamma[j]+3, fVogammaH2plus_node[l].getADValue(6*j+3));
                mat_data->addMatrixValue(mapping_Vogamma[l]+3, mapping_Vogamma[j]+4, fVogammaH2plus_node[l].getADValue(6*j+4));
                mat_data->addMatrixValue(mapping_Vogamma[l]+3, mapping_Vogamma[j]+5, fVogammaH2plus_node[l].getADValue(6*j+5));
                //get dfVogammaH/dx
                mat_data->addMatrixValue(mapping_Vogamma[l]+4, mapping_Vogamma[j], fVogammaH_node[l].getADValue(6*j));
                mat_data->addMatrixValue(mapping_Vogamma[l]+4, mapping_Vogamma[j]+1, fVogammaH_node[l].getADValue(6*j+1));
                mat_data->addMatrixValue(mapping_Vogamma[l]+4, mapping_Vogamma[j]+2, fVogammaH_node[l].getADValue(6*j+2));
                mat_data->addMatrixValue(mapping_Vogamma[l]+4, mapping_Vogamma[j]+3, fVogammaH_node[l].getADValue(6*j+3));
                mat_data->addMatrixValue(mapping_Vogamma[l]+4, mapping_Vogamma[j]+4, fVogammaH_node[l].getADValue(6*j+4));
                mat_data->addMatrixValue(mapping_Vogamma[l]+4, mapping_Vogamma[j]+5, fVogammaH_node[l].getADValue(6*j+5));
                //get dfVogammaH2/dx
                mat_data->addMatrixValue(mapping_Vogamma[l]+5, mapping_Vogamma[j], fVogammaH2_node[l].getADValue(6*j));
                mat_data->addMatrixValue(mapping_Vogamma[l]+5, mapping_Vogamma[j]+1, fVogammaH2_node[l].getADValue(6*j+1));
                mat_data->addMatrixValue(mapping_Vogamma[l]+5, mapping_Vogamma[j]+2, fVogammaH2_node[l].getADValue(6*j+2));
                mat_data->addMatrixValue(mapping_Vogamma[l]+5, mapping_Vogamma[j]+3, fVogammaH2_node[l].getADValue(6*j+3));
                mat_data->addMatrixValue(mapping_Vogamma[l]+5, mapping_Vogamma[j]+4, fVogammaH2_node[l].getADValue(6*j+4));
                mat_data->addMatrixValue(mapping_Vogamma[l]+5, mapping_Vogamma[j]+5, fVogammaH2_node[l].getADValue(6*j+5));
            }
            vec_data->addVectorValue(mapping_Vogamma[l], -fVogamma_node[l].getValue());
            vec_data->addVectorValue(mapping_Vogamma[l]+1, -fVogammaplus_node[l].getValue());
            vec_data->addVectorValue(mapping_Vogamma[l]+2, -fVogammaHplus_node[l].getValue());
            vec_data->addVectorValue(mapping_Vogamma[l]+3, -fVogammaH2plus_node[l].getValue());
            vec_data->addVectorValue(mapping_Vogamma[l]+4, -fVogammaH_node[l].getValue());
            vec_data->addVectorValue(mapping_Vogamma[l]+5, -fVogammaH2_node[l].getValue());
        }
    }
    mat_data->assemble();
}


//the below code will not run

/**
 * @brief NonlinearPoisson::applyNewtonPBCOnPatch: Vector
 * @param patch
 * @param time
 * @param dt
 * @param component_name
 */
void NonlinearPoisson::applyNewtonPBCOnPatch(hier::Patch<NDIM> & patch,
                                          const double  time,
                                          const double  dt,
                                          const string& component_name){
    /** 取出本地PatchGeometry **/
    tbox::Pointer< hier::PatchGeometry<NDIM> > patch_geo = patch.getPatchGeometry();
    /** get PatchTopology **/
    tbox::Pointer< hier::PatchTopology<NDIM> > patch_top = patch.getPatchTopology();

    /** (形函数, 积分器, 单元)管理器 **/
    tbox::Pointer<ElementManager<NDIM> >
            ele_manager = ElementManager<NDIM>::getManager();

    /** 取出本地Patch的结点坐标数组 **/
    tbox::Pointer<pdat::NodeData<NDIM, double> > node_coord = patch_geo->getNodeCoordinates();

    /** setup the source target **/
    tbox::Pointer<bound_condition> bound_source = new bound_condition();

    /** get node values **/
    tbox::Pointer<pdat::VectorData<NDIM, double> >
            V_soln = patch.getPatchData(d_poisson_sol_cur_id);

    int num_nodes = patch.getNumberOfNodes(1);

    int * dof_map_semi = d_dof_semi_info->getDOFMapping(patch, hier::EntityUtilities::NODE);

    tbox::Array<int> face_node_ext, face_node_idx;
    patch_top->getFaceAdjacencyNodes(face_node_ext, face_node_idx);

    /** assign boundary condition to right side vector **/
    /** obtain the right side vector **/
    tbox::Pointer<pdat::VectorData<NDIM,double> >
            fval_data = patch.getPatchData(d_poisson_newton_fval_id);
    double * fval_val = fval_data->getPointer();
    MPI_Barrier(MPI_COMM_WORLD);

    int num_bc = EBoundassignment.getSize();
    for(int i_bc=0; i_bc<num_bc; ++i_bc){
        if (patch_geo->hasEntitySet(EBoundassignment[i_bc].bc_no, hier::EntityUtilities::NODE)){
            const tbox::Array<int> &entity_idx = patch_geo->getEntityIndicesInSet(EBoundassignment[i_bc].bc_no,
                                                 hier::EntityUtilities::NODE,
                                                 num_nodes);
            int size = entity_idx.getSize();

            if(EBoundassignment[i_bc].bc_type == "ohmic"){
                for (int i = 0; i < size; ++i) {
                    int index = dof_map_semi[entity_idx[i]];
                    fval_val[index] = 0.0;
                    fval_val[index+1] = 0.0;
                    fval_val[index+2] = 0.0;
                }
            }else if (EBoundassignment[i_bc].bc_type == "gate-simplified"){
                /** load fixed volt value **/
                double Vapp = bound_source->source_assign(EBoundassignment[i_bc].sc_type,time, dt,EBoundassignment[i_bc].sc_coe);
                double WORKFUNC = EBoundassignment[i_bc].sc_coe[3];
                double dins = EBoundassignment[i_bc].sc_coe[4]*m;
                double epsilon_ins = EBoundassignment[i_bc].sc_coe[5];
                double local_Qt = EBoundassignment[i_bc].sc_coe[6];
                double Qins = local_Qt*pow(m,-2);
                if(time<=1e-9){
                    Qins = 0;
                    //Vapp = 0;
                }
                double value_rhs = (e*Qins/eps0 + epsilon_ins * (Vapp - WORKFUNC)/dins);
                tbox::Pointer<BaseElement<NDIM> > ele_quadrangle = ele_manager->getElement("LinearQuadrangle");
                int num_faces = patch.getNumberOfEntities(hier::EntityUtilities::FACE, 1);
                int num_face = EBoundassignment[i_bc].bc_face.getSize();
                for (int i_f = 0; i_f < num_face; ++i_f){
                    if(patch_geo->hasEntitySet(EBoundassignment[i_bc].bc_face[i_f], hier::EntityUtilities::FACE)){
                        //cout<<"get the gate-simplified face:"<<EBoundassignment[i_bc].bc_face[i_f]<<endl;
                        const tbox::Array<int> & face_id = patch_geo->getEntityIndicesInSet(EBoundassignment[i_bc].bc_face[i_f],
                                                           hier::EntityUtilities::FACE,
                                                           num_faces);
                        int size_f = face_id.size();
                        for (int i_e = 0; i_e < size_f; ++i_e){
                            int n_fvertex = face_node_ext[face_id[i_e] + 1] - face_node_ext[face_id[i_e]];
                            tbox::Array<hier::DoubleVector<NDIM> > face_vertex(n_fvertex);
                            tbox::Array<int> node_n(n_fvertex);
                            tbox::Array<int> mapping_semi(n_fvertex);
                            tbox::Array<double> V_node(n_fvertex);
                            tbox::Array<double> value(n_fvertex);

                            /** inform for node **/
                            for (int i1 = 0, j = face_node_ext[face_id[i_e]]; i1 < n_fvertex; ++i1, ++j) {
                                node_n[i1] = face_node_idx[j];
                                mapping_semi[i1] = dof_map_semi[face_node_idx[j]];
                                V_node[i1] = (*V_soln)(mapping_semi[i1]);
                                for(int k = 0; k< NDIM; ++k){
                                    face_vertex[i1][k] = (*node_coord)(k, face_node_idx[j])*m;
                                }
                            }

                            /** element vector RHS **/
                            tbox::Pointer<tbox::Vector<double> > ele_vec = new tbox::Vector<double>();
                            ele_vec->resize(n_fvertex);

                            for (int i_n = 0; i_n < n_fvertex; ++i_n){
                                value[i_n] = - value_rhs + epsilon_ins*V_node[i_n]/dins;
                            }
                            /** building element vector **/
                            ele_quadrangle->buildCVFEMRhs(face_vertex, value, dt, time, ele_vec);

                            /** construct fval vector **/
                            for(int i2= 0; i2 < n_fvertex; ++i2)
                                fval_val[mapping_semi[i2]] += (*ele_vec)[i2];
                        }
                    }else{
                        //cout<<"can't get the gate-simplified face:"<<EBoundassignment[i_bc].bc_face[i_f]<<endl;
                    }
                }
            }else if (EBoundassignment[i_bc].bc_type == "gate-diode"){
                string filename = "diode_voltage";
                std::fstream input(filename.c_str());
                if(!input) {
                    TBOX_ERROR("fail to open file: "<< filename << endl);
                }
                string line;
                getline(input,line);
                //line.clear();
                double diode_voltage = 0;
                stringstream ss(line);
                ss>>diode_voltage;
                input.close();
                //cout<<"diode_voltage:"<<diode_voltage<<endl;
                /** load fixed volt value **/
                double dins = EBoundassignment[i_bc].sc_coe[4]*m;
                double epsilon_ins = EBoundassignment[i_bc].sc_coe[5];
                double local_Qt = EBoundassignment[i_bc].sc_coe[6];
                double Qins = local_Qt*pow(m,-2);
                if(time<=1e-9){
                    Qins = 0;
                    diode_voltage = 0;
                }
                double value_rhs = (e*Qins/eps0 + epsilon_ins * (diode_voltage)/dins);
                tbox::Pointer<BaseElement<NDIM> > ele_quadrangle = ele_manager->getElement("LinearQuadrangle");
                int num_faces = patch.getNumberOfEntities(hier::EntityUtilities::FACE, 1);
                int num_face = EBoundassignment[i_bc].bc_face.getSize();
                for (int i_f = 0; i_f < num_face; ++i_f){
                    if(patch_geo->hasEntitySet(EBoundassignment[i_bc].bc_face[i_f], hier::EntityUtilities::FACE)){
                        //cout<<"get the gate-simplified face:"<<EBoundassignment[i_bc].bc_face[i_f]<<endl;
                        const tbox::Array<int> & face_id = patch_geo->getEntityIndicesInSet(EBoundassignment[i_bc].bc_face[i_f],
                                                           hier::EntityUtilities::FACE,
                                                           num_faces);
                        int size_f = face_id.size();
                        for (int i_e = 0; i_e < size_f; ++i_e){
                            int n_fvertex = face_node_ext[face_id[i_e] + 1] - face_node_ext[face_id[i_e]];
                            tbox::Array<hier::DoubleVector<NDIM> > face_vertex(n_fvertex);
                            tbox::Array<int> node_n(n_fvertex);
                            tbox::Array<int> mapping_semi(n_fvertex);
                            tbox::Array<double> V_node(n_fvertex);
                            tbox::Array<double> value(n_fvertex);

                            /** inform for node **/
                            for (int i1 = 0, j = face_node_ext[face_id[i_e]]; i1 < n_fvertex; ++i1, ++j) {
                                node_n[i1] = face_node_idx[j];
                                mapping_semi[i1] = dof_map_semi[face_node_idx[j]];
                                V_node[i1] = (*V_soln)(mapping_semi[i1]);
                                for(int k = 0; k< NDIM; ++k){
                                    face_vertex[i1][k] = (*node_coord)(k, face_node_idx[j])*m;
                                }
                            }

                            /** element vector RHS **/
                            tbox::Pointer<tbox::Vector<double> > ele_vec = new tbox::Vector<double>();
                            ele_vec->resize(n_fvertex);

                            for (int i_n = 0; i_n < n_fvertex; ++i_n){
                                value[i_n] = - value_rhs + epsilon_ins*V_node[i_n]/dins;
                            }
                            /** building element vector **/
                            ele_quadrangle->buildCVFEMRhs(face_vertex, value, dt, time, ele_vec);

                            /** construct fval vector **/
                            for(int i2= 0; i2 < n_fvertex; ++i2)
                                fval_val[mapping_semi[i2]] += (*ele_vec)[i2];
                        }
                    }else{
                        cout<<"can't get the gate-diode face:"<<EBoundassignment[i_bc].bc_face[i_f]<<endl;
                    }
                }
            }
        }
    }
}

/**
 * @brief NonlinearPoisson::applyPoissonBCOnPatch: vector
 * @param patch
 * @param time
 * @param dt
 * @param component_name
 */
void NonlinearPoisson::applyPoissonBCOnPatch(hier::Patch<NDIM> & patch,
                    const double  time,
                    const double  dt,
                    const string& component_name){
    /** 取出本地PatchGeometry **/
    tbox::Pointer< hier::PatchGeometry<NDIM> > patch_geo = patch.getPatchGeometry();
    /** get PatchTopology **/
    tbox::Pointer< hier::PatchTopology<NDIM> > patch_top = patch.getPatchTopology();

    /** obtain the right side vector **/
    tbox::Pointer<pdat::VectorData<NDIM,double> >
            fval_data = patch.getPatchData(d_poisson_fval_id);

    /** (形函数, 积分器, 单元)管理器 **/
    tbox::Pointer<ElementManager<NDIM> >
            ele_manager = ElementManager<NDIM>::getManager();

    /** 取出本地Patch的结点坐标数组 **/
    tbox::Pointer<pdat::NodeData<NDIM, double> > node_coord = patch_geo->getNodeCoordinates();

    /** setup the source target **/
    tbox::Pointer<bound_condition> bound_source = new bound_condition();

    /** get node values **/
    tbox::Pointer<pdat::VectorData<NDIM, double> >
            V_soln = patch.getPatchData(d_poisson_sol_cur_id);

    int num_nodes = patch.getNumberOfNodes(1);

    int * dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);
    int * dof_map_semi = d_dof_semi_info->getDOFMapping(patch, hier::EntityUtilities::NODE);

    tbox::Array<int> face_node_ext, face_node_idx;
    patch_top->getFaceAdjacencyNodes(face_node_ext, face_node_idx);

    /** assign boundary condition to right side vector **/
    double * fval_val = fval_data->getPointer();

    int num_bc = EBoundassignment.getSize();
    for(int i_bc=0; i_bc<num_bc; ++i_bc){
        if (patch_geo->hasEntitySet(EBoundassignment[i_bc].bc_no, hier::EntityUtilities::NODE)){
            const tbox::Array<int> &entity_idx = patch_geo->getEntityIndicesInSet(EBoundassignment[i_bc].bc_no,
                                                 hier::EntityUtilities::NODE,
                                                 num_nodes);
            int size = entity_idx.getSize();

            if(EBoundassignment[i_bc].bc_type == "ohmic"){
                for (int i = 0; i < size; ++i) {
                    int index = dof_map[entity_idx[i]];
                    fval_val[index] = 0.0;
                }
            }else if (EBoundassignment[i_bc].bc_type == "schottky"){
                for (int i = 0; i < size; ++i) {
                    int index = dof_map[entity_idx[i]];
                    fval_val[index] = 0.0;
                }
            }else if (EBoundassignment[i_bc].bc_type == "gate-simplified"){
                /** load fixed volt value **/
                double dins = EBoundassignment[i_bc].sc_coe[4]*m;
                double epsilon_ins = EBoundassignment[i_bc].sc_coe[5];
                double Vapp = bound_source->source_assign(EBoundassignment[i_bc].sc_type,time, dt,EBoundassignment[i_bc].sc_coe);
                double WORKFUNC = EBoundassignment[i_bc].sc_coe[3];
                double local_Qt = EBoundassignment[i_bc].sc_coe[6];
                double Qins = local_Qt*pow(m,-2);
                if(int(time/dt)==0){
                    Qins = 0;
                    //Vapp = 0;
                }
                double value_rhs = (e*Qins/eps0 + epsilon_ins * (Vapp - WORKFUNC)/dins);
                tbox::Pointer<BaseElement<NDIM> > ele_quadrangle = ele_manager->getElement("LinearQuadrangle");
                int num_faces = patch.getNumberOfEntities(hier::EntityUtilities::FACE, 1);
                int num_face = EBoundassignment[i_bc].bc_face.getSize();
                for (int i_f = 0; i_f < num_face; ++i_f){
                    if(patch_geo->hasEntitySet(EBoundassignment[i_bc].bc_face[i_f], hier::EntityUtilities::FACE)){
                        //cout<<"get the gate-simplified face:"<<EBoundassignment[i_bc].bc_face[i_f]<<endl;
                        const tbox::Array<int> & face_id = patch_geo->getEntityIndicesInSet(EBoundassignment[i_bc].bc_face[i_f],
                                                           hier::EntityUtilities::FACE,
                                                           num_faces);
                        int size_f = face_id.size();
                        for (int i_e = 0; i_e < size_f; ++i_e){
                            int n_fvertex = face_node_ext[face_id[i_e] + 1] - face_node_ext[face_id[i_e]];
                            tbox::Array<hier::DoubleVector<NDIM> > face_vertex(n_fvertex);
                            tbox::Array<int> node_n(n_fvertex);
                            tbox::Array<int> mapping(n_fvertex);
                            tbox::Array<int> mapping_semi(n_fvertex);
                            tbox::Array<double> V_node(n_fvertex);
                            tbox::Array<double> value(n_fvertex);

                            /** inform for node **/
                            for (int i1 = 0, j = face_node_ext[face_id[i_e]]; i1 < n_fvertex; ++i1, ++j) {
                                node_n[i1] = face_node_idx[j];
                                mapping[i1] = dof_map[face_node_idx[j]];
                                mapping_semi[i1] = dof_map_semi[face_node_idx[j]];
                                V_node[i1] = (*V_soln)(mapping_semi[i1]);
                                for(int k = 0; k< NDIM; ++k){
                                    face_vertex[i1][k] = (*node_coord)(k, face_node_idx[j])*m;
                                }
                            }

                            /** element vector RHS **/
                            tbox::Pointer<tbox::Vector<double> > ele_vec = new tbox::Vector<double>();
                            ele_vec->resize(n_fvertex);

                            for (int i_n = 0; i_n < n_fvertex; ++i_n){
                                value[i_n] = - value_rhs + epsilon_ins*V_node[i_n]/dins;
                            }
                            /** building element vector **/
                            ele_quadrangle->buildCVFEMRhs(face_vertex, value, dt, time, ele_vec);

                            /** construct fval vector **/
                            for(int i2= 0; i2 < n_fvertex; ++i2)
                                fval_val[mapping[i2]] += (*ele_vec)[i2];
                        }
                    }else{
                        cout<<"can't get the gate-simplified face:"<<EBoundassignment[i_bc].bc_face[i_f]<<endl;
                    }
                }
            }
        }
    }
}

/**
 * @brief NonlinearPoisson::applyCarnBCOnPatch: Vector
 * @param patch
 * @param time
 * @param dt
 * @param component_name
 */
void NonlinearPoisson::applyCarnBCOnPatch(hier::Patch<NDIM> & patch,
                      const double  time,
                      const double  dt,
                      const string& component_name){
    /** 取出本地PatchGeometry **/
    tbox::Pointer< hier::PatchGeometry<NDIM> > patch_geo = patch.getPatchGeometry();
    /** Get Local PatchTopology **/
    tbox::Pointer< hier::PatchTopology<NDIM> > patch_top = patch.getPatchTopology();
    /** obtain the current solution vector **/
    tbox::Pointer<pdat::VectorData<NDIM,double> >
      fval_data = patch.getPatchData(d_carn_fval_id);


    tbox::Pointer<pdat::VectorData<NDIM,double> >
      V_soln = patch.getPatchData(d_poisson_sol_cur_id);
    tbox::Pointer<pdat::NodeData<NDIM, double> >
      dop = patch.getPatchData(d_dopant_id);

    /** quasi fermi level **/
    tbox::Pointer<pdat::VectorData<NDIM,double> >
      Fn_soln = patch.getPatchData(d_fermin_sol_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
      Fh_soln = patch.getPatchData(d_fermih_sol_id);

    /** face - entity relationship **/
    tbox::Pointer<pdat::FaceData<NDIM, int> >
      face_mat = patch.getPatchData(d_face_mat_id);
    tbox::Pointer<pdat::NodeData<NDIM, int> >
      node_mat = patch.getPatchData(d_node_mat_id);
    /** heterojunction - node relationship **/
    tbox::Pointer<pdat::NodeData<NDIM, int> >
      heter_node = patch.getPatchData(d_heter_node_id);

    /** 取出本地Patch的结点坐标数组 **/
    tbox::Pointer<pdat::NodeData<NDIM, double> >
      node_coord = patch_geo->getNodeCoordinates();

    /** 面——单元邻接关系数据 **/
    tbox::Array<int>face_cell_ext,face_cell_idx;
    patch_top->getFaceAdjacencyCells(face_cell_ext,face_cell_idx);
    /** 单元节点邻接关系数据 **/
    tbox::Array<int>face_node_ext,face_node_idx;
    patch_top->getFaceAdjacencyNodes(face_node_ext,face_node_idx);

    int num_nodes = patch.getNumberOfNodes(1);
    //int num_faces = patch.getNumberOfEntities(hier::EntityUtilities::FACE,1);

    tbox::Pointer<Variables_Functions> Var_func = new Variables_Functions();

    int * dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);
    int * dof_map_semi = d_dof_semi_info->getDOFMapping(patch, hier::EntityUtilities::NODE);

    /** 取出矩阵向量的核心数据结构的指针 **/
    double *fval_val = fval_data->getPointer();

    /** material list **/
    tbox::Pointer<MaterialManager<NDIM> >
            mat_manager = MaterialManager<NDIM>::getManager();

    /** setup the source target **/
    tbox::Pointer<bound_condition> bound_source = new bound_condition();

    double tmp_val;

    /** boundary condition source **/
    int num_bc = EBoundassignment.getSize();
    tbox::Array<double> bound_coe;
    bound_coe.resizeArray(9);

    for(int i_bc=0; i_bc<num_bc; ++i_bc){
        if (EBoundassignment[i_bc].bc_type == "ohmic"){
            if (patch_geo->hasEntitySet(EBoundassignment[i_bc].bc_no, hier::EntityUtilities::NODE)){
                const tbox::Array<int> &entity_idx =
                       patch_geo->getEntityIndicesInSet(EBoundassignment[i_bc].bc_no,
                                             hier::EntityUtilities::NODE,
                                             num_nodes);
                double app_volt = bound_source->source_assign(EBoundassignment[i_bc].sc_type,
                                                              time,
                                                              dt,
                                                              EBoundassignment[i_bc].sc_coe);
                string contact_type = "Ohmic_contact";
                int size = entity_idx.getSize();
                for (int i = 0; i < size; ++i) {
                    int index = dof_map[entity_idx[i]];
                    int index_s = dof_map_semi[entity_idx[i]];

                    /** setup material target **/
                    int dom_no;
                    dom_no = 0;
                    for(int i_m = 0; i_m < material_list.size(); ++i_m){
                        if(material_list[i_m].mat_type == "semiconductor"){
                            if((*node_mat)(i_m,entity_idx[i])==1){
                                dom_no = i_m;
                            }
                        }
                    }

                    tbox::Pointer<BaseMaterial<NDIM> > mat_dom =
                            mat_manager->getMaterial(material_list[dom_no].mat_name);

                    bound_coe[0] = app_volt;
                    bound_coe[1] = T;
                    bound_coe[2] = mat_dom->effectiveAffinity(T_amb);
                    bound_coe[3] = mat_dom->bandGap(T_amb);
                    bound_coe[4] = (*Fn_soln)(index);
                    bound_coe[5] = (*Fh_soln)(index);
                    bound_coe[6] = (*V_soln)(index_s);
                    bound_coe[7] = (*dop)(0,entity_idx[i]);
                    bound_coe[8] = mat_dom->intrinsicDensityNi(T_amb)*pow(m,-3);

                    tmp_val = bound_source->carn_bound_value(contact_type, car_bound_mode, time, dt, kb, e, bound_coe);
                    fval_val[index] = tmp_val;
                }
            }
        }
    }
}

/**
 * @brief NonlinearPoisson::setCarhPhysicalBC: Matrix
 * @param patch
 * @param time
 * @param dt
 * @param component_name
 */
void NonlinearPoisson::setCarhPhysicalBC(hier::Patch<NDIM> & patch,
                            const double  time,
                            const double  dt,
                            const string& component_name){
    /** 取出本地PatchGeometry **/
    tbox::Pointer< hier::PatchGeometry<NDIM> > patch_geo = patch.getPatchGeometry();
    /** Get Local PatchTopology **/
    tbox::Pointer< hier::PatchTopology<NDIM> > patch_top = patch.getPatchTopology();
    /** 获取网格片矩阵对象 **/
    tbox::Pointer<pdat::CSRMatrixData<NDIM,double> >
      mat_data = patch.getPatchData(d_carh_mat_id);

    /** face - entity relationship **/
    tbox::Pointer<pdat::FaceData<NDIM, int> >
      face_mat= patch.getPatchData(d_face_mat_id);

    /** 取出本地Patch的结点坐标数组 **/
    tbox::Pointer<pdat::NodeData<NDIM, double> >
      node_coord = patch_geo->getNodeCoordinates();
    /** 面——单元邻接关系数据 **/
    tbox::Array<int>face_cell_ext,face_cell_idx;
    patch_top->getFaceAdjacencyCells(face_cell_ext,face_cell_idx);
    /** 单元节点邻接关系数据 **/
    tbox::Array<int>face_node_ext,face_node_idx;
    patch_top->getFaceAdjacencyNodes(face_node_ext,face_node_idx);

    int num_nodes = patch.getNumberOfNodes(1);
    //int num_faces = patch.getNumberOfEntities(hier::EntityUtilities::FACE,1);

    /** material list **/
    tbox::Pointer<MaterialManager<NDIM> >
            mat_manager = MaterialManager<NDIM>::getManager();

    tbox::Pointer<Variables_Functions> Var_func = new Variables_Functions();

    int * dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);
    //int * dof_map_semi = d_dof_semi_info->getDOFMapping(patch, hier::EntityUtilities::NODE);

    /** 取出矩阵向量的核心数据结构的指针 **/
    int * row_start = mat_data->getRowStartPointer();
    int * col_idx   = mat_data->getColumnIndicesPointer();
    double * mat_val    = mat_data->getValuePointer();

    /** boundary condition source **/
    int num_bc = EBoundassignment.getSize();

    for(int i_bc=0; i_bc<num_bc; ++i_bc){
        if (EBoundassignment[i_bc].bc_type == "ohmic"){
            if (patch_geo->hasEntitySet(EBoundassignment[i_bc].bc_no, hier::EntityUtilities::NODE)){
                const tbox::Array<int> &entity_idx = patch_geo->getEntityIndicesInSet(EBoundassignment[i_bc].bc_no,
                                                     hier::EntityUtilities::NODE,
                                                     num_nodes);
                int size = entity_idx.getSize();

                for (int i = 0; i < size; ++i) {
                    int index = dof_map[entity_idx[i]];
                    for (int j = row_start[index]; j < row_start[index+1]; ++j){
                        if (col_idx[j] == index){
                            mat_val[j] = 1.0;   /**< 对角线元素 */
                        }
                        else {
                            mat_val[j] = 0.0;   /**< 非对角线元素 */
                        }
                    }
                }
            }
        }
    }
}

/**
 * @brief NonlinearPoisson::applyCarhBCOnPatch: Vector
 * @param patch
 * @param time
 * @param dt
 * @param component_name
 */
void NonlinearPoisson::applyCarhBCOnPatch(hier::Patch<NDIM> & patch,
                      const double  time,
                      const double  dt,
                      const string& component_name){
    /** 取出本地PatchGeometry **/
    tbox::Pointer< hier::PatchGeometry<NDIM> > patch_geo = patch.getPatchGeometry();
    /** Get Local PatchTopology **/
    tbox::Pointer< hier::PatchTopology<NDIM> > patch_top = patch.getPatchTopology();
    /** obtain the current solution vector **/
    tbox::Pointer<pdat::VectorData<NDIM,double> >
      fval_data = patch.getPatchData(d_carh_fval_id);

    tbox::Pointer<pdat::VectorData<NDIM,double> >
      V_soln = patch.getPatchData(d_poisson_sol_cur_id);

    tbox::Pointer<pdat::NodeData<NDIM, double> >
      dop = patch.getPatchData(d_dopant_id);

    /** quasi fermi level **/
    tbox::Pointer<pdat::VectorData<NDIM,double> >
      Fn_soln = patch.getPatchData(d_fermin_sol_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
      Fh_soln = patch.getPatchData(d_fermih_sol_id);

    /** face - entity relationship **/
    tbox::Pointer<pdat::FaceData<NDIM, int> >
      face_mat = patch.getPatchData(d_face_mat_id);
    tbox::Pointer<pdat::NodeData<NDIM, int> >
      node_mat = patch.getPatchData(d_node_mat_id);

    /** 取出本地Patch的结点坐标数组 **/
    tbox::Pointer<pdat::NodeData<NDIM, double> >
      node_coord = patch_geo->getNodeCoordinates();

    /** 面——单元邻接关系数据 **/
    tbox::Array<int>face_cell_ext,face_cell_idx;
    patch_top->getFaceAdjacencyCells(face_cell_ext,face_cell_idx);
    /** 单元节点邻接关系数据 **/
    tbox::Array<int>face_node_ext,face_node_idx;
    patch_top->getFaceAdjacencyNodes(face_node_ext,face_node_idx);

    int num_nodes = patch.getNumberOfNodes(1);
    //int num_faces = patch.getNumberOfEntities(hier::EntityUtilities::FACE,1);

    int * dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);
    int * dof_map_semi = d_dof_semi_info->getDOFMapping(patch, hier::EntityUtilities::NODE);

    /** 取出矩阵向量的核心数据结构的指针 **/
    double * fval_val = fval_data->getPointer();
    /** setup the source target **/
    tbox::Pointer<bound_condition> bound_source = new bound_condition();

    /** material list **/
    tbox::Pointer<MaterialManager<NDIM> >
            mat_manager = MaterialManager<NDIM>::getManager();

    tbox::Pointer<Variables_Functions> Var_func = new Variables_Functions();

    /** boundary condition source **/
    double tmp_val;
    int num_bc = EBoundassignment.getSize();

    tbox::Array<double> bound_coe;
    bound_coe.resizeArray(9);

    for(int i_bc=0; i_bc<num_bc; ++i_bc){
        if (EBoundassignment[i_bc].bc_type == "ohmic"){
            if (patch_geo->hasEntitySet(EBoundassignment[i_bc].bc_no, hier::EntityUtilities::NODE)){

                const tbox::Array<int> &entity_idx =
                       patch_geo->getEntityIndicesInSet(EBoundassignment[i_bc].bc_no,
                                             hier::EntityUtilities::NODE,
                                             num_nodes);

                double app_volt = bound_source->source_assign(EBoundassignment[i_bc].sc_type,
                                                              time,
                                                              dt,
                                                              EBoundassignment[i_bc].sc_coe);
                string contact_type = "Ohmic_contact";
                int size = entity_idx.getSize();
                for (int i = 0; i < size; ++i) {
                    int index = dof_map[entity_idx[i]];
                    int index_s = dof_map_semi[entity_idx[i]];

                    /** setup material target **/
                    int dom_no;
                    dom_no = 0;
                    for(int i_m = 0; i_m < material_list.size(); ++i_m){
                        if(material_list[i_m].mat_type == "semiconductor"){
                            if((*node_mat)(i_m,entity_idx[i])==1){
                                dom_no = i_m;
                            }
                        }
                    }

                    tbox::Pointer<BaseMaterial<NDIM> > mat_dom =
                            mat_manager->getMaterial(material_list[dom_no].mat_name);

                    bound_coe[0] = app_volt;
                    bound_coe[1] = T;
                    bound_coe[2] = mat_dom->effectiveAffinity(T_amb);
                    bound_coe[3] = mat_dom->bandGap(T_amb);
                    bound_coe[4] = (*Fn_soln)(index);
                    bound_coe[5] = (*Fh_soln)(index);
                    bound_coe[6] = (*V_soln)(index_s);
                    bound_coe[7] = (*dop)(0,entity_idx[i]);
                    bound_coe[8] = mat_dom->intrinsicDensityNi(T_amb)*pow(m,-3);

                    tmp_val = bound_source->carh_bound_value(contact_type, car_bound_mode, time, dt, kb, e, bound_coe);
                    fval_val[index] = tmp_val;
                }
            }
        }
    }
}


/**
* @brief NonlinearPoisson::applyConstraint2NP: apply constraint to hole and electron
* @param patch
* @param time
* @param dt
* @param component_name
*/
void NonlinearPoisson::applyConstraint2NP(hier::Patch<NDIM>& patch,
                                        const double time,
                                        const double dt,
                                        const string& component_name){
   /** 取出本地PatchTopology **/
   tbox::Pointer< hier::PatchTopology<NDIM> > patch_top = patch.getPatchTopology();
   tbox::Pointer<pdat::NodeData<NDIM, int> >
     semi_node = patch.getPatchData(d_semi_node_id);
   tbox::Pointer<pdat::NodeData<NDIM, int> >
     sio2_node = patch.getPatchData(d_sio2_node_id);
   /** get node values **/
   tbox::Pointer<pdat::VectorData<NDIM, double> >
     N_soln = patch.getPatchData(d_carn_sol_cur_id);
   tbox::Pointer<pdat::VectorData<NDIM, double> >
     H_soln = patch.getPatchData(d_carh_sol_cur_id);
   /** get node values **/
   tbox::Pointer<pdat::VectorData<NDIM, double> >
     V_soln = patch.getPatchData(d_poisson_sol_cur_id);

   tbox::Pointer<pdat::NodeData<NDIM, int> >
     heter_node = patch.getPatchData(d_heter_node_id);
   tbox::Pointer<pdat::NodeData<NDIM, int> >
     hnode_mat = patch.getPatchData(d_hnode_mat_id);

   tbox::Pointer< pdat::NodeData<NDIM, int> >
           node_mat = patch.getPatchData(d_node_mat_id);

   int * dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);
   int * dof_map_semi = d_dof_semi_info->getDOFMapping(patch, hier::EntityUtilities::NODE);

   /** 获取本地单元(边, 面)周围结点的索引关系 **/
   tbox::Array<int> can_extent, can_indices;
   patch_top->getCellAdjacencyNodes(can_extent, can_indices);
   int num_nodes = patch.getNumberOfNodes(1);

   /** material list **/
   tbox::Pointer<MaterialManager<NDIM> >
           mat_manager = MaterialManager<NDIM>::getManager();

   /** setup target for class 'Variables_Functions' **/
   tbox::Pointer<Variables_Functions> Var_func = new Variables_Functions();

   double scaling_constraint_value = constraint_value*pow(m,-3);

   /** traversing nodes **/
   for(int i = 0; i < num_nodes; ++i){
       if(((*semi_node)(0,i) == 1)||((*sio2_node)(0,i)== 1)){
           int index_s = dof_map_semi[i];
           int index = dof_map[i];
           #if 0
           if(((*N_soln)(index)<=0)&&index==27){
               cout<<"(*N_soln)(index)"<<(*N_soln)(index)<<endl;
               cout<<"N_soln->getPointer()[index]"<<N_soln->getPointer()[index]<<endl;
               cout<<"index:"<<index<<endl;
           }
           #endif
           #if 1
           /** get material for cell **/
           tbox::Pointer<BaseMaterial<NDIM> > mat_cell;
           if((*semi_node)(0,i) == 1) mat_cell = mat_manager->getMaterial("matSilicon");
           if((*sio2_node)(0,i) == 1) mat_cell = mat_manager->getMaterial("matSiO2");
           double ni_tmp = mat_cell->intrinsicDensityNi(T_amb)*pow(m,-3);
           (*N_soln)(index) = Var_func->Constraint((*N_soln)(index),pow(ni_tmp,2)/(*H_soln)(index));
           (*H_soln)(index) = Var_func->Constraint((*H_soln)(index),pow(ni_tmp,2)/(*N_soln)(index));
           (*V_soln)(index_s+1) = Var_func->Constraint((*V_soln)(index_s+1),pow(ni_tmp,2)/(*V_soln)(index_s+2));
           (*V_soln)(index_s+2) = Var_func->Constraint((*V_soln)(index_s+2),pow(ni_tmp,2)/(*V_soln)(index_s+1));
           #else
           (*N_soln)(index) = Var_func->Constraint((*N_soln)(index),scaling_constraint_value);
           (*H_soln)(index) = Var_func->Constraint((*H_soln)(index),scaling_constraint_value);
           (*V_soln)(index_s+1) = Var_func->Constraint((*V_soln)(index_s+1),scaling_constraint_value);
           (*V_soln)(index_s+2) = Var_func->Constraint((*V_soln)(index_s+2),scaling_constraint_value);
           #endif
       }
   }
}


/**
 * @brief NonlinearPoisson::buildCarhFXOnPatch: CV-FEM-SG for hole continuity equation
 * @param patch
 * @param time
 * @param dt
 * @param component_name
 */
void NonlinearPoisson::buildCarhFXOnPatch(hier::Patch<NDIM> & patch,
                            const double  time,
                            const double  dt,
                            const string& component_name){
    /** (形函数, 积分器, 单元)管理器 **/
    tbox::Pointer<ElementManager<NDIM> >
      ele_manager = ElementManager<NDIM>::getManager();
    tbox::Pointer<BaseElement<NDIM> >
      ele = ele_manager->getElement(d_element_type);

    /** Get Local PatchGeometry **/
    tbox::Pointer< hier::PatchGeometry<NDIM> > patch_geo = patch.getPatchGeometry();
    /** Get Local PatchTopology **/
    tbox::Pointer< hier::PatchTopology<NDIM> > patch_top = patch.getPatchTopology();
    /** 取出本地Patch的结点坐标数组 **/
    tbox::Pointer<pdat::NodeData<NDIM, double> >
      node_coord = patch_geo->getNodeCoordinates();
    /** get data of RHS **/
    tbox::Pointer<pdat::VectorData<NDIM,double> >
      vec_data = patch.getPatchData(d_carh_fval_id);

    /** get node values **/
    tbox::Pointer<pdat::VectorData<NDIM,double> >
      N_data = patch.getPatchData(d_carn_sol_cur_id);
    tbox::Pointer<pdat::VectorData<NDIM,double> >
      H_data = patch.getPatchData(d_carh_sol_cur_id);

    /** cell - entity relationship **/
    tbox::Pointer<pdat::CellData<NDIM, int> >
      cell_mat = patch.getPatchData(d_cell_mat_id);
    /** cell - heter node relationship **/
    tbox::Pointer<pdat::CellData<NDIM, int> >
      cell_hnode = patch.getPatchData(d_cell_hnode_id);
    /** heterojunction - node relationship **/
    tbox::Pointer<pdat::CellData<NDIM, int> >
      heter_node = patch.getPatchData(d_cell_hnode_id);
    tbox::Pointer<pdat::CellData<NDIM, int> >
      semiinsbc_node = patch.getPatchData(d_semiinsbc_node_id);

    /** get carh plot values **/
    tbox::Pointer<pdat::NodeData<NDIM, double> >
      carh_plot_data = patch.getPatchData(d_carh_plot_id);

    /** 获取本地单元(边, 面)周围结点的索引关系 **/
    tbox::Array<int> can_extent, can_indices;
    patch_top->getCellAdjacencyNodes(can_extent, can_indices);

    int num_cells = patch.getNumberOfCells(1);

    int * dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);
    int * dof_map_semi = d_dof_semi_info->getDOFMapping(patch, hier::EntityUtilities::NODE);

    /** material list **/
    tbox::Pointer<MaterialManager<NDIM> >
            mat_manager = MaterialManager<NDIM>::getManager();

    tbox::Pointer<Variables_Functions> Var_func = new Variables_Functions();

    vec_data->fillAll(0.0);

    /** Traversing Unit **/
    for(int i = 0; i < num_cells; ++i){

        int n_vertex = can_extent[i+1] - can_extent[i];
        tbox::Array<hier::DoubleVector<NDIM> > vertex(n_vertex);
        tbox::Array<int> node_n(n_vertex);
        tbox::Array<int> mapping(n_vertex);
        tbox::Array<int> mapping_semi(n_vertex);
        tbox::Array<double> T_node(n_vertex);
        tbox::Array<double> N_node(n_vertex);
        tbox::Array<double> H_node(n_vertex);
        tbox::Array<double> RGrate(n_vertex);
        tbox::Array<double> carh_node(n_vertex);

        /** inform for node **/
        for (int i1 = 0, j = can_extent[i]; i1 < n_vertex; ++i1, ++j) {
            mapping[i1] = dof_map[can_indices[j]];
            mapping_semi[i1] = dof_map_semi[can_indices[j]];
            for(int k = 0; k< NDIM; ++k){
              vertex[i1][k] = (*node_coord)(k, can_indices[j])*m;
            }
            node_n[i1] = can_indices[j];
            //cout<< can_indices[j]<<endl;
            carh_node[i1] = (*carh_plot_data)(0, can_indices[j])*pow(m,-3);
        }

        /** get material for cell **/
        tbox::Pointer<BaseMaterial<NDIM> > mat_cell =
                mat_manager->getMaterial(material_list[(*cell_mat)(0,i)].mat_name);

        for (int i_n = 0; i_n < n_vertex; ++i_n){
            T_node[i_n]=T;
            N_node[i_n]=(*N_data)(mapping[i_n]);
            H_node[i_n]=(*H_data)(mapping[i_n]);
            if((material_list[(*cell_mat)(0, i)].mat_type == "dielectric")&&((*semiinsbc_node)(0,node_n[i_n])==1)){
                N_node[i_n]=(*N_data)(mapping[i_n]+1);
                H_node[i_n]=(*H_data)(mapping[i_n]+1);
            }
        }

        double T_cell = Var_func->node2ElemVariable(vertex,T_node,dt,time);

        /** semi-variables for RG rate **/
        double e_trap = mat_cell->electronTrap(T_cell);
        bool fermi_st = false;

        tbox::Pointer<tbox::Vector<double> > ele_vec = new tbox::Vector<double>();
        ele_vec->resize(n_vertex);

        /** two ways to calculate the recombination & generation rate
         1. get the cell carrier density the calculate the cell rate
         2. get the node rate and then calculate the element rate **/
        for (int i_n = 0; i_n < n_vertex; ++i_n){
            double ni_tmp = mat_cell->intrinsicDensityNi(T_amb)*pow(m,-3);
            RGrate[i_n] = -1 * Var_func->recombGenerationSRH(N_node[i_n], H_node[i_n], 0.0, 0.0, T_node[i_n], kb, e, ni_tmp, e_trap, fermi_st);
        }
        /** building element vector **/
        ele->buildCVFEMRhs(vertex, RGrate, dt, time,  ele_vec);

        if(transient_state){
            //if transient state is true, the below will run, else the below will not run
            tbox::Pointer<tbox::Vector<double> > cv_vol = new tbox::Vector<double>();
            cv_vol->resize(n_vertex);
            tbox::Array<double> value(n_vertex);
            for (int i_n = 0; i_n < n_vertex; ++i_n){
                value[i_n] = 1;
            }
            ele->buildCVFEMRhs(vertex, value, dt, time, cv_vol);

            for (int i_n = 0; i_n < n_vertex; ++i_n){
                (*ele_vec)[i_n] = dt * (*ele_vec)[i_n] + (*cv_vol)[i_n] * carh_node[i_n];
            }
        }

        for(int i2= 0; i2 < n_vertex; ++i2){
            vec_data->addVectorValue(mapping[i2], (*ele_vec)[i2]);
        }
    }
}

/**
 * @brief check the node order in element
 **/
void NonlinearPoisson::checkOrder(hier::Patch<NDIM> & patch,
                const double  time,
                const double  dt,
                const string& component_name){
    /** (形函数, 积分器, 单元)管理器 **/
    tbox::Pointer<ElementManager<NDIM> >
      ele_manager = ElementManager<NDIM>::getManager();
    tbox::Pointer<BaseElement<NDIM> > ele = ele_manager->getElement(d_element_type);

    /** 取出本地PatchGeometry **/
    tbox::Pointer< hier::PatchGeometry<NDIM> > patch_geo = patch.getPatchGeometry();
    /** 取出本地PatchTopology **/
    tbox::Pointer< hier::PatchTopology<NDIM> > patch_top = patch.getPatchTopology();
    /** 取出本地Patch的结点坐标数组 **/
    tbox::Pointer<pdat::NodeData<NDIM, double> >
            node_coord = patch_geo->getNodeCoordinates();

    /** 获取本地单元(边, 面)周围结点的索引关系 **/
    tbox::Array<int> can_extent, can_indices;
    patch_top->getCellAdjacencyNodes(can_extent, can_indices);
    /** 面——单元邻接关系数据 **/
    tbox::Array<int>face_cell_ext,face_cell_idx;
    patch_top->getFaceAdjacencyCells(face_cell_ext,face_cell_idx);
    /** 单元节点邻接关系数据 **/
    tbox::Array<int>face_node_ext,face_node_idx;
    patch_top->getFaceAdjacencyNodes(face_node_ext,face_node_idx);

    //int num_faces = patch.getNumberOfFaces(1);
    int num_cells = patch.getNumberOfCells(1);

    int * dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);
    int * dof_map_semi = d_dof_semi_info->getDOFMapping(patch, hier::EntityUtilities::NODE);

    //int num_edge = cell_edge;
    //int num_ver = cell_point;
    /** Traversing Unit **/
    for(int i = 0; i < num_cells; ++i){
        int n_vertex = can_extent[i+1] - can_extent[i];
        tbox::Array<hier::DoubleVector<NDIM> > vertex(n_vertex);
        tbox::Array<int> mapping(n_vertex);
        tbox::Array<int> mapping_semi(n_vertex);
        tbox::Array<int> node_n(n_vertex);

        /** 取出单元结点坐标，以及填写映射值 **/
        for (int i1 = 0, j = can_extent[i]; i1 < n_vertex; ++i1, ++j) {
            mapping[i1] = dof_map[can_indices[j]];
            mapping_semi[i1] = dof_map_semi[can_indices[j]];
            for(int k = 0; k< NDIM; ++k){
                vertex[i1][k] = (*node_coord)(k, can_indices[j])*m;
            }
            node_n[i1] = can_indices[j];
        }

        /** 取出单元结点坐标，判断单元节点顺序是否正确 **/
        if(vertex[0][2]<vertex[4][2]||vertex[1][2]<vertex[5][2]||vertex[2][2]<vertex[6][2]||vertex[3][2]<vertex[7][2]){
            cout<<i<<" element has wrong z order"<<endl;
            for(int i=0; i<8; i++){
                cout<<vertex[i][0]<<" "<<vertex[i][1]<<" "<<vertex[i][2]<<endl;
            }
        }else if(vertex[3][1]<vertex[0][1]||vertex[2][1]<vertex[1][1]||vertex[7][1]<vertex[4][1]||vertex[6][1]<vertex[5][1]){
            cout<<i<<" element has wrong y order"<<endl;
            for(int i=0; i<8; i++){
                cout<<vertex[i][0]<<" "<<vertex[i][1]<<" "<<vertex[i][2]<<endl;
            }
        }else if(vertex[1][0]<vertex[0][0]||vertex[2][0]<vertex[3][0]||vertex[5][0]<vertex[4][0]||vertex[6][0]<vertex[7][0]){
            cout<<i<<" element has wrong x order"<<endl;
            for(int i=0; i<8; i++){
                cout<<vertex[i][0]<<" "<<vertex[i][1]<<" "<<vertex[i][2]<<endl;
            }
        }else{}
    }
}

/**
 * @brief NonlinearPoisson::transferCarhSoln: transfer current solution to plot
 * @param patch
 * @param time
 * @param dt
 * @param component_name
 */
void NonlinearPoisson::transferCarhSoln(hier::Patch<NDIM> & patch,
                                      const double  time,
                                      const double  dt,
                                      const string& component_name)
{
    tbox::Pointer<pdat::NodeData<NDIM, double> >
      plot_data = patch.getPatchData(d_carh_plot_id);
    tbox::Pointer<pdat::NodeData<NDIM, double> >
      log_plot_data = patch.getPatchData(d_logh_plot_id);
    tbox::Pointer<pdat::NodeData<NDIM, double> >
      fermi_plot_data = patch.getPatchData(d_fermih_plot_id);

    tbox::Pointer<pdat::NodeData<NDIM, double> >
      carh_sol_1_data = patch.getPatchData(d_carh_sol_1_id);
    tbox::Pointer<pdat::NodeData<NDIM, double> >
      fh_sol_1_data = patch.getPatchData(d_fh_sol_1_id);

    tbox::Pointer<pdat::VectorData<NDIM, double> >
      solution = patch.getPatchData(d_carh_sol_cur_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
      log_solution = patch.getPatchData(d_logh_sol_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
      fermi_solution = patch.getPatchData(d_fermih_sol_id);

    tbox::Pointer<pdat::NodeData<NDIM, int> >
      semi_node = patch.getPatchData(d_semi_node_id);
    tbox::Pointer<pdat::NodeData<NDIM, int> >
      sio2_node = patch.getPatchData(d_sio2_node_id);

    int * dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);

    int num_nodes = patch.getNumberOfNodes(1);

    for(int i = 0; i < num_nodes; ++i){
        if(((*semi_node)(0,i) == 1)||((*sio2_node)(0,i)== 1)){
            int index = dof_map[i];
            plot_data->getPointer()[i] = solution->getPointer()[index]/pow(m,-3);
            log_plot_data->getPointer()[i] = log_solution->getPointer()[index];
            fermi_plot_data->getPointer()[i] = fermi_solution->getPointer()[index];

            (*carh_sol_1_data)(0,i) = (*solution)(index);
            (*fh_sol_1_data)(0,i) = (*fermi_solution)(index);

        }
    }
}


/**
 * @brief NonlinearPoisson::getIV: I-V data
 * @param patch
 * @param time
 * @param dt
 * @param component_name
 */
void NonlinearPoisson::getIV(hier::Patch<NDIM> &patch,
                                      const double time,
                                      const double dt,
                                      const string & component_name){
    /** (形函数, 积分器, 单元)管理器 **/
    tbox::Pointer<ElementManager<NDIM> >
      ele_manager = ElementManager<NDIM>::getManager();
    tbox::Pointer<BaseElement<NDIM> >
      ele = ele_manager->getElement(d_element_type);

    /** 取出本地PatchGeometry **/
    tbox::Pointer< hier::PatchGeometry<NDIM> > patch_geo = patch.getPatchGeometry();
    /** 取出本地PatchTopology **/
    tbox::Pointer< hier::PatchTopology<NDIM> > patch_top = patch.getPatchTopology();
    /** 取出本地Patch的结点坐标数组 **/
    tbox::Pointer<pdat::NodeData<NDIM, double> >
      node_coord = patch_geo->getNodeCoordinates();

    /** get node values **/
    tbox::Pointer<pdat::VectorData<NDIM, double> >
      V_soln = patch.getPatchData(d_poisson_sol_cur_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
      N_soln = patch.getPatchData(d_carn_sol_cur_id);
    tbox::Pointer<pdat::VectorData<NDIM, double> >
      H_soln = patch.getPatchData(d_carh_sol_cur_id);

    /** 取出本地cotact id值 **/
    tbox::Pointer<pdat::NodeData<NDIM, int> >
      contact_node_data = patch.getPatchData(d_contact_node_id);

    /** current density vector **/
    tbox::Pointer<pdat::CellData<NDIM,double> >
      J_data = patch.getPatchData(d_J_vec_id);

    tbox::Pointer<Variables_Functions>
            Var_func = new Variables_Functions();

    /** 获取本地单元(边, 面)周围结点的索引关系 **/
    tbox::Array<int> can_extent, can_indices;
    patch_top->getCellAdjacencyNodes(can_extent, can_indices);

    /** setup the source target **/
    tbox::Pointer<bound_condition> bound_source = new bound_condition();

    //int num_faces = patch.getNumberOfFaces(0);
    int num_cells = patch.getNumberOfCells(0);

    /** 取出自由度映射信息 **/
    int * dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);
    int * dof_map_semi = d_dof_semi_info->getDOFMapping(patch, hier::EntityUtilities::NODE);

    /** boundary condition source **/
    int num_bc = EBoundassignment.getSize();

    /** cell - entity relationship **/
    tbox::Pointer<pdat::CellData<NDIM, int> >
            cell_mat = patch.getPatchData(d_cell_mat_id);

    /** material list **/
    tbox::Pointer<MaterialManager<NDIM> >
            mat_manager = MaterialManager<NDIM>::getManager();

    tbox::Array<double> S_coe;
    S_coe.resizeArray(num_bc);
    tbox::Array<double> volt_coe;
    volt_coe.resizeArray(num_bc);
    tbox::Array<double> current_coe;
    current_coe.resizeArray(num_bc);
    for(int i_bc=0; i_bc < num_bc; ++i_bc){
        current_coe[i_bc] = 0;
        S_coe[i_bc] = 0;
    }

    for(int i_bc = 0; i_bc < num_bc; ++i_bc){
        /** load fixed volt value **/
        double app_volt = bound_source->source_assign(EBoundassignment[i_bc].sc_type,
                                                      time, dt,
                                                      EBoundassignment[i_bc].sc_coe);
        volt_coe[i_bc] = app_volt;  //get the app_volt
    }

    for(int i=0; i<num_cells; i++){
        int n_vertex = can_extent[i+1] - can_extent[i];
        tbox::Array<hier::DoubleVector<NDIM> > vertex(n_vertex);
        tbox::Array<int> mapping(n_vertex);
        tbox::Array<int> mapping_semi(n_vertex);
        tbox::Array<int> node_n(n_vertex);
        tbox::Array<int> contact_node(n_vertex);
        tbox::Array<double> T_node(n_vertex);
        tbox::Array<double> V_node(n_vertex);
        tbox::Array<double> N_node(n_vertex);
        tbox::Array<double> H_node(n_vertex);

        /** get cell node coord **/
        for (int i1 = 0, j = can_extent[i]; i1 < n_vertex; ++i1, ++j) {
            mapping[i1] = dof_map[can_indices[j]];
            mapping_semi[i1] = dof_map_semi[can_indices[j]];
            for(int k = 0; k< NDIM; ++k){
                vertex[i1][k] = (*node_coord)(k, can_indices[j])*m;
            }
            node_n[i1] = can_indices[j];
            contact_node[i1] = (*contact_node_data)(0, node_n[i1]);
        }

        for (int i_n = 0; i_n < n_vertex; ++i_n){
            T_node[i_n]=T;
            V_node[i_n]=(*V_soln)(mapping_semi[i_n]);
            N_node[i_n]=(*N_soln)(mapping[i_n]);
            H_node[i_n]=(*H_soln)(mapping[i_n]);
        }

        /// 初始化contact面单元的面积和外法向量
        tbox::Pointer<tbox::Array<double> > norm_vec = new tbox::Array<double>(3);
        tbox::Pointer<tbox::Array<double> > face_area = new tbox::Array<double>(1);
        (*face_area)[0]=0.0;
        for (int i = 0; i < 3; ++i) {
          (*norm_vec)[i] = 0.0;
        }
    }

    #if 0
    for(int i_bc = 0; i_bc < num_bc; ++i_bc){
        /** output the current and voltage information **/
        cout<<"volt"<<i_bc<<":"<<volt_coe[i_bc]<<endl;
        cout<<"current"<<i_bc<<":"<<current_coe[i_bc]<<endl;
        cout<<"S"<<i_bc<<":"<<S_coe[i_bc]<<endl;
        Var_func->DataOutput("current",i_bc,current_coe[i_bc]);
        Var_func->DataOutput("voltage",i_bc,volt_coe[i_bc]);
    }
    #endif
}


