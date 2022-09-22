//
// 文件名:     NonlinearPoissonLevel.C
// 软件包:     JAUMIN 
// 版权　:     北京应用物理与计算数学研究所 
// 版本号:     $Revision: 0 $
// 修改　:     $Date: Mon Apr 13 11:08:15 2015 $
// 描述　:     求解非线性Poisson方程的网格层算法类
// 类别　:     %Internal File% ( Don't delete this line )
//
#define USE_JAUMIN

#include "Pointer.h"
#include "LinearSolverManager.h"
#include "NonlinearPoissonLevel.h" 
#include "DOFInfo.h"
#include "NonlinearPoisson.h"
#include "GridGeometry.h"
#include "Global_variable.h"
#include "Bound_condition.h"
#include "Variables_Functions.h"

#include <Utilities.h>

//#include "poisson_function.h"
//#include "nssolver.h"
#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#include<fstream>
#endif

/*************************************************************************
 * 构造函数.
 *************************************************************************/
NonlinearPoissonLevel::NonlinearPoissonLevel(const string& object_name,
    tbox::Pointer<NonlinearPoisson> strategy,
    tbox::Pointer<tbox::Database> input_db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
  assert(!object_name.empty());
  assert(strategy.getPointer() != NULL);
#endif
  /** 创建网格片积分算法 **/
  d_patch_strategy = strategy;
  d_object_name = object_name;

  /** create linear solver **/
  d_linear_solver_db = input_db->getDatabase ("LinearSolver");
  d_linear_solver_manager = solv::LinearSolverManager<NDIM>::getManager();

  d_linear_solver_LU =  d_linear_solver_manager->lookupLinearSolver(
              d_linear_solver_db->getDatabase("linearSolver_LU")->getString("solver_name"));
  d_linear_solver_PCG  =  d_linear_solver_manager->lookupLinearSolver(
              d_linear_solver_db->getDatabase("linearSolver_PCG")->getString("solver_name"));

  /** nonlinear iteration parameters **/
  tol_cp = input_db->getDatabase("Iteration_param")->getDouble("tol_cp");
  tol_non = input_db->getDatabase("Iteration_param")->getDouble("tol_non");
  tol_semi = input_db->getDatabase("Iteration_param")->getDouble("tol_semi");
  tol_poi = input_db->getDatabase("Iteration_param")->getDouble("tol_poi");
  error_cp = input_db->getDatabase("Iteration_param")->getDouble("error_cp");
  error_non = input_db->getDatabase("Iteration_param")->getDouble("error_non");
  error_semi = input_db->getDatabase("Iteration_param")->getDouble("error_semi");
  error_poi = input_db->getDatabase("Iteration_param")->getDouble("error_poi");
  iter_cp_max = input_db->getDatabase("Iteration_param")->getDouble("iter_cp_max");
  iter_non_max = input_db->getDatabase("Iteration_param")->getDouble("iter_non_max");
  iter_semi_max = input_db->getDatabase("Iteration_param")->getDouble("iter_semi_max");
  iter_poi_max = input_db->getDatabase("Iteration_param")->getDouble("iter_poi_max");
  enable_selfDop = input_db->getDatabase("selfDop")->getBool("enable_selfDop");
  save_bas = input_db->getBool("save_bas");
  enable_tid = input_db->getBool("enable_tid");
  get_InitialSol = input_db->getBool("get_InitialSol");
  debug_flag = input_db->getString("debug_flag");
  save_PoissonSolcur = input_db->getBool("save_PoissonSolcur");
  get_PoissonSolcur = input_db->getBool("get_PoissonSolcur");

  /** read boundary condition information **/
  tbox::Pointer<bound_condition> bc_manager = new bound_condition();
  EBoundassignment = bc_manager->getEBCassignment(input_db);
  num_bc = EBoundassignment.getSize();

}

/** ***********************************************************************
 * 析构函数.
 ********************************************************************** **/
NonlinearPoissonLevel::~NonlinearPoissonLevel() {
}

/** ***********************************************************************
 * 初始化该积分算法: 创建所有计算需要的积分构件.
 *
 * 该函数创建了1个内存构件, 1个初始化构件, 4个数值构件. 
 * 这些构件所操作的数据片, 
 * 由函数 d_patch_strategy->initializeComponent() 指定.
 *
 *********************************************************************** **/
void NonlinearPoissonLevel::initializeLevelIntegrator(
        tbox::Pointer<algs::IntegratorComponentManager<NDIM> > manager)
{
  /** 初值构件: 管理当前值数据片的内存以及初始化 **/
  d_init_intc = new algs::InitializeIntegratorComponent<NDIM>("INIT", 
                                                              d_patch_strategy,
                                                              manager);

  d_dop_intc = new algs::InitializeIntegratorComponent<NDIM>("DOP_INIT",
                                                                    d_patch_strategy,
                                                                    manager);

  d_tid_intc = new algs::InitializeIntegratorComponent<NDIM>("TID_INIT",
                                                                    d_patch_strategy,
                                                                    manager);

  d_init_dof_intc = new algs::InitializeIntegratorComponent<NDIM>("DOF_INIT",
                                                                    d_patch_strategy,
                                                                    manager);
  /**8.11**/

  d_init_sol_intc = new algs::InitializeIntegratorComponent<NDIM>("INIT_SOL",
                                                              d_patch_strategy,
                                                              manager);

  /** initial solution vector **/
  d_solvec_init_intc = new algs::NumericalIntegratorComponent<NDIM>("SOLINIT",
                                                                    d_patch_strategy,
                                                                    manager);
  /** check the node order in element**/
  d_checkorder_intc = new algs::NumericalIntegratorComponent<NDIM>("CHECKORDER",
                                                                   d_patch_strategy,
                                                                   manager);
  /** 内存构件: alloc/dealloc memory for vectors **/
  d_alloc_vec = new algs::MemoryIntegratorComponent<NDIM>("AllocVec",
                                                           d_patch_strategy,
                                                           manager);
  /** 内存构件: alloc/dealloc memory for matrixs **/
  d_alloc_CCmat = new algs::MemoryIntegratorComponent<NDIM>("AllocCCMat",
                                                            d_patch_strategy,
                                                            manager);

  d_alloc_Pmat = new algs::MemoryIntegratorComponent<NDIM>("AllocPMat",
                                                           d_patch_strategy,
                                                           manager);

  d_alloc_Newton_Pmat = new algs::MemoryIntegratorComponent<NDIM>("AllocNewtonPMat",
                                                           d_patch_strategy,
                                                           manager);

  d_alloc_Cmat = new algs::MemoryIntegratorComponent<NDIM>("AllocCMat",
                                                           d_patch_strategy,
                                                           manager);
  /** 内存构件: alloc/dealloc memory for right side vectors **/
  d_alloc_CCrhs = new algs::MemoryIntegratorComponent<NDIM>("AllocCCRhs",
                                                            d_patch_strategy,
                                                            manager);
  d_alloc_Prhs = new algs::MemoryIntegratorComponent<NDIM>("AllocPRhs",
                                                           d_patch_strategy,
                                                           manager);
  d_alloc_Newton_Prhs = new algs::MemoryIntegratorComponent<NDIM>("AllocNewtonPRhs",
                                                           d_patch_strategy,
                                                           manager);

  d_alloc_Crhs = new algs::MemoryIntegratorComponent<NDIM>("AllocCRhs",
                                                           d_patch_strategy,
                                                           manager);

  /** 数值构件 assign boundary condition to solution vector **/
  d_CCapp_bc  = new algs::NumericalIntegratorComponent<NDIM>("CCAPPBC",
                                                             d_patch_strategy,
                                                             manager);
  d_Papp_bc   = new algs::NumericalIntegratorComponent<NDIM>("PAppBC",
                                                            d_patch_strategy,
                                                            manager);
  d_Newton_Papp_bc   = new algs::NumericalIntegratorComponent<NDIM>("NewtonPAppBC",
                                                            d_patch_strategy,
                                                            manager);
  d_Napp_bc   = new algs::NumericalIntegratorComponent<NDIM>("NAppBC",
                                                            d_patch_strategy,
                                                            manager);
  d_Happ_bc   = new algs::NumericalIntegratorComponent<NDIM>("HAppBC",
                                                            d_patch_strategy,
                                                            manager);
  /** 数值构件 construct of mass matrix **/
  d_CCmat_intc = new algs::NumericalIntegratorComponent<NDIM>("CCMAT",
                                                            d_patch_strategy,
                                                            manager);
  d_Pmat_intc = new algs::NumericalIntegratorComponent<NDIM>("PMAT",
                                                            d_patch_strategy,
                                                            manager);
  d_Newton_Pmat_intc = new algs::NumericalIntegratorComponent<NDIM>("NewtonPMAT",
                                                            d_patch_strategy,
                                                            manager);
  d_adtest_intc = new algs::NumericalIntegratorComponent<NDIM>("ADTEST",
                                                            d_patch_strategy,
                                                            manager);
  d_CVmat_intc = new algs::NumericalIntegratorComponent<NDIM>("CVMAT",
                                                            d_patch_strategy,
                                                            manager);
  /** 数值构件 construction of right side vectors **/
  d_CCfx_intc = new algs::NumericalIntegratorComponent<NDIM>("CCFX",
                                                           d_patch_strategy,
                                                           manager);
  d_Pfx_intc = new algs::NumericalIntegratorComponent<NDIM>("PFX",
                                                           d_patch_strategy,
                                                           manager);
  d_Newton_Pfx_intc = new algs::NumericalIntegratorComponent<NDIM>("NewtonPFX",
                                                           d_patch_strategy,
                                                           manager);
  d_Nfx_intc = new algs::NumericalIntegratorComponent<NDIM>("NFX",
                                                             d_patch_strategy,
                                                             manager);
  d_Hfx_intc = new algs::NumericalIntegratorComponent<NDIM>("HFX",
                                                             d_patch_strategy,
                                                             manager);
  /** 数值构件 assign boundary to matrix **/
  d_CCphybc_intc = new algs::NumericalIntegratorComponent<NDIM>("CCPHYBC",
                                                              d_patch_strategy,
                                                              manager);
  d_Pphybc_intc = new algs::NumericalIntegratorComponent<NDIM>("PPHYBC",
                                                              d_patch_strategy,
                                                              manager);
  d_Newton_Pphybc_intc = new algs::NumericalIntegratorComponent<NDIM>("NewtonPPHYBC",
                                                              d_patch_strategy,
                                                              manager);
  d_Nphybc_intc = new algs::NumericalIntegratorComponent<NDIM>("NPHYBC",
                                                              d_patch_strategy,
                                                              manager);
  d_Hphybc_intc = new algs::NumericalIntegratorComponent<NDIM>("HPHYBC",
                                                              d_patch_strategy,
                                                              manager);
  /** 数值构件 Joule heating, E field, Current density, Fermi level
      Fermi intrinsic && Temperature dependent semi variables **/
  d_Cdens_intc = new algs::NumericalIntegratorComponent<NDIM>("CURRDENS",
		                                                      d_patch_strategy,
															  manager);
  d_Flevel_intc = new algs::NumericalIntegratorComponent<NDIM>("FLEVEL",
		                                                      d_patch_strategy,
															  manager);
  d_NPconstraint_intc = new algs::NumericalIntegratorComponent<NDIM>("NPCONS",
                                                                     d_patch_strategy,
                                                                     manager);
  d_Fintrin_intc= new algs::NumericalIntegratorComponent<NDIM>("FINTRIN",
                                                              d_patch_strategy,
		                                                      manager);
  d_equilibrium_intc = new algs::NumericalIntegratorComponent<NDIM>("EQUILIBRIUM",
                                                                    d_patch_strategy,
                                                                    manager);

  /** 数值构件 solution post process for step (output) **/
  d_sol_process_intc = new algs::NumericalIntegratorComponent<NDIM>("SOLPROCESS",
                                                                    d_patch_strategy,
                                                                    manager);
  /** 数值构件 solution transport for draw figure **/
  d_Ppost_intc = new algs::NumericalIntegratorComponent<NDIM>("PPOST",
                                                            d_patch_strategy,
                                                            manager);
  d_Npost_intc = new algs::NumericalIntegratorComponent<NDIM>("NPOST",
                                                            d_patch_strategy,
                                                            manager);
  d_Hpost_intc = new algs::NumericalIntegratorComponent<NDIM>("HPOST",
                                                            d_patch_strategy,
                                                            manager);
  d_E_intc = new algs::NumericalIntegratorComponent<NDIM>("calculate_E",
                                                           d_patch_strategy,
                                                           manager);
  d_J_intc = new algs::NumericalIntegratorComponent<NDIM>("calculate_J",
                                                           d_patch_strategy,
                                                           manager);
  d_GetIV_intc = new algs::NumericalIntegratorComponent<NDIM>("get_IV",
                                                               d_patch_strategy,
                                                               manager);

  /** 数值构件 solution transport for nonlinear iteration **/
  d_Psolution_intc = new algs::NumericalIntegratorComponent<NDIM>("PSOLTRANS",
                                                            d_patch_strategy,
                                                            manager);
  d_GummelPsolution_intc = new algs::NumericalIntegratorComponent<NDIM>("GummelPSOLTRANS",
                                                            d_patch_strategy,
                                                            manager);
  d_Semisolution_intc = new algs::NumericalIntegratorComponent<NDIM>("SEMISOLTRANS",
                                                              d_patch_strategy,
                                                              manager);
  d_GummelSemisolution_intc = new algs::NumericalIntegratorComponent<NDIM>("GummelSEMISOLTRANS",
                                                              d_patch_strategy,
                                                              manager);
  /** 归约构件 searching the maximum change in solutions **/
  d_Dopreduction_intc = new algs::ReductionIntegratorComponent<NDIM>("MAX_DOP",
                                                                   MPI_MAX,
                                                                   d_patch_strategy,
                                                                   manager);
  d_Preduction_intc = new algs::ReductionIntegratorComponent<NDIM>("POISSON",
                                                                   MPI_MAX,
                                                                   d_patch_strategy,
                                                                   manager);
  d_Poireduction_intc = new algs::ReductionIntegratorComponent<NDIM>("GummelPOISSON",
                                                                   MPI_MAX,
                                                                   d_patch_strategy,
                                                                   manager);
  d_Smeireduction_intc = new algs::ReductionIntegratorComponent<NDIM>("SEMI",
                                                                      MPI_MAX,
                                                                      d_patch_strategy,
                                                                      manager);
  d_GummelSmeireduction_intc = new algs::ReductionIntegratorComponent<NDIM>("GummelSEMI",
                                                                      MPI_MAX,
                                                                      d_patch_strategy,
                                                                      manager);
  d_getI_intc = new algs::ReductionIntegratorComponent<NDIM>("get_I",
                                                             MPI_SUM,
                                                             d_patch_strategy,
                                                             manager);
  d_saveSol_intc = new algs::ReductionIntegratorComponent<NDIM>("saveSol",
                                                                      MPI_MAX,
                                                                      d_patch_strategy,
                                                                      manager);
  /** 复制构件: 接受数值解，把数值解从新值数据片复制到当前数据片,该类型数据片只能有一个 **/
  d_rebalance_intc = new algs::RebalanceIntegratorComponent<NDIM>(
      "INIT_REBALANCE", d_patch_strategy, manager);

  d_dt_intc = new algs::DtIntegratorComponent<NDIM>("Dt", d_patch_strategy, manager);

}

void NonlinearPoissonLevel::initializeLevelData(const tbox::Pointer< hier::BasePatchLevel<NDIM> > level, 
                                                const double init_data_time,
                                                const bool   initial_time)
{
    d_init_sol_intc->initializeLevelData(level,
                                         init_data_time,
                                         initial_time);
    if(!enable_selfDop){
        /** 通过在GriddingAlgorithm中设置"init_on_proc_zero=TRUE",
            读入的网格只在0号进程形成1个网格片，因此数据的初始化是串行进行的 **/
        d_init_intc->initializeLevelData(level,
                                       init_data_time,
                                       initial_time);
    }
    if(!enable_selfDop||get_InitialSol||get_PoissonSolcur){
        /**  将网格及数据一并并行分发 **/
        d_rebalance_intc->rebalance(level);
    }

    if(enable_selfDop){
        d_dop_intc->initializeLevelData(level,
                                                 init_data_time,
                                                 initial_time);
    }

    d_init_dof_intc->initializeLevelData(level,
                                         init_data_time,
                                         initial_time);

    if(enable_tid){
        d_tid_intc->initializeLevelData(level,
                                             init_data_time,
                                             initial_time);
    }
}

/** ***********************************************************************
 * 获取网格层上的计算时间步长.
 *
 * @note: 求解问题如果是静态的, 则该函数将不被调用, 如果是时间发展的, 则该函数调
 * 用有限元构件的函数getLevelDt(),该构件又进一步自动调用
 * d_fem_strategy->getCellDt(), 完成时间步长的求解.
 *
 ********************************************************************** **/
double NonlinearPoissonLevel::getLevelDt(
    const tbox::Pointer<hier::BasePatchLevel<NDIM> > level,
    const double dt_time, const bool initial_time, const int flag_last_dt,
    const double last_dt) {
#ifdef DEBUG_CHECK_ASSERTIONS
  assert(!level.isNull());
#endif
  return (d_dt_intc->getLevelDt(level, dt_time, initial_time, flag_last_dt,
                                last_dt, false));
}

/** ***********************************************************************
 * 向前积分一个时间步.
 *
 * 注解: 该函数调用了有限元构件对象的5个函数，
 * 分别完成计算矩阵右端项, 设置边界条件, 求解线性系统的计算, 误差计算.
 *
 ********************************************************************** **/
int NonlinearPoissonLevel::advanceLevel(
              const tbox::Pointer< hier::BasePatchLevel<NDIM> > level, 
              const double current_time,
              const double predict_dt,
              const double max_dt,
              const double min_dt,
              const bool   first_step,
              const int    step_number,
              double&      actual_dt)
{
#ifdef DEBUG_CHECK_ASSERTIONS
  assert(!level.isNull());
#endif

  actual_dt = predict_dt;

  const tbox::Pointer< hier::PatchLevel<NDIM> > patch_level = level;

  /** memory assignment **/
  d_alloc_vec->allocatePatchData(patch_level, current_time);

  #if 0
  /** check the node order in element**/
  tbox::pout << "---------check the node order in element---------"<< std::endl
             << "-------------------------------------------------"<< std::endl;
  d_checkorder_intc->computing(patch_level,current_time, actual_dt);
  #endif

  /** ***********initialization - calculate the linear poisson equation*********** **/
  tbox::pout << "-----start initialization of electro-poisson-----"<< std::endl
             << "-------------------------------------------------"<< std::endl;

  if(first_step){
      double MAX_DOP;
      d_Dopreduction_intc->reduction(&MAX_DOP, 1, patch_level, current_time, actual_dt, true);
      //cout<<"MAX_DOP:"<<MAX_DOP<<endl;
      d_patch_strategy->parameter.setmaxdop(MAX_DOP);
  }

  /** initial solution vector **/
  d_solvec_init_intc->computing(patch_level,current_time, actual_dt);

  /** ********************** main loop ********************** **/
  tbox::pout << "***** start electro-thermal nonlinear process *****"<< std::endl;

  int iter_num_cp, iter_num_semi, iter_num_poi;

  iter_num_cp = 0;
  error_cp = 1.0;

  /** number of iterations **/
  // iter_num_cp += 1;

  /** ********** loop of Poisson-Transport equations *************** **/
  tbox::pout << "***** start semi_iteration of electron-transport! *****"<< std::endl;

  iter_num_semi = 0; error_semi = 1.0;

  if(debug_flag=="debug0"){
      iter_poi_max = 0;
      iter_semi_max = 0;
  }if(debug_flag=="debug1"){
      iter_poi_max = 1;
      iter_semi_max = 1;
  }

  if(current_time<-100){
      while(error_semi > tol_semi && iter_num_semi < iter_semi_max){//iter_semi_max
          //cout<<"dt:"<<actual_dt<<endl;

          /** number of iterations **/
          iter_num_semi += 1;

          if(tbox::MPI::getRank()==0) tbox::pout<<"iter_num_semi = "<<iter_num_semi<<endl;

          tbox::pout << "***** start poisson_iteration of electro-poisson *****"<< std::endl;

          iter_num_poi = 0;
          error_poi = 1.0;

          while(error_poi > tol_poi && iter_num_poi < iter_poi_max){//iter_poi_max
              if(get_InitialSol&&iter_num_semi==1){
                  break;
              }

              /** number of iterations **/
              iter_num_poi +=1;

              //Gummel interative method
              //cout<<"allocatePathData"<<endl;
              d_alloc_Pmat->allocatePatchData(patch_level, current_time);
              d_alloc_Prhs->allocatePatchData(patch_level, current_time);
              //cout<<"end of allocatePathData"<<endl;
              d_Pmat_intc->computing(patch_level, current_time, actual_dt);
              d_Pfx_intc->computing(patch_level, current_time, actual_dt);
              d_Pphybc_intc->computing(patch_level, current_time, actual_dt);
              //d_Papp_bc->computing(patch_level, current_time, actual_dt);
              int P_matrix_id = d_patch_strategy->getPoissonMatrixID();
              int P_fx_id = d_patch_strategy->getPoissonFxID();
              int P_sol_id = d_patch_strategy->getPoissonSolID();
              d_linear_solver_LU->setMatrix(P_matrix_id);
              d_linear_solver_LU->setRHS(P_fx_id);
              d_linear_solver_LU->solve(first_step, P_sol_id, patch_level, d_linear_solver_db->getDatabase("linearSolver_LU"));
              MatrixEquationOutput(patch_level,P_matrix_id,P_fx_id,"poi_matrix","poi_fx");
              double poisson_error[1];
              d_Poireduction_intc->reduction(&poisson_error[0], 1, patch_level, current_time, actual_dt, true);
              error_poi = poisson_error[0];
              d_alloc_Pmat->deallocatePatchData(patch_level);
              d_alloc_Prhs->deallocatePatchData(patch_level);
              d_GummelPsolution_intc->computing(patch_level,current_time, actual_dt);
              if(tbox::MPI::getRank()==0){
                  tbox::pout<<"step: "<<iter_num_poi<<endl;
                  tbox::pout<<"potential_error: "<<error_poi<<endl;
              }
          }
          d_alloc_Cmat->allocatePatchData(patch_level, current_time);
          d_alloc_Crhs->allocatePatchData(patch_level, current_time);
          d_CVmat_intc->computing(patch_level, current_time, actual_dt);
          d_Nfx_intc-> computing(patch_level,current_time, actual_dt);
          d_Nphybc_intc-> computing(patch_level,current_time, actual_dt);
          //d_Hfx_intc-> computing(patch_level,current_time, actual_dt);
          //d_Napp_bc->computing(patch_level,current_time, actual_dt);
          //d_Hphybc_intc-> computing(patch_level,current_time, actual_dt);
          //d_Happ_bc->computing(patch_level,current_time, actual_dt);
          int n_mat_id = d_patch_strategy->getCarnMatrixID();
          int n_fx_id = d_patch_strategy->getCarnFxID();
          int n_sol_id = d_patch_strategy->getCarnSolID();
          d_linear_solver_LU->setMatrix(n_mat_id);
          d_linear_solver_LU->setRHS(n_fx_id);
          d_linear_solver_LU->solve(first_step, n_sol_id, patch_level, d_linear_solver_db->getDatabase("linearSolver_LU"));
          int h_mat_id = d_patch_strategy->getCarhMatrixID();
          int h_fx_id = d_patch_strategy->getCarhFxID();
          int h_sol_id = d_patch_strategy->getCarhSolID();
          d_linear_solver_LU->setMatrix(h_mat_id);
          d_linear_solver_LU->setRHS(h_fx_id);
          d_linear_solver_LU->solve(first_step, h_sol_id, patch_level, d_linear_solver_db->getDatabase("linearSolver_LU"));
          #if 1
          MatrixEquationOutput(patch_level,n_mat_id,n_fx_id,"n_matrix","n_fx");
          MatrixEquationOutput(patch_level,h_mat_id,h_fx_id,"h_matrix","h_fx");
          #endif
          d_alloc_Cmat->deallocatePatchData(patch_level);
          d_alloc_Crhs->deallocatePatchData(patch_level);
          d_NPconstraint_intc->computing(patch_level,current_time, actual_dt);
          d_GummelSemisolution_intc->computing(patch_level,current_time, actual_dt); // poi = n,p
          d_Flevel_intc->computing(patch_level,current_time, actual_dt); //n,p = poi
          double gummel_semi_error[1];
          d_GummelSmeireduction_intc->reduction(&gummel_semi_error[0], 1, patch_level, current_time, actual_dt, true);
          error_semi = gummel_semi_error[0];
          tbox::pout<<"concentration_error: "<<error_semi<<endl;
          d_Semisolution_intc->computing(patch_level,current_time, actual_dt);
          //fermin_sol_scr fermih_sol_scr carn_stetch_data carh_stetch_data poisson_sol_scr
          error_Output("error_poi", iter_num_poi, error_poi,current_time, actual_dt);
          error_Output("error_semi", iter_num_poi, error_semi,current_time, actual_dt);
      }
  }else{
      while(error_semi > tol_semi && iter_num_semi < iter_semi_max){ //iter_semi_max
          //cout<<"dt:"<<actual_dt<<endl;

          /** number of iterations **/
          iter_num_semi += 1;

          if(tbox::MPI::getRank()==0) tbox::pout<<"iter_num_semi = "<<iter_num_semi<<endl;

          tbox::pout << "***** start poisson_iteration of electro-poisson *****"<< std::endl;

          iter_num_poi = 0;
          error_poi = 1.0;

          while(error_poi > tol_poi && iter_num_poi < iter_poi_max){ //iter_poi_max
              /** number of iterations **/
              iter_num_poi +=1;
              #if 1
              //Newton interative method
              /** allocate Poisson stiff matrix (In the Newton method the matrix and rhs should updated every step) **/
              d_alloc_Newton_Pmat->allocatePatchData(patch_level, current_time);
              d_alloc_Newton_Prhs->allocatePatchData(patch_level, current_time);

              //d_GetIV_intc->computing(patch_level, current_time, actual_dt);

              /** setup poisson matrix **/
              //d_adtest_intc->computing(patch_level, current_time, actual_dt);
              d_Newton_Pmat_intc->computing(patch_level, current_time, actual_dt);
              /** computing Poisson RHS **/
              //d_Newton_Pfx_intc->computing(patch_level, current_time, actual_dt);
              /** assign boundary condtion **/
              d_Newton_Pphybc_intc-> computing(patch_level, current_time, actual_dt);
              //d_Newton_Papp_bc->computing(patch_level, current_time, actual_dt);

              /** get id for matrix, right side vector, and solution vector successively **/
              int matrix_id = d_patch_strategy->getPoissonNewtonMatrixID();
              int fx_id = d_patch_strategy->getPoissonNewtonFxID();
              int sol_id = d_patch_strategy->getPoissonNewtonSolID();

              #if 1
              MatrixEquationOutput(patch_level,matrix_id,fx_id,"newton_matrix","newton_fx");
              #endif

              /** setup solver **/
              d_linear_solver_LU->setMatrix(matrix_id);
              d_linear_solver_LU->setRHS(fx_id);
              d_linear_solver_LU->solve(first_step, sol_id, patch_level, d_linear_solver_db->getDatabase("linearSolver_LU"));

              /** calculate the maximum error **/
              double poi_error[1];
              d_Preduction_intc->reduction(&poi_error[0], 1, patch_level, current_time, actual_dt, true);
              error_poi = poi_error[0];

              /** finding maximum potential variation (based on the scaled variable) **/
              double semi_error[1];
              d_Smeireduction_intc->reduction(&semi_error[0], 1, patch_level, current_time, actual_dt, true);
              error_semi = semi_error[0];

              /** data output **/
              tbox::pout<<"step: "<<iter_num_poi<<endl;
              tbox::pout<<"potential_error: "<<error_poi<<endl;
              tbox::pout<<"concentration_error: "<<error_semi<<endl;

              /** deallocate memory for poisson **/
              d_alloc_Newton_Pmat->deallocatePatchData(patch_level);
              d_alloc_Newton_Prhs->deallocatePatchData(patch_level);

              /** if(iter_num_poi == 1){
                  error_poi = 0.0;
              } **/

              /** solution transport for poisson equation (potential) **/
              d_Psolution_intc->computing(patch_level,current_time, actual_dt);

              /** update quasi-fermi level based on the updated carrier density **/
              d_Flevel_intc->computing(patch_level,current_time, actual_dt);

              /** apply the constraint on the updated carrier density **/
              d_NPconstraint_intc->computing(patch_level,current_time, actual_dt);

              /** solution transport for semiconductor equations **/
              d_Semisolution_intc->computing(patch_level,current_time, actual_dt);
              //fermin_sol_scr fermih_sol_scr carn_stetch_data carh_stetch_data poisson_sol_scr

              /** if(iter_num_semi == 1){
                  error_semi = 0.0;
              } **/
              error_Output("error_poi", iter_num_poi, error_poi,current_time, actual_dt);
              error_Output("error_semi", iter_num_poi, error_semi,current_time, actual_dt);
              #endif
          }
      }
  }

  /** ************save PoissonSolcur************* **/
  if(save_PoissonSolcur){
      //int PoissonSolcur_id = d_patch_strategy->getPoissonSolcurID();
      //stringstream PoissonSolcur_filename;
      //PoissonSolcur_filename<<"PoissonSolcur_"<<int(current_time/actual_dt);
      //SolutionOutput(patch_level,PoissonSolcur_id,PoissonSolcur_filename.str());
      //double flag[num_bc];
      //d_saveSol_intc->reduction(flag, num_bc, patch_level, current_time, actual_dt);
  }
  /** ************后处理计算 (计算数值解与解析解误差)************* **/
  d_sol_process_intc->computing(patch_level,current_time, actual_dt); //logn_sol_data logh_sol_data
  d_Ppost_intc->computing(patch_level, current_time, actual_dt, false); //plot_data sol_1_data (possion)
  d_Npost_intc->computing(patch_level, current_time, actual_dt, false); //carn_sol_1_data fn_sol_1_data
  //d_Hpost_intc->computing(patch_level, current_time, actual_dt, false); //carh_sol_1_data fh_sol_1_data
  d_E_intc->computing(patch_level, current_time, actual_dt); //E data
  d_J_intc->computing(patch_level, current_time, actual_dt); //J_data
  //d_GetIV_intc->computing(patch_level, current_time, actual_dt); //get I-V data
  double I[num_bc];
  d_getI_intc->reduction(I, num_bc, patch_level, current_time, actual_dt);
  IV_Output(I, num_bc, EBoundassignment, current_time, actual_dt);
  d_alloc_vec->deallocatePatchData(patch_level);
  return(0);
}





