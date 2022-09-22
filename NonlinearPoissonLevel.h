//
// 文件名:     NonlinearPoissonLevel.h
// 软件包:     JAUMIN 
// 版权　:     北京应用物理与计算数学研究所 
// 版本号:     $Revision: 0 $
// 修改　:     $Date: Sun Apr 12 11:14:26 2015 $
// 描述　:     求解非线性Poisson方程的网格层算法类
// 类别　:     %Internal File% ( Don't delete this line )
//


#ifndef included_NonlinearPoissonLevel
#define included_NonlinearPoissonLevel

#include "TimeIntegratorLevelStrategy.h"
#include "PatchLevel.h"
#include "InitializeIntegratorComponent.h"
#include "NumericalIntegratorComponent.h"
#include "DtIntegratorComponent.h"
#include "ReductionIntegratorComponent.h"
#include "CopyIntegratorComponent.h"
#include "RebalanceIntegratorComponent.h"

#include "StandardComponentPatchStrategy.h"
#include "MemoryIntegratorComponent.h"
#include "LinearSolverManager.h"


#include "BaseLinearSolver.h"
#include "JVector.h"
#include "NSAbstractFunctions.h"
#include "J_BaseNonlinearSolver.h"
#include "J_KINSOLSolver.h"
#include "J_SNESSolver.h"

#include "NonlinearPoisson.h"
#include "Global_variable.h"
#include "Variables_Functions.h"

using namespace JAUMIN;


/**
 * @brief 该类从网格层时间积分算法策略类 algs::TimeIntegratorLevelStrategy 派生,
 * 实现有限元方法的流程控制.
 */   

class NonlinearPoissonLevel : 
  public algs::TimeIntegratorLevelStrategy<NDIM>
{
public:
  /**
   * @brief 构造函数.
   * @param object_name   输入参数, 字符串, 表示对象名称.
   * @param strategy      输入参数, 指针, 网格片策略类派生类对象.
   */     
  NonlinearPoissonLevel(
      const string& object_name,
      tbox::Pointer<NonlinearPoisson> strategy,
      tbox::Pointer<tbox::Database> input_db);

  /**
   * @brief 析构函数.
   */
  virtual ~NonlinearPoissonLevel();

  ///@name 重载基类algs::TimeIntegratorLevelStrategy<NDIM>的函数
  //@{
  /**
   * @brief 初始化该积分算法: 创建所有计算需要的积分构件.
   * @param manager 输入参数, 指针, 指向积分构件管理器.
   */
  void initializeLevelIntegrator(
        tbox::Pointer<algs::IntegratorComponentManager<NDIM> > manager);

  /**
   * @brief 初始化指定网格层的数据片.
   *
   * 具体地，该函数完成以下操作：
   * - 若输入参数 initial_time 为真，则根据初始条件，
   *   为当前网格层上的所有数据片 <uval,current> 赋初值;
   * - 若输入参数 initial_time 为假（此时刚完成负载调整），
   *   则将数据片 <uval,current> 的值从旧网格层复制到新网格层.
   *
   * @param level 输入参数, 指针, 指向待初始化网格层.
   * @param init_data_time 输入参数, 双精度浮点型, 表示初始化的时刻.
   * @param initial_time 输入参数, 逻辑型, 真值表示当前时刻为计算的初始时刻.
   *
   * @note
   * 该函数调用了1个初值构件对象，该对象又进一步自动调用函数
   * FEM::initializePatchData(), 逐个网格片地完成初始时刻的数据初始化.
   */
  void initializeLevelData(const tbox::Pointer< hier::BasePatchLevel<NDIM> > level, 
                           const double init_data_time,
                           const bool   initial_time = true);


  /**
       * @brief 返回指定网格层的时间步长.
       *
       * @param level 输入参数, 指针, 指向网格层.
       * @param dt_time 输入参数, 双精度浮点型, 表示计算时间步长的当前时刻.
       * @param initial_time 输入参数, 逻辑型, 真值表示当前时刻为计算的初始时刻.
       * @param flag_last_dt 输入参数, 整型, 表示上个时间步积分返回的状态.
       * @param last_dt 输入参数, 双精度浮点型, 表示上个时间步长.
       * @return 双精度浮点型, 表示网格层的时间步长.
       *
       * @note
       * 该函数调用了1个时间步长构件对象，该对象又进一步自动调用函数
       * FEM::getPatchDt(), 逐个网格片地计算稳定时间步长.
       */
    double getLevelDt(const tbox::Pointer<hier::BasePatchLevel<NDIM> > level,
                        const double dt_time, const bool initial_time,
                        const int flag_last_dt, const double last_dt);


  /**
   * @brief 网格层向前积分一个时间步. 
   *
   * @param level 输入参数, 指针, 指向待积分的网格层.
   * @param current_time 输入参数, 双精度浮点型, 表示时间步的起始时刻.
   * @param predict_dt 输入参数, 双精度浮点型, 表示为该时间步预测的时间步长.
   * @param max_dt 输入参数, 双精度浮点型, 表示时间步允许的最大时间步长.
   * @param min_dt 输入参数, 双精度浮点型, 表示时间步允许的最小时间步长.
   * @param first_step 输入参数, 逻辑型, 真值当前为重构后或时间步序列的第1步.
   * @param step_number, 输入参数, 整型, 表示积分步数.
   * @param actual_dt 输出参数, 双精度浮点型, 表示时间步实际采用的时间步长.
   * @return 整型, 表示该时间步积分的状态. 
   *
   */ 
  int advanceLevel(const tbox::Pointer< hier::BasePatchLevel<NDIM> > level, 
                   const double current_time, 
                   const double predict_dt,
                   const double max_dt,
                   const double min_dt,
                   const bool   first_step, 
                   const int    step_number,
                   double &     actual_dt);

  //print matrix equation
  void MatrixEquationOutput(tbox::Pointer< hier::PatchLevel<NDIM> > patch_level,
                              int matrix_id,
                              int fx_id,
                              string matrix_name,
                              string fx_name){
       //print left matrix
       #if 1
       tbox::Pointer<JPSOL::JMatrix<NDIM,double> > d_matrix;
       d_matrix = new JPSOL::JMatrix<NDIM,double>(patch_level,matrix_id);
       d_matrix->printMatrix(matrix_name,false,false);
       //print right vector
       tbox::Pointer<JPSOL::JVector<NDIM,double> > d_fx;
       d_fx = new JPSOL::JVector<NDIM,double>(patch_level,fx_id);
       d_fx->printVector(fx_name,false,false);
       #endif
  };

  //print solution
  void SolutionOutput(tbox::Pointer< hier::PatchLevel<NDIM> > patch_level,
                              int fx_id,
                              string fx_name){
       //print vector
       tbox::Pointer<JPSOL::JVector<NDIM,double> > d_fx;
       d_fx = new JPSOL::JVector<NDIM,double>(patch_level,fx_id);
       d_fx->printVector(fx_name,false,false);
  };

  void IV_Output(double* I,int num_bc,tbox::Array<electrobc> EBoundassignment,
                 const double time, const double dt){
      if(tbox::MPI::getRank()==0){
          tbox::Pointer<Variables_Functions> Var_func = new Variables_Functions();
          /** setup the source target **/
          tbox::Pointer<bound_condition> bound_source = new bound_condition();
          for(int i_bc=0; i_bc<num_bc; i_bc++){
              double app_volt = bound_source->source_assign(EBoundassignment[i_bc].sc_type,
                                                            time, dt,
                                                            EBoundassignment[i_bc].sc_coe);
              Var_func->DataOutput("current",i_bc,I[i_bc],time,dt);
              Var_func->DataOutput("voltage",i_bc,app_volt,time,dt);
          }
      }
  }

  void error_Output(string outputname, int iter, double error,
                    const double time, const double dt){
      if(tbox::MPI::getRank()==0){
          tbox::Pointer<Variables_Functions> Var_func = new Variables_Functions();
          Var_func->ErrorOutput(outputname,iter,error,time,dt);
      }
  }

  //@}
  tbox::Pointer<NonlinearPoisson> get_patch_strategy(){ return d_patch_strategy; }
  // poisson equation
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > get_Pphybc_intc(){ return d_Pphybc_intc; }
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > get_Pfx_intc(){ return d_Pfx_intc; }
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > get_Pmat_intc(){ return d_Pmat_intc; }
  tbox::Pointer<JPSOL::JVector<NDIM,double> > get_Pstetch_vector(){ return poisson_stetch; }
  tbox::Pointer<JPSOL::JMatrix<NDIM,double> > get_Pmatrix(){ return poisson_mat; }
  // post simulation
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > get_Ppost_intc(){ return d_Ppost_intc; }

 private:

  /*!@brief 对象名称. */
  string d_object_name; 
 // string d_nonlinear_solver_name;

  //bool d_uses_preconditioner;
  //bool d_uses_explicit_jacobian;

  /*!@brief 网格片策略类 */
  tbox::Pointer<NonlinearPoisson> d_patch_strategy;

  /* memory component*/
  tbox::Pointer<algs::MemoryIntegratorComponent<NDIM> > d_alloc_n_c_vec;
  tbox::Pointer<algs::MemoryIntegratorComponent<NDIM> > d_alloc_vec;
  tbox::Pointer<algs::MemoryIntegratorComponent<NDIM> > d_alloc_CCmat;
  tbox::Pointer<algs::MemoryIntegratorComponent<NDIM> > d_alloc_Pmat;
  tbox::Pointer<algs::MemoryIntegratorComponent<NDIM> > d_alloc_Newton_Pmat;
  tbox::Pointer<algs::MemoryIntegratorComponent<NDIM> > d_alloc_Cmat;
  tbox::Pointer<algs::MemoryIntegratorComponent<NDIM> > d_alloc_CCrhs;
  tbox::Pointer<algs::MemoryIntegratorComponent<NDIM> > d_alloc_Prhs;
  tbox::Pointer<algs::MemoryIntegratorComponent<NDIM> > d_alloc_Newton_Prhs;
  tbox::Pointer<algs::MemoryIntegratorComponent<NDIM> > d_alloc_Crhs;
  /*!@brief 初始构件 */
  tbox::Pointer<algs::InitializeIntegratorComponent<NDIM> > d_init_intc;
  tbox::Pointer<algs::InitializeIntegratorComponent<NDIM> > d_dop_intc;
  tbox::Pointer<algs::InitializeIntegratorComponent<NDIM> > d_tid_intc;
  tbox::Pointer<algs::InitializeIntegratorComponent<NDIM> > d_init_dof_intc;
  /**8.11**/
  tbox::Pointer<algs::InitializeIntegratorComponent<NDIM> > d_init_sol_intc;

  /*initialize solution vector*/
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_solvec_init_intc;
  /** check the node order in element**/
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_checkorder_intc;

  /*!@brief 数值构件，设定边值 */
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_CCapp_bc;
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_Papp_bc;
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_Newton_Papp_bc;
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_Napp_bc;
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_Happ_bc;
  /*!@brief 数值构件，组装矩阵 */
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_CCmat_intc;
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_Pmat_intc;
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_Newton_Pmat_intc;
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_adtest_intc;
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_Nmat_intc;
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_Hmat_intc;
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_CVmat_intc;
  /*!@brief 数值构件，计算非线性残差F(x) */
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_CCfx_intc;
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_Pfx_intc;
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_Newton_Pfx_intc;
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_Nfx_intc;
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_Hfx_intc;
  /*!@brief 数值构件，设置边界条件 */
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_CCphybc_intc;
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_Pphybc_intc;
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_Newton_Pphybc_intc;
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_Nphybc_intc;
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_Hphybc_intc;
  /* @brief joule heat, e-field, current density*/
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_Jheat_intc;
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_Efield_intc;
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_Cdens_intc;

  /* @brief semiconductor variables*/
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_Semivar_intc;
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_SemivarInit_intc;
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_Flevel_intc;
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_NPconstraint_intc;
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_Fintrin_intc;
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_equilibrium_intc;
  /* @brief solution scaled*/
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_sol_scal_intc;
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_sol_antiscal_intc;
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_sol_anticar_intc;
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_sol_carnew_intc;
  /* @ post solution process*/
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_sol_process_intc;
  /* solution transport */
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_Psolution_intc;
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_GummelPsolution_intc;
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_Semisolution_intc;
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_GummelSemisolution_intc;
  /*!@brief 数值构件，计算矩阵乘向量 */
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_mv_intc;
  /*!@brief 数值构件，后处理计算 */
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_Ppost_intc;
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_Npost_intc;
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_Hpost_intc;

  /*!@brief 电场数值构件，后处理计算 */
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_E_intc;
  /*!@brief 电流密度数值构件，后处理计算 */
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_J_intc;
  /*!@brief 获取电流电压数值构件，后处理计算 */
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_GetIV_intc;

  /*!@brief 规约构件 */
  tbox::Pointer<algs::ReductionIntegratorComponent<NDIM> > d_Preduction_intc;
  tbox::Pointer<algs::ReductionIntegratorComponent<NDIM> > d_Dopreduction_intc;
  tbox::Pointer<algs::ReductionIntegratorComponent<NDIM> > d_Smeireduction_intc;
  tbox::Pointer<algs::ReductionIntegratorComponent<NDIM> > d_GummelSmeireduction_intc;
  tbox::Pointer<algs::ReductionIntegratorComponent<NDIM> > d_getI_intc;
  tbox::Pointer<algs::ReductionIntegratorComponent<NDIM> > d_Poireduction_intc;
  tbox::Pointer<algs::ReductionIntegratorComponent<NDIM> > d_saveSol_intc;

  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_poison_comm_intc;
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_cceq_comm_intc;

  /*!@brief 网格重分布构件  */
  tbox::Pointer<algs::RebalanceIntegratorComponent<NDIM> > d_rebalance_intc;
  tbox::Pointer<algs::DtIntegratorComponent<NDIM> > d_dt_intc;

  /* linear solver */
  tbox::Pointer<solv::LinearSolverManager<NDIM> > d_linear_solver_manager;
  tbox::Pointer<solv::BaseLinearSolver<NDIM> > d_linear_solver_LU;
  tbox::Pointer<solv::BaseLinearSolver<NDIM> > d_linear_solver_PCG;
  tbox::Pointer<tbox::Database> d_linear_solver_db;

  /// 非线性系统的解向量
  // poisson equation
  tbox::Pointer<JPSOL::JMatrix<NDIM,double> > poisson_mat;
  tbox::Pointer<JPSOL::JVector<NDIM,double> > poisson_sol;
  tbox::Pointer<JPSOL::JVector<NDIM,double> > poisson_fval;
  tbox::Pointer<JPSOL::JVector<NDIM,double> > poisson_stetch;
  tbox::Pointer<JPSOL::JVector<NDIM,double> > poisson_sol_pre;


  /// 非线性解法器基类
//  NSSlover * d_nonlinear_solver;
//  JPSOL::NSAbstractFunctions<NDIM, double> * d_nonlinear_function;

  /// error for nonlinear iterations
  double error_cp, error_non, error_semi, error_poi, tol_cp, tol_non, tol_semi, tol_poi;
  int iter_cp_max, iter_non_max, iter_semi_max, iter_poi_max;
  /** boundary assignment **/
  tbox::Array<electrobc> EBoundassignment;
  int num_bc;
  bool enable_selfDop;
  bool save_bas;
  bool enable_tid;
  bool get_InitialSol;
  bool save_PoissonSolcur;
  bool get_PoissonSolcur;
  string debug_flag;
};

#endif
