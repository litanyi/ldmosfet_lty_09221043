//
// 文件名:     NonlinearPoisson.h
// 软件包:     JAUMIN 
// 版权　:     北京应用物理与计算数学研究所 
// 版本号:     $Revision: 0 $
// 修改　:     $Date: Wed Apr 15 08:47:03 2015 $
// 描述　:     
// 类别　:     %Internal File% ( Don't delete this line )
//


#ifndef included_NonlinearPoisson
#define included_NonlinearPoisson

#include "StandardComponentPatchStrategy.h"
#include "JaVisDataWriter.h"
#include "DOFInfo.h"
#include "RestartManager.h"

#include "Bound_condition.h"
#include "material_database.h"
#include "BaseMaterial.h"
#include "Global_variable.h"
#include "MaterialManager.h"
#include "physical_unit.h"

#include "BaseElement.h"

using namespace JAUMIN;
using namespace PhysicalUnit;

using PhysicalUnit::cm;
using PhysicalUnit::nm;
using PhysicalUnit::m;
using PhysicalUnit::s;
using PhysicalUnit::V;
using PhysicalUnit::C;
using PhysicalUnit::K;
using PhysicalUnit::g;
using PhysicalUnit::A;
using PhysicalUnit::eV;
using PhysicalUnit::kb;
using PhysicalUnit::e;
using PhysicalUnit::eps0;
using PhysicalUnit::mu0;
using PhysicalUnit::Pa;
using PhysicalUnit::max_dop_ptr;

class NonlinearPoisson
: public algs::StandardComponentPatchStrategy<NDIM>, tbox::Serializable
{
public:

  /*! @brief 构造函数.
   * @param object_name          输入参数, 字符串, 表示对象名称.
   */
  NonlinearPoisson(const string& object_name,
                   bool is_from_restart,
                   tbox::Pointer<tbox::Database> input_db);

  /*!
   * @brief 析构函数.
   */
  virtual ~NonlinearPoisson();

  /// @name 重载基类 algs::StandardComponentPatchStrategy<NDIM> 的函数:
  // @{
      
  /*!
   * @brief 初始化指定的积分构件.
   *
   * 注册待填充的数据片或待调度内存空间的数据片到积分构件.
   *
   * @param component 输入参数, 指针, 指向待初始化的积分构件对象.
   */
  void initializeComponent(algs::IntegratorComponent<NDIM> * component) const;

  /**
   * @brief 初始化数据片（支持有限元初值构件）.
   *
   * @param patch          输入参数, 网格片类, 表示网格片.
   * @param time           输入参数, 双精度浮点型, 表示初始化的时刻.
   * @param initial_time   输入参数, 逻辑型, 真值表示当前时刻为计算的初始时刻.
   * @param component_name 输入参数, 字符串, 当前调用该函数的初值构件之名称.
   */
  void initializePatchData(hier::Patch<NDIM> & patch,
                           const double  time, 
                           const bool    initial_time, 
                           const string& component_name); 

  /**
   * @brief get the cell - entity relationship
   **/
  void getCellEntityRelation(hier::Patch<NDIM> & patch,
                           const string& component_name);

  /**
   * @brief get the cell - material relationship
   **/
  void getCellMatRelation(hier::Patch<NDIM> & patch,
                           const string& component_name);

  /**
   * @brief get the cell - entity relationship: (num_entity, num_nodes)
   **/
  void getNodeEntityRelation(hier::Patch<NDIM> &patch,
  		                   const string& component_name);

  /**
   * @brief get the cell - entity relationship: (num_mat, num_nodes)
   **/
  void getNodeMatRelation(hier::Patch<NDIM> &patch,
                           const string& component_name);

  /**
   * @brief initialize quasi fermi level
   **/
  void initialQuasiFermilevel(hier::Patch<NDIM> & patch,
                        const double  time,
                        const double    initial_time,
                        const string& component_name);

  /**
   * @brief calculate Quasi Fermi level
   */
  void calculateQuasiFermilevel(hier::Patch<NDIM> &patch,
                        const double time,
                        const double dt,
                        const string & component_name);

  /**
   * @brief apply constraint to hole and electron
   */
  void applyConstraint2NP(hier::Patch<NDIM> &patch,
                          const double time,
                          const double dt,
                          const string & component_name);

  /**
   * @brief initialize the solution vectors
   **/
  void solVecInitialization(hier::Patch<NDIM> & patch,
                          const double  time,
                          const double  dt,
                          const string& component_name);

  /**
   * @brief check the node order in element
   **/
  void checkOrder(hier::Patch<NDIM> & patch,
                  const double  time,
                  const double  dt,
                  const string& component_name);

  /**
   * @brief get the neutral status for semiconductor domains
   **/
  void neutral_status(hier::Patch<NDIM> & patch,
                      const double  time,
                      const double dt,
                      const string& component_name);

  /**
   * @brief get the neutral status for semiconductor domains
   **/
  void potential_bound(hier::Patch<NDIM> & patch,
                       const double  time,
                       const double dt,
                       const string& component_name);


  /**
   * @brief reduceOnPatch: get the system level value
   * @param vector
   * @param len
   * @param patch
   * @param time
   * @param dt
   * @param component_name
   */
  void reduceOnPatch(double* vector, int len,
                          hier::Patch<NDIM>& patch,
                          const double time, const double dt,
                          const string& component_name);

  /**
   * @brief PotentialMaxError: get the maximum change in potential solution
   * @param vector
   * @param patch
   * @param time
   * @param dt
   * @param component_name
   */
  void PotentialMaxError(double* vector,
                          hier::Patch<NDIM>& patch,
                          const double time, const double dt,
                          const string& component_name);

  /**
   * @brief NonlinearPoisson::PoissonMaxError: to find the maximum change in potential
   * @param vector
   * @param patch
   * @param time
   * @param dt
   * @param component_name
   */
  void PoissonMaxError(double* vector,
                                      hier::Patch<NDIM>& patch,
                                      const double time, const double dt,
                                      const string& component_name);

  /**
   * @brief SemiMaxError: get the maximum change in quasi fermi level
   * @param vector
   * @param len
   * @param patch
   * @param time
   * @param dt
   * @param component_name
   */
  void SemiMaxError(double* vector, int len,
                          hier::Patch<NDIM>& patch,
                          const double time, const double dt,
                          const string& component_name);

  /**
   * @brief NonlinearPoisson::MaxDop: to find the maximum doping
   * @param vector
   * @param len
   * @param patch
   * @param time
   * @param dt
   * @param component_name
   */
  void MaxDop(double* vector, int len,
                                      hier::Patch<NDIM>& patch,
                                      const double time, const double dt,
                                      const string& component_name);

  /**
   * @brief NonlinearPoisson::saveSol
   * @param vector
   * @param len
   * @param patch
   * @param time
   * @param dt
   * @param component_name
   */
  void saveSol(double* vector, int len,
                                      hier::Patch<NDIM>& patch,
                                      const double time, const double dt,
                                      const string& component_name);

  /**
   * @brief NonlinearPoisson::GummelSemiMaxError: to find the maximum change in quasi fermi level
   * @param vector
   * @param len
   * @param patch
   * @param time
   * @param dt
   * @param component_name
   */
  void GummelSemiMaxError(double* vector, int len,
                                      hier::Patch<NDIM>& patch,
                                      const double time, const double dt,
                                      const string& component_name);

  /**
   * @brief NonlinearPoisson::getI: to find the current
   * @param vector
   * @param len
   * @param patch
   * @param time
   * @param dt
   * @param component_name
   */
  void getI(double* vector, int len,
            hier::Patch<NDIM>& patch,
            const double time, const double dt,
            const string& component_name);

  /**
   * @brief transportPoissonSol: from current to stetch
   * @param patch
   * @param time
   * @param dt
   * @param component_name
   */
  void transportPoissonSol(hier::Patch<NDIM>& patch,
                           const double time,
  						   const double dt,
                           const string& component_name);

  /**
   * @brief transportSemiSol: from current to stetch
   * @param patch
   * @param time
   * @param dt
   * @param component_name
   */
  void transportSemiSol(hier::Patch<NDIM>& patch,
                           const double time,
  						   const double dt,
                           const string& component_name);

  /**
   * @brief NonlinearPoisson::transportGummelSemiSol: fermi level-from current to stetch
   * @param patch
   * @param time
   * @param dt
   * @param component_name
   */
  void transportGummelSemiSol(hier::Patch<NDIM>& patch,
                                          const double time,
                                          const double dt,
                                          const string& component_name);

  /**
   * @brief transportSystemSol: from current to stetch (temperature)
   * @param patch
   * @param time
   * @param dt
   * @param component_name
   */
  void transportSystemSol(hier::Patch<NDIM>& patch,
                           const double time,
  						   const double dt,
                           const string& component_name);

  /**
   * @brief stepSolPostProcess: log(carrier_density)
   * @param patch
   * @param time
   * @param dt
   * @param component_name
   */
  void stepSolPostProcess(hier::Patch<NDIM>& patch,
                           const double time,
                           const double dt,
                           const string& component_name);

/** *************************************************************************
 * build up the left matrix and right vector
 ************************************************************************* **/
  /**
   * @brief buildCurConserveMatOnPatch: current continuity (metal region)
   * @param patch
   * @param time
   * @param dt
   * @param component_name
   */
  void buildCurConserveMatOnPatch(hier::Patch<NDIM> & patch,
                             const double  time,
                             const double  dt,
                             const string& component_name);

  /**
   * @brief buildCurConserveFXOnPatch: current continuity (metal region)
   * @param patch
   * @param time
   * @param dt
   * @param component_name
   */
  void buildCurConserveFXOnPatch(hier::Patch<NDIM> & patch,
                             const double  time,
                             const double  dt,
                             const string& component_name);

  /**
   * @brief setCurConservePhysicalBC: matrix
   * @param patch
   * @param time
   * @param dt
   * @param component_name
   */
  void setCurConservePhysicalBC(hier::Patch<NDIM> & patch,
                             const double  time,
                             const double  dt,
                             const string& component_name);

  /**
   * @brief applyCurConserveBCOnPatch: vector
   * @param patch
   * @param time
   * @param dt
   * @param component_name
   */
  void applyCurConserveBCOnPatch(hier::Patch<NDIM> & patch,
                             const double  time,
                             const double  dt,
                             const string& component_name);
  /**
   * @brief buildPoissonMatOnPatch: for initialization
   * @param patch
   * @param time
   * @param dt
   * @param component_name
   */
  void buildPoissonMatOnPatch(hier::Patch<NDIM> & patch,
                            const double  time,
                            const double  dt,
                            const string& component_name);

  /**
   * @brief buildPoissonFXOnPatch: for initialization
   * @param patch
   * @param time
   * @param dt
   * @param component_name
   */
  void buildPoissonFXOnPatch(hier::Patch<NDIM> & patch,
                            const double  time,
                            const double  dt,
                            const string& component_name);

  /**
   * @brief setPoissonPhysicalBC: matrix
   * @param patch
   * @param time
   * @param dt
   * @param component_name
   */
  void setPoissonPhysicalBC(hier::Patch<NDIM> & patch,
                            const double  time,
                            const double  dt,
                            const string& component_name);

  /**
   * @brief applyPoissonBCOnPatch: vector
   * @param patch
   * @param time
   * @param dt
   * @param component_name
   */
  void applyPoissonBCOnPatch(hier::Patch<NDIM> & patch,
                      const double  time,
                      const double  dt,
                      const string& component_name);

  void adTest(hier::Patch<NDIM> & patch,
              const double  time,
              const double  dt,
              const string& component_name);

  /**
   * @brief NonlinearPoisson::buildNewtonVogammaMatOnPatch: Newton-Raphson
   * @param patch
   * @param time
   * @param dt
   * @param component_name
   */
  void buildNewtonVogammaMatOnPatch(hier::Patch<NDIM> & patch,
                                            const double  time,
                                            const double  dt,
                                            const string& component_name);

  /**
   * @brief buildNewtonPMatOnPatch: Raphson-Newton
   * @param patch
   * @param time
   * @param dt
   * @param component_name
   */
  void buildNewtonPMatOnPatch(hier::Patch<NDIM> & patch,
                              const double  time,
                              const double  dt,
                              const string& component_name);

  /**
   * @brief buildNewtonPFXOnPatch: Raphson-Newton
   * @param patch
   * @param time
   * @param dt
   * @param component_name
   */
  void buildNewtonPFXOnPatch(hier::Patch<NDIM> & patch,
                             const double  time,
                             const double  dt,
                             const string& component_name);

  /**
   * @brief setNewtonPPhysicalBC: matrix
   * @param patch
   * @param time
   * @param dt
   * @param component_name
   */
  void setNewtonPPhysicalBC(hier::Patch<NDIM> & patch,
                            const double  time,
                            const double  dt,
                            const string& component_name);

  /**
   * @brief applyNewtonPBCOnPatch: vector
   * @param patch
   * @param time
   * @param dt
   * @param component_name
   */
  void applyNewtonPBCOnPatch(hier::Patch<NDIM> & patch,
                             const double  time,
                             const double  dt,
                             const string& component_name);

  /**
   * @brief transferPoissonSoln: transfer solution to plot data
   * @param patch
   * @param time
   * @param dt
   * @param component_name
   */
  void transferPoissonSoln(hier::Patch<NDIM> & patch,
                    const double  time,
                    const double  dt,
                    const string& component_name);

  /**
   * @brief NonlinearPoisson::transportGummelPoissonSol: potential-from current to stetch
   * @param patch
   * @param time
   * @param dt
   * @param component_name
   */
  void transportGummelPoissonSol(hier::Patch<NDIM>& patch,
                                          const double time,
                                          const double dt,
                                          const string& component_name);

  /**
   * @brief buildCarMatOnPatch: CV-FEM-SG for carrier continuity equation
   * @param patch
   * @param time
   * @param dt
   * @param component_name
   */
  void buildCarMatOnPatch(hier::Patch<NDIM> & patch,
                            const double  time,
                            const double  dt,
                            const string& component_name);

  /**
   * @brief buildCarnFXOnPatch: rhs for electron continuity equation
   * @param patch
   * @param time
   * @param dt
   * @param component_name
   */
  void buildCarnFXOnPatch(hier::Patch<NDIM> & patch,
                            const double  time,
                            const double  dt,
                            const string& component_name);

  /**
   * @brief setCarnPhysicalBC: matrix
   * @param patch
   * @param time
   * @param dt
   * @param component_name
   */
  void setCarnPhysicalBC(hier::Patch<NDIM> & patch,
                            const double  time,
                            const double  dt,
                            const string& component_name);

  /**
   * @brief applyCarnBCOnPatch: vector
   * @param patch
   * @param time
   * @param dt
   * @param component_name
   */
  void applyCarnBCOnPatch(hier::Patch<NDIM> & patch,
                      const double  time,
                      const double  dt,
                      const string& component_name);

  /**
   * @brief transferCarnSoln: transfer solution to plot data
   * @param patch
   * @param time
   * @param dt
   * @param component_name
   */
  void transferCarnSoln(hier::Patch<NDIM> & patch,
                    const double  time,
                    const double  dt,
                    const string& component_name);

  /**
   * @brief buildCarhFXOnPatch: rhs for hole continuity equation
   * @param patch
   * @param time
   * @param dt
   * @param component_name
   */
  void buildCarhFXOnPatch(hier::Patch<NDIM> & patch,
                            const double  time,
                            const double  dt,
                            const string& component_name);

  /**
   * @brief setCarhPhysicalBC: matrix
   * @param patch
   * @param time
   * @param dt
   * @param component_name
   */
  void setCarhPhysicalBC(hier::Patch<NDIM> & patch,
                            const double  time,
                            const double  dt,
                            const string& component_name);

  /**
   * @brief applyCarhBCOnPatch: vector
   * @param patch
   * @param time
   * @param dt
   * @param component_name
   */
  void applyCarhBCOnPatch(hier::Patch<NDIM> & patch,
                      const double  time,
                      const double  dt,
                      const string& component_name);

  /**
   * @brief transferCarhSoln: transfer solution to plot data
   * @param patch
   * @param time
   * @param dt
   * @param component_name
   */
  void transferCarhSoln(hier::Patch<NDIM> & patch,
                    const double  time,
                    const double  dt,
                    const string& component_name);


  /**
   * @brief solutionAntiScal: antiscaling for solution
   * @param patch
   * @param time
   * @param dt
   * @param component_name
   */
  void solutionAntiScal(hier::Patch<NDIM> & patch,
                              const double  time,
                              const double  dt,
                              const string& component_name);

  /**
   * @brief solutionScal: scaling solution
   * @param patch
   * @param time
   * @param dt
   * @param component_name
   */
  void solutionScal(hier::Patch<NDIM> & patch,
                          const double  time,
                          const double  dt,
                          const string& component_name);

  /**
   * @brief carrierAntiScal: anti-scaling carrier solution
   * @param patch
   * @param time
   * @param dt
   * @param component_name
   */
  void carrierAntiScal(hier::Patch<NDIM> & patch,
                          const double  time,
                          const double  dt,
                          const string& component_name);

  /*************************************************************************
   * calculate the electric field
   ************************************************************************/
  void calculateEx(hier::Patch<NDIM>& patch,
                       const double time,
                       const double dt,
                       const string& component_name);

  /*************************************************************************
   * calculate the electric current density
   ************************************************************************/
  void calculateElemJ(hier::Patch<NDIM> &patch,
                       const double time,
                       const double dt,
                       const string & component_name);

  /*************************************************************************
   * get I-V data
   ************************************************************************/
  void getIV(hier::Patch<NDIM> &patch,
               const double time,
               const double dt,
               const string & component_name);

  /** 
   * @brief 注册模型变量, 完成用户变量的注册.
   * 
   */
  void registerModelVariable();
   
  /*!
   * @brief 支撑指定名称的有限元数值构件, 在单个网格片上完成后数值计算.
   *
   * @param patch          输入参数, 网格片类, 表示网格片.
   * @param time           输入参数, 双精度浮点型, 表示当前时刻.
   * @param dt             输入参数, 双精度浮点型, 表示时间步长.
   * @param initial_time   输入参数, BOOL型, 是否为初始时刻.
   * @param component_name 输入参数, 字符串, 表示数值构件的名称.
   *
   */
  void computeOnPatch(hier::Patch<NDIM>& patch,
                      const double  time,
                      const double  dt,
                      const bool    initial_time,
                      const string& component_name);

  /** 
   * @brief 注册可视化数据.
   */
  void registerPlotData(tbox::Pointer<appu::JaVisDataWriter<NDIM> > javis_data_writer);

  /**
   * @brief getPatchDt
   * @param patch
   * @param time
   * @param initial_time
   * @param flag_last_dt
   * @param last_dt
   * @param component_name
   * @return
   */
  double getPatchDt(hier::Patch<NDIM>& patch, const double time,
                    const bool initial_time, const int flag_last_dt,
                    const double last_dt, const string& component_name);

  /** 
   * @brief 从输入数据库中读取参数, 并设置到计算流程中.
   * 
   * @param input_db 输入参数, 指针, 指向输入数据库.
   */
   void setParameter(tbox::Pointer<tbox::Database> input_db);

   /** 
    * @brief 从输入数据库中读如参数。
    * 
    * @param db 输入参数，指针，指向输入数据库。
    */
   void getFromInput(tbox::Pointer<tbox::Database> db);

   /** 
    * @brief 从重启动数据库中读如参数。
    * 
    * @param db 输入参数，指针，指向输入数据库。
    */
   void getFromRestart(tbox::Pointer<tbox::Database> db);

   /** 
    * @brief 将数据写入到重启动数据库。
    * 
    * @param db 输入参数，指针，指向重启动数据库。
    */
   void putToDatabase(tbox::Pointer<tbox::Database> db);
   /** 
    * @brief  获取矩阵id 
    * 
    * @return 整型，矩阵id
    */
   void setPoissonFXDataID(int fval_id)
   {
     if(fval_id >= 0) d_poisson_fval_id = fval_id;
   }

   /******************************************************************************
    get the id of matrix and solutions return int variable
    */
   // current conservation equation
   int getCurConserveMatID()
   {return d_curconserve_mat_id;}
   int getCurConserveFvalID()
   {return d_curconserve_fval_id;}
   int getCurConserveSolID()
   {return d_curconserve_sol_id;}

   // Poisson equation
   int getPoissonMatrixID()
   {return d_poisson_mat_id;}
   int getPoissonFxID()
   {return d_poisson_fval_id;}
   int getPoissonSolID()
   {return d_poisson_sol_id;}
   int getPoissonSolcurID()
   {return d_poisson_sol_cur_id;}
   int getPoissonStetchID()
   {return d_poisson_stetch_id;}
   int getSolID()
   {return d_sol_id;}

   int getPoissonNewtonMatrixID()
   {return d_poisson_newton_mat_id;}
   int getPoissonNewtonSolID()
   {return d_poisson_newton_sol_id;}   
   int getPoissonNewtonFxID()
   {return d_poisson_newton_fval_id;}


   int getVogammaNewtonMatrixID()
   {return d_Vogamma_newton_mat_id;}
   int getVogammaNewtonSolID()
   {return d_Vogamma_newton_sol_id;}
   int getVogammaNewtonFxID()
   {return d_Vogamma_newton_fval_id;}

   int getVodeltaNewtonMatrixID()
   {return d_Vodelta_newton_mat_id;}
   int getVodeltaNewtonSolID()
   {return d_Vodelta_newton_sol_id;}
   int getVodeltaNewtonFxID()
   {return d_Vodelta_newton_fval_id;}

   // Carrier transport equation: electron
   int getCarnMatrixID()
   {return d_carn_mat_id;}
   int getCarnSolID()
   {return d_carn_sol_cur_id;}
   int getCarnStetchID()
   {return d_carn_stetch_id;}
   int getCarnFxID()
   {return d_carn_fval_id;}
   int getFerminSolID()
   {return d_fermin_sol_id;}

   // Carrier transport equation: hole
   int getCarhMatrixID()
   {return d_carh_mat_id;}
   int getCarhSolID()
   {return d_carh_sol_cur_id;}
   int getCarhStetchID()
   {return d_carh_stetch_id;}
   int getCarhFxID()
   {return d_carh_fval_id;}
   int getFermihSolID()
   {return d_fermih_sol_id;}

/******************************************************************/

   /** 
    * @brief  获取右端项id 
    * 
    * @return 整型，右端项id
    */
   void setPoissonMVDataID(int soln_id, int prod_id, int vec_id)
   {
     if(soln_id >= 0) d_poisson_stetch_id = soln_id;
   }

   tbox::Pointer<tbox::Vector<double> > getRhs(int n_vertex,bool save_bas,int cell,
                                                  tbox::Pointer<BaseElement<NDIM> > ele,
                                                  tbox::Pointer<pdat::CellData<NDIM, double> > bas_rhs,
                                                  tbox::Array<double> value,
                                                  tbox::Array<hier::DoubleVector<NDIM> > vertex,
                                                  const double  time,
                                                  const double  dt){
       tbox::Pointer<tbox::Vector<double> > cv_vol = new tbox::Vector<double>();
       cv_vol->resize(n_vertex);
       /** building element vector **/
       if(save_bas){
           for(int i_v = 0; i_v < cell_point; ++i_v){
               (*cv_vol)[i_v] = (*bas_rhs)(i_v,cell)*value[i_v];
           }
       }else{
           ele->buildCVFEMRhs(vertex, value, dt, time, cv_vol);
       }
       return cv_vol;
   }

   tbox::Array<tbox::Array<double> > getInteg(int n_vertex,bool save_bas,int cell,
                                                  tbox::Pointer<BaseElement<NDIM> > ele,
                                                  tbox::Pointer<pdat::CellData<NDIM, double> > bas_integ,
                                                  tbox::Array<hier::DoubleVector<NDIM> > vertex,
                                                  const double  time,
                                                  const double  dt){
       tbox::Array<tbox::Array<double> > r_integ;
       if(save_bas){
           r_integ.resizeArray(cell_edge);
           for (int i_row = 0; i_row < cell_edge; ++i_row){
               r_integ[i_row].resizeArray(cell_point);
           }
           for(int i_v = 0; i_v < cell_point; ++i_v){
               for(int i_e = 0; i_e < cell_edge; ++i_e){
                   r_integ[i_e][i_v] = ((*bas_integ)((i_e*cell_point+i_v),cell));
               }
           }
       }else{
           r_integ = ele->buildCVFEMInteg(vertex, dt, time);
       }
       return r_integ;
   }

   void SolRead(hier::Patch<NDIM> & patch,const std::string filename){
       tbox::Pointer<pdat::NodeData<NDIM, double> >
               sol_data = patch.getPatchData(d_sol_id);

       std::fstream input(filename);

       if(!input) {
           TBOX_ERROR("fail to open file: "<< filename << endl);
       }

       string line;

       getline(input,line);

       line.clear();

       while( getline(input, line) )
       {
           stringstream buf(line);
           int node_id;
           double volt_node;
           double n_node;
           double p_node;

           //read and convert data
           buf >> node_id;
           buf >> volt_node;
           buf >> n_node;
           buf >> p_node;

           (*sol_data)(0,node_id) = volt_node;
           (*sol_data)(1,node_id) = n_node;
           (*sol_data)(2,node_id) = p_node;
           //cout<<"sol_data:"<<(*sol_data)(0,node_id)<<" "<<(*sol_data)(1,node_id)<<" "<<(*sol_data)(2,node_id)<<endl;
       }
       input.close();
   }

   void PoissonSolcurRead(hier::Patch<NDIM> & patch,const std::string filename){
       tbox::Pointer<pdat::NodeData<NDIM, double> >
               sol_data = patch.getPatchData(d_sol_id);

       std::fstream input(filename);

       if(!input) {
           TBOX_ERROR("fail to open file: "<< filename << endl);
       }

       string line;

       getline(input,line);

       line.clear();

       while( getline(input, line) )
       {
           stringstream buf(line);
           int node_id;
           double node_data;

           //read and convert data
           buf >> node_id;
           buf >> node_data;

           (*sol_data)(node_id-3*(node_id/3),node_id/3) = node_data;
           //cout<<"sol_data:"<<(*sol_data)(0,node_id)<<" "<<(*sol_data)(1,node_id)<<" "<<(*sol_data)(2,node_id)<<endl;
       }
       input.close();
   }

   /** 
    * @brief 获取自由度信息
    * 
    * 
    * @return 指针，指向自由度信息
    */
   tbox::Pointer<solv::DOFInfo<NDIM> > getDOFInfo()
   {
     return d_dof_info;
   }

   class Parameters{
       public:
       double MAX_DOP;
       void setmaxdop(double value){MAX_DOP = value;}
       double getmaxdop(){return MAX_DOP;}
   };
   Parameters parameter;

 private:

  /*!@brief 对象名.  */
  string d_object_name;

  /// 自由度信息
  tbox::Pointer<solv::DOFInfo<NDIM> > d_dof_info;
  tbox::Pointer<solv::DOFInfo<NDIM> > d_dof_semi_info;
  tbox::Pointer<solv::DOFInfo<NDIM> > d_dof_Vogamma_info;

  /// 有限元计算的形函数类型, 单元类型, 积分器类型. 
  string d_shape_func_type;
  string d_element_type;
  string d_integrator_type;
  tbox::Array<string> d_constraint_types;
  tbox::Array<int>    d_constraint_marks;
  
  /** 对应于变量的id, 通过id应户可以获取变量数据 **/

  /** 8.11 sol**/
  int d_sol_id;
  /** current conservation equation **/
  int d_curconserve_mat_id;
  int d_curconserve_fval_id;
  int d_curconserve_sol_id;

  /** poisson equation **/
  int d_poisson_mat_id;
  int d_poisson_fval_id;
  int d_poisson_sol_id;
  int d_poisson_stetch_id;
  int d_poisson_sol_cur_id;
  int d_poisson_plot_id;

  int d_poisson_newton_mat_id;
  int d_poisson_newton_fval_id;
  int d_poisson_newton_sol_id;

  /** carrier transport-electron **/
  int d_carn_mat_id;
  int d_carn_fval_id;
  int d_carn_stetch_id;
  int d_carn_sol_cur_id;
  int d_logn_sol_id;
  int d_carn_plot_id;
  int d_logn_plot_id;

  int d_fermin_sol_id;
  int d_fermin_stetch_id;
  int d_fermin_plot_id;

  int d_nc_id;

  /** carrier transport-hole **/
  int d_carh_mat_id;
  int d_carh_fval_id;
  int d_carh_stetch_id;
  int d_carh_sol_cur_id;
  int d_logh_sol_id;
  int d_carh_plot_id;
  int d_logh_plot_id;

  int d_fermih_sol_id;
  int d_fermih_stetch_id;
  int d_fermih_plot_id;

  int d_nv_id;

  /** unscaled variable **/
  int d_poisson_unscal_id;
  int d_carn_unscal_id;
  int d_carh_unscal_id;
  int d_fermin_unscal_id;
  int d_fermih_unscal_id;

  //Vogamma-x
  int d_Vogamma_newton_mat_id;
  int d_Vogamma_newton_fval_id;
  int d_Vogamma_newton_sol_id;

  int d_Vogamma_sol_cur_id;

  int d_Vogamma_plot_id;
  int d_Vogammaplus_plot_id;
  int d_VogammaHplus_plot_id;
  int d_VogammaH2plus_plot_id;
  int d_VogammaH_plot_id;
  int d_VogammaH2_plot_id;


  int d_Vogamma_id;
  int d_VogammaPlus_id;
  int d_VogammaHPlus_id;
  int d_VogammaH2Plus_id;
  int d_VogammaH_id;
  int d_VogammaH2_id;

  int d_delta_Vogamma_id;
  int d_delta_VogammaPlus_id;
  int d_delta_VogammaHPlus_id;
  int d_delta_VogammaH2Plus_id;
  int d_delta_VogammaH_id;
  int d_delta_VogammaH2_id;

  //Vodelta-x
  int d_Vodelta_newton_mat_id;
  int d_Vodelta_newton_fval_id;
  int d_Vodelta_newton_sol_id;

  int d_Vodelta_id;
  int d_VodeltaPlus_id;
  int d_VodeltaHPlus_id;
  int d_VodeltaH2Plus_id;
  int d_VodeltaH_id;
  int d_VodeltaH2_id;

  int d_delta_Vodelta_id;
  int d_delta_VodeltaPlus_id;
  int d_delta_VodeltaHPlus_id;
  int d_delta_VodeltaH2Plus_id;
  int d_delta_VodeltaH_id;
  int d_delta_VodeltaH2_id;

  int d_E_vec_id;
  int d_E_cell_id;
  int d_Jn_vec_id;
  int d_Jp_vec_id;
  int d_J_vec_id;
  int d_J_edge_id;

  /** cell to entity no. relationship **/
  int d_cell_entity_id;
  int d_cell_mat_id;
  int d_node_entity_id;
  int d_node_mat_id;
  int d_semi_node_id;
  int d_metal_node_id;
  int d_sio2_node_id;
  int d_siclass_node_id;
  int d_contact_node_id;

  int d_interface_node_id;
  int d_heter_node_id;
  int d_hnode_entity_id;
  int d_hnode_mat_id;
  int d_cell_hnode_id;

  int d_semiinsbc_node_id;

  int d_face_entity_id;
  int d_face_mat_id;

  /** dopant **/
  int d_dopant_id;
  int d_maxdoping_id;

  int d_ND_id;
  int d_NA_id;

  /** polarized charge **/
  int d_polarizedp_id;

  /** plot node variables: for debug **/
  int d_node_var_plot_id;

  /** carrier density scalor **/
  double ni_o;
  double thermo_volt_o;

  /** polarization charge density **/
  double polarized_p;

  /** time step **/
  double Dt;

  /** solution storage **/
  int d_carn_sol_1_id;
  int d_carh_sol_1_id;
  int d_fn_sol_1_id;
  int d_fh_sol_1_id;
  int d_fn_sol_0_id;
  int d_fh_sol_0_id;
  int d_phi_sol_0_id;
  int d_phi_sol_1_id;

  //bas vector
  int d_bas_integ_id;
  int d_bas_rhs_id;

  /** material list **/
  tbox::Array<material_pro> material_list;
  int n_mat;
  int n_entity;

  /** boundary assignment **/
  tbox::Array<electrobc> EBoundassignment;
  tbox::Array<carrierbc> CBoundassignment;
  /** dopant list **/
  tbox::Array<double> dopant_list;
  tbox::Array<doplist> Dopassignment;

  /** entity-mat relation **/
  /** dopant list **/
  tbox::Array<int> entity_mat_list;

  tbox::Pointer<tbox::Database> GetPointor_db;

  /** cell_point and cell_edge information for semi**/
  int cell_point;
  int cell_edge;
  /** dopant_file_name**/
  string dopant_file_name;
  string sol_file_name;
  int num_sol_file_name;
  string PoissonSolcur_filename;
  int num_PoissonSolcur_filename;

  bool transient_state;
  double constraint_value;

  bool car_bound_mode;
  bool enable_selfDop;
  bool save_bas;
  bool save_carrier;
  bool save_sol;
  bool enable_tid;
  bool get_InitialSol;
  bool get_PoissonSolcur;

  double Not;
  string DD_type;

  /* unit scaling factor*/
  bool enable_unit_scalling;
  double T;
};
#endif
