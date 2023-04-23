#include "include_common.h"


namespace xj_dy_ns
{
  class MPC_controller
  {
  private:
    /* data */
    //这里的A和B是离散化后的A和B
    Eigen::MatrixXd A_;
    Eigen::MatrixXd B_;
    int u_size_;
    int A_dof_;
  public:
    MPC_controller(/* args */);
    MPC_controller(Eigen::MatrixXd A,Eigen::MatrixXd B,double T);
    MPC_controller(Eigen::MatrixXd Lambda_d,Eigen::MatrixXd D_d,Eigen::MatrixXd K_d,double T);

    void init(Eigen::MatrixXd A,Eigen::MatrixXd B,double T);    //输入的是线性A和B。计算出来给内部赋值的事离散的A和B
    void impedance_init_A_B(Eigen::MatrixXd Lambda_d,Eigen::MatrixXd D_d,Eigen::MatrixXd K_d,double T);//直接用阻抗参数计算
    void robot_M_init_A_B(Eigen::MatrixXd Lambda_now,double T);//直接用惯性力计算

    ~MPC_controller();

    Eigen::VectorXd u_MPC_cal_no_constrain(double q,double r,int N,Eigen::VectorXd X);
    Eigen::VectorXd u_MPC_cartesian_cal_no_constrain(Eigen::Vector4d q,Eigen::Vector2d r,int N,Eigen::VectorXd X);

    Eigen::MatrixXd K_MPC_cal_no_constrain(Eigen::MatrixXd Q,Eigen::MatrixXd R,int N);
  };
    
};