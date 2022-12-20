
#include "impedance_controller.h"
namespace xj_dy_ns
{


    ImpedanceController::ImpedanceController(/* args */)
    {

    }
    ImpedanceController::ImpedanceController(int dof)
    {
        init_param(dof);
    }

    /**
     * @brief 初始化成员参数矩阵的大小（阶数），
     * 
     * @param dof 
     * @return true 
     * @return false 
     */
    bool ImpedanceController::init_param(int dof)
    {
        this->DOF_ = dof;
    }
    /**
     * @brief 通过完全输入参数来计算机器人的阻抗控制的输出
     * 
     * @param Lanmbda_d 期望笛卡尔空间下的惯量
     * @param D_d 笛卡尔空间下的期望阻尼
     * @param K_d 笛卡尔空间下的期望刚度
     * @param inv_jacobe 雅可比的逆矩阵，因为如果是冗余自由度的机械臂，需要在动力学中计算伪逆
     * @param M_q 在关节空间下描述的机器人本体惯量
     * @param jacobe 末端在世界坐标系下的雅可比矩阵
     * @param d_jacobe 雅可比矩阵的微分
     * @param x_d 笛卡尔空间下的期望位置
     * @param dx_d 笛卡尔空间下的期望速度
     * @param ddx_d 笛卡尔空间下的期望加速度
     * @param x 笛卡尔空间下的工具坐标系当前位置
     * @param dx 笛卡尔空间下工具坐标系当前速度
     * @param F_ext 笛卡尔空间下工具坐标系当前所受外力，方向为受力方向，不是六维力传感器直接读出来的数据，这两个差一个负号
     * @param tor_C 关节空间下的科氏力离心力补偿力矩，不是真实方向！！！
     * @param tor_g 关节空间下的重力矩的补偿力矩，不是真实方向！！！
     * @return Eigen::VectorXd 
     */
    Eigen::VectorXd ImpedanceController::tau_impedance_cal(Eigen::Matrix<double,6,6> Lanmbda_d,
                                Eigen::Matrix<double,6,6> D_d,
                                Eigen::Matrix<double,6,6> K_d,
                                Eigen::Matrix<double,Eigen::Dynamic,6> inv_jacobe,
                                Eigen::MatrixXd M_q,
                                Eigen::Matrix<double,6,Eigen::Dynamic> jacobe,
                                Eigen::Matrix<double,6,Eigen::Dynamic> d_jacobe,
                                Eigen::Matrix<double,6,1> x_d,
                                Eigen::Matrix<double,6,1> dx_d,
                                Eigen::Matrix<double,6,1> ddx_d,
                                Eigen::Matrix<double,6,1> x,
                                Eigen::Matrix<double,6,1> dx,
                                Eigen::Matrix<double,6,1> F_ext,
                                Eigen::VectorXd tor_C,
                                Eigen::VectorXd tor_g
                                )
    {

    }






    
    ImpedanceController::~ImpedanceController()
    {
    }


    
}//namespace xj_dy_ns
