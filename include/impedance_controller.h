/**
 * @file impedance_controller.h
 * @author 王虓 (2569209027@qq.com)
 * @brief 阻抗控制器
 * @version 0.1
 * @date 2022-12-16
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#ifndef IMPEDANCE_CONTROLLER_H
#define IMPEDANCE_CONTROLLER_H

#include "include_common.h"
#include "math_wx.h"


namespace xj_dy_ns
{
    class ImpedanceController
    {
    private:
        Eigen::Matrix<double,6,6> Lambda_d_,Dd_,Kd;//笛卡尔阻抗参数
        Eigen::MatrixXd Mq;
        int DOF_;
    public:
        ImpedanceController(int dof);
        ImpedanceController();

        bool init_param(int dof);
        void set_Mq(Eigen::MatrixXd Mq);

        //静态成员函数可以不用生成对象就调用
        static Eigen::VectorXd tau_impedance_cal(Eigen::Matrix<double,6,6> Lanmbda_d,
                                        Eigen::Matrix<double,6,6> D_d,
                                        Eigen::Matrix<double,6,6> K_d,
                                        Eigen::Matrix<double,Eigen::Dynamic,6> inv_jacobe,
                                        Eigen::MatrixXd M_q,
                                        Eigen::Matrix<double,6,Eigen::Dynamic> jacobe,
                                        Eigen::Matrix<double,6,Eigen::Dynamic> d_jacobe,
                                        Eigen::Matrix<double,6,1> x_err,
                                        Eigen::Matrix<double,6,1> dx_d,
                                        Eigen::Matrix<double,6,1> ddx_d,
                                        Eigen::Matrix<double,6,1> dx,
                                        Eigen::VectorXd dq,
                                        Eigen::Matrix<double,6,1> F_ext,
                                        Eigen::VectorXd tor_C,
                                        Eigen::VectorXd tor_g
                                        );
        
        //计算阻抗控制计算出来的力矩
        static Eigen::VectorXd tau_impedance_cal(Eigen::Matrix<double,6,6> Lanmbda_d,
                                        Eigen::Matrix<double,6,6> D_d,
                                        Eigen::Matrix<double,6,6> K_d,
                                        Eigen::Matrix<double,Eigen::Dynamic,6> inv_jacobe,
                                        Eigen::MatrixXd M_q,
                                        Eigen::Matrix<double,6,Eigen::Dynamic> jacobe,
                                        Eigen::Matrix<double,6,Eigen::Dynamic> d_jacobe,
                                        Eigen::Matrix<double,6,1> x_err,
                                        Eigen::Matrix<double,6,1> dx_d,
                                        Eigen::Matrix<double,6,1> ddx_d,
                                        Eigen::Matrix<double,6,1> dx,
                                        Eigen::VectorXd dq,
                                        Eigen::Matrix<double,6,1> F_ext,
                                        Eigen::VectorXd tor__C_G
                                        );

        //计算阻抗控制计算出来的力矩，加入期望力
        static Eigen::VectorXd tau_impedance_cal(Eigen::Matrix<double,6,6> Lanmbda_d,
                                        Eigen::Matrix<double,6,6> D_d,
                                        Eigen::Matrix<double,6,6> K_d,
                                        Eigen::Matrix<double,Eigen::Dynamic,6> inv_jacobe,
                                        Eigen::MatrixXd M_q,
                                        Eigen::Matrix<double,6,Eigen::Dynamic> jacobe,
                                        Eigen::Matrix<double,6,Eigen::Dynamic> d_jacobe,
                                        Eigen::Matrix<double,6,1> x_err,
                                        Eigen::Matrix<double,6,1> dx_d,
                                        Eigen::Matrix<double,6,1> ddx_d,
                                        Eigen::Matrix<double,6,1> dx,
                                        Eigen::VectorXd dq,
                                        Eigen::Matrix<double,6,1> F_ext,
                                        Eigen::Matrix<double,6,1> F_d,
                                        Eigen::VectorXd tor__C_G,
                                        Eigen::Matrix<double,6,1> &xr_err,
                                        Eigen::Matrix<double,6,1> &dxr_err,
                                        double dt
                                        );



        static Eigen::VectorXd tau_impedance_cal_axis(Eigen::Matrix<double,6,6> Lanmbda_d,
                                        Eigen::Matrix<double,6,6> D_d,
                                        Eigen::Matrix<double,6,6> K_d,
                                        Eigen::MatrixXd M_q,
                                        Eigen::Matrix<double,6,Eigen::Dynamic> jacobe,
                                        Eigen::Matrix<double,6,Eigen::Dynamic> d_jacobe,
                                        Eigen::Matrix4d T_d,
                                        Eigen::Matrix<double,3,1> dx_d,
                                        Eigen::Matrix3d dR_d,
                                        Eigen::Matrix<double,3,1> ddx_d,
                                        Eigen::Matrix3d ddR_d,
                                        Eigen::Matrix<double,6,1> F_ext,
                                        Eigen::Matrix<double,6,1> F_d,
                                        Eigen::VectorXd tor__C_G,
                                        Eigen::Matrix<double,3,1> dx,
                                        Eigen::VectorXd dq,
                                        Eigen::Matrix4d T_now,
                                        Eigen::Matrix<double,6,1> &xr_err,
                                        Eigen::Matrix<double,6,1> &dxr_err,
                                        Eigen::Ref<Eigen::MatrixXd> derr_new,
                                        double dt
                                        );

        //这个是目标减去当前位姿态来返回笛卡尔空间下的位置误差和轴角误差
        static Eigen::Matrix<double,6,1> x_err_cal(Eigen::Matrix4d T_d,Eigen::Matrix4d T_now);
        
        //混在一起的轴角转化为旋转矩阵
        static Eigen::Matrix3d axis2rot(Eigen::Vector3d n_theta);

        ~ImpedanceController();
    };   

    
}

#endif
