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
        Eigen::VectorXd tau_impedance_cal(Eigen::Matrix<double,6,6> Lanmbda_d,
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
                                        );
        ~ImpedanceController();
    };
    

    
}

#endif
