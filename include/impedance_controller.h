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
        Eigen::MatrixXd Md_,Dd_,Kd_;//关节阻抗控制参数
        Eigen::Matrix<double,6,6> Mxd_,Dxd_,Kxd;//笛卡尔阻抗参数
        int DOF_;
    public:
        ImpedanceController(int dof);
        ImpedanceController();

        bool init_param(int dof);
        ~ImpedanceController();
    };
    

    
}

#endif
