/**
 * @file OVI_Controller.h
 * @author 李长骏 提出者，理论推导和公式实现
 * @author 王虓 (2569209027@qq.com) 代码编写
 * @author 赵飞 西安交通大学机器人所 指导老师
 * @brief 这个代码是一个结合用最优控制实现变阻抗控制的一个库，名字叫OVIC(optimal variable impedance controller)
 * @version 0.1
 * @date 2022-12-14
 * 
 * @copyright Copyright (c) 2022 
 * 
 */
#ifndef OVI_CONTROLLER_H
#define OVI_CONTROLLER_H

#include "include_common.h"

namespace xj_dy_ns
{
    /**
     * @brief 看文件开头的描述
     * 
     */
    class OVI_Controller
    {
    private:
        /* data */
    public:
        OVI_Controller();
        bool init_param();//初始化传入参数
        static std::vector<Eigen::MatrixXd> P_matrix_cal(double dt,
                                                        std::vector<Eigen::MatrixXd> Md,
                                                        std::vector<Eigen::MatrixXd> Dd,
                                                        std::vector<Eigen::MatrixXd> Kd,
                                                        Eigen::MatrixXd S,
                                                        Eigen::MatrixXd Q,
                                                        Eigen::MatrixXd R
                                                        );
        static Eigen::VectorXd OVI_tor_cal(Eigen::VectorXd epsil,
                                            Eigen::VectorXd d_epsil,
                                            Eigen::MatrixXd P,
                                            double dt);

        ~OVI_Controller();
    };
    
}








#endif
