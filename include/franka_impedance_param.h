#ifndef FRANKA_IMPEDANCE_PARAM_H 
#define FRANKA_IMPEDANCE_PARAM_H

#include <iostream>
#include <Eigen/Eigen>
#include <ros/ros.h>

namespace xj_dy_ns
{

    


    class franka_impedance_param
    {
    private:
        /* data */
    public:
        franka_impedance_param();
        ~franka_impedance_param();
        static void set_cartation_param6X6_(const double& x,
                                    const double& r,
                                    Eigen::Matrix<double,6,6>& matrix_6X6);//设置成一样
        static void set_cartation_param6X6_(const Eigen::Vector3d& x,
                                        const Eigen::Vector3d& r,
                                        Eigen::Matrix<double,6,6>& matrix_6X6);
        static void set_cartation_impedance_(const Eigen::Vector3d & Kx,
                                        const Eigen::Vector3d& Kr,
                                        Eigen::Matrix<double,6,6>& matrix_K,
                                        Eigen::Matrix<double,6,6>& matrix_damping);//根据刚度设置阻尼
        static void set_cartation_impedance_(const double& Kx,
                                        const double& Kr,
                                        Eigen::Matrix<double,6,6>& matrix_K,
                                        Eigen::Matrix<double,6,6>& matrix_damping);//根据刚度设置阻尼
        static void set_cartation_impedance_(const double& Kx,
                                        const double& Kr,
                                        Eigen::Matrix<double,6,6>& matrix_K,
                                        Eigen::Matrix<double,6,6>& matrix_damping,
                                        const Eigen::Matrix<double,6,6>&  Lambda_d);//根据刚度设置阻尼
        static void set_cartation_impedance_(const double& Kx,
                                        const double& Kr,
                                        Eigen::MatrixXd& matrix_K,
                                        Eigen::MatrixXd& matrix_damping,
                                        const Eigen::MatrixXd&  Lambda_d);//根据刚度设置阻尼

        static void set_cartation_param6X6_(const double& x,
                                        Eigen::Matrix<double,6,6>& matrix_6X6);
        static void set_cartation_param6X6_(const double& x,
                                        Eigen::MatrixXd& matrix_6X6);
        static void set_cartation_impedance_(const double& x,
                                        Eigen::Matrix<double,6,6>& matrix_K,
                                        Eigen::Matrix<double,6,6>& matrix_D);
        static void set_cartation_impedance_(const double& x,
                                        Eigen::MatrixXd& matrix_K,
                                        Eigen::MatrixXd& matrix_D);
        static void set_joint_impedance_(const double& x,
                                        Eigen::MatrixXd& matrix_K,
                                        Eigen::MatrixXd& matrix_D);
    };


}



#endif