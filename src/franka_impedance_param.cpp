#include "franka_impedance_param.h"

namespace xj_dy_ns
{

    
    franka_impedance_param::franka_impedance_param()
    {
        ;
    }
    
    franka_impedance_param::~franka_impedance_param()
    {
        ;
    }


       /**
    * @brief 通过向量设置刚度矩阵
    * @param x 迪卡尔坐标下线性刚度向量或者阻尼
    * @param r 迪卡尔坐标系下的旋转刚度向量或者阻尼
    * @param matrix_6X6 要修改的刚度活阻尼矩阵
    */
    void franka_impedance_param::set_cartation_param6X6_(const Eigen::Vector3d& x,
                                        const Eigen::Vector3d& r,
                                        Eigen::Matrix<double,6,6>& matrix_6X6)
    {
        matrix_6X6 = Eigen::Matrix<double,6,6>::Identity();
        for (int i = 0; i < 3; i++)
        {
            matrix_6X6(i,i) = x(i);
            matrix_6X6(i+3,i+3) = r(i);
        }
    }

       /**
    * @brief 通过向量设置刚度矩阵
    * @param x 迪卡尔坐标下线性刚度向量或者阻尼
    * @param matrix_6X6 要修改的刚度或者阻尼矩阵
    */
    void franka_impedance_param::set_cartation_param6X6_(const double& x,
                                        Eigen::Matrix<double,6,6>& matrix_6X6)
    {
        matrix_6X6 = Eigen::Matrix<double,6,6>::Identity()*x;
    }

       /**
    * @brief 通过向量设置刚度矩阵
    * @param x 迪卡尔坐标下线性刚度向量或者阻尼
    * @param matrix_K 要修改的刚度
    * @param matrix_D 
    */
    void franka_impedance_param::set_cartation_impedance_(const double& x,
                                        Eigen::Matrix<double,6,6>& matrix_K,
                                        Eigen::Matrix<double,6,6>& matrix_D)
    {
        set_cartation_param6X6_(x,matrix_K);
        matrix_D = matrix_K.cwiseSqrt() *2*0.707;
    }


   /**
    * @brief 通过数值设置刚度或阻尼矩阵,所有线性参数都相等,所有旋转参数也都相等
    * @param x 迪卡尔坐标下线性数值刚度或者阻尼
    * @param r 迪卡尔坐标系下的旋转刚度或者阻尼
    * @param matrix_6X6 要修改的刚度或者阻尼矩阵
    */
    void franka_impedance_param::set_cartation_param6X6_(const double& x,
                                const double& r,
                                Eigen::Matrix<double,6,6>& matrix_6X6)//刚度设置成一样
    {
        matrix_6X6.topLeftCorner(3,3) = Eigen::Matrix3d::Identity()*x;
        matrix_6X6.bottomRightCorner(3,3) = Eigen::Matrix3d::Identity()*r;
    }


   /**
    * @brief 通过刚度来设置阻尼参数,使用最佳阻尼比0.707来设置
    * @param x 迪卡尔坐标下线性数值刚度
    * @param r 迪卡尔坐标系下的旋转刚度
    * @param matrix_K 要修改的刚度矩阵
    * @param matrix_damping 要修改的阻尼矩阵
    */
    void franka_impedance_param::set_cartation_impedance_(const Eigen::Vector3d& Kx,
                                        const Eigen::Vector3d& Kr,
                                        Eigen::Matrix<double,6,6>& matrix_K,
                                        Eigen::Matrix<double,6,6>& matrix_damping)//根据刚度设置阻尼
    {
        set_cartation_param6X6_(Kx,Kr,matrix_K);
        // Damping ratio = 0.707    2*0.707
        matrix_damping.setIdentity();
        matrix_damping = matrix_K.cwiseSqrt() *2*0.707;//开方
    }

   /**
    * @brief 通过刚度来设置阻尼参数,使用最佳阻尼比0.707来设置
    * @param x 迪卡尔坐标下线性刚度向量
    * @param r 迪卡尔坐标系下的旋转刚度向量
    * @param matrix_K 要修改的刚度矩阵
    * @param matrix_damping 要修改的阻尼矩阵
    */
    void franka_impedance_param::set_cartation_impedance_(const double& Kx,
                                        const double& Kr,
                                        Eigen::Matrix<double,6,6>& matrix_K,
                                        Eigen::Matrix<double,6,6>& matrix_damping)//根据刚度设置阻尼
    {
        matrix_K.setIdentity();
        matrix_K.topLeftCorner(3,3) = Eigen::Matrix3d::Identity()*Kx;
        matrix_K.bottomRightCorner(3,3) = Eigen::Matrix3d::Identity()*Kr;
        matrix_damping.setIdentity();
        matrix_damping = 2*0.707*matrix_K.cwiseSqrt();
    }

    /**
     * @brief joint刚度
     * 
     * @param x 
     * @param matrix_K 
     * @param matrix_D 
     */
    void franka_impedance_param::set_joint_impedance_(const double& x,
                                        Eigen::MatrixXd& matrix_K,
                                        Eigen::MatrixXd& matrix_D)
    {
        matrix_K = x * matrix_K.setIdentity();
        matrix_D = 2*0.707*matrix_K.cwiseSqrt();
    }

}