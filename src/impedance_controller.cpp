
#include "impedance_controller.h"
/// @brief 
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
        return true;
    }


    /**
     * @brief 通过完全输入参数来计算机器人的阻抗控制的输出，这是一个静态成员函数，可以直接调用不用生成对象
     * 
     * @param Lanmbda_d 期望笛卡尔空间下的惯量
     * @param D_d 笛卡尔空间下的期望阻尼
     * @param K_d 笛卡尔空间下的期望刚度
     * @param inv_jacobe 雅可比的逆矩阵，因为如果是冗余自由度的机械臂，需要在动力学中计算伪逆
     * @param M_q 在关节空间下描述的机器人本体惯量
     * @param jacobe 末端在世界坐标系下的雅可比矩阵
     * @param d_jacobe 雅可比矩阵的微分
     * @param x_err 笛卡尔空间下的误差
     * @param dx_d 笛卡尔空间下的期望速度
     * @param ddx_d 笛卡尔空间下的期望加速度
     * @param dx 笛卡尔空间下工具坐标系当前速度
     * @param dq 关节速度
     * @param F_ext 笛卡尔空间下工具坐标系当前所受外力，方向为受力方向，不是六维力传感器直接读出来的数据，这两个差一个负号
     * @param tor_C 关节空间下的科氏力离心力补偿力矩，不是真实方向！！！也就是惯性力的反方向，或者说是与加速度方向一致
     * @param tor_g 关节空间下的重力矩的补偿力矩，不是真实方向！！！这个是重力的反方向。
     * @return Eigen::VectorXd tau_imp_cmd计算出补偿力矩的大小
     */
    Eigen::VectorXd ImpedanceController::tau_impedance_cal(Eigen::Matrix<double,6,6> Lanmbda_d,
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
                                )
    {
        Eigen::VectorXd tau_imp_cmd;
        tau_imp_cmd.setZero(tor_C.rows());//初始化补偿力矩的大小
        Eigen::MatrixXd k1 = M_q*inv_jacobe;
        Eigen::MatrixXd k2 = k1*Lanmbda_d.inverse();

        tau_imp_cmd = k1*(ddx_d-d_jacobe*dq) 
        + tor_C 
        + tor_g 
        + k2*(D_d*(dx_d-dx) 
        + K_d*(x_err))
        + (k2 - jacobe.transpose())*F_ext;
        return tau_imp_cmd;
    }

    /**
     * @brief 通过完全输入参数来计算机器人的阻抗控制的输出，这是一个静态成员函数，可以直接调用不用生成对象
     * 
     * @param Lanmbda_d 期望笛卡尔空间下的惯量
     * @param D_d 笛卡尔空间下的期望阻尼
     * @param K_d 笛卡尔空间下的期望刚度
     * @param inv_jacobe 雅可比的逆矩阵，因为如果是冗余自由度的机械臂，需要在动力学中计算伪逆
     * @param M_q 在关节空间下描述的机器人本体惯量
     * @param jacobe 末端在世界坐标系下的雅可比矩阵
     * @param d_jacobe 雅可比矩阵的微分
     * @param x_err 笛卡尔空间下的误差
     * @param dx_d 笛卡尔空间下的期望速度
     * @param ddx_d 笛卡尔空间下的期望加速度
     * @param dx 笛卡尔空间下工具坐标系当前速度
     * @param dq 关节速度
     * @param F_ext 笛卡尔空间下工具坐标系当前所受外力，方向为受力方向，不是六维力传感器直接读出来的数据，这两个差一个负号
     * @param tor__C_G 关节空间下的科氏力离心力和重力的补偿力矩，不是真实方向！！！也就是惯性力的反方向，或者说是与加速度方向一致
     * @return Eigen::VectorXd tau_imp_cmd计算出补偿力矩的大小
     */
    Eigen::VectorXd ImpedanceController::tau_impedance_cal(Eigen::Matrix<double,6,6> Lanmbda_d,
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
                                )
    {
        //测试显示输入参数程序

        // printf("\033[1;32;40m cmd_tor =   \n");
        // std::cout<<"Lanmbda_d ="<<Lanmbda_d<<std::endl;
        // std::cout<<"D_d ="<<D_d<<std::endl;
        // std::cout<<"K_d ="<<K_d<<std::endl;
        // std::cout<<"inv_jacobe ="<<inv_jacobe<<std::endl;
        // std::cout<<"jacobe ="<<jacobe<<std::endl;
        // std::cout<<"d_jacobe ="<<d_jacobe<<std::endl;
        // std::cout<<"x_err ="<<x_err<<std::endl;
        // std::cout<<"dx_d ="<<dx_d<<std::endl;
        // std::cout<<"ddx_d ="<<ddx_d<<std::endl;
        // std::cout<<"dx ="<<dx<<std::endl;
        // std::cout<<"dq ="<<dq<<std::endl;
        // std::cout<<"F_ext ="<<F_ext<<std::endl;
        // std::cout<<"tor__C_G ="<<tor__C_G<<std::endl;
        // printf(" \033[0m \n");
        


        Eigen::VectorXd tau_imp_cmd;
        tau_imp_cmd.setZero(tor__C_G.rows());//初始化补偿力矩的大小
        Eigen::MatrixXd k1 = M_q*inv_jacobe;
        Eigen::MatrixXd k2 = k1*Lanmbda_d.inverse();

        tau_imp_cmd = k1*(ddx_d-d_jacobe*dq) 
        + tor__C_G
        + k2*(D_d*(dx_d-dx)
        + K_d*(x_err))
        + (k2 - jacobe.transpose())*F_ext;
        return tau_imp_cmd;
    }


    /**
     * @brief 通过目标位姿和当前位姿，来计算误差
     * 
     * @param T_d 目标其次矩阵
     * @param T_now 当前齐次矩阵
     * @return Eigen::Matrix<double,6,1> 返回的是误差向量。
     */
    Eigen::Matrix<double,6,1> ImpedanceController::x_err_cal(Eigen::Matrix4d T_d,Eigen::Matrix4d T_now)
    {
        Eigen::Matrix3d R_err = T_d.topLeftCorner(3,3)*((T_now.topLeftCorner(3,3)).transpose());
        double xita = acos((R_err.trace()-1)/2);//计算角度
        Eigen::Vector3d n = Eigen::Vector3d::Zero();
        n(0)=R_err(2,1)-R_err(1,2);
        n(1)=R_err(0,2)-R_err(2,0);
        n(2)=R_err(1,0)-R_err(0,1);

        n = n/(2*sin(xita));
        Eigen::Matrix<double,6,1> err = Eigen::Matrix<double,6,1>::Zero();
        err.bottomRows(3) = n*xita;
        err.topRows(3) = T_d.topRightCorner(3,1) - T_now.topRightCorner(3,1);
        //调试输出
        // printf("\033[1;31;40m err =   \n");
        // std::cout<<err<<std::endl;
        // printf("   \033[0m \n");
        return err;
    }

    
    ImpedanceController::~ImpedanceController()
    {
    }


    
}//namespace xj_dy_ns
