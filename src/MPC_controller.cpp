

#include "MPC_controller.h"

namespace xj_dy_ns
{
    /**
     * @brief 无参构造函数
     * 
     */
    MPC_controller::MPC_controller()
    {
        ;
    }

    MPC_controller::MPC_controller(Eigen::MatrixXd A,Eigen::MatrixXd B,double T)
    {
        init(A,B,T);
    }

    /**
     * @brief Construct a new mpc controller::mpc controller object
     * 
     * @param Lambda_d 期望惯量
     * @param D_d 期望阻尼
     * @param K_d 期望刚度
     * @param T 控制周期
     */
    MPC_controller::MPC_controller(Eigen::MatrixXd Lambda_d,Eigen::MatrixXd D_d,Eigen::MatrixXd K_d,double T)
    {
        impedance_init_A_B(Lambda_d,D_d,K_d,T);
    }

    /**
     * @brief 使用线性状态方程初始化
     * 
     * @param A 线性状态方程A
     * @param B 线性状态方程B
     * @param T 控制周期
     */
    void MPC_controller::init(Eigen::MatrixXd A,Eigen::MatrixXd B,double T)
    {
        u_size_=B.cols();
        int dof_A=A.cols();
        A_dof_=dof_A;
        A_=Eigen::MatrixXd::Identity(dof_A,dof_A)+T*A;
        B_=T*B;
    }

    /**
     * @brief 
     * 
     * @param Lambda_d 期望惯量
     * @param D_d 期望阻尼
     * @param K_d 期望刚度
     * @param T 控制周期
     */
    void MPC_controller::impedance_init_A_B(Eigen::MatrixXd Lambda_d,Eigen::MatrixXd D_d,Eigen::MatrixXd K_d,double T)
    {
        int dof = Lambda_d.cols();
        Eigen::MatrixXd A,B;
        A.setZero(2*dof,2*dof);
        B.setZero(2*dof,dof);
        A.topRightCorner(dof,dof).setIdentity();
        Eigen::MatrixXd Lambda_inv=Lambda_d.inverse();
        A.bottomLeftCorner(dof,dof)=-Lambda_inv*K_d;
        A.bottomRightCorner(dof,dof)=-Lambda_inv*D_d;
        B.bottomRows(dof)=-Lambda_inv;
        init(A,B,T);
    }

    /**
     * @brief 直接用惯量阵计算MPC参数
     * 
     * @param Lambda_now 
     * @param T 
     */
    void MPC_controller::robot_M_init_A_B(Eigen::MatrixXd Lambda_now,double T)
    {
        int dof = Lambda_now.cols();
        Eigen::MatrixXd A,B;
        A.setZero(2*dof,2*dof);
        B.setZero(2*dof,dof);
        A.topRightCorner(dof,dof).setIdentity();
        Eigen::MatrixXd Lambda_inv=Lambda_now.inverse();
        B.bottomRows(dof)=-Lambda_inv;
        init(A,B,T);
    }

    
    /**
     * @brief 
     * 
     * @param q 恒定状态权重参数
     * @param r 恒定输出权重参数
     * @param N 预测步长
     * @param X 当前状态
     * @return Eigen::VectorXd 
     */
    Eigen::VectorXd MPC_controller::u_MPC_cal_no_constrain(double q,double r,int N,Eigen::VectorXd X)
    {
        //状态权重
        Eigen::MatrixXd Q=q*Eigen::MatrixXd::Identity(A_dof_,A_dof_);
        // Q.bottomRightCorner(A_dof_/2,A_dof_/2).setZero();
        //输出权重
        Eigen::MatrixXd R=r*Eigen::MatrixXd::Identity(u_size_,u_size_);
        //求解K_MPC
        Eigen::MatrixXd K_MPC=K_MPC_cal_no_constrain(Q,R,N);
        //求解输出
        // printf("\033[1;31;40m =   \n");
        // std::cout<<"计算K_MPC = "<< K_MPC<<std::endl;
        // std::cout<<"计算K_MPC 行= "<< K_MPC.rows()<<"列="<<K_MPC.cols()<<std::endl;
        // std::cout<<"显示Q= "<< Q<<"R="<<R<<std::endl;

        // printf(" \033[0m \n");
        Eigen::VectorXd u=-K_MPC*X;
        return u;
    }

    /**
     * @brief 无约束MPC求解 K_MPC
     * 
     * @param Q 
     * @param R 
     * @param N 
     * @return Eigen::MatrixXd 
     */
    Eigen::MatrixXd MPC_controller::K_MPC_cal_no_constrain(Eigen::MatrixXd Q,Eigen::MatrixXd R,int N)
    {
        //初始化大小
        Eigen::MatrixXd F=Eigen::MatrixXd::Zero(N*A_dof_,A_dof_);
        Eigen::MatrixXd Phi=Eigen::MatrixXd::Zero(N*A_dof_,N*this->u_size_);
        Eigen::MatrixXd Q_mathcal=Eigen::MatrixXd::Identity(N*A_dof_,N*A_dof_);
        Eigen::MatrixXd R_mathcal=Eigen::MatrixXd::Identity(N*u_size_,N*u_size_);
        Eigen::MatrixXd A_n=Eigen::MatrixXd::Identity(A_dof_,A_dof_);
        Eigen::MatrixXd AB=B_;
        for (int i = 0; i < N; i++)
        {
            //定义花写Q和R
            Q_mathcal.block(i*A_dof_,i*A_dof_,A_dof_,A_dof_)=Q;
            R_mathcal.block(i*u_size_,i*u_size_,u_size_,u_size_)=R;
            //求解Phi
            for (int j = 0; j < N-i; j++)
            {
                Phi.block(i*A_dof_+j*A_dof_,u_size_*j,A_dof_,u_size_)=AB;
            }
            AB=A_*AB;

            //求解F
            A_n=A_*A_n;
            F.block(i*A_dof_,0,A_dof_,A_dof_)=A_n;
        }
        // printf("\033[1;31;40m =   \n");
        // std::cout<<"计算K_MPC = "<< K_MPC<<std::endl;
        // std::cout<<"计算Phi 行= "<< Phi.rows()<<"列="<<Phi.cols()<<std::endl;
        // std::cout<<"计算Phi= "<< Phi<<std::endl<<"F="<<F<<std::endl;
        // std::cout<<"计算A_= "<< A_<<std::endl;

        // std::cout<<"计算F 行= "<< F.rows()<<"列="<<F.cols()<<std::endl;
        // std::cout<<"A_dof_= "<< A_dof_<<std::endl;

        // printf(" \033[0m \n");
        Eigen::MatrixXd tmp1=(Phi.transpose()*Q_mathcal*Phi+R_mathcal);
        tmp1=tmp1.inverse();
        Eigen::MatrixXd tmp2=(Phi.transpose()*Q_mathcal*F);

        Eigen::MatrixXd K_MPC=Eigen::MatrixXd::Identity(u_size_,N*A_dof_)*tmp1*tmp2;
        return K_MPC;
    }

    

};