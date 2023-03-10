/**
 * @file robot_dynamic.h
 * @author 王虓 (2569209027@qq.com)
 * @author 赵飞 西安交通大学机器人所指导老师
 * @brief 代码用来计算串联机械臂的运动学和动力学，内有雅可比矩阵计算，坐标变换，正向运动学，牛顿欧拉法迭代，拉格朗日法计算惯性力。
 *          可操作度计算以及可操作度优化力矩
 * @version 0.1
 * @date 2022-12-14
 * 
 * @copyright Copyright (c) 2022 西安交通大学机器人所 指导老师：赵飞
 * 
 */
#include "include_common.h"



namespace xj_dy_ns
{
    class Robot_dynamic
    {
    private:
        /* data */
    
    Eigen::Matrix<double,Eigen::Dynamic,7> DH_table;    //DH未换算单位的表
    Eigen::Matrix<double,Eigen::Dynamic,7> DH_table_red;//DH表换算成可以直接计算的表
    Eigen::Matrix<double,Eigen::Dynamic,1> q_now;
    Eigen::Matrix<double,Eigen::Dynamic,1> gq;
    std::vector<Eigen::Matrix<double,4,4>> T_;
    std::vector<Eigen::Matrix<double,3,1>> Pc;//质心向量
    std::vector<double> m_;
    std::vector<Eigen::Matrix<double,3,3> > Ic_;//转动惯量
    Eigen::Matrix<double,3,1> Pc_eff_;//末端执行器质心向量
    double m_eff_;//末端执行器质量
    Eigen::Matrix<double,3,3> Ic_eff_;//末端执行器转动惯量在质心处的表达

    

    Eigen::Matrix<double,Eigen::Dynamic,1> G_;//补偿重力的关节力矩值
    std::vector<Eigen::Matrix<double,3,1>> F_T_EE;
    Eigen::Matrix<double,4,4> T_tool_;
    
    std::vector<Eigen::Matrix<double,4,4>> _0T_i;//世界坐标系下
    Eigen::Matrix<double,4,4> _0T_tool;//世界坐标系下的末端齐次变换矩阵

    std::vector<Eigen::Matrix<double,3,1>> v_;//第i个坐标系原点的速度
    std::vector<Eigen::Matrix<double,3,1>> w_;
    std::vector<Eigen::Matrix<double,3,1>> a_;
    std::vector<Eigen::Matrix<double,3,1>> dw_;
    Eigen::Matrix<double,Eigen::Dynamic,1> dq_;//关节速度
    Eigen::Matrix<double,Eigen::Dynamic,1> dq_last_;//关节速度

    Eigen::Matrix<double,Eigen::Dynamic,1> ddq_;//关节加速度
    Eigen::Matrix<double,Eigen::Dynamic,1> ddq_last_;//上一次关节加速度

    std::vector<Eigen::Matrix<double,6,Eigen::Dynamic>> jacobi_ci_;
    Eigen::Matrix<double,6,Eigen::Dynamic> jacobi_;//末端雅可比矩阵
    Eigen::Matrix<double,6,Eigen::Dynamic> d_jacobi_;//末端雅可比矩阵的导数
    

    Eigen::Matrix<double,Eigen::Dynamic,1> tor_CpM_neton_;
    Eigen::Matrix<double,Eigen::Dynamic,1> tor_CpM_neton_last_;
    std::vector<double> motor_I_; //电机转子的惯量

    Eigen::MatrixXd M_q_;
    Eigen::VectorXd f_mu_;//摩擦系数
    Eigen::VectorXd f_s_;//静摩擦力
    Eigen::VectorXd tor_friction_;
    Eigen::VectorXd tor_CpG_neton_;//科氏力离心力和重力矩方向
    double manipulabilityIndex_position_;
    Eigen::Matrix<double,Eigen::Dynamic,6> jacobe_pse_inv_;
    Eigen::Matrix<double,6,6> Lambda_now_;
    Eigen::VectorXd limit_max;
    Eigen::VectorXd limit_min;
    Eigen::MatrixXd nullspace_jacobi_;

    Eigen::Matrix<double,6,1> F_ext_;


    bool init_size_flag=false;
    int DOF_ = 0;
    public:
        Robot_dynamic();//需要调用resize_param函数来初始化大小
        Robot_dynamic(const std::string& url,const int& DOF);
        Robot_dynamic(int dof);//自带resize_param初始化参数大小
        ~Robot_dynamic();
        bool read_dynamic(const std::string& url,const int& DOF);
        void read_dynamic_type_Gong(Eigen::VectorXd param,int dof);
        void read_dynamic_type_franka(Eigen::MatrixXd param,int dof);
        void set_DH_table(Eigen::Matrix<double,Eigen::Dynamic,7> dh_table,int dof);

        void set_q_now(Eigen::Matrix<double,Eigen::Dynamic,1> q_now);
        Eigen::Matrix<double,4,4> Ti_cal(const int n,float theta);
        std::vector<Eigen::Matrix<double,3,1>> neton_iter(
            int i,
            Eigen::Matrix<double,3,1> F_ip1,
            Eigen::Matrix<double,3,1> T_ip1,
            Eigen::Matrix<double,3,1> F_i,
            Eigen::Matrix<double,3,1> T_i);//牛顿欧拉法计算这一次的力
        void T_cal();
        std::vector<Eigen::Matrix4d> T_cal(Eigen::VectorXd q);
        Eigen::Matrix<double,3,3> get_R(int i);
        void tor_gravity_cal();//计算每个关节的重力项
        Eigen::Matrix<double,Eigen::Dynamic,1> get_G_();//获取重力的关节力矩值
        Eigen::VectorXd tor_gravity_and_Cq_cal();//获取使用牛顿欧拉法迭代计算重力+科氏力离心力
        double get_qi_now(int i);
        Eigen::VectorXd get_q_now();
        void jacobi_cal();
        void djacobe_cal();
        Eigen::Matrix<double,3,3> get_0R(int i);
        std::vector<Eigen::Matrix3d> get_0R(std::vector<Eigen::Matrix4d> T_all);
        void set_tool(Eigen::Matrix<double,4,4> T_tool);
        bool joint_is_rev(int i);
        std::vector<Eigen::Matrix<double,3,1>>  vel_iter(
            int ip1,
            Eigen::Matrix<double,3,1> vi,
            Eigen::Matrix<double,3,1> wi,
            double theta_dot_ip1
            );    //迭代计算i+1个坐标系的线速度和角速度
        
        std::vector<Eigen::Matrix<double,3,1>> a_iter(
            int ip1,
            Eigen::Matrix<double,3,1> ai,
            Eigen::Matrix<double,3,1> dwi,
            Eigen::Matrix<double,3,1> wi,
            double dtheta_ip1,
            double ddtheta_ip1);      //迭代计算加速度
        void a_cal();//计算加速度和角加速度
        std::vector<Eigen::Matrix<double,6,1>> a_cal(Eigen::VectorXd dq,Eigen::VectorXd ddq);//计算通过传入参数计算质心的加速度，多用于忽悠牛顿和欧拉

        void vel_cal();//更新计算速度
        std::vector<Eigen::Matrix<double,6,1>> vel_cal(Eigen::VectorXd dq);//通过dq计算速度，用于忽悠牛顿和欧拉

        Eigen::Matrix<double,6,1> get_vel_w_jacobe();
        std::vector<Eigen::Matrix<double,6,Eigen::Dynamic>> get_Pc_jacobe();
        Eigen::Matrix<double,6,1> get_vel_w_iter();

        void set_dq_now(Eigen::Matrix<double,Eigen::Dynamic,1> dq);//设置当前关节速度
        void ddq_cal(ros::Duration period);//通过差分法计算当前关节加速度
        void set_last_dq();//设置lastdq
        
        std::vector<Eigen::Matrix<double,3,1>> get_a_dw_now(int i);//获取第i个坐标系的线加速度和角加速度在第i个坐标系下的表达
        Eigen::Matrix<double,3,1> get_a_now(int i);//获取第i个坐标系的线加速度在第i个坐标系下的表达
        Eigen::Matrix<double,3,1> get_dw_now(int i);//获取第i个坐标系的线角速度在第i个坐标系下的表达
        std::vector<Eigen::Matrix<double,3,1>> get_i_M_C_cal(int i);
        std::vector<Eigen::Matrix<double,3,1>> get_i_M_C_cal(int i,Eigen::Matrix<double,6,1> a_dw_Pci,Eigen::Matrix<double,6,1> v_w_Pci);


        Eigen::Matrix<double,3,1> a_Pc_cal(int i);//计算第i根连杆的加速度，角加速度不用计算
        std::vector<Eigen::Matrix<double,3,1>> a_Pc_cal(std::vector<Eigen::Matrix<double,6,1>> vel_frame,
                                        std::vector<Eigen::Matrix<double,6,1>> a_dw_frame);//计算全部的质心线加速度，角加速度不用计算

        void tor_M_C_neton_cal_();
        Eigen::VectorXd tor_M_C_neton_cal_(Eigen::VectorXd dq,
                                            Eigen::VectorXd ddq);//通过给定参数计算牛顿欧拉法，用来忽悠牛顿欧拉

        Eigen::VectorXd get_tor_CpM_neton_();//获取惯性力和离心力
        Eigen::VectorXd get_tor_CpG_neton_();//获取

        Eigen::Matrix<double,Eigen::Dynamic,1> tor_filter(Eigen::VectorXd tor_,
        double hz,
        double period,
        Eigen::VectorXd* tor_last);

        static Eigen::Matrix<double,Eigen::Dynamic,1> tor_filter2(Eigen::VectorXd& tor_,
        double hz,
        double period,
        Eigen::VectorXd& tor_last);

        

        Eigen::Matrix<double,Eigen::Dynamic,1>* get_tor_CpM_neton_last_ptr();
        Eigen::Matrix<double,4,4> get_0Ti(int i);
        void resize_param(int dof);
        void jacobe_ci_cal();
        Eigen::Matrix<double,6,Eigen::Dynamic> jacobe_cal(Eigen::VectorXd q);

        Eigen::Matrix<double,6,Eigen::Dynamic> get_jacobe_tool();
        Eigen::Matrix<double,6,Eigen::Dynamic> get_djacobe_tool();


        bool get_joint_isrevolutor(int i);//判断是否是旋转关节
        Eigen::MatrixXd M_q_cal_Lagrange();//用拉格朗日方法求解惯性矩阵
        void updata_cal();//更新内部参数
        void updata_cal(Eigen::VectorXd q,
                        Eigen::VectorXd dq);//更新内部参数
        


        
        Eigen::MatrixXd get_matrix_Mq();
        void set_friction_param(Eigen::VectorXd f_param);
        void set_friction_param(Eigen::VectorXd f_s,Eigen::VectorXd f_mu);
        Eigen::VectorXd friction_cal();//内部参数计算
        Eigen::VectorXd friction_cal(Eigen::VectorXd dq);//外部参数计算
        double manipulabilityIndex_position();//计算内部可操作度
        double manipulabilityIndex_position(Eigen::VectorXd q);//计算内部可操作度
        Eigen::VectorXd manipulabilityOptimization_tor_cal(const Eigen::VectorXd& q, const double k_0);
        Eigen::VectorXd manipulabilityOptimization_tor_cal_2(const Eigen::VectorXd& q, const double k_0);//讲关节限位更新到可操作度权重里

        static Eigen::Matrix<double,Eigen::Dynamic,6> pseudo_inverse_jacobe_cal(Eigen::Matrix<double,6,Eigen::Dynamic> jacobi);

        static Eigen::Matrix<double,6,6> Lambda_now_cal(Eigen::MatrixXd Mq,Eigen::Matrix<double,6,-1> jacobi);
        Eigen::Matrix<double,6,6> get_Lambda_now();
        Eigen::Matrix<double,6,6> Lambda_now_cal();
        Eigen::Matrix<double,Eigen::Dynamic,6> get_pseudo_inverse_jacobe();
        Eigen::Matrix4d get_T_0tool();//获得末端坐标系到世界坐标系的位姿变换矩阵
        Eigen::VectorXd get_dq_now();

        static Eigen::MatrixXd nullspace_matrix_jacobi_cal(Eigen::MatrixXd jacobi_p_inv,Eigen::MatrixXd jacobi);//计算零空间矩阵
        static Eigen::MatrixXd nullspace_matrix_Nd_tau_cal(Eigen::MatrixXd jacobi,Eigen::MatrixXd Lambda,Eigen::MatrixXd Mq);
        Eigen::MatrixXd get_nullspace_matrix();//获得内部零空间矩阵

        Eigen::VectorXd limit_optimiza_tor_cal(double kd);//内部计算关节力矩限制

        void set_endeffector_dynamic_param(Eigen::Matrix<double,3,1> Pc_eff,
                                            double m_eff,
                                            Eigen::Matrix<double,3,3> Ic_eff);//设置末端执行器的动力学参数
        static Eigen::Matrix<double,6,1> F_ext_cal_by_tau_ext(Eigen::MatrixXd jacobi_r_inv,Eigen::VectorXd tau_ext);
        Eigen::Matrix<double,6,1> F_ext_cal_inner(Eigen::VectorXd tau_measure);


    //         Eigen::Matrix<double,3,1> Pc_eff_;//末端执行器质心向量
    // double m_eff_;//末端执行器质量
    // Eigen::Matrix<double,3,3> Ic_eff_;//末端执行器转动惯量在质心处的表达


    };
    
    
}
