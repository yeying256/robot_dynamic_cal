#include <ros/ros.h>
#include <Eigen/Eigen>
#include <fstream>
#include <sstream>



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
    Eigen::Matrix<double,Eigen::Dynamic,1> G_;//补偿重力的关节力矩值
    std::vector<Eigen::Matrix<double,3,1>> F_T_EE;
    Eigen::Matrix<double,4,4> T_tool_;
    
    std::vector<Eigen::Matrix<double,4,4>> _0T_i;//世界坐标系下
    Eigen::Matrix<double,4,4> _0T_tool;//世界坐标系下的末端齐次变换矩阵

    std::vector<Eigen::Matrix<double,3,1>> v_;
    std::vector<Eigen::Matrix<double,3,1>> w_;
    std::vector<Eigen::Matrix<double,3,1>> a_;
    std::vector<Eigen::Matrix<double,3,1>> dw_;
    Eigen::Matrix<double,Eigen::Dynamic,1> dq_;//关节速度
    Eigen::Matrix<double,Eigen::Dynamic,1> dq_last_;//关节速度

    Eigen::Matrix<double,Eigen::Dynamic,1> ddq_;//关节加速度
    Eigen::Matrix<double,Eigen::Dynamic,1> ddq_last_;//上一次关节加速度

    std::vector<Eigen::Matrix<double,6,Eigen::Dynamic>> jacobi_ci_;
    Eigen::Matrix<double,6,Eigen::Dynamic> jacobi_;//雅可比矩阵
    std::vector<Eigen::Matrix<double,3,3> > Ic_;//转动惯量
    Eigen::Matrix<double,Eigen::Dynamic,1> tor_CpM_neton_;
    Eigen::Matrix<double,Eigen::Dynamic,1> tor_CpM_neton_last_;
    std::vector<double> motor_I_; //电机转子的惯量

    Eigen::MatrixXd M_q_;

    bool init_size_flag=false;
    int DOF_ = 0;
    public:
        Robot_dynamic();//需要调用resize_param函数来初始化大小
        Robot_dynamic(const std::string& url,const int& DOF);
        Robot_dynamic(int dof);//自带resize_param初始化参数大小
        ~Robot_dynamic();
        bool read_dynamic(const std::string& url,const int& DOF);
        void read_dynamic_type_Gong(Eigen::VectorXd param,int dof);
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
        Eigen::Matrix<double,3,3> get_R(int i);
        void tor_gravity_cal();//计算每个关节的重力项
        Eigen::Matrix<double,Eigen::Dynamic,1> get_G_();//获取补偿重力的关节力矩值
        double get_qi_now(int i);
        void jacobi_cal();
        Eigen::Matrix<double,3,3> get_0R(int i);
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

        void tor_M_C_neton_cal_();
        std::vector<Eigen::Matrix<double,6,1>> tor_M_C_neton_cal_(Eigen::Matrix<double,6,1>dq,
                                                                    Eigen::Matrix<double,6,1>ddq);

        Eigen::Matrix<double,Eigen::Dynamic,1> get_tor_CpM_neton_();
        Eigen::Matrix<double,Eigen::Dynamic,1> tor_filter(Eigen::Matrix<double,Eigen::Dynamic,1> tor_,
        double hz,
        double period,
        Eigen::Matrix<double,Eigen::Dynamic,1>* tor_last);
        

        Eigen::Matrix<double,Eigen::Dynamic,1>* get_tor_CpM_neton_last_ptr();
        Eigen::Matrix<double,4,4> get_0Ti(int i);
        void resize_param(int dof);
        void jacobe_ci_cal();
        bool get_joint_isrevolutor(int i);//判断是否是旋转关节
        Eigen::MatrixXd M_q_cal_Lagrange();//用拉格朗日方法求解惯性矩阵
        void updata_cal();

    };
    
    
}