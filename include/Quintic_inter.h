#pragma once
#include <Eigen/Eigen>
#include <math.h>

namespace xj_dy_ns
{
    class Quintic_inter
    {
    private:
        Eigen::Matrix<double,6,Eigen::Dynamic> ai;//列的数量由插值的参数数量决定
        double deltaT_;
        int val_num_;//被插值的参数数量
        Eigen::Matrix<double,6,Eigen::Dynamic> state_;
        //顺序：每一行都是一个同一个待插值的变量的不同状态，从上到下分别是这个变量的 x0 dx0 ddx0 x1 dx1 ddx1 状态
        //每一列都是不同带插值的变量的6个状态。
        Eigen::Matrix<double,Eigen::Dynamic,6> A_;
        bool flag_cal = false;
    public:
        Quintic_inter();
        Quintic_inter(double deltaT);
        Quintic_inter(double deltaT,Eigen::Matrix<double,6,Eigen::Dynamic> state_);
        Quintic_inter(double deltaT,Eigen::VectorXd x0,Eigen::VectorXd x1);
        void ai_cal();//多项式参数计算
        void set_deltaT(double deltaT);
        void set_start_and_end_state(Eigen::Matrix<double,6,Eigen::Dynamic> state_);
        void set_start_and_end_state(Eigen::VectorXd x0,Eigen::VectorXd x1);
        Eigen::VectorXd q_now_cal(double t_now);
        ~Quintic_inter();
    };
    

    
}