#include <Quintic_inter.h>


namespace xj_dy_ns
{
    /**
     * @brief 默认构造函数，啥也不是
    */
    Quintic_inter::Quintic_inter()
    {;    }
    
    /**
     * @brief 默认析构函数，啥也不是
    */
    Quintic_inter::~Quintic_inter()
    {;    }

    /**
     * @brief 带deltaT的构造函数
     * @param deltaT 插值持续时间
    */
    Quintic_inter::Quintic_inter(double deltaT)
    {
        set_deltaT(deltaT);
    }

    /**
     * @brief 带参数和始末状态的的构造函数
     * @param deltaT 插值持续时间
     * @param state_ 始末状态 顺序：每一行都是一个同一个待插值的变量的不同状态，从上到下分别是这个变量的 x0 dx0 ddx0 x1 dx1 ddx1 状态每一列都是不同带插值的变量的6个状态。
    */
    Quintic_inter::Quintic_inter(double deltaT,Eigen::Matrix<double,6,Eigen::Dynamic> state)
    {
        this->deltaT_=deltaT;
        set_start_and_end_state(state);
    }

    /**
     * @brief 带参数和始末位置状态的的构造函数，始末速度和加速度默认为0
     * @param deltaT 插值持续时间
     * @param x0 不同参数的初始位置向量
     * @param x1 不同参数的结束位置向量
    */
    Quintic_inter::Quintic_inter(double deltaT,Eigen::VectorXd x0,Eigen::VectorXd x1)
    {
        this->deltaT_=deltaT;
        set_start_and_end_state(x0,x1);
    }

    /**
     * @brief 设置deltaT
     * @param deltaT 插值持续时间
    */
    void Quintic_inter::set_deltaT(double deltaT)
    {this->deltaT_=deltaT;}

    void Quintic_inter::set_start_and_end_state(Eigen::Matrix<double,6,Eigen::Dynamic> state)
    {
        this->state_=state;
        this->val_num_=state.cols();
        ai_cal();
    }

    void Quintic_inter::set_start_and_end_state(Eigen::VectorXd x0,Eigen::VectorXd x1)
    {
        val_num_= x0.size();
        this->state_.setZero(6,val_num_);
        this->state_.row(0) = x0.transpose();
        state_.row(3) = x1.transpose();
        ai_cal();
    }
    Eigen::VectorXd Quintic_inter::q_now_cal(double t_now)
    {
        Eigen::VectorXd x_now;
        Eigen::Matrix<double,6,1> pow_t_;

        x_now.setZero(this->val_num_);
        // std::cout<<"state_="<<state_<<std::endl;

        if (!this->flag_cal)
        {
            printf("\033[1;31;40m  五次插值有bug,初始化参数不全 \033[0m \n");
            return x_now;
        }
        pow_t_<<
        1,t_now,pow(t_now,2),pow(t_now,3),pow(t_now,4),pow(t_now,5);
        // for (int i = 0; i < this->val_num_; i++)
        // {
        //     x_now(i) = A_*;
        // }
        x_now = A_*pow_t_;

        if (t_now>=this->deltaT_)
        {
        x_now = last_x_;
        }
        else
        {
        x_now = A_*pow_t_;
        }
        last_x_=x_now;
        return x_now;
        
        
        
    }

    void Quintic_inter::ai_cal()//多项式参数计算
    {
        //6行 n列
        Eigen::Matrix<double,6,6> matrix_delta_t;
        matrix_delta_t<<
        1,0,0,0,0,0,
        0,1,0,0,0,0,
        0,0,2,0,0,0,
        1,deltaT_,pow(deltaT_,2),pow(deltaT_,3),pow(deltaT_,4),pow(deltaT_,5),
        0,1,2*deltaT_,3*pow(deltaT_,2),4*pow(deltaT_,3),5*pow(deltaT_,4),
        0,0,2,6*deltaT_,12*pow(deltaT_,2),20*pow(deltaT_,3);
        this->A_=(matrix_delta_t.inverse() * state_).transpose();

        flag_cal = true;
    }

}