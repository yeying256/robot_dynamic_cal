#include "PID_controller.h"



namespace xj_dy_ns
{
    PID_controller::PID_controller(/* args */)
    {
        m_last_err = 0;
        m_i_err = 0;
        m_d_err = 0;
        max_output_ = 0.1;
    }
    
    PID_controller::~PID_controller()
    {
    }

    /**
     * @brief 最大输出默认0.1
     * 
     * @param p 
     * @param i 
     * @param d 
     */
    void PID_controller::init(double p, double i, double d)
    {
        m_p=p;
        m_i=i;
        m_d=d;
        m_last_err = 0;
        m_i_err = 0;
        m_d_err = 0;
        max_output_ = 0.1;
    }

    PID_controller::PID_controller(double p, double i, double d ,double max_output):m_p(p), m_i(i), m_d(d),max_output_(max_output)
    {
        m_last_err = 0;
        m_i_err = 0;
        m_d_err = 0;
    }


    /**
     * @brief 加入了maxoutput的初始化函数
     * 
     * @param p 
     * @param i 
     * @param d 
     * @param max_output 限制最大输出
     */
    void PID_controller::init(double p, double i, double d,double max_output)
    {
        m_p=p;
        m_i=i;
        m_d=d;
        m_last_err = 0;
        m_i_err = 0;
        m_d_err = 0;
        max_output_ = max_output;
    }


    /**
     * @brief pid计算
     * 
     * @param err 
     * @return double 输出pid的计算值 
     */
    double PID_controller::PID(double err)
    {
        m_i_err += m_last_err;
        m_d_err = err - m_last_err;
        double k = m_p*err + m_i*m_i_err + m_d*m_d_err;
        m_last_err = err;
        // if(fabs(k) < 0.0001) k = fabs(k) / k * 0.01;
        if(fabs(k) > max_output_ ) k = fabs(k) / k * max_output_;
        return k;
    }

    /**
     * @brief Construct a new pid controller::pid controller object 构造函数初始化pid的值
     * 
     * @param p 
     * @param i 
     * @param d 
     */
    PID_controller::PID_controller(double p, double i, double d):m_p(p), m_i(i), m_d(d)
    {
        m_last_err = 0;
        m_i_err = 0;
        m_d_err = 0;
    }
}