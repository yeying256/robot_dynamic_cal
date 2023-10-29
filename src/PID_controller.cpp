#include "PID_controller.h"



namespace xj_dy_ns
{
    PID_controller::PID_controller(/* args */)
    {
        m_last_err = 0;
        m_i_err = 0;
        m_d_err = 0;
    }
    
    PID_controller::~PID_controller()
    {
    }

    void PID_controller::init(double p, double i, double d)
    {
        m_p=p;
        m_i=i;
        m_d=d;
        m_last_err = 0;
        m_i_err = 0;
        m_d_err = 0;
    }

    double PID_controller::PID(double err)
    {
        m_i_err += m_last_err;
        m_d_err = err - m_last_err;
        double k = m_p*err + m_i*m_i_err + m_d*m_d_err;
        m_last_err = err;
        // if(fabs(k) < 0.0001) k = fabs(k) / k * 0.01;
        if(fabs(k) > 0.1 ) k = fabs(k) / k * 0.1;
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