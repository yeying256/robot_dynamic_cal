#ifndef PID_CONTROLLER_H
#define PID_CONTROLLER_H
#include "include_common.h"


namespace xj_dy_ns
{
    class PID_controller
    {
    private:
        /* data */
        private:
        double m_p;
        double m_i;
        double m_d;
        double m_last_err;
        double m_i_err;
        double m_d_err;
        double max_output_;
    public:
        PID_controller(/* args */);
        void init(double p, double i, double d);
        void init(double p, double i, double d,double max_output);
        ~PID_controller();
        PID_controller(double p, double i, double d);
        PID_controller(double p, double i, double d ,double max_output);

        double PID(double err);
    };
    

    
}



#endif