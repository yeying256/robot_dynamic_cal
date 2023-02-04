#ifndef MATH_H
#define MATH_H

#include <iostream>
#include <Eigen/Eigen>

namespace xj_dy_ns
{
    class math_wx
    {
    private:
        
    public:
        math_wx(/* args */);
        ~math_wx();

        //龙格库塔法迭代计算
        static Eigen::VectorXd Runge_Kutta_itr( Eigen::MatrixXd A,
                                                Eigen::VectorXd B,
                                                Eigen::VectorXd X,
                                                double dt
                                                );
    };
    
}



#endif