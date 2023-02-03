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
                                                Eigen::MatrixXd B,
                                                Eigen::VectorXd X,
                                                double dt
                                                );
    };
    
}