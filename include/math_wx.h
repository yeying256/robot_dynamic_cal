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

    static void Cal_Axis_Angle_Jo_g(Eigen::MatrixXd R,
                            Eigen::MatrixXd dq,
                            Eigen::MatrixXd J,
                            Eigen::MatrixXd dJ,
                            Eigen::MatrixXd Rd,
                            Eigen::MatrixXd d_Rd,
                            Eigen::MatrixXd d2_Rd,
                            Eigen::Ref<Eigen::MatrixXd> out_Jo,
                            Eigen::Ref<Eigen::MatrixXd> out_d_Jo,
                            Eigen::Ref<Eigen::MatrixXd> out_g,
                            Eigen::Ref<Eigen::MatrixXd> out_d_g,
                            Eigen::Ref<Eigen::MatrixXd> fix_Jparam);

    static Eigen::MatrixXd R_to_so3_separate(Eigen::MatrixXd R);

    static double GCW_acos(double x);
    static Eigen::Vector3d GCW_cross(Eigen::Vector3d a, Eigen::Vector3d b);

    static Eigen::Matrix3d return_S(Eigen::Vector3d a);




    };
    
}



#endif