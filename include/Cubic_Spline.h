#ifndef CUBIC_SPLINE_H
#define CUBIC_SPLINE_H
#include "include_common.h"

namespace xj_dy_ns
{
    class Cubic_Spline
    {
    private:
        // 每个时间段 比如0.1 0.1 0.1…… 0.1
        Eigen::VectorXd h_;

        // 数据的自由度，有几组数据需要差值
        int Dof_,num_point_,num_time_;
        Eigen::Matrix<double,2,Eigen::Dynamic> bound_;
        
        // 插值点 
        Eigen::MatrixXd point_;

        // 边界加速度 如果是单纯离散点就是代表弯矩
        Eigen::MatrixXd M_;
        // 维度和M一致，在方程组中是AX=b 中的b
        Eigen::MatrixXd d_;

        // 中间参数 mui = h_i/( h_i + h_{i+1} )这个的数量比时间段少1
        Eigen::VectorXd mu_;
        // 中间参数 lambda_i = 1-mu_i 比时间段少一
        Eigen::VectorXd lambda_;





    public:
        // 边界条件类型
        enum BOUND_type{
            ACC = 0,
            VELOCITY,
            CYCLE
        }bound_type_;

        // 构造函数
        Cubic_Spline(/* args */);
        ~Cubic_Spline();

        /**
         * @brief 初始化三次样条插值
         * 
         * @param point 插值点，每一行代表一组数据（x0,x1……xn），每一列代表这组数据的一个点（xk，yk……）
         * @param h 每一段时间间隔
         * @param bound_type 边界条件类型，0是加速度 1 是速度 2是周期
         * @param bound 边界条件，每一行是不同组数据的边界条件，第一列是初始的边界条件，第二列是结束的边界条件
         */
        void init(Eigen::MatrixXd point,
                Eigen::VectorXd h,
                BOUND_type bound_type,
                Eigen::Matrix<double,2,Eigen::Dynamic> bound);

        /**
         * @brief 返回数值
         * 
         * @param t 输入时间
         * @return Eigen::ArrayXd 
         */
        Eigen::ArrayXd cal(double t);
        


    };
    

    
    ;
}


#endif