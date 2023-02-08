#include "math_wx.h"

namespace xj_dy_ns
{


    math_wx::math_wx(/* args */)
    {
    }
    
    math_wx::~math_wx()
    {
    }

    /**
     * @brief 龙格库塔法 \dot(X) = AX + B 输入第n次的X，返回第n+1的X
     * 
     * @param A 状态方程的A
     * @param B 状态方程的B
     * @param X 上一次的状态量
     * @param dt 步长
     * @return Eigen::VectorXd 返回微分方程当前计算的状态量
     */
    Eigen::VectorXd math_wx::Runge_Kutta_itr(Eigen::MatrixXd A,
                                            Eigen::VectorXd B,
                                            Eigen::VectorXd X,
                                            double dt
                                            )
    {
        Eigen::VectorXd K1,K2,K3,K4,X_now = Eigen::VectorXd::Zero(X.size());
        K1=A*X+B;
        K2=A*(X+0.5*dt*K1)+B;   
        K3=A*(X+0.5*dt*K2)+B;   
        K4=A*(X+dt*K3)+B;       
        X_now=X+dt/6*(K1+2*K2+2*K3+K4);
        return X_now;
    }
}