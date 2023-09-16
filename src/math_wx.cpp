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


    /**
     * @brief 计算重要参数
     * 
     * @param R 
     * @param dq 
     * @param J 
     * @param dJ 
     * @param Rd 
     * @param d_Rd 
     * @param d2_Rd 
     * @param out_Jo 是3Xdof
     * @param out_d_Jo 是3Xdof
     * @param out_g 是3X1
     * @param out_d_g 是3X1
     * @param fix_Jparam 修正矩阵 旋转雅可比修正矩阵，左乘雅可比等于角度雅可比
     */
    void math_wx::Cal_Axis_Angle_Jo_g(Eigen::MatrixXd R,
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
                            Eigen::Ref<Eigen::MatrixXd> fix_Jparam)

    {
    // J 和 dJ 是 3×7 的， 只有角度部分的
    #pragma region 
        // Eigen::MatrixXd R(3, 3);
        // R = T.block<3, 3>(0, 0);

        Eigen::MatrixXd d_R(3, 3);
        d_R = return_S(J.block(0, 0, 3, 7) * dq) * R;

        // Eigen::MatrixXd Rd(3, 3);
        // Rd = Td.block<3, 3>(0, 0);

        // Eigen::MatrixXd d_Rd(3, 3);
        // d_Rd = d_Td.block<3, 3>(0, 0);

        // Eigen::MatrixXd d2_Rd(3, 3);
        // d2_Rd = d2_Td.block<3, 3>(0, 0);

        Eigen::MatrixXd so3(4, 1);
        so3 = R_to_so3_separate(R * Rd.transpose());

        Eigen::MatrixXd l(3, 1);
        l = so3.block<3, 1>(0, 0);//轴角的轴
        

        double theta = so3(3, 0);//轴角的角
        

        //Eigen::MatrixXd dJ(6, 7);
        // dJ = get_dJ(L, R_coordinate, robot_joint_parameters, q, dq, 7, L_end);
        //dJ = J_d_d2.D_Value * 1.0;

        Eigen::MatrixXd S_Jdq(3, 3);
        S_Jdq = return_S(J.block(0, 0, 3, 7) * dq);//计算角速度叉乘算子

        Eigen::MatrixXd Jo_2(3, 7);
        Eigen::MatrixXd GAMMA_2(3, 1);

        #pragma region theta == 0
        // if (theta == 0) { //这部分没验证过，因为根本进不来
        //   ROS_INFO("jinlaile");
        //   Eigen::Matrix<double, 3, 1> e1;
        //   Eigen::Matrix<double, 3, 1> e2;
        //   Eigen::Matrix<double, 3, 1> e3;
        //   e1 << 1, 0, 0;
        //   e2 << 0, 1, 0;
        //   e3 << 0, 0, 1;

        //   Jo_2 = -e1 * e3.transpose() * return_S(R * Rd.transpose() * e2) * J.block<3, 7>(3, 0) -
        //          e2 * e1.transpose() * return_S(R * Rd.transpose() * e3) * J.block<3, 7>(3, 0) -
        //          e3 * e2.transpose() * return_S(R * Rd.transpose() * e1) * J.block<3, 7>(3, 0);

        //   GAMMA_2 =
        //       -e1 * e3.transpose() * return_S(R * Rd.transpose() * e2) * dJ.block<3, 7>(3, 0) * dq -
        //       e2 * e1.transpose() * return_S(R * Rd.transpose() * e3) * dJ.block<3, 7>(3, 0) * dq -
        //       e3 * e2.transpose() * return_S(R * Rd.transpose() * e1) * dJ.block<3, 7>(3, 0) * dq +
        //       e1 * e3.transpose() * S_Jdq * S_Jdq * R * Rd.transpose() * e2 +
        //       e2 * e1.transpose() * S_Jdq * S_Jdq * R * Rd.transpose() * e3 +
        //       e3 * e2.transpose() * S_Jdq * S_Jdq * R * Rd.transpose() * e1 +
        //       e1 * e3.transpose() * (2 * d_R * d_Rd.transpose() + R * d2_Rd.transpose()) * e2 +
        //       e2 * e1.transpose() * (2 * d_R * d_Rd.transpose() + R * d2_Rd.transpose()) * e3 +
        //       e3 * e2.transpose() * (2 * d_R * d_Rd.transpose() + R * d2_Rd.transpose()) * e1 -
        //       e1 * e3.transpose() * (d_R * Rd.transpose() + R * d_Rd.transpose()) *
        //           (d_R * Rd.transpose() + R * d_Rd.transpose()) * e2 -
        //       e2 * e1.transpose() * (d_R * Rd.transpose() + R * d_Rd.transpose()) *
        //           (d_R * Rd.transpose() + R * d_Rd.transpose()) * e3 -
        //       e3 * e2.transpose() * (d_R * Rd.transpose() + R * d_Rd.transpose()) *
        //           (d_R * Rd.transpose() + R * d_Rd.transpose()) * e1;

        //   M_cart = (combine_matrix(J.block<3, 7>(0, 0), Jo_2) * B.inverse() * J.transpose()).inverse();

        //   tau_d =
        //       pseudoInverse(M_cart * combine_matrix(J.block<3, 7>(0, 0), Jo_2) * B.inverse()) *
        //       (fext - fd - K_cart * combine_matrix(T.block<3, 1>(0, 3) - Td.block<3, 1>(0, 3), so3) -
        //        B_cart * combine_matrix(J.block<3, 7>(0, 0) * dq - d_Td.block<3, 1>(0, 3),
        //                                return_S_vector(d_R * Rd.transpose() + R * d_Rd.transpose())) -   // 这里的return_S_vector 开始编错了 时候验证的 好吧一直没进来过， 就没验证过
        //        M_cart * combine_matrix(dJ.block<3, 7>(0, 0) * dq - d2_Td.block<3, 1>(0, 3), GAMMA_2) +   // 下面就没用到过  return_S_vector 也就没验证过
        //        M_cart * combine_matrix(J.block<3, 7>(0, 0), Jo_2) * B.inverse() * tau_C -
        //        M_cart * combine_matrix(J.block<3, 7>(0, 0), Jo_2) * B.inverse() * J.transpose() * fext);

        //   tau_d << 0, 0, 0, 0, 0, 0, 0;
        //   //return limit_tau_d(tau_d);
        // }
        #pragma endregion 

        Eigen::MatrixXd d_l(3, 1);
        d_l = (Rd * d_R.transpose() + d_Rd * R.transpose() + d_R * Rd.transpose() +
            R * d_Rd.transpose()) *
            so3.block<3, 1>(0, 0) / 2.0 / (1 - cos(so3(3, 0)));
        

        Eigen::Matrix<double, 3, 1> beta;
        beta << 1, 1, 1;
        for (int i = 0; i < 3; i++) {
        if (abs(l(i, 0)) > 0.4) {               //设置一个垂直的向量
            beta(i, 0) = (l(0, 0) + l(1, 0) + l(2, 0) - l(i, 0)) / (-l(i, 0));
            beta = beta * (1.0 / beta.norm());  //单位化处理
            break;                              //处理完一个向量，就垂直了
        }

        }
        Eigen::Matrix<double, 3, 1> alpha;
        alpha = GCW_cross(l, beta);

        Eigen::Matrix<double, 1, 1> d_theta_tmp;
        d_theta_tmp =
            alpha.transpose() * (d_R * Rd.transpose() + R * d_Rd.transpose()) * beta / cos(theta);
        double d_theta = d_theta_tmp(0, 0);

        fix_Jparam = ( theta/2/(1-cos(theta)) )*(Rd*R.transpose()*return_S(l) - return_S(R*Rd.transpose()*l))
                - 1/cos(theta)*l*alpha.transpose()*return_S(R*Rd.transpose()*beta);

        
        out_Jo = fix_Jparam*J;

        // out_Jo = ( theta/2/(1-cos(theta)) )*(Rd*R.transpose()*return_S(l)*J - return_S(R*Rd.transpose()*l)*J)
        //         - 1/cos(theta)*l*alpha.transpose()*return_S(R*Rd.transpose()*beta)*J;

        out_g = ( theta/2/(1-cos(theta)) )*(R*d_Rd.transpose()*l+d_Rd*R.transpose()*l) 
                + 1/cos(theta)*l*alpha.transpose()*R*d_Rd.transpose()*beta;
        //经过验证 out_Jo 和 out_g 是对的

        Eigen::Matrix<double,3,1> w;
        w = -(return_S(l)).completeOrthogonalDecomposition().pseudoInverse()*d_l;
        
        Eigen::Matrix<double,3,1> d_alpha;
        Eigen::Matrix<double,3,1> d_beta;

        d_beta = return_S(w)*beta;
        d_alpha = return_S(w)*alpha;

        out_d_Jo =  (d_theta/2/(1-cos(theta)) - theta*sin(theta)*d_theta/2/(1-cos(theta))/(1-cos(theta))   )
                    *( Rd*R.transpose()*return_S(l)*J - return_S(R*Rd.transpose()*l)*J  )
                +( theta/2/(1-cos(theta)) )
                    *(   d_Rd*R.transpose()*return_S(l)*J + Rd*d_R.transpose()*return_S(l)*J 
                    + Rd*R.transpose()*return_S(d_l)*J + Rd*R.transpose()*return_S(l)*dJ
                    - return_S(d_R*Rd.transpose()*l)*J - return_S(R*d_Rd.transpose()*l)*J
                    - return_S(R*Rd.transpose()*d_l)*J - return_S(R*Rd.transpose()*l)*dJ   )
                -1/cos(theta)*d_l*alpha.transpose()*return_S(R*Rd.transpose()*beta)*J
                -1/cos(theta)*l*d_alpha.transpose()*return_S(R*Rd.transpose()*beta)*J
                -(sin(theta)*d_theta/cos(theta)/cos(theta))*l*alpha.transpose()*return_S(R*Rd.transpose()*beta)*J
                -1/cos(theta)*l*alpha.transpose()*(
                    return_S(d_R*Rd.transpose()*beta)*J + return_S(R*d_Rd.transpose()*beta)*J
                    +return_S(R*Rd.transpose()*d_beta)*J + return_S(R*Rd.transpose()*beta)*dJ)  ;

        out_d_g = (d_theta/2/(1-cos(theta)) - theta*sin(theta)*d_theta/2/(1-cos(theta))/(1-cos(theta))   )
                    *(R*d_Rd.transpose()*l+d_Rd*R.transpose()*l)
                +( theta/2/(1-cos(theta)) )
                    *(d_R*d_Rd.transpose()*l+R*d2_Rd.transpose()*l+R*d_Rd.transpose()*d_l
                    +d2_Rd*R.transpose()*l+d_Rd*d_R.transpose()*l+d_Rd*R.transpose()*d_l)
                +1/cos(theta)*d_l*alpha.transpose()*R*d_Rd.transpose()*beta
                +1/cos(theta)*l*d_alpha.transpose()*R*d_Rd.transpose()*beta
                +(sin(theta)*d_theta/cos(theta)/cos(theta))*l*alpha.transpose()*R*d_Rd.transpose()*beta
                +1/cos(theta)*l*alpha.transpose()*(d_R*d_Rd.transpose()*beta+R*d2_Rd.transpose()*beta+R*d_Rd.transpose()*d_beta);

        // out_data.block<3,1>(0,0) = l;
        // out_data.block<3,1>(0,1) = d_l;
        // out_data.block<3,1>(0,2) = beta;
        // out_data.block<3,1>(0,3) = d_beta;
        // out_data.block<3,1>(0,4) = alpha;
        // out_data.block<3,1>(0,5) = d_alpha;
        #pragma endregion
    }

    Eigen::MatrixXd math_wx::R_to_so3_separate(Eigen::MatrixXd R) {
        Eigen::MatrixXd so3(4, 1);
        so3(3, 0) = GCW_acos((R(0, 0) + R(1, 1) + R(2, 2) - 1) /
                            2.0);  // acos 直接是 0 到 180度  就是说如果 R 穿过 I 的时候，
                                    // 前后角度都是正的， 但是 轴线方向会直接 变 反向

        Eigen::MatrixXd temp(3, 1);
        temp << R(2, 1) - R(1, 2), R(0, 2) - R(2, 0), R(1, 0) - R(0, 1);
        if (abs(sin(so3(3, 0))) < 1e-6) {
            so3.block<3, 1>(0, 0) << 1, 0, 0;
            so3(3, 0) = 0;
        } else if (so3(3, 0) < 0) {
            so3.block<3, 1>(0, 0) = -temp * (1.0 / temp.norm());//norm返回的是二范数
        } else {
            so3.block<3, 1>(0, 0) = temp * (1.0 / temp.norm());
        }
        // so3.block<3, 1>(0, 0) = so3.block<3, 1>(0, 0) / so3.block<3, 1>(0, 0).norm();
        return so3;
    }


    double math_wx::GCW_acos(double x) {
        if (x > 1) {
            return acos(1);
        } else if (x < -1) {
            return acos(-1);
        } else {
            return acos(x);
        }
    }


    Eigen::Vector3d math_wx::GCW_cross(Eigen::Vector3d a, Eigen::Vector3d b) {
        return return_S(a) * b;
    }


    /**
     * @brief 计算叉乘算子
     * 
     * @param a 
     * @return Eigen::Matrix3d 
     */
    Eigen::Matrix3d math_wx::return_S(Eigen::Vector3d a) {
        Eigen::Matrix3d b;
        b << 0, -a(2), a(1), a(2), 0, -a(0), -a(1), a(0), 0;
        return b;
    }

}