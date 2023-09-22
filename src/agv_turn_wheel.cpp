/**
 * @file agv_turn_wheel.cpp
 * @author wangxiao (2569209027@qq.com)
 * @brief 
 * @version 0.1
 * @date 2023-09-20
 * 
 * @copyright Copyright (c) 2023
 * 
 */


#include "agv_turn_wheel.h"


namespace xj_dy_ns
{
    agv_turn_wheel::agv_turn_wheel(/* args */)
    {
    }
    
    agv_turn_wheel::~agv_turn_wheel()
    {
    }

    agv_turn_wheel::agv_turn_wheel(double lr,double lf,double B,double R)
    {
        init(lr, lf, B, R);
    }

    /**
     * @brief 
     *  车体坐标系在车的正中心，左下角轮子为2号主动轮，右上角为1号主动轮。且两个轮子的偏置都朝右，与X轴平行。
		↑y
		丨   1
		丨  /|
		丨 / | lf
		丨/φ |
--------＋-------->x
	   /丨  lr
	L /	丨
	 /	丨
	/	丨
  2	    丨

     * @param lr Lright x轴方向轮子偏转中心到车体中心的距离
     * @param lf Lfrount y轴方向的轮子偏转中心到车体的距离
     * @param B 偏置距离，轮子滚动轴心至轮系与车体固定点的水平距离
     * @param R 轮子滚动半径
     * @return true 
     * @return false 
     */
    bool agv_turn_wheel::init(double lr,double lf,double B,double R)
    {
        double l =sqrt(pow(lr,2)+pow(lr,2));
        this->SinPhi_=lf/l;
        this->CosPhi_=lr/l;
        this->B_=B;
        this->R_=R;
        return true;
    }

    void agv_turn_wheel::Steer_Wheel_Kinematics(Eigen::Ref<Eigen::Vector2d> Wheel_Drive_Omega,
	                        Eigen::Ref<Eigen::Vector2d> Wheel_Steer_Omega,

							Eigen::Vector2d Wheel_Steer_Alpha,
							Eigen::Vector3d xyw
							)
    {
        Wheel_Drive_Omega(0)=(-sinf(Wheel_Steer_Alpha[0])*xyw(0) + cosf(Wheel_Steer_Alpha[0])*xyw(1) + xyw(2) *L_ * (SinPhi_*sinf(Wheel_Steer_Alpha[0]) + CosPhi_*cosf(Wheel_Steer_Alpha[0])))/(float)R_;
        Wheel_Steer_Omega[0]=(-cosf(Wheel_Steer_Alpha[0])*xyw(0) - sinf(Wheel_Steer_Alpha[0])*xyw(1) + xyw(2) * (B_ - (L_ * (CosPhi_*sinf(Wheel_Steer_Alpha[0]) - SinPhi_*cosf(Wheel_Steer_Alpha[0])))))/(float)B_;
        Wheel_Drive_Omega[1]=(-sinf(Wheel_Steer_Alpha[1])*xyw(0) + cosf(Wheel_Steer_Alpha[1])*xyw(1) - xyw(2) *       L_ * (SinPhi_*sinf(Wheel_Steer_Alpha[1]) + CosPhi_*cosf(Wheel_Steer_Alpha[1])))/(float)R_;
        Wheel_Steer_Omega[1]=(-cosf(Wheel_Steer_Alpha[1])*xyw(0) - sinf(Wheel_Steer_Alpha[1])*xyw(1) + xyw(2) * (B_ - (L_ * (SinPhi_*cosf(Wheel_Steer_Alpha[1]) - CosPhi_*sinf(Wheel_Steer_Alpha[1])))))/(float)R_;
    }



}
