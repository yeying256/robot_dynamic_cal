#ifndef AGV_TURN_WHEEL_H 
#define AGV_TURN_WHEEL_H

#include <iostream>
#include <Eigen/Eigen>
#include <ros/ros.h>

namespace xj_dy_ns
{
    class agv_turn_wheel
    {
    private:
    double SinPhi_,CosPhi_,B_,L_,R_;
        

//           lr: 1
//   lf: 1

//   #轮子与车体固定点到轮子滚动中心的垂直距离
//   B: 1

//   #轮子的半径
//   R: 1
    public:
        agv_turn_wheel();
        agv_turn_wheel(double lr,double lf,double B,double R);
        bool init(double lr,double lf,double B,double R);


        void Steer_Wheel_Kinematics(Eigen::Ref<Eigen::Vector2d> Wheel_Drive_Omega,
	                        Eigen::Ref<Eigen::Vector2d> Wheel_Steer_Omega,

							Eigen::Vector2d Wheel_Steer_Alpha,
							Eigen::Vector3d xyw
							);

        ~agv_turn_wheel();
    };





}


#endif