#ifndef AGV_TURN_WHEEL_H 
#define AGV_TURN_WHEEL_H

#include <iostream>
#include <Eigen/Eigen>
#include <ros/ros.h>

#include <tf/transform_broadcaster.h>
#include <nav_msgs/Odometry.h>

namespace xj_dy_ns
{
    class agv_turn_wheel
    {
    private:
    double SinPhi_,CosPhi_,B_,L_,R_,Phi_;
    Eigen::Vector2d Wheel_Drive_Omega_,Wheel_Steer_Omega_,Wheel_Steer_Alpha_;
    Eigen::Vector3d odom_,Vxyw_now_;
    ros::Duration period_;
        

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
							Eigen::Vector3d Vxyw_cmd
							);
        
        //双舵轮对角布置运动学正解                    
        Eigen::Vector3d Steer_Wheel_forward(Eigen::Vector2d Wheel_Drive_Omega,
                                    Eigen::Vector2d Wheel_Steer_Omega,
                                    Eigen::Vector2d Wheel_Steer_Alpha);
        

        void update(Eigen::Vector2d Wheel_Drive_Omega,
                                    Eigen::Vector2d Wheel_Steer_Omega,
                                    Eigen::Vector2d Wheel_Steer_Alpha,
                                    const ros::Duration& period);
        void odom_updata();
        ~agv_turn_wheel();

        void tf_odom_trans(nav_msgs::Odometry &odom,
                            geometry_msgs::TransformStamped &odom_trans,
                            ros::Time time_now,
                            std::string frame_id,
                            std::string robot_frame_id);
        
    };





}


#endif