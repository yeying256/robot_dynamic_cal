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
    Eigen::Vector2d Wheel_Drive_Omega_,Wheel_Steer_Omega_,Wheel_Steer_Alpha_,wheel_rad_last_,Alpha_last_;
    Eigen::Vector3d odom_,Vxyw_now_,d_xyrz_; //d_xyrz_微分量
    ros::Duration period_;
    int num_wheel_;


    /*
    Drive_Omega驱动轮角速度
    Steer_Omega舵的角速度
    Alpha轮子与x轴正方向夹角
    vx，vy，omega 给定的指令速度
    B 偏置距离
    L 偏心距离
    R 滚动半径
    Phi 某一个轮子与车体坐标系原点连线相对于x轴正向的夹角
    */
   //小写的都是向量
    Eigen::VectorXd b_,l_,r_,phi_,alpha_init_;//这是一个向量

    

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
        bool init(Eigen::VectorXd phi_,Eigen::VectorXd B,Eigen::VectorXd R,Eigen::VectorXd L,Eigen::VectorXd alpha_init);


        void Inverse_Kinematics_new(Eigen::Ref<Eigen::VectorXd> Drive,
	                        Eigen::Ref<Eigen::VectorXd> Turn,
							Eigen::VectorXd Alpha,
							Eigen::Vector3d Vxyw_cmd
							);
        void Inverse_Kinematics_new2(Eigen::Ref<Eigen::VectorXd> Drive,
	                        Eigen::Ref<Eigen::VectorXd> Turn,
							Eigen::VectorXd Alpha,
							Eigen::Vector3d Vxyw_cmd
							);


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

        void update(Eigen::Vector2d wheel_rad_now,
                                    Eigen::Vector2d Alpha,
                                    const ros::Duration& period);
        // 使用反馈的速度来计算里程计
        void odom_updata();
        // 使用反馈的位置微小位移来计算里程计
        void d_odom_updata();

        ~agv_turn_wheel();

        void tf_odom_trans(nav_msgs::Odometry &odom,
                            geometry_msgs::TransformStamped &odom_trans,
                            ros::Time time_now,
                            std::string frame_id,
                            std::string robot_frame_id);
        Eigen::Vector3d Forward_Kinematics2(Eigen::Vector2d Drive_Omega,
                                    Eigen::Vector2d turn_Omega,
                                    Eigen::Vector2d Alpha);//只适用于两个主动舵轮的正运动学，如果多于这个数，必须要重算。
        
    };





}


#endif