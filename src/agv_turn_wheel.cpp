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

// #define OLD_VISION


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
        L_ =l;
        this->SinPhi_=lf/l;
        this->CosPhi_=lr/l;
        this->B_=B;
        this->R_=R;
        Phi_=atan(lf/lr);
        return true;
    }

    /**
     * @brief 
     * 		↑y
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
     * 所有的参数都是与x轴的夹角，逆时针为正，赋值小写的都是向量
     * @param phi
     * @param B 
     * @param R 
     * @param L 
     * @param alpha_init 轮子初始状态与x轴的夹角
     * @return true 
     * @return false 
     */
    bool agv_turn_wheel::init(Eigen::VectorXd phi,
    Eigen::VectorXd B,
    Eigen::VectorXd R,
    Eigen::VectorXd L,
    Eigen::VectorXd alpha_init)
    {
        
        this->num_wheel_=phi.size();//舵轮的个数
        this->phi_=phi;
        this->b_=B;
        this->l_=L;
        this->r_=R;
        this->alpha_init_=alpha_init;
        return true;
    }

   /**
    * @brief 
    * 
    * @param Drive 行走电机的指令速度
    * @param Turn 转向轮子的指令速度
    * @param Alpha 当前相对初始位置的转向轮子的角位置
    * @param Vxyw_cmd 指令速度
    */
    void agv_turn_wheel::Inverse_Kinematics_new(Eigen::Ref<Eigen::VectorXd> Drive,
	                        Eigen::Ref<Eigen::VectorXd> Turn,
							Eigen::VectorXd Alpha,
							Eigen::Vector3d Vxyw_cmd
							)
    {
        if (Alpha.size()!=alpha_init_.size())
        {
            std::cout<<"逆向运动学中Alpha和alpha_init_长度不相等"<<std::endl;
            return;
        }
        
        Eigen::VectorXd alpha_now = Alpha + alpha_init_;
        // std::cout<<"检测是否运行到此处1"<<std::endl;

        for (int i = 0; i < this->num_wheel_; i++)
        {
            Eigen::Vector2d Psi;
            Psi<<cos(alpha_now[i]),sin(alpha_now[i]);//向量，轮子当前所处的角度
            double Steer_Omega_Vx=-sin(alpha_now[i])*Vxyw_cmd[0];//vx
            double Steer_Omega_Vy= cos(alpha_now[i])*Vxyw_cmd[1];//vy
            // std::cout<<"检测是否运行到此处2"<<std::endl;
        
            if(Psi.dot(Vxyw_cmd.topRows(2))<0)//输入车体中心速度的向量，内积＜0则反转轮子自旋的角度。
            {
            Steer_Omega_Vx=-Steer_Omega_Vx;
            Steer_Omega_Vy=-Steer_Omega_Vy;
            }
            // std::cout<<"检测是否运行到此处3"<<std::endl;

            // Drive(0) = ((-sin(Alpha[0])*Vxyw_cmd(0)) + cos(Alpha[0])*Vxyw_cmd(1) + Vxyw_cmd(2)*    L*cos(Phi-Alpha[0]))  /R;
            Drive[i] = ( cos(alpha_now[i])*Vxyw_cmd[0]  + sin(alpha_now[i])*Vxyw_cmd[1]  + Vxyw_cmd[2] *l_[i]*sin(alpha_now[i]-phi_[i]))/r_[i];

            // Drive[i] = ( cos(alpha_now[i])*Vxyw_cmd[0]  + sin(alpha_now[i])*Vxyw_cmd[1]  + Vxyw_cmd[2] *l_[i]*sin(phi_[i]-alpha_now[i]))/r_[i];
            Turn[i] = ( Steer_Omega_Vx + Steer_Omega_Vy + abs(Vxyw_cmd[2])* ( b_[i]+l_[i]*cos(alpha_now[i]-phi_[i]) ) )/b_[i];
            // std::cout<<"检测是否运行到此处4"<<std::endl;
        // Turn(0)  = ( Steer_Omega_1_Vx+Steer_Omega_1_Vy + abs(Vxyw_cmd(2))*(B-(L*sin(Alpha[0]-Phi))))/B;
            // 

        }
        
    }

    /**
     * @brief 双舵轮机器人移动底盘逆向运动学
     * 
     * @param Wheel_Drive_Omega 传出参数，行走电机角速度
     * @param Wheel_Steer_Omega 传出参数，转向电机角速度
     * @param Wheel_Steer_Alpha 传入参数，转向电机的当前角度
     * @param Vxyw_cmd 传入参数，机器人指令速度
     */
    void agv_turn_wheel::Steer_Wheel_Kinematics(Eigen::Ref<Eigen::Vector2d> Drive,
	                        Eigen::Ref<Eigen::Vector2d> Turn,

							Eigen::Vector2d Alpha,
							Eigen::Vector3d Vxyw_cmd
							)
    {

        double B=this->B_;
        double L=this->L_;
        double Phi = this->Phi_;
        double R=R_;
        

        #ifndef OLD_VISION

        

        Drive(0) = ((-sin(Alpha[0])*Vxyw_cmd(0)) + cos(Alpha[0])*Vxyw_cmd(1) + Vxyw_cmd(2)*    L*cos(Phi-Alpha[0]))  /R;
        Drive(1) = ((-sin(Alpha[1])*Vxyw_cmd(0)) + cos(Alpha[1])*Vxyw_cmd(1) - Vxyw_cmd(2)*    L*cos(Phi-Alpha[1]))  /R;


        double Steer_Omega_1_Vx,Steer_Omega_1_Vy,Steer_Omega_2_Vx,Steer_Omega_2_Vy,Steer_Omega_1_Omega,Steer_Omega_2_Omega;
        Steer_Omega_1_Vx=(-cos(Alpha[0])*Vxyw_cmd(0));
        Steer_Omega_1_Vy=- sin(Alpha[0])*Vxyw_cmd(1);
        Steer_Omega_2_Vx=(-cos(Alpha[1])*Vxyw_cmd(0));
        Steer_Omega_2_Vy=- sin(Alpha[1])*Vxyw_cmd(1);

        std::vector<Eigen::Vector2d> vector_alpha;
        vector_alpha.resize(2);
        vector_alpha[0]<< -sin(Alpha[0]),cos(Alpha[0]);
        vector_alpha[1]<< -sin(Alpha[1]),cos(Alpha[1]);




        // if(vector_alpha[0].dot(Vxyw_cmd.topRows(2))<0)
        // {
        //     Steer_Omega_1_Vx=-Steer_Omega_1_Vx;
        //     Steer_Omega_1_Vy=-Steer_Omega_1_Vy;

        // }

        // if (vector_alpha[1].dot(Vxyw_cmd.topRows(2))<0)
        // {
        //     Steer_Omega_2_Vx=-Steer_Omega_2_Vx;
        //     Steer_Omega_2_Vy=-Steer_Omega_2_Vy;
        // }

        Turn(0)  = ( Steer_Omega_1_Vx+Steer_Omega_1_Vy + Vxyw_cmd(2)*(B-(L*sin(Alpha[0]-Phi))))/B;
        Turn(1)  = ( Steer_Omega_2_Vx+Steer_Omega_2_Vy + Vxyw_cmd(2)*(B+(L*sin(Alpha[1]-Phi))))/B;

        

        #else




        // std::cout<<"\033[1;36;40m "<<"Wheel_Steer_Alpha="<<Wheel_Steer_Alpha<<"\033[0m "<<std::endl;
        // std::cout<<"\033[1;36;40m "<<"Wheel_Steer_Alpha[1]="<<Wheel_Steer_Alpha[1]<<"\033[0m "<<std::endl;
        // std::cout<<"\033[1;36;40m "<<"Wheel_Steer_Alpha(1)="<<Wheel_Steer_Alpha(1)<<"\033[0m "<<std::endl;


        // std::cout<<"\033[1;36;40m "<<"xyw="<<xyw<<"\033[0m "<<std::endl;
        // std::cout<<"\033[1;36;40m "<<"Wheel_Drive_Omega="<<Wheel_Drive_Omega<<"\033[0m "<<std::endl;
        // std::cout<<"\033[1;36;40m "<<"Wheel_Steer_Omega="<<Wheel_Steer_Omega<<"\033[0m "<<std::endl;


        // std::cout<<"\033[1;36;40m "<<"(float)R_="<<(float)R_<<"\033[0m "<<std::endl;
        // std::cout<<"\033[1;36;40m "<<"xyw(2) *L_ * (SinPhi_*sinf(Wheel_Steer_Alpha[0]) + CosPhi_*cosf(Wheel_Steer_Alpha[0])) ="<<xyw(2) *L_ * (SinPhi_*sinf(Wheel_Steer_Alpha[0]) + CosPhi_*cosf(Wheel_Steer_Alpha[0]))<<"\033[0m "<<std::endl;
        // std::cout<<"\033[1;36;40m "<<"L_="<<L_<<"\033[0m "<<std::endl;




        double Drive_Omega_1 = ((-sin(Alpha_1)*Vx) + cos(Alpha_1)*Vy + Omega*    L*cos(Phi-Alpha_1))  /R;
        double Steer_Omega_1 = ((-cos(Alpha_1)*Vx) - sin(Alpha_1)*Vy + Omega*(B-(L*sin(Alpha_1-Phi))))/B;
        double Drive_Omega_2 = ((-sin(Alpha_2)*Vx) + cos(Alpha_2)*Vy - Omega*    L*cos(Phi-Alpha_2))  /R;
        double Steer_Omega_2 = ((-cos(Alpha_2)*Vx) - sin(Alpha_2)*Vy + Omega*(B+(L*sin(Alpha_2-Phi))))/B;


        Wheel_Drive_Omega(0) = Drive_Omega_1;
        Wheel_Drive_Omega(1) = Drive_Omega_2;


        Wheel_Steer_Omega[0] = Steer_Omega_1;
        Wheel_Steer_Omega[1] = Steer_Omega_2;

        #endif

    }




    /**
     * @brief 正向运动学
     * 
     * @param Wheel_Drive_Omega 两个驱动轮的角速度 
     * @param Wheel_Steer_Omega 两个转向电机的角速度
     * @param Wheel_Steer_Alpha 两个转向电机的角度
     * @return Eigen::Vector3d 整体机器人的vx vy w
     */
    Eigen::Vector3d agv_turn_wheel::Steer_Wheel_forward(Eigen::Vector2d Wheel_Drive_Omega,
                                    Eigen::Vector2d Wheel_Steer_Omega,
                                    Eigen::Vector2d Wheel_Steer_Alpha)
    {   
        double Vx,Vy,Omega;
        double R=R_;
        double B=B_;
        double Alpha_1 = Wheel_Steer_Alpha[0];
        double Alpha_2 = Wheel_Steer_Alpha[1];
        double L=L_;
        double Phi = Phi_;
        double Drive_Omega_1 = Wheel_Drive_Omega[0];
        double Drive_Omega_2 = Wheel_Drive_Omega[1];
        double Steer_Omega_1 = Wheel_Steer_Omega[0];
        double Steer_Omega_2 = Wheel_Steer_Omega[1];





        Vx    = ((R*(B*cos(Alpha_1)+L*sin(Phi))*(Drive_Omega_1*sin(Alpha_1)+Drive_Omega_2*sin(Alpha_2)))-(R*Drive_Omega_1*sin(Alpha_1)*B*(cos(Alpha_1)+cos(Alpha_2)))+(B*(B*cos(Alpha_1)+L*sin(Phi))*(Steer_Omega_1*cos(Alpha_1)+Steer_Omega_2*cos(Alpha_2)))-(Steer_Omega_1*B*cos(Alpha_1)*B*(cos(Alpha_1)+cos(Alpha_2))))/(B*(cos(Alpha_2)-cos(Alpha_1))-2*L*sin(Phi));
    
        Omega = (2*Vx + R*(Drive_Omega_1*sin(Alpha_1)+Drive_Omega_2*sin(Alpha_2)) + B*(Steer_Omega_1*cos(Alpha_1)+Steer_Omega_2*cos(Alpha_2)))/(B*(cos(Alpha_1)+cos(Alpha_2)));
	
	    Vy    =  R*Drive_Omega_1*cos(Alpha_1)+Omega*(B*sin(Alpha_1)-L*cos(Phi))-Steer_Omega_1*B*sin(Alpha_1);

        Eigen::Vector3d vxyw;
        vxyw<< Vx,Vy,Omega;
        return vxyw;
    }

    /**
     * @brief 更新当前参数
     * 
     * @param Wheel_Drive_Omega 
     * @param Wheel_Steer_Omega 
     * @param Wheel_Steer_Alpha 
     */
    void agv_turn_wheel::update(Eigen::Vector2d Wheel_Drive_Omega,
                            Eigen::Vector2d Wheel_Steer_Omega,
                            Eigen::Vector2d Wheel_Steer_Alpha,
                            const ros::Duration& period)
    {
        Wheel_Drive_Omega_=Wheel_Drive_Omega;
        Wheel_Steer_Omega_=Wheel_Steer_Omega;
        Wheel_Steer_Alpha_=Wheel_Steer_Alpha;
        
        Vxyw_now_=Steer_Wheel_forward(Wheel_Drive_Omega_,Wheel_Steer_Omega_,Wheel_Steer_Alpha_);
        this->period_=period;
        odom_updata();

    }

    /**
     * @brief 使用内部参数更新里程计
     * 
     */
    void agv_turn_wheel::odom_updata()
    {
        Eigen::Matrix3d A;
        A<<cos(odom_(2)),-sin(odom_(2)),0,
            sin(odom_(2)),cos(odom_(2)),0,
            0,0,1;
        odom_=odom_+A*Vxyw_now_*period_.toSec();
        
    }

    /**
     * @brief 
     * 
     * @param odom 要发布的odom话题
     * @param odom_trans 要发布的tf变换
     * @param time_now 当前的时间
     * @param frame_id 里程计frame名称
     * @param robot_frame_id 机器人变换后的frame名称
     */
    void agv_turn_wheel::tf_odom_trans(nav_msgs::Odometry &odom,
                    geometry_msgs::TransformStamped &odom_trans,
                    ros::Time time_now,
                    std::string frame_id,
                    std::string robot_frame_id)
    {

        //tf信息加载

        //当前系统时间戳
        odom_trans.header.stamp = time_now;
        //里程计坐标系
        odom_trans.header.frame_id = frame_id;
        //机器人坐标系
        odom_trans.child_frame_id = robot_frame_id;
        //当前tf变换的位置部分
        odom_trans.transform.translation.x = this->odom_(0);
        odom_trans.transform.translation.y = this->odom_(1);
        odom_trans.transform.translation.z = 0;
        //当前tf变换的旋转部分
        geometry_msgs::Quaternion odom_quat = tf::createQuaternionMsgFromYaw(this->odom_(2));//转化成四元数
        odom_trans.transform.rotation = odom_quat;

        //odom信息更新
        //当前系统时间戳
        odom.header.stamp = time_now;
        //里程计坐标系
        odom.header.frame_id = frame_id;
        //机器人坐标系
        odom.child_frame_id = robot_frame_id;
        //设置位置部分
        odom.pose.pose.position.x=odom_(0);
        odom.pose.pose.position.y=odom_(1);
        odom.pose.pose.position.z=0.0;
        //设置姿态部分
        odom.pose.pose.orientation = odom_quat;
        //设置线速度
        odom.twist.twist.linear.x = Vxyw_now_(0);
        odom.twist.twist.linear.y = Vxyw_now_(1);
        odom.twist.twist.linear.z = 0.0;
        //设置角速度
        odom.twist.twist.angular.x = 0.0;
        odom.twist.twist.angular.y = 0.0;
        odom.twist.twist.angular.z = Vxyw_now_(2);
    }



}
