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

#define OLD_VISION


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
     * @brief 双舵轮机器人移动底盘逆向运动学
     * 
     * @param Wheel_Drive_Omega 传出参数，行走电机角速度
     * @param Wheel_Steer_Omega 传出参数，转向电机角速度
     * @param Wheel_Steer_Alpha 传入参数，转向电机的当前角度
     * @param xyw 传入参数，机器人指令速度
     */
    void agv_turn_wheel::Steer_Wheel_Kinematics(Eigen::Ref<Eigen::Vector2d> Wheel_Drive_Omega,
	                        Eigen::Ref<Eigen::Vector2d> Wheel_Steer_Omega,

							Eigen::Vector2d Wheel_Steer_Alpha,
							Eigen::Vector3d xyw
							)
    {

        double Vx = xyw(0);
        double Vy = xyw(1);
        double Omega = xyw(2);
        double B=this->B_;
        double L=this->L_;
        double Phi = this->Phi_;
        double R=R_;


        


        double Alpha_1 = Wheel_Steer_Alpha[0];
        double Alpha_2 = Wheel_Steer_Alpha[1];
        

        #ifndef OLD_VISION

        while(Alpha_1>2*M_PI)
        {
            Alpha_1-=2*M_PI;
        }

        while(Alpha_1<-2*M_PI)
        {
            Alpha_1+=2*M_PI;
        }

        while(Alpha_2>2*M_PI)
        {
            Alpha_2-=2*M_PI;
        }

        while(Alpha_2<-2*M_PI)
        {
            Alpha_2+=2*M_PI;
        }


        
        Wheel_Drive_Omega(0)=(-sinf(Wheel_Steer_Alpha[0])*xyw(0) 
        + cosf(Wheel_Steer_Alpha[0])*xyw(1) 
        + xyw(2) *L_ * (SinPhi_*sinf(Wheel_Steer_Alpha[0]) + CosPhi_*cosf(Wheel_Steer_Alpha[0])))/(float)R_;
        // Wheel_Steer_Omega[0]=(-cosf(Wheel_Steer_Alpha[0])*xyw(0) 
        // - sinf(Wheel_Steer_Alpha[0])*xyw(1) 
        // + xyw(2) * (B_ - (L_ * (CosPhi_*sinf(Wheel_Steer_Alpha[0]) - SinPhi_*cosf(Wheel_Steer_Alpha[0])))))/(float)B_;
        Wheel_Drive_Omega[1]=(-sinf(Wheel_Steer_Alpha[1])*xyw(0) 
        + cosf(Wheel_Steer_Alpha[1])*xyw(1) 
        - xyw(2) * L_ * (SinPhi_*sinf(Wheel_Steer_Alpha[1]) + CosPhi_*cosf(Wheel_Steer_Alpha[1])))/(float)R_;
        // Wheel_Steer_Omega[1]=(-cosf(Wheel_Steer_Alpha[1])*xyw(0) 
        // - sinf(Wheel_Steer_Alpha[1])*xyw(1) 
        // + xyw(2) * (B_ - (L_ * (SinPhi_*cosf(Wheel_Steer_Alpha[1]) - CosPhi_*sinf(Wheel_Steer_Alpha[1])))))/(float)B_;

        double Steer_Omega_1_Vx,Steer_Omega_1_Vy,Steer_Omega_2_Vx,Steer_Omega_2_Vy,Steer_Omega_1_Omega,Steer_Omega_2_Omega;
        Steer_Omega_1_Vx=(-cos(Alpha_1)*Vx);
        Steer_Omega_1_Vy=(-sin(Alpha_1)*Vy);
        Steer_Omega_2_Vx=(-cos(Alpha_2)*Vx);
        Steer_Omega_2_Vy=(-sin(Alpha_2)*Vy);

        Steer_Omega_1_Omega = fabs(Omega)*(B-(L*sin(Alpha_1-Phi)));
        Steer_Omega_2_Omega = fabs(Omega)*(B+(L*sin(Alpha_2-Phi)));



        if((xyw(0) >0)&&(0<Alpha_1)&&(Alpha_1<3.1415926535))
            Steer_Omega_1_Vx=-Steer_Omega_1_Vx;
        else if((Vx<0)&&(3.1415926535<Alpha_1))
            Steer_Omega_1_Vx=-Steer_Omega_1_Vx;
        if((Vy>0)&&(1.5707963267948966192313216916398<Alpha_1)&&(Alpha_1<4.7123889803846898576939650749193))
            Steer_Omega_1_Vy=-Steer_Omega_1_Vy;
        else if((Vy<0)&&((Alpha_1<1.5707963267948966192313216916398)||(Alpha_1>4.7123889803846898576939650749193)))
            Steer_Omega_1_Vy=-Steer_Omega_1_Vy;

        if((Vx>0)&&(0<Alpha_2)&&(Alpha_2<3.1415926535))
            Steer_Omega_2_Vx=-Steer_Omega_2_Vx;
        else if((Vx<0)&&(3.1415926535<Alpha_2))
            Steer_Omega_2_Vx=-Steer_Omega_2_Vx;
        if((Vy>0)&&(1.5707963267948966192313216916398<Alpha_2)&&(Alpha_2<4.7123889803846898576939650749193))
            Steer_Omega_2_Vy=-Steer_Omega_2_Vy;
        else if((Vy<0)&&((Alpha_2<1.5707963267948966192313216916398)||(Alpha_2>4.7123889803846898576939650749193)))
            Steer_Omega_2_Vy=-Steer_Omega_2_Vy;


        



        
        
        double Steer_Omega_1 = (Steer_Omega_1_Vx + Steer_Omega_1_Vy + Steer_Omega_1_Omega)/B;
        double Steer_Omega_2 = (Steer_Omega_2_Vx + Steer_Omega_2_Vy + Steer_Omega_2_Omega)/B;


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



}
