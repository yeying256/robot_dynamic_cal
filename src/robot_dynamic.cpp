#include "robot_dynamic.h"



namespace xj_dy_ns
{
    /**
     * @brief 构造函数啥也不是，需要手动重置参数大小，调用resize_param函数
    */
    Robot_dynamic::Robot_dynamic(/* args */)
    {
        ;
    }

    /**
     * @brief 构造函数，根据dh表文件地址初始化动力学参数和运动学参数
     * @param dof 自由度数
    */
    Robot_dynamic::Robot_dynamic(const std::string& url,const int& DOF)
    {
        read_dynamic(url,DOF);
        this->resize_param(DOF);
    }

    /**
     * @brief 构造函数，仅初始化参数大小
     * @param dof 自由度数
    */
    Robot_dynamic::Robot_dynamic(int dof)
    {
        this->resize_param(dof);
    }
    
    Robot_dynamic::~Robot_dynamic()
    {
        ;
    }

    /**
     * @brief 初始化参数大小
     * @param param 西西里亚诺那本书定义的参数
     * @param dof 自由度数
    */
    void Robot_dynamic::resize_param(int dof)
    {
        this->DOF_=dof;
        this->a_.resize(dof);
        this->ddq_.setZero(dof);
        this->ddq_last_.setZero(dof);
        m_.resize(dof);
        Pc.resize(dof);
        Ic_.resize(dof);
        motor_I_.resize(dof);
        this->dq_.setZero(dof);
        this->dq_last_.setZero(dof);
        this->G_.setZero(dof);
        this->v_.resize(dof);
        this->w_.resize(dof);
        this->T_.resize(dof);
        this->jacobi_ci_.resize(dof);
        for (int i = 0; i < dof; i++)
        {
            jacobi_ci_[i].setZero(6,dof);
        }
        tor_CpM_neton_= Eigen::VectorXd::Zero(dof);
        DH_table_red.setZero(dof,7);
        DH_table.setZero(dof,7);
        _0T_i.resize(dof);
        M_q_.setZero(dof,dof);
        q_now.resize(dof);
        this->tor_CpM_neton_.resize(dof);
        this->tor_CpM_neton_last_.resize(dof);


        // Eigen::Matrix<double,DOF,7> DH_table;
        DH_table =  Eigen::Matrix<double,Eigen::Dynamic,7>::Zero(dof,7);
        // Eigen::Matrix<double,DOF,7> DH_table_red;
        DH_table_red = Eigen::Matrix<double,Eigen::Dynamic,7>::Zero(dof,7);
        // Eigen::Matrix<double,DOF,1> q_now;
        q_now = Eigen::VectorXd::Zero(dof);
        // Eigen::Matrix<double,DOF,1> gq;
        gq = Eigen::VectorXd::Zero(dof);
        // Eigen::Matrix<double,DOF,1> G_;
        G_=Eigen::VectorXd::Zero(dof);
        // Eigen::Matrix<double,DOF,1> dq_;//关节速度
        dq_=Eigen::VectorXd::Zero(dof);
        // Eigen::Matrix<double,DOF,1> dq_last_;//关节速度
        dq_last_ = Eigen::VectorXd::Zero(dof);
        // dq_last_.resize(DOF);
        // Eigen::Matrix<double,DOF,1> ddq_;//关节加速度
        ddq_= Eigen::VectorXd::Zero(dof);
        // Eigen::Matrix<double,DOF,1> ddq_last_;//上一次关节加速度
        ddq_last_= Eigen::VectorXd::Zero(dof);
        // Eigen::Matrix<double,6,DOF> jacobi_;//雅可比矩阵
        // jacobi_.resize(6,DOF);
        jacobi_ =Eigen::Matrix<double, 6, -1>::Zero(6,dof);
        // Eigen::Matrix<double,DOF,1> tor_CpM_neton_;
        tor_CpM_neton_= Eigen::VectorXd::Zero(dof);
        // Eigen::Matrix<double,DOF,1> tor_CpM_neton_last_;
        tor_CpM_neton_last_= Eigen::VectorXd::Zero(dof);
    //初始化矩阵大小

        // q_now <<0,0,0,0,0,0;
        T_.resize(dof);
        Pc.resize(dof);
        m_.resize(dof);
        F_T_EE.resize(2);
        T_tool_=Eigen::Matrix<double,4,4>::Identity();
        this->v_.resize(dof);
        this->dw_.resize(dof);
        this->w_.resize(dof);
        this->a_.resize(dof);
        this->Ic_.resize(dof);


        init_size_flag = true;

        this->f_mu_.setZero(dof);//摩擦系数
        this->f_s_.setZero(dof);//静摩擦力
        this->tor_friction_.setZero(dof);//摩擦力



        
    }

    /**
     * @brief 龚陈威根据西西里亚诺那本书定义的参数。
     * @param param 西西里亚诺那本书定义的参数
     * @param dof 自由度数
    */
    void Robot_dynamic::read_dynamic_type_Gong(Eigen::VectorXd param,int dof)
    {
        if (param.size() != dof*11+dof*2)
        {
            printf("\033[1;31;40m 参数数量不对\033[0m \n");
            return;
        }
        //初始化一些参数的大小
        
        this->DOF_=dof;
        for (int i = 0; i < dof; i++)                           //自由度数
        {
            this->m_[i] = param(i*11 + 0);                      //设置m
            this->Pc[i] = param.block<3,1>(i*11+1,0)/m_[i];     //设置质心向量
            this->Ic_[i](0,0) = param(i*11+4);
            this->Ic_[i](1,1) = param(i*11+7);
            this->Ic_[i](2,2) = param(i*11+9);

            this->Ic_[i](0,1) = -param(i*11+5);
            this->Ic_[i](1,0) = -param(i*11+5);

            this->Ic_[i](0,2) = -param(i*11+6);
            this->Ic_[i](2,0) = -param(i*11+6);

            this->Ic_[i](1,2) = -param(i*11+8);
            this->Ic_[i](2,1) = -param(i*11+8);

            this->motor_I_[i] = param(i*11 + 10);               //设置转子惯量

            this->f_s_(i) = param(11*dof + i-1);
            this->f_mu_(i) = param(11*dof + i-1+7);

        }
    }

    void Robot_dynamic::set_DH_table(Eigen::Matrix<double,Eigen::Dynamic,7> dh_table,int dof)
    {
        this->DH_table = dh_table;
        
        for (int i = 0; i < dof; i++)
        {
            // DH_table_red(i,0) = DH_table(i,0)/1000.0f;//ai-1换算成m
            // DH_table_red(i,2) = DH_table(i,2)/1000.0f;//di换算成m
            // DH_table_red(i,1) = DH_table(i,1)*M_PI/180;//alphai-1换算成弧度
            for (int j = 0; j < 7; j++)
            {
                if (j == 0 || j==2)
                {
                    DH_table_red(i,j) = DH_table(i,j)/1000.0f;//ai-1换算成m
                }else if(j==4)//条件变量
                {DH_table_red(i,j) = DH_table(i,j);
                }
                else
                {//换算成弧度
                    DH_table_red(i,j) = DH_table(i,j)*M_PI/180.0f;
                }
                std::cout <<DH_table_red(i,j)<<",";
            }
            std::cout<<std::endl;
        }

    }
    


    /**
     * @brief 根据文件更新内部所有参数
     * @param url 文件路径
     * @param DOF 自由度数目
     * @return 是否成功读取了文件
    */
    bool Robot_dynamic::read_dynamic(const std::string& url,const int& DOF)
    {
        DOF_=DOF;
        // Eigen::Matrix<double,DOF,7> DH_table;
        DH_table =  Eigen::Matrix<double,Eigen::Dynamic,7>::Zero(DOF,7);
        // Eigen::Matrix<double,DOF,7> DH_table_red;
        DH_table_red = Eigen::Matrix<double,Eigen::Dynamic,7>::Zero(DOF,7);
        // Eigen::Matrix<double,DOF,1> q_now;
        q_now = Eigen::VectorXd::Zero(DOF);
        // Eigen::Matrix<double,DOF,1> gq;
        gq = Eigen::VectorXd::Zero(DOF);
        // Eigen::Matrix<double,DOF,1> G_;
        G_=Eigen::VectorXd::Zero(DOF);
        // Eigen::Matrix<double,DOF,1> dq_;//关节速度
        dq_=Eigen::VectorXd::Zero(DOF);
        // Eigen::Matrix<double,DOF,1> dq_last_;//关节速度
        dq_last_ = Eigen::VectorXd::Zero(DOF);
        // dq_last_.resize(DOF);
        // Eigen::Matrix<double,DOF,1> ddq_;//关节加速度
        ddq_= Eigen::VectorXd::Zero(DOF);
        // Eigen::Matrix<double,DOF,1> ddq_last_;//上一次关节加速度
        ddq_last_= Eigen::VectorXd::Zero(DOF);
        // Eigen::Matrix<double,6,DOF> jacobi_;//雅可比矩阵
        // jacobi_.resize(6,DOF);
        jacobi_ =Eigen::Matrix<double, 6, -1>::Zero(6,DOF);
        // Eigen::Matrix<double,DOF,1> tor_CpM_neton_;
        tor_CpM_neton_= Eigen::VectorXd::Zero(DOF);
        // Eigen::Matrix<double,DOF,1> tor_CpM_neton_last_;
        tor_CpM_neton_last_= Eigen::VectorXd::Zero(DOF);
    //初始化矩阵大小

        // q_now <<0,0,0,0,0,0;
        T_.resize(DOF);
        Pc.resize(DOF);
        m_.resize(DOF);
        F_T_EE.resize(2);
        T_tool_=Eigen::Matrix<double,4,4>::Identity();
        this->v_.resize(DOF);
        this->dw_.resize(DOF);
        this->w_.resize(DOF);
        this->a_.resize(DOF);
        this->Ic_.resize(DOF);

        for (int i = 0; i < 2; i++)
        {
        F_T_EE.at(i)<<0,0,0;
        }


        std::ifstream inFile;
        inFile.open(url.c_str());
        if (!inFile.is_open())
        {
            printf("\033[1;31;40m 打开链接 \033[0m \n");//红色
            std::exit(1);
        }

        std::string lineStr;
        printf("\033[1;32;40m测试开始\033[0m  \n");//绿色
        std::vector<std::vector<std::string> > strArray;

        while (getline(inFile, lineStr)) 
        {
            std::cout << "\033[1;32;40m"<< lineStr <<"\033[0m"<< std::endl;
            std::stringstream ss(lineStr);
            std::string str;
            std::vector<std::string> lineArray;
            while (getline(ss, str, ',')){
                lineArray.push_back(str);   //每一行的数组
            }
            strArray.push_back(lineArray);  //包括所有的行
        }

        for (size_t i = 0; i < strArray.size(); i++)
        {
            for (size_t j = 0; j < strArray.at(i).size(); j++)
            {
                std::cout<<strArray.at(i).at(j)<<std::endl;
            }
        }

        int row = DOF;
        int col = 7;
        std::cout<<"以下是度为单位,mm为单位"<<std::endl;
        for (int i = 0; i < row; i++)//从第16行开始 一共6行
        {
            for (int j = 0; j < col; j++)//从第1列开始，一共7列，分别是
            //ai-1，alphai-1，di，thetai，
            //isrevolutor，min，max
            {
                DH_table(i,j)=atof(strArray.at(i+DOF*2+4).at(j+1).c_str());
                std::cout<<DH_table(i,j)<<",";      //DH 参数
            }
            std::cout<<std::endl;//显示参数
        }

        std::cout<<"下面是m为单位和弧度为单位的DH_table"<<std::endl;
        for (int i = 0; i < row; i++)
        {
            // DH_table_red(i,0) = DH_table(i,0)/1000.0f;//ai-1换算成m
            // DH_table_red(i,2) = DH_table(i,2)/1000.0f;//di换算成m
            // DH_table_red(i,1) = DH_table(i,1)*M_PI/180;//alphai-1换算成弧度
            for (int j = 0; j < col; j++)
            {
                if (j == 0 || j==2)
                {
                    DH_table_red(i,j) = DH_table(i,j)/1000.0f;//ai-1换算成m
                }else if(j==4)//条件变量
                {DH_table_red(i,j) = DH_table(i,j);
                }
                else
                {//换算成弧度
                    DH_table_red(i,j) = DH_table(i,j)*M_PI/180.0f;
                }
                std::cout <<DH_table_red(i,j)<<",";
            }
            std::cout<<std::endl;
        }
        
        //下面是质心参数
        std::cout<<"下面是质心参数"<<std::endl;
        row = DOF;
        col = 3;
        for (int i = 0; i < row; i++)
        {
            for (int j = 0; j < col; j++)
            {
                Pc.at(i)(j) = atof(strArray.at(i+3+DOF).at(j+1).c_str())/1000.0f;//换算单位为m
                std::cout<<Pc.at(i)(j)<<",";
            }
            m_.at(i) = atof(strArray.at(i+3+DOF).at(4).c_str());
            std::cout<<m_.at(i)<<",";
            std::cout<<std::endl;
        }

        //下面是惯量矩阵参数
        std::cout<<"下面是惯量矩阵参数"<<std::endl;
        row = DOF;//6行数据
        col = 6;//6列数据
        for (int i = 0; i < DOF; i++)//每一行都是每一个惯量矩阵
        {
            std::vector<double> data_tmp6;
            data_tmp6.clear();
            for (int j = 0; j < col; j++)
            {
                data_tmp6.push_back(atof(strArray.at(i+2).at(j+1).c_str())/1000000.0f);//换算单位为mm2*kg;
            }
            Ic_[i]=Eigen::Matrix<double,3,3>::Zero();

            Ic_[i](0,1)=data_tmp6.at(1);        //XY
            Ic_[i](0,2)=data_tmp6.at(2);        //XZ
            Ic_[i](1,2)=data_tmp6.at(4);        //YZ
            Ic_[i] = Ic_[i] + Ic_[i].transpose().eval();       //加上对称项 eval()可以消除别名问题
            
            Ic_[i](0,0)=data_tmp6[0];           //XX
            Ic_[i](1,1)=data_tmp6[3];           //YY
            Ic_[i](2,2)=data_tmp6[5];           //ZZ
            std::cout<<"惯量参数I("<<i<<")="<<std::endl<<Ic_.at(i)<<","<<std::endl;
        }
        return true;
    }//采集参数

    bool Robot_dynamic::get_joint_isrevolutor(int i)
    {
        if (((uint8_t)DH_table_red(i,4))==1)
        {return true;}
        else
        {return false;}
        
    }


    /** 
     * @brief 函数简要说明-计算n-1到n的齐次坐标变换矩阵(MDH)
     *        //ai-1，alphai-1，di，thetai，
            //isrevolutor，min，max
     * @param n    n计算n-1到n的齐次坐标变换矩阵中的n,如要计算56的变换矩阵，就写6
     * (但是c语言数组是从0开始的，所以平常说的5到6变换矩阵是第六个矩阵，但是从0开始，所以填写5)
     * @param theta        变量的值 @see CTest
     *
     * @return 返回说明
     *     -<em>Eigen::Matrix<double,4,4></em> 4X4的矩阵
     *     -<em>true</em> succeed
     */
    Eigen::Matrix<double,4,4> Robot_dynamic::Ti_cal(const int n,float theta)
    {
        Eigen::Matrix<double,4,4> T;
        float a =this->DH_table_red(n,0);
        float alpha =this->DH_table_red(n,1);
        float d =this->DH_table_red(n,2);
        float thetad = this->DH_table_red(n,3);
        uint8_t isrevolutor =(uint8_t)this->DH_table_red(n,4);
        if (isrevolutor == 1)//计算转动关节
        {
            T<<cos(theta),-sin(theta),0,a,
                sin(theta)*cos(alpha),cos(theta)*cos(alpha),-sin(alpha),-sin(alpha)*d,
                sin(theta)*sin(alpha),cos(theta)*sin(alpha),cos(alpha),cos(alpha)*d,
                0,0,0,1;

            // std::cout<<"fuck!!!!!!!!!!!!! cos(theta) = "<<cos(theta)<<"T(1)="<<T(1)<<std::endl;



        }else//计算移动关节
        {//移动关节的话传进来是变量d
            d = theta;
            T<<cos(thetad),-sin(thetad),0,a,
                sin(thetad)*cos(alpha),cos(thetad)*cos(alpha),-sin(alpha),-sin(alpha)*d,
                sin(thetad)*sin(alpha),cos(thetad)*sin(alpha),cos(alpha),cos(alpha)*d,
                0,0,0,1;
        }
        // cout<<"T"<<n<<"="<<T<<endl;
        this->T_[n] = T;
        return T;
    }

        /** 
     * @brief 函数简要说明-更新关节角度，单位rad
     * @param q_now  类型：Eigen::Matrix<double,DOF,1>类型的
     * @return 返回说明无返回值
     */
    void Robot_dynamic::set_q_now(Eigen::Matrix<double,Eigen::Dynamic,1> q_now)
    {
        this->q_now = q_now;
    }


    /** 
     * @brief 函数简要说明-更新所有齐次坐标变换矩阵(MDH),以及0到n的其次坐标变换矩阵
     *        //ai-1，alphai-1，di，thetai，
            //isrevolutor，min，max
     */
    void Robot_dynamic::T_cal()
    {
        for (int i = 0; i < DOF_; i++)
        {
            //Eigen::Matrix<double,4,4> Robot_dynamic::Ti_cal(const int n,float theta)
            this->T_.at(i) = this->Ti_cal(i,this->q_now(i));
            Eigen::Matrix<double,4,4> temp_0_Ti;
            temp_0_Ti.setIdentity();
            for (int j = 0; j < i; j++)
            {
                temp_0_Ti=temp_0_Ti*T_[j];//坐标变换矩阵
            }
            this->_0T_i[i] =temp_0_Ti * T_[i];//计算了从0坐标系到i坐标系的变换矩阵
            
        // cout<<"T_ = "<< T_.at(i)<<endl;

        }
    }

    /** 
     * @brief 函数简要说明-更新这一个坐标系到上一个坐标系的旋转矩阵
     * @param i  int型，这一个坐标系
     * @return  返回Eigen::Matrix<double,3,3>的3X3的旋转矩阵
     */
    Eigen::Matrix<double,3,3> Robot_dynamic::get_R(int i)
    {
        // cout<<"T_"<<i<<T_.at(i)<<endl;
        if (i>=DOF_)
        {
             Eigen::Matrix<double,3,3> tmp;
             tmp<<1,0,0,0,1,0,0,0,1;
             return tmp;
        }
        else
        {
        // cout<<"get_T("<<DOF-i<<") = "<<endl<<this->T_.at(i)<<endl;
        return (this->T_.at(i).topLeftCorner(3,3));
        }
    }
    /** 
     * @brief 函数简要说明-更新这一个坐标系到世界坐标系下的旋转矩阵,此函数需要_0_Ti计算后才能得到
     * @param i  int型
     * @return  返回Eigen::Matrix<double,3,3>的3X3的旋转矩阵
     */
    Eigen::Matrix<double,3,3> Robot_dynamic::get_0R(int n)
    {
        // Eigen::Matrix<double,3,3> _0Rn;
        // _0Rn<<1,0,0,
        //         0,1,0,
        //         0,0,1;
        // for (int i = 0; i <= n; i++)//从0算到n
        // {
        //     _0Rn = _0Rn*get_R(i);
        // }
        Eigen::Matrix<double,3,3> _0Rn;
        _0Rn =  _0T_i[n].topLeftCorner(3,3);

        return _0Rn;
    }

    /** 
     * @brief 函数简要说明-使用牛顿欧拉法计算力和扭矩，在第i个坐标系下表达
     * @param i 第几个坐标系下的受力，
     * @param F_ip1 第i+1坐标系受到的力在第i+1坐标系下的表达
     * @param T_ip1 第i+1坐标系受到的力矩在第i+1坐标系下的表达
     * @param F_i 第i坐标系受到的力在第i坐标系下的表达
     * @param T_i 第i坐标系受到的力在第i坐标系下的表达 以上三个参数均为3X1向量
     * @return  vector<Eigen::Matrix<double,3,1>>，返回力和扭矩，第0个vector参数是力，第1个vector是力矩，
     * 力和力矩均为3X1的矩阵，在第i个坐标系下表达
     */
    std::vector<Eigen::Matrix<double,3,1>> Robot_dynamic::neton_iter(
        int i,
        Eigen::Matrix<double,3,1> F_ip1,
        Eigen::Matrix<double,3,1> T_ip1,
        Eigen::Matrix<double,3,1> F_i,
        Eigen::Matrix<double,3,1> T_i)
    {
        std::vector<Eigen::Matrix<double,3,1>> F_T;//第0个是F，第1个是T
        
        Eigen::Matrix<double,3,1> Fi_i;//返回的在i坐标系下的 合 力
        Eigen::Matrix<double,3,1> Ti_i;//返回的在i坐标系下的 合 力矩
        //如果在计算最后一个力，后面没有其次坐标旋转矩阵了，那么是旋转矩阵直接取0  
        if (i==DOF_-1)//若传入的是计算最后一个坐标系下的力
        {
            Fi_i = F_i;//受到的力
            //受到的合力矩 = T_i + 第i个坐标系下表达的质心叉乘受力
            Ti_i = T_i + Pc.at(i).cross(F_i);
        }
        else if (i>=0 && i<DOF_-1)
        {
            Eigen::Matrix<double,4,1> Pc_4;//i+1的坐标系下的质心在i+1坐标系下的齐次坐标表示
            Pc_4<< (Pc.at(i+1)),
            1.0;
            // Eigen::Matrix<double,3,1> temp2 = (T_.at(i+1)*Pc_4).topLeftCorner(3,1);
            Eigen::Matrix<double,3,1> temp2 = (T_.at(i+1)).topRightCorner(3,1);

            //第i+1个质心在i这个坐标系下的表示
            Fi_i = get_R(i+1)*F_ip1 + F_i;
            //在第i个坐标系下受力 = i+1坐标系下的受力做旋转变换变换到i坐标系下 + 这个坐标系下的受力
            Ti_i = get_R(i+1)*T_ip1 + T_i+Pc.at(i).cross(F_i)
            +temp2.cross(get_R(i+1)*F_ip1);
            // + temp2.cross(get_R(i+1)*F_ip1);
            //在第i个坐标系下受力矩 = i+1坐标系下的受力矩做旋转变换变换到i坐标系下 + i坐标系下的受力，
            // + 上一个坐标系受到的力
        }
        F_T.resize(2);
        F_T.at(0) = Fi_i;
        F_T.at(1) = Ti_i;
        return F_T;
        
    }

    /** 
     * @brief 函数简要说明-更新对象中的重力补偿
     * 没有输入参数和返回值
     */
    void Robot_dynamic::tor_gravity_cal()
    {
        Eigen::Matrix<double,3,1> Ti_tmp,Fi_tmp,Tip1_tmp,Fip1_tmp;
        Eigen::Matrix<double,3,1> g_w(0,0,-9.8);    //世界坐标系下的g
        Eigen::Matrix<double,3,3> R_0_g,R_0_g_init;            //世界坐标系到当前计算的旋转矩阵
        std::vector<Eigen::Matrix<double,3,1>> F_T; 
        F_T.resize(2);
        R_0_g_init<<1,0,0,0,1,0,0,0,1;
        for (int i = 0; i < DOF_; i++)//循环DOF次
        {
            // R_0_g<<1,0,0,0,1,0,0,0,1;
            R_0_g = R_0_g_init;

            for (int j = 0;j<(DOF_-i);j++)//计算0到坐标系的旋转矩阵
            {
                R_0_g = R_0_g*get_R(j);//右乘相对坐标系旋转矩阵
            }
             //cout<<"get_R("<<DOF-i<<") = "<<endl<<R_0_g<<endl;
            Fi_tmp = m_.at(DOF_-i-1) * R_0_g.transpose()*g_w;//mg在第DOF-i个坐标系(本坐标系)下的表达
            //Ti_tmp = this->Pc.at(DOF-i-1).cross(Fi_tmp);
            Ti_tmp <<0,0,0;
            if (i==0)//最后一个坐标系下
            {
                F_T =neton_iter(DOF_-i-1,this->F_T_EE.at(0),this->F_T_EE.at(1),Fi_tmp,Ti_tmp);

            }else if(i>0 && i<DOF_)//从倒数第二个坐标系开始到第一个坐标系，
            {
            F_T = neton_iter(DOF_-i-1,F_T.at(0),F_T.at(1),Fi_tmp,Ti_tmp);
            }
            // Eigen::Matrix<double,3,1> zi(0,0,1);
            // Eigen::Matrix<double,3,1>(0,0,1);
            G_(DOF_-i-1) = F_T.at(1).transpose() * Eigen::Matrix<double,3,1>(0,0,1);
            // cout<<Eigen::Matrix<double,3,1>(0,0,1)<<endl;
        }
    }



    /** 
     * @brief 获取补偿重力的关节力矩值
     * @return  返回结构体中的补偿重力的关节力矩值
     */
    Eigen::Matrix<double,Eigen::Dynamic,1> Robot_dynamic::get_G_()//获取补偿重力的关节力矩值
    {
        for (int i = 0; i < DOF_; i++)
        {
            /* code */
            if(fabs(G_(i))<0.0001)
            G_(i) = 0;
        }
        
        
        return this->G_;
    }

    /** 
     * @brief 返回数组中第i个关节的关节位置
     * @param i 返回数组中第i个关节的关节位置
     * @return 无
     */
    double Robot_dynamic::get_qi_now(int i)
    {
        return this->q_now(i);
    }

    /** 
     * @brief 通过几何法计算雅可比矩阵
     * @return  无
     */
    void Robot_dynamic::jacobi_cal()
    {
        Eigen::Matrix<double,6,Eigen::Dynamic> Jacobi;
        Jacobi.resize(6,DOF_);
        Eigen::Matrix<double,3,1> zi;//因为用DH变换法，所以容易的到每个转轴都是001
        zi<<0,0,1;
        
        //计算前三行,前三行是计算世界坐标系下线速度的
        //思路：只有此轴运动，其他轴不动计算此轴对末端速度的影响
        for (int i = 0; i < DOF_; i++)//每一列循环
        {
            if (joint_is_rev(i))//如果是旋转关节
            {
            Eigen::Matrix<double,4,4> Ti_tool= Eigen::Matrix<double,4,4>::Identity();
            for (int j = i+1; j <DOF_ ; j++)//Ti都是第i-1个坐标系变换到第i个坐标系的变换矩阵
            {
                Ti_tool =Ti_tool*this->T_.at(j);
            }
            Ti_tool = Ti_tool*this->T_tool_;//乘上末端的工具坐标系变换
            Eigen::Matrix<double,3,1> _iPi_tool= Ti_tool.topRightCorner(3,1);
            Eigen::Matrix<double,3,1> J02_i = this->get_0R(i)*zi.cross(_iPi_tool);//雅可比矩阵1到3行，第i列
            Jacobi.block<3,1>(0,i) = J02_i;
            //////////////////上面是雅可比矩阵前三行的计算
            //////////////////下面是雅可比矩阵后三行的计算
            Jacobi.block<3,1>(3,i) = this->get_0R(i) *zi;
            }//
            else//移动关节
            {
                Eigen::Matrix<double,6,1> J_line_i=Eigen::Matrix<double,6,1>::Zero();
                J_line_i.block<3,1>(0,i) = this->get_0R(i) *zi;//方向就是z轴的方向
                Jacobi.block<6,1>(0,i) = J_line_i;
            }
        }//列循环结束

        this->jacobi_=Jacobi;
    }

    /** 
     * @brief 设置工具坐标系
     * @param T_tool 工具刀尖处相对于最后一个坐标系的相对位置
     * @return  无
     */
    void Robot_dynamic::set_tool(Eigen::Matrix<double,4,4> T_tool)
    {
        this->T_tool_ = T_tool;
    }

    /** 
     * @brief 判断是否是旋转关节，DH表从0列开始，第4列是是否为旋转关节，1是 0否
     * @param i 判断是哪一个关节，从0开始
     * @return  true 是，false 否
     */
    bool Robot_dynamic::joint_is_rev(int i)//第四列是是否是1
    {
        if (this->DH_table(i,4) == 1)
        {
            return true;
        }else
        {
            return false;
        }
    }

    /** 
     * @brief 用迭代法计算速度
     * @param ip1 判断是返回第几个坐标系的速度，将来返回的是i+1个坐标系原点的线速度和角速度
     * @param vi 第i个坐标系原点的线速度
     * @param wi 第i个坐标系原点的角速度
     * @param theta_dot_ip1 第i+1个轴的角速度
     * @return  std::vector<Eigen::Matrix<double,3,1>>类型的第i+1个坐标系的速度和角速度。速度为第一项，角速度为第二项
     * 在i+1坐标系下表达
     */
    std::vector<Eigen::Matrix<double,3,1>>  Robot_dynamic::vel_iter(
            int ip1,
            Eigen::Matrix<double,3,1> vi,
            Eigen::Matrix<double,3,1> wi,
            double theta_dot_ip1
            )
    {
        std::vector<Eigen::Matrix<double,3,1>> v_w_ip1;
        Eigen::Matrix<double,3,1> v_ip1;
        Eigen::Matrix<double,3,1> w_ip1;
        Eigen::Matrix<double,3,1> z_ip1;//i+1轴在i+1坐标系下的单位向量表示
        v_w_ip1.resize(2);
        z_ip1<<0,0,1;
        
        w_ip1 = get_R(ip1).transpose()*wi + theta_dot_ip1*z_ip1;
        Eigen::Vector3d temp=T_.at(ip1).topRightCorner(3,1);
        v_ip1 = get_R(ip1).transpose()*(vi + wi.cross(temp));

        v_w_ip1.at(0) = v_ip1;
        v_w_ip1.at(1) = w_ip1;
        return v_w_ip1;
    }

    /** 
     * @brief 用迭代法计算每个坐标系的线速度和角速度，并更新到对象中去
     * @return  无返回值
     */
    void Robot_dynamic::vel_cal()
    {
        std::vector<Eigen::Matrix<double,3,1>> v_w_ip1;
        v_w_ip1.resize(2);
        for (int i = 0; i < DOF_; i++)
        {
            if (i==0)
            {
                v_w_ip1 = vel_iter(i,Eigen::Vector3d::Zero(),Eigen::Vector3d::Zero(),this->dq_(i));
            }else
            {   
                v_w_ip1 = vel_iter(i,v_w_ip1.at(0),v_w_ip1.at(1),this->dq_(i));
            }
            this->v_.at(i) = v_w_ip1.at(0);
            this->w_.at(i) = v_w_ip1.at(1);
        }
    }

    /** 
     * @brief 用迭代法计算每个坐标系的线速度和角速度，并更新到对象中去，通过dq计算速度，用于忽悠牛顿和欧拉
     * @param dq 
     * @return 返回每个坐标系原点的速度
     */
    std::vector<Eigen::Matrix<double,6,1>> Robot_dynamic::vel_cal(Eigen::VectorXd dq)//通过dq计算速度，用于忽悠牛顿和欧拉
    {
        std::vector<Eigen::Matrix<double,6,1>> v_w_all;
        v_w_all.resize(DOF_);


        std::vector<Eigen::Matrix<double,3,1>> v_w_ip1;
        v_w_ip1.resize(2);
        for (int i = 0; i < DOF_; i++)
        {
            if (i==0)
            {
                // v_w_ip1 = vel_iter(i,Eigen::Vector3d::Zero(),Eigen::Vector3d::Zero(),this->dq_(i));
                v_w_ip1 = vel_iter(i,Eigen::Vector3d::Zero(),Eigen::Vector3d::Zero(),dq(i));
            }else
            {   
                v_w_ip1 = vel_iter(i,v_w_ip1.at(0),v_w_ip1.at(1),dq(i));
                // v_w_ip1 = vel_iter(i,v_w_ip1.at(0),v_w_ip1.at(1),this->dq_(i));
            }
            // this->v_.at(i) = v_w_ip1.at(0);
            // this->w_.at(i) = v_w_ip1.at(1);
            v_w_all.at(i).topLeftCorner(3,1) = v_w_ip1.at(0);
            v_w_all.at(i).bottomLeftCorner(3,1) = v_w_ip1.at(1);
        }
        return v_w_all;

    }

    /**
     * @brief 获取末端雅可比矩阵计算出来的速度
     * 
     * @return Eigen::Matrix<double,6,1> 
     */
    Eigen::Matrix<double,6,1> Robot_dynamic::get_vel_w_jacobe()
    {
        return this->jacobi_*this->dq_;
    }

    /**
     * @brief 返回质心的雅可比矩阵，一共有DOF_个
     * 
     * @return std::vector<Eigen::Matrix<double,6,Eigen::Dynamic>> 
     */
    std::vector<Eigen::Matrix<double,6,Eigen::Dynamic>> Robot_dynamic::get_Pc_jacobe()
    {
        return this->jacobi_ci_;
    }



    Eigen::Matrix<double,6,1> Robot_dynamic::get_vel_w_iter()
    {
        Eigen::Matrix<double,6,1> v_w;
        v_w.topLeftCorner(3,1) = this->v_.at(DOF_-1);
        v_w.bottomLeftCorner(3,1) = this->w_.at(DOF_-1);
        Eigen::Matrix<double,6,6> E_R_R = Eigen::Matrix<double,6,6>::Zero();//转换矩阵
        E_R_R.topLeftCorner(3,3) = get_0R(DOF_-1);
        E_R_R.bottomRightCorner(3,3) = get_0R(DOF_-1);
        v_w = E_R_R*v_w;
        // std::cout<<"w_="<<this->w_.at(DOF-1)<<std::endl;
        // std::cout<<"E_R_R="<<E_R_R<<std::endl;

        return v_w;
    }

    /** 
     * @brief 设置成员的当前关节速度
     * @return  无返回值
     */
    void Robot_dynamic::set_dq_now(Eigen::Matrix<double,Eigen::Dynamic,1> dq)//设置当前关节速度
    {

        this->dq_ = this->tor_filter(dq,30,0.008,&(this->dq_last_));
        
    }

    /** 
     * @brief 用差分计算加速度
     * @param period 运行周期
     * @return  无返回值
     */
    void Robot_dynamic::ddq_cal(ros::Duration period)
    {
        Eigen::Matrix<double, 6, 1> ddq_now = (this->dq_-this->dq_last_)/period.toSec();
        // this->ddq_ = Eigen::Matrix<double, DOF, 1>::Zero();
        // ROS_INFO("period = %f s",period.toSec());
        this->ddq_ = this->tor_filter(ddq_now,5.0,0.008,&(this->ddq_last_));
    }

    /** 
     * @brief 保存一次的角速度
     * @return  无返回值
     */
    void Robot_dynamic::set_last_dq()
    {
        this->dq_last_ = this->dq_;
    }

     /** 
     * @brief 计算i+1的加速度和角加速度
     * @param ip1 第i+1个坐标系计算
     * @param ai 上一个线速度
     * @param dwi 上一个角加速度
     * @param wi 上一个角速度
     * @param dtheta_ip1 这一个轴的关节速度
     * @param ddtheta_ip1 这一个轴的关节角加速度
     * @return  返回i+1的加速度和角加速度
     */
    std::vector<Eigen::Matrix<double,3,1>> Robot_dynamic::a_iter(
            int ip1,
            Eigen::Matrix<double,3,1> ai,
            Eigen::Matrix<double,3,1> dwi,
            Eigen::Matrix<double,3,1> wi,
            double dtheta_ip1,
            double ddtheta_ip1)      //迭代计算加速度
    {
        std::vector<Eigen::Matrix<double,3,1>> a_w_ip1; //输出的结果
        Eigen::Matrix<double,3,1> dw_ip1;               //下一个角速度
        Eigen::Matrix<double,3,1> a_ip1;                //下一个线加速度
        Eigen::Matrix<double,3,1> z_ip1(0,0,1);         //i+1的z轴单位向量
        Eigen::Matrix<double,3,1> i_p_i_ip1;            //i到i+1坐标系向量在第i坐标系下的表达
        i_p_i_ip1 = this->T_.at(ip1).topRightCorner(3,1);//是3X1列向量
        a_w_ip1.resize(2);                              //

        dw_ip1 =this->get_R(ip1).transpose() * dwi      //角加速度 上一个角加速度左乘从i到i+1坐标系下的变换旋转矩阵
        + (this->get_R(ip1).transpose() * wi).cross(dtheta_ip1 * z_ip1) //
        + ddtheta_ip1*z_ip1;
        a_ip1 = this->get_R(ip1).transpose() *(ai + dwi.cross(i_p_i_ip1) + wi.cross(wi.cross(i_p_i_ip1)));
        a_w_ip1.at(0) = a_ip1;
        a_w_ip1.at(1) = dw_ip1;

        return a_w_ip1;
    }

    /** 
     * @brief 更新的加速度和角加速度，此函数必须运行在计算速度和角速度函数之后，因为要用到速度和角速度
     * @return  无返回值
     */
    void Robot_dynamic::a_cal()//计算加速度和角加速度
    {
        Eigen::Matrix<double,3,1> ai=Eigen::Matrix<double,3,1>::Zero(),
                                    dwi=Eigen::Matrix<double,3,1>::Zero(),
                                    wi=Eigen::Matrix<double,3,1>::Zero();
        std::vector<Eigen::Matrix<double,3,1>> a_w_ip1;//第一个是线加速度，第二个是角加速度
        a_w_ip1.resize(2);
        

        for (int ip1 = 0; ip1 < DOF_; ip1++)
        {
            if (ip1==0)
            {
                a_w_ip1 = a_iter(ip1,ai,dwi,wi,this->dq_[ip1],this->ddq_[ip1]);
            }
            else
            {
                a_w_ip1 = a_iter(ip1,a_w_ip1.at(0),a_w_ip1.at(1),this->w_.at(ip1-1),this->dq_[ip1],this->ddq_[ip1]);
            }
            this->a_.at(ip1) = a_w_ip1.at(0);
            this->dw_.at(ip1) = a_w_ip1.at(1);
        }
    }


    /**
     * @brief 计算此条件下的的加速度和角加速度，用来忽悠牛顿和欧拉
     * 
     * @param dq 关节速度
     * @param ddq 关节加速度
     * @return std::vector<Eigen::Matrix<double,6,1>> 自由度数目的vector 每个成员前三行是线加速度，后三行是角加速度，描述的均为坐标系原点
     */
    std::vector<Eigen::Matrix<double,6,1>> Robot_dynamic::a_cal(Eigen::VectorXd dq,Eigen::VectorXd ddq)
    {
        Eigen::Matrix<double,3,1> ai=Eigen::Matrix<double,3,1>::Zero(),
                                    dwi=Eigen::Matrix<double,3,1>::Zero(),
                                    wi=Eigen::Matrix<double,3,1>::Zero();
        std::vector<Eigen::Matrix<double,3,1>> a_w_ip1;//第一个是线加速度，第二个是角加速度
        a_w_ip1.resize(2);

        std::vector<Eigen::Matrix<double,6,1>> v_w_frame;
        v_w_frame = this->vel_cal(dq);//计算新的速度

        std::vector<Eigen::Matrix<double,6,1>> a_dw_all;//要返回的坐标系原点的线加速度和角加速度
        a_dw_all.resize(DOF_);

        for (int ip1 = 0; ip1 < DOF_; ip1++)
        {
            if (ip1==0)
            {
                // a_w_ip1 = a_iter(ip1,ai,dwi,wi,this->dq_[ip1],this->ddq_[ip1]);
                a_w_ip1 = a_iter(ip1,ai,dwi,wi,dq[ip1],ddq[ip1]);

            }
            else
            {
                // a_w_ip1 = a_iter(ip1,a_w_ip1.at(0),a_w_ip1.at(1),this->w_.at(ip1-1),this->dq_[ip1],this->ddq_[ip1]);
                a_w_ip1 = a_iter(ip1,a_w_ip1.at(0),a_w_ip1.at(1),v_w_frame.at(ip1-1).bottomLeftCorner(3,1),dq[ip1],ddq[ip1]);

            }
            a_dw_all.at(ip1).topLeftCorner(3,1) = a_w_ip1.at(0);
            a_dw_all.at(ip1).bottomLeftCorner(3,1) = a_w_ip1.at(1);
        }
        return a_dw_all;
    }
    
    /** 
     * @brief 返回a和w
     * @param i 输入第几个坐标系
     * @return  返回v和w
     */
    std::vector<Eigen::Matrix<double,3,1>> Robot_dynamic::get_a_dw_now(int i)
    {
        std::vector<Eigen::Matrix<double,3,1>> a_dw_now;
        a_dw_now.resize(2);
        a_dw_now.at(0) = this->a_.at(i);
        a_dw_now.at(1) = this->dw_.at(i);
        return a_dw_now;
    }

    /** 
     * @brief 返回第i个坐标系的线加速度
     * @param i 输入第几个坐标系
     * @return  返回v和w
     */
    Eigen::Matrix<double,3,1> Robot_dynamic::get_a_now(int i)//获取第i个坐标系的线加速度在第i个坐标系下的表达
    {
        return this->a_.at(i);
    }

    /** 
     * @brief 返回第i的dw(角加速度)
     * @param i 输入第几个坐标系
     * @return  返回v和w
     */
    Eigen::Matrix<double,3,1> Robot_dynamic::get_dw_now(int i)//获取第i个坐标系的线角速度在第i个坐标系下的表达
    {
        return this->dw_.at(i);
    }

    /** 
     * @brief 计算惯性力+科氏力用牛顿迭代法，必须计算完速度和加速度才能进行这一个计算
     * @param i 返回哪一个坐标系的力和力矩
     * @return  返回 第i个坐标系下的力和力矩，vector第一个是受力，第二个是力矩，
     */
    std::vector<Eigen::Matrix<double,3,1>> Robot_dynamic::get_i_M_C_cal(int i)
    {
        std::vector<Eigen::Matrix<double,3,1>> F_T_M_C;
        F_T_M_C.resize(2);
        F_T_M_C.at(0) = this->m_.at(i) * a_Pc_cal(i);
        F_T_M_C.at(1) = Ic_[i]*dw_[i]
        + w_[i].cross(Ic_[i] *w_[i]);
        
        return F_T_M_C;
    }


    /** 
     * @brief 计算惯性力+科氏力用牛顿迭代法，必须计算完速度和加速度才能进行这一个计算,可以直接指定一些数据进行单独计算
     * @param i 返回哪一个坐标系的力和力矩
     * @param a_dw_Pci 第i个质心的线加速度和角加速度
     * @param v_w_Pci 第i个质心的线速度和角速度
     * @return  返回 第i个坐标系下的力和力矩，vector第一个是受力，第二个是力矩，
     */
    std::vector<Eigen::Matrix<double,3,1>> Robot_dynamic::get_i_M_C_cal(int i,
    Eigen::Matrix<double,6,1> a_dw_Pci,
    Eigen::Matrix<double,6,1> v_w_Pci
    )
    {
        //把分块提前赋好值 要不然会出bug
        Eigen::Vector3d a_pci = a_dw_Pci.topRows(3);
        Eigen::Vector3d dw_pci = a_dw_Pci.bottomRows(3);
        Eigen::Vector3d v_pci = v_w_Pci.topRows(3);//这里的v没有用到，所以传进来的是坐标系原点的线速度对旋转关节也没关系
        Eigen::Vector3d w_pci = v_w_Pci.bottomRows(3);

        std::vector<Eigen::Matrix<double,3,1>> F_T_M_C;//力和力矩，第一项是力，第二项是力矩
        F_T_M_C.resize(2);
        F_T_M_C.at(0) = this->m_.at(i) * a_pci;
        F_T_M_C.at(1) = Ic_[i]*dw_pci
        + w_pci.cross(Ic_[i] *w_pci);
        return F_T_M_C;
    }


    /** 
     * @brief 计算第i个质心在第i个坐标系下的加速度,通过私有参数加速度和角加速度以及角速度计算。
     * @param i 返回哪个杆件的加速度
     * @return  返回第i个质心在第i个坐标系下的加速度
     */
    Eigen::Matrix<double,3,1> Robot_dynamic::a_Pc_cal(int i)
    {
        Eigen::Matrix<double,3,1> a_pc_i;
        a_pc_i = this->a_.at(i) 
        + this->dw_.at(i).cross(this->Pc.at(i)) 
        + this->w_.at(i).cross(this->w_[i].cross(this->Pc[i])); 
        return a_pc_i;
    }

    /**
     * @brief 通过坐标系下的速度和加速度计算，计算所有的 质心用的线加速度 
     * 
     * @param vel_frame DOF_个坐标系下的速度向量，第i个vector是第i个坐标系的信息，成员前三行是线速度，后三行是角速度
     * @param a_dw_frame DOF_个坐标系下的加速度向量，第i个vector是第i个坐标系的信息，成员前三行是线加速度，后三行是角加速度
     * @return std::vector<Eigen::Matrix<double,3,1>> 返回DOF_个连杆质心的线加速度
     */
    std::vector<Eigen::Matrix<double,3,1>> Robot_dynamic::a_Pc_cal(std::vector<Eigen::Matrix<double,6,1>> vel_frame,
                                        std::vector<Eigen::Matrix<double,6,1>> a_dw_frame)
    {

        std::vector<Eigen::Matrix<double,3,1>> a_pc;
        a_pc.resize(DOF_);
        
        for (int i = 0; i < DOF_; i++)
        {
            Eigen::Vector3d ai = a_dw_frame[i].topRows(3);      //获取第i坐标系线加速度
            Eigen::Vector3d dwi = a_dw_frame[i].bottomRows(3);  //获取第i坐标系角加速度
            Eigen::Vector3d wi = vel_frame[i].bottomRows(3);    //获取第i坐标系角速度
            a_pc[i] = ai
            + dwi.cross(this->Pc.at(i))
            + wi.cross(wi.cross(this->Pc.at(i)));
        }
        return a_pc;
        // Eigen::Matrix<double,3,1> a_pc_i;
        // a_pc_i = this->a_.at(i) 
        // + this->dw_.at(i).cross(this->Pc.at(i)) 
        // + this->w_.at(i).cross(this->w_[i].cross(this->Pc[i])); 
        // return a_pc_i;
    }

    /** 
     * @brief 通过牛顿欧拉法迭代更新内部的科氏力+离心力
     * @return  更新内部资源，无返回值
     */
    void Robot_dynamic::tor_M_C_neton_cal_()
    {
        Eigen::Matrix<double,3,1> Ti_tmp,Fi_tmp,Tip1_tmp,Fip1_tmp;

        Eigen::Matrix<double,3,3> R_0_g ,R_0_g_init;            //世界坐标系到当前计算的旋转矩阵
        std::vector<Eigen::Matrix<double,3,1>> F_T;             //计算临时的迭代过来的力和力矩
        F_T.resize(2);
        // R_0_g_init<<1,0,0,0,1,0,0,0,1;
        for (int i = 0; i < DOF_; i++)//循环DOF次
        {
            std::vector<Eigen::Matrix<double,3,1> >  F_T__M_C_i= this->get_i_M_C_cal(DOF_-i-1);
            // R_0_g<<1,0,0,0,1,0,0,0,1;
            if (i==0)//最后一个坐标系下
            {
                F_T =neton_iter(DOF_-i-1,
                Eigen::Matrix<double,3,1>::Zero(),
                Eigen::Matrix<double,3,1>::Zero(),
                F_T__M_C_i.at(0),
                F_T__M_C_i.at(1));
            }else if(i>0 && i<DOF_)//从倒数第二个坐标系开始到第一个坐标系，
            {
            F_T = neton_iter(DOF_-i-1,
            F_T.at(0),
            F_T.at(1),
            F_T__M_C_i.at(0),
            F_T__M_C_i.at(1));
            }
            // Eigen::Matrix<double,3,1> zi(0,0,1);
            // Eigen::Matrix<double,3,1>(0,0,1);
            this->tor_CpM_neton_(DOF_-i-1) = F_T.at(1).transpose() * Eigen::Matrix<double,3,1>(0,0,1);
            // cout<<Eigen::Matrix<double,3,1>(0,0,1)<<endl;
        }
    }


    Eigen::VectorXd Robot_dynamic::tor_M_C_neton_cal_(Eigen::VectorXd dq,
                                                                    Eigen::VectorXd ddq)
    {
        Eigen::Matrix<double,3,1> Ti_tmp,Fi_tmp,Tip1_tmp,Fip1_tmp;
        Eigen::Matrix<double,3,3> R_0_g ,R_0_g_init;            //世界坐标系到当前计算的旋转矩阵
        std::vector<Eigen::Matrix<double,3,1>> F_T;             //计算临时的迭代过来的力和力矩
        F_T.resize(2);
        // R_0_g_init<<1,0,0,0,1,0,0,0,1;
        std::vector<Eigen::Matrix<double,6,1>> vel_frame = this->vel_cal(dq);//计算出每个坐标原点的线速度和角速度
        std::vector<Eigen::Matrix<double,6,1>> a_dw_frame = this->a_cal(dq,ddq);
        std::vector<Eigen::Matrix<double,6,1>> a_dw_Pc = a_dw_frame;//先初始化，此时角加速度是没问题的，但是线加速度需要进行处理
        Eigen::VectorXd tor_CpM_neton;
        tor_CpM_neton.setZero(DOF_);
        
        std::vector<Eigen::Matrix<double,3,1>> a_Pc = this->a_Pc_cal(vel_frame,a_dw_frame);
        for (int i = 0; i < DOF_; i++)
        {
            a_dw_Pc[i].topRows(3) = a_Pc[i];
        }
        

        for (int i = 0; i < DOF_; i++)//循环DOF次
        {
            //这里需要用的是质心速度和质心加速度
            std::vector<Eigen::Matrix<double,3,1> >  F_T__M_C_i= this->get_i_M_C_cal(DOF_-i-1,a_dw_Pc[DOF_-i-1],vel_frame[DOF_-i-1]);
            
            // R_0_g<<1,0,0,0,1,0,0,0,1;
            if (i==0)//最后一个坐标系下
            {
                F_T =neton_iter(DOF_-i-1,
                Eigen::Matrix<double,3,1>::Zero(),
                Eigen::Matrix<double,3,1>::Zero(),
                F_T__M_C_i.at(0),
                F_T__M_C_i.at(1));
            }else if(i>0 && i<DOF_)//从倒数第二个坐标系开始到第一个坐标系，
            {
            F_T = neton_iter(DOF_-i-1,
            F_T.at(0),
            F_T.at(1),
            F_T__M_C_i.at(0),
            F_T__M_C_i.at(1));
            }
            // Eigen::Matrix<double,3,1> zi(0,0,1);
            // Eigen::Matrix<double,3,1>(0,0,1);
            // this->tor_CpM_neton_(DOF_-i-1) = F_T.at(1).transpose() * Eigen::Matrix<double,3,1>(0,0,1);
            tor_CpM_neton(DOF_-i-1) = F_T.at(1).transpose() * Eigen::Matrix<double,3,1>(0,0,1);

            
            // cout<<Eigen::Matrix<double,3,1>(0,0,1)<<endl;
        }
        return tor_CpM_neton;
    }

    /** 
     * @brief 返回对象中的用牛顿欧拉法计算出来的科氏力和离心力的和
     * @return 返回对象中的用牛顿欧拉法计算出来的科氏力和离心力的和
     */
    Eigen::Matrix<double,Eigen::Dynamic,1> Robot_dynamic::get_tor_CpM_neton_()
    {
        return this->tor_CpM_neton_;
    }

    /** 
     * @brief 返回滤波后的力矩向量
     * @param tor_ 
     * @param hz
     * @param period
     * @param tor_last
     * @return 返回对象中的用牛顿欧拉法计算出来的科氏力和离心力的和
     */
    Eigen::Matrix<double,Eigen::Dynamic,1> Robot_dynamic::tor_filter(Eigen::Matrix<double,Eigen::Dynamic,1> tor_,
                                                            double hz,
                                                            double period,
                                                            Eigen::Matrix<double,Eigen::Dynamic,1>* tor_last)
    {
        Eigen::Matrix<double,Eigen::Dynamic,1> tor_out;
        double alpha =2*3.14*hz*period;
        alpha = alpha/(1+alpha);
        tor_out = tor_*alpha + (1-alpha)*(*tor_last);
        *tor_last = tor_out;
        return tor_out;
    }

    /** 
     * @brief 返回上一次的科氏力和离心力的补偿力矩 的指针
     * @return 返回上一次的科氏力和离心力的补偿力矩 的指针
     */
    Eigen::Matrix<double,Eigen::Dynamic,1>* Robot_dynamic::get_tor_CpM_neton_last_ptr()
    {
        return &(this->tor_CpM_neton_last_);
    }

    /**
     * @brief 返回0坐标系到第i个坐标系的齐次矩阵
     * @param i 第i个
     * @return 返回一个4X4的矩阵
     */
    Eigen::Matrix<double,4,4> Robot_dynamic::get_0Ti(int i)
    {
        return this->_0T_i[i];
    }

    /**
     * @brief 计算质心雅可比矩阵，更新内部的质心雅可比矩阵
     * @remark 0.002ms就可以算一次
    */
    void Robot_dynamic::jacobe_ci_cal()
    {
        //计算质心的位置
        std::vector<Eigen::Vector4d> pci_w;
        pci_w.resize(DOF_);
        for (int i = 0; i < DOF_; i++)
        {
            pci_w[i].setZero();                 //置零
            pci_w[i](3)=1;                      //最后一个元素变成1
            pci_w[i].topRows(3) = this->Pc[i];  //加入前三行元素，此时为局部坐标系下的位置向量
            pci_w[i]=this->_0T_i[i] * pci_w[i]; //计算出世界坐标系下的质心向量，注意前面3行是向量，最后一行是1
        }
        for (int i = 0; i < DOF_; i++)
        {
            for (int j = 0; j <= i; j++)
            {
                if (this->get_joint_isrevolutor(i))//旋转关节
                {
                    Eigen::Vector3d Pi_j = pci_w[i].topRightCorner(3,1)-_0T_i[j].topRightCorner(3,1);
                    jacobi_ci_[i].block<3,1>(0,j) = _0T_i[j].block<3,1>(0,2).cross(Pi_j); 
                    jacobi_ci_[i].block<3,1>(3,j) = _0T_i[j].block<3,1>(0,2);
                }
                else//移动关节
                {
                    jacobi_ci_[i].block<3,1>(0,j) = _0T_i[j].block<3,1>(0,2);
                }
            }
        }
    }

    /**
     * @brief 用拉格朗日方法更新内部的惯性矩阵，并返回
     * @return 返回惯性矩阵
    */
    Eigen::MatrixXd Robot_dynamic::M_q_cal_Lagrange()//用拉格朗日方法求解惯性矩阵
    {
        Eigen::MatrixXd M_q;//惯性矩阵
        M_q.setZero(DOF_,DOF_);

        if (!this->init_size_flag)
        {
            printf("\033[1;31;40 出错啦！！！您还没初始化呢！！ \033[0m \n");
            return M_q;
        }
        
        
        Eigen::Matrix<double,3,Eigen::Dynamic> Jci_x;//雅可比矩阵前三行
        Eigen::Matrix<double,3,Eigen::Dynamic> Jci_w;//雅可比矩阵后三行

        Jci_x.setZero(3,DOF_);
        Jci_w.setZero(3,DOF_);

        for (int i = 0; i < DOF_; i++)
        {
            Jci_x = jacobi_ci_[i].topRows(3);//雅可比矩阵前三行
            Jci_w = jacobi_ci_[i].bottomRows(3);//雅可比矩阵后三行
            M_q = M_q+this->m_[i]*(Jci_x.transpose()*Jci_x) + Jci_w.transpose() * get_0R(i) * Ic_[i] *get_0R(i).transpose()*Jci_w;
        }
        M_q_=M_q;//将计算出来的惯性矩阵复制给成员变量
        return M_q;
    }

    /**
     * @brief 每次必须进行的计算，将计算结果储存到对象中，减少重复计算
    */
    void Robot_dynamic::updata_cal()
    {
        T_cal();//运动学建模
        jacobe_ci_cal();
        jacobi_cal();
        tor_M_C_neton_cal_();
        M_q_cal_Lagrange();
    }

    Eigen::MatrixXd Robot_dynamic::get_matrix_Mq()
    {
        return this->M_q_;
    }

    /**
     * @brief 设置摩擦参数到对象中
     * 
     * @param f_param 前dof个是静摩擦参数，后dof个是mu
     */
    void Robot_dynamic::set_friction_param(Eigen::VectorXd f_param)
    {
        this->f_s_ = f_param.topRows(DOF_);
        this->f_mu_= f_param.bottomRows(DOF_);
    }

    /**
     * @brief 设置摩擦参数到对象中
     * 
     * @param f_s 静摩擦参数
     * @param f_mu mu
     */
    void Robot_dynamic::set_friction_param(Eigen::VectorXd f_s,Eigen::VectorXd f_mu)
    {
        this->f_s_ = f_s;
        this->f_mu_= f_mu;
    }

    /**
     * @brief 计算摩擦力，更新内部参数，如果想补偿，请加负号
     * 
     * @return Eigen::VectorXd 摩擦力正负号代表方向，不是补偿力矩！！！补偿请加负号
     */
    Eigen::VectorXd Robot_dynamic::friction_cal()                     
    {
        this->tor_friction_ = -(this->f_mu_.cwiseProduct(this->dq_) +  dq_.cwiseSign()*this->f_s_);
        return tor_friction_;
    }

    /**
     * @brief 通过外参计算摩擦力,，如果想补偿，请加负号
     * 
     * @param dq 传进去的关节角速度
     * @return Eigen::VectorXd 摩擦力正负号代表方向，不是补偿力矩！！！补偿请加负号
     */
    Eigen::VectorXd Robot_dynamic::friction_cal(Eigen::VectorXd dq)  
    {
        Eigen::VectorXd tor_friction;
        tor_friction = -(this->f_mu_.cwiseProduct(dq) +  dq.cwiseSign()*this->f_s_);
        return tor_friction;
    }

}//namespace xj_dy_ns 