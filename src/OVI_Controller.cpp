#include "OVI_Controller.h"



namespace xj_dy_ns
{

    OVI_Controller::OVI_Controller(/* args */)
    {
    }
    
    OVI_Controller::~OVI_Controller()
    {
    }

    std::vector<Eigen::MatrixXd> OVI_Controller::P_matrix_cal(double dt,
                                                        std::vector<Eigen::MatrixXd> Md,
                                                        std::vector<Eigen::MatrixXd> Dd,
                                                        std::vector<Eigen::MatrixXd> Kd,
                                                        Eigen::MatrixXd S,
                                                        Eigen::MatrixXd Q,
                                                        Eigen::MatrixXd R
                                                        )
    {
        // printf("\033[1;31;40m  \n");
        // std::cout<<"这里是OVI_Controller::P_matrix_cal的开始部分="<<std::endl;
        // printf(" \033[0m \n");
        int len = Md.size();//返回序列长度
        int dof = Md[0].rows();//返回行数获取自由度
        std::vector<Eigen::MatrixXd> P;//创建P矩阵
        if (len!=Dd.size() || len!=Kd.size())
        {
            printf("\033[1;31;40m  \n");
            std::cout<<"这里是OVI_Controller::P_matrix_cal函数,这里的输入参数Md,Dd,Kd序列长度不一致";
            printf(" \033[0m \n");
        }
        P.resize(len);
        P.back()=S;
        for (int i = len-1; i >0; i--)
        {
            // printf("\033[1;31;40m  \n");
            // std::cout<<"这里是OVI_Controller::P_matrix_cal的循环代码 此时i="<<i<<std::endl;
            // printf(" \033[0m \n");
            Eigen::MatrixXd Md_tmp,Dd_tmp,Kd_tmp,P_tmp;
            P_tmp = P[i];
            Md_tmp = Md[i];
            Dd_tmp = Dd[i];
            Kd_tmp = Kd[i];
            Eigen::MatrixXd A;
            Eigen::MatrixXd B;
            A.setZero(2*dof,2*dof);
            B.setZero(2*dof,dof);
            A.topLeftCorner(dof,dof).setZero();
            A.topRightCorner(dof,dof).setIdentity();
            A.bottomLeftCorner(dof,dof) = -Md_tmp.inverse() * Kd_tmp;
            A.bottomRightCorner(dof,dof) = -Md_tmp.inverse() *Dd_tmp;
            B.topRows(dof).setZero();
            B.bottomRows(dof).setIdentity();

            // printf("\033[1;31;40m  \n");
            // std::cout<<"这里是OVI_Controller::P_matrix_cal的循环代码 此时计算完了A和B"<<std::endl;
            // printf(" \033[0m \n");

            Eigen::MatrixXd X1= -P_tmp*A -A.transpose()*P_tmp + P_tmp*B*R.inverse()*B.transpose()*P_tmp-Q;
            // printf("\033[1;31;40m  \n");
            // std::cout<<"这里是OVI_Controller::P_matrix_cal的循环代码 计算X1"<<std::endl;
            // printf(" \033[0m \n");
            Eigen::MatrixXd P_tmp_X1=P_tmp+0.5*dt*X1;
            // printf("\033[1;31;40m  \n");
            // std::cout<<"这里是OVI_Controller::P_matrix_cal的循环代码 计算P_tmp_X1"<<std::endl;
            // printf(" \033[0m \n");
            Eigen::MatrixXd X2 = -P_tmp_X1*A - A.transpose()*P_tmp_X1 +P_tmp_X1*B*R.inverse()*B.transpose()*P_tmp_X1 -Q;
            Eigen::MatrixXd P_tmp_X2 =P_tmp+0.5*dt*X2;
            // printf("\033[1;31;40m  \n");
            // std::cout<<"这里是OVI_Controller::P_matrix_cal的循环代码 计算P_tmp_X2"<<std::endl;
            // printf(" \033[0m \n");
            Eigen::MatrixXd X3 = -P_tmp_X2*A - A.transpose()*P_tmp_X2 + P_tmp_X2*B*R.inverse()*B.transpose()*P_tmp_X2 - Q;
            Eigen::MatrixXd P_tmp_X3 =(P_tmp+dt*X3);
            // printf("\033[1;31;40m  \n");
            // std::cout<<"这里是OVI_Controller::P_matrix_cal的循环代码 计算P_tmp_X3"<<std::endl;
            // printf(" \033[0m \n");
            Eigen::MatrixXd X4 = -P_tmp_X3*A - A.transpose()*P_tmp_X3 + P_tmp_X3*B*R.inverse()*B.transpose()*P_tmp_X3 - Q;
            // printf("\033[1;31;40m  \n");
            // std::cout<<"这里是OVI_Controller::P_matrix_cal的循环代码 计算X4"<<std::endl;
            // printf(" \033[0m \n");
            Eigen::MatrixXd pl = P_tmp - dt/6 *(X1 + 2*X2 + 2*X3 + X4);
            P[i-1]=pl;
            // printf("\033[1;31;40m  \n");
            // std::cout<<"这里是OVI_Controller::P_matrix_cal的循环代码 计算完了pl"<<std::endl;
            // printf(" \033[0m \n");
        }
        return P;
    }
    Eigen::VectorXd OVI_Controller::OVI_tor_cal(Eigen::VectorXd epsil,
                                            Eigen::VectorXd d_epsil,
                                            Eigen::MatrixXd P,
                                            double dt)
    {

        int dof = epsil.rows();
        // printf("\033[1;31;40m  \n");
        // std::cout<<"这里是OVI_Controller::OVI_tor_cal 计算了dof="<<dof<<std::endl;
        // printf(" \033[0m \n");
        Eigen::MatrixXd R = dt*Eigen::MatrixXd::Identity(dof,dof);
        
        // std::cout<<"\033[1;31;40m "<<"这里是OVI_Controller::OVI_tor_cal 测试1 R="<<R<<"\033[1;31;40m "<<std::endl;
        Eigen::MatrixXd B;
        B.setZero(2*dof,dof);
        B.topRows(dof).setZero();
        B.bottomRows(dof).setIdentity();
        Eigen::VectorXd x;
        x.setZero(2*dof);
        x.topRows(dof) = epsil;
        x.bottomRows(dof) = d_epsil;
        Eigen::MatrixXd K= R.inverse()*B.transpose()*P;
        // std::cout<<"\033[1;31;40m "<<"这里是OVI_Controller::OVI_tor_cal 测试2 K="<<K<<"\033[1;31;40m "<<std::endl;
        Eigen::VectorXd miu = -K*x;
        return miu;
    }
    

}