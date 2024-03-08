#include "Cubic_Spline.h"

namespace xj_dy_ns
{

    Cubic_Spline::Cubic_Spline(/* args */)
    {
    }
    
    Cubic_Spline::~Cubic_Spline()
    {
    }

    void Cubic_Spline::init(Eigen::MatrixXd point,
        Eigen::VectorXd h,
        BOUND_type bound_type,
        Eigen::Matrix<double,2,Eigen::Dynamic> bound)
    {
        // 先做数据检查
        // 数据的自由度
        Dof_ = point.cols();
        // 点的个数
        this->num_point_ = point.rows();
        this->point_ = point;
        // 时间段的个数
        num_time_ = num_point_-1;
        // 重塑M的特征，是一个dof列的向量组 AM=d,求解M
        M_.resize(num_point_,Dof_);
        d_.resize(num_point_,Dof_);

        
        if (num_time_ == h.size())
        {
            this->h_ = h;
            time_all_ =0.0;
            for (size_t i = 0; i < h_.size(); i++)
            {
                time_all_+=h_(i);
            }
            
        }else
        {
            std::cout<<"Cubic_Spline:: error!! number of time is error"<<std::endl;
            return;
        }
        // 判断边界条件是否正常
        if (bound.cols() != Dof_)
        {
            std::cout<<"Cubic_Spline:: error!! number of dof of bound is error"<<std::endl;
            return;
        }else
        {
            bound_ = bound;
        }

        // 计算d1到dn-1 这里的下标比较混乱，小心一点
        for (int i = 0; i < num_point_-2; i++)
        {
            d_.row(i+1) = (6.0/(h_(i)+h_(i+1))) 
            * ( (point_.row(i+2) - point_.row(i+1))/h_(i+1) 
                - (point_.row(i+1) - point_.row(i))/h_(i) );
        }//还有两个，在后面搞

        // 判断时间段个数是否正常
        mu_.resize(num_time_-1);
        lambda_.resize(num_time_-1);

        // 计算mu 和lanmbda 这个的数量比时间段还要少1
        for (int i = 0; i < num_time_-1; i++)
        {
            mu_(i) = h_(i)/(h_(i)+h_(i+1));
            lambda_(i) = 1-mu_(i);
        }
        switch (bound_type)
        {
        case VELOCITY:
            {
                /*
                [2  1  0       0     ... 0
                 u1 2  lambda1
                 0  u2 2      lambad2... 0
                 ...
                 0  0          u_n-1 2 lambda_n-1 
                 0  0                1  2
                ] * M = d
                主对角元素为2 右上次对角线元素是1 lambda 1..... lambda_n-1
                左下次对角线元素是u1 。。。 un-1 1
                
                */
               Eigen::MatrixXd A = 2 * Eigen::MatrixXd::Identity(num_point_,num_point_);
               A(0,1) = 1;
               A(num_time_,num_time_-1) = 1;
               for (size_t i = 0; i < num_time_-1; i++)
               {
                // 对角线右上方
                A(i+1,i+2) = lambda_(i);
                // 对角线左下方
                A(i+1,i) = mu_(i);
               }

                // 计算d0和 dn
                d_.row(0) = 6.0/h_(0)*((point_.row(1)-point_.row(0))/h_(0) - bound_.row(0));
                d_.row(num_point_-1) = 6.0/h_(num_time_-1) 
                    * (bound_.row(1) - (point_.row(num_point_-1)- point_.row(num_point_-2))/h_(num_time_-1));

                // 对于稠密矩阵，使用常规的求解器 求解
                M_ = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(d_);
                
               
            }
            break;

        case ACC:
        break;
        
        default:
            break;
        }
        ;
    }

    Eigen::VectorXd Cubic_Spline::cal(double t)
    {
        // 区间左侧时间和区间右侧时间
        double t_l = 0.0,t_r = 0.0;
        
        Eigen::ArrayXd out;
        out.resize(Dof_);

        // 区间左侧点标号
        int index_l = 0;
        // 判断是否在区间中
        if (t>=time_all_)
        {
            // t_r = time_all_;
            // t_l = time_all_ - h_(h_.size()-1);
            return this->point_.bottomRows(1).transpose();
        }
        else if (t<=0)
        {
            return this->point_.row(0).transpose();

        }
        else
        {
            for (int i = 0; i < num_time_; i++)
            {
                t_r+=this->h_(i);
                if (t<t_r)
                {
                    index_l = i;
                    break;
                }
                t_l+=this->h_(i);
            }
        }
        
        
        
            // 检查过了
            Eigen::MatrixXd result = pow((t_r-t),3)/(6.0*h_(index_l))*this->M_.row(index_l) 
            +pow((t-t_l),3)/(6.0*h_(index_l))*M_.row(index_l+1)
            // +(point_.row(index_l)-pow(h_(index_l),2)/6.0*M_.row(index_l))
            +( point_.row(index_l)- pow(h_(index_l),2)/6.0 *M_.row(index_l) )*(t_r - t)/h_(index_l)
            +(point_.row(index_l+1)-pow(h_(index_l),2)/6.0 *M_.row(index_l+1)) * (t - t_l)/h_(index_l); 
            // std::cout<<"t_l = "<<t_l<<",,,t_r = "<<t_r<<std::endl ;
        
        return result.transpose();
    }

    Cubic_Spline::Cubic_Spline(Eigen::MatrixXd point,
        Eigen::VectorXd h,
        BOUND_type bound_type,
        Eigen::Matrix<double,2,Eigen::Dynamic> bound)
    {
        init(point,h,bound_type,bound);
    }

    Eigen::MatrixXd Cubic_Spline::cal(Eigen::VectorXd t)
    {
        Eigen::MatrixXd result;
        // std::cout<<"Dof_ = "<<Dof_<<std::endl;
        result.resize(t.size(),this->Dof_);
        // std::cout<<"result.rows() = "<<result.rows()<<std::endl;
        for (int i = 0; i < result.rows(); i++)
        {
            result.row(i) = cal(t(i)).transpose();
            // std::cout<< "t "<<i<<"=" <<t(i)<<std::endl<<"cal(t(i)) = "<<cal(t(i))<<std::endl;
        }
        return result;
    }





}