#include "mean_filter.h"


namespace xj_dy_ns
{
    MeanFilter::MeanFilter()
    {
        ;
    }
    MeanFilter::MeanFilter(int window_size,int dof) {
        window_size_ = window_size;
        dof_ = dof;
        sum_.setZero(dof);
        index_ = 0;

        buffer_.resize(window_size_);
        for (int i = 0; i < window_size_; i++) {
            buffer_[i].setZero(dof);
        }
    }

    void MeanFilter::init(int window_size,int dof)
    {
        window_size_ = window_size;
        dof_ = dof;
        sum_.setZero(dof);
        index_ = 0;

        buffer_.resize(window_size_);
        for (int i = 0; i < window_size_; i++) {
            buffer_[i].setZero(dof);
        }
    }

    Eigen::VectorXd MeanFilter::Filter(Eigen::VectorXd x) 
    {
        sum_ = sum_ - buffer_[index_] + x;
        buffer_[index_] = x;
        index_ = (index_ + 1) % window_size_;
        return sum_ / window_size_;
    }


}