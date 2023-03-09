#include <iostream>
#include <vector>
#include <Eigen/Eigen>
namespace xj_dy_ns
{

    class MeanFilter {
    public:
        MeanFilter();
        MeanFilter(int window_size,int dof);
        void init(int window_size,int dof);

        Eigen::VectorXd Filter(Eigen::VectorXd x);

    private:
        int window_size_; // 滤波器窗口大小
        Eigen::VectorXd sum_; // 滤波器窗口内所有数据的总和
        int index_; // 当前数据在滤波器窗口中的位置
        int dof_;
        std::vector<Eigen::VectorXd> buffer_; // 滤波器窗口缓存
    };


}