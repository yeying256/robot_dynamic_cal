
#include "impedance_controller.h"
namespace xj_dy_ns
{


    ImpedanceController::ImpedanceController(/* args */)
    {

    }
    ImpedanceController::ImpedanceController(int dof)
    {
        init_param(dof);
    }

    /**
     * @brief 初始化成员参数矩阵的大小（阶数），
     * 
     * @param dof 
     * @return true 
     * @return false 
     */
    bool ImpedanceController::init_param(int dof)
    {
        this->DOF_ = dof;
        this->Md_.setIdentity(dof,dof);
        this->Kd_.setIdentity(dof,dof);
        this->Dd_.setIdentity(dof,dof);
        
    }


    
    ImpedanceController::~ImpedanceController()
    {
    }


    
}//namespace xj_dy_ns
