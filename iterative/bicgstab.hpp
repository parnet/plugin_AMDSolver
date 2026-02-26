#pragma once

#include <rocalution.hpp>

#include "lib_algebra/operator/linear_solver/external_solvers/external_solvers.h"



#include "../iterative/iterative_solver.hpp"
#include "../solver.hpp"
#include "../local_data.hpp"


namespace ug{

template <typename TAlgebra>
class AMD_BiCGStab : public IExternalSolver <TAlgebra>,
            public AMD_IterativeLinearSolver {
public:
    using matrix_t = CPUAlgebra::matrix_type;
    using vector_t = CPUAlgebra::vector_type;

    AMD_BiCGStab() = default;

    ~AMD_BiCGStab() override {
        clear();

        };

    [[nodiscard]] const char * double_name() const override {  return "AMD_BiCGStab"; };

    void double_init(const matrix_t &mat) override ;

    bool double_apply(vector_t &c, const vector_t &d) override ;

    void set_preconditioner(SmartPtr<AMD_Solver> precond){ // todo move to baseclass
        sp_precond_ = precond;
        solver_.SetPreconditioner(sp_precond_->get_solver()); // add in init?
    }


    amd_solver& get_solver() override {
        return  solver_;
    }
private:

    void clear() {
        solver_.Clear();
        sp_precond_->get_solver().Clear();
    }
private:
    LocalData data;

    rocalution::BiCGStab<
        rocalution::LocalMatrix<double>,
        rocalution::LocalVector<double>,
        double
    > solver_;

    SmartPtr<AMD_Solver> sp_precond_;
};

template <typename TAlgebra>
void AMD_BiCGStab<TAlgebra>::double_init(const matrix_t &mat){

        /*mat.copy_crs(size_t &numRows, size_t &numCols,
                std::vector<value_type> &argValues, std::vector<int> &argRowStart,
                std::vector<int> &argColInd)*/
        solver_.Clear();
        data.init(mat);
        solver_.SetOperator(data.A_);

        //solver_.InitMaxIter(1000);
        //solver_.InitTol(1e-14, 1e-10, 1e+10); // absolut, relative, divergence
        this->solver_.InitTol(conv_check_->tolerance_absolute, conv_check_->tolerance_reduction, conv_check_->tolerance_divergence);
        this->solver_.InitMaxIter(conv_check_->iterations_max);
        this->solver_.SetResidualNorm(conv_check_->resnorm_);
        //solver_.SetResidualNorm(2); // L2 norm
        data.sync_operator();
        solver_.Build();
}

template <typename TAlgebra>
bool AMD_BiCGStab<TAlgebra>::double_apply(vector_t &c, const vector_t &d){
    data.double_pre_apply(c,d);
    // Solve A * x = b
    solver_.Solve(data.b_, &data.x_);
    data.double_post_apply(c,d); // result in c
    return true;
}


}