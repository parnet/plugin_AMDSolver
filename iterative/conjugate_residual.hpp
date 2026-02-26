#pragma once

#include <rocalution.hpp>

#include "lib_algebra/operator/linear_solver/external_solvers/external_solvers.h"



#include "../iterative//iterative_solver.hpp"
#include "../solver.hpp"
#include "../local_data.hpp"


namespace ug {
template <typename TAlgebra>
class AMD_ConjugateResidual : public IExternalSolver <TAlgebra>,
                          public AMD_IterativeLinearSolver {
public:
    using matrix_t = CPUAlgebra::matrix_type;
    using vector_t = CPUAlgebra::vector_type;


    AMD_ConjugateResidual() = default;

    virtual ~AMD_ConjugateResidual(){
        clear();
    };

    [[nodiscard]] const char * double_name() const override {  return "AMD_BiCGStab"; };

    void double_init(const matrix_t &mat) override ;

    bool double_apply(vector_t &c, const vector_t &d) override ;

    void set_preconditioner(SmartPtr<AMD_Solver> precond){
        sp_precond_ = precond;
        solver_.SetPreconditioner(sp_precond_->get_solver());
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
}