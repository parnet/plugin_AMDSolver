#pragma once

#include <rocalution.hpp>

#include "lib_algebra/operator/linear_solver/external_solvers/external_solvers.h"

#include "../local_data.hpp"
#include "../preconditioner/preconditioner.hpp"
#include "../iterative/iterative_solver.hpp"
#include "../solver.hpp"

#include "../iterative//iterative_solver.hpp"

#include <vector>

namespace ug{

template <typename TAlgebra>
class AMD_RugeStuebenAMG : public IExternalSolver <TAlgebra>,
    public AMD_IterativeLinearSolver {
public:
    using matrix_t = CPUAlgebra::matrix_type;
    using vector_t = CPUAlgebra::vector_type;

    AMD_RugeStuebenAMG() = default;

    ~AMD_RugeStuebenAMG() override {
        clear();

        };

    const char * double_name() const override {  return "AMD_BiCGStab"; };

    void double_init(const matrix_t &mat) override ;

    bool double_apply(vector_t &c, const vector_t &d) override ;


    amd_solver& get_solver() {
        return solver_;
    }


    void set_preconditioner(SmartPtr<AMD_Preconditioner> precond){
        sp_precond_ = precond;
    }

    void set_smoother(SmartPtr<AMD_IterativeLinearSolver> smoother){
        sp_smoother_ = smoother;
        //solver_.SetSmoother(&sp_smoother_->get_solver())
    }
    void set_solver(SmartPtr<AMD_Solver> solver){
        sp_solver_ = solver;
        solver_.SetPreconditioner(sp_solver_->get_solver());
    }

    void set_coarsest_level(int coarse_size){
        solver_.SetCoarsestLevel(coarse_size);
    }

    /// base amg


    void set_strength_threashold(float eps){
        solver_.SetStrengthThreshold(eps);
    }

    void set_coarsening_strategy(const char * strategy){
        if (strcmp(strategy,"greedy")  == 0){
            solver_.SetCoarseningStrategy(rocalution::CoarseningStrategy::Greedy);
        } else if (strcmp("pmis",strategy) == 0){
            solver_.SetCoarseningStrategy(rocalution::CoarseningStrategy::PMIS);
        }
    }

    void set_interpolation_type(const char * type){
        if (strcmp(type,"direct") == 0 ){
            solver_.SetInterpolationType(rocalution::InterpolationType::Direct);
        } else if (strcmp(type,"extpi") == 0){
            solver_.SetInterpolationType(rocalution::InterpolationType::ExtPI);
        }
    }

    void set_interpolation_ff1_limit(bool ff1){
        solver_.SetInterpolationFF1Limit(ff1);
    }

    void re_build_numeric(){
        solver_.ReBuildNumeric();
    }

    void set_cycle(const char * name) {
        if (strcmp(name,"V")==0){
            solver_.SetCycle(rocalution::_cycle::Vcycle);
        } else if (strcmp(name,"F")==0) {
            solver_.SetCycle(rocalution::_cycle::Fcycle);
        } else if (strcmp(name, "W")==0) {
            solver_.SetCycle(rocalution::_cycle::Wcycle);
        } else if (strcmp(name,"K")==0) {
            solver_.SetCycle(rocalution::_cycle::Kcycle);
        }
    }

    void set_smooth_num(int pre, int post){
        solver_.SetSmootherPreIter(pre);
        solver_.SetSmootherPreIter(post);
    }

private:

    void clear() {
        solver_.Clear();
        sp_precond_->get_solver().Clear();
    }
private:
    LocalData data;

    rocalution::RugeStuebenAMG<
        rocalution::LocalMatrix<double>,
        rocalution::LocalVector<double>,
        double
    > solver_;

    SmartPtr<AMD_Preconditioner> sp_precond_;
    SmartPtr<AMD_IterativeLinearSolver> sp_smoother_;
    SmartPtr<AMD_Solver> sp_solver_;
};

template <typename TAlgebra>
void AMD_RugeStuebenAMG<TAlgebra>::double_init(const matrix_t &mat){


    /*mat.copy_crs(size_t &numRows, size_t &numCols,
            std::vector<value_type> &argValues, std::vector<int> &argRowStart,
            std::vector<int> &argColInd)*/
    solver_.Clear();

    data.init(mat);


    solver_.SetOperator(data.A_);

    //sp_precond_->get_preconditioner().SetOperator(data.A_);

    //solver_.SetPreconditioner(sp_precond_->get_preconditioner());
    solver_.InitMaxIter(1000);
    solver_.InitTol(1e-14, 1e-10, 1e+10); // absolut, relative, divergence
    solver_.SetResidualNorm(2); // L2 norm
    data.sync_operator();
    solver_.Build();

}

template <typename TAlgebra>
bool AMD_RugeStuebenAMG<TAlgebra>::double_apply(vector_t &c, const vector_t &d){

    data.double_pre_apply(c,d);

    // Solve A * x = b
    solver_.Solve(data.b_, &data.x_);

    data.double_post_apply(c,d); // result in c
    return true;
}


}