#pragma once

#include <rocalution.hpp>

#include "../iteration_control.hpp"
#include "../solver.hpp"

namespace ug {
class AMD_IterativeLinearSolver : public AMD_Solver{
public:
    using amd_iterativelinearsolver = rocalution::IterativeLinearSolver<
    rocalution::LocalMatrix<double>,
    rocalution::LocalVector<double>,
    double
>;

    ~AMD_IterativeLinearSolver() override = default;


    void set_conv_check(SmartPtr<AMD_IterationControl > conv_check) {
        conv_check_ = conv_check;
    }


    void set_preconditioner(SmartPtr<AMD_Solver> solver) {
        preconditioner_ = solver;
    };

    /*
    void record_residual_history() {
        RecordResidualHistory();
    };

    void record_history_to_file(std::string filename) {
        RecordHistory();
    };

    void get_iteration_count() {
        GetIterationCount();
    };

    void get_current_residual() {
        GetCurrentResidual();
    };

    void get_solver_status() {
        GetSolverStatus();
    };

    void get_amax_residual_index() {
        GetAmaxResidualIndex();
    };*/


    SmartPtr<AMD_Solver> preconditioner_;
    SmartPtr<AMD_IterationControl> conv_check_;

};

}
