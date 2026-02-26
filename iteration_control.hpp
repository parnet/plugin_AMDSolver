#pragma once

namespace ug {

class AMD_IterationControl {
public:
    /*void init(double abs_tol, double rel_tol, double div_tol, int max_iter) {
        tolerance_absolute = abs_tol;
        tolerance_reduction = rel_tol;
        tolerance_divergence = div_tol;

        iterations_max = max_iter;
    };*/

    void init(double abs_tol, double rel_tol, double div_tol, int min_iter, int max_iter) {
        tolerance_absolute = abs_tol;
        tolerance_reduction = rel_tol;
        tolerance_divergence = div_tol;

        iterations_min = min_iter;
        iterations_max = max_iter;
    };

    void init_tol(double abs_tol, double rel_tol, double div_tol) {
        tolerance_absolute = abs_tol;
        tolerance_reduction = rel_tol;
        tolerance_divergence = div_tol;
    };

    void init_min_iter(int min_iter) {
        iterations_min = min_iter;
    };

    void init_max_iter(int max_iter) {
        iterations_max = max_iter;
    };

    void set_residual_norm(int resnorm) {
        resnorm_ = resnorm;
    };

    void set_verbose(bool verbose) {
        verbose_ = verbose;
    };



    int resnorm_ = 2;

    int iterations_max = 100;
    int iterations_min = 0;

    double tolerance_absolute = 1e-14;
    double tolerance_reduction = 1e-10;
    double tolerance_divergence = 1e10;
    bool verbose_ = false;
};

}