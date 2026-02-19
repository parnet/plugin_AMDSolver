#pragma once

#include <rocalution.hpp>

#include "lib_algebra/operator/linear_solver/external_solvers/external_solvers.h"
#include <vector>

namespace ug{

template <typename TAlgebra>
class AMD_BiCGStab : public IExternalSolver <TAlgebra> {
public:
    using matrix_t = CPUAlgebra::matrix_type;
    using vector_t = CPUAlgebra::vector_type;

    AMD_BiCGStab(){

        };

    virtual ~AMD_BiCGStab(){
        clear();
        rocalution::stop_rocalution();
        };

    const char * double_name() const override {  return "AMD_BiCGStab"; };

    void double_init(const matrix_t &mat) override ;

    bool double_apply(vector_t &c, const vector_t &d) override ;
private:

    void clear() {
        solver_.Clear();
        precond_.Clear();
        A_.Clear();
        b_.Clear();
        x_.Clear();
    }
private:

    int n_ = 0;
    size_t stn_   = 0;
    int nnz_ = 0;

    std::vector<double> argValues;
    std::vector<int> argRowStart;
    std::vector<int> argColInd;

    rocalution::LocalMatrix<double> A_;
    rocalution::LocalVector<double> b_;
    rocalution::LocalVector<double> x_;

    rocalution::BiCGStab<
        rocalution::LocalMatrix<double>,
        rocalution::LocalVector<double>,
        double
    > solver_;

    rocalution::ILU<
        rocalution::LocalMatrix<double>,
        rocalution::LocalVector<double>,
        double
    > precond_; // todo ::SetLevel with k

    /*Jacobi<LocalVector<double>, LocalMatrix<double>, double> precond;
    precond.Init(A_gpu);
    precond.SetMaxIter(1);

     // simple Jacobi*/
     // ILUT
    // IC

    /*         rocalution::AMG< rocalution::LocalMatrix<double>,
                             rocalution::LocalVector<double>,
                             double > amg;

            amg.SetOperator(A_);

            // Typical AMG settings obs! 3-10 times matrix memory
            //  as solver:
            amg.SetCoarseningStrategy(rocalution::AMG_RS);
            amg.SetInterpolationStrategy(rocalution::AMG_PMIS);
            amg.SetRelaxationType(rocalution::AMG_FSAI);
            amg.Build();
            // as preconditioner
            solver.SetPreconditioner(amg_precond);
            solver.Build()

    */
};

template <typename TAlgebra>
void AMD_BiCGStab<TAlgebra>::double_init(const matrix_t &mat){


    /*mat.copy_crs(size_t &numRows, size_t &numCols,
            std::vector<value_type> &argValues, std::vector<int> &argRowStart,
            std::vector<int> &argColInd)*/
    solver_.Clear();

    mat.copy_crs(stn_, stn_, argValues, argRowStart, argColInd);
    n_ = stn_;
    nnz_ = argValues.size();


    A_.CopyFromHostCSR(
            argRowStart.data(),
            argColInd.data(),
            argValues.data(),
            "SystemMatrix",
             nnz_,
             n_, n_
    );

    A_.MoveToAccelerator();
    solver_.SetOperator(A_);
    precond_.Set(0);
    precond_.SetOperator(A_);
    solver_.SetPreconditioner(precond_);
    solver_.InitMaxIter(1000);
    solver_.InitTol(1e-14, 1e-10, 1e+2); // absolut, relative, divergence
    solver_.SetResidualNorm(2); // L2 norm
    solver_.Build();

}

template <typename TAlgebra>
bool AMD_BiCGStab<TAlgebra>::double_apply(vector_t &c, const vector_t &d){

    if (c.size() != n_ || d.size() != n_) {
        std::cerr << "Size mismatch\n";
        return false;
    }


    const size_t Nd = d.size();
    auto* internal_d = new double[Nd];
    for (size_t i = 0; i < Nd; ++i)
        internal_d[i] = d[i];

    const size_t Nc = d.size();
    auto* internal_c = new double[Nc];
    for (size_t i = 0; i < Nc; ++i)
        internal_c[i] = c[i];

    b_.Allocate("d", n_);
    x_.Allocate("c", n_);
    b_.CopyFromHostData(internal_d);
    x_.CopyFromHostData(internal_c);

    b_.MoveToAccelerator();
    x_.MoveToAccelerator();

    // Solve A * x = b
    solver_.Solve(b_, &x_);

    x_.MoveToHost();
    x_.CopyToHostData(internal_c);

    //double val = 0.0;
    for (size_t i = 0; i < Nc; ++i) {
        c[i] = internal_c[i];
      //  val += internal_c[i];
    }
    //std::cout << "--------------------------------------------- val=" << val << std::endl;

    return true;
}


}