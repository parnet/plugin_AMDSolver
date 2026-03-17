#pragma once


#include <rocalution/rocalution.hpp>

#include "lib_algebra/operator/linear_solver/external_solvers/external_solvers.h"

#include "../local_data.hpp"


namespace ug {

template <typename TAlgebra>
class AMD_QR : public IExternalSolver <TAlgebra>,
    public AMD_Solver {
public:
    using matrix_t = CPUAlgebra::matrix_type;
    using vector_t = CPUAlgebra::vector_type;

    AMD_QR() = default;

    ~AMD_QR() override {
        clear();
        };

    amd_solver& get_solver() override {
        return solver_;
    };

    [[nodiscard]] const char * double_name() const override {  return "AMD_QR"; };

    void double_init(const matrix_t &mat) override ;

    bool double_apply(vector_t &c, const vector_t &d) override ;
private:

    void clear() {
        data.clear();
        solver_.Clear();
    }
private:
    LocalData data;

    rocalution::QR<
        rocalution::LocalMatrix<double>,
        rocalution::LocalVector<double>,
        double
    > solver_;

};
// todo set conv_check, abs, rel, divergence, iteration, norm(1,2,3)
template <typename TAlgebra>
void AMD_QR<TAlgebra>::double_init(const matrix_t &mat){
    solver_.Clear();

    data.init(mat);

    solver_.SetOperator(data.A_);

    data.sync_operator();

    solver_.Build();
}

template <typename TAlgebra>
bool AMD_QR<TAlgebra>::double_apply(vector_t &c, const vector_t &d){

    data.double_pre_apply(c,d);
    // Solve A * x = b
    solver_.Solve(data.b_, &data.x_);

    data.double_post_apply(c,d);
    //std::cout << "--------------------------------------------- val=" << val << std::endl;

    return true;
}


}