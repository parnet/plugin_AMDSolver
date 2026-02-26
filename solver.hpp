#pragma once


namespace ug {

class AMD_Solver {
public:
    using amd_solver = rocalution::Solver<
    rocalution::LocalMatrix<double>,
    rocalution::LocalVector<double>,
    double
>;

    virtual ~AMD_Solver() = default;

    [[nodiscard]] virtual const char * double_name() const = 0;
    [[nodiscard]] virtual amd_solver& get_solver() = 0;

};
}