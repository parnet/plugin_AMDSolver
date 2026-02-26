#pragma once


#include <rocalution.hpp>

#include <vector>

namespace ug{

class AMD_Preconditioner : public AMD_Solver {
public:
    AMD_Preconditioner() = default;

    ~AMD_Preconditioner() override = default;



    using amd_preconditioner = rocalution::Preconditioner<
        rocalution::LocalMatrix<double>,
        rocalution::LocalVector<double>,
        double
    >;

    virtual amd_preconditioner& get_preconditioner() = 0;
};



class AMD_Jacobi : public AMD_Preconditioner {
public:
    using matrix_t = CPUAlgebra::matrix_type;
    using vector_t = CPUAlgebra::vector_type;

    AMD_Jacobi() = default;

    ~AMD_Jacobi() override = default;

    [[nodiscard]] const char * double_name() const override {  return "AMD_Jacobi"; };

    rocalution::Jacobi<
        rocalution::LocalMatrix<double>,
        rocalution::LocalVector<double>,
        double
    > preconditioner_;

    amd_preconditioner& get_preconditioner() override {
        return preconditioner_;
    };

    amd_solver& get_solver() override {
        return preconditioner_;
    };
};




class AMD_GaussSeidel : AMD_Preconditioner {
public:
    using matrix_t = CPUAlgebra::matrix_type;
    using vector_t = CPUAlgebra::vector_type;

    AMD_GaussSeidel() = default ;

    ~AMD_GaussSeidel() override = default;

    const char * double_name() const override {  return "AMD_GaussSeidel"; };

    rocalution::GS<
        rocalution::LocalMatrix<double>,
        rocalution::LocalVector<double>,
        double
    > preconditioner_;

    amd_preconditioner& get_preconditioner() override {
        return preconditioner_;
    };

    amd_solver& get_solver() override {
        return preconditioner_;
    };
};


class AMD_SymmetricGaussSeidel : AMD_Preconditioner {
public:
    using matrix_t = CPUAlgebra::matrix_type;
    using vector_t = CPUAlgebra::vector_type;

    AMD_SymmetricGaussSeidel() =default;

    ~AMD_SymmetricGaussSeidel() override = default;

    [[nodiscard]] const char * double_name() const override {  return "AMD_SymmetricGaussSeidel"; };

    rocalution::SGS<
        rocalution::LocalMatrix<double>,
        rocalution::LocalVector<double>,
        double
    > preconditioner_;

    amd_preconditioner& get_preconditioner() override {
        return preconditioner_;
    };

    amd_solver& get_solver() override {
        return preconditioner_;
    };
};

class AMD_IncompleteLU : AMD_Preconditioner {
public:
    using matrix_t = CPUAlgebra::matrix_type;
    using vector_t = CPUAlgebra::vector_type;

    AMD_IncompleteLU() = default;

    ~AMD_IncompleteLU() override = default;

    void set(int p, bool level = true){ // based on power, level
        preconditioner_.Set(p,level);
    }

    [[nodiscard]]
    const char * double_name() const override {  return "AMD_IncompleteLU"; };
private:
    rocalution::ILU<
        rocalution::LocalMatrix<double>,
        rocalution::LocalVector<double>,
        double
    > preconditioner_;

    amd_preconditioner& get_preconditioner() override {
        return preconditioner_;
    };

    amd_solver& get_solver() override {
        return preconditioner_;
    };
};


class AMD_ItIncompleteLU0 : AMD_Preconditioner {
public:
    using matrix_t = CPUAlgebra::matrix_type;
    using vector_t = CPUAlgebra::vector_type;

    AMD_ItIncompleteLU0() = default;

    void set_algorithm(int alg){
        if(0 < alg && alg <= 4) {
           preconditioner_.SetAlgorithm(static_cast<rocalution::_itilu0_alg>(alg));
        }
    }

    void SetOptions(int option){
        if(0 < option && option < 64) {
            preconditioner_.SetOptions(option);
        }
    }

    void SetMaxIter(int maxiter){
        preconditioner_.SetMaxIter(maxiter);
    }

    void SetTolerance(double tol){
        preconditioner_.SetTolerance(tol);
    }

    ~AMD_ItIncompleteLU0() override = default;


    [[nodiscard]] const char * double_name() const override {  return "AMD_ItIncompleteLU0"; };

    rocalution::ItILU0<
        rocalution::LocalMatrix<double>,
        rocalution::LocalVector<double>,
        double
    > preconditioner_;

    amd_preconditioner& get_preconditioner() override {
        return preconditioner_;
    };

    amd_solver& get_solver() override {
        return preconditioner_;
    };
};


class AMD_ILUT : AMD_Preconditioner {
public:
    using matrix_t = CPUAlgebra::matrix_type;
    using vector_t = CPUAlgebra::vector_type;

    AMD_ILUT() = default;

    void set(double t){
        preconditioner_.Set(t);
    }
    void set(double t, int maxrow) {
        preconditioner_.Set(t,maxrow);
    }

    ~AMD_ILUT() override = default;

    [[nodiscard]] const char * double_name() const override {  return "AMD_ILUT"; };

    rocalution::ILUT<
        rocalution::LocalMatrix<double>,
        rocalution::LocalVector<double>,
        double
    > preconditioner_;

    amd_preconditioner& get_preconditioner() override {
        return preconditioner_;
    };

    amd_solver& get_solver() override {
        return preconditioner_;
    };
};


class AMD_IC : AMD_Preconditioner {
public:
    using matrix_t = CPUAlgebra::matrix_type;
    using vector_t = CPUAlgebra::vector_type;

    AMD_IC() = default;

    ~AMD_IC() override = default;


    [[nodiscard]] const char * double_name() const override {  return "AMD_IC"; };

    rocalution::IC<
        rocalution::LocalMatrix<double>,
        rocalution::LocalVector<double>,
        double
    > preconditioner_;

    amd_preconditioner& get_preconditioner() override {
        return preconditioner_;
    };

    amd_solver& get_solver() override {
        return preconditioner_;
    };
};


class AMD_VariablePreconditioner : AMD_Preconditioner {
public:
    using matrix_t = CPUAlgebra::matrix_type;
    using vector_t = CPUAlgebra::vector_type;
    using amd_preconditioner_ptr = rocalution::Preconditioner<
        rocalution::LocalMatrix<double>,
        rocalution::LocalVector<double>,
        double> * ;

    AMD_VariablePreconditioner() = default;

    void add_preconditioner(SmartPtr<AMD_Preconditioner> precond){
        preconditioner_list.push_back(precond);
    }

    ~AMD_VariablePreconditioner() override {
        if (preconditioner_ptr_list != nullptr) {
            delete[] preconditioner_ptr_list;
        }
    };


    [[nodiscard]] const char * double_name() const override {  return "AMD_VariablePreconditioner"; };

    std::vector<SmartPtr<AMD_Preconditioner>> preconditioner_list;

    rocalution::VariablePreconditioner<
        rocalution::LocalMatrix<double>,
        rocalution::LocalVector<double>,
        double
    > preconditioner_;

     amd_preconditioner_ptr *preconditioner_ptr_list = nullptr;

    amd_preconditioner& get_preconditioner() override {
        if (preconditioner_ptr_list != nullptr) {
            delete[] preconditioner_ptr_list;
        }
        preconditioner_ptr_list = new amd_preconditioner_ptr[preconditioner_list.size()];
        for (size_t i = 0; i < preconditioner_list.size(); ++i) {
            preconditioner_ptr_list[i] = &preconditioner_list[i]->get_preconditioner();
        }
        return preconditioner_;
    };

    amd_solver& get_solver() override {
        return preconditioner_;
    };
};

}