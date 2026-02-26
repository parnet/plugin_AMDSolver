#pragma once

#include <rocalution.hpp>

namespace ug{

class LocalData {
public:
    using matrix_t = ug::CPUAlgebra::matrix_type;
    using vector_t = ug::CPUAlgebra::vector_type;

    void clear() {
        A_.Clear();
        b_.Clear();
        x_.Clear();
    }

    void init(const matrix_t &mat) {
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
        A_.MoveToAcceleratorAsync();
    }

    void double_pre_apply(vector_t &c, const vector_t &d){
        if (c.size() != n_ || d.size() != n_) {
            std::cerr << "Size mismatch\n";
        }

        const size_t Nd = d.size();
        internal_d = new double[Nd];
        for (size_t i = 0; i < Nd; ++i)
            internal_d[i] = d[i];

        const size_t Nc = d.size();
        internal_c = new double[Nc];
        for (size_t i = 0; i < Nc; ++i)
            internal_c[i] = c[i];

        b_.Allocate("d", n_);
        x_.Allocate("c", n_);
        b_.CopyFromHostData(internal_d);
        x_.CopyFromHostData(internal_c);

        b_.MoveToAcceleratorAsync();
        x_.MoveToAcceleratorAsync();
        b_.Sync();
        x_.Sync();
    }

    void double_post_apply(vector_t &c, const vector_t &d){
        const size_t Nc = c.size();
        x_.MoveToHost();
        x_.CopyToHostData(internal_c);

        //double val = 0.0;
        for (size_t i = 0; i < Nc; ++i) {
            c[i] = internal_c[i];
        }
        delete[] internal_d;
        delete[] internal_c;
    }

    void sync_operator(){
        A_.Sync();
    }


    int n_ = 0;
    size_t stn_   = 0;
    int nnz_ = 0;

    //double* c;
    double* internal_c;
    double* internal_d;

    std::vector<double> argValues;
    std::vector<int> argRowStart;
    std::vector<int> argColInd;

    rocalution::LocalMatrix<double> A_;
    rocalution::LocalVector<double> b_;
    rocalution::LocalVector<double> x_;

};
}