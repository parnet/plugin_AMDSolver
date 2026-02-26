#pragma once

namespace ug{
class AMD_Rocalution {
public:

    static AMD_Rocalution& instance() {
        static AMD_Rocalution instance;  // Created on first use
        return instance;
    }

    void init(){
        if(!initialized){
            rocalution::set_device_rocalution(0);
            int rank = 0;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            rocalution::init_rocalution(rank);
            std::cout << "get_arch_rocalution: " << rocalution::get_arch_rocalution() << std::endl;
            initialized = true;
        }
    }

    void finish(){
        if(initialized){
            rocalution::stop_rocalution();
            initialized = false;
        }
    }


private:
    AMD_Rocalution() = default;

    ~AMD_Rocalution() {
        if (initialized) {
            rocalution::stop_rocalution();
        }
    }

    bool initialized = false;
};


// lua methods to call

void rocalution_init(){
    AMD_Rocalution::instance().init();
    }

void rocalution_stop(){
    AMD_Rocalution::instance().finish();
    }

}