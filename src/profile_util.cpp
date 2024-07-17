/*! \file time_util.cpp
 *  \brief Get timing
 */

#include "profile_util.h"

namespace profiling_util {
    std::string __when(){
        auto log_time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
        std::string whenbuff=std::ctime(&log_time);
        whenbuff.erase(std::find(whenbuff.begin(), whenbuff.end(), '\n'), whenbuff.end());
        return whenbuff;
    }
    /*
    void _place_cmd(const std::string cmd)
    {
        auto status = std::system(cmd.c_str());
        // need to add a throw/catch
    };
    void _place_long_lived_cmd(const std::string cmd, const float t)
    {
        while () 
        {
            auto status = std::system(cmd.c_str());
            sleep(t);
        }
    };
    */
    template std::tuple<double, double, double, double>get_stats(std::vector<double> input);
    template std::tuple<float, float, float, float>get_stats(std::vector<float> input);

    #ifdef _MPI
    static bool _PU_USING_MPI=true;
    #endif
    #ifdef _HIP
    static bool _PU_USING_HIP=true;
    #endif
    #ifdef _CUDA
    static bool _PU_USING_CUDA=true;
    #endif
    #ifdef _OPENMP
    static bool _PU_USING_OPENMP=true;
    #endif
}

