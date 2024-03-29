/*! \file time_util.cpp
 *  \brief Get timing
 */

#include "profile_util.h"

/// get the time taken to do some comptue 
namespace profiling_util {
    template <typename T>
    inline
    std::basic_ostream<T> &operator<<(std::basic_ostream<T> &os, const Timer &t) {
        os << ns_time(t.get());
        return os;
    }

    std::string ReportTimeTaken(
        Timer &t, 
        const std::string &function, 
        const std::string &line_num)
    {
        std::string new_ref = "@"+function+" L"+line_num;
        std::ostringstream report;
        report <<"Time taken between : " << new_ref << " - " << t.get_ref() << " : " << ns_time(t.get());
        return report.str();
    }

    float GetTimeTaken(
        Timer &t, 
        const std::string &function, 
        const std::string &line_num)
    {
        return static_cast<float>((t.get()));
    }

#if defined(_GPU)
    std::string ReportTimeTakenOnDevice(
        Timer &t, 
        const std::string &function, 
        const std::string &line_num)
    {
        std::string new_ref = "@"+function+" L"+line_num;
        std::ostringstream report;
        if (t.get_use_device()) report << "Time taken on device between : " ;
        else report << "NO DEVICE to measure : ";
        auto ns = ns_time(t.get_on_device());
        report << t.get_device_swap_info();
        report << new_ref << " - " << t.get_ref() << " : " << ns;
        return report.str();
    }

    float GetTimeTakenOnDevice(
        Timer &t, 
        const std::string &function, 
        const std::string &line_num)
    {
        return t.get_on_device();
    }
#endif

} 

