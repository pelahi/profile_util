/*! \file mem_util.cpp
 *  \brief Get memory 
 */

#include "profile_util.h"

/// get the memory use looking as the /proc/self/status file 
namespace profiling_util {

    // return the memory usage;
    memory_usage get_memory_usage() {
        memory_usage usage;

        const char *stat_file = "/proc/self/status";
        std::ifstream f(stat_file);
        if (!f.is_open()) {
            std::cerr << "Couldn't open " << stat_file << " for memory usage reading" <<std::endl;
        }
        for (std::string line; std::getline(f, line); ) {
            auto start = line.find("VmSize:");
            if (start != std::string::npos) {
                std::istringstream is(line.substr(start + 7));
                is >> usage.vm.current;
                continue;
            }
            start = line.find("VmPeak:");
            if (start != std::string::npos) {
                std::istringstream is(line.substr(start + 7));
                is >> usage.vm.peak;
                continue;
            }
            start = line.find("VmRSS:");
            if (start != std::string::npos) {
                std::istringstream is(line.substr(start + 6));
                is >> usage.rss.current;
                continue;
            }
            start = line.find("VmHWM:");
            if (start != std::string::npos) {
                std::istringstream is(line.substr(start + 6));
                is >> usage.rss.peak;
                continue;
            }
        }

        // all values above are in kB
        usage.vm.current *= 1024;
        usage.vm.peak *= 1024;
        usage.rss.current *= 1024;
        usage.rss.peak *= 1024;
        return usage;
    }

    std::string ReportMemUsage(
        const std::string &function, 
        const std::string &line_num
        )
    {
        std::string report;
        memory_usage mem;
        std::tie(report, mem) = GetMemUsage(function, line_num);
        return report;
    }

    //report usage along with change relative to another sampling of memory
    std::string ReportMemUsage(
        const memory_usage &prior_mem_usage,
        const std::string &function, 
        const std::string &line_num
        )
    {
        std::string report;
        memory_usage mem;
        std::tie(report, mem) = GetMemUsage(prior_mem_usage, function, line_num);
        return report;
    }

    std::tuple<std::string, memory_usage> GetMemUsage(
        const std::string &function, 
        const std::string &line_num
        )
    {
        auto memory_usage = get_memory_usage();
        std::ostringstream memory_report;
        auto append_memory_stats = [&memory_report](const char *name, const memory_stats &stats) {
            memory_report << name << " current/peak: " << memory_amount(stats.current) << " / " << memory_amount(stats.peak);
        };
        memory_report << "Memory report @ " << function << " L"<<line_num <<" : ";
        append_memory_stats("VM", memory_usage.vm);
        memory_report << "; ";
        append_memory_stats("RSS", memory_usage.rss);
        return std::make_tuple(memory_report.str(), memory_usage);
    }

    //report usage along with change relative to another sampling of memory
    std::tuple<std::string, memory_usage> GetMemUsage(
        const memory_usage &prior_mem_usage,
        const std::string &function, 
        const std::string &line_num
        )
    {
        auto memory_usage = get_memory_usage();
        memory_usage.vm.change = memory_usage.vm.current - prior_mem_usage.vm.current;
        memory_usage.rss.change = memory_usage.rss.current - prior_mem_usage.rss.current;
        std::ostringstream memory_report;
        auto append_memory_stats = [&memory_report](const char *name, const memory_stats &stats) {
            memory_report << name << " current/peak/change : " << memory_amount(stats.current) << " / " << memory_amount(stats.peak)<< " / "<< memory_amount(stats.change);
        };
        memory_report << "Memory report @ " << function << " L"<<line_num <<" : ";
        append_memory_stats("VM", memory_usage.vm);
        memory_report << "; ";
        append_memory_stats("RSS", memory_usage.rss);
        return std::make_tuple(memory_report.str(), memory_usage);
    }

} 

