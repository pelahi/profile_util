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

    // execute a command on the command line
    std::string exec_sys_cmd(std::string cmd) 
    {
        std::string result;
        if (cmd.size()==0) return result;
        std::array<char, 128> buffer;
        std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd.c_str(), "r"), pclose);
        if (!pipe) {
            throw std::runtime_error("popen() failed!");
        }
        while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
            result += buffer.data();
        }
        return result;
    }

    sys_memory_stats get_system_memory()
    {
        auto text = exec_sys_cmd("free | head -n 2 | tail -n 1");
        std::size_t doubleSpace = text.find("  ");
        while (doubleSpace != std::string::npos)
        {
            text.erase(doubleSpace, 1);
            doubleSpace = text.find("  ");
        }
        std::vector<std::string> words{};
        size_t pos = 0;
        std::string space_delimiter = " ";
        sys_memory_stats sysmem;
        // order of data should follow the following headers
        // std::vector<std::string> headers = {"total", "used", "free", "shared", "buff/cache", "available"};
        while ((pos = text.find(space_delimiter)) != std::string::npos) {
            if (text.substr(0, pos) != space_delimiter) words.push_back(text.substr(0, pos));
            text.erase(0, pos + space_delimiter.length());
        }
        words.push_back(text.substr(0,text.size()-1));
        // all memory reported back in kilobytes
        sysmem.total = std::stoul(words[1])*1024;
        sysmem.used = std::stoul(words[2])*1024;
        sysmem.free = std::stoul(words[3])*1024;
        sysmem.shared = std::stoul(words[4])*1024;
        sysmem.cache = std::stoul(words[5])*1024;
        sysmem.avail = std::stoul(words[6])*1024;
        return sysmem;
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

    std::string ReportSystemMem(
        const std::string &function, 
        const std::string &line_num
        )
    {
        std::string report;
        sys_memory_stats mem;
        std::tie(report, mem) = GetSystemMem(function, line_num);
        return report;
    }

    //report usage along with change relative to another sampling of memory
    std::string ReportSystemMem(
        const sys_memory_stats &prior_mem_usage,
        const std::string &function, 
        const std::string &line_num
        )
    {
        std::string report;
        sys_memory_stats mem;
        std::tie(report, mem) = GetSystemMem(prior_mem_usage, function, line_num);
        return report;
    }

    std::tuple<std::string, sys_memory_stats> GetSystemMem(
        const std::string &function, 
        const std::string &line_num
        )
    {
        auto sys_mem = get_system_memory();
        std::ostringstream memory_report;
        auto append_memory_stats = [&memory_report](const char *name, const size_t &stat) {
            memory_report << name << ": " << memory_amount(stat);
        };
        memory_report << "Memory report @ " << function << " L"<<line_num <<" : ";
        append_memory_stats("Total ", sys_mem.total);memory_report << "; ";
        append_memory_stats("Used  ", sys_mem.used);memory_report << "; ";
        append_memory_stats("Free  ", sys_mem.free);memory_report << "; ";
        append_memory_stats("Shared", sys_mem.shared);memory_report << "; ";
        append_memory_stats("Cache ", sys_mem.cache);memory_report << "; ";
        append_memory_stats("Avail ", sys_mem.avail);memory_report << "; ";
        return std::make_tuple(memory_report.str(), sys_mem);
    }

    //report usage along with change relative to another sampling of memory
    std::tuple<std::string, sys_memory_stats> GetSystemMem(
        const sys_memory_stats &prior_mem_usage,
        const std::string &function, 
        const std::string &line_num
        )
    {
        auto sys_mem = get_system_memory();
        std::ostringstream memory_report;
        auto append_memory_stats = [&memory_report](const char *name, const size_t stat, const size_t diff) {
            memory_report << name << ": current/change" << memory_amount(stat)<<" / "<< memory_amount(diff);
        };
        memory_report << "Memory report @ " << function << " L"<<line_num <<" : ";
        append_memory_stats("Total ", sys_mem.total, sys_mem.total-prior_mem_usage.total);memory_report << "; ";
        append_memory_stats("Used  ", sys_mem.used, sys_mem.used-prior_mem_usage.used);memory_report << "; ";
        append_memory_stats("Free  ", sys_mem.free, sys_mem.free-prior_mem_usage.free);memory_report << "; ";
        append_memory_stats("Shared", sys_mem.shared, sys_mem.shared-prior_mem_usage.shared);memory_report << "; ";
        append_memory_stats("Cache ", sys_mem.cache, sys_mem.cache-prior_mem_usage.cache);memory_report << "; ";
        append_memory_stats("Avail ", sys_mem.avail, sys_mem.avail-prior_mem_usage.avail);memory_report << "; ";
        return std::make_tuple(memory_report.str(), sys_mem);
    }
} 

