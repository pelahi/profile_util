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

    void profiling_util::StateSampler::_launch()
    {
        std::string cmd;
        std::string s_cpu;
        s_cpu = " ps -p " + std::to_string(pid) + " -o %cpu | tail -n 1";
        
        cmd =  profiling_util::StateSampler::_set_sampling(s_cpu, cpu_usage_fname); 
        (*threads).emplace_back(std::thread(&profiling_util::StateSampler::_place_long_lived_cmd, this, cmd, sample_time));
#ifdef _GPU
        if (use_device) 
        {
            // usage cmd 
            std::string s_gpu;
            s_gpu = std::string(pu_gpuMonitorCmd) + std::string(" ") + std::string(pu_gpu_usage_request) + std::string(" ") + std::string(pu_gpu_formating);
            cmd = _set_sampling(s_gpu, gpu_usage_fname);
            // run the gpu commands as a specific process
            // (*threads).emplace_back(std::thread(_place_cmd, cmd));
            (*threads).emplace_back(std::thread(&profiling_util::StateSampler::_place_long_lived_cmd, this, cmd, sample_time));

            // energy cmd 
            s_gpu = std::string(pu_gpuMonitorCmd) + " " + std::string(pu_gpu_energy_request) + " " + std::string(pu_gpu_formating);
            cmd = _set_sampling(s_gpu, gpu_energy_fname);
            // run the gpu commands as a specific process
            (*threads).emplace_back(std::thread(&profiling_util::StateSampler::_place_long_lived_cmd, this, cmd, sample_time));
        }
#endif
    }
    // profiling_util::StateSampler::StateSampler(const std::string &f, const std::string &l, float samples_per_sec, bool _use_device) : profiling_util::Timer(f,l,_use_device), threads(0)
    profiling_util::StateSampler::StateSampler(const std::string &f, const std::string &l, float _sample_time_in_sec, bool _use_device) : profiling_util::Timer::Timer(f,l,_use_device)
    {
        pid = getpid();
        // set the random seed based on curent time 
        std::srand(static_cast<unsigned>(std::time(nullptr)));
        id = std::rand();
        sample_time = _sample_time_in_sec*1000.0;//convert to micro seconds
        use_device = _use_device;

        cpu_usage_fname = ".sampler.cpu_usage." + std::to_string(id) + ".txt";
#ifdef _GPU
        gpu_usage_fname = ".sampler.gpu_usage."+std::to_string(id)+".txt";
        gpu_energy_fname = ".sampler.gpu_energy."+std::to_string(id)+".txt";
#endif
        // allocate the thread vector
        threads = new std::vector<std::thread>;
        _launch();

    }
    profiling_util::StateSampler::~StateSampler()
    {
        Pause();
        delete threads;
        // and remove files
        std::filesystem::remove(cpu_usage_fname);
        std::filesystem::remove(cpu_energy_fname);
#ifdef _GPU
        std::filesystem::remove(gpu_usage_fname);
        std::filesystem::remove(gpu_energy_fname);
#endif

    }
    void profiling_util::StateSampler::Pause()
    {
        std::unique_lock<std::mutex> lock(mtx);
        stopFlag = true;
        cv.notify_all();
        for (auto &t: *threads) t.join();
        (*threads).clear();
    }
    void profiling_util::StateSampler::Restart()
    {
        stopFlag = false;
        cv.notify_all();
        _launch();
    }
    std::string ReportCPUUsage(profiling_util::StateSampler &s, 
        const std::string &function, 
        const std::string &line_num)
    {
        // std::cout<< " going to try getting some stuff without pausing "<<std::endl;
        s.Pause();
        std::ifstream file(s.GetCPUUsageFname());
        std::vector<double> content;
        std::string line;
        while (std::getline(file, line)) {
            content.push_back(std::stof(line));
        }
        auto [ave, std, max, min] = get_stats(content);
        std::string new_ref = "@"+function+" L"+line_num;
        std::ostringstream report;
        report <<"CPU usage statistics taken between : " << new_ref << " - " << s.get_ref() << " : ";
        report <<" [ave,std,min,max] = [ "<<ave<<", "<<std<<", "<<min<<", "<<max<<" ]";
        s.Restart();
        return report.str();
    }
#ifdef _GPU
    std::string ReportGPUUsage(profiling_util::StateSampler &s, 
        const std::string &function, 
        const std::string &line_num, 
        int gid)
    {
        // std::cout<< " going to try getting some stuff without pausing "<<std::endl;
        s.Pause();
        std::ifstream file(s.GetGPUUsageFname());
        std::vector<double> content;
        std::string line;
        while (std::getline(file, line)) {
            content.push_back(std::stof(line));
        }
        auto [ave, std, max, min] = get_stats(content);
        std::string new_ref = "@"+function+" L"+line_num;
        std::ostringstream report;
        report <<"GPU usage statistics taken between : " << new_ref << " - " << s.get_ref() << " : ";
        report <<" [ave,std,min,max] = [ "<<ave<<", "<<std<<", "<<min<<", "<<max<<" ]";
        s.Restart();
        return report.str();
    }
    std::string ReportGPUEnergy(profiling_util::StateSampler &s, 
        const std::string &function, 
        const std::string &line_num,
        int gid)
    {
        // std::cout<< " going to try getting some stuff without pausing "<<std::endl;
        s.Pause();
        std::ifstream file(s.GetGPUEnergyFname());
        std::vector<double> content;
        std::string line;
        while (std::getline(file, line)) {
            content.push_back(std::stof(line));
        }
        auto [ave, std, max, min] = get_stats(content);
        std::string new_ref = "@"+function+" L"+line_num;
        std::ostringstream report;
        report <<"GPU energy statistics taken between : " << new_ref << " - " << s.get_ref() << " : ";
        report <<" [ave,std,min,max] = [ "<<ave<<", "<<std<<", "<<min<<", "<<max<<" ]";
        s.Restart();
        return report.str();
    }
#endif

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

