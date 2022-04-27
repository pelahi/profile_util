/*! \file profile_util.h
 *  \brief this file contains all function prototypes of the code
 */

#ifndef _PROFILE_UTIL
#define _PROFILE_UTIL

#include <cstring>
#include <string>
#include <tuple>
#include <ostream>
#include <sstream>
#include <iostream>
#include <chrono>
#include <fstream>
#include <iomanip>

#include <sched.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/sysinfo.h>

#ifdef _MPI
#include <mpi.h>
#endif 

#ifdef _OPENMP 
#include <omp.h>
#endif

namespace profiling_util {

    /// function that converts the mask of thread affinity to human readable string 
    void cpuset_to_cstr(cpu_set_t *mask, char *str);
    /// reports the parallelAPI 
    /// @return string of MPI comm size and OpenMP version and max threads for given rank
    /// \todo needs to be generalized to report parallel API of code and not library
    std::string ReportParallelAPI();
    /// reports binding of MPI comm world and each ranks thread affinity 
    /// @return string of MPI comm rank and thread core affinity 
    std::string ReportBinding();
    /// reports thread affinity within a given scope, thus depends if called within OMP region 
    /// @param func function where called in code, useful to provide __func__ and __LINE
    /// @param line code line number where called
    /// @return string of thread core affinity 
    std::string ReportThreadAffinity(std::string func, std::string line);
#ifdef _MPI
    /// reports thread affinity within a given scope, thus depends if called within OMP region, MPI aware
    /// @param func function where called in code, useful to provide __func__ and __LINE
    /// @param line code line number where called
    /// @param comm MPI communicator
    /// @return string of MPI comm rank and thread core affinity 
    std::string MPIReportThreadAffinity(std::string func, std::string line, MPI_Comm &comm)
#endif

    namespace detail {

        template <int N, typename T>
        struct _fixed {
            T _val;
        };

        template <typename T, int N, typename VT>
        inline
        std::basic_ostream<T> &operator<<(std::basic_ostream<T> &os, detail::_fixed<N, VT> v)
        {
            os << std::setprecision(N) << std::fixed << v._val;
            return os;
        }

    } // namespace detail

    ///
    /// Sent to a stream object, this manipulator will print the given value with a
    /// precision of N decimal places.
    ///
    /// @param v The value to send to the stream
    ///
    template <int N, typename T>
    inline
    detail::_fixed<N, T> fixed(T v) {
        return {v};
    }

    namespace detail {

        struct _memory_amount {
            std::size_t _val;
        };

        struct _microseconds_amount {
            std::chrono::microseconds::rep _val;
        };

        template <typename T>
        inline
        std::basic_ostream<T> &operator<<(std::basic_ostream<T> &os, const detail::_memory_amount &m)
        {

            if (m._val < 1024) {
                os << m._val << " [B]";
                return os;
            }

            float v = m._val / 1024.;
            const char *suffix = " [KiB]";

            if (v > 1024) {
                v /= 1024;
                suffix = " [MiB]";
            }
            if (v > 1024) {
                v /= 1024;
                suffix = " [GiB]";
            }
            if (v > 1024) {
                v /= 1024;
                suffix = " [TiB]";
            }
            if (v > 1024) {
                v /= 1024;
                suffix = " [PiB]";
            }
            if (v > 1024) {
                v /= 1024;
                suffix = " [EiB]";
            }
            // that should be enough...

            os << fixed<3>(v) << suffix;
            return os;
        }

        template <typename T>
        inline
        std::basic_ostream<T> &operator<<(std::basic_ostream<T> &os, const detail::_microseconds_amount &t)
        {
            auto time = t._val;
            if (time < 1000) {
                os << time << " [us]";
                return os;
            }

            time /= 1000;
            if (time < 1000) {
                os << time << " [ms]";
                return os;
            }

            float ftime = time / 1000.f;
            const char *prefix = " [s]";
            if (ftime > 60) {
                ftime /= 60;
                prefix = " [min]";
                if (ftime > 60) {
                    ftime /= 60;
                    prefix = " [h]";
                    if (ftime > 24) {
                        ftime /= 24;
                        prefix = " [d]";
                    }
                }
            }
            // that should be enough...

            os << fixed<3>(ftime) << prefix;
            return os;
        }

    } // namespace detail

    ///
    /// Sent to a stream object, this manipulator will print the given amount of
    /// memory using the correct suffix and 3 decimal places.
    ///
    /// @param v The value to send to the stream
    ///
    inline
    detail::_memory_amount memory_amount(std::size_t amount) {
        return {amount};
    }

    ///
    /// Sent to a stream object, this manipulator will print the given amount of
    /// nanoseconds using the correct suffix and 3 decimal places.
    ///
    /// @param v The value to send to the stream
    ///
    inline
    detail::_microseconds_amount us_time(std::chrono::microseconds::rep amount) {
        return {amount};
    }

    struct memory_stats {
        std::size_t current;
        std::size_t peak;
        std::size_t change;
    };

    struct memory_usage {
        memory_stats vm;
        memory_stats rss;
    };

    ///get memory usage
    memory_usage get_memory_usage();
    ///report memory usage from within a specific function/scope
    ///usage would be from within a function use 
    ///auto l=std::to_string(__LINE__); auto f = __func__; GetMemUsage(f,l);
    std::string ReportMemUsage(const std::string &f, const std::string &l);
    /// like above but also reports change relative to another sampling of memory 
    std::string ReportMemUsage(const memory_usage &prior_mem_use, const std::string &f, const std::string &l);
    /// like ReportMemUsage but also returns the mem usage 
    std::tuple<std::string, memory_usage> GetMemUsage(const std::string &f, const std::string &l);
    std::tuple<std::string, memory_usage> GetMemUsage(const memory_usage &prior_mem_use, const std::string &f, const std::string &l);

    /// Timer class. 
    /// In code create an instance of time and then just a mantter of 
    /// creating an instance and then reporting it. 
    class Timer {

    public:

        using clock = std::chrono::high_resolution_clock;
        using duration = typename std::chrono::microseconds::rep;
        

        /*!
         * Returns the number of milliseconds elapsed since the reference time
         * of the timer
         *
         * @return The time elapsed since the creation of the timer, in [us]
         */
        inline
        duration get() const {
            return std::chrono::duration_cast<std::chrono::microseconds>(clock::now() - tref).count();
        }

        /*!
         * Returns the number of milliseconds elapsed since the creation
         * of the timer
         *
         * @return The time elapsed since the creation of the timer, in [us]
         */
        inline
        duration get_creation() const {
            return std::chrono::duration_cast<std::chrono::microseconds>(clock::now() - t0).count();
        }

        void set_ref(const std::string &new_ref)
        {
            ref = new_ref;
            t0 = clock::now();
        };
        std::string get_ref() const 
        {
            return ref;
        };

        Timer(const std::string &f, const std::string &l) {
            ref="@"+f+" L"+l;
            t0 = clock::now();
            tref = t0;
        }

    private:
        clock::time_point t0;
        clock::time_point tref;
        std::string ref;
    };

    /// get the time taken between some reference time (which defaults to creation of timer )
    /// and current call
    std::string ReportTimeTaken(const Timer &t, const std::string &f, const std::string &l);
    float GetTimeTaken(const Timer &t, const std::string &f, const std::string &l);
}

/// MPI helper routines
//{@
#ifdef _MPI 
//#define MPIOnly0 if (ThisTask == 0)
std::string MPICallingRank(int task);
#endif
//@}
/// \def utility definitions 
//@{
#define _where_calling_from "@"<<__func__<<" L"<<std::to_string(__LINE__)
//@}
/// \defgroup LogAffinity
/// Log thread affinity and parallelism either to std or an ostream
//@{
#define LogParallelAPI() std::cout<<_where_calling_from<<"\n"<<profiling_util::ReportParallelAPI()<<std::endl;
#define LogBinding() std::cout<<_where_calling_from<<"\n"<<profiling_util::ReportBinding()<<std::endl;
#define LogThreadAffinity() printf("%s \n", profiling_util::ReportThreadAffinity(__func__, std::to_string(__LINE__)).c_str());
#define LoggerThreadAffinity(logger) logger<<profiling_util::ReportThreadAffinity(__func__, std::to_string(__LINE__))<<std::endl;
#ifdef _MPI
#define MPILog0ThreadAffinity() if(ThisTask == 0) std::cout<<profiling_util::ReportThreadAffinity(__func__, std::to_string(__LINE__))<<std::endl;
#define MPILogger0ThreadAffinity(logger) if(ThisTask == 0)logger<<profiling_util::ReportThreadAffinity(__func__, std::to_string(__LINE__))<<std::endl;
#define MPILogThreadAffinity(comm) std::cout<<profiling_util::MPIReportThreadAffinity(__func__, std::to_string(__LINE__), comm)<<std::endl;
#define MPILoggerThreadAffinity(logger,comm) logger<<profiling_util::MPIReportThreadAffinity(__func__, std::to_string(__LINE__), comm)<<std::endl;
#define MPILog0ParallelAPI() if(ThisTask==0) std::cout<<_where_calling_from<<"\n"<<profiling_util::ReportParallelAPI()<<std::endl;
#define MPILog0Binding() if (ThisTask == 0) std::cout<<_where_calling_from<<"\n"<<profiling_util::ReportBinding()<<std::endl;
#endif
//@}

/// \defgroup LogMem
/// Log memory usage either to std or an ostream
//@{
#define LogMemUsage() std::cout<<profiling_util::ReportMemUsage(__func__, std::to_string(__LINE__))<<std::endl;
#define LoggerMemUsage(logger) logger<<profiling_util::ReportMemUsage(__func__, std::to_string(__LINE__))<<std::endl;
#ifdef _MPI
#define MPILogMemUsage() std::cout<<profiling_util::MPICallingRank(ThisTask)<<profiling_util::ReportMemUsage(__func__, std::to_string(__LINE__))<<std::endl;
#define MPILoggerMemUsage(logger) logger<<profiling_util::MPICallingRank(ThisTask)<<profiling_util::ReportMemUsage(__func__, std::to_string(__LINE__))<<std::endl;
#endif
//@}


/// \defgroup LogTime
/// Log time taken either to std or an ostream
//@{
#define LogTimeTaken(timer) std::cout<<profiling_util::ReportTimeTaken(timer, __func__, std::to_string(__LINE__))<<std::endl;
#define LoggerTimeTaken(logger,timer) logger<<profiling_util::ReportTimeTaken(timer,__func__, std::to_string(__LINE__))<<std::endl;
#ifdef _MPI
#define MPILogTimeTaken(timer) std::cout<<profiling_util::MPICallingRank(ThisTask)<<profiling_util::ReportTimeTaken(timer, __func__, std::to_string(__LINE__))<<std::endl;
#define MPILoggerTimeTaken(logger,timer) logger<<profiling_util::MPICallingRank(ThisTask)<<profiling_util::ReportTimeTaken(timer,__func__, std::to_string(__LINE__))<<std::endl;
#endif 
#define NewTimer() profiling_util::Timer(__func__, std::to_string(__LINE__));
//@}

/// \defgroup C_naming
/// Extern C interface
//@{
extern "C" {
    /// \defgropu LogAffinity_C
    //@{
    int report_parallel_api(char *str);
    #define log_parallel_api() printf("@%s L%d %s\n", __func__, __LINE__, profiling_util::ReportParallelAPI().c_str());
    int report_binding(char *str);
    #define log_binding() printf("@%s L%d %s\n", __func__, __LINE__, profiling_util::ReportBinding().c_str());
    int report_thread_affinity(char *str, char *f, int l);
    #define log_thread_affinity() printf("%s\n", profiling_util::ReportThreadAffinity(__func__, std::to_string(__LINE__)).c_str());
    #ifdef _MPI
        #define mpi_log0_thread_affinity() if(ThisTask == 0) printf("%s\n", ReportThreadAffinity(__func__, std::to_string(__LINE__)).c_str());
        #define mpi_log_thread_affinity() printf("%s\n", profiling_util::MPIReportThreadAffinity(__func__, std::to_string(__LINE__)).c_str());
        #define mpi_log0_parallel_api() if(ThisTask==0) printf("@%s L%d %s\n", __func__, __LINE__, profiling_util::ReportParallelAPI().c_str());
        #define mpi_log0_binding() if (ThisTask == 0) printf("@%s L%d %s\n", __func__, __LINE__, profiling_util::ReportBinding().c_str());
    #endif
    //@}

    /// \defgroup LogMem_C
    /// Log memory usage either to std or an ostream
    //@{
    #define log_mem_usage() printf("%s \n", profiling_util::ReportMemUsage(__func__, std::to_string(__LINE__)).c_str());
    #ifdef _MPI
    #define mpi_log_mem_usage() printf("%s %s \n", profiling_util::MPICallingRank(ThisTask).c_str(), profiling_util::ReportMemUsage(__func__, std::to_string(__LINE__)).c_str());
    #endif
    //@}

    /// \defgroup LogTime_C
    /// \todo Still implementing a C friendly structure to record timing. 
    //@{
    struct timer_c {
        double t0, tref;
        char ref[2000], where[2000];
        timer_c(){
            t0 = tref = 0;
            memset(ref, 0, sizeof(ref));
            memset(where, 0, sizeof(where));
        }
        void set_timer_ref(double _t0, char _ref[2000]) {
            t0 = _t0;
            strcpy(ref,_ref);
        }
    };
    void report_time_taken(char *str, timer_c &t);

    //@}
}
//@}

#endif
