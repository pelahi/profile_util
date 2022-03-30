/*! \file proto.h
 *  \brief this file contains all function prototypes of the code
 */

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

namespace profiling_util {

    void cpuset_to_cstr(cpu_set_t *mask, char *str);
    std::string ReportThreadAffinity();

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
