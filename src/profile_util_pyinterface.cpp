/*! \file python_interface.cpp
 *  \brief Define the interface to the profile util classes and calls
*/

//possibly other stuff to include
// cfg['extra_link_args'] = ['...']
// cfg['extra_compile_args'] = ['...']
// cfg['libraries'] = ['...']

<%
cfg['compiler_args'] = ['-std=c++20', '-O3']
cfg['include_dirs'] = ['../include/']
cfg['sources'] = ['profile_util.cpp', 'mem_util.cpp', 'thread_affinity.cpp', 'time_util.cpp']
cfg['dependencies'] = ['profile_util.h', 'profile_util_gpu_api.h']
cfg['parallel'] = True
setup_pybind11(cfg)
%>

#include "profile_util.h"
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace std;
using namespace profiling_util;

/// @brief python module for profile_util
/// @param  profile_util name of python library to import
/// @param  m python module class
PYBIND11_MODULE(profile_util, m) {
    py::class_<Timer>(m, "Timer")
        .def(py::init<const std::string&, const std::string&, const std::string&, bool>())
        .def("get_use_device", &Timer::get_use_device)
        .def("get", &Timer::get)
        .def("get_creation", &Timer::get_creation)
#if defined(_GPU)
        .def("get_on_device", &Timer::get_on_device)
        .def("get_device_swap_info", &Timer::get_device_swap_info)
#endif
        .def("set_ref", &Timer::set_ref)
        .def("get_ref", &Timer::get_ref);
}

