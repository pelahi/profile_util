/*! \file profile_util_pyinterface.cpp
 *  \brief Define the interface to the profile util classes and calls
*/

//possibly other stuff to include
// cfg['extra_link_args'] = ['...']
// cfg['extra_compile_args'] = ['...']
// cfg['libraries'] = ['...']
#ifndef _USING_SETUP
<%
import subprocess, shutil
def my_preprocess_script():
    # Example sed command

    shutil.copyfile("get_revision.cpp.in", "pybind11_git_revision.cpp")
    #sed_command = ['sed', '-i', '', 's//new_text/g', "pybind11_git_revision.cpp"]
    #subprocess.run(sed_command, check=True)

cfg['preprocess_script'] = my_preprocess_script
cfg['compiler_args'] = ['-std=c++20', '-O3']
cfg['include_dirs'] = ['../include/']
cfg['sources'] = ['profile_util.cpp', 'mem_util.cpp', 'thread_affinity_util.cpp', 'time_util.cpp',  'pybind11_git_revision.cpp']
cfg['dependencies'] = ['../include/profile_util.h', '../include/profile_util_gpu.h']
cfg['parallel'] = True
setup_pybind11(cfg)
%>
#endif

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

    /// @defgroup Thread_affinity
        //@{
    m.def("cpuset_to_cstr", &cpuset_to_cstr);
    m.def("MPICallingRank", &MPICallingRank);
    m.def("ReportParallelAPI", &ReportParallelAPI);
    m.def("ReportBinding", &ReportBinding);
    m.def("ReportThreadAffinity", &ReportThreadAffinity);
#ifdef _MPI
    m.def("MPIReportThreadAffinity", &MPIReportThreadAffinity);
#endif
    //@}
}

