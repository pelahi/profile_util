#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
#include <profile_util.h>

#ifdef USEOPENMP
#include <omp.h>
#endif

int main() {
    LogParallelAPI();
    LogBinding();
  #ifdef USEOPENMP
    int nthreads = omp_get_max_threads();
  #endif

    auto time_mem = NewTimer();
    std::vector<int> x_int, y_int;
    std::vector<float> x_float, y_float;
    std::vector<double> x_double, y_double;
    const unsigned int Nentries = 24.0*1024.0*1024.0*1024.0/8.0/6.0;
    

    std::cout<<"Vectorization test running with "<<Nentries<<" requiring "<<Nentries*2*(sizeof(int)+sizeof(float)+sizeof(double))/1024./1024./1024.<<"GB"<<std::endl;
    x_int.resize(Nentries);
    y_int.resize(Nentries);
    x_float.resize(Nentries);
    y_float.resize(Nentries);
    x_double.resize(Nentries);
    y_double.resize(Nentries);
    LogMemUsage();
    LogTimeTaken(time_mem);

    auto time_sillycalc = NewTimer();

    #ifdef USEOPENMP
    #pragma omp parallel for \
    default(none) \
    shared(x_int, y_int, x_float, y_float, x_double, y_double, Nentries) \
    schedule(static) if (nthreads > 1)
    #endif
    for (auto i=0u; i<Nentries; i++) {
      x_int[i] = i;
      auto temp = x_int[i];
      y_int[i] = temp+temp*pow(temp,2) + temp/(temp+1);
    }

    #ifdef USEOPENMP
    #pragma omp parallel for \
    default(none) shared(x_float, y_float, Nentries) \
    schedule(static) if (nthreads > 1)
    #endif
    for (auto i=0u; i<Nentries; i++) {
      x_float[i] = i;
      auto tempf = x_float[i];
      y_float[i] = tempf+tempf*pow(tempf,2) + tempf/(tempf+1);
    }

    #ifdef USEOPENMP
    #pragma omp parallel for \
    default(none) shared(x_double, y_double, Nentries) \
    schedule(static) if (nthreads > 1)
    #endif
    for (auto i=0u; i<Nentries; i++) {
      x_double[i] = i;
      auto tempd = x_double[i];
      y_double[i] = tempd+tempd*pow(tempd,2) + tempd/(tempd+1);
    }
    LogTimeTaken(time_sillycalc);
    time_mem = NewTimer();
    x_int.clear();
    x_float.clear();
    x_double.clear();
    y_int.clear();
    y_float.clear();
    y_double.clear();
    x_int.shrink_to_fit();
    x_float.shrink_to_fit();
    x_double.shrink_to_fit();
    y_int.shrink_to_fit();
    y_float.shrink_to_fit();
    y_double.shrink_to_fit();
    LogMemUsage();
    LogTimeTaken(time_mem);
}
