#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
#include <algorithm>
#include <random>
#include <tuple>
#include <profile_util.h>

#ifdef USEOPENMP
#include <omp.h>
#endif

#include <mpi.h>

int ThisTask, NProcs;

template<class T> std::vector<T> allocate_and_init_vector(unsigned long long N)
{
    std::vector<T> v(N);
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 gen(rd()); // seed the generator
    std::normal_distribution<> distr(0,1);
    auto t_init = NewTimer();
    #pragma omp parallel \
    default(none) shared(v,N) firstprivate(gen,distr) LOGGING() \
    if (N>10000)
    {
        #pragma omp critical
        LogThreadAffinity();
        #pragma omp for 
        for (auto &x:v)
        {
            x = distr(gen);
        }
    }
    LogTimeTaken(t_init);
    LogMemUsage();
    return v;
}


template<class T> T vector_sq_and_sum(std::vector<T> &v)
{
    T sum = 0;
    auto t1 = NewTimer();
    #pragma omp parallel \
    default(none) shared(v,sum) LOGGING() \
    if (v.size()>1000)
    {
        #pragma omp critical
        LogThreadAffinity();
        #pragma omp for reduction(+:sum) nowait 
        for (auto &x:v) 
        {
            //x = x*x;
            sum += x*x;
        }
    }
    LogTimeTaken(t1);
    t1 = NewTimer();
    T sum_serial = 0;
    for (auto &x:v) 
    {
        sum_serial += x*x;
    }
    LogTimeTaken(t1);
    return sum;
}

template<class T> void recursive_vector(std::vector<T> &v)
{
    T sum = 0;
    if (v.size() < 1000000) {
        sum = vector_sq_and_sum(v);
    }
    else 
    {
        std::random_device rd; // obtain a random number from hardware
        std::mt19937 gen(rd()); // seed the generator
        std::uniform_int_distribution<> distr(0,v.size()); // define the range
        auto split = distr(gen);
        std::vector<T> left, right;
        // std::cout<<"@"<<__func__<<" L"<<__LINE__<<" spliting vector of size "<<v.size()<<" at "<<split<<" at level "<<omp_get_level()<<std::endl;
        #pragma omp parallel \
        default(none) shared(v, split, left, right) LOGGING() 
        {
            #pragma omp critical
            LogThreadAffinity();
            #pragma omp single
            {
                // #pragma taskgroup 
                {
                    #pragma task 
                    {
                        #pragma omp critical
                        LogThreadAffinity();
                        std::string s = "Left split of " + std::to_string(split);
                        printf("%s\n",s.c_str());
                        std::copy(v.begin(), v.begin() + split, std::back_inserter(left));
                        recursive_vector(left);
                    }
                    #pragma task 
                    {
                        std::copy(v.begin() + split, v.end(), std::back_inserter(right));
                        recursive_vector(right);
                    }
                    #pragma taskwait 
                }
            }
            #pragma omp for 
            for (auto i=0;i<split;i++) v[i]=left[i];

            #pragma omp for 
            for (auto i=split;i<v.size();i++) v[i]=right[i-split];
        }
    }
}

// allocate mem and logs time taken and memory usage
void allocate_mem(unsigned long long Nentries, std::vector<int> &x_int, std::vector<int> &y_int, 
    std::vector<float> &x_float, std::vector<float> &y_float, 
    std::vector<double> &x_double, std::vector<double> &y_double) 
{
    auto time_mem = NewTimer();
    Log()<<"Vectorization test running with "<<Nentries<<" requiring "<<Nentries*2*(sizeof(int)+sizeof(float)+sizeof(double))/1024./1024./1024.<<"GB"<<std::endl;
    x_int.resize(Nentries);
    y_int.resize(Nentries);
    x_float.resize(Nentries);
    y_float.resize(Nentries);
    x_double.resize(Nentries);
    y_double.resize(Nentries);
    LogMemUsage();
    LogTimeTaken(time_mem);
}

void deallocate_mem( std::vector<int> &x_int, std::vector<int> &y_int, 
    std::vector<float> &x_float, std::vector<float> &y_float, 
    std::vector<double> &x_double, std::vector<double> &y_double) 
{
    auto time_mem = NewTimer();
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

// tests openmp vectorization, logs time taken
void vector_vectorization_test(unsigned long long Nentries, 
    std::vector<int> &x_int, std::vector<int> &y_int, 
    std::vector<float> &x_float, std::vector<float> &y_float, 
    std::vector<double> &x_double, std::vector<double> &y_double) 
{
    #ifdef USEOPENMP
    int nthreads = omp_get_max_threads();
    #endif
    auto time_sillycalc = NewTimer();
    #ifdef USEOPENMP
    #pragma omp parallel for \
    default(none) \
    shared(x_int, y_int, x_float, y_float, x_double, y_double, Nentries, std::cout) \
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
}

#define LogMPITest() if (ThisTask==0) std::cout<<" running "<<mpifunc<< " test"<<std::endl;
#define LogMPIBroadcaster() if (ThisTask == itask) std::cout<<ThisTask<<" running "<<mpifunc<<" broadcasting "<<sendsize<<" GB"<<std::endl;
#define LogMPISender() if (ThisTask == itask) std::cout<<ThisTask<<" running "<<mpifunc<<" sending "<<sendsize<<" GB"<<std::endl;
#define LogMPIReceiver() if (ThisTask == itask) std::cout<<ThisTask<<" running "<<mpifunc<<std::endl;
#define LogMPIAllComm() if (ThisTask == 0) std::cout<<ThisTask<<" running "<<mpifunc<<" all "<<sendsize<<" GB"<<std::endl;


int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &NProcs);
    MPI_Comm_rank(comm, &ThisTask);
    MPISetLoggingComm(comm);

    auto start = std::chrono::system_clock::now();
    std::time_t start_time = std::chrono::system_clock::to_time_t(start);
    if (ThisTask==0) std::cout << "Starting job at " << std::ctime(&start_time);
    auto maxgb=0.8;
    if (argc == 2) maxgb = atof(argv[1]);
    
    LogParallelAPI();
    LogBinding();

    std::vector<int> x_int, y_int;
    std::vector<float> x_float, y_float;
    std::vector<double> x_double, y_double;
    const unsigned long long Nentries = 24.0*1024.0*1024.0*1024.0/8.0/6.0/NProcs/16.0;
    //allocate, test vectorization and deallocate
    //functions showcase use of logging time taken and mem usage
    allocate_mem(Nentries, x_int, y_int, x_float, y_float, x_double, y_double);
    vector_vectorization_test(Nentries, x_int, y_int, x_float, y_float, x_double, y_double);
    deallocate_mem(x_int, y_int, x_float, y_float, x_double, y_double);

    //allocate mem and init vector using random numbers 
    unsigned long long N=100000000;
    x_double = allocate_and_init_vector<double>(N);
    //recursive call highlights thread affinity reporting
    auto t1 = NewTimer();
    recursive_vector(x_double);
    LogTimeTaken(t1);

    auto end = std::chrono::system_clock::now();
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);
    if (ThisTask==0) std::cout << "Ending job at " << std::ctime(&end_time);
}
