#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
#include <algorithm>
#include <random>
#include <tuple>
#include <thread>
#include <profile_util.h>


#ifdef USEOPENMP
#include <omp.h>
#endif

#include <mpi.h>

// if want to try running code but not do any actual communication
// #define TURNOFFMPI

int ThisTask, NProcs;

#define LogMPITest() if (ThisTask==0) std::cout<<" running "<<mpifunc<< " test"<<std::endl;
#define LogMPIBroadcaster() if (ThisTask == itask) std::cout<<ThisTask<<" running "<<mpifunc<<" broadcasting "<<sendsize<<" GB"<<std::endl;
#define LogMPISender() if (ThisTask == itask) std::cout<<ThisTask<<" running "<<mpifunc<<" sending "<<sendsize<<" GB"<<std::endl;
#define LogMPIReceiver() if (ThisTask == itask) std::cout<<ThisTask<<" running "<<mpifunc<<std::endl;
#define LogMPIAllComm() if (ThisTask == 0) std::cout<<ThisTask<<" running "<<mpifunc<<" all "<<sendsize<<" GB"<<std::endl;

struct Options
{
    // what times of communication  to do
    bool igather = true, ireduce = true, iscatter = true;
    bool ibcast = false, isendrecv = true;
    // max message size in GB
    double maxgb = 1.0;
};

std::tuple<int,
    std::vector<MPI_Comm> ,
    std::vector<std::string> ,
    std::vector<int>, 
    std::vector<int>, 
    std::vector<int>>
    MPIAllocateComms()
{
    // number of comms is 2, 4, 8, ... till MPI_COMM_WORLD;
    int numcoms = std::floor(log(static_cast<double>(NProcs))/log(2.0))+1;
    int numcomsold = numcoms;
    std::vector<MPI_Comm> mpi_comms(numcoms);
    std::vector<std::string> mpi_comms_name(numcoms);
    std::vector<int> ThisLocalTask(numcoms), NProcsLocal(numcoms), NLocalComms(numcoms);

    mpi_comms[0] = MPI_COMM_WORLD;
    ThisLocalTask[0] = ThisTask;
    NProcsLocal[0] = NProcs;
    NLocalComms[0] = 1;
    mpi_comms_name[0] = "world";
    for (auto i=1;i<=numcomsold;i++) 
    {
        NLocalComms[i] = NProcs/pow(2,i);
        if (NLocalComms[i] < 2) 
        {
            numcoms = i-1;
            break;
        }
        auto ThisLocalCommFlag = ThisTask % NLocalComms[i];
        mpi_comms_name[i] = std::to_string(pow(2,i));
        MPI_Comm_split(MPI_COMM_WORLD, ThisLocalCommFlag, ThisTask, &mpi_comms[i]);
        MPI_Comm_rank(mpi_comms[i], &ThisLocalTask[i]);
        MPI_Comm_size(mpi_comms[i], &NProcsLocal[i]);
    }
    ThisLocalTask.resize(numcoms);
    NLocalComms.resize(numcoms);
    NProcsLocal.resize(numcoms);
    mpi_comms.resize(numcoms);
    mpi_comms_name.resize(numcoms);
    for (auto i=0;i<=numcoms;i++) 
    {
        if (ThisTask==0) std::cout<<" MPI communicator "<<mpi_comms_name[i]<<" has size of "<<NProcsLocal[i]<<" and there are "<<NLocalComms[i]<<" communicators"<<std::endl;
    }
    MPI_Barrier(mpi_comms[0]);
    return std::make_tuple(numcoms,
        std::move(mpi_comms),
        std::move(mpi_comms_name), 
        std::move(ThisLocalTask), 
        std::move(NProcsLocal), 
        std::move(NLocalComms)
        );
}
void MPIFreeComms(std::vector<MPI_Comm> &mpi_comms, std::vector<std::string> &mpi_comms_name){
    for (auto i=1;i<=mpi_comms.size();i++) {
        if (ThisTask==0) std::cout<<"Freeing "<<mpi_comms_name[i]<<std::endl;
        MPI_Comm_free(&mpi_comms[i]);
    }
}

std::vector<unsigned long long> MPISetSize(double maxgb) 
{
    std::vector<unsigned long long> sizeofsends(4);
    sizeofsends[0] = 1024.0*1024.0*1024.0*maxgb/8;
    for (auto i=1;i<sizeofsends.size();i++) sizeofsends[i] = sizeofsends[i-1]/8;
    std::sort(sizeofsends.begin(),sizeofsends.end());
    if (ThisTask==0) for (auto &x: sizeofsends) std::cout<<"Messages of "<<x*sizeof(double)/1024./1024./1024.<<" GB"<<std::endl;
    MPI_Barrier(MPI_COMM_WORLD);
    return sizeofsends;
}

void MPIReportTimeStats(profiling_util::Timer time1, std::string f, std::string l)
{
    std::vector<float> times(NProcs);
    auto p = times.data();
    auto time_taken = profiling_util::GetTimeTaken(time1, f, l);
    MPI_Gather(&time_taken, 1, MPI_FLOAT, p, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
    if (ThisTask == 0) {
        auto ave = 0.0, std = 0.0;
        for (auto &t:times) 
        {
            ave += t;
            std += t*t;
        }
        ave /= static_cast<float>(NProcs);
        std = sqrt((std - ave*ave*NProcs)/(NProcs-1.0));
        std::cout<<f<<" time taken stats is "<<ave<<" +/- "<<std<<std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
}


void MPITestBcast(Options &opt) 
{
    MPI_Status status;
    std::string mpifunc;
    auto[numcoms, mpi_comms, mpi_comms_name, ThisLocalTask, NProcsLocal, NLocalComms] = MPIAllocateComms();
    std::vector<double> data, allreducesum;
    
    double * p1 = nullptr, *p2 = nullptr;
    auto  sizeofsends = MPISetSize(opt.maxgb);
   
    // run broadcast 
    mpifunc = "Bcast";
    LogMPITest();
    for (auto itask=0;itask<NProcs;itask++) 
    {
        auto time1 = NewTimer();
        for (auto i=0;i<sizeofsends.size();i++) {
            auto time2 = NewTimer();
            auto sendsize = sizeofsends[i]*sizeof(double)/1024./1024./1024.;
            LogMPIBroadcaster();
            data.resize(sizeofsends[i]);
            p1 = data.data();
#ifdef TURNOFFMPI
#else 
            MPI_Bcast(p1, sizeofsends[i], MPI_DOUBLE, itask, mpi_comms[0]);
#endif
            if (ThisTask==itask) LogTimeTaken(time2);
        }
        if (ThisTask==itask) LogTimeTaken(time1);
        MPI_Barrier(mpi_comms[0]);
    }
    p1 = p2 = nullptr;
    data.clear();
    data.shrink_to_fit();
    MPIFreeComms(mpi_comms, mpi_comms_name);
}

void MPITestSendRecv(Options &opt) 
{

}
void MPITestAllGather(Options &opt) 
{

}

void MPITestAllScatter(Options &opt) 
{

}

void MPITestAllReduce(Options &opt) 
{
    MPI_Status status;
    std::string mpifunc;
    auto[numcoms, mpi_comms, mpi_comms_name, ThisLocalTask, NProcsLocal, NLocalComms] = MPIAllocateComms();
    std::vector<double> data, allreducesum;
    
    double * p1 = nullptr, *p2 = nullptr;
    auto  sizeofsends = MPISetSize(opt.maxgb);

    // now allreduce 
    mpifunc = "allreduce";
    LogMPITest();
    for (auto i=0;i<sizeofsends.size();i++) 
    {
        auto sendsize = sizeofsends[i]*sizeof(double)/1024./1024./1024.;
        LogMPIAllComm();
        data.resize(sizeofsends[i]);
        allreducesum.resize(sizeofsends[i]);
        if (ThisTask==0) LogMemUsage();
        for (auto &d:data) d = pow(2.0,ThisTask);
        p1 = data.data();
        p2 = allreducesum.data();
        auto time1 = NewTimer();
        for (auto j=1;j<=mpi_comms.size();j++) 
        {
            auto time2 = NewTimer();
            if (ThisLocalTask[j] == 0) std::cout<<ThisTask<<" / "<<ThisLocalTask[j]<<" : Communicating using comm "<<mpi_comms_name[j]<<std::endl;
#ifdef TURNOFFMPI
#else 
            MPI_Allreduce(p1, p2, sizeofsends[i], MPI_DOUBLE, MPI_SUM, mpi_comms[j]);
#endif
            MPIReportTimeStats(time2, __func__, std::to_string(__LINE__));
        }
        {
            auto time2 = NewTimer();
            if (ThisLocalTask[0] == 0) std::cout<<ThisTask<<" / "<<ThisLocalTask[0]<<" : Communicating using comm "<<mpi_comms_name[0]<<std::endl;
#ifdef TURNOFFMPI
            std::this_thread::sleep_for(std::chrono::seconds(10));
            for (auto &d:allreducesum) d = pow(2.0,ThisTask);
#else 
            MPI_Allreduce(p1, p2, sizeofsends[i], MPI_DOUBLE, MPI_SUM, mpi_comms[0]);
#endif
            MPIReportTimeStats(time2, __func__, std::to_string(__LINE__));
        }
        if (ThisTask==0) LogTimeTaken(time1);
    }
    data.clear();
    data.shrink_to_fit();
    allreducesum.clear();
    allreducesum.shrink_to_fit();
    MPIFreeComms(mpi_comms, mpi_comms_name);
}

void MPIRunTests(Options &opt)
{
    if (opt.igather) MPITestAllGather(opt);
    if (opt.iscatter) MPITestAllScatter(opt);
    if (opt.ireduce) MPITestAllReduce(opt);
    if (opt.ibcast) MPITestBcast(opt);
    if (opt.isendrecv) MPITestSendRecv(opt);
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &NProcs);
    MPI_Comm_rank(comm, &ThisTask);
    Options opt;

    auto start = std::chrono::system_clock::now();
    std::time_t start_time = std::chrono::system_clock::to_time_t(start);
    if (ThisTask==0) std::cout << "Starting job at " << std::ctime(&start_time);
    if (argc == 2) opt.maxgb = atof(argv[1]);
    
    if (ThisTask==0) LogParallelAPI();
    MPI_Barrier(MPI_COMM_WORLD);
    LogBinding();
    MPI_Barrier(MPI_COMM_WORLD);
    MPIRunTests(opt);

    auto end = std::chrono::system_clock::now();
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);
    if (ThisTask==0) std::cout << "Ending job at " << std::ctime(&end_time);
}
