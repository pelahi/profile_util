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

int Niter=1;
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

    for (auto i=0;i<=numcomsold;i++) 
    {
        NLocalComms[i] = NProcs/pow(2,i+1);
        if (NLocalComms[i] < 2) 
        {
            numcoms = i+1;
            break;
        }
        auto ThisLocalCommFlag = ThisTask % NLocalComms[i];
        MPI_Comm_split(MPI_COMM_WORLD, ThisLocalCommFlag, ThisTask, &mpi_comms[i]);
        MPI_Comm_rank(mpi_comms[i], &ThisLocalTask[i]);
        MPI_Comm_size(mpi_comms[i], &NProcsLocal[i]);
        int tasktag = ThisTask;
        MPI_Bcast(&tasktag, 1, MPI_INTEGER, 0, mpi_comms[i]);
        mpi_comms_name[i] = "Comm tag " + std::to_string(static_cast<int>(pow(2,i+1)))+" world rank " + std::to_string(tasktag);
    }
    mpi_comms[numcoms-1] = MPI_COMM_WORLD;
    ThisLocalTask[numcoms-1] = ThisTask;
    NProcsLocal[numcoms-1] = NProcs;
    NLocalComms[numcoms-1] = 1;
    mpi_comms_name[numcoms-1] = "world";
    ThisLocalTask.resize(numcoms);
    NLocalComms.resize(numcoms);
    NProcsLocal.resize(numcoms);
    mpi_comms.resize(numcoms);
    mpi_comms_name.resize(numcoms);
    for (auto i=0;i<numcoms;i++) 
    {
        if (ThisLocalTask[i]==0) std::cout<<" MPI communicator "<<mpi_comms_name[i]<<" has size of "<<NProcsLocal[i]<<" and there are "<<NLocalComms[i]<<" communicators"<<std::endl;
    }
    MPI_Barrier(mpi_comms[numcoms-1]);
    return std::make_tuple(numcoms,
        std::move(mpi_comms),
        std::move(mpi_comms_name), 
        std::move(ThisLocalTask), 
        std::move(NProcsLocal), 
        std::move(NLocalComms)
        );
}

void MPIFreeComms(std::vector<MPI_Comm> &mpi_comms, std::vector<std::string> &mpi_comms_name){
    for (auto i=0;i<mpi_comms.size()-1;i++) {
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

std::vector<float> MPIGatherTimeStats(profiling_util::Timer time1, std::string f, std::string l)
{
    std::vector<float> times(NProcs);
    auto p = times.data();
    auto time_taken = profiling_util::GetTimeTaken(time1, f, l);
    MPI_Gather(&time_taken, 1, MPI_FLOAT, p, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
    return times;
}

std::tuple<float, float, float, float> TimeStats(std::vector<float> times) 
{
    auto ave = 0.0, std = 0.0;
    auto mint = times[0];
    auto maxt = times[0];
    for (auto &t:times)
    {
        ave += t;
        std += t*t;
        mint = std::min(mint,t);
        maxt = std::max(maxt,t);
    }
    float n = times.size();
    ave /= n;
    if (n>1) std = sqrt((std - ave*ave*n)/(n-1.0));
    else std = 0;
    return std::make_tuple(ave, std, mint, maxt);
}

void MPIReportTimeStats(profiling_util::Timer time1, std::string commname, std::string f, std::string l)
{
    auto times = MPIGatherTimeStats(time1, f, l);
    if (ThisTask == 0) {
        auto[ave, std, mint, maxt] = TimeStats(times);
        std::cout<<"MPI "<<commname<<" comm "<<f<<"@"<<l<<": time taken stats is "<<ave<<" +/- "<<std<<" with [min,max]=["<<mint<<","<<maxt<<"]"<<std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

void MPIReportTimeStats(std::vector<float> times, std::string commname, std::string f, std::string l)
{
    auto[ave, std, mint, maxt] = TimeStats(times);
    if (ThisTask==0) std::cout<<"MPI "<<commname<<" comm "<<f<<"@"<<l<<": time taken stats is "<<ave<<" +/- "<<std<<" with [min,max]=["<<mint<<","<<maxt<<"]"<<std::endl;
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
    MPI_Status status;
    std::string mpifunc;
    auto[numcoms, mpi_comms, mpi_comms_name, ThisLocalTask, NProcsLocal, NLocalComms] = MPIAllocateComms();
    std::vector<double> senddata, receivedata;
    
    double * p1 = nullptr, *p2 = nullptr;
    auto  sizeofsends = MPISetSize(opt.maxgb);

    // now allreduce 
    mpifunc = "sendrecv";
    LogMPITest();
    for (auto i=0;i<sizeofsends.size();i++) 
    {
        auto sendsize = sizeofsends[i]*sizeof(double)/1024./1024./1024.;
        LogMPIAllComm();
        senddata.resize(sizeofsends[i]);
        receivedata.resize(sizeofsends[i]);
        if (ThisTask==0) LogMemUsage();
        for (auto &d:senddata) d = pow(2.0,ThisTask);
        p1 = senddata.data();
        p2 = receivedata.data();
        auto time1 = NewTimer();
        for (auto j=0;j<mpi_comms.size();j++)
        {
#ifdef TURNOFFMPI
#else
            if (ThisLocalTask[j] == 0) std::cout<<ThisTask<<" / "<<ThisLocalTask[j]<<" : Communicating using comm "<<mpi_comms_name[j]<<std::endl;
            std::vector<float> times;
            for (auto iter=0;iter<Niter;iter++) {
                auto time2 = NewTimer();
                std::vector<MPI_Request> sendreqs, recvreqs;
                for (auto isend=0;isend<NProcsLocal[j];isend++) 
                {
                    if (isend != ThisLocalTask[j]) 
                    {
                        MPI_Request request;
                        int tag = isend*NProcsLocal[j]+ThisLocalTask[j];
                        MPI_Isend(p1, sizeofsends[i], MPI_DOUBLE, isend, tag, mpi_comms[j], &request);
                        sendreqs.push_back(request);
                    }
                }
                for (auto irecv=0;irecv<NProcsLocal[j];irecv++) 
                {
                    if (irecv != ThisLocalTask[j]) 
                    {
                        MPI_Request request;
                        int tag = ThisLocalTask[j]*NProcsLocal[j]+irecv;
                        MPI_Irecv(p2, sizeofsends[i], MPI_DOUBLE, irecv, tag, mpi_comms[j], &request);
                        recvreqs.push_back(request);
                    }
                }
                MPI_Waitall(recvreqs.size(), recvreqs.data(), MPI_STATUSES_IGNORE);
                auto times_tmp = MPIGatherTimeStats(time2, __func__, std::to_string(__LINE__));
                times.insert(times.end(), times_tmp.begin(), times_tmp.end());
            }
            MPI_Barrier(MPI_COMM_WORLD);
            MPIReportTimeStats(times, mpi_comms_name[j], __func__, std::to_string(__LINE__));
#endif
        }
        if (ThisTask==0) LogTimeTaken(time1);
    }
    senddata.clear();
    senddata.shrink_to_fit();
    receivedata.clear();
    receivedata.shrink_to_fit();
    MPIFreeComms(mpi_comms, mpi_comms_name);

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
        for (auto j=0;j<mpi_comms.size();j++) 
        {
#ifdef TURNOFFMPI
#else
            if (ThisLocalTask[j] == 0) std::cout<<ThisTask<<" / "<<ThisLocalTask[j]<<" : Communicating using comm "<<mpi_comms_name[j]<<std::endl;
            std::vector<float> times;
            for (auto iter=0;iter<Niter;iter++) {
                auto time2 = NewTimer();
                MPI_Allreduce(p1, p2, sizeofsends[i], MPI_DOUBLE, MPI_SUM, mpi_comms[j]);
                auto times_tmp = MPIGatherTimeStats(time2, __func__, std::to_string(__LINE__));
                times.insert(times.end(), times_tmp.begin(), times_tmp.end());
            }
            MPI_Barrier(MPI_COMM_WORLD);
            MPIReportTimeStats(times, mpi_comms_name[j], __func__, std::to_string(__LINE__));
#endif
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
    if (argc >= 2) opt.maxgb = atof(argv[1]);
    if (argc == 3) Niter = atof(argv[2]);
    
    if (ThisTask==0) LogParallelAPI();
    MPI_Barrier(MPI_COMM_WORLD);
    LogBinding();
    MPI_Barrier(MPI_COMM_WORLD);
    MPIRunTests(opt);

    auto end = std::chrono::system_clock::now();
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);
    if (ThisTask==0) std::cout << "Ending job at " << std::ctime(&end_time);
}
