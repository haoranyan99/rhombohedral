// Util/MPI/mpi.h
#pragma once

#include <iostream>
#include <stdexcept>
#include <string>

#ifdef USE_MPI
  #include <mpi.h>
#endif

namespace rgmpi {


// Return true if MPI has been initialized (USE_MPI enabled)
inline bool inited()
{
#ifdef USE_MPI
    int flag = 0;
    MPI_Initialized(&flag);
    return (flag != 0);
#else
    return false;
#endif
}


inline void barrier()
{
#ifdef USE_MPI
    if (inited()) MPI_Barrier(MPI_COMM_WORLD);
#endif
}

inline void finalize()
{
#ifdef USE_MPI
    MPI_Finalize();
#endif
}

// Safe rank/size getter: if MPI not inited -> (0,1)
inline void rank_size(int& rank, int& nprocs)
{
#ifdef USE_MPI
    if (inited()) {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
        return;
    }
#endif
    rank   = 0;
    nprocs = 1;
}





inline std::pair<int,int> block_1d_int(int n/*task_num*/, int r/*rank*/, int p/*proc_num*/) {
    if (p <= 1) return {0, n};
    const int base = n / p;
    const int rem  = n % p;
    const int my_n = base + (r < rem ? 1 : 0);
    const int my_s = r * base + std::min(r, rem);
    return {my_s/*task_start_idx*/, my_s + my_n/*task_end_idx*/};
}

// Abort all ranks if MPI is inited; otherwise throw.
// NOTE: This function does NOT call MPI_Init for you.
[[noreturn]] inline void abort_all(const std::string& msg, int errcode = 1)
{
#ifdef USE_MPI
    if (inited()) {
        int rank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        std::cerr << "\n[MPI_ABORT rank " << rank << "] " << msg << "\n";
        MPI_Abort(MPI_COMM_WORLD, errcode);
    }
#endif
    throw std::runtime_error(msg);
}


inline void reduce_sum(double local, double& global)
{
    #ifdef USE_MPI
        if (inited()) {
            MPI_Reduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            return;
        }
    #endif
    global = local;  
}
inline void reduce_sum(long long local, long long& global)
{
    #ifdef USE_MPI
        if (inited()) {
            MPI_Reduce(&local, &global, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
            return;
        }
    #endif
    global = local;
}


inline void allreduce_sum(double local, double& global)
{
    #ifdef USE_MPI
        if (inited()) {
            MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            return;
        }
    #endif
    global = local;
}
inline void allreduce_sum(long long local, long long& global)
{
    #ifdef USE_MPI
        if (inited()) {
            MPI_Allreduce(&local, &global, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
            return;
        }
    #endif
    global = local;
}

inline void allreduce_min(double local, double& global)
{
    #ifdef USE_MPI
        if (inited()) {
            MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
            return;
        }
    #endif
    global = local;
}

inline void allreduce_max(double local, double& global)
{
    #ifdef USE_MPI
        if (inited()) {
            MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            return;
        }
    #endif
    global = local;
}










} // namespace rgmpi
