#pragma once

void launch_gpu_kernel(int kernel_type, int blocksize, int threadsperblock,
    int Nentries, int Nloop, const float *din, float *dout);
