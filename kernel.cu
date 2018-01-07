
#include <cuda.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <iostream>
#include <algorithm>

#include "d:\Dokumente\OVGU\GPU\cudaSample\solution\src\cuda_util.h"

cudaError_t addWithCuda(int *c, const int *a, const int *b, unsigned int size);

__global__ void addKernel(int *c, const int *a, const int *b)
{
    int i = threadIdx.x;
    c[i] = a[i] + b[i];
}

__global__ void diffuseKernel(float *des, float *des_fin, const int res, const float diff, const float dt)	//diffusion l�sst sich mit einem Faltungsfilter simulieren
{
	//ermittlung der Position
	int2 id;
	id.x = blockIdx.x * blockDim.x + threadIdx.x;
	id.y = blockIdx.y * blockDim.y + threadIdx.y;

	if (id.x >= res || id.y >= res)
		return;
	float sum = 0.f;
	if (id.x > 0)
	{
		sum += des_fin[id.y * res + id.x - 1] - des_fin[id.y * res + id.x];
	}
	if (id.x < res - 1)
	{
		sum += des_fin[id.y * res + id.x + 1] - des_fin[id.y * res + id.x];
	}
	if (id.y > 0)
	{
		sum += des_fin[(id.y - 1) * res + id.x] - des_fin[id.y * res + id.x];
	}
	if (id.y < res - 1)
	{
		sum += des_fin[(id.y + 1) * res + id.x] - des_fin[id.y * res + id.x];
	}
	des[id.y * res + id.x] = des_fin[id.y * res + id.x] + sum * diff * dt;
}

int main()
{
	const int res = 500;		//image size, resolution per axis
	const float diff = 0.2;	//diffusion speed
	const float dt = 0.1;	//virtual time between to frames
	const int frames = 100;	//amount of frames to render
	//const float src = 1;	//denisty in source field;	TODO:Change to field with production speed per tile
	
	cudaError_t cudaStatus = fluidSimulation(res, diff, dt, frames);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "addWithCuda failed!");
        return 1;
    }
    // cudaDeviceReset must be called before exiting in order for profiling and
    // tracing tools such as Nsight and Visual Profiler to show complete traces.
    cudaStatus = cudaDeviceReset();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceReset failed!");
        return 1;
    }

    return 0;
}
cudaError_t fluidSimulation(const int res, const float diff, const float dt, const int frames)
{
	size_t pixel = res * res;

	float *dev_des;		//field on Device with informatioon about density
	float *dev_des_fin;	//last Completed rendert density field
	float2 *dev_vel;	//field with velocity information		TODO:aproximate Velocity and decress the resolution

	float *des_start;	//denisty distribution at start
	des_start = (float*)malloc(pixel * sizeof(float));
	for (size_t i = (res / 5) * 2; i < (res / 5) * 3; ++i)
		for (size_t j = (res / 5) * 2; i < (res / 5) * 3; ++j)
			des_start[i + j * res] = 1;

	int deviceCount = 0;
	cudaGetDeviceCount(&deviceCount);
	if (0 == deviceCount) {
		std::cerr << "No CUDA device found." << std::endl;
	}
	cudaDeviceProp devProp;
	cudaGetDeviceProperties(&devProp, 0);
	printDeviceProps(devProp);

	cudaError_t cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) 
	{
		std::cerr << "cudaSetDevice failed!" << std::endl;
		goto End;
	}

	cudaStatus = cudaMalloc((void**)&dev_des, pixel * sizeof(float));
	if(cudaStatus == cudaSuccess)
		cudaStatus = cudaMalloc((void**)&dev_des_fin, pixel * sizeof(float));
	if (cudaStatus == cudaSuccess)
		cudaStatus = cudaMalloc((void**)&dev_vel, pixel * sizeof(float2));
	if (cudaStatus != cudaSuccess)
	{
		std::cerr << "cudaMalloc failed!" << std::endl;
		goto End;
	}
	//TODO:Generate Velocity field
	cudaStatus = cudaMemset(dev_des, 0, pixel * sizeof(float));
	if (cudaStatus == cudaSuccess)
		cudaStatus = cudaMemcpy(dev_des_fin, des_start, pixel, cudaMemcpyHostToDevice);
	if (cudaStatus == cudaSuccess)
		cudaStatus = cudaMemset(dev_vel, 0, pixel * sizeof(float2));
	if (cudaStatus != cudaSuccess)
	{
		std::cerr << "initialisation failed!" << std::endl;
		goto End;
	}
	const int MAX_THREADS_PER_BLOCK = 1024;
	size_t blocks, threadsPerBlock;
	threadsPerBlock = std::min((int)pixel, MAX_THREADS_PER_BLOCK);
	blocks = pixel / MAX_THREADS_PER_BLOCK;
	if (pixel % MAX_THREADS_PER_BLOCK != 0)
		blocks++;
	std::cout << "need " << blocks << " blocks with max " << threadsPerBlock << "per block" << std::endl;
	diffuseKernel<<< blocks, threadsPerBlock >>>(dev_des, dev_des_fin, res, diff, dt);

End:
	cudaFree(dev_des);
	cudaFree(dev_des_fin);
	cudaFree(dev_vel);
}
cudaError_t addWithCuda(int *c, const int *a, const int *b, unsigned int size)
{
    int *dev_a = 0;
    int *dev_b = 0;
    int *dev_c = 0;
    cudaError_t cudaStatus;

    // Choose which GPU to run on, change this on a multi-GPU system.
    cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
        goto Error;
    }

    // Allocate GPU buffers for three vectors (two input, one output)    .
    cudaStatus = cudaMalloc((void**)&dev_c, size * sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    cudaStatus = cudaMalloc((void**)&dev_a, size * sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    cudaStatus = cudaMalloc((void**)&dev_b, size * sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    // Copy input vectors from host memory to GPU buffers.
    cudaStatus = cudaMemcpy(dev_a, a, size * sizeof(int), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

    cudaStatus = cudaMemcpy(dev_b, b, size * sizeof(int), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

    // Launch a kernel on the GPU with one thread for each element.
    addKernel<<<1, size>>>(dev_c, dev_a, dev_b);

    // Check for any errors launching the kernel
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "addKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
        goto Error;
    }
    
    // cudaDeviceSynchronize waits for the kernel to finish, and returns
    // any errors encountered during the launch.
    cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
        goto Error;
    }

    // Copy output vector from GPU buffer to host memory.
    cudaStatus = cudaMemcpy(c, dev_c, size * sizeof(int), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

Error:
    cudaFree(dev_c);
    cudaFree(dev_a);
    cudaFree(dev_b);
    
    return cudaStatus;
}
