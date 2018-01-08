
#include <cuda.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <atomic>
#include <thread>

#include "d:\Dokumente\OVGU\GPU\cudaSample\solution\src\cuda_util.h"

cudaError_t fluidSimulation(const int res, const float diff, const float dt, const int frames);
std::atomic_bool copyDevDesFin(false);
__global__ void addVelKernel(float* vel_x, float* vel_y, float* vel_src_x, float* vel_src_y, const float dt, const int res)
{
	int2 id;
	id.x= blockIdx.x * blockDim.x + threadIdx.x;
	id.y= blockIdx.y * blockDim.y + threadIdx.y;
	int i = id.x + id.y * res;
	vel_x[i] += dt * vel_src_x[i];
	if (vel_x[i] > 1.f)
		vel_x[i] = 1.f;
	if (vel_x[i] < -1.f)
		vel_x[i] = -1.f;
	vel_y[i] += dt * vel_src_y[i];
	if (vel_y[i] > 1.f)
		vel_y[i] = 1.f;
	if (vel_y[i] < -1.f)
		vel_y[i] = -1.f;
}
__global__ void advectKernel(float *des, float* des_fin, float* vel_x, float* vel_y, const int res, const float dt)	//berchnen wie viel prozent von welcher Zelle nach dt in der aktuellen Zelle landet
{
	int2 id;
	id.x = blockIdx.x * blockDim.x + threadIdx.x;
	id.y = blockIdx.y * blockDim.y + threadIdx.y;
	if (id.x >= res || id.y >= res)
		return;
	
	float2 d;		//travelle way from the now center Particel
	d.x = - dt * vel_x[id.y * res + id.x];	//det * vel < 0.5 !!
	d.y = - dt * vel_y[id.y * res + id.x];
	//boundary check, wall perfect reflaction
	if ((id.x == 0 && d.x < 0.f) || (id.x == res-1 && d.x > 0.f))
		d.x = 0.f;
	if ((id.y == 0 && d.y < 0.f) || (id.y == res - 1 && d.y > 0.f))
		d.y = 0.f;
	int dx = 1, dy = 1;
	if (d.x < 0.f)
	{
		dx = -1;
		d.x = -d.x;
	}
	if (d.y < 0.f)
	{
		dy = -1;
		d.y = -d.y;
	}
	des[id.x + id.y * res] = des_fin[id.x + id.y * res] * (1.f - d.x) * (1.f - d.y)
		+ des_fin[id.x + dx + id.y * res] * d.x * (1.f - d.y)
		+ des_fin[id.x + (id.y + dy) * res] * (1.f - d.x) * d.y
		+ des_fin[id.x + dx + (id.y + dy) * res] * d.x * d.y;
	if (id.x == 0 || id.x == res - 1 || id.y == 0 || id.y == res - 1)
		des[id.x + id.y * res] = des_fin[id.x + id.y * res] * d.x * d.y;
}
__global__ void diffuseKernel(float *des, float *des_fin, const int res, const float diff, const float dt)	//diffusion lässt sich mit einem Faltungsfilter simulieren
{
	//ermittlung der Position
	int2 id;
	id.x = blockIdx.x * blockDim.x + threadIdx.x;
	id.y = blockIdx.y * blockDim.y + threadIdx.y;
	//float subtraction = Eviel IMP
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
	//if (des[id.y * res + id.x] < 0.f)
	//	des[id.y * res + id.x] = 0.f;
	if (des[id.y * res + id.x] > 1.f)
		des[id.y * res + id.x] = 1.f;
}

int main()
{
	const int res = 150;		//image size, resolution per axis
	const float diff = 0.4f;	//diffusion speed
	const float visc = 0.7f;	//viscosity 
	const float dt = 0.3;	//virtual time between to frames
	const int frames = 20;	//amount of frames to render
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
void safeFrame(int num, float* picture, float* dev_picture, const int res)	//very slow
{
	cudaError_t cudaStatus = cudaMemcpy(picture, dev_picture, res * res * sizeof(float), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess)
	{
		std::cerr << "copy frame from device failed!" << std::endl;
	}
	copyDevDesFin = false;
	FILE *output;
	std::stringstream ss;
	ss << "frame" << num << ".pnm";
	output = fopen(ss.str().c_str(), "w");
	ss.str("");
	ss << "P5 " << res << ' ' << res << " 255 ";
	fprintf(output, ss.str().c_str());
	for (int i = 0; i < res*res; ++i)
	{
		if (picture[i] < 0.f)
		{
			std::cerr << "ERROR" << std::endl;
			return;
		}
		fprintf(output, "%c", ((int)(picture[i] * 255.f) == '\n' ? (int)(picture[i] * 255.f) + 1 : (int)(picture[i] * 255.f)));
	}
	fclose(output);
}

cudaError_t fluidSimulation(const int res, const float diff, const float dt, const int frames)
{
	size_t pixel = res * res;

	float *dev_des;		//field on Device with informatioon about density
	float *dev_des_fin;	//last Completed rendert density field
	float *dev_vel_x;	//field with velocity information		TODO:aproximate Velocity and decress the resolution
	float *dev_vel_x_fin;	//use float insted of float2 to reuse diffuseKernel
	float *dev_vel_y;
	float *dev_vel_y_fin;
	float *des_start;	//denisty distribution at start
	float *dev_des_src;	//particel sources
	float *dev_vel_src_x;	//velocity sources
	float *dev_vel_src_y;
	float *vel_src_x;
	float *vel_src_y;
	vel_src_x = (float*)malloc(pixel * sizeof(float));
	vel_src_y = (float*)malloc(pixel * sizeof(float));
	int center = res / 2;
	for(size_t j = 0; j < res; ++j)
		for (size_t i = 0; i < res; ++i)
		{
			if (std::abs((int)i - center) < 20 && std::abs((int)j - center) < 20)
				vel_src_y[j * res + i] = -1.f;
			else
				vel_src_y[j * res + i] = 0.f;
			vel_src_x[j * res + i] = 0.f;
		}
	des_start = (float*)malloc(pixel * sizeof(float));
	int min = (res / 4)*2;
	int max = (res / 4)*3;
	std::cout << "border for img " << min << " " << max << std::endl;
	for (size_t j = 0; j < res; ++j)
		for (size_t i = 0; i < res; ++i)
		{
			if (i >= min && i <= max && j >= min && j <= max)
				des_start[i + j * res] = 1.f;
			else
				des_start[i + j * res] = 0.f;

		}

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
	}

	cudaStatus = cudaMalloc((void**)&dev_des, pixel * sizeof(float));
	if(cudaStatus == cudaSuccess)
		cudaStatus = cudaMalloc((void**)&dev_des_fin, pixel * sizeof(float));
	if (cudaStatus == cudaSuccess)
		cudaStatus = cudaMalloc((void**)&dev_vel_x, pixel * sizeof(float));
	if (cudaStatus == cudaSuccess)
		cudaStatus = cudaMalloc((void**)&dev_vel_y, pixel * sizeof(float));
	if (cudaStatus == cudaSuccess)
		cudaStatus = cudaMalloc((void**)&dev_vel_x_fin, pixel * sizeof(float));
	if (cudaStatus == cudaSuccess)
		cudaStatus = cudaMalloc((void**)&dev_vel_y_fin, pixel * sizeof(float));
	if (cudaStatus == cudaSuccess)
		cudaStatus = cudaMalloc((void**)&dev_vel_src_x, pixel * sizeof(float));
	if (cudaStatus == cudaSuccess)
		cudaStatus = cudaMalloc((void**)&dev_vel_src_y, pixel * sizeof(float));
	if (cudaStatus != cudaSuccess)
	{
		std::cerr << "cudaMalloc failed!" << std::endl;
	}
	//TODO:Generate Velocity field
	cudaStatus = cudaMemcpy(dev_des_fin, des_start, pixel * sizeof(float), cudaMemcpyHostToDevice);
	if (cudaStatus == cudaSuccess)
		cudaStatus = cudaMemcpy(dev_vel_src_x, vel_src_x, pixel * sizeof(float), cudaMemcpyHostToDevice);
	if (cudaStatus == cudaSuccess)
		cudaStatus = cudaMemcpy(dev_vel_src_y, vel_src_y, pixel * sizeof(float), cudaMemcpyHostToDevice);
	if (cudaStatus == cudaSuccess)
		cudaStatus = cudaMemset(dev_vel_x_fin, 0, pixel * sizeof(float));
	if (cudaStatus == cudaSuccess)
		cudaStatus = cudaMemset(dev_vel_y_fin, 0, pixel * sizeof(float));
	if (cudaStatus != cudaSuccess)
	{
		std::cerr << "initialisation failed!" << std::endl;
	}
	const int MAX_THREADS_PER_BLOCK = 32;	//per dim
	unsigned int blocks, threadsPerBlock;
	threadsPerBlock = std::min((int)res, MAX_THREADS_PER_BLOCK);
	blocks = res / MAX_THREADS_PER_BLOCK;
	if (res % MAX_THREADS_PER_BLOCK != 0)
		blocks++;
	std::cout << "need " << blocks << " per dim blocks with max " << threadsPerBlock << "per block per dim" << std::endl;
	float* picture = (float*)malloc(pixel * sizeof(float));
	dim3 blockSize = dim3(blocks, blocks);
	dim3 threadSize = dim3(threadsPerBlock, threadsPerBlock);
	std::thread safePicThread;
	for (size_t frame = 0; frame < frames * 40; ++frame)
	{
		addVelKernel << <blockSize, threadSize>> > (dev_vel_x_fin, dev_vel_y_fin, dev_vel_src_x, dev_vel_src_y, dt, res);
		cudaStatus = cudaDeviceSynchronize();
		
		std::cout << "strat " << frame << std::endl;
		if (frame % 40 == 0)
		{
			if(safePicThread.joinable())
				safePicThread.join();
			copyDevDesFin = true;
			safePicThread = std::thread(safeFrame, frame, picture, dev_des_fin, res);
		}

		advectKernel << <blockSize, threadSize >> > (dev_des, dev_des_fin, dev_vel_x_fin, dev_vel_y_fin, res, dt);
		//TODO: more realistic velocity
		advectKernel <<<blockSize, threadSize>>> (dev_vel_x, dev_vel_x_fin, dev_vel_x_fin, dev_vel_y_fin, res, dt);
		advectKernel <<<blockSize, threadSize>>> (dev_vel_y, dev_vel_y_fin, dev_vel_x_fin, dev_vel_y_fin, res, dt);

		cudaStatus = cudaDeviceSynchronize();
		if (cudaStatus != cudaSuccess) {
			std::cerr << "advect Kernel failed!" << std::endl;
		}
		std::swap(dev_des, dev_des_fin);
		std::swap(dev_vel_x, dev_vel_x_fin);
		std::swap(dev_vel_y, dev_vel_y_fin);

		diffuseKernel << <blockSize, threadSize >> > (dev_des, dev_des_fin, res, diff, dt);
		diffuseKernel << <blockSize, threadSize >> > (dev_vel_x, dev_vel_x_fin, res, diff, dt);
		diffuseKernel << <blockSize, threadSize >> > (dev_vel_y, dev_vel_y_fin, res, diff, dt);
		cudaStatus = cudaDeviceSynchronize();
		if (cudaStatus != cudaSuccess) {
			std::cerr << "diffuseKernel failed!" << std::endl;
		}

		while (copyDevDesFin);
		std::swap(dev_des, dev_des_fin);
		std::swap(dev_vel_x, dev_vel_x_fin);
		std::swap(dev_vel_y, dev_vel_y_fin);
	}
	safePicThread.join();
	cudaFree(dev_des);
	cudaFree(dev_des_fin);
	cudaFree(dev_vel_x);
	cudaFree(dev_vel_y);
	cudaFree(dev_vel_x_fin);
	cudaFree(dev_vel_y_fin);

	return cudaStatus;
}