
#include <cuda.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <atomic>
#include <thread>
#include <chrono>
#include <intrin.h>

#include "d:\Dokumente\OVGU\GPU\cudaSample\solution\src\cuda_util.h"
std::atomic_bool copyDevDesFin = false;
int main()
{

    return 0;
}

void safeFrame(int num, float* picture, float* dev_picture, const int2 res)	//very slow
{
	cudaError_t cudaStatus = cudaMemcpy(picture, dev_picture, res.x * res.y * sizeof(float), cudaMemcpyDeviceToHost);
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
	ss << "P5 " << res.x << ' ' << res.y << " 255 ";
	fprintf(output, ss.str().c_str());
	for (int i = 0; i < res.x*res.y; ++i)
	{
		if (picture[i] < 0.f)
		{
			std::cerr << "ERROR" << std::endl;
			return;
		}
		fprintf(output, "%c", picture[i] == 0 ? (int)0 : (int)255);
	}
	fclose(output);
}

cudaError_t fluidSimulation(const int2 res, const int n, const float diff, const float dt, const int frames, float* pos_x, float* pos_y)
{
	const int pixel = res.x * res.y;
	float *dev_vel_x;
	float *dev_vel_y;
	float *dev_pos_x;
	float *dev_pos_y;
	float *dev_dvl_x;
	float *dev_dvl_y;
	uint8_t *dev_pic;

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
	cudaMalloc((void**)dev_vel_x, n * sizeof(float));
	cudaMalloc((void**)dev_vel_y, n * sizeof(float));
	cudaMalloc((void**)dev_pos_x, n * sizeof(float));
	cudaMalloc((void**)dev_pos_y, n * sizeof(float));
	cudaMalloc((void**)dev_dvl_x, n * sizeof(float));
	cudaMalloc((void**)dev_dvl_y, n * sizeof(float));
	cudaMalloc((void**)dev_pic, pixel * sizeof(uint8_t));

	cudaMemset(dev_vel_x, 0, n * sizeof(float));
	cudaMemset(dev_vel_y, 0, n * sizeof(float));
	cudaMemset(dev_dvl_x, 0, n * sizeof(float));
	cudaMemset(dev_dvl_y, 0, n * sizeof(float));
	cudaMemset(dev_pic, 0, n * sizeof(float));

	cudaMemcpy(dev_pos_x, pos_x, n * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_pos_y, pos_y, n * sizeof(float), cudaMemcpyHostToDevice);

	const int MAX_THREADS_PER_BLOCK = 1024;	//per dim
	unsigned int blocks, threadsPerBlock;
	threadsPerBlock = std::min(n, MAX_THREADS_PER_BLOCK);
	blocks = n / MAX_THREADS_PER_BLOCK;
	if (n % MAX_THREADS_PER_BLOCK != 0)
		blocks++;
	std::cout << "need " << blocks << " per dim blocks with max " << threadsPerBlock << "per block per dim" << std::endl;
	float* picture = (float*)malloc(pixel * sizeof(float));
	dim3 blockSize = dim3(blocks, blocks);
	dim3 threadSize = dim3(threadsPerBlock, threadsPerBlock);


	std::thread safePicThread;
	const int STEPS_BETWEEN_FRAMES = 30;
	float *dev_p, *dev_diff;		//field to save vel diff and presuare temp
	for (size_t frame = 0; frame <= frames * STEPS_BETWEEN_FRAMES; ++frame)
	{
		//std::cout << "strat " << frame << std::endl;
		if (frame % STEPS_BETWEEN_FRAMES == 0)
		{
			if (safePicThread.joinable())
			{
				safePicThread.join();
				copyDevDesFin = true;
				safePicThread = std::thread(safeFrame, frame, picture, dev_pic, res);
				while (copyDevDesFin);
			}
			else
				__debugbreak();
		}



		cudaStatus = cudaDeviceSynchronize();
		if (cudaStatus != cudaSuccess) {
			std::cerr << "diffuseKernel failed!" << std::endl;
		}
		cudaStatus = cudaDeviceSynchronize();
	}
	safePicThread.join();
	cudaFree(dev_vel_x);
	cudaFree(dev_vel_y);
	cudaFree(dev_pos_x);
	cudaFree(dev_pos_y);
	cudaFree(dev_dvl_x);
	cudaFree(dev_dvl_y);
	cudaFree(dev_pic);

	return cudaStatus;
}