
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
#define N 100
#define AMOUNT_SM 20
__global__ void calculateMovmentKernel(unsigned int *eInG, uint4* areal, const float H, const float dt, const float visc, const float d, const float g, float *pos_x, float *pos_y, float *vel_x, float *vel_y, float *posN_x, float *posN_y, float *velN_x, float *velN_y);
__global__ void sumbissionKernel(const int2 res, const int2 fields, uint4 *dev_areal, uint32_t *dev_eInG, const float H, const float dt, const float visc, const float d, const float g, float *pos_x, float *pos_y, float *vel_x, float *vel_y, float *posN_x, float *posN_y, float *velN_x, float *velN_y)
{
	const int h = h + ((int)h) < h ? 1 : 0;
	const int2 size = { res.x / fields.x, res.y / fields.y };
	__shared__ uint32_t group[N];			//array partikel -> group // 0x00 & group = gruppe 1
	__shared__ uint4 areal[AMOUNT_SM];		//(x,y) top left corner, w = width, z = height
	__shared__ uint32_t eInG[AMOUNT_SM];//particel / group
	if (threadIdx.x < AMOUNT_SM)
	{
		eInG[threadIdx.x] = 0;
		if (threadIdx.x % fields.x == 0)	//left wall
		{
			areal[threadIdx.x].x = 0;
			areal[threadIdx.x].w = size.x + h;
		}
		else if (threadIdx.x % fields.x == fields.x - 1)
		{
			areal[threadIdx.x].x = res.x - areal[threadIdx.x - 1].x - areal[threadIdx.x - 1].w + h;
			areal[threadIdx.x].w = res.x - areal[threadIdx.x].x;
		}
		else
		{
			areal[threadIdx.x].x = areal[threadIdx.x].x + areal[threadIdx.x].w - h;
			areal[threadIdx.x].w = size.x + h + h;
		}

		if (threadIdx.y / fields.x == 0)	//left wall
		{
			areal[threadIdx.x].y = 0;
			areal[threadIdx.x].z = size.y + h;
		}
		else if (threadIdx.x / fields.x == fields.y - 1)
		{
			areal[threadIdx.x].y = res.y - areal[threadIdx.x - 1].y - areal[threadIdx.x - 1].z + h;
			areal[threadIdx.x].z = res.y - areal[threadIdx.x].y;
		}
		else
		{
			areal[threadIdx.x].y = areal[threadIdx.x].x + areal[threadIdx.x].z - h;
			areal[threadIdx.x].z = size.y + h + h;
		}
	}
	__syncthreads();

	for (int i = 0; i * blockDim.x + threadIdx.x < N; ++i)
	{
		unsigned int x = pos_x[i * blockDim.x + threadIdx.x] / (size.x + h);	//x = x cordinate von min group
		unsigned int y = pos_y[i * blockDim.x + threadIdx.x] / (size.y + h);
		unsigned int gNr = x + y * fields.x;
		group[i * blockDim.x + threadIdx.x] = 0x00000000 | (0x00000001 << gNr);
		atomicAdd(eInG + gNr, 1);
		if(pos_x[i * blockDim.x + threadIdx.x] - x * (size.x + h) >= size.x && x < fields.x - 1)	//im geteilten bereich zwischen zwei Boxen horizontal
		{
			group[i * blockDim.x + threadIdx.x] |= (0x00000001 << (gNr + 1));
			atomicAdd(eInG + gNr + 1, 1);
			if (pos_y[i * blockDim.x + threadIdx.x] - y * (size.y + h) >= size.y && y < fields.y - 1)	//im geteilten breich zwischen vier Boxen 
			{
				group[i * blockDim.x + threadIdx.x] |= (0x00000001 << (gNr + fields.x));
				group[i * blockDim.x + threadIdx.x] |= (0x00000001 << (gNr + fields.x + 1));
				eInG[gNr + fields.x] ++;
				atomicAdd(eInG + gNr + fields.x + 1, 1);
			}
		}
		else if (pos_y[i * blockDim.x + threadIdx.x] - y * (size.y + h) >= size.y && y < fields.y - 1)	//im geteilten breich zwischen zewi Boxen horizontal
		{
			group[i * blockDim.x + threadIdx.x] |= (0x00000001 << (gNr + fields.x));
			atomicAdd(eInG + gNr + fields.x, 1);
		}
	}
	__syncthreads();

	if (threadIdx.x < AMOUNT_SM)
	{
		dev_areal[threadIdx.x] = areal[threadIdx.x];
		dev_eInG[threadIdx.x] = eInG[threadIdx.x];
	}
	unsigned int maxE = 0;
	if (threadIdx.x == 0)
	{
		for (int i = 0; i < AMOUNT_SM; ++i)
			if (eInG[i] > maxE)
				maxE = eInG[i] > maxE;
	}
	__syncthreads();
	dim3 blocks(fields.x, fields.y);
	calculateMovmentKernel << <1024, blocks, 5 * maxE * sizeof(float) + maxE * sizeof(uint32_t)>> > (dev_eInG, dev_areal, H, dt, visc, d, g, pos_x, pos_y, vel_x, vel_y, posN_x, posN_y, velN_x, velN_y);
}
__global__ void calculateMovmentKernel(unsigned int *eInG, uint4* areal, const float H, const float dt, const float visc, const float d, const float g, float *pos_x, float *pos_y, float *vel_x, float *vel_y, float *posN_x, float *posN_y, float *velN_x, float *velN_y)
{
	extern __shared__ float *shared;
	const unsigned int eInA= eInG[blockIdx.x];
	float *pres;
	float *s_pos_x = shared + eInA;
	float *s_pos_y = shared + 2 * eInA;
	float *s_vel_x = shared + 3 * eInA;
	float *s_vel_y = shared + 4 * eInA;
	uint32_t *id = (uint32_t*)(shared + 5 * eInA);
	const float min = H / 30.f;
	const float max = 100.f;
	uint4 a = areal[blockIdx.x];
	__shared__ unsigned int pos;
	if (threadIdx.x == 0)
		pos = -1;
	for (unsigned int i = 0; i * blockDim.x + threadIdx.x < N; ++i)
	{
		if(pos_x[i * blockDim.x + threadIdx.x] > a.x && pos_x[i * blockDim.x + threadIdx.x] < (a.x + a.w))
			if (pos_y[i * blockDim.x + threadIdx.x] > a.y && pos_y[i * blockDim.x + threadIdx.x] < (a.y + a.z))
			{
				unsigned int p = atomicAdd(&pos, 1);
				s_pos_x[p] = pos_x[i * blockDim.x + threadIdx.x];
				s_pos_y[p] = pos_y[i * blockDim.x + threadIdx.x];
				s_vel_x[p] = vel_x[i * blockDim.x + threadIdx.x];
				s_vel_y[p] = vel_y[i * blockDim.x + threadIdx.x];
				id[p] = i * blockDim.x + threadIdx.x;
			}
	}
	__syncthreads();
	float dxSq;
	float2 dx;
	for (unsigned int i = 0; i * blockDim.x + threadIdx.x <= pos; ++i)
	{
		pres[i * blockDim.x + threadIdx.x] = 0.f;
		for (unsigned int j = 0; j <= pos; ++j)
		{
			dx.x = s_pos_x[j] - s_pos_x[i * blockDim.x + threadIdx.x];
			dx.y = s_pos_y[j] - s_pos_y[i * blockDim.x + threadIdx.x];
			dxSq = dx.x*dx.x + dx.y*dx.y;
			if (dxSq < H*H)
			{
				if (dxSq < min*min)
					pres[i * blockDim.x + threadIdx.x] += 1.f / (min*min);
				else
					pres[i * blockDim.x + threadIdx.x] += 1.f / dxSq;
			}
		}
	}
	__syncthreads();
	float2 dv = { 0.f, g };
	float absDx;
	for (unsigned int i = 0; i * blockDim.x + threadIdx.x <= pos; ++i)
	{
		int k = i * blockDim.x + threadIdx.x;
		if (s_pos_x[k] > a.x + a.w - H && ! (blockIdx.x == gridDim.x - 1)
			|| s_pos_y[k] > a.y + a.z - H && ! (blockIdx.y == gridDim.y - 1))	//if in border area and not end of screen
				continue;
		for (unsigned int j = 0; j <= pos; ++j)
		{
			dx.x = s_pos_x[j] - s_pos_x[k];
			dx.y = s_pos_y[j] - s_pos_y[k];
			dxSq = dx.x*dx.x + dx.y*dx.y;
			if (dxSq < min*min)
			{
				//calculate dVel from pressuar
				dv.x -= dx.x * (pres[j] + pres[k]) / min * d;
				dv.y -= dx.y * (pres[j] + pres[k]) / min * d;
			}
			else
			{ 
				absDx = std::sqrt(dxSq);
				//calculate dVel from pressuar
				dv.x -= dx.x * (pres[j] + pres[k]) / absDx * d;
				dv.y -= dx.y * (pres[j] + pres[k]) / absDx * d;
				//alculate dVel from visc
				dx = { -dx.y, dx.x };		//rotate 90°
				float v1_x = dx.x * (dx.x * s_vel_x[k] + dx.y * s_vel_y[k]) / dxSq;	//projection from vel[k] on ortogonal to dx
				float v1_y = dx.y * (dx.x * s_vel_x[k] + dx.y * s_vel_y[k]) / dxSq;
				float v2_x = dx.x * (dx.x * s_vel_x[j] + dx.y * s_vel_y[j]) / dxSq;
				float v2_y = dx.y * (dx.x * s_vel_x[j] + dx.y * s_vel_y[j]) / dxSq;
				float2 dvel = {v1_x - v2_x, v1_y - v2_y};
				dv.x += dvel.x * visc / absDx;
				dv.y += dvel.y * visc / absDx;
			}
		}
		float dvSq = dv.x*dv.x + dv.y*dv.y;
		if (dvSq > max*max)		//max speed
		{
			float c = max / dvSq;
			c = std::sqrt(c);
			dv.x *= c;
			dv.y *= c;
		}
		velN_x[id[k]] = s_vel_x[k] +(dv.x * dt);
		velN_y[id[k]] = s_vel_y[k] + (dv.y * dt);
		posN_x[id[k]] = s_pos_x[k] + (s_vel_x[k] * dt);
		posN_y[id[k]] = s_pos_y[k] + (s_vel_x[k] * dt);
	}
}
cudaError_t fluidSimulation(const int2 res, const int2 fields, const float r, const float dt, const float visc, const float d, const float g, const int frames, float* pos_x, float* pos_y);
int main()
{
	const int2 res = { 1000, 800 };
	const int2 fields = {5, 4};
	const float r = 1.f;
	const float dt = 0.01f;
	const float visc = 0.2f;
	const float d = 0.5f;
	const float g = 9.8f;
	const int frames = 100;
	float pos_x[N];
	float pos_y[N];
	float x = 0.f;
	float y = 0.f;
	for (int i = 0; i < N; ++i)
	{
		pos_x[i] = x + (i % 3 == 0 ? 0.5f*r : 0);
		pos_y[i] = y;
		x += 1.8f*r;
		if (x >= res.x)
		{
			x = 0.1f;
			y += 1.8f;
		}
	}
	fluidSimulation(res, fields, r, dt, visc, d, g, frames, pos_x, pos_y);
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

cudaError_t fluidSimulation(const int2 res, const int2 fields, const float r, const float dt, const float visc, const float d, const float g, const int frames, float* pos_x, float* pos_y)
{
	const float H = 2.f * r;
	const int pixel = res.x * res.y;
	float *dev_vel_x;
	float *dev_vel_y;
	float *dev_pos_x;
	float *dev_pos_y;
	float *dev_velN_x;
	float *dev_velN_y;
	float *dev_posN_x;
	float *dev_posN_y;
	uint32_t *dev_eInG;
	uint4 *dev_areal;
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
	cudaMalloc((void**)dev_vel_x, N * sizeof(float));
	cudaMalloc((void**)dev_vel_y, N * sizeof(float));
	cudaMalloc((void**)dev_pos_x, N * sizeof(float));
	cudaMalloc((void**)dev_pos_y, N * sizeof(float));
	cudaMalloc((void**)dev_velN_x, N * sizeof(float));
	cudaMalloc((void**)dev_velN_y, N * sizeof(float));
	cudaMalloc((void**)dev_posN_x, N * sizeof(float));
	cudaMalloc((void**)dev_posN_y, N * sizeof(float));
	cudaMalloc((void**)dev_pic, pixel * sizeof(uint8_t));
	cudaMalloc((void**)dev_areal, AMOUNT_SM * sizeof(uint4));

	cudaMemcpy(dev_pos_x, pos_x, N * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_pos_y, pos_y, N * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemset(dev_vel_x, 0, N * sizeof(float));
	cudaMemset(dev_vel_y, 0, N * sizeof(float));
	cudaMemset(dev_pic, 0, N * sizeof(float));

	cudaMemcpy(dev_pos_x, pos_x, N * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_pos_y, pos_y, N * sizeof(float), cudaMemcpyHostToDevice);

	const int MAX_THREADS_PER_BLOCK = 1024;	//per dim
	
	float* picture = (float*)malloc(pixel * sizeof(float));
	dim3 blockSize = dim3(5, 2);
	dim3 threadSize = dim3(1024);


	std::thread safePicThread;
	const int STEPS_BETWEEN_FRAMES = 30;
	float *dev_p, *dev_diff;		//field to save vel diff and presuare temp
	for (size_t frame = 0; frame <= frames * STEPS_BETWEEN_FRAMES; ++frame)
	{
		//std::cout << "strat " << frame << std::endl;
#ifdef DRAW_PIC
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
#endif
		//calculate 1 frame
		sumbissionKernel << <MAX_THREADS_PER_BLOCK, 1 >> >
			(res, fields, dev_areal, dev_eInG, H, dt, visc, d, g,
				dev_pos_x, dev_pos_y, dev_vel_x, dev_vel_y, dev_posN_x, dev_posN_y, dev_velN_x, dev_velN_y);

		cudaStatus = cudaDeviceSynchronize();
		if (cudaStatus != cudaSuccess) {
			std::cerr << "diffuseKernel failed!" << std::endl;
		}
		std::swap(dev_posN_x, dev_pos_x);
		std::swap(dev_posN_y, dev_pos_y);
		std::swap(dev_velN_x, dev_vel_x);
		std::swap(dev_velN_y, dev_vel_y);
	}
	safePicThread.join();
	cudaFree(dev_vel_x);
	cudaFree(dev_vel_y);
	cudaFree(dev_pos_x);
	cudaFree(dev_pos_y);
	cudaFree(dev_velN_x);
	cudaFree(dev_velN_y);
	cudaFree(dev_posN_x);
	cudaFree(dev_posN_y);
	cudaFree(dev_areal);
	cudaFree(dev_pic);

	return cudaStatus;
}