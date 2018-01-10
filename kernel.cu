
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
__global__ void addSrcKernel(float* vel_x, float* vel_y, float* des, float* vel_src_x, float* vel_src_y, float* des_src, const float dt, const int res)
{
	int2 id;
	id.x= blockIdx.x * blockDim.x + threadIdx.x;
	id.y= blockIdx.y * blockDim.y + threadIdx.y;
	int i = id.x + id.y * res;
	if (i >= res*res)
		return;

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
	
	des[i] += dt * des_src[i];
	if (des[i] > 1.f)
		des[i] = 1.f;
	if (des[i] < 0.f)
		des[i] = 0.f;
}
__global__ void velDiffKernel(float* vel_x, float* vel_y, float* diff,const int res)
{
	int2 id;
	id.x = blockIdx.x * blockDim.x + threadIdx.x;
	id.y = blockIdx.y * blockDim.y + threadIdx.y;
	if (id.x >= res || id.y >= res)
		return;
	int i = id.x + id.y * res;

	diff[i] = -0.5f * (std::abs(vel_x[i - 1] - vel_x[i + 1]) + std::abs(vel_y[i - res] - vel_y[i + res]));
}
__global__ void presuerKernel(float* diff, float* p, const int res)
{
	int2 id;
	id.x = blockIdx.x * blockDim.x + threadIdx.x;
	id.y = blockIdx.y * blockDim.y + threadIdx.y;
	if (id.x >= res || id.y >= res)
		return;
	int i = id.x + id.y * res;

	p[i] = (diff[i] + diff[i + 1] + diff[i - 1] + diff[i + res] + diff[i - res]) / 5.f;
}
__global__ void pressVelKernel(float* vel_x, float*  vel_y, float* p, const int res, const float dt)
{
	int2 id;
	id.x = blockIdx.x * blockDim.x + threadIdx.x;
	id.y = blockIdx.y * blockDim.y + threadIdx.y;
	if (id.x >= res || id.y >= res)
		return;
	int i = id.x + id.y * res;
	vel_x[i] += (p[i + 1] - p[i - 1]) * dt;		//IMP
	if (vel_x[i] > 1.f)
		vel_x[i] = 1.f;
	if (vel_x[i] < -1.f)
		vel_x[i] = -1.f;

	vel_y[i] += (p[i + res] - p[i - res]) * dt;
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
	d.y = -dt * vel_y[id.y * res + id.x];

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
	//TODO denisty constant
	
	if (id.x == 0 || id.x == res - 1 || id.y == 0 || id.y == res - 1)	//boundarey, no fluid leg
	{
		des[id.x + id.y * res] = des_fin[id.x + id.y * res] * (1.f - d.x) * (1.f - d.y);
		if ((id.x != 0 || dx != -1) && (id.x != res - 1 || dx != 1))
			des[id.x + id.y * res] += des_fin[id.x + dx + id.y * res] * d.x * (1.f - d.y);
		//else
		//	des[id.x + id.y * res] -= des_fin[id.x + id.y * res] * d.x * (1.f - d.y);
		if ((id.y != 0 || dy != -1) && (id.y != res - 1 || dy != 1))
			des[id.x + id.y * res] += des_fin[id.x + (id.y + dy) * res] * (1.f - d.x) * d.y;
		//else
		//	des[id.x + id.y * res] -= des_fin[id.x + id.y * res] * (1.f - d.x) * d.y;
		if ((id.y != 0 || dy != -1) && (id.y != res - 1 || dy != 1) && (id.x != 0 || dx != -1) && (id.x != res - 1 || dx != 1))
			des[id.x + id.y * res] += des_fin[id.x + dx + (id.y + dy) * res] * d.x * d.y;
		//else
		//	des[id.x + id.y * res] -= des_fin[id.x + id.y * res] * d.x * d.y;
	}
	else
	{
		des[id.x + id.y * res] = des_fin[id.x + id.y * res] * (1.f - d.x) * (1.f - d.y)
			+ des_fin[id.x + dx + id.y * res] * d.x * (1.f - d.y)
			+ des_fin[id.x + (id.y + dy) * res] * (1.f - d.x) * d.y
			+ des_fin[id.x + dx + (id.y + dy) * res] * d.x * d.y;
	}
}
__global__ void diffuseKernel(float *des, float *des_fin, const int res, const float diff, const float dt)	//diffusion l�sst sich mit einem Faltungsfilter simulieren
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
	const int2 res = { 500, 300 };//pic resolution
	const int n = 500; //amount of Particle
	const float p = 0.3f; //how mch room is filled with fluid
	const float dt = 0.3f;		//time between animation steps
	const float visc = 0.8f;	//visosity		//F=-visc*dv.x/d(p1, p2).x
	const float g = 0.9f;		//gravity
	const float r = std::sqrt(p / (float) n);		//particle radius
	const float sf = 0.5f;		//value define distance between surface and particle
    return 0;
}
float temperatur(float x1, float y1, float x2, float y2)	//calculate temp on p1 when p2 is a heat source
{
	return 0.f;
}
void bounc(float2 *f, float pos_x, float pos_y, float f_pos_x, float f_pos_y)
{
	float dx = std::abs(pos_x - f_pos_x), dy = std::abs(pos_y - f_pos_y);
	float force = 1.f / (dx*dx + dy*dy);
	force /= (dx + dy);
	(*f).x += force * dx;
	(*f).y += force * dy;
}
void sphFluidSimulationCPU(const float2 res, const int n, const float dt, const float visc, const float g, const float r, const float sf)
{
	float *par_vel_x;		//velocity from Particle
	float *par_vel_y;
	float *par_vel_x_tmp;
	float *par_vel_y_tmp;
	float *par_pos_x;		//position particle
	float *par_pos_y;
	float *par_pos_x_tmp;
	float *par_pos_y_tmp;

	//calculate force on Particle
	for (int i = 0; i < n; ++i)
	{
		float2 f = {0.f, g};
		int amountN;		//amount nightbour
		int *nei;			//nighbour			partice in interactionrad bzw d(p1, p2) < 2*r
		for (int k = 0; k < amountN; ++k)
		{
			f.x += visc * (par_vel_x[nei[k]] - par_vel_x[i]) * std::abs(par_pos_y[nei[k]] - par_pos_y[i]);
			f.y += visc * (par_vel_y[nei[k]] - par_vel_y[i]) * std::abs(par_pos_x[nei[k]] - par_pos_x[i]);
			bounc(&f, par_pos_x[i], par_pos_y[i], par_pos_x[nei[k]], par_pos_y[nei[k]]);					//calculate bounce force
		}
		f.x *= dt;	//f * t / m = dv, m = 1 => f * t = dv
		f.y *= dt;
		par_vel_x_tmp[i] = par_vel_x[i] + f.x;
		par_vel_y_tmp[i] = par_vel_y[i] + f.y;

		par_pos_x_tmp[i] = par_vel_x_tmp[i] * dt;
		par_pos_y_tmp[i] = par_vel_y_tmp[i] * dt;
	}
	std::swap(par_pos_x_tmp, par_pos_x);
	std::swap(par_pos_x_tmp, par_pos_x);
	std::swap(par_pos_x_tmp, par_pos_x);
	std::swap(par_pos_x_tmp, par_pos_x);
	float temp;
	float step = 1.f/(float)std::min(res.x, res.y);
	std::stringstream ss;
	ss.str("P1 ");
	ss << res.x << ' ' << res.y;
	for (int i = 0; i < res.y; ++i)
	{
		for (int j = 0; j < res.x; ++j)
		{
			temp = 0.f;
			for (int p = 0; p < n; ++p)
			{
				temp += temperatur((float)i*step, (float)j*step, par_pos_x[p], par_pos_y[p]);
			}
			if (temp >= sf)
				ss << " 1";
			else
				ss << " 0";
		}
	}
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
		fprintf(output, "%c", (255 - (int)(picture[i] * 255.f) == '\n' ? 255 - (int)(picture[i] * 255.f) + 1 : 255 - (int)(picture[i] * 255.f)));
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
	float *des_src;
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
			if (std::abs((int)j - center) < 2 && std::abs((int)i - center) < 10)
				vel_src_y[j * res + i] = -0.8f;
			else
				vel_src_y[j * res + i] = 0.f;
			vel_src_x[j * res + i] = 0.f;
		}
	
	des_start = (float*)malloc(pixel * sizeof(float));
	int min = (res / 15)*7;
	int max = (res / 15)*8;
	std::cout << "border for img " << min << " " << max << std::endl;
	for (size_t j = 0; j < res; ++j)
		for (size_t i = 0; i < res; ++i)
		{
			if (std::abs((int)j - center) < 2 && std::abs((int)i - center) < 10)
				des_start[i + j * res] = .1f;
			else
				des_start[i + j * res] = 0.f;

		}

	des_src = des_start;

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
	if (cudaStatus == cudaSuccess)
		cudaStatus = cudaMalloc((void**)&dev_des_src, pixel * sizeof(float));
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
	if (cudaStatus == cudaSuccess)
		cudaStatus = cudaMemcpy(dev_des_src, des_src, pixel * sizeof(float), cudaMemcpyHostToDevice);
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
	const int STEPS_BETWEEN_FRAMES = 30;
	float *dev_p, *dev_diff;		//field to save vel diff and presuare temp
	for (size_t frame = 0; frame <= frames * STEPS_BETWEEN_FRAMES; ++frame)
	{
		//std::cout << "strat " << frame << std::endl;
		if (frame % STEPS_BETWEEN_FRAMES == 0)
		{
			if(safePicThread.joinable())
				safePicThread.join();
			copyDevDesFin = true;
			safePicThread = std::thread(safeFrame, frame, picture, dev_des_fin, res);
		}
		diffuseKernel << <blockSize, threadSize >> > (dev_des, dev_des_fin, res, diff, dt);
		diffuseKernel << <blockSize, threadSize >> > (dev_vel_x, dev_vel_x_fin, res, diff, dt);
		diffuseKernel << <blockSize, threadSize >> > (dev_vel_y, dev_vel_y_fin, res, diff, dt);
		cudaStatus = cudaDeviceSynchronize();
		if (cudaStatus != cudaSuccess) {
			std::cerr << "diffuseKernel failed!" << std::endl;
		}
		
		while (copyDevDesFin);	//wait for pic generation
		std::swap(dev_des, dev_des_fin);
		std::swap(dev_vel_x, dev_vel_x_fin);
		std::swap(dev_vel_y, dev_vel_y_fin);

		//calculate very aproxed Velocety evolution over time
		//1. Calculate vel diff
		//2. Calculate preser from this
		//3. Change vel						TODO: make more effizient
		dev_p = dev_vel_y;				//dev_vel_y..Buffer field to temp save pressuar
		dev_diff = dev_vel_x;			//dev_vel_x..Buffer Field
		velDiffKernel << <blockSize, threadSize >> > (dev_vel_x_fin, dev_vel_y_fin, dev_diff, res);
		cudaStatus = cudaDeviceSynchronize();
		presuerKernel << <blockSize, threadSize >> > (dev_diff, dev_p, res);
		cudaStatus = cudaDeviceSynchronize();
		pressVelKernel << <blockSize, threadSize >> > (dev_vel_x_fin, dev_vel_y_fin, dev_p, res, dt);
		cudaStatus = cudaDeviceSynchronize();
		if (cudaStatus != cudaSuccess) {
			std::cerr << "presuerCalculation failed!" << std::endl;
		}
		advectKernel << <blockSize, threadSize >> > (dev_des, dev_des_fin, dev_vel_x_fin, dev_vel_y_fin, res, dt);
		advectKernel <<<blockSize, threadSize>>> (dev_vel_x, dev_vel_x_fin, dev_vel_x_fin, dev_vel_y_fin, res, dt);
		advectKernel <<<blockSize, threadSize>>> (dev_vel_y, dev_vel_y_fin, dev_vel_x_fin, dev_vel_y_fin, res, dt);

		cudaStatus = cudaDeviceSynchronize();
		if (cudaStatus != cudaSuccess) {
			std::cerr << "advect Kernel failed!" << std::endl;
		}

		std::swap(dev_des, dev_des_fin);
		std::swap(dev_vel_x, dev_vel_x_fin);
		std::swap(dev_vel_y, dev_vel_y_fin);

		addSrcKernel << <blockSize, threadSize >> > (dev_vel_x_fin, dev_vel_y_fin, dev_des_fin, dev_vel_src_x, dev_vel_src_y, dev_des_src, dt, res);
		cudaStatus = cudaDeviceSynchronize();
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