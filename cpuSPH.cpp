#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <vector>
#include <memory>
#include <windows.h>
#include <WinUser.h>
#include <math.h>
#include <intrin.h>
#include <atomic>
#include <chrono>

std::atomic_bool drawing = false;
struct float2 {
	void operator+=(const float2& add)
	{
		this->x += add.x;
		this->y += add.y;
	}
	void operator-=(const float2& add)
	{
		this->x -= add.x;
		this->y -= add.y;
	}
	void operator*= (const float& vec)
	{
		*this = *this * vec;
	}
	float2 operator*(const float& f)
	{
		return { this->x * f, this->y * f };
	}
	float2 operator-(const float2& m)
	{
		return { this->x - m.x, this->y - m.y };
	}
	float2 operator+(const float2& add)
	{
	return { this->x + add.x, this->y + add.y };
	}
	float operator*(const float2& vec)
	{
		return this->x * vec.x + this->y * vec.y;
	}
	bool operator<(const float& num)
	{
		return (this->sq()) < num;
	}
	bool operator== (const float2& vec)
	{
		if (vec.x == this->x && vec.y == this->y)
			return true;
		return false;
	}
	float sq()
	{
		return this->x * this->x + this->y * this->y;
	}
	float2 abs()
	{
		return { std::abs(this->x), std::abs(this->y) };
	}
	float x;
	float y;
};

HINSTANCE hInst;
LRESULT CALLBACK WindProcedure(HWND hWnd, UINT Msg, WPARAM wParam, LPARAM lParam);
const int res[] = { 500, 250 };
const float a = res[0] * res[1];
#define n 1000//amount of Particle
const float p = 0.3f; //how mch room is filled with fluid
const float visc = 0.0001f;	//visosity		//F=-visc*dv.x/d(p1, p2).x
const float g = 0.9f;		//gravity
const float r = 4.f;// std::sqrt(p * a / (float)n);		//particle radius
const float h = 2.f * r;
const float min = h / 5.f;
const int frameTimeMs = 10;
const float dt = 0.01f;		//time between animation steps
const float d = 0.0001f; //federkonstante
const float pre = 0.4f;

BYTE *pic;
size_t bytePerLine;
float2 pos[n];
float2 posP[n];
float2 posN[n];
float2 vel[n];
float2 velP[n];
float2 velN[n];
float2 dVel[n];
float pres[n];
byte   map[(n/8 + 1) * n]; //bit map for neighbor
INT WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance,
	LPSTR lpCmdLine, int nCmdShow)
{
	float x = 0.f;
	float y = 0.f;
	for (int i = 0; i < n; ++i)
	{
		y += h;
		if (y >= res[1] - 1)
		{
			y = 0.f;
			x += h;
		}
		pos[i] = { x + (i%2 == 0 ? 0.f : 0.5f), y };
		vel[i] = { 0.f, 0.f };
	}
	WNDCLASSEX  WndCls;
	static char szAppName[] = "BitmapIntro";
	MSG         Msg;

	bytePerLine = res[0] / 8 + 1;
	bytePerLine += 4 - bytePerLine % 4;		//ro size = x * 4byte

	pic = (BYTE*)malloc(bytePerLine * res[1] * sizeof(byte));
	memset(pic, 0, bytePerLine * res[1] * sizeof(byte));

	hInst = hInstance;
	WndCls.cbSize = sizeof(WndCls);
	WndCls.style = CS_OWNDC | CS_VREDRAW | CS_HREDRAW;
	WndCls.lpfnWndProc = WindProcedure;
	WndCls.cbClsExtra = 0;
	WndCls.cbWndExtra = 0;
	WndCls.hInstance = hInst;
	WndCls.hIcon = LoadIcon(NULL, IDI_APPLICATION);
	WndCls.hCursor = LoadCursor(NULL, IDC_ARROW);
	WndCls.hbrBackground = (HBRUSH)GetStockObject(WHITE_BRUSH);
	WndCls.lpszMenuName = NULL;
	WndCls.lpszClassName = szAppName;
	WndCls.hIconSm = LoadIcon(hInstance, IDI_APPLICATION);
	RegisterClassEx(&WndCls);
	CreateWindowEx(WS_EX_OVERLAPPEDWINDOW,
		szAppName,
		"Bitmaps Fundamentals",
		WS_OVERLAPPEDWINDOW | WS_VISIBLE,
		CW_USEDEFAULT,
		CW_USEDEFAULT,
		res[0] + 50,
		res[1] + 60,
		NULL,
		NULL,
		hInstance,
		NULL);

	while (GetMessage(&Msg, NULL, 0, 0))
	{
		TranslateMessage(&Msg);
		DispatchMessage(&Msg);
	}

	return static_cast<int>(Msg.wParam);
}
void getNearst() //fill adiazentz matrix
{
	const int w = n / 8 + 1;
	memset(map, 0, w * n);
	float2 dx;
	for (int i = 0; i < n; ++i)
	{
		pres[i] = 0.f;
		for (int j = 0; j < n; ++j)
		{
			if (i == j) continue;
			dx = (pos[i] - pos[j]);
			if (dx.sq() < h*h)
			{
				if (dx.sq() == 0)
					continue;
				map[w * i + j / 8] |= 0x80 >> (j % 8);
				if (dx.sq() < min*min)
					pres[i] += 1 / (min*min);
				else
					pres[i] += 1 / dx.sq();	//equivalent zu normierter vector / länge
			}
		}
		if (pres[i] * d > 100)
			pres[i] = 100.f;
	}
}
float W(float p) //p = d^2 / h^2
{
	if (p < 1.f)
	{
		return 0.2f * (6.f * p*p - 3 * p+ 1);
	}
	else 0.f;
}
float2 deltaVel(int id, float2 *pos, float2 *vel)
{
	const int w = n / 8 + 1;
	dVel[id] = { 0.f, g };
	float absDx;
	float2 dx;
	float2 vq;
	for (int i = 0; i < n; ++i)
	{
		if (isnan(pos[i].x) || isnan(pos[i].y))
		{
			std::stringstream ss;
			ss << pos[i].x << " : " << pos[i].y;
			MessageBox(NULL, ss.str().c_str(), "ERROR", MB_ICONERROR);
			__debugbreak();
			continue;
		}
		if (!(pos[id] == pos[id]))
		{
			if (!(pos[id].x == pos[id].x))
				__debugbreak();
			else if (!(pos[id].y == pos[id].y))
				__debugbreak();
			else
				__debugbreak();
			continue;
		}
		dx = pos[i] - pos[id];
		if (isnan(dx.y) || isnan(dx.x))
			__debugbreak();
		if (map[id * w + i / 8] & 0x80 >> (i % 8))
		{

			if (dx.sq() == 0)
			{
				continue;
				//__debugbreak();
			}
			if (id == i)
				__debugbreak();
			absDx = std::abs(std::sqrt(dx.sq()));
			if (isnan(absDx))
				__debugbreak();
			if (absDx < min)
			{
				if (absDx == 0)
					continue;
				dx *= (min / absDx);
				absDx = min;
			}
			{
				dVel[id] -= dx * ((pres[i] + pres[id]) / absDx) * d;
			}
			float2 odx = { dx.y, -dx.x };//ortogonal zu dx
			if (odx.sq() == 0 || absDx == 0)
				__debugbreak();
			float2 v1 = odx * (odx*vel[id] / odx.sq());
			float2 v2 = odx * (odx*vel[i] / odx.sq());
			float2 dv = v2 - v1;
			if (!(dv == dv))
				__debugbreak();
			dVel[id] += dv * (visc / absDx);
		}
	}
	if (!(dVel[id] == dVel[id]))
		__debugbreak();
	if (dVel[id].sq() > 10000)
	{
		__debugbreak();
		float c = 10000 / dVel[id].sq();
		c = std::sqrt(c);
		dVel[id] *= c;
	}
	return dVel[id];
}
float2 deltaPos(int id, float2 *pos, float2 *vel)
{
	const int w = n / 8 + 1;
	float2 dx = { 0.f, 0.f};
	for (int i = 0; i < n; ++i)
	{
		if (map[id * w + i / 8] & 0x80 >> (1 % 8))
		{
			dx += (vel[i] - vel[id])*W((pos[i] - pos[id]).sq() / (h*h));
		}
	}
	dx *= 0.5f;
	dx += vel[id];
	if (isnan(dx.x) || isnan(dx.y))
		__debugbreak();
	return dx;
}
void ceckPos(float2 *pos)
{
	for(int id = 0; id < n; ++id)
		for (int i = id + 1; i < n; ++i)
		{
			if (isnan(pos[i].x) || isnan(pos[i].y))
				__debugbreak();
			float2 dx = pos[i] - pos[id];
			if (dx.sq() < min*min)
			{
				if (dx.sq() > 0)
				{
					float2 dvel[2];		//vel welche durch kollision verändert wird
					dvel[0] = dx * ((dx * vel[id]) / dx.sq());
					dvel[1] = dx * ((dx * vel[i]) / dx.sq());
					vel[id] -= dvel[0];
					vel[i] -= dvel[1];
					dvel[0] = (dvel[1] + dvel[0]) * 0.5f;
					vel[id] += dvel[0];
					vel[i] += dvel[0];
				}
				else
				{
					float2 vq = vel[i] + vel[id];
					vq *= 0.5f;
					vel[id] = vq;
					vel[i] = vq;
				}
				//no position correction
				/*dx *= (min / absDx);
				absDx = std::sqrt(dx.sq());
				if (isnan(absDx))
					__debugbreak();
				if (absDx < 0)
					__debugbreak();
				pos[i] += (dx * 0.5f);
				pos[id] -= (dx * 0.5f);
				id = 0;*/
			}
			if (isnan(pos[i].x) || isnan(pos[i].y) || isnan(pos[id].x) || isnan(pos[id].y))
				__debugbreak();
		}
}
void calculatehalf()
{
	for (int i = 0; i < n; ++i)
	{
		velP[i] = vel[i] + deltaVel(i, pos, vel) * (dt * 0.5f);
		if (!(velP == velP))
			__debugbreak();
		posP[i] = pos[i] + deltaPos(i, pos, vel) * dt * 0.5f;
		if (!(posP == posP))
			__debugbreak();
	}
	ceckPos(posP);
}
void aproximateTimeStep()
{
	for (int i = 0; i < n; ++i)
	{
		velN[i] = vel[i] + deltaVel(i, posP, velP) * dt;
		if (!(velN == velN))
			__debugbreak();
		posN[i] = pos[i] + deltaPos(i, posP, velP) * dt;
		if (!(posN == posN))
			__debugbreak();
	}
	ceckPos(posN);
	std::swap(posN, pos);
	std::swap(velN, vel);
}
void boundaryCheck()
{
	const float bD = 0.0f;
	const int w = n / 8 + 1;
	for (int i = 0; i < n; ++i)
	{
		if (pos[i].x < 0)
		{
			pos[i].x = 0.f;
			vel[i].x = -vel[i].x * bD;
			for (int j = 0; j < n; ++j)
				if (map[w*j + j / 8] & 0x80 >> (j % 8))
					vel[j].x = vel[i].x;
		}
		else if (pos[i].x >= res[0] - 1)
		{
			pos[i].x = res[0] - 1;
			vel[i].x = -vel[i].x * bD;
			for (int j = 0; j < n; ++j)
				if (map[w*j + j / 8] & 0x80 >> (j % 8))
					vel[j].x = vel[i].x;
		}
		if (pos[i].y >= res[1] - 1)
		{
			pos[i].y = res[1] - 1;
			vel[i].y = -vel[i].y * bD;
			for (int j = 0; j < n; ++j)
				if (map[w*j + j / 8] & 0x80 >> (j % 8))
					vel[j].y = vel[i].y;
		}
		else if (pos[i].y < 0)
		{
			pos[i].y = 0;
			vel[i].y = -vel[i].y * bD;
			for (int j = 0; j < n; ++j)
				if (map[w*j + j / 8] & 0x80 >> (j % 8))
					vel[j].y = vel[i].y;
		}
		if (!(pos[i] == pos[i]))
			__debugbreak();
	}
}
void renderNewPic(HWND hWnd, int loops)	//flip each bit
{
	std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < loops; ++i)
	{
		getNearst();

		calculatehalf();

		boundaryCheck();

		aproximateTimeStep();

		boundaryCheck();
	}
	std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
	std::stringstream ss;
	ss << time_span.count() << " seconds";
	//MessageBox(hWnd, ss.str().c_str(), "Zeit für 100", MB_ICONINFORMATION);
	memset(pic, 0, bytePerLine * res[1] * sizeof(BYTE));
	for (int i = 0; i < n; ++i)
	{
		int x = pos[i].x;
		int y = pos[i].y;
		pic[y * bytePerLine + x / 8] |= 0x80 >> (x % 8);
	}
}
int count = 0;
BITMAPINFO *bmi;
LRESULT CALLBACK WindProcedure(HWND hWnd, UINT Msg,
	WPARAM wParam, LPARAM lParam)
{	
	switch (Msg)
	{
	case WM_CREATE: {
		SetTimer(hWnd, NULL, frameTimeMs, (TIMERPROC)NULL);
		bmi = (BITMAPINFO*)malloc(sizeof(BITMAPINFO) + sizeof(RGBQUAD));
		memset(bmi, 0, sizeof(BITMAPINFO) + sizeof(RGBQUAD));
		RGBQUAD color[2];
		color[1].rgbReserved = 0x00;
		color[1].rgbBlue = 0x00;
		color[1].rgbGreen = 0xFF;
		color[1].rgbRed = 0xFF;

		color[0].rgbReserved = 0x00;
		color[0].rgbBlue = 0x00;
		color[0].rgbGreen = 0x00;
		color[0].rgbRed = 0x00;

		bmi->bmiHeader.biWidth = res[0];
		bmi->bmiHeader.biHeight = -res[1];
		bmi->bmiHeader.biBitCount = 1;
		bmi->bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
		bmi->bmiHeader.biPlanes = 1;
		bmi->bmiHeader.biClrUsed = 2;
		bmi->bmiColors[0] = color[0];
		bmi->bmiColors[1] = color[1];
		}
		break;
	case WM_DESTROY:
		KillTimer(hWnd, WM_TIMER);
		PostQuitMessage(WM_QUIT);
		break;
	case WM_PAINT: {
			HDC hDC, MemDCExercising;
			PAINTSTRUCT Ps;
			HBITMAP bmpExercising;
			hDC = BeginPaint(hWnd, &Ps);
			// Load the bitmap from the resource
			bmpExercising = CreateDIBitmap(hDC, &bmi->bmiHeader, CBM_INIT, pic, bmi, DIB_RGB_COLORS);
			// Create a memory device compatible with the above DC variable
			MemDCExercising = CreateCompatibleDC(hDC);
			// Select the new bitmap
			SelectObject(MemDCExercising, bmpExercising);

			// Copy the bits from the memory DC into the current dc
			StretchBlt(hDC, 10, 10, res[0], res[1], MemDCExercising, 0, 0, res[0], res[1], SRCCOPY);

			// Restore the old bitmap
			DeleteDC(MemDCExercising);
			DeleteObject(bmpExercising);
			EndPaint(hWnd, &Ps); 
			drawing = false;
		}
		break;
	case WM_TIMER: {
		if (!drawing)
		{
			drawing = true;
			renderNewPic(hWnd, 100);
			InvalidateRgn(hWnd, NULL, FALSE);
		}
	}
	default:
		return DefWindowProc(hWnd, Msg, wParam, lParam);
	}
	return 0;
}
