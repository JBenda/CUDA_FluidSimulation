#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <iostream>
#include <vector>
#include <memory>
#include <windows.h>
#include <WinUser.h>
#include "resource.h"

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
		return {this->x + add.x, this->y + add.y};
	}
	bool operator<(const float& num)
	{
		return (this->sq()) < num;
	}
	float sq()
	{
		return this->x * this->x + this->y * this->y;
	}
	float x;
	float y;
};

HINSTANCE hInst;
LRESULT CALLBACK WindProcedure(HWND hWnd, UINT Msg, WPARAM wParam, LPARAM lParam);
const int res[] = { 200, 100 };
const float a = res[0] * res[1];
#define n 2 //amount of Particle
const float p = 0.3f; //how mch room is filled with fluid
const float visc = 0.8f;	//visosity		//F=-visc*dv.x/d(p1, p2).x
const float g = 0.9f;		//gravity
const float r = 1;// std::sqrt(p * a / (float)n);		//particle radius
const int frameTimeMs = 25;
const float dt = 0.025f;		//time between animation steps
BYTE *pic;
size_t bytePerLine;
float2 pos[n];
float2 posN[n];
float2 vel[n];
float2 velN[n];
INT WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance,
	LPSTR lpCmdLine, int nCmdShow)
{
	for (int i = 0; i < n; ++i)
	{
		pos[i] = {(float)(i % res[0])*40, float(i / res[0])};
		vel[i] = {i==0? 20.f: 1.f, 0.f};
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
		res[0] * 8 + 50,
		res[1] * 8 + 60,
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

void checkBound(int id)
{
	if (posN[id].x < 0)
	{
		posN[id].x = 0;
		velN[id].x = 0.f;
	}
	if (posN[id].x >= res[0])
	{
		posN[id].x = res[0] - 1;
		velN[id].x = 0.f;
	}
	if (posN[id].y < 0)
	{
		posN[id].y = 0;
		velN[id].y = 0.f;
	}
	if (posN[id].y >= res[1])
	{
		posN[id].y = res[1] - 1;
		velN[id].y = 0.f;
	}
}
void physik(int id)
{
	float2 force = { 0.f, 0.f };
	force.y = 0.2f;
	float absF;
	for (int i = 0; i < n; ++i)
	{
		if (i == id) continue;
		float2 d;
		d = pos[id] - pos[i];
		if (std::abs(d.x )< 1 && d.sq() > 0)	//in radius
		{
			if((d.x > 0 && vel[i].x > 0) || (d.x < 0 && vel[i].x < 0))
				force.x = (vel[i].x  - vel[id].x);
		}
	}
	velN[id] = vel[id] + force;
	posN[id] = pos[id] + vel[id] * dt;
	checkBound(id);
}
void renderNewPic()	//flip each bit
{
	for(int i = 0; i < n; ++i)
		physik(i);
	std::swap(vel, velN);
	std::swap(pos, posN);
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
			StretchBlt(hDC, 10, 10, res[0] * 8, res[1] * 8, MemDCExercising, 0, 0, res[0], res[1], SRCCOPY);

			// Restore the old bitmap
			DeleteDC(MemDCExercising);
			DeleteObject(bmpExercising);
			EndPaint(hWnd, &Ps); 
		}
		break;
	case WM_TIMER: {
		renderNewPic();
		InvalidateRgn(hWnd, NULL, FALSE);
	}
	default:
		return DefWindowProc(hWnd, Msg, wParam, lParam);
	}
	return 0;
}