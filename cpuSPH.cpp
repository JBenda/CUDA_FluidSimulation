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
	float x;
	float y;
};

HINSTANCE hInst;
LRESULT CALLBACK WindProcedure(HWND hWnd, UINT Msg, WPARAM wParam, LPARAM lParam);
const int res[] = { 100, 50 };
const float a = res[0] * res[1];
#define n 100//amount of Particle
const float p = 0.3f; //how much room is filled with fluid
const float visc = 0.8f;	//visosity		//F=-visc*dv.x/d(p1, p2).x
const float g = 0.9f;		//gravity
const float r = 3;// std::sqrt(p * a / (float)n);		//particle radius
const int frameTimeMs = 20;
const float dt = 0.02f;		//time between animation steps
BYTE *pic;
size_t bytePerLine;
float2 pos[n];
float2 vel[n];
INT WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance,
	LPSTR lpCmdLine, int nCmdShow)
{
	for (int i = 0; i < n; ++i)
	{
		pos[i] = { (float)(i * 6 / res[1]), (float)(i % (res[1] / 6) * 6) };
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
void ballCollision(int b1, int b2)
{
	float2 dv = vel[b1] - vel[b2];
	float2 dx = pos[b1] - pos[b2];
	float normSq = dx.sq();
	float norm = std::sqrt(normSq);
	if (norm - 2.f * r< dt * std::abs(dv*dx) / norm)
	{
		float2 a = dx * ((dv*dx) / normSq) + float2({0.f, g});
		vel[b1] = vel[b1] - a * dt;
		vel[b2] = vel[b2] + a * dt;
	}
}
void wallCollision(int id)
{
	if (pos[id].x - r < 0 || pos[id].x + r > res[0])
		vel[id].x = vel[id].x * (-1.f);
	if (pos[id].y - r < 0 || pos[id].x + r > res[1])
		vel[id].y = vel[id].y * (-1.f);
}
void calculateCollision()
{
	for (int i = 0; i < n; ++i)
	{
		for (int j = i + 1; j < n; ++j)
		{
			ballCollision(i, j);
		}
		wallCollision(i);
	}
}
void moveStep()
{
	for (int i = 0; i < n; ++i)
	{
		pos[i] += vel[i] * dt;
	}
}
void renderNewPic(HWND hWnd)
{
	calculateCollision();
	moveStep();
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
			StretchBlt(hDC, 10, 10, res[0] * 8, res[1] * 8, MemDCExercising, 0, 0, res[0], res[1], SRCCOPY);

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
			renderNewPic(hWnd);
			InvalidateRgn(hWnd, NULL, FALSE);
		}
	}
	default:
		return DefWindowProc(hWnd, Msg, wParam, lParam);
	}
	return 0;
}