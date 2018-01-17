#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <iostream>
#include <vector>
#include <memory>
#include <windows.h>
#include <WinUser.h>
#include <math.h>

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
		return {this->x + add.x, this->y + add.y};
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
#define n 6 //amount of Particle
const float p = 0.3f; //how mch room is filled with fluid
const float visc = 0.8f;	//visosity		//F=-visc*dv.x/d(p1, p2).x
const float g = 0.9f;		//gravity
const float r = 3;// std::sqrt(p * a / (float)n);		//particle radius
const int frameTimeMs = 20;
const float dt = 0.02f;		//time between animation steps
const float roh0 = 1.f;
BYTE *pic;
size_t bytePerLine;
float2 pos[n];
float2 posN[n];
float2 vel[n];
float2 velN[n];
float2 dVel[n];
float  rho[n];
float  rhoN[n];
byte   map[(n/8 + 1) * (n/8 + 1)]; //bit map for neighbor
INT WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance,
	LPSTR lpCmdLine, int nCmdShow)
{
	for (int i = 0; i < n; ++i)
	{
		pos[i] = {(float)(i % 2)*10 + 10, (float) i / 2 * 5};
		rho[i] = 1.f;
	}
	vel[0] = {0.f, 0.f};
	vel[1] = { -2.f ,0.f };
	vel[2] = { 5.f, 0.f };
	vel[3] = { -2.f, 0.f };
	vel[4] = { 3.f, 0.f };
	vel[5] = { -4.f, 0.f};
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
float deltaRho(int id)
{
	if(vel[id].x != 0 && vel[id].y != 0)
		return -rho[id] * (dVel[id].x / vel[id].x + dVel[id].y / vel[id].y);
	//if(dVel[id].x != 0 || dVel[id].y != 0)
	//	int m = MessageBox(NULL, (LPCSTR)"Mist", (LPCSTR)"ERROR", MB_ICONWARNING);
	return 0.f;
}
float getVisc(int id1, int id2)
{
	float2 dVel = vel[id1] - vel[id2];
	float2 dPos = pos[id1] - pos[id2];
	if (dVel.x * dPos.x < 0 && dVel.y * dPos.y < 0)
	{
		return dVel * dPos / (dPos * dPos);
	}
	return 0.f;
}
float W(int id1, int id2)
{
	float2 d = pos[id1] - pos[id2];
	float p = d.sq() / (r * r);
	if (p < 0.5)
	{
		return 40 / (7 * 3.14f) * (6 * p*p*p - 6 * p*p + 1);
	}
	else if (p < 1.f)
	{
		p = 1.f - p;
		return 40 / (7 * 3.14f) * 2 * p*p*p;
	}
	else
		return 0.f;
}
float deltaW(int id1, int id2)
{
	float2 d = pos[id1] - pos[id2];
	float p = d.sq() / (r * r);
	if (p < 0.5)
	{
		return 240.f / (7.f * 3.14f) * (3 * p*p - 2 * p);
	}
	else if (p < 1.f)
	{
		p = 1.f - p;
		return -240.f / (7.f* 3.14f) * p*p;
	}
	else
		return 0.f;
}
const float c = 250.f;
float cacPresRho(int id)
{
	float r = c*c * (rho[id] - roh0);
	if (rho[id] / roh0 < 0.9f)
		return 0.f;
	else return r;
}
float2 deltaVel(int id)
{
	const int w = n / 8 + 1; //map width
	dVel[id] = { 0.f, 0.f};
	float pr = cacPresRho(id);
	for (int i = 0; i < n; ++i)
	{
		if (map[id * w + i / 8] & (0x01 << (i % 8)))
		{
			if(vel[i].x != 0)
				dVel[id].x = rho[id] * (pr + cacPresRho(i) + getVisc(id, i)) * deltaW(id, i) / vel[i].x;
			if(vel[i].y != 0)
				dVel[id].y = rho[id] * (pr + cacPresRho(i) + getVisc(id, i)) * deltaW(id, i) / vel[i].y;
		}
	}
	return dVel[id];
}
float2 deltaPos(int id)
{
	const int w = n / 8 + 1;
	float2 d = { 0.f, 0.f };
	for (int i = 0; i < n; ++i)
	{
		if (map[id * w + i / 8] & (0x01 << (i % 8)))
		{
			d += (vel[i] - vel[id]) * (W(id, i) * 2.f * rho[i] / rho[id]);
		}
	}
	d *= 0.5f;
	return d + vel[id];
}
void getNearst() //fill adiazentz matrix
{
	const int w = n / 8 + 1;
	memset(map, 0, w * w);
	float2 d;
	for(int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
		{
			if (i == j) continue;
			d = (pos[i] - pos[j]);
			if (d.sq() < r*r)
			{
				map[w * i + j / 8] |= 0x01 << (j % 8);
			}
		}
}
void calculatehalf()
{
	for (int i = 0; i < n; ++i)
	{
		velN[i] = vel[i] + deltaVel(i) * dt * 0.5f;
		if (!(velN[i] == velN[i]))
		{
			int m = MessageBox(NULL, (LPCSTR)"Mist", (LPCSTR)"ERROR", MB_ICONWARNING);
		}
		rhoN[i] = rho[i] + deltaRho(i) * dt * 0.5f;
		if (!(rhoN[i] == rhoN[i]))
		{
			int m = MessageBox(NULL, (LPCSTR)"Mist", (LPCSTR)"ERROR", MB_ICONWARNING);
		}
		posN[i] = pos[i] + deltaPos(i) * dt * 0.5f;
		if (!(posN[i] == posN[i]))
		{
			int m = MessageBox(NULL, (LPCSTR)"Mist", (LPCSTR)"ERROR", MB_ICONWARNING);
		}
	}
}
void aproximateTimeStep()
{
	for (int i = 0; i < n; ++i)
	{
		velN[i] = vel[i] + deltaVel(i) * dt;
		rhoN[i] = rho[i] + deltaRho(i) * dt;
		posN[i] = pos[i] + deltaPos(i) * dt;
	}
}
void renderNewPic()	//flip each bit
{
	getNearst();
	
	calculatehalf();
	std::swap(vel, velN);
	std::swap(rho, rhoN);
	std::swap(pos, posN);

	aproximateTimeStep();
	std::swap(vel, velN);
	std::swap(pos, posN);
	std::swap(rho, rhoN);

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