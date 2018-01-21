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
const int res[] = { 500, 250 };
const float a = res[0] * res[1];
<<<<<<< HEAD
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
=======
#define n 1000//amount of Particle
const float p = 0.3f; //how mch room is filled with fluid
const float visc = 0.8f;	//visosity		//F=-visc*dv.x/d(p1, p2).x
const float g = 0.9f;		//gravity
const float r = 1.6f;// std::sqrt(p * a / (float)n);		//particle radius
const float h = 2 * r;
const float min = h / 5.f;
const int frameTimeMs = 10;
const float dt = 0.01f;		//time between animation steps
const float d = 0.1f; //federkonstante
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
float2 pres[n];
byte   map[(n/8 + 1) * n]; //bit map for neighbor
>>>>>>> b67cda9cc322cb7634b9ddee8863ca8301d7bf60
INT WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance,
	LPSTR lpCmdLine, int nCmdShow)
{
	float x = 0.f;
	float y = 0.f;
	for (int i = 0; i < n; ++i)
	{
<<<<<<< HEAD
		pos[i] = { (float)(i * 6 / res[1]), (float)(i % (res[1] / 6) * 6) };
=======
		y += 3.f;
		if (y >= res[1] - 1)
		{
			y = 0.f;
			x += 3.f;
		}
		pos[i] = { x, y };
>>>>>>> b67cda9cc322cb7634b9ddee8863ca8301d7bf60
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
<<<<<<< HEAD
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
=======

void getNearst() //fill adiazentz matrix
{
	const int w = n / 8 + 1;
	memset(map, 0, w * n);
	float2 d;
	for (int i = 0; i < n; ++i)
	{
		pres[i] = {0.f, 0.f};
		for (int j = 0; j < n; ++j)
		{
			if (i == j) continue;
			d = (pos[i] - pos[j]);
			if (d.sq() < h*h)
			{
				if (d.sq() == 0)
					__debugbreak();
				map[w * i + j / 8] |= 0x80 >> (j % 8);
				float c = std::sqrt(d.sq());
				c *= 0.5f;
				float y = std::sqrt(r*r - c*c);
				y *= 2.f;
				std::swap(d.x, d.y);
				d.x = -d.x;
				d *= 2.f;
				pres[i] += d;
			}
		}
		if (pres[i].x > h || pres[i].y > h)
		{
			float c = h / (pres[i].x > pres[i].y ? pres[i].x : pres[i].y);
			pres[i].x *= c;
			pres[i].y *= c;
		}
		if (pres[i].x > 2.f*h || pres[i].y > 2.f*h)
			__debugbreak();
	}
}
float distance(int id, int g) //distance between point[id] and g, point[g] e g && vel[g] || g
{
	return 1.f;
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
			if (id == i)
				__debugbreak();
			absDx = std::abs(std::sqrt(dx.sq()));
			if (isnan(absDx))
				__debugbreak();
			if (absDx < min * 0.5f)
			{
				__debugbreak();
			}
			else
			{
				dVel[id] += (dx - (dx * (h / absDx))) * 0.01f;	//bounce
				float absP = std::sqrt(pres[i].sq());
				if(absP > 0)
					dVel[id] -= dx * (dx * pres[i] / (absP * absDx) / absDx) * 0.9f;
			}
			//dVel[id] += vel[i] * (visc / distance(id, i))
>>>>>>> b67cda9cc322cb7634b9ddee8863ca8301d7bf60
		}
		wallCollision(i);
	}
<<<<<<< HEAD
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
=======
	if (!(dVel[id] == dVel[id]))
		__debugbreak();
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
			float absDx = std::sqrt(dx.sq());
			if (absDx < min)
			{
				float2 dvel[2];		//vel welche durch kollision verändert wird
				dvel[0] = dx * ((dx * vel[id]) / dx.sq());
				dvel[1] = dx * ((dx * vel[i]) / dx.sq());
				vel[id] -= dvel[0];
				vel[i] -= dvel[1];
				dvel[0] = (dvel[1] + dvel[0]) * 0.5f;
				vel[id] += dvel[0];
				vel[i] += dvel[0];
				dx *= (min / absDx);
				absDx = std::sqrt(dx.sq());
				if (isnan(absDx))
					__debugbreak();
				if (absDx < 0)
					__debugbreak();
				pos[i] += (dx * 0.5f);
				pos[id] -= (dx * 0.5f);
				id = 0;
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
	const float bD = 0.f;
	const int w = n / 8 + 1;
	for (int i = 0; i < n; ++i)
	{
		if (pos[i].x < 0)
		{
			//pos[i].x = 0.f;
			vel[i].x = -vel[i].x * bD;
			for (int j = 0; j < n; ++j)
				if (map[w*j + j / 8] & 0x80 >> (j % 8))
					vel[j].x = vel[i].x;
		}
		else if (pos[i].x >= res[0] - 1)
		{
			//pos[i].x = res[0] - 1;
			vel[i].x = -vel[i].x * bD;
			for (int j = 0; j < n; ++j)
				if (map[w*j + j / 8] & 0x80 >> (j % 8))
					vel[j].x = vel[i].x;
		}
		if (pos[i].y >= res[1] - 1)
		{
			//pos[i].y = res[1] - 1;
			vel[i].y = -vel[i].y * bD;
			for (int j = 0; j < n; ++j)
				if (map[w*j + j / 8] & 0x80 >> (j % 8))
					vel[j].y = vel[i].y;
		}
		else if (pos[i].y < 0)
		{
			//pos[i].y = 0;
			vel[i].y = -vel[i].y * bD;
			for (int j = 0; j < n; ++j)
				if (map[w*j + j / 8] & 0x80 >> (j % 8))
					vel[j].y = vel[i].y;
		}
		if (!(pos[i] == pos[i]))
			__debugbreak();
	}
}
void renderNewPic(HWND hWnd)	//flip each bit
{
	getNearst();

	calculatehalf();

	aproximateTimeStep();

	boundaryCheck();

>>>>>>> b67cda9cc322cb7634b9ddee8863ca8301d7bf60
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
			renderNewPic(hWnd);
			InvalidateRgn(hWnd, NULL, FALSE);
		}
	}
	default:
		return DefWindowProc(hWnd, Msg, wParam, lParam);
	}
	return 0;
}