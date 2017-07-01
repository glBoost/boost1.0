#include <GL/glut.h>
#include <vector>
#include <iostream>
#include <stdio.h>
#include <memory.h>
#include <algorithm>
#include <cmath>
#include <time.h>

using namespace std;

/*  屏幕大小  */
#define Width 640
#define Height 480

// 向量类 
struct Vector {
	float u, v, w;
	Vector(float _u, float _v, float _w) :u(_u), v(_v), w(_w) {}
	Vector() :u(0.0), v(0.0), w(0.0) {}
	float operator*(const Vector &d)  //点积 
	{
		return u*d.u + v*d.v + w*d.w;
	}
	Vector operator*(float f)	//数乘
	{
		return Vector(u*f, v*f, w*f);
	}
	Vector operator-(const Vector &d)	//向量减法
	{
		return Vector(u - d.u, v - d.v, w - d.w);
	}
	Vector operator+(const Vector &d)	//向量加法
	{
		return Vector(u + d.u, v + d.v, w + d.w);
	}
	float length()  //模长 
	{
		return sqrt(u*u + v*v + w*w);
	}
	void Normalize()  //单位化 
	{
		float l = length();
		u /= l;
		v /= l;
		w /= l;
		return;
	}
};

struct Vector4D {
	float a[4];
	Vector4D() {
		memset(a, 0, sizeof(a));
	}
	Vector4D(float x, float y, float z, float w) {
		a[0] = x;
		a[1] = y;
		a[2] = z;
		a[3] = w;
	}
	Vector4D(Vector v) {
		a[0] = v.u;
		a[1] = v.v;
		a[2] = v.w;
		a[3] = 1;
	}
	float& operator[](int i) {
		return a[i];
	}
};	

//矩阵类
struct Matrix {
	float a[4][4];
	Matrix() {
		memset(a, 0, sizeof(a));
	}
	//平移
	void ModelT(float x, float y, float z) {
		a[0][0] = 1;
		a[1][1] = 1;
		a[2][2] = 1;
		a[3][3] = 1;
		a[0][3] = x;
		a[1][3] = y;
		a[2][3] = z;
	}
	//旋转
	void ModelR(float angle_x, float angle_y, float angle_z) {

	}
	//缩放
	void ModelS(float x, float y, float z) {

	}
	float* operator[](int x) {
		return a[x];
	}
	void ViewT(float x, float y, float z) {
		a[0][0] = 1;
		a[1][1] = 1;
		a[2][2] = 1;
		a[3][3] = 1;
		a[0][3] = -x;
		a[1][3] = -y;
		a[2][3] = -z;
	}
	void ViewR(Vector u, Vector v, Vector n) {
		a[0][0] = u.u;
		a[0][1] = u.v;
		a[0][2] = u.w;
		a[1][0] = v.u;
		a[1][1] = v.v;
		a[1][2] = v.w;
		a[2][0] = n.u;
		a[2][1] = n.v;
		a[2][2] = n.w;
		a[3][3] = 1;
	}
	Matrix operator*(Matrix b) {
		Matrix result;
		for(int i = 0; i < 4; ++i)
			for (int j = 0; j < 4; ++j)
			{
				float sum = 0;
				for (int k = 0; k < 4; ++k)
				{
					sum += a[i][k] * b.a[k][j];
				}
				result[i][j] = sum;
			}
		return result;
	}
};

Vector4D operator*(Matrix m, Vector4D v) {
	Vector4D result;
	for (int i = 0; i < 4; ++i)
	{
		float sum = 0;
		for (int j = 0; j < 4; ++j)
		{
			m.a[i][j] * v.a[j];
		}
		result[i] = sum;
	}
}

Vector pathway[12][1000];
Vector trace[1000];
float cameraR;

struct Color {
	int r, g, b;
	Color(int _r, int _g, int _b) :r(_r), g(_g), b(_b) {}
	Color() :r(0), g(0), b(0) {}
	Color operator*(float f)
	{
		return Color(r*f, g*f, b*f);
	}
	Color operator-(const Color &d)
	{
		return Color(r - d.r, g - d.g, b - d.b);
	}
	Color operator+(const Color &d)
	{
		return Color(r + d.r, g + d.g, b + d.b);
	}
};

/*
// 球类 
struct Sphere {
	float r;
	Vector o;
	Color color;
} s[4];

//光源位置,可手动调整
Vector light(300, 300, 100);

//一些随便设的参数 
#define MinWeight 0.0001
#define Wr 0.3
#define Wt 0.2
#define Ks 0.5
#define Kt 0.2
#define Ka 0.3
#define Kd 0.6

void RayTracing(Vector start, Vector direction, float weight, Color* color) {
	if (weight < MinWeight)
	{
		*color = Color(0, 0, 0);
	}
	else
	{
		//计算光线与所有物体的交点中离start最近的点;
		float t = 10000.0;
		int index = -1;

		for (int i = 0; i < 4; ++i)
		{
			Vector v = start - s[i].o;
			float a = direction*v, b = v*v - s[i].r * s[i].r;
			float delta = a*a - b;
			if (delta < 0)
			{
				continue;
			}
			else
			{
				float tmpt = -a - sqrt(delta);
				if (tmpt < 1e-6)
					continue;
				if (tmpt < t)
				{
					t = tmpt;
					index = i;
				}
			}
		}
		if (index == -1)
			*color = Color(0, 0, 0);
		else
		{
			Vector p = start + direction * t;	//交点位置
			Vector ld = p - light;	//光源光线入射方向
			ld.Normalize();
			Vector N = p - s[index].o;	//法线方向
			N.Normalize();
			Vector R = N * (2 * -(direction*N)) + direction;	//反射方向
			Vector lr = N * (2 * -(ld*N)) + ld;
			Vector V = direction * -1;	//视点方向

			int flag = 0;
			for (int i = 0; i < 4; ++i)
			{
				Vector v = p - s[i].o;
				Vector dir = ld * -1;
				float a = dir*v, b = v*v - s[i].r * s[i].r;
				float delta = a*a - b;
				if (delta < 0)
				{
					continue;
				}
				else
				{
					float tmpt1 = -a + sqrt(delta);
					float tmpt2 = -a - sqrt(delta);
					if (tmpt1 > 1e-3 || tmpt2 > 1e-3)
					{
						flag = 1;
						break;
					}
				}
			}

			//在交点处用局部光照模型计算出的光强;
			float I = 0;
			if (flag)
				I = Ka;
			else
				I = Ka + Kd * (max((double)-(N*ld), 0.0)) + Ks * (V*lr)*(V*lr)*(V*lr)*(V*lr); //不知道为啥镜面反射的效果很奇怪所以去掉了 
			if (I > 1)	//防止叠加的光强过大 
				I = 1;
			Color Ilocal = s[index].color * I;	//原颜色*光强得到显示亮度 

			Color *Ir = new Color();
			RayTracing(p, R, weight*Wr, Ir);	//反射
												//Vector T;
			Color *It = new Color();
			//RayTracing(p, T, weight*Wt, It);	//折射，暂时不算（因为没有折射率。。。） 
			*color = Ilocal + *Ir * Ks + *It * Kt;
			delete Ir;	//回收内存 
			delete It;
		}
	}
}


void myDisplay(void)
{
	glBegin(GL_POINTS);
	for (int i = 0; i < Width; ++i)
	{
		for (int j = 0; j < Height; ++j)
		{
			Color *c = new Color();
			RayTracing(Vector(i, j, 1000), Vector(0, 0, -1), 1, c);
			if (c->r > 255)
				c->r = 255;
			if (c->g > 255)
				c->g = 255;
			if (c->b > 255)
				c->b = 255;
			glColor3ub(c->r, c->g, c->b);
			glVertex2d(i, j);
			delete c;
		}
	}
	glEnd();
	glFlush();
}
*/

void myDisplay() {
	Vector tmp[12][1000];
	memcpy(tmp, pathway, sizeof(tmp));
	Vector camera;
}

void transform(Vector view, Vector4D Vertex) {

	Matrix ModelViewProject, Project, View, Model;
	Matrix ViewR, ViewT, ModelS, ModelR, ModelT;
	
	/*Project：投影变换矩阵
	View：视图变换矩阵
	Model：模型变换矩阵*/
	ModelViewProject = Project * View * Model;

	/*S：缩放 R：旋转 T：平移*/
	ModelViewProject = Project
		* ViewR * ViewT
		* ModelS * ModelR * ModelT;
	
	Vector4D gl_position = ModelViewProject * Vertex;
}

void InitEnvironment()
{
	glClearColor(0, 0, 0, 0);
	glClear(GL_COLOR_BUFFER_BIT);
	glPointSize(1);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluOrtho2D(0, 640, 480, 0);
}

//不用鼠标 
void OnMouse(int button, int state, int x, int y)
{

	if (button == GLUT_LEFT_BUTTON&&state == GLUT_DOWN)	 	//
	{

	}
	if (button == GLUT_RIGHT_BUTTON&&state == GLUT_DOWN)	//
	{

	}
}

//还不知道怎么保持比例，没用 
void myReshape(int w, int h)
{
	gluOrtho2D(0, w, h, 0);
	glViewport(0, Width, Height, 0);//改变显示区域，起始位置为客户端窗口左下角（非坐标原点）
	glLoadIdentity();
}


int main(int argc, char *argv[])
{
	srand(time(0));

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);
	glutInitWindowPosition(100, 100);
	glutInitWindowSize(Width, Height);
	glutCreateWindow("RayTracing");

	InitEnvironment();

	/*
	//初始化球，可手动调整参数 
	s[0].o = Vector(150, 200, 0);
	s[0].r = 70;
	s[0].color = Color(255, 0, 0);
	s[1].o = Vector(270, 200, 150);
	s[1].r = 30;
	s[1].color = Color(0, 255, 0);
	s[2].o = Vector(400, 200, 50);
	s[2].r = 100;
	s[2].color = Color(0, 0, 255);
	s[3].o = Vector(530, 200, 300);
	s[3].r = 70;
	s[3].color = Color(100, 100, 100);
	*/

	glutMouseFunc(&OnMouse);
	glutDisplayFunc(&myDisplay);
	//glutReshapeFunc(&myReshape);

	glutMainLoop();

	return 0;

}