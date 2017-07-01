#include <GL/glut.h>
#include <vector>
#include <iostream>
#include <stdio.h>
#include <memory.h>
#include <algorithm>
#include <cmath>
#include <time.h>

using namespace std;

/*  ��Ļ��С  */
#define Width 640
#define Height 480

// ������ 
struct Vector {
	float u, v, w;
	Vector(float _u, float _v, float _w) :u(_u), v(_v), w(_w) {}
	Vector() :u(0.0), v(0.0), w(0.0) {}
	float operator*(const Vector &d)  //��� 
	{
		return u*d.u + v*d.v + w*d.w;
	}
	Vector operator*(float f)	//����
	{
		return Vector(u*f, v*f, w*f);
	}
	Vector operator-(const Vector &d)	//��������
	{
		return Vector(u - d.u, v - d.v, w - d.w);
	}
	Vector operator+(const Vector &d)	//�����ӷ�
	{
		return Vector(u + d.u, v + d.v, w + d.w);
	}
	float length()  //ģ�� 
	{
		return sqrt(u*u + v*v + w*w);
	}
	void Normalize()  //��λ�� 
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

//������
struct Matrix {
	float a[4][4];
	Matrix() {
		memset(a, 0, sizeof(a));
	}
	//ƽ��
	void ModelT(float x, float y, float z) {
		a[0][0] = 1;
		a[1][1] = 1;
		a[2][2] = 1;
		a[3][3] = 1;
		a[0][3] = x;
		a[1][3] = y;
		a[2][3] = z;
	}
	//��ת
	void ModelR(float angle_x, float angle_y, float angle_z) {

	}
	//����
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
// ���� 
struct Sphere {
	float r;
	Vector o;
	Color color;
} s[4];

//��Դλ��,���ֶ�����
Vector light(300, 300, 100);

//һЩ�����Ĳ��� 
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
		//�����������������Ľ�������start����ĵ�;
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
			Vector p = start + direction * t;	//����λ��
			Vector ld = p - light;	//��Դ�������䷽��
			ld.Normalize();
			Vector N = p - s[index].o;	//���߷���
			N.Normalize();
			Vector R = N * (2 * -(direction*N)) + direction;	//���䷽��
			Vector lr = N * (2 * -(ld*N)) + ld;
			Vector V = direction * -1;	//�ӵ㷽��

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

			//�ڽ��㴦�þֲ�����ģ�ͼ�����Ĺ�ǿ;
			float I = 0;
			if (flag)
				I = Ka;
			else
				I = Ka + Kd * (max((double)-(N*ld), 0.0)) + Ks * (V*lr)*(V*lr)*(V*lr)*(V*lr); //��֪��Ϊɶ���淴���Ч�����������ȥ���� 
			if (I > 1)	//��ֹ���ӵĹ�ǿ���� 
				I = 1;
			Color Ilocal = s[index].color * I;	//ԭ��ɫ*��ǿ�õ���ʾ���� 

			Color *Ir = new Color();
			RayTracing(p, R, weight*Wr, Ir);	//����
												//Vector T;
			Color *It = new Color();
			//RayTracing(p, T, weight*Wt, It);	//���䣬��ʱ���㣨��Ϊû�������ʡ������� 
			*color = Ilocal + *Ir * Ks + *It * Kt;
			delete Ir;	//�����ڴ� 
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
	
	/*Project��ͶӰ�任����
	View����ͼ�任����
	Model��ģ�ͱ任����*/
	ModelViewProject = Project * View * Model;

	/*S������ R����ת T��ƽ��*/
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

//������� 
void OnMouse(int button, int state, int x, int y)
{

	if (button == GLUT_LEFT_BUTTON&&state == GLUT_DOWN)	 	//
	{

	}
	if (button == GLUT_RIGHT_BUTTON&&state == GLUT_DOWN)	//
	{

	}
}

//����֪����ô���ֱ�����û�� 
void myReshape(int w, int h)
{
	gluOrtho2D(0, w, h, 0);
	glViewport(0, Width, Height, 0);//�ı���ʾ������ʼλ��Ϊ�ͻ��˴������½ǣ�������ԭ�㣩
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
	//��ʼ���򣬿��ֶ��������� 
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