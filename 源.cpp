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
#define WIDTH 640
#define HEIGHT 480
#define PI 3.14159265

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
	Vector cross(Vector V) {
		Vector result;
		result.u = v*V.w - w*V.v;
		result.v = w*V.u - u*V.w;
		result.w = u*V.v - v*V.u;
		return result;
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
	//w归一化
	void identity() {
		if (fabs(a[3]) < 1e-6)
		{
			a[0] *= 10000;
			a[1] *= 10000;
			a[2] *= 0;
			return;
		}
		a[0] /= a[3];
		a[1] /= a[3];
		a[2] /= a[3];
		return;
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
		memset(a, 0, sizeof(a));
		a[0][0] = 1;
		a[1][1] = 1;
		a[2][2] = 1;
		a[3][3] = 1;
		a[0][3] = x;
		a[1][3] = y;
		a[2][3] = z;
	}
	//旋转
	void ModelRx(float angle_x) {
		memset(a, 0, sizeof(a));
		float s = sin(angle_x * PI / 180);
		float c = cos(angle_x * PI / 180);
		a[0][0] = 1;
		a[1][1] = a[2][2] = c;
		a[2][1] = s;
		a[1][2] = -s;
		a[3][3] = 1;
		return;
	}
	void ModelRy(float angle_y) {
		memset(a, 0, sizeof(a));
		float s = sin(angle_y * PI / 180);
		float c = cos(angle_y * PI / 180);
		a[1][1] = 1;
		a[2][2] = a[0][0] = c;
		a[0][2] = s;
		a[2][0] = -s;
		a[3][3] = 1;
		return;
	}
	void ModelRz(float angle_z) {
		memset(a, 0, sizeof(a));
		float s = sin(angle_z * PI / 180);
		float c = cos(angle_z * PI / 180);
		a[2][2] = 1;
		a[0][0] = a[1][1] = c;
		a[1][0] = s;
		a[0][1] = -s;
		a[3][3] = 1;
		return;
	}
	//缩放
	void ModelS(float fx, float fy, float fz) {
		memset(a, 0, sizeof(a));
		a[0][0] = fx;
		a[1][1] = fy;
		a[2][2] = fz;
		a[3][3] = 1;
	}
	float* operator[](int x) {
		return a[x];
	}
	void ViewT(float x, float y, float z) {
		memset(a, 0, sizeof(a));
		a[0][0] = 1;
		a[1][1] = 1;
		a[2][2] = 1;
		a[3][3] = 1;
		a[0][3] = -x;
		a[1][3] = -y;
		a[2][3] = -z;
	}
	void ViewR(Vector u, Vector v, Vector n) {
		memset(a, 0, sizeof(a));
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
	void Project(float z_near, float z_far, int width, int height) {
		memset(a, 0, sizeof(a));
		a[0][0] = z_near / (width / 2);
		a[1][1] = z_near / (height / 2);
		a[2][2] = -(z_far + z_near) / (z_far - z_near);
		a[2][3] = 2 * z_far*z_near / (z_far - z_near);
		a[3][2] = 1;
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
	friend ostream& operator<<(ostream& o, const Matrix &m) {
		for (int i = 0; i < 4; ++i)
		{
			for (int j = 0; j < 4; ++j)
			{
				o << m.a[i][j] << ' ';
			}
			o << endl;
		}
		return o;
	}
};

Vector4D operator*(Matrix m, Vector4D v) {
	Vector4D result;
	for (int i = 0; i < 4; ++i)
	{
		float sum = 0;
		for (int j = 0; j < 4; ++j)
		{
			sum += m.a[i][j] * v.a[j];
		}
		result[i] = sum;
	}
	return result;
}

/////////////////////////////////////////
// B样条建模
//	     by 文敬瑄
/////////////////////////////////////////

#define X 0
#define Y 1
#define Z 2
#define MM 20
typedef float Vector1[3];

Vector1 origin[12];
Vector1 cir[12][500];
Vector1 points[9];
int n = 8;
int k = 7;
float tt[16];
int count = 50;
Vector1 V[500];

void Init(void)
{

	for (int q = 0; q < 12; q++) {
		origin[q][X] = 0.025*cos(PI / 6 * q);
		origin[q][Y] = 0.025*sin(PI / 6 * q);
		origin[q][Z] = 0;
	}

	points[0][X] = 0.5f;
	points[0][Y] = 0.5f;
	points[0][Z] = 0.5f;

	points[1][X] = -0.5f;
	points[1][Y] = 0.5f;
	points[1][Z] = 0.5f;

	points[2][X] = -0.5f;
	points[2][Y] = -0.5f;
	points[2][Z] = 0.5f;

	points[3][X] = 0.5f;
	points[3][Y] = -0.5f;
	points[3][Z] = 0.5f;

	points[4][X] = 0.5f;
	points[4][Y] = -0.5f;
	points[4][Z] = -0.5f;

	points[5][X] = -0.5f;
	points[5][Y] = -0.5f;
	points[5][Z] = -0.5f;

	points[6][X] = -0.5f;
	points[6][Y] = 0.5f;
	points[6][Z] = -0.5f;

	points[7][X] = 0.5f;
	points[7][Y] = 0.5f;
	points[7][Z] = -0.5f;

	points[8][X] = 0.5f;
	points[8][Y] = 0.5f;
	points[8][Z] = 0.5f;


	for (int i = 0; i < 16; i++) {
		tt[i] = 1.0 / 15 * i;
	}


	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glColor3f(1.0f, 1.0f, 0.0f);
	// 把着色模式设置为单调着色
	glShadeModel(GL_FLAT); //glShadeModel(GL_SMOOTH);
						   // 把顺时针环绕的多边形设为正面，这与默认是相反的，因为我们使用的是三角形扇
	glFrontFace(GL_CW);
	glOrtho(-1.0f, 1.0f, -1.0f, 1.0f, -1.0f, 1.0f);
}


void Line(Vector1 a, Vector1 b) {
	glLineWidth(2);
	glBegin(GL_LINES);
	glColor3f(0.2, 0.4, 0.7);
	glPolygonMode(GL_FRONT, GL_FILL);
	glVertex3f(a[X], a[Y], a[Z]);
	glVertex3f(b[X], b[Y], b[Z]);
	glEnd();
	glFlush();
}

void Bspline(Vector1 P[], float T[], int k, float t, int j, Vector1 V) {
	int i, r, temp, temp1;
	Vector1 Q[MM];
	float lamta;
	temp = j - k + 1;
	for (i = 0; i < k; i++) {
		Q[i][X] = P[temp + i][X];
		Q[i][Y] = P[temp + i][Y];
		Q[i][Z] = P[temp + i][Z];
	}
	for (r = 1; r < k; r++) {
		for (i = j; i >= temp + r; i--) {
			lamta = (t - T[i]) / (T[i + k - r] - T[i]);
			temp1 = i - temp;
			Q[temp1][X] = lamta*Q[temp1][X] + (1.0 - lamta)*Q[temp1 - 1][X];
			Q[temp1][Y] = lamta*Q[temp1][Y] + (1.0 - lamta)*Q[temp1 - 1][Y];
			Q[temp1][Z] = lamta*Q[temp1][Z] + (1.0 - lamta)*Q[temp1 - 1][Z];
		}
	}
	V[X] = Q[k - 1][X];
	V[Y] = Q[k - 1][Y];
	V[Z] = Q[k - 1][Z];
}

void build_ground(Vector1 cur, Vector1 last, int t) {
	Vector1 c[20];
	Vector1 p1;
	p1[X] = cur[X] - last[X];
	p1[Y] = cur[Y] - last[Y];
	p1[Z] = cur[Z] - last[Z];
	float len = sqrt(p1[X] * p1[X] + p1[Y] * p1[Y] + p1[Z] * p1[Z]);
	p1[X] /= len;
	p1[Y] /= len;
	p1[Z] /= len;
	float sa, ca, sb, cb;
	sb = p1[X];
	cb = sqrt(1 - sb*sb);
	ca = p1[Z] / cb;
	sa = -p1[Y] / cb;
	float mat[6] = { cb, 0,sa*sb,ca,-ca*sb,sa };

	for (int i = 0; i < 12; i++) {
		cir[i][t][X] = origin[i][X] * mat[0] + origin[i][Y] * mat[1] + cur[X];
		cir[i][t][Y] = origin[i][X] * mat[2] + origin[i][Y] * mat[3] + cur[Y];
		cir[i][t][Z] = origin[i][X] * mat[4] + origin[i][Y] * mat[5] + cur[Z];
	}
}

void DisplayBspine(Vector1 P[], float T[], int n, int k, int count) {
	/*for (int q = 0; q < 12;q++) {
	origin[q][X] = 0.07*cos(PI / 6 * q);
	origin[q][Y] = 0.07*sin(PI / 6 * q);
	origin[q][Z] = 0;
	}*/
	Init();
	int i, j;
	float delta, t;
	delta = (T[n + 1] - T[k - 1]) / count;
	t = T[k - 1];
	j = k - 1;
	Bspline(P, T, k, t, j, V[0]);
	for (i = 1; i <= count; i++) {
		t = t + delta;
		while (t > T[j + 1])
			j++;
		Bspline(P, T, k, t, j, V[i]);
		build_ground(V[i], V[i - 1], i);
		for (int l = 0; l < 12; l++) {
			Line(cir[l][i], cir[(l + 1) % 12][i]);
			if (i > 1)
				Line(cir[l][i], cir[l][i - 1]);
		}
		Line(V[i], V[i - 1]);
	}

}

//////////////////////////////////////////
//	B样条结束
//////////////////////////////////////////
Vector pathway[12][1000];
Vector trace[1000];
float cameraR;
int angle;

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

Matrix view_transform(Vector camera, Vector x, Vector y, Vector z) {

#if 0
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
#endif

	Matrix viewT;
	viewT.ViewT(camera.u, camera.v, camera.w);
	Matrix viewR;
	viewR.ViewR(x, y, z);
	//cout << "???" << endl << viewT << viewR << "!!!" << endl;		//调试输出
	Matrix ViewTransform = viewR * viewT;

	return ViewTransform;
}

Matrix project_transform(float z_near, float z_far, int width, int height)
{
	Matrix Project;
	Project.Project(z_near, z_far, width, height);
	return Project;
}

void myDisplay() {
	int index = 0;
	//while (index < 10)
	//{
		glClear(GL_COLOR_BUFFER_BIT);
		Vector4D tmp[12][1000];							//用于坐标变换的四元组
		for (int i = 0; i < 12; ++i)
			for (int j = index; j < index + 50; ++j)
			{
				tmp[i][j] = Vector4D(pathway[i][j]);
				//cout << tmp[i][j].a[0] << ' ' << tmp[i][j].a[1] << ' ' << tmp[i][j].a[2] << ' ' << tmp[i][j].a[3] << endl;
			}


		Vector y0 = pathway[0][index] - trace[index];	//法向
		Vector z = trace[index + 1] - trace[index];		//相机方向
		z.Normalize();
		y0.Normalize();
		Vector x0 = y0.cross(z);
		x0.Normalize();
		Vector y = x0 * cos(angle*PI / 180) + y0 *sin(angle*PI / 180);
		Vector camera = y * 8 + trace[index];			//相机位置
		Vector x = y.cross(z);
		x.Normalize();
		y.Normalize();
		y = y * -1;

		cout << "camera: " << camera.u << ' ' << camera.v << ' ' << camera.w << endl;
		cout << "x: " << x.u << ' ' << x.v << ' ' << x.w << endl;
		cout << "y: " << y.u << ' ' << y.v << ' ' << y.w << endl;
		cout << "z: " << z.u << ' ' << z.v << ' ' << z.w << endl;

		Matrix ViewTransform = view_transform(camera, x, y, z);
		Matrix ProjectTransform = project_transform(20, 200, WIDTH, HEIGHT);

		cout << ViewTransform << ProjectTransform;

		for (int i = 0; i < 12; ++i)
			for (int j = index; j < index + 50; ++j)
			{
				//cout << tmp[i][j].a[0] << ' ' << tmp[i][j].a[1] << ' ' << tmp[i][j].a[2] << ' ' << tmp[i][j].a[3] << endl;
				tmp[i][j] = ViewTransform * tmp[i][j];	//相机坐标系变换没问题！！！
				//cout << tmp[i][j].a[0] << ' ' << tmp[i][j].a[1] << ' ' << tmp[i][j].a[2] << ' ' << tmp[i][j].a[3] << endl;
				tmp[i][j] = ProjectTransform * tmp[i][j];
				tmp[i][j].identity();
				//cout << tmp[i][j].a[0] << ' ' << tmp[i][j].a[1] << ' ' << tmp[i][j].a[2] << ' ' << tmp[i][j].a[3] << endl;
			}

		glColor3b(127, 0, 0);
		for (int i = 0; i < 12; ++i)
		{
			for (int k = index; k < index + 49; ++k)
			{
				int j = (i + 1) % 12;
#if 0
				cout << '(' << tmp[i][k].a[0] << ',' << tmp[i][k].a[1] << ')' << ' ';
				cout << '(' << tmp[i][k + 1].a[0] << ',' << tmp[i][k + 1].a[1] << ')' << ' ';
				cout << '(' << tmp[j][k + 1].a[0] << ',' << tmp[j][k + 1].a[1] << ')' << ' ';
				cout << '(' << tmp[j][k].a[0] << ',' << tmp[j][k].a[1] << ')' << ' ';
				cout << endl << endl;
#endif
				//纹理映射
				glColor3b(127, 127, 127);
				glBegin(GL_POLYGON);
				glVertex2f(tmp[i][k].a[0], tmp[i][k].a[1]);
				glVertex2f(tmp[i][k + 1].a[0], tmp[i][k + 1].a[1]);
				glVertex2f(tmp[j][k + 1].a[0], tmp[j][k + 1].a[1]);
				glVertex2f(tmp[j][k].a[0], tmp[j][k].a[1]);
				glEnd();

				glColor3b(0, 0, 0);
				glBegin(GL_LINE_LOOP);
				glVertex2f(tmp[i][k].a[0], tmp[i][k].a[1]);
				glVertex2f(tmp[i][k + 1].a[0], tmp[i][k + 1].a[1]);
				glVertex2f(tmp[j][k + 1].a[0], tmp[j][k + 1].a[1]);
				glVertex2f(tmp[j][k].a[0], tmp[j][k].a[1]);
				glEnd();
			}

		}
		glFlush();
	//	index++;
	//}
}

void InitEnvironment()
{
	glClearColor(0, 0, 0, 0);
	glClear(GL_COLOR_BUFFER_BIT);
	//glPointSize(1);
	glLineWidth(2);
	//glMatrixMode(GL_MODELVIEW);
	//glLoadIdentity();
	//gluOrtho2D(0, 640, 480, 0);
	angle = 100;
	for (int i = 0; i < 12; ++i)
	{
		trace[i] = Vector(0, i, 0);
		for (int j = 0; j < 200; ++j)
		{
			pathway[i][j] = Vector(10 * sin(i *PI / 6), j, 10 * cos(i *PI / 6));
		}
	}
}

void SpecialKeys(int key, int x, int y)
{
	if (key == GLUT_KEY_LEFT)
		angle += 5;

	if (key == GLUT_KEY_RIGHT)
		angle -= 5;
	if (angle > 360)
		angle -= 360;
	if (angle < 0)
		angle += 360;
	myDisplay();
}

//不用鼠标 
void OnMouse(int button, int state, int x, int y)
{

	if (button == GLUT_LEFT_BUTTON&&state == GLUT_DOWN)	 	//
	{
		return;
	}
	if (button == GLUT_RIGHT_BUTTON&&state == GLUT_DOWN)	//
	{
		return;
	}
}

//还不知道怎么保持比例，没用 

void myReshape(int w, int h)
{

}



int main(int argc, char *argv[])
{
	//srand(time(0));

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);
	glutInitWindowPosition(100, 100);
	glutInitWindowSize(WIDTH, HEIGHT);
	
	glutCreateWindow("Boost v1.0");

	InitEnvironment();

	glutMouseFunc(&OnMouse);
	glutDisplayFunc(&myDisplay);
	glutReshapeFunc(&myReshape);
	glutSpecialFunc(SpecialKeys);

	glutMainLoop();

	return 0;

}