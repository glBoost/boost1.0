#include <GL/glut.h>
#include <vector>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <algorithm>
#include <cmath>
#include <time.h>
using namespace std;

/*  ��Ļ��С  */
#define WIDTH 640
#define HEIGHT 480

#define PI 3.14159265
#define BMP_Header_Length 54	//24λBMP�ļ�ͷ��С
bool game_over = false;
bool game_start = true;
GLuint tex_id = 0;
GLuint tex[10];
GLuint start_pic = 0, over_pic = 0;

/////////////////////////////////////////////
// ��ά������ �� ��Ҫ���ڵ��ʾ����������  //
/////////////////////////////////////////////
struct Vector {
	float u, v, w;
	Vector(float _u, float _v, float _w) :u(_u), v(_v), w(_w) {}
	Vector() :u(0.0), v(0.0), w(0.0) {}

	//�������
	float operator*(const Vector &d)
	{
		return u*d.u + v*d.v + w*d.w;
	}

	//��������
	Vector operator*(float f)
	{
		return Vector(u*f, v*f, w*f);
	}

	//��������
	Vector operator-(const Vector &d)
	{
		return Vector(u - d.u, v - d.v, w - d.w);
	}

	//�����ӷ�
	Vector operator+(const Vector &d)
	{
		return Vector(u + d.u, v + d.v, w + d.w);
	}

	//�������
	Vector cross(Vector V) {
		Vector result;
		result.u = v*V.w - w*V.v;
		result.v = w*V.u - u*V.w;
		result.w = u*V.v - v*V.u;
		return result;
	}

	//����ģ��
	float length()
	{
		return sqrt(u*u + v*v + w*w);
	}

	//������λ��
	void Normalize()
	{
		float l = length();
		u /= l;
		v /= l;
		w /= l;
		return;
	}
};

//////////////////////////////////////
//	��ά�����ࣺ��Ҫ��������任	//
//////////////////////////////////////
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

	//������ά����������ά����
	Vector4D(Vector v) {
		a[0] = v.u;
		a[1] = v.v;
		a[2] = v.w;
		a[3] = 1;
	}

	//[]�������װ
	float& operator[](int i) {
		return a[i];
	}

	//��һ����ת��ΪͶӰ��Ļ����
	void identity() {
		//z = 0�Ĵ����z < 0�Ĳü�
		if (a[3] < 1e-6)
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

//////////////////////////////////////
//	��ά�����ࣺ��Ҫ��������任	//
//////////////////////////////////////
struct Matrix {
	float a[4][4];
	Matrix() {
		memset(a, 0, sizeof(a));
	}

	//ƽ�Ʊ任
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

	//�����������ת�任
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

	//���ű任
	void ModelS(float fx, float fy, float fz) {
		memset(a, 0, sizeof(a));
		a[0][0] = fx;
		a[1][1] = fy;
		a[2][2] = fz;
		a[3][3] = 1;
	}

	//���������
	float* operator[](int x) {
		return a[x];
	}

	//�۲�����ϵƽ�Ʊ任
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

	//��������ϵת��Ϊ�۲�����ϵ
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

	//ͶӰ�任
	void Project(float z_near, float z_far, int width, int height) {
		memset(a, 0, sizeof(a));
		a[0][0] = z_near / (width / 2);
		a[1][1] = z_near / (height / 2);
		a[2][2] = -(z_far + z_near) / (z_far - z_near);
		a[2][3] = 2 * z_far*z_near / (z_far - z_near);
		a[3][2] = 1;
	}

	//����˷�
	Matrix operator*(Matrix b) {
		Matrix result;
		for (int i = 0; i < 4; ++i)
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

	//�������������debug��
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

//����������ĳ˷�
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

//////////////////////////////////
//	����ʱ�����������õ����ݽṹ
//////////////////////////////////
struct Rec
{
	float dis;				//�������ĵ�����ľ���
	Vector4D vertices[4];	//���ε��ĸ�����
	int tex_id;				//���ε�����
};
Rec rec[2000];	//���������Ҫ���ľ���

				//�Ӵ�С����
int cmpRec(const void* p1, const void* p2)
{
	return (*(Rec*)p2).dis > (*(Rec*)p1).dis ? 1 : -1;
}

//��С��������
int cmpInt(const void* p1, const void* p2)
{
	return (*(int*)p2) < (*(int*)p1) ? 1 : -1;
}
////////////////////////////////////
//	B������Ҫ������
////////////////////////////////////
#define X 0
#define Y 1
#define Z 2
#define MM 20
#define random(x) (rand()%x)

typedef float Vector1[3];
GLfloat xRot = 0;
GLfloat yRot = 0;
GLfloat zRot = 0;

Vector1 origin[12];
Vector1 cir[12][2000];
Vector1 points[11];
int N = 10;
int K = 4;
float tt[100];
int seg_count = 2000;
Vector1 V[50000];
int ran[2][2000];

struct Barrier
{
	int b_line, b_index;
	Vector4D vertices[8];
}barriers[2000];

///////////////////////////////
//	����任��Ҫ������
////////////////////////////////
Vector pathway[12][2000];	//12������
Vector trace[2000];			//������
float cameraR;				//�������Բ����İ뾶
int angle;					//�������λ����Ĭ��λ��ƫת�ĽǶ�
float angle_a;
float fangle;

							//////////////////////////////
							//	ͼ�������Ҫ������
							//////////////////////////////
int index = 0;					//��ǰλ������Ļ�׼�±�
float position = 0.0;			//�����ʵ��λ��
float speed = 0.04;				//����ƶ����ٶ�
float offset = 0.0;				//λ�ó�����׼�Ĳ���
float Radio = 0.015;				//����뾶
int b_start = 0, b_end = 0, c_end = 0;	//��Ҫ�����ϰ���Ŀ�ʼ���ꡢ�������ꡢ��Ҫ����·���ϰ���Ľ�������
Vector4D tmp[12][2000];			//��ʱ4ά����������任
Barrier b_tmp[30];				//��ʱ�ϰ�������任

								////////////////////////////////////
								//	���еķ���
								////////////////////////////////////
void make_barriers();
void Init(void);
void Bspline(Vector1 P[], float T[], int k, float t, int j, Vector1 V);
void build_ground(Vector1 cur, Vector1 last, int t);
void DisplayBspine(Vector1 P[], float T[], int n, int k, int seg_seg_count);
Matrix view_transform(Vector camera, Vector x, Vector y, Vector z);
Matrix project_transform(float z_near, float z_far, int width, int height);
int power_of_two(int n);
GLuint load_texture(const char* file_name);
bool hit_barrier(int s, int e);
void GameScreen();
void GameEnd();
void Game();
void InitEnvironment();
void GameInit();
void myDisplay();
void SpecialKeys(int key, int x, int y);
void OnMouse(int button, int state, int x, int y);
void Update();
int main(int argc, char *argv[]);

/////////////////////////////////////////
// B������ģ
//	     by �ľ��u
/////////////////////////////////////////

//��������ϰ�
void make_barriers() {
	srand((int)time(0));
	for (int i = 0; i < (seg_count / 6); i++) {
		ran[0][i] = random(12);
		ran[1][i] = random(seg_count - 15) + 15;
	}
	//���ϰ���index����
	qsort(ran[1], seg_count / 6, sizeof(ran[1][0]), cmpInt);
}

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
	points[1][Z] = 0.2f;

	points[2][X] = -0.5f;
	points[2][Y] = -0.5f;
	points[2][Z] = 1.3f;

	points[3][X] = 0.5f;
	points[3][Y] = -0.5f;
	points[3][Z] = -0.5f;

	points[4][X] = -0.5f;
	points[4][Y] = -0.7f;
	points[4][Z] = -0.5f;

	points[5][X] = -0.5f;
	points[5][Y] = 0.5f;
	points[5][Z] = -0.4f;

	points[6][X] = 0.1f;
	points[6][Y] = 1.2f;
	points[6][Z] = -0.9f;

	points[7][X] = 0.75f;
	points[7][Y] = 0.4f;
	points[7][Z] = -0.5f;

	points[8][X] = 1.2f;
	points[8][Y] = 0.9f;
	points[8][Z] = 0.0f;

	points[9][X] = 0.75f;
	points[9][Y] = 0.5f;
	points[9][Z] = 0.6f;

	points[10][X] = 0.5f;
	points[10][Y] = 0.5f;
	points[10][Z] = 0.5f;

	for (int i = 0; i < (N + K - 5); i++) {
		tt[i + 3] = 1.0 / (N + K - 6) * i;
	}

	for (int i = 0; i < 3; i++) {
		tt[i] = 0;
		tt[N + K - i] = 1;
	}
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
	float mat[6] = { cb, 0, sa*sb, ca, -ca*sb, sa };

	for (int i = 0; i < 12; i++) {
		cir[i][t][X] = origin[i][X] * mat[0] + origin[i][Y] * mat[1] + cur[X];
		cir[i][t][Y] = origin[i][X] * mat[2] + origin[i][Y] * mat[3] + cur[Y];
		cir[i][t][Z] = origin[i][X] * mat[4] + origin[i][Y] * mat[5] + cur[Z];
	}


	if (t == 1) {
		for (int i = 0; i < 12; i++) {
			cir[i][0][X] = cir[i][t][X] - cur[X] + last[X];
			cir[i][0][Y] = cir[i][t][Y] - cur[Y] + last[Y];
			cir[i][0][Z] = cir[i][t][Z] - cur[Z] + last[Z];
		}
	}
}

void DisplayBspine(Vector1 P[], float T[], int n, int k, int seg_seg_count) {
	Init();
	int i, j;
	float delta, t;
	delta = (T[n + 1] - T[k - 1]) / seg_count;
	t = T[k - 1];
	j = k - 1;
	Bspline(P, T, k, t, j, V[0]);
	for (i = 1; i <= seg_count; i++) {
		t = t + delta;
		while (t > T[j + 1])
			j++;
		Bspline(P, T, k, t, j, V[i]);
		build_ground(V[i], V[i - 1], i);
	}
}

#if 0
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
#endif

//�ӽǱ任����ļ���
Matrix view_transform(Vector camera, Vector x, Vector y, Vector z) {

	Matrix viewT;
	viewT.ViewT(camera.u, camera.v, camera.w);	//���ƽ�Ƶ�ԭ��
	Matrix viewR;
	viewR.ViewR(x, y, z);						//����ϵ��һ
	Matrix ViewTransform = viewR * viewT;
	return ViewTransform;
}

//ͶӰ�任����ļ���
Matrix project_transform(float z_near, float z_far, int width, int height)
{
	Matrix Project;
	Project.Project(z_near, z_far, width, height);
	return Project;
}

//���һ�������Ƿ�Ϊ2�������η�������ǣ�����1�����򷵻�0
int power_of_two(int n)
{
	if (n <= 0)
		return 0;
	return (n & (n - 1)) == 0;
}

/////////////////////////////////////////////////
//	��ȡһ��BMP�ļ���Ϊ����
//	���ʧ�ܣ�����0������ɹ�������������
//////////////////////////////////////////////////
GLuint load_texture(const char* file_name)
{
	GLint width, height, total_bytes;
	GLubyte* pixels = NULL;
	GLuint last_texture_ID, texture_ID = 0;

	// ���ļ������ʧ�ܣ�����
	FILE* pFile = fopen(file_name, "rb");
	if (pFile == 0)
		return 0;

	// ��ȡ�ļ���ͼ��Ŀ�Ⱥ͸߶�
	fseek(pFile, 0x0012, SEEK_SET);
	fread(&width, 4, 1, pFile);
	fread(&height, 4, 1, pFile);
	fseek(pFile, BMP_Header_Length, SEEK_SET);

	// ����ÿ��������ռ�ֽ����������ݴ����ݼ����������ֽ���
	{
		GLint line_bytes = width * 3;
		while (line_bytes % 4 != 0)
			++line_bytes;
		total_bytes = line_bytes * height;
	}

	// �����������ֽ��������ڴ�
	pixels = (GLubyte*)malloc(total_bytes);
	if (pixels == 0)
	{
		fclose(pFile);
		return 0;
	}

	// ��ȡ��������
	if (fread(pixels, total_bytes, 1, pFile) <= 0)
	{
		free(pixels);
		fclose(pFile);
		return 0;
	}

	// �ھɰ汾��OpenGL��
	// ���ͼ��Ŀ�Ⱥ͸߶Ȳ��ǵ������η�������Ҫ��������
	// ���ﲢû�м��OpenGL�汾�����ڶ԰汾�����ԵĿ��ǣ����ɰ汾����
	// ���⣬�����Ǿɰ汾�����°汾��
	// ��ͼ��Ŀ�Ⱥ͸߶ȳ�����ǰOpenGLʵ����֧�ֵ����ֵʱ��ҲҪ��������
	/*
	{
	GLint max;
	glGetIntegerv(GL_MAX_TEXTURE_SIZE, &max);
	if (!power_of_two(width)
	|| !power_of_two(height)
	|| width > max
	|| height > max)
	{
	const GLint new_width = 256;
	const GLint new_height = 256; // �涨���ź��µĴ�СΪ�߳���������
	GLint new_line_bytes, new_total_bytes;
	GLubyte* new_pixels = 0;

	// ����ÿ����Ҫ���ֽ��������ֽ���
	new_line_bytes = new_width * 3;
	while (new_line_bytes % 4 != 0)
	++new_line_bytes;
	new_total_bytes = new_line_bytes * new_height;

	// �����ڴ�
	new_pixels = (GLubyte*)malloc(new_total_bytes);
	if (new_pixels == 0)
	{
	free(pixels);
	fclose(pFile);
	return 0;
	}

	// ������������
	gluScaleImage(GL_RGB,
	width, height, GL_UNSIGNED_BYTE, pixels,
	new_width, new_height, GL_UNSIGNED_BYTE, new_pixels);

	// �ͷ�ԭ�����������ݣ���pixelsָ���µ��������ݣ�����������width��height
	free(pixels);
	pixels = new_pixels;
	width = new_width;
	height = new_height;
	}
	}
	*/

	// ����һ���µ�������
	glGenTextures(1, &texture_ID);
	if (texture_ID == 0)
	{
		free(pixels);
		fclose(pFile);
		return 0;
	}

	// ���µ������������������������
	glBindTexture(GL_TEXTURE_2D, texture_ID);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0,
		GL_BGR_EXT, GL_UNSIGNED_BYTE, pixels);

	free(pixels);
	return texture_ID;
}

bool hit_barrier(int s, int e)
{
	for (int i = s; i < e; i++)
	{
		if (barriers[i].b_index == index)
		{
			if (barriers[i].b_line == ((angle + 260) % 360) / 30)
				return true;
			if (barriers[i].b_line == ((angle + 280) % 360) / 30)
				return true;
		}
		if (barriers[i].b_index > index)
			return false;
	}
}

void GameScreen() {
	glClear(GL_COLOR_BUFFER_BIT);
	glBindTexture(GL_TEXTURE_2D, start_pic);
	glBegin(GL_QUADS);
	glTexCoord2f(0.0f, 0.0f); glVertex2f(-1.0f, -1.0f);
	glTexCoord2f(0.0f, 1.0f); glVertex2f(-1.0f, 1.0f);
	glTexCoord2f(1.0f, 1.0f); glVertex2f(1.0f, 1.0f);
	glTexCoord2f(1.0f, 0.0f); glVertex2f(1.0f, -1.0f);
	glEnd();
	glFlush();
	glutSwapBuffers();
}

void GameEnd() {
	glClear(GL_COLOR_BUFFER_BIT);
	glBindTexture(GL_TEXTURE_2D, over_pic);
	glBegin(GL_QUADS);
	glTexCoord2f(0.0f, 0.0f); glVertex2f(-1.0f, -1.0f);
	glTexCoord2f(0.0f, 1.0f); glVertex2f(-1.0f, 1.0f);
	glTexCoord2f(1.0f, 1.0f); glVertex2f(1.0f, 1.0f);
	glTexCoord2f(1.0f, 0.0f); glVertex2f(1.0f, -1.0f);
	glEnd();
	glFlush();
	glutSwapBuffers();
}

void Game() {
	//����
	glClear(GL_COLOR_BUFFER_BIT);
	glBindTexture(GL_TEXTURE_2D, tex_id);

	//������ʼ��
	for (int i = 0; i < 12; ++i)
		for (int j = index; j < index + 40; ++j)
		{
			tmp[i][j] = Vector4D(pathway[i][j]);
		}

	//�������λ�ü��������ϵ��������
	Vector y0 = pathway[0][index] - trace[index];	//������ת����ϵ��y��
	Vector z = trace[index + 1] - trace[index];		//�������	
	y0.Normalize();
	Vector x0 = y0.cross(z);						//��ת����ϵ��x������
	x0.Normalize();
	Vector y = x0 * cos(angle*PI / 180) + y0 *sin(angle*PI / 180);	//��ʵ���λ���������������ڵķ���
	y.Normalize();
	Vector camera = y * Radio + trace[index];			//���λ�� = ����λ�� + �뾶
	camera = camera + z * offset;						//					  + ƫ��
	z.Normalize();
	Vector x = y.cross(z);							//�������ϵ��x��
	x.Normalize();
	y = y * -1;										//�������ϵ��y��

													/*cout << "camera: " << camera.u << ' ' << camera.v << ' ' << camera.w << endl;
													cout << "x: " << x.u << ' ' << x.v << ' ' << x.w << endl;
													cout << "y: " << y.u << ' ' << y.v << ' ' << y.w << endl;
													cout << "z: " << z.u << ' ' << z.v << ' ' << z.w << endl;*/

													//����任����
	Matrix ViewTransform = view_transform(camera, x, y, z);		//�������ϵ�任����
	Matrix ProjectTransform = project_transform(100, 500, WIDTH, HEIGHT);	//ͶӰ�任����

																			/*cout << ViewTransform << ProjectTransform;*/

																			//�ҵ���Ҫ�����ϰ���
	while (barriers[b_start].b_index < index)
		b_start++;
	while (barriers[b_end].b_index <= index + 39)
		b_end++;
	while (barriers[c_end].b_index <= index + 59)
		c_end++;

	if (hit_barrier(b_start, b_end))
		game_over = true;

	//�����ܱڵ���������ϵ��ת��
	for (int i = 0; i < 12; ++i)
		for (int j = index; j < index + 40; ++j)
			tmp[i][j] = ViewTransform * tmp[i][j];
	//�ɼ��ϰ��ﶥ����������ϵ��ת��
	for (int i = b_start; i < b_end; ++i)
	{
		b_tmp[i - b_start].b_line = barriers[i].b_line;
		b_tmp[i - b_start].b_index = barriers[i].b_index;
		for (int j = 0; j < 8; ++j)
			b_tmp[i - b_start].vertices[j] = ViewTransform * barriers[i].vertices[j];
	}

	int r_num = 0;
	//�����㷨�����ݾ����Զ���Ѿ�������
	for (int i = 0; i < 12; ++i)
		for (int j = index; j < index + 39; ++j)
		{
			//����ܱ��Ƿ�Ϳ��ɫ
			bool is_color = false;
			for (int k = b_start; k < c_end; ++k)
				if (i == barriers[k].b_line)
					if (j<barriers[k].b_index && (j + 21)>barriers[k].b_index)
						is_color = true;

			if (is_color)
				rec[r_num].tex_id = 1 + i / 2;
			else
			{
				if (j < 500)
					rec[r_num].tex_id = 0;
				else if (j < 600)
					rec[r_num].tex_id = 7;
				else
					rec[r_num].tex_id = 8;
			}
				

			rec[r_num].vertices[0] = tmp[i][j];
			rec[r_num].vertices[1] = tmp[i][j + 1];
			rec[r_num].vertices[2] = tmp[(i + 1) % 12][j + 1];
			rec[r_num].vertices[3] = tmp[(i + 1) % 12][j];
			Vector4D* p1 = &tmp[i][j];
			Vector4D* p2 = &tmp[(i + 1) % 12][j + 1];
			rec[r_num++].dis = (p1->a[0] + p2->a[0])*(p1->a[0] + p2->a[0]) + (p1->a[1] + p2->a[1])*(p1->a[1] + p2->a[1]) + (p1->a[2] + p2->a[2])*(p1->a[2] + p2->a[2]);
		}
	for (int i = 0; i < b_end - b_start; ++i)
	{
		int line_no = b_tmp[i].b_line;
		for (int j = 0; j < 3; ++j)
		{
			rec[r_num].tex_id = 1 + line_no / 2;
			rec[r_num].vertices[0] = b_tmp[i].vertices[j];
			rec[r_num].vertices[1] = b_tmp[i].vertices[j + 4];
			rec[r_num].vertices[2] = b_tmp[i].vertices[(j + 3) % 4 + 4];
			rec[r_num].vertices[3] = b_tmp[i].vertices[(j + 3) % 4];
			Vector4D* p1 = &b_tmp[i].vertices[j];
			Vector4D* p2 = &b_tmp[i].vertices[(j + 3) % 4 + 4];
			rec[r_num++].dis = (p1->a[0] + p2->a[0])*(p1->a[0] + p2->a[0]) + (p1->a[1] + p2->a[1])*(p1->a[1] + p2->a[1]) + (p1->a[2] + p2->a[2])*(p1->a[2] + p2->a[2]);
		}
		rec[r_num].tex_id = 1 + line_no / 2;
		rec[r_num].vertices[0] = b_tmp[i].vertices[4];
		rec[r_num].vertices[1] = b_tmp[i].vertices[5];
		rec[r_num].vertices[2] = b_tmp[i].vertices[6];
		rec[r_num].vertices[3] = b_tmp[i].vertices[7];
		Vector4D* p1 = &b_tmp[i].vertices[4];
		Vector4D* p2 = &b_tmp[i].vertices[6];
		rec[r_num++].dis = (p1->a[0] + p2->a[0])*(p1->a[0] + p2->a[0]) + (p1->a[1] + p2->a[1])*(p1->a[1] + p2->a[1]) + (p1->a[2] + p2->a[2])*(p1->a[2] + p2->a[2]);
	}
	qsort(rec, r_num, sizeof(rec[0]), cmpRec);

	//�ܱڵ����ϰ��ﶥ��ͶӰ�任����Ļ
	for (int r = 0; r < r_num; ++r)
		for (int i = 0; i < 4; ++i)
		{
			rec[r].vertices[i] = ProjectTransform * rec[r].vertices[i];
			rec[r].vertices[i].identity();
		}

	glColor3b(127, 0, 0);

	//����
	for (int r = 0; r < r_num; ++r)
	{
		//����ӳ��
		glBindTexture(GL_TEXTURE_2D, tex[rec[r].tex_id]);
		glBegin(GL_QUADS);
		glTexCoord2f(0.0f, 0.0f); glVertex2f(rec[r].vertices[0].a[0], rec[r].vertices[0].a[1]);
		glTexCoord2f(0.0f, 1.0f); glVertex2f(rec[r].vertices[1].a[0], rec[r].vertices[1].a[1]);
		glTexCoord2f(1.0f, 1.0f); glVertex2f(rec[r].vertices[2].a[0], rec[r].vertices[2].a[1]);
		glTexCoord2f(1.0f, 0.0f); glVertex2f(rec[r].vertices[3].a[0], rec[r].vertices[3].a[1]);
		glEnd();
	}
	glFlush();
	glutSwapBuffers();
}

void GameInit() {
	position = 0;
	speed = 0.04;
	b_start = b_end = c_end = 0;
	angle = 0;
	fangle = 0;
	angle_a = 0;
	make_barriers();													//�����ϰ�

																		//����b�������߳�ʼ���ܵ�����
	for (int j = 0; j < seg_count; ++j)
	{
		for (int i = 0; i < 12; ++i) {
			//	pathway[i][j] = Vector(10 * sin(i *PI / 6), j, 10 * cos(i *PI / 6));
			pathway[i][j] = Vector(cir[i][j][X], cir[i][j][Y], cir[i][j][Z]);
		}
		trace[j] = Vector(V[j][X], V[j][Y], V[j][Z]);
		//trace[j] = Vector(0, j, 0);
	}

	//ת���ϰ��ﶥ��
	for (int u = 0; u < (seg_count / 6); u++) {
		barriers[u].b_line = ran[0][u];
		barriers[u].b_index = ran[1][u];
		int bl = barriers[u].b_line;
		int bi = barriers[u].b_index;
		barriers[u].vertices[0] = Vector4D(pathway[bl][bi]);
		barriers[u].vertices[1] = Vector4D(pathway[(bl + 1) % 12][bi]);
		barriers[u].vertices[2] = Vector4D(pathway[(bl + 1) % 12][bi + 1]);
		barriers[u].vertices[3] = Vector4D(pathway[bl][bi + 1]);
		//Vector b_height = (pathway[(bl + 1) % 12][bi] - pathway[bl][bi]).cross(pathway[bl][bi + 1] - pathway[bl][bi]);
		Vector b_height = (Vector(V[bi][X], V[bi][Y], V[bi][Z]) - ((pathway[(bl + 1) % 12][bi] + pathway[bl][bi])*0.5f))*0.5f;
		//Vector b_height = (Vector(V[bi][X], V[bi][Y], V[bi][Z]) - pathway[(bl + 1) % 12][bi])*0.3f;
		//Vector b_height = (pathway[bl][bi + 1] - pathway[bl][bi]).cross(pathway[(bl + 1) % 12][bi] - pathway[bl][bi]);
		barriers[u].vertices[4] = Vector4D(pathway[bl][bi] + b_height);
		barriers[u].vertices[5] = Vector4D(pathway[(bl + 1) % 12][bi] + b_height);
		barriers[u].vertices[6] = Vector4D(pathway[(bl + 1) % 12][bi + 1] + b_height);
		barriers[u].vertices[7] = Vector4D(pathway[bl][bi + 1] + b_height);
	}
}

void myDisplay() {
	if (game_start)
		GameScreen();
	else if (!game_over)
		Game();
	else
		GameEnd();
}

void InitEnvironment()
{
	glClearColor(0, 0, 0, 0);
	glClear(GL_COLOR_BUFFER_BIT);
	glEnable(GL_TEXTURE_2D);											//��������

																		//��ȡ����
	tex[0] = load_texture("wenli2.bmp");
	tex[1] = load_texture("wenli2_1.bmp");
	tex[2] = load_texture("wenli2_2.bmp");
	tex[3] = load_texture("wenli2_3.bmp");
	tex[4] = load_texture("wenli2_4.bmp");
	tex[5] = load_texture("wenli2_5.bmp");
	tex[6] = load_texture("wenli2_6.bmp");
	tex[7] = load_texture("wenli3.bmp");
	tex[8] = load_texture("wenli1.bmp");
	start_pic = load_texture("start.bmp");
	over_pic = load_texture("over.bmp");
	//cout << tex_id << ' ' << start_pic << ' ' << over_pic << endl;
	game_start = true;
	game_over = true;

	glLineWidth(2);
	angle = 0;															//��ʼ���ĽǶ�
	DisplayBspine(points, tt, N, K, seg_count);							//b�����γɹܵ�
	GameInit();
	return;
}

//�����¼�
void SpecialKeys(int key, int x, int y)
{
	if (key == GLUT_KEY_LEFT)
		fangle += 5;
	if (key == GLUT_KEY_RIGHT)
		fangle -= 5;
	if (angle >= 360)
		fangle -= 360;
	if (angle < 0)
		fangle += 360;

	myDisplay();
}

//��� 
void OnMouse(int button, int state, int x, int y)
{
	if (button == GLUT_LEFT_BUTTON&&state == GLUT_DOWN)	 	//
	{
		if (game_start)
		{

			if (x > 250 && x < 390 && y > 250 && y < 300)
			{
				GameInit();
				game_start = false;
				game_over = false;
			}
			else if (x > 250 && x < 390 && y > 350 && y < 400)
			{
				exit(0);
			}
		}
		else if (game_over)
		{
			if (x > 90 && x < 240 && y > 400 && y < 450)
			{
				GameInit();
				game_over = false;
			}
			else if (x > 400 && x < 550 && y > 400 && y < 450)
			{
				game_start = true;
				game_over = true;
			}
		}
		else
		{
			angle_a = -(float)((x - 320) / 200.0);
		}
		return;
	}
	if (button == GLUT_RIGHT_BUTTON&&state == GLUT_DOWN)	//
	{
		return;
	}
	if (button == GLUT_LEFT_BUTTON && state == GLUT_UP)
	{
		if (!game_start && !game_over)
		{
			angle_a = 0.0;
		}
	}
}

//ˢ�º���
void Update()
{
	if (game_start)
		GameScreen();
	else if (position < 1500 && !game_over)
	{
		/*if (position < 500)
			tex_id = tex1;
		else if (position < 600)
			tex_id = tex3;
		else
			tex_id = tex2;*/
		tex_id = tex[0];
		position += speed;
		fangle += angle_a;
		if (fangle >= 360)
			fangle -= 360;
		if (fangle < 0)
			fangle += 360;
		angle = (int)fangle;
		index = (int)position;
		offset = position - (float)index;
		glutPostRedisplay();
		//Sleep(10);
	}
	else
		GameEnd();
}

int main(int argc, char *argv[])
{
	//srand(time(0));

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
	glutInitWindowPosition(100, 100);
	glutInitWindowSize(WIDTH, HEIGHT);

	glutCreateWindow("Boost v1.0");

	InitEnvironment();

	glutIdleFunc(&Update);				//CPU��ʱ����
	glutMouseFunc(&OnMouse);			//��꺯��
	glutDisplayFunc(&myDisplay);		//��ʾ����
	//glutReshapeFunc(&myReshape);		//���ڸı亯��
	glutSpecialFunc(SpecialKeys);		//���̺���
	glutMainLoop();

	return 0;
}