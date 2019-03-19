#include <stdbool.h>
#include <mgl2/mgl_cf.h>
#include <mgl2/glut.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <GL/glut.h>

#include "mgl-glue.h"

#define N 100
#define FPS 30
#define FRAME_PERIOD ((int) floor(1e9/FPS))

HMGL gr;
HMDT dat, du, dv;

double zz[N*N], uu[N*N], vv[N*N], shift = 0.0;
GLenum doubleBuffer = GL_FALSE;
GLint windW = 1000, windH = 500;

int win1;
int frames = 0, updates = 0;
struct timespec later = {0, 0};

void
draw()
{
	mgl_clf(gr);
	//printf("draw called\n");

	//mgl_set_size(gr, windW, windH);

	mgl_subplot(gr, 2, 1, 0, "");
	mgl_title(gr, "Electric potential \\phi", "", 5.0);
	//mgl_dens(gr, (HCDT) dat, 0, 0);
	mgl_contf(gr, (HCDT) dat, "", "");
	mgl_axis_grid(gr, "xy", "", "");

	mgl_subplot(gr, 2, 1, 1, "");

	mgl_title(gr, "Electric field \\textbf{E}", "", 5.0);
//	mgl_axis(gr, "xy", "", "");
//	mgl_box(gr);
	mgl_axis_grid(gr, "xy", "", "");
	//mgl_plot(gr, (HCDT) dat, "b", "");
	//mgl_contf(gr, (HCDT) dat, "", "");
	//mgl_colorbar(gr, ">");
	//mgl_dens(gr, (HCDT) dat, 0, 0);
	//mgl_grad(gr, (HCDT) dat, 0, 0);
	mgl_set_meshnum(gr, 20);
	mgl_vect_2d(gr,du,dv,"b2", "");
	//mgl_vect_2d(gr,du,dv,"BbcyrR2", "");
	//mgl_dew_2d(gr,du,dv,0,"");
	//mgl_flow_2d(gr,du,dv,"","");

	mgl_finish(gr);
	glFlush();
	//glFinish();

	frames++;
}

void
idle()
{
	int i, j;
	double x, y;
	struct timespec now;
	//printf("idle called\n");
	for(i=0; i<N; i++)
	{
		for(j=0; j<N; j++)
		{
			x = (double) i / N * 3.0;
			y = (double) j / N * 3.0;
			zz[i*N + j] = sin(shift) + cos(x) * cos(y);
			uu[i*N + j] = 0.6*sin(2*x+shift)*sin(3*y) + 0.4*cos(3*(x*y));
			vv[i*N + j] = 0.6*cos(2*x+shift)*cos(3*y) + 0.4*cos(3*(x*y));
			//cos((double) i/N * 10.0 + shift);
		}
	}

	updates++;
	shift += 0.05;
	mgl_data_set_double(dat, zz, N, N, 1);
	mgl_data_set_double(du, uu, N, N, 1);
	mgl_data_set_double(dv, vv, N, N, 1);

	clock_gettime(CLOCK_MONOTONIC, &now);

	if((now.tv_sec > later.tv_sec) ||
		((now.tv_sec > later.tv_sec) && (now.tv_nsec >= later.tv_nsec)))
	{
		later.tv_sec = now.tv_sec + 1;
		//later.tv_nsec = now.tv_nsec;

		fprintf(stderr, "FPS = %d, UPS = %d\n", frames, updates);
		frames = 0;
		updates = 0;
	}

	glutSetWindow(win1);
	glutPostRedisplay();
}

void timer(int value)
{
	glutTimerFunc(1000.0/FPS, timer, 0);
	idle();
}

void
Reshape(int width, int height)
{
	windW = width;
	windH = height;

	fprintf(stderr, "Window size is now %dx%d\n", width, height);

	glViewport(0, 0, windW, windH);

	// FIXME: Neither this one
	mgl_set_size(gr, windW, windH);
}

/* ARGSUSED1 */
static void
Key(unsigned char key, int x, int y)
{
	switch (key) {
	case 'q':
	case 27:
		exit(0);
	}
}

int main(int argc, char *argv[])
{
	int i;
	GLenum type;

	glutInitWindowSize(windW, windH);
	glutInit(&argc, argv);

	type = GLUT_RGB;
	type |= (doubleBuffer) ? GLUT_DOUBLE : GLUT_SINGLE;
	glutInitDisplayMode(type);

	glDisable(GL_DITHER);

	win1 = glutCreateWindow("plot particles");
	glutPositionWindow(5, 5);
	glutReshapeFunc(Reshape);
	glutKeyboardFunc(Key);
	glutDisplayFunc(draw);
	glClearColor(0.0, 0.0, 0.0, 0.0);

	mgl_set_num_thr(1);

	dat = mgl_create_data_size(N, N, 1);
	du = mgl_create_data_size(N, N, 1);
	dv = mgl_create_data_size(N, N, 1);

	gr = mgl_create_graph_gl();

	//gr = canvas_create();

	// FIXME: This doesn't set the size properly
	mgl_set_size(gr, windW, windH);

	fprintf(stderr, "Canvas size %dx%d, window size %dx%d\n",
		mgl_get_width(gr),
		mgl_get_height(gr),
		windW, windH);

	mgl_set_font_size(gr, 2.0);

	//glutIdleFunc(idle);
	timer(0);

	glutMainLoop();

	return 0;
}
