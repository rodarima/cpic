#include <stdbool.h>
#include <mgl2/mgl_cf.h>
#include <mgl2/glut.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <GL/glut.h>

/* FIXME: If the title is called, nothing is shown on the plot */
#define SHOW_TITLE 0

#define N 100
#define FPS 30
#define FRAME_PERIOD ((int) floor(1e9/FPS))

HMGL gr;
HMDT dat, du, dv;
double zz[N*N], uu[N*N], vv[N*N], shift = 0.0;
GLenum doubleBuffer = GL_FALSE;
GLint windW = 800, windH = 800;
int win1;
int frames = 0;
struct timespec later = {0, 0};
struct timespec next_frame = {0, 0};

void
draw()
{
	mgl_clf(gr);
	//printf("draw called\n");

#if SHOW_TITLE
	mgl_title(gr, "This title", "", 10);
#endif
	mgl_axis(gr, "xy", "", "");
	mgl_box(gr);
	//mgl_title(gr, "This title", "#", 2.0);
	//mgl_axis_grid(gr, "xy", "", "");
	//mgl_plot(gr, (HCDT) dat, "b", "");
	//mgl_contf(gr, (HCDT) dat, "", "");
	//mgl_colorbar(gr, ">");
	//mgl_dens(gr, (HCDT) dat, 0, 0);
	//mgl_grad(gr, (HCDT) dat, 0, 0);
	mgl_set_meshnum(gr, 30);
	mgl_vect_2d(gr,du,dv,"b2", "");
	//mgl_vect_2d(gr,du,dv,"BbcyrR2", "");
	//mgl_dew_2d(gr,du,dv,0,"");
	//mgl_flow_2d(gr,du,dv,"","");
	mgl_finish(gr);
	glFlush();
	//glFinish();

	frames++;
	clock_gettime(CLOCK_MONOTONIC, &next_frame);
	next_frame.tv_nsec += FRAME_PERIOD;
	if(next_frame.tv_nsec >= 1000000000)
	{
		next_frame.tv_nsec -= 1000000000;
		next_frame.tv_sec += 1;
	}
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

		fprintf(stderr, "FPS = %d\n", frames);
		frames = 0;
	}

	/* Before calling the redraw function, we can sleep here */
	//clock_nanosleep(CLOCK_MONOTONIC, TIMER_ABSTIME, &next_frame, NULL);

	glutSetWindow(win1);
	glutPostRedisplay();
}

void timer(int value)
{
	//glutPostRedisplay();
	glutTimerFunc(1000.0/FPS, timer, 0);
	idle();
}

void
reload(void *p)
{
	//printf("reload called\n");
}

void
Reshape(int width, int height)
{
//	windW = width;
//	windH = height;
//
//	glViewport(0, 0, windW, windH);
//
//	glMatrixMode(GL_PROJECTION);
//	glLoadIdentity();
//	gluOrtho2D(-0.5, windW + 0.5, -0.5, windH + 0.5);
//	glMatrixMode(GL_MODELVIEW);

	mgl_clf(gr);
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
	mgl_set_font_size(gr, 2.0);
	//mgl_load_font(gr, "DejaVu Sans", NULL);
	//mgl_create_graph_glut(draw, "title", (void *) dat, reload);

	//glutIdleFunc(idle);
	timer(0);

	glutMainLoop();

	return 0;
}
