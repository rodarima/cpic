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
#define FPS 60
#define FRAME_PERIOD ((int) floor(1e9/FPS))

HMGL gr;
HMDT dat;
double yy[N], shift = 0.0;
GLenum doubleBuffer = GL_FALSE;
GLint windW = 800, windH = 300;
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
	mgl_plot(gr, (HCDT) dat, "b", "");
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
	int i;
	struct timespec now;
	//printf("idle called\n");
	for(i=0; i<N; i++)
	{
		//yy[i] = ((double) rand()) / RAND_MAX;
		yy[i] = cos((double) i/N * 10.0 + shift);
	}
	shift += 0.1;
	mgl_data_set_double(dat, yy, N, 1, 1);

	clock_gettime(CLOCK_MONOTONIC, &now);

	if((now.tv_sec >= later.tv_sec) && (now.tv_nsec >= later.tv_nsec))
	{
		later.tv_sec = now.tv_sec + 1;
		later.tv_nsec = now.tv_nsec;

		fprintf(stderr, "FPS = %d\n", frames);
		frames = 0;
	}

	/* Before calling the redraw function, we can sleep here */
	clock_nanosleep(CLOCK_MONOTONIC, TIMER_ABSTIME, &next_frame, NULL);

	glutSetWindow(win1);
	glutPostRedisplay();
}

void
reload(void *p)
{
	//printf("reload called\n");
}

void
Reshape(int width, int height)
{
	//windW = width;
	//windH = height;

	//glViewport(0, 0, windW, windH);

	//glMatrixMode(GL_PROJECTION);
	//glLoadIdentity();
	//gluOrtho2D(-0.5, windW + 0.5, -0.5, windH + 0.5);
	//glMatrixMode(GL_MODELVIEW);

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

	dat = mgl_create_data_size(N, 1, 1);

	gr = mgl_create_graph_gl();
	//mgl_create_graph_glut(draw, "title", (void *) dat, reload);

	glutIdleFunc(idle);

	glutMainLoop();

	return 0;
}
