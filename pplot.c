#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <GL/glut.h>

#include "specie.h"

GLenum doubleBuffer;
GLint windW = 500, windH = 500;

#define MAX_HIST 30
#define TAG 'p'
#define MAX_LINE 256

particle_t *particles[MAX_HIST];
int hist = 0;

int nparticles;
int shape;
float dx, dt;
int play = 1;
int clear = 0;
int plotting = 1;

static int
get_line(char *buf, size_t n)
{
	while(getline(&buf, &n, stdin))
	{
		if(buf[0] == TAG)
			return 0;

		printf("%s", buf);
		fflush(stdout);
	}
	return -1;
}

static void
Init(void)
{
	int i;
	particle_t *p;
	char buf[MAX_LINE];

	if(get_line(buf, MAX_LINE))
	{
		fprintf(stderr, "EOF before read config\n");
		exit(1);
	}

	sscanf(buf, "p %d %d %f %f", &nparticles, &shape, &dx, &dt);
	printf("Number of particles=%d, shape=%d\n", nparticles, shape);

	for(i = 0; i < MAX_HIST; i++)
		particles[i] = calloc(sizeof(particle_t), nparticles);

	glClearColor(0.0, 0.0, 0.0, 0.0);

	glDisable(GL_DITHER);
}

void
Reshape(int width, int height)
{
	windW = width;
	windH = height;

	glViewport(0, 0, windW, windH);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(-0.5, windW + 0.5, -0.5, windH + 0.5);
	glMatrixMode(GL_MODELVIEW);
}

/* ARGSUSED1 */
static void
Key(unsigned char key, int x, int y)
{
	switch (key) {
	case 'p':
		play = !play;
		break;
	case 'c':
		clear = !clear;
		break;
	case 'n':
		plotting = !plotting;
		break;
	case 'q':
	case 27:
		exit(0);
	}
}

void
Idle(void)
{
	float t, x, u;
	int i, j, ret;
	particle_t *p;
	char buf[MAX_LINE];

	if(!play)
	{
		glutPostRedisplay();
		return;
	}

	hist = (hist + 1) % MAX_HIST;

	for(i = 0; i < nparticles; i++)
	{
		/* Copy old particle into particles1 */
		//memcpy(&particles1[prev_hist][i], &particles[hist][i], sizeof(particle_t));

		if(get_line(buf, MAX_LINE))
		{
			fprintf(stderr, "EOF reached\n");
			exit(0);
		}

		ret = sscanf(buf, "p %d %f %f",
				&j, &x, &u);

		if(ret != 3)
		{
			printf("ret=%d i=%d j=%d, exitting\n", ret, i, j);
			exit(1);
		}

		p = &particles[hist][j];
		p->x = x;
		p->u = u;
	}

	glutPostRedisplay();
}

static int
get_curve(float *xx, float *yy, int *segment, int pi)
{
	particle_t *p;
	int from_hist = (hist + 1) % MAX_HIST;
	float x0 = 0, y0 = 0, x1, y1;
	int i, si = 0, sl = 0;

	for(i = 0; i < MAX_HIST; i++, from_hist = (from_hist + 1) % MAX_HIST, sl++)
	{
		p = &particles[from_hist][pi];

		x1 = (p->x / (dx * shape)) * windW;
		y1 = (p->u / (8.0 * 3e8)) * windH;

		/* Center y, as u goes from about -c to +c */
		y1 += windH / 2.0;

		if(i > 0 && fabs(x1 - x0) > windW/4) // Wrap around the screen
		{
			/* Create another segment */
			segment[si++] = sl;
			sl = 0;
		}

//		/* For plotting purposes, we need at least one pixel drawn */
//		if(i > 0 && fabs(x0 - x1) < 1.0)
//		{
//			x1 = x0 + 1.0;
//		}

		xx[i] = x1;
		yy[i] = y1;

		x0 = x1;
		y0 = y1;
	}

	segment[si++] = sl;

	return si;
}

static void
plot_particle(int pi)
{
	float x[MAX_HIST], y[MAX_HIST];
	int segments[MAX_HIST];
	int i, j, k=0, n;
	float cm, ch, cl;
	float h = 1.0, l = 0.3;

	n = get_curve(x, y, segments, pi);

	glLineWidth(2.0);

	for(i = 0; i < n; i++)
	{

		glBegin(GL_LINE_STRIP);
		for(j = 0; j<segments[i]; j++)
		{
			cm = (float) k / MAX_HIST;
			if(cm < 0.3) cm = 0.3;
			ch = cm * h;
			cl = cm * l;
			if(pi % 2)
				glColor3f(ch, cl, cl);
			else
				glColor3f(cl, ch, cl);
			glVertex2f(x[k], y[k]);
			k++;
		}
		glEnd();
	}
}


static void
plot_particles()
{
	int i;

	if(clear)
		glClear(GL_COLOR_BUFFER_BIT);

	for(i=0; i<nparticles; i++)
	{
		plot_particle(i);
	}
}

void
Display(void)
{
	if(!plotting)
		return;

	plot_particles();

	if (doubleBuffer) {
		glutSwapBuffers();
	} else {
		glFlush();
	}
}

void
Visible(int state)
{
	if (state == GLUT_VISIBLE) {
		glutIdleFunc(Idle);
	} else {
		glutIdleFunc(NULL);
	}
}

int
main(int argc, char **argv)
{
	GLenum type;

	glutInitWindowSize(windW, windH);
	glutInit(&argc, argv);

	type = GLUT_RGB;
	type |= (doubleBuffer) ? GLUT_DOUBLE : GLUT_SINGLE;
	glutInitDisplayMode(type);
	glutCreateWindow("plot");
	//glutCreateWindow("two streams");
	//glutFullScreen();

	Init();

	glutReshapeFunc(Reshape);
	glutKeyboardFunc(Key);
	glutVisibilityFunc(Visible);
	glutDisplayFunc(Display);
	glutMainLoop();
	return 0;             /* ANSI C requires main to return int. */
}
