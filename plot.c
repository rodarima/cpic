#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <GL/glut.h>
#include <fftw3.h>
#include <unistd.h>

#include "specie.h"

GLenum doubleBuffer = GL_FALSE;
GLint windW = 800, windH = 300;
int win1, win2, win3, win4;

#define MAX_HIST 10
#define MAX_LINE 256
#define MAX_V 20.0
#define MAX_POS 819
#define DT 5.0e-4

particle_t *particles[MAX_HIST];
int hist = 0;

double pos_vec[MAX_POS];
double pos_fft[MAX_POS];
int pos_i = 0;
int redraw_fft = 0, recompute_fft = 0;
double maxfft = 0.0;
double last_freq_peak = 0.0;
double freq_peak = 0.0;
double max_freq_peak = 0.0;
int cursor_fft = 0;

int nparticles;
int nnodes;
double dx, dt;
int play = 1;
int clear = 0;
int plotting = 1;
int grid = 1;
int iter = 0;

int arg_freq = 0, arg_energy = 0, arg_particles = 0, arg_field = 0;

double maxEE = 0, maxTE = 0, maxKE = 0;
double minEE = 1e30, minTE = 1e30, minKE = 1e30;
double maxEE2 = 0, maxTE2 = 0, maxKE2 = 0;
double minEE2 = 1e30, minTE2 = 1e30, minKE2 = 1e30;
int cursor_x;
double TE0, EE0, KE0;
double TE=0, EE=0, KE=0;

#define ENERGY_HIST 1024

double EV[ENERGY_HIST];

double *gE, *grho, *gphi;
double maxphi = -1e10;
double maxrho = -1e10;
double maxE = -1e10;
double minphi = 1e10;
double minrho = 1e10;
double minE = 1e10;

static int
get_line(char *buf, size_t n, int prefix)
{
	while(getline(&buf, &n, stdin))
	{
		if(buf[0] == prefix)
			return 0;

		printf("%s", buf);
		fflush(stdout);
	}
	return -1;
}

static void
init_particles(void)
{
	int i;
	particle_t *p;
	char buf[MAX_LINE];

	if(get_line(buf, MAX_LINE, 'p'))
	{
		fprintf(stderr, "EOF before read config\n");
		exit(1);
	}

	sscanf(buf, "p %d %d %le %le", &nparticles, &nnodes, &dx, &dt);
	printf("Number of particles=%d, nnodes=%d\n", nparticles, nnodes);

	for(i = 0; i < MAX_HIST; i++)
		particles[i] = calloc(sizeof(particle_t), nparticles);


	grho = malloc(nnodes*sizeof(double));
	gphi = malloc(nnodes*sizeof(double));
	gE = malloc(nnodes*sizeof(double));

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
		grid = 1;
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
idle_particles(char *line)
{
	float t, x, u;
	int i, j, ret;
	particle_t *p;
	char buf[MAX_LINE];
	size_t n = MAX_LINE;


	if(!play)
	{
		//glutPostRedisplay();
		return;
	}

	hist = (hist + 1) % MAX_HIST;
	pos_i = (pos_i + 1) % MAX_POS;
	if(pos_i == 0)
	{
		//fprintf(stderr, "Completed FFT reached\n");
		recompute_fft = 1;
	}

	iter++;

	for(i = 0; i < nparticles; i++)
	{
		/* Copy old particle into particles1 */
		//memcpy(&particles1[prev_hist][i], &particles[hist][i], sizeof(particle_t));

		if(!fgets(buf, n, stdin))
		{
			fprintf(stderr, "EOF reached\n");
			exit(0);
		}

		ret = sscanf(buf, "%d %f %f",
				&j, &x, &u);

		if(ret != 3)
		{
			printf("ret=%d i=%d j=%d, exitting\n", ret, i, j);
			exit(1);
		}

		p = &particles[hist][j];
		p->x = x;
		p->u = u;

		if(i == 0)
		{
			pos_vec[pos_i] = u;
			j = (pos_i + MAX_POS-1) % MAX_POS;
			if(pos_vec[j] <= 0 && pos_vec[pos_i] > 0)
			{
				printf("loop iterations %d\n", iter);
				iter = 0;
			}
		}
	}

	if(arg_particles)
	{
		glutSetWindow(win1);
		glutPostRedisplay();
	}
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

		x1 = (p->x / (dx * nnodes)) * windW;
		y1 = (p->u / (MAX_V)) * windH;

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
	float h = 1.0, l = 0.1;

	n = get_curve(x, y, segments, pi);


	/* Erase a bit of the tail */
	glLineWidth(5.0);
	glColor3f(0,0,0);

	for(i = 0; i < n; i++)
	{
		glBegin(GL_LINE_STRIP);
		for(j = 0; j<segments[i]; j++)
		{
			glVertex2f(x[k], y[k]);
			k++;
			if(k >= 2) break;
		}
		glEnd();
		if(k >= 2) break;
	}

	k = 0;
	glLineWidth(2.0);

	for(i = 0; i < n; i++)
	{

		glBegin(GL_LINE_STRIP);
		for(j = 0; j<segments[i]; j++)
		{
			//cm = (float) k / MAX_HIST;
			//if(cm < 0.5) cm = 0.5;
			if(k < 2)
				cm = 0.3;
			else
				cm = 1.0;

			//fprintf(stderr, "k=%d, MAX_HIST=%d\n", k, MAX_HIST);

			ch = cm * h;
			cl = cm * l;
			if(pi % 2)
				glColor3f(ch, cl, cl);
			else
				glColor3f(cl, ch, cl);
			glVertex2f(x[k], y[k]);
			glVertex2f(x[k]-1, y[k]);
			k++;
		}
		//exit(0);
		glEnd();
	}
}

static void
plot_grid()
{
	int i;
	double x;
	double gray = 0.1;

	glColor3f(gray, gray, gray);

	for(i=0; i<nnodes; i++)
	{
		x = ((double)i) / nnodes * windW;

		glBegin(GL_LINES);
		glVertex2f(x, -1);
		glVertex2f(x, windH);
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
display_particles(void)
{
	if(!plotting)
		return;

	//fprintf(stderr, "Plotting particles\n");

	if(grid)
	{
		plot_grid();
		grid = 0;
	}
	plot_particles();
	//printf("iteration %d\n", iter);

	if (doubleBuffer) {
		glutSwapBuffers();
	} else {
		glFlush();
	}
}


void
idle_energy(char *line)
{
	int ret;

	if(!play)
	{
		//glutPostRedisplay();
		return;
	}

	TE0 = TE;
	EE0 = EE;
	KE0 = KE;

	ret = sscanf(line, "e %le %le %le",
			&TE, &EE, &KE);

	if(ret != 3)
	{
		printf("ret=%d, exitting\n", ret);
		exit(1);
	}

	glutSetWindow(win2);
	glutPostRedisplay();
}

static void
plot()
{
	char buf[MAX_LINE];
	int i;
	double x, yTE, yTE0, yEE, yEE0, yKE, yKE0;
	double lTE = (TE);
	double lEE = (EE);
	double lKE = (KE);
	double lTE0 = (TE0);
	double lEE0 = (EE0);
	double lKE0 = (KE0);
	int dirty = 0;

	//fprintf(stderr, "TE=%e\n", TE);
	//fprintf(stderr, "EE=%e\n", EE);

	x = (double) cursor_x;

	cursor_x = (cursor_x + 1) % windW;

	if(lTE > maxTE)
	{
		maxTE = 2.0 * lTE;
		dirty = 1;
	}

	//if(lEE > maxEE) maxEE = lEE;
	//if(lKE > maxKE) maxKE = lKE;

	//if(lTE < minTE) minTE = lTE;
	//if(lEE < minEE) minEE = lEE;
	//if(lKE < minKE) minKE = lKE;

	//if(maxTE == 0.0) maxTE = 1e-10;
	//if(maxEE == 0.0) maxEE = 1e-10;
	//if(maxKE == 0.0) maxKE = 1e-10;
	minTE = 0.0;
	minKE = 0.0;
	minEE = 0.0;

	maxEE = maxTE;
	maxKE = maxTE;

	yTE  = (lTE  - minTE) / (maxTE - minTE) * windH;
	yTE0 = (lTE0 - minTE) / (maxTE - minTE) * windH;
	yEE  = (lEE  - minEE) / (maxEE - minEE) * windH;
	yEE0 = (lEE0 - minEE) / (maxEE - minEE) * windH;
	yKE  = (lKE  - minKE) / (maxKE - minKE) * windH;
	yKE0 = (lKE0 - minKE) / (maxKE - minKE) * windH;


	//fprintf(stderr, "x=%e yTE0=%e\n", x, yTE0);
	//fprintf(stderr, "x=%e yTE=%e\n", x-1, yTE);

	glLineWidth(1.0);

	/* Clear previous segment */
	if(dirty)
		glColor3f(0.5, 0.5, 0.5);
	else
		glColor3f(0.0, 0.0, 0.0);
	glBegin(GL_LINES);
	glVertex2f(x, 0);
	glVertex2f(x, windH);
	glEnd();

	/* Draw total energy */
	glColor3f(1.0, 0.0, 0.0);
	glBegin(GL_LINES);
	glVertex2f(x-1, yTE0);
	glVertex2f(x, yTE);
	glEnd();

	/* Draw electrostatic energy */
	glColor3f(1.0, 1.0, 0.0);
	glBegin(GL_LINES);
	glVertex2f(x-1, yEE0);
	glVertex2f(x, yEE);
	glEnd();

	/* Draw kinetic energy */
	glColor3f(0.0, 1.0, 0.0);
	glBegin(GL_LINES);
	glVertex2f(x-1, yKE0);
	glVertex2f(x, yKE);
	glEnd();


}

void
display_energy(void)
{
	if(!play)
		return;
	//fprintf(stderr, "Plotting energy\n");
	plot();

	if (doubleBuffer) {
		glutSwapBuffers();
	} else {
		glFlush();
	}
}

static void
plot_fft()
{
	int x;
//	glClear(GL_COLOR_BUFFER_BIT);

	glLineWidth(1.0);

	glColor3f(0., 0., 0.);
	glBegin(GL_LINES);
	glVertex2f((double) cursor_fft, -1.0);
	glVertex2f((double) cursor_fft, windH+1);
	glEnd();

	glColor3f(1., 1., 0.);
	glBegin(GL_LINES);
	glVertex2f((double) cursor_fft - 1, (last_freq_peak/max_freq_peak) * windH);
	glVertex2f((double) cursor_fft, (freq_peak/max_freq_peak) * windH);
	glEnd();


	//fprintf(stderr, "Freq y = %f\n", (freq_peak/max_freq_peak) * windH);

//	glBegin(GL_LINES);
//	//glVertex2f(10.0, 20.0);
//	//glVertex2f(50.0, 100.0);
//	for(x = 0; x < windW; x++)
//	{
//		glVertex2f((double) x, -1.0);
//		glVertex2f((double) x, pos_fft[x]/maxfft * windH);
//	}
//	glEnd();

	cursor_fft = (cursor_fft + 1) % windW;
}

void
display_fft(void)
{
	if(!play)
		return;
	//fprintf(stderr, "Plotting fft\n");
	plot_fft();
	//fprintf(stderr, "Done plotting\n");


	if (doubleBuffer) {
		glutSwapBuffers();
	} else {
		glFlush();
	}
}

void
idle_fft()
{
	int i;
	//fprintf(stderr, "idle fft\n");

	if(!play)
		return;

	if(!recompute_fft)
		return;

	recompute_fft = 0;
	redraw_fft = 1;
	//fprintf(stderr, "idle fft runs\n");

	/* Compute the FFT and plot it */
	fftw_complex signal[MAX_POS];
	fftw_complex result[MAX_POS];

	fftw_plan plan = fftw_plan_dft_1d(MAX_POS, signal, result,
			FFTW_FORWARD, FFTW_ESTIMATE);

	for (i = 0; i < MAX_POS; ++i)
	{
		signal[i][0] = pos_vec[i];
		signal[i][1] = 0.0;
	}

	fftw_execute(plan);

	last_freq_peak = freq_peak;

	maxfft = 0.0;
	for (i = 1; i < MAX_POS/2; i++)
	{
		pos_fft[i] = sqrt(result[i][0] * result[i][0] +
			  result[i][1] * result[i][1]);

		pos_fft[i] = (pos_fft[i]);

		if(pos_fft[i] > maxfft)
		{
			maxfft = pos_fft[i];
			/* FIXME: The FFT is symmetric */
			freq_peak = ((double) i * 2.0) / (MAX_POS * DT);
			//fprintf(stderr, "maxfft = %e\n", maxfft);
			if (max_freq_peak < freq_peak) max_freq_peak = 3.0 * freq_peak;
		}
	}

	fprintf(stderr, "Frequency peak at %e Hz (max %e Hz)\n", freq_peak, max_freq_peak);
//	maxfft = 0.0;
//	for (i = 1; i < windW; i++)
//	{
//		pos_fft[i] = sqrt(result[i][0] * result[i][0] +
//			  result[i][1] * result[i][1]);
//
//		pos_fft[i] = (pos_fft[i]);
//
//		if(pos_fft[i] > maxfft)
//		{
//			maxfft = 1.5 * pos_fft[i];
//			//fprintf(stderr, "maxfft = %e\n", maxfft);
//		}
//	}
	//fprintf(stderr, "Finally set maxfft = %e\n", maxfft);

	fftw_destroy_plan(plan);

	/* Mark for redisplay */
	glutSetWindow(win3);
	glutPostRedisplay();
}

void
idle_field(char *line)
{
	int i, j, ret;
	particle_t *p;
	char buf[MAX_LINE];
	size_t n = MAX_LINE;
	double rho, phi, E;

	if(!play)
	{
		//glutPostRedisplay();
		return;
	}

//	maxphi = -1e10;
//	maxrho = -1e10;
//	maxE = -1e10;

	for(i = 0; i < nnodes; i++)
	{
		/* Copy old particle into particles1 */
		//memcpy(&particles1[prev_hist][i], &particles[hist][i], sizeof(particle_t));

		if(!fgets(buf, n, stdin))
		{
			fprintf(stderr, "EOF reached\n");
			exit(0);
		}

		ret = sscanf(buf, "%lf %lf %lf",
				&rho, &phi, &E);

		if(ret != 3)
		{
			printf("ret=%d i=%d, exitting\n", ret, i);
			exit(1);
		}

		if(maxphi < fabs(phi)) maxphi = fabs(phi);
		if(maxrho < fabs(rho)) maxrho = fabs(rho);
		if(maxE   < fabs(E))   maxE   = fabs(E);

		grho[i] = rho;
		gphi[i] = phi;
		gE[i] = E;
	}

	if(arg_field)
	{
		glutSetWindow(win4);
		glutPostRedisplay();
	}
}

static void
display_field()
{
	int i;
	double x, x0;
	double gray = 0.1;

	//glutSetWindow(win4);
	//glClearColor(0.0, 0.0, 0.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT);
	glColor3f(1.0, 1.0, 0.0);

	//glBegin(GL_LINES);
	glBegin(GL_LINE_STRIP);
	for(i=0; i<nnodes; i++)
	{
		x = ((double)i) / (nnodes) * windW;
		glVertex2f(x, grho[i]/maxrho * windH/2 + windH/2);
	}
	glEnd();
	glColor3f(1.0, 0.0, 0.0);
	glBegin(GL_LINE_STRIP);
	for(i=0; i<nnodes; i++)
	{
		x = ((double)i) / (nnodes) * windW;
		glVertex2f(x, (gphi[i]/maxphi) * windH/2 + windH/2);
//		glVertex2f(x+2, (gphi[i] - minphi)/(maxphi - minphi) * windH);
	}
	glEnd();
	glColor3f(0.0, 1.0, 0.0);
	glBegin(GL_LINE_STRIP);
	for(i=0; i<nnodes; i++)
	{
		x = ((double)i) / (nnodes) * windW;
		glVertex2f(x, gE[i]/maxE * windH/2 + windH/2);
//		glVertex2f(x+2, (gE[i] - minE)/(maxE - minE) * windH);
	}
	glEnd();
	glColor3f(.4, .4, .4);
	glBegin(GL_LINES);
	glVertex2f(0, windH/2);
	glVertex2f(windW, windH/2);
	glEnd();

}

void
idle()
{
	//fprintf(stderr, "idle\n");
	if(!play)
		return;

	char buf[MAX_LINE];
	size_t n = MAX_LINE;

	if(!fgets(buf, n, stdin))
		return;

	//fprintf(stderr, "Loading info\n");

	if(buf[0] == 'p')
		idle_particles(buf);
	else if(buf[0] == 'e' && arg_energy)
		idle_energy(buf);
	else if(buf[0] == 'f' && arg_field)
		idle_field(buf);

	if(arg_freq)
		idle_fft();
}

void
idle_redisplay()
{
	//fprintf(stderr, "idle redisplay\n");
	if(plotting)
	{
		//fprintf(stderr, "Calling redisplay for particles\n");
		glutPostRedisplay();
	}
}


void
visible_energy(int state)
{
	if (state == GLUT_VISIBLE) {
		glutIdleFunc(idle);
	} else {
		glutIdleFunc(NULL);
	}
}

void
visible_particles(int state)
{
	if (state == GLUT_VISIBLE) {
		glutIdleFunc(idle);
	} else {
		glutIdleFunc(NULL);
	}
}

int
parse_args(int argc, char *argv[])
{
	int opt;
	int any = 0;

	while((opt = getopt(argc, argv, "epfF")) != -1)
	{
		switch(opt)
		{
			case 'e':
				arg_energy = 1;
				any = 1;
				break;
			case 'p':
				arg_particles = 1;
				any = 1;
				break;
			case 'f':
				arg_field = 1;
				any = 1;
				break;
			case 'F':
				arg_freq = 1;
				any = 1;
				break;
			default:
				fprintf(stderr, "Usage: %s [-e] [-p] [-f] [-F]\n", argv[0]);
				exit(EXIT_FAILURE);
		}
	}

	if(!any)
	{
		fprintf(stderr, "Please choose at least one plot: -e, -p or -f\n");
		exit(EXIT_FAILURE);
	}

	return 0;
}

int
main(int argc, char **argv)
{
	int show_energy = 1;
	int winy = 20;
	GLenum type;

	parse_args(argc, argv);

	glutInitWindowSize(windW, windH);
	glutInit(&argc, argv);

	type = GLUT_RGB;
	type |= (doubleBuffer) ? GLUT_DOUBLE : GLUT_SINGLE;
	glutInitDisplayMode(type);

	init_particles();

	if(arg_particles)
	{
		win1 = glutCreateWindow("plot");
		glutPositionWindow(5, winy);
		winy += windH + 30;
		glutReshapeFunc(Reshape);
		glutKeyboardFunc(Key);
		glutDisplayFunc(display_particles);
		glClearColor(0.0, 0.0, 0.0, 0.0);
		glDisable(GL_DITHER);
	}

	if(arg_field)
	{
		win4 = glutCreateWindow("plot fields");
		glutPositionWindow(5, winy);
		winy += windH + 30;
		glutReshapeFunc(Reshape);
		glutKeyboardFunc(Key);
		//glutVisibilityFunc(visible_energy);
		glutDisplayFunc(display_field);
	}

	if (arg_energy)
	{
		win2 = glutCreateWindow("plot energy");
		glutPositionWindow(5, winy);
		winy += windH + 30;
		glutReshapeFunc(Reshape);
		glutKeyboardFunc(Key);
		glutDisplayFunc(display_energy);
	}


	if(arg_freq)
	{
		win3 = glutCreateWindow("plot frequency");
		glutPositionWindow(5, winy);
		winy += windH + 30;
		glutReshapeFunc(Reshape);
		glutKeyboardFunc(Key);
		//glutVisibilityFunc(visible_energy);
		glutDisplayFunc(display_fft);
	}


	glutIdleFunc(idle);

	glutMainLoop();
	return 0;             /* ANSI C requires main to return int. */
}
