#include "plot.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <GLFW/glfw3.h>
#include <fftw3.h>
#include <unistd.h>

#include "mat.h"
#include "sim.h"
#include "specie.h"
#include "log.h"
#include "config.h"

GLenum doubleBuffer = GL_FALSE;
GLint windW = 800, windH = 300;
int win1, win2, win3, win4;

//#define MAX_HIST 10
#define MAX_LINE 256
#define MAX_POS 819
#define MAX_LOOP 5

particle_t *particles[MAX_HIST];
int hist = 0;

double pos_vec[MAX_POS];
double pos_fft[MAX_POS];
int pos_i = 0;
int maxloops = 0;
int redraw_fft = 0, recompute_fft = 0;
double maxfft = 0.0;
double last_freq_peak = 0.0;
double freq_peak = 0.0;
double max_freq_peak = 0.0;
int cursor_fft = 0;

int dim;
int nparticles;
int nblocks, blocksize, nnodes;
double L, dx, dt;
int play = 1;
int clear = 0;
int plotting = 1;
int fast = 0;
int grid = 1;
int iter = 0;
int loop_n = 0;
double maxfps = 0.0;
double maxv = 0.0;
double frame_dt = 0.0;

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
double trigger_factor = 0.0;

config_t conf;

int
parse_config(plot_t *plot, config_t *conf)
{
	int size[2];
	const char *video_fname = NULL;

	/* First set all direct configuration variables */
	config_lookup_array_int(conf, "plot.size", size, 2);
	plot->w = size[0];
	plot->h = size[1];

	config_lookup_float(conf, "plot.max_fps", &plot->maxfps);
	config_lookup_float(conf, "plot.max_velocity", &plot->maxv);
	config_lookup_int(conf, "plot.max_loops", &plot->maxloops);
	config_lookup_float(conf, "plot.trigger_factor", &plot->trigger_factor);

	config_lookup_string(conf, "plot.video_file", &video_fname);

	if(video_fname)
	{
		plot->video_enabled = 1;
		plot->video = video_init(plot->w, plot->h, video_fname);
	}

	/* Then compute the rest */
	nnodes = nblocks * blocksize;
	dx = L / nnodes;

	if(maxfps > 0)
		frame_dt = 1/maxfps;

	fprintf(stderr, "maxfps set to %f, frame_dt = %e\n", maxfps, frame_dt);

	return 0;
}

int
read_config(config_t *conf)
{
	char buf[MAX_LINE];
	size_t n = MAX_LINE;
	FILE *f;

	if(!fgets(buf, n, stdin))
	{
		fprintf(stderr, "EOF before reading config file line\n");
		return 1;
	}

	/* Remove newline if any */
	buf[strcspn(buf, "\n")] = 0;

	f = fopen(buf, "r");

	if(!f)
	{
		perror("fopen");
		return 1;
	}

	config_init(conf);

	/* Read the configuration from stdin */
	if(config_read(conf, f) == CONFIG_FALSE)
	{
		err("Configuration read failed\n");
		return 1;
	}

	fclose(f);

	return 0;
}

int
plot_redraw(plot_t *plot)
{
	HMGL gr;
	int i, j;
	specie_t *s;
	particle_t *p;
	sim_t *sim;
	double tot_p;

	sim = plot->sim;

	gr = plot->gr;

	mgl_clf(gr);

	mgl_set_font_size(gr, 3.0);
	mgl_set_mark_size(gr, 0.35);

	/* Update energies */
	mgl_data_roll(plot->EE, 'x', 1);
	mgl_data_roll(plot->KE, 'x', 1);
	mgl_data_roll(plot->TE, 'x', 1);
	mgl_data_roll(plot->pE, 'x', 1);
	mgl_data_roll(plot->P[X], 'x', 1);
	mgl_data_roll(plot->P[Y], 'x', 1);

	mgl_data_put_val(plot->EE, sim->energy_electrostatic, MAX_HIST-1, 0, 0);
	mgl_data_put_val(plot->KE, sim->energy_kinetic, MAX_HIST-1, 0, 0);
	mgl_data_put_val(plot->TE,
			sim->energy_kinetic+sim->energy_electrostatic,
			MAX_HIST-1, 0, 0);

	tot_p = sqrt(sim->total_momentum[X]*sim->total_momentum[X]+
			sim->total_momentum[Y]*sim->total_momentum[Y]);

	mgl_data_put_val(plot->pE, 1000 + 100*sim->species[0].particles[0].E[X], MAX_HIST-1, 0, 0);
	//mgl_data_put_val(plot->P[X], sim->total_momentum[X], MAX_HIST-1, 0, 0);
	mgl_data_put_val(plot->P[X], tot_p, MAX_HIST-1, 0, 0);
	mgl_data_put_val(plot->P[Y], sim->total_momentum[Y], MAX_HIST-1, 0, 0);


	mgl_subplot(gr, 2, 2, 0, "");
	mgl_set_ranges(gr, 0, 1, 0,
			mgl_data_max(plot->TE), 0, 1);
	mgl_title(gr, "Energy", "", 5.0);
	mgl_plot(gr, plot->EE, "r", "");
	mgl_plot(gr, plot->KE, "b", "");
	mgl_plot(gr, plot->TE, "k", "");
	mgl_plot(gr, plot->P[X], "g", "");
	mgl_plot(gr, plot->P[Y], "G", "");
//	mgl_plot(gr, plot->pE, "g", "");
	mgl_add_legend(gr, "Total", "k");
	mgl_add_legend(gr, "Kinetic", "b");
	mgl_add_legend(gr, "Electrostatic", "r");
	mgl_legend(gr, 0, "#", "");
//	mgl_axis(gr, "y", "", "");


//	mgl_subplot(gr, 2, 2, 0, "");
//	mgl_title(gr, "Particle x-y space", "", 5.0);
//	mgl_plot_xy(gr, plot->x, plot->y, "#s ", "");
//	//mgl_axis_grid(gr, "xy", "", "");
//	mgl_axis(gr, "xy", "", "");
//	mgl_set_ranges(gr, 0.0, 64.0, 0.0, 64.0, -1, 1);
//	mgl_label(gr, 'x', "Position x", 0.0, "");
//	mgl_label(gr, 'y', "Position y", 0.0, "");


	mgl_subplot(gr, 2, 2, 1, "");
	mgl_set_ranges(gr, 0.0, sim->nnodes[X], 0.0, sim->nnodes[Y],
			mgl_data_min(plot->rho), mgl_data_max(plot->rho));
//	mgl_set_ranges(gr, 0.0, sim->nnodes[X]-1, 0.0, sim->nnodes[Y]-1,
//			-100, 100);
	mgl_set_ticks(gr, 'x', 1, 1, 0);
	mgl_set_ticks(gr, 'y', 1, 1, 0);
	mgl_title(gr, "Charge density \\rho", "", 5.0);
//	mgl_colorbar(gr, "<");
	mgl_contf(gr, (HCDT) plot->rho, "", "");
//	mgl_surf(gr, (HCDT) plot->rho, "b", "");
	mgl_axis_grid(gr, "xy", "", "");
	mgl_axis(gr, "xy", "", "");
//	mgl_boxs(gr, (HCDT) plot->rho, "wk", "");

	mgl_subplot(gr, 2, 2, 2, "");
	mgl_set_ranges(gr, 0.0, sim->nnodes[X], 0.0, sim->nnodes[Y],
			mgl_data_min(plot->phi), mgl_data_max(plot->phi));
	mgl_title(gr, "Electric potential \\phi", "", 5.0);
	mgl_contf(gr, (HCDT) plot->phi, "", "");
//	mgl_boxs(gr, (HCDT) plot->phi, "wk", "");
//	mgl_colorbar(gr, ">");
	mgl_axis_grid(gr, "xy", "", "");
	mgl_axis(gr, "xy", "", "");

//	mgl_subplot(gr, 2, 2, 3, "");
//	mgl_set_ranges(gr, 0.0, 64.0, -10, 10,
//			mgl_data_min(plot->E0), mgl_data_max(plot->E0));
//	mgl_title(gr, "Electric field E_x", "", 5.0);
//	mgl_contf(gr, (HCDT) plot->E0, "", "");
//	mgl_axis_grid(gr, "xy", "", "");
//	mgl_axis(gr, "xy", "", "");
//	mgl_colorbar(gr, ">");

//	mgl_subplot(gr, 2, 2, 3, "");
//	mgl_title(gr, "Electric field E", "", 5.0);
//	mgl_axis(gr, "xy", "", "");
//	mgl_set_meshnum(gr, 30);
//	mgl_set_ranges(gr, 0.0, sim->nnodes[X], 0.0, sim->nnodes[Y],
//			mgl_data_min(plot->E[X]),
//			mgl_data_max(plot->E[Y]));
//	mgl_vect_2d(gr, plot->E[X], plot->E[Y], "b2", "");

	mgl_subplot(gr, 2, 2, 3, "");
	mgl_title(gr, "Particle x-y space", "", 5.0);
	//mgl_axis_grid(gr, "xy", "", "");
	mgl_set_ranges(gr, 0.0, sim->L[X], 0.0, sim->L[Y], -1, 1);
	mgl_axis(gr, "xy", "", "");

	for(j=0; j<sim->nspecies; j++)
	{
		s = &plot->sim->species[j];
		for(i = 0; i<s->nparticles; i++)
		{
			p = &s->particles[i];

			mgl_data_set_value(plot->x, p->x[0], i, 0, 0);
			mgl_data_set_value(plot->v, p->u[0], i, 0, 0);
			mgl_data_set_value(plot->y, p->x[1], i, 0, 0);

			//mgl_mark(gr, p->x[0], p->u[0], 0.0, "o");
		}
		mgl_plot_xy(gr, plot->x, plot->y, "#s ", "");
	}
	mgl_label(gr, 'x', "Position x", 0.0, "");
	mgl_label(gr, 'y', "Position y", 0.0, "");



	mgl_finish(gr);
	//glFlush();
	return 0;
}

plot_t *
plot_init(sim_t *sim)
{
	plot_t *plot;

	plot = calloc(sizeof(plot_t), 1);

	plot->sim = sim;
	plot->wait = 1;

	if(parse_config(plot, sim->conf))
		return NULL;

	return plot;
}

static void
key(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	plot_t *plot;

	plot = glfwGetWindowUserPointer(window);

	if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
		glfwSetWindowShouldClose(window, GLFW_TRUE);

	if (key == GLFW_KEY_SPACE && action == GLFW_PRESS)
		plot->paused = !plot->paused;

	if (key == GLFW_KEY_W && action == GLFW_PRESS)
		plot->wait = !plot->wait;
}



/* Executed in plotter thread, should never return */
void *
plot_loop(void *p)
{
	plot_t *plot;
	sim_t *sim;
	GLFWwindow* window;

	sim = (sim_t *) p;

	plot = plot_init(sim);

	if(!glfwInit())
		return NULL;

	glfwWindowHint(GLFW_FLOATING, GLFW_TRUE);

	window = glfwCreateWindow(plot->w, plot->h,
			"plot particles", NULL, NULL);

	glfwSetWindowTitle(window, "plot particles");


	glfwSetWindowUserPointer(window, plot);

	if(!window)
		return NULL;

	glfwSetKeyCallback(window, key);
	glfwMakeContextCurrent(window);
	glfwSwapInterval(1);

	int mx = sim->nnodes[X];
	int my = sim->nnodes[Y];
	int mz = sim->nnodes[Z];
	int np = sim->species[0].nparticles;

	mgl_set_num_thr(1);
	plot->gr = mgl_create_graph_gl();

	plot->rho = mgl_create_data();
	plot->phi = mgl_create_data();
	plot->E[X] = mgl_create_data();
	plot->E[Y] = mgl_create_data();
	plot->x = mgl_create_data_size(np, 1, 1);
	plot->y = mgl_create_data_size(np, 1, 1);
	plot->v = mgl_create_data_size(np, 1, 1);

	plot->EE = mgl_create_data_size(MAX_HIST, 1, 1);
	plot->KE = mgl_create_data_size(MAX_HIST, 1, 1);
	plot->TE = mgl_create_data_size(MAX_HIST, 1, 1);
	plot->pE = mgl_create_data_size(MAX_HIST, 1, 1);
	plot->P[X] = mgl_create_data_size(MAX_HIST, 1, 1);
	plot->P[Y] = mgl_create_data_size(MAX_HIST, 1, 1);

	mgl_data_put_val(plot->EE, 0.0, -1, -1, -1);
	mgl_data_put_val(plot->KE, 0.0, -1, -1, -1);
	mgl_data_put_val(plot->TE, 0.0, -1, -1, -1);
	mgl_data_put_val(plot->pE, 0.0, -1, -1, -1);
	mgl_data_put_val(plot->P[X], 0.0, -1, -1, -1);
	mgl_data_put_val(plot->P[Y], 0.0, -1, -1, -1);

	pthread_mutex_lock(&sim->lock);
	mgl_data_link(plot->phi, sim->field->phi->data, mx, my, mz);
	mgl_data_link(plot->rho, sim->field->rho->data, mx, my, mz);
	mgl_data_link(plot->E[X], sim->field->E[X]->data, mx, my, mz);
	mgl_data_link(plot->E[Y], sim->field->E[Y]->data, mx, my, mz);
	pthread_mutex_unlock(&sim->lock);

	mgl_set_font_size(plot->gr, 2.0);
	//mgl_set_quality(plot->gr, 4);

	while (!glfwWindowShouldClose(window))
	{
		glfwMakeContextCurrent(window);
		glfwGetFramebufferSize(window, &plot->w, &plot->h);

		if(plot->wait)
		{
			glViewport(0, 0, plot->w, plot->h);
			//glClear(GL_COLOR_BUFFER_BIT);

			mgl_set_size(plot->gr, plot->w, plot->h);
		}

		pthread_mutex_lock(&sim->lock);
		while(sim->run == 1)
			pthread_cond_wait(&sim->signal, &sim->lock);

		if(plot->wait)
		{

			plot_redraw(plot);

	//		glfwMakeContextCurrent(window);

			//plot_energy(plot);

		}

		while(plot->paused)
		{
			glfwPollEvents();
		}

		sim->run = 1;
		pthread_cond_signal(&sim->signal);
		pthread_mutex_unlock(&sim->lock);

		if(plot->video_enabled)
		{
			video_update(plot->video);
		}

		if(plot->wait)
		{
			glfwSwapBuffers(window);
		}

		glfwPollEvents();
	}

	if(plot->video_enabled)
	{
		video_finish(plot->video);
	}

	glfwTerminate();
	return NULL;
}

/* Executed in the simulator thread. Just create the plotter thread and exit */
int
plot_thread_init(sim_t *sim)
{
	if(pthread_create(&sim->plot_thread, NULL,
				plot_loop, (void *) sim) != 0)
	{
		perror("pthread_create");
		return 1;
	}

	return 0;
}
