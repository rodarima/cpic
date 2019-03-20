#include <stdbool.h>
#include <mgl2/mgl_cf.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <GLFW/glfw3.h>

#include "mgl-glue.h"

#define N 100
#define FPS 30
#define FRAME_PERIOD ((int) floor(1e9/FPS))

HMGL gr;
HMDT dat, du, dv;

double zz[N*N], uu[N*N], vv[N*N], shift = 0.0;

int frames = 0, updates = 0;
struct timespec later = {0, 0};

void
plot()
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
update_data()
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
}

static void
key(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
		glfwSetWindowShouldClose(window, GLFW_TRUE);
}

int main(int argc, char *argv[])
{
	int width = 1000, height = 500;

	if(!glfwInit())
		return 1;

	GLFWwindow* window = glfwCreateWindow(width, height,
			"plot particles", NULL, NULL);

	if(!window)
		return 1;

	glfwSetKeyCallback(window, key);
	glfwMakeContextCurrent(window);
	glfwSwapInterval(1);



	mgl_set_num_thr(1);

	dat = mgl_create_data_size(N, N, 1);
	du = mgl_create_data_size(N, N, 1);
	dv = mgl_create_data_size(N, N, 1);

	gr = mgl_create_graph_gl();

	mgl_set_font_size(gr, 2.0);




	while (!glfwWindowShouldClose(window))
	{

		update_data();

		glfwGetFramebufferSize(window, &width, &height);

		glViewport(0, 0, width, height);
		glClear(GL_COLOR_BUFFER_BIT);

		mgl_set_size(gr, width, height);

		plot();

		glfwSwapBuffers(window);
		glfwPollEvents();
	}

	glfwTerminate();

	return 0;
}
