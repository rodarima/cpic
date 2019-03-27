#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <GLFW/glfw3.h>

#define N 100

#define WIDTH 500
#define HEIGHT 500

int
main(int argc, char *argv[])
{
	int i;
	char* cmd;
	FILE* f;
	int *buf;

	// start ffmpeg telling it to expect raw rgba 720p-60hz frames
	// -i - tells it to read frames from stdin
	asprintf(&cmd, "/usr/bin/ffmpeg -r 60 -f rawvideo -pix_fmt rgba -s %dx%d -i - "
		"-threads 0 -preset fast -y -pix_fmt yuv420p -crf 21 -vf vflip output.mp4",
		WIDTH, HEIGHT);


	fprintf(stderr, "Trying to run: %s\n", cmd);
	f = popen(cmd, "w");

	if(!f)
	{
		perror("popen");
		return 1;
	}

	buf = malloc(WIDTH * HEIGHT * sizeof(int));

	if(!buf)
	{
		perror("malloc");
		return 1;
	}

	if(!glfwInit())
		return 1;

	GLFWwindow* window = glfwCreateWindow(WIDTH, HEIGHT,
			"recording opengl test", NULL, NULL);

	glfwHideWindow(window);
	glfwMakeContextCurrent(window);

	glViewport(0, 0, WIDTH, HEIGHT);

	glMatrixMode( GL_PROJECTION );
	glLoadIdentity();

	glOrtho(0.0, WIDTH, 0.0, HEIGHT, 0.0, 1.0);

	for(i=0; i<N; i++)
	{

		glClearColor(0, 0, 0, 0);
		glClear(GL_COLOR_BUFFER_BIT);
		glPointSize(10);
		glLineWidth(2.5);
		glColor3f(1.0, 0.0, 0.0);

		glBegin(GL_LINES);
		glVertex2f(10.0, 10.0);
		glVertex2f(50.0, 50.0 + 5.0 * i);
		glEnd();

		//glFlush();
		glFinish();
		//glutSwapBuffers();
		//glfwSwapBuffers(window);
		//glfwPollEvents();

		/* Send to ffmpeg */
		glReadPixels(0, 0, WIDTH, HEIGHT, GL_RGBA, GL_UNSIGNED_BYTE, buf);

		//memset(buf, 0xffffffff, 20*HEIGHT*sizeof(int));

		fwrite(buf, WIDTH*HEIGHT*sizeof(int), 1, f);
	}

	glfwTerminate();

	pclose(f);

	return 0;
}
