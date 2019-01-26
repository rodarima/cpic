
/* Copyright (c) Mark J. Kilgard, 1994. */

/**
 * (c) Copyright 1993, Silicon Graphics, Inc.
 * ALL RIGHTS RESERVED 
 * Permission to use, copy, modify, and distribute this software for 
 * any purpose and without fee is hereby granted, provided that the above
 * copyright notice appear in all copies and that both the copyright notice
 * and this permission notice appear in supporting documentation, and that 
 * the name of Silicon Graphics, Inc. not be used in advertising
 * or publicity pertaining to distribution of the software without specific,
 * written prior permission. 
 *
 * THE MATERIAL EMBODIED ON THIS SOFTWARE IS PROVIDED TO YOU "AS-IS"
 * AND WITHOUT WARRANTY OF ANY KIND, EXPRESS, IMPLIED OR OTHERWISE,
 * INCLUDING WITHOUT LIMITATION, ANY WARRANTY OF MERCHANTABILITY OR
 * FITNESS FOR A PARTICULAR PURPOSE.  IN NO EVENT SHALL SILICON
 * GRAPHICS, INC.  BE LIABLE TO YOU OR ANYONE ELSE FOR ANY DIRECT,
 * SPECIAL, INCIDENTAL, INDIRECT OR CONSEQUENTIAL DAMAGES OF ANY
 * KIND, OR ANY DAMAGES WHATSOEVER, INCLUDING WITHOUT LIMITATION,
 * LOSS OF PROFIT, LOSS OF USE, SAVINGS OR REVENUE, OR THE CLAIMS OF
 * THIRD PARTIES, WHETHER OR NOT SILICON GRAPHICS, INC.  HAS BEEN
 * ADVISED OF THE POSSIBILITY OF SUCH LOSS, HOWEVER CAUSED AND ON
 * ANY THEORY OF LIABILITY, ARISING OUT OF OR IN CONNECTION WITH THE
 * POSSESSION, USE OR PERFORMANCE OF THIS SOFTWARE.
 * 
 * US Government Users Restricted Rights 
 * Use, duplication, or disclosure by the Government is subject to
 * restrictions set forth in FAR 52.227.19(c)(2) or subparagraph
 * (c)(1)(ii) of the Rights in Technical Data and Computer Software
 * clause at DFARS 252.227-7013 and/or in similar or successor
 * clauses in the FAR or the DOD or NASA FAR Supplement.
 * Unpublished-- rights reserved under the copyright laws of the
 * United States.  Contractor/manufacturer is Silicon Graphics,
 * Inc., 2011 N.  Shoreline Blvd., Mountain View, CA 94039-7311.
 *
 * OpenGL(TM) is a trademark of Silicon Graphics, Inc.
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <GL/glut.h>

#include "specie.h"

GLenum doubleBuffer;
GLint windW = 500, windH = 500;

particle_t *particles0;
particle_t *particles1;
int nparticles;
int shape;
float dx, dt;
int play = 1;

static void
Init(void)
{
	int i;
	particle_t *p;

	scanf("%d %d %f %f", &nparticles, &shape, &dx, &dt);
	printf("Number of particles=%d, shape=%d\n", nparticles, shape);

	particles0 = calloc(sizeof(particle_t), nparticles);
	particles1 = calloc(sizeof(particle_t), nparticles);

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
	case 'q':
	case 27:
		exit(0);
	}
}

void
Idle(void)
{
	float t;
	int i, j, ret;
	particle_t *p;

	if(!play)
	{
		glutPostRedisplay();
		return;
	}

	for(i = 0; i < nparticles; i++)
	{
		p = &particles0[i];
		ret = scanf("%d %f %f",
				&j, &p->x, &p->u);
		if(ret == EOF)
		{
			printf("EOF reached\n");
			exit(0);
		}

		if((ret != 3) || (j != i))
		{
			printf("ret=%d i=%d j=%d, exitting\n", ret, i, j);
			exit(1);
		}
	}

	glutPostRedisplay();
}

static void
plot_particle(int i)
{
	particle_t *p0 = &particles0[i];
	particle_t *p1 = &particles1[i];

	float x0, y0;
	float x1, y1, width;

	x0 = (p0->x / (dx * shape)) * windW;
	y0 = (p0->u / (4.0 * 3e8)) * windH;
	//y0 = ((float) i / (float) nparticles) * windH;
	//x0 += windW / 2.0;
	y0 += windH / 2.0;

	x1 = (p1->x / (dx * shape)) * windW;
	y1 = (p1->u / (4.0 * 3e8)) * windH;
	//y0 = ((float) i / (float) nparticles) * windH;
	//x0 += windW / 2.0;
	y1 += windH / 2.0;

	if(fabs(x0 - x1) > dx)
		x1 = x0;


//	if (!(x >= 0.0 && x < windW && y >= 0.0 && y < windH))
//	{
//		printf("Particle out of bounds x=%f(%d) y=%f(%d)\n",
//				x, windW, y, windH);
//	}

	glLineWidth(2.0);
	if(i % 2)
		glColor3f(1.0, 0.3, 0.3);
	else
		glColor3f(0.3, 1.0, 0.3);

	glBegin(GL_LINES);
	glVertex2f(x0, y0);
	glVertex2f(x1, y1);
//	glVertex2f(x+1, y);
//	glVertex2f(x+1, y+1);
	glEnd();
}


static void
plot_particles()
{
	int i;
	glClear(GL_COLOR_BUFFER_BIT);

	for(i=0; i<nparticles; i++)
	{
		plot_particle(i);
		memcpy(&particles1[i], &particles0[i], sizeof(particle_t));
	}
}

void
Display(void)
{
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
	glutCreateWindow("two streams");
	//glutFullScreen();

	Init();

	glutReshapeFunc(Reshape);
	glutKeyboardFunc(Key);
	glutVisibilityFunc(Visible);
	glutDisplayFunc(Display);
	glutMainLoop();
	return 0;             /* ANSI C requires main to return int. */
}
