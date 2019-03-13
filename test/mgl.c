#include <stdbool.h>
#include <mgl2/mgl_cf.h>
#include <mgl2/glut.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>

#define N 100

double yy[N];

int
draw(HMGL gr, void *p)
{
	printf("draw called\n");
	HMDT dat = (HMDT) p;
	int i;

	for(i=0; i<N; i++)
	{
		yy[i] = ((double) rand()) / RAND_MAX;
	}
	mgl_data_set_double(dat, yy, N, 1, 1);
	mgl_plot(gr, (HCDT) dat, "", "");
	return 1;
}

void
reload(void *p)
{
	printf("reload called\n");
}

int main()
{
	HMDT dat;
	HMGL gr;
	int i;

	dat = mgl_create_data_size(N, 1, 1);

	for(i=0; i<N; i++)
	{
		yy[i] = ((double) rand()) / RAND_MAX;
	}

	mgl_data_set_double(dat, yy, N, 1, 1);
	gr = mgl_create_graph_glut(draw, "title", (void *) dat, reload);

	return 0;
}

//#include <mgl2/Fl_MathGL.h>
//int main(int argc,char **argv)
//{
//    mglFLTK gr("test");
//    mgl_fltk_thr();
//
//    mglPoint pnt;   // yours data
//    for(int i=0;i<10;i++)   // do calculation
//    {
//#if defined(WIN32) || defined(_MSC_VER) || defined(__BORLANDC__)
//        Sleep(1000);
//#else
//        sleep(1);           // which can be very long
//#endif
//        pnt = mglPoint(2*mgl_rnd()-1,2*mgl_rnd()-1);
//        gr.Clf();           // make new drawing
//        gr.Line(mglPoint(),pnt,"Ar2");
//        char str[10] = "i=0";   str[3] = '0'+i;
//        gr.Puts(mglPoint(),"");
//        gr.Update();        // update window
//    }
//    return 0;   // finish calculations and close the window
//}
