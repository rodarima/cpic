#include <mgl2/opengl.h>
#include <mgl2/glut.h>
#include <mgl2/canvas.h>
#include "mgl-glue.h"

//-----------------------------------------------------------------------------
/// Base class for windows containing MathGL graphics
class CanvasGLUT : public mglGraph
{
	public:

	~CanvasGLUT() {}
	CanvasGLUT() : mglGraph(1, 1600, 800) {}
	HMGL get_gr() { return gr; }

};

HMGL
canvas_create()
{
	CanvasGLUT *c = new CanvasGLUT();
	return (HMGL) c->get_gr();
}


void
canvas_set_size(HMGL gr, int w, int h)
{
	CanvasGLUT *g = dynamic_cast<CanvasGLUT *>(gr);
	//mgl_set_size(gr, w, h);
	//g->SetSize(w, h, true);
	return;
//	CanvasGLUT *g = dynamic_cast<CanvasGLUT *>(gr);
//	//g->set_size(w, h);
//	//g->SetSize(w, h, true);
//	g->mglCanvas::SetSize(1,1,false);
//	g->SetSize(1,1,true);
//	fprintf(stderr, "Size set to %dx%d\n", w, h);
}


//int
//main(int argc, char *argv[])
//{
//	CanvasGLUT *gr = new CanvasGLUT();
//
//	gr->set_size(100, 50);
//	return 0;
//}
