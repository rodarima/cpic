#include "video.h"

#include <stdio.h>
#include <string.h>
#include <GL/gl.h>

video_t *
video_init(int w, int h, const char *fname)
{
	char *cmd;
	video_t *v;

	v = safe_malloc(sizeof(*v));

	v->w = w;
	v->h = h;

	// start ffmpeg telling it to expect raw rgba 720p-60hz frames
	// -i - tells it to read frames from stdin
	asprintf(&cmd, "/usr/bin/ffmpeg -r 60 -f rawvideo "
			"-pix_fmt rgba -s %dx%d -i - "
			"-y -crf 20 -vf vflip,format=pix_fmts=yuv420p %s",
			w, h, fname);

	v->pipe = popen(cmd, "w");

	if(!v->pipe)
	{
		perror("popen");
		return NULL;
	}

	free(cmd);

	v->buffer = safe_malloc(w * h * sizeof(int));

	return v;
}

int
video_update(video_t *v)
{
	/* Send to ffmpeg */
	glReadPixels(0, 0, v->w, v->h, GL_RGBA, GL_UNSIGNED_BYTE, v->buffer);

	fwrite(v->buffer, v->w * v->h * sizeof(int), 1, v->pipe);

	return 0;
}

int
video_finish(video_t *v)
{
	pclose(v->pipe);
	free(v->buffer);
	free(v);

	return 0;
}
