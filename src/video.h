#pragma once

#define _GNU_SOURCE
#include <stdio.h>

struct video {
	int w, h;
	FILE *pipe;
	unsigned char *buffer;
};

typedef struct video video_t;

video_t *
video_init(int w, int h, const char *fname);

int
video_update(video_t *v);

int
video_finish(video_t *v);
