// Simple demo of a 2D line plot.
//
// Copyright (C) 2011  Alan W. Irwin
//
// This file is part of PLplot.
//
// PLplot is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published
// by the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// PLplot is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with PLplot; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
//
//

//#include "plcdemos.h"
#include <plplot/plplot.h>
#include <math.h>
#include <unistd.h>

#define NSIZE    101
#define DT 1e-2

int
data(int n, PLFLT *x, PLFLT *y, PLFLT t)
{
	int i;

	for(i=0; i<n; i++)
	{
		y[i] = cos(i*DT + t);
	}

	return 0;
}

int
main( int argc, char *argv[] )
{
	PLFLT t = 0.0;
	PLFLT x[NSIZE], y[NSIZE];
	PLFLT xmin = 0., xmax = 1., ymin = -1., ymax = 1.;
	int i;

	// Prepare data to be plotted.
	for ( i = 0; i < NSIZE; i++ )
	{
		x[i] = (PLFLT) i/NSIZE;
		y[i] = ymax * x[i] * x[i];
	}

	// Parse and process command line arguments
	plparseopts( &argc, argv, PL_PARSE_FULL );

	// Initialize plplot
	plinit();

	// Create a labelled box to hold the plot.
	plenv( xmin, xmax, ymin, ymax, 0, 0 );

	while(1)
	{
		//pllab( "x", "y=100 x#u2#d", "Simple PLplot demo of a 2D line plot" );
		data(NSIZE, x, y, t);
		t += DT;
		// Plot the data that was prepared above.
		plline( NSIZE, x, y );
		usleep(100000);
		plclear();
	}


	// Close PLplot library
	plend();

	exit( 0 );
}
