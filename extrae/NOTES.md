Thu, 13 Jun 2019 10:02:29 +0200

Problems observed when trying to run extrae with:

- nanos6 + TAMPI
- FFTW3 with threads

## Using TAMPI as dynamic library: `-ltampi-c`

Using `LD_PRELOAD=${EXTRAE_HOME}/lib/libnanostrace.so` I can run extrae, but I
only obtain one CPU in paraver. I have `NANOS6_EXTRAE_AS_THREADS=0` but the same
happens with `NANOS6_EXTRAE_AS_THREADS=1`.

## Using TAMPI as static library: `$(TAMPI_HOME)/lib/libtampi-c.a`

	LD_PRELOAD		AS_THREADS	TAMPI		RESULT

	libnanostrace.so	1		static		Crash 1
	libnanostrace.so	0		static		Crash 1
	libnanosmpitrace.so	1		static		Crash 1
	libnanosmpitrace.so	0		static		Crash 1
	libmpitrace.so		1		static		Crash 1
	libmpitrace.so		0		static		Crash 1
	libnanostrace.so	1		dynamic		Only 1 process
	libnanostrace.so	0		dynamic		Only 1 process
	libnanosmpitrace.so	1		dynamic		Crash 1
	libnanosmpitrace.so	0		dynamic		Crash 1
	libmpitrace.so		1		dynamic		Crash 1
	libmpitrace.so		0		dynamic		Crash 1

The column `AS_THREADS` means the environment variable
`NANOS6_EXTRAE_AS_THREADS`.

Crash 1:

	*** The MPI_Init_thread() function was called before MPI_INIT was invoked.
	*** This is disallowed by the MPI standard.
	*** Your MPI job will now abort.
	[s06r2b52:425577] Local abort before MPI_INIT completed completed
	successfully, but am not able to aggregate error messages, and not able to
	guarantee that all other processes were killed!

Only 1 process:

	When I load the trace in paraver, only one process is shown. If I use
	NANOS6_EXTRAE_AS_THREADS=1, then multiple threads are shown, but only
	for one process.

