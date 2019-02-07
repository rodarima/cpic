#pragma GCC visibility push(default)
#pragma GCC visibility pop
#pragma GCC visibility push(default)
#pragma GCC visibility pop
#pragma GCC visibility push(default)
#pragma GCC visibility pop
#pragma GCC visibility push(default)
#pragma GCC visibility pop
#pragma GCC visibility push(default)
#pragma GCC visibility pop
#pragma GCC visibility push(default)
#pragma GCC visibility push(default)
#pragma GCC visibility pop
#pragma GCC visibility push(default)
#pragma GCC visibility pop
#pragma GCC visibility push(default)
#pragma GCC visibility pop
#pragma GCC visibility push(default)
#pragma GCC visibility pop
#pragma GCC visibility push(default)
#pragma GCC visibility pop
#pragma GCC visibility push(default)
#pragma GCC visibility pop
#pragma GCC visibility push(default)
#pragma GCC visibility pop
#pragma GCC visibility push(default)
#pragma GCC visibility pop
static __inline unsigned int __bswap_32(unsigned int __bsx)
{
  return __builtin_bswap32(__bsx);
}
typedef unsigned long int __uint64_t;
static __inline __uint64_t __bswap_64(__uint64_t __bsx)
{
  return __builtin_bswap64(__bsx);
}
struct _IO_FILE_plus;
extern struct _IO_FILE_plus _IO_2_1_stdin_;
extern struct _IO_FILE_plus _IO_2_1_stdout_;
extern struct _IO_FILE_plus _IO_2_1_stderr_;
struct _IO_FILE;
extern struct _IO_FILE *stdin;
extern struct _IO_FILE *stdout;
extern struct _IO_FILE *stderr;
extern int sys_nerr;
extern const char *const sys_errlist[];
struct specie;
typedef struct specie specie_t;
typedef unsigned long int size_t;
extern void *malloc(size_t __size) __attribute__((__nothrow__)) __attribute__((__leaf__)) __attribute__((__malloc__));
struct mcc_struct_anon_6;
typedef struct mcc_struct_anon_6 mat_t;
struct particle;
struct block;
typedef struct block block_t;
struct  specie
{
  int dim;
  int *shape;
  float q;
  float m;
  mat_t *E;
  mat_t *B;
  mat_t *J;
  mat_t *rho;
  float dt;
  float t;
  float dx;
  float C;
  float e0;
  int nparticles;
  struct particle *particles;
  int nblocks;
  int blocksize;
  int nnodes;
  block_t *blocks;
};
struct  particle
{
  int i;
  float x;
  float u;
  float E;
  float B;
  float J;
  struct particle *next;
  struct particle *prev;
};
typedef struct particle particle_t;
mat_t *mat_init(int dim, int *shape, float v);
specie_t *specie_alloc(int dim, int *shape, int nparticles)
{
  specie_t *s;
  s = malloc(sizeof(specie_t));
  (*s).dim = dim;
  (*s).nparticles = nparticles;
  (*s).particles = malloc(nparticles * sizeof(particle_t));
  (*s).E = mat_init(dim, shape, 0.00000000000000000000000000000000000000000000000000000e+00);
  (*s).B = mat_init(dim, shape, 0.00000000000000000000000000000000000000000000000000000e+00);
  (*s).J = mat_init(dim, shape, 0.00000000000000000000000000000000000000000000000000000e+00);
  (*s).rho = mat_init(dim, shape, 0.00000000000000000000000000000000000000000000000000000e+00);
  return s;
}
struct  mcc_struct_anon_6
{
  float *data;
  int *shape;
  int dim;
  int size;
};
int field_init(mat_t *f)
{
  int i;
  for (i = 0; i < (*f).size; i++)
    {
      ((float *)(*f).data)[i] = 0.00000000000000000000000000000000000000000000000000000e+00;
    }
  return 0;
}
extern int rand(void) __attribute__((__nothrow__)) __attribute__((__leaf__));
int particles_init(specie_t *s)
{
  int i;
  particle_t *p;
  int total_nodes = (*s).blocksize * (*s).nblocks;
  for (i = 0; i < (*s).nparticles; i++)
    {
      p = &(*s).particles[i];
      (*p).i = i;
      (*p).x = (float)rand() / 2147483647 * total_nodes * (*s).dx;
      (*p).u = 2.00000000000000000000000000000000000000000000000000000e+00 * (i % 2 - 5.00000000000000000000000000000000000000000000000000000e-01) * 5.00000000000000000000000000000000000000000000000000000e-01 * (*s).C;
      (*p).E = 0.00000000000000000000000000000000000000000000000000000e+00;
      (*p).J = 0.00000000000000000000000000000000000000000000000000000e+00;
    }
  return 0;
}
extern int printf(const char *__restrict __format, ...);
int blocks_init(specie_t *s);
specie_t *specie_init()
{
  specie_t *s;
  int dim = 1;
  int shape = 400;
  int nparticles = 100000;
  s = specie_alloc(dim, &shape, nparticles);
  (*s).shape = malloc(sizeof(int) * 1);
  (*s).shape[0] = shape;
  (*s).t = 0.00000000000000000000000000000000000000000000000000000e+00;
  (*s).C = 2.99792458000000000000000000000000000000000000000000000e+08;
  (*s).q =  -1.60217662000000009010239718003012684178391545553893783e-19 * 1.00000000000000002081668171172168513294309377670288086e-03;
  (*s).m = 9.10938355999999975159749255595945385778904388621038995e-31;
  (*s).e0 = 8.85000000000000045813426548913611588420558007328509120e-12;
  (*s).dt = 1.00000000000000002092256083012847267532663408928783610e-08;
  (*s).dx = 1.00000000000000000000000000000000000000000000000000000e+01 * ((*s).C / 2.00000000000000000000000000000000000000000000000000000e+00) * (*s).dt * 1.00000000000000000000000000000000000000000000000000000e+02;
  (*s).nblocks = shape;
  (*s).blocksize = 10000;
  (*s).nnodes = (*s).nblocks * (*s).blocksize;
  printf("%d %d %10.e %10.e\n", nparticles, (*s).nnodes, (*s).dx, (*s).dt);
  particles_init(s);
  blocks_init(s);
  return s;
}
void specie_step(specie_t *s)
{
  (*s).t += (*s).dt;
}
int specie_print(specie_t *s)
{
  int i;
  particle_t *p;
  for (i = 0; i < (*s).nparticles; i++)
    {
      p = &(*s).particles[i];
      printf("%d %10.3e %10.3e\n", (*p).i, (*p).x, (*p).u);
    }
  return 0;
}
