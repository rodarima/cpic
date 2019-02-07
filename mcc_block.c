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
extern int signgam;
enum mcc_enum_anon_4
{
  _IEEE_ =  -1,
  _SVID_ = 0,
  _XOPEN_ = 1,
  _POSIX_ = 2,
  _ISOC_ = 3
};
typedef enum mcc_enum_anon_4 _LIB_VERSION_TYPE;
extern _LIB_VERSION_TYPE _LIB_VERSION;
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
struct block;
typedef struct block block_t;
struct particle;
typedef struct particle particle_t;
typedef unsigned long int size_t;
struct mcc_struct_anon_6;
typedef struct mcc_struct_anon_6 mat_t;
struct  block
{
  int i;
  float x;
  mat_t *E;
  mat_t *J;
  float rE;
  float lE;
  float rJ;
  float lJ;
  particle_t *particles;
  particle_t *left;
  particle_t *right;
};
mat_t *vec_init(int size, float v);
int block_init(block_t *b, size_t blocksize, particle_t *particles, size_t n)
{
  (*b).E = vec_init(blocksize, 0.00000000000000000000000000000000000000000000000000000e+00);
  (*b).J = vec_init(blocksize, 0.00000000000000000000000000000000000000000000000000000e+00);
  return 0;
}
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
void block_add_particle(block_t *b, particle_t *p)
{
  do
    {
      if ((*b).particles)
        {
          (*p).prev = (*(*b).particles).prev;
          (*(*(*b).particles).prev).next = p;
          (*(*b).particles).prev = p;
          (*p).next = (void *)0;
        }
      else
        {
          (*b).particles = p;
          (*(*b).particles).prev = (*b).particles;
          (*(*b).particles).next = (void *)0;
        }
    }
  while (0);
}
struct specie;
typedef struct specie specie_t;
extern void *calloc(size_t __nmemb, size_t __size) __attribute__((__nothrow__)) __attribute__((__leaf__)) __attribute__((__malloc__));
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
extern float ceilf(float __x) __attribute__((__nothrow__)) __attribute__((__leaf__)) __attribute__((__const__));
extern double floor(double __x) __attribute__((__nothrow__)) __attribute__((__leaf__)) __attribute__((__const__));
int blocks_init(specie_t *s)
{
  block_t *blocks;
  size_t nparticles_block;
  size_t i;
  size_t j;
  block_t *b;
  particle_t *p;
  blocks = calloc((*s).nblocks, sizeof(block_t));
  nparticles_block = (size_t)ceilf((float)(*s).nparticles / (float)(*s).nblocks);
  for ((i = 0, j = 0); i < (*s).nblocks; i++)
    {
      b = &blocks[i];
      (*b).i = i;
      (*b).E = vec_init((*s).blocksize, 0.00000000000000000000000000000000000000000000000000000e+00);
      (*b).J = vec_init((*s).blocksize, 0.00000000000000000000000000000000000000000000000000000e+00);
      (*b).x = i * (*s).dx * (*s).blocksize;
    }
  for (i = 0; i < (*s).nparticles; i++)
    {
      p = &(*s).particles[i];
      j = (int)floor((*p).x / ((*s).dx * (*s).blocksize));
      block_add_particle(&blocks[j], p);
    }
  (*s).blocks = blocks;
  return 0;
}
extern int printf(const char *__restrict __format, ...);
void block_print(block_t *block)
{
  particle_t *p;
  size_t i;
  p = (*block).particles;
  for ((p = (*block).particles, i = 0); p; (p = (*p).next, i++))
    {
      printf("%zu %10.3e %10.3e\n", i, (*p).x, (*p).u);
    }
}
void blocks_print(block_t *blocks, size_t n)
{
  size_t i;
  for (i = 0; i < n; i++)
    {
      ;
      block_print(&blocks[i]);
    }
}
int block_particle_J(specie_t *s, block_t *b)
{
  particle_t *p;
  for (p = (*b).particles; p; p = (*p).next)
    {
      (*p).J = (*s).q * (*p).u / (*s).dt;
    }
  return 0;
}
struct  mcc_struct_anon_6
{
  float *data;
  int *shape;
  int dim;
  int size;
};
extern void *memset(void *__s, int __c, size_t __n) __attribute__((__nothrow__)) __attribute__((__leaf__)) __attribute__((__nonnull__(1)));
extern void exit(int __status) __attribute__((__nothrow__)) __attribute__((__leaf__)) __attribute__((__noreturn__));
extern void __assert_fail(const char *__assertion, const char *__file, unsigned int __line, const char *__function) __attribute__((__nothrow__)) __attribute__((__leaf__)) __attribute__((__noreturn__));
int block_field_J(specie_t *s, block_t *b)
{
  particle_t *p;
  float px;
  float deltax;
  int j0;
  float w1;
  float w0;
  int j1;
  static const char __MERCURIUM_PRETTY_FUNCTION__[41L] = "int block_field_J(specie_t *, block_t *)";
  float *J = (*(*b).J).data;
  int size = (*(*b).J).size;
  float dx = (*s).dx;
  float x0 = (*b).x;
  float x1 = (*b).x + dx * size;
  float xhalf = (x0 + x1) / 2.00000000000000000000000000000000000000000000000000000e+00;
  memset(J, 0, sizeof(float) * size);
  (*b).rJ = 0.00000000000000000000000000000000000000000000000000000e+00;
  ;
  for (p = (*b).particles; p; p = (*p).next)
    {
      px = (*p).x;
      deltax = px - x0;
      if (!(x0 <= px && px <= x1))
        {
          ;
          exit(1);
        }
      j0 = (int)floor(deltax / (*s).dx);
      j0 >= 0 ? (void)0 : __assert_fail("j0 >= 0", "block.c", 144, __MERCURIUM_PRETTY_FUNCTION__);
      j0 < size ? (void)0 : __assert_fail("j0 < size", "block.c", 145, __MERCURIUM_PRETTY_FUNCTION__);
      w1 = deltax / (*s).dx;
      w0 = 1.00000000000000000000000000000000000000000000000000000e+00 - w1;
      if (px >= x1 - dx)
        {
          J[j0] += w0 * (*p).J;
          (*b).rJ += w1 * (*p).J;
          ;
        }
      else
        {
          j1 = (j0 + 1) % size;
          j1 < size ? (void)0 : __assert_fail("j1 < size", "block.c", 162, __MERCURIUM_PRETTY_FUNCTION__);
          J[j0] += w0 * (*p).J;
          J[j1] += w1 * (*p).J;
          ;
        }
    }
  return 0;
}
int block_comm_field_J(block_t *dst, block_t *left)
{
  (*(*dst).J).data[0] += (*left).rJ;
  (*left).rJ = 0.00000000000000000000000000000000000000000000000000000e+00;
  return 0;
}
int block_field_E(specie_t *s, block_t *b)
{
  int i;
  int size = (*(*b).E).size;
  float *E = (*(*b).E).data;
  float *J = (*(*b).J).data;
  float coef =  -(*s).dt / (*s).e0;
  for (i = 0; i < size; i++)
    {
      E[i] += coef * J[i];
      ;
    }
  return 0;
}
int block_comm_field_E(block_t *dst, block_t *right)
{
  (*dst).rE = (*(*right).E).data[0];
  return 0;
}
int block_particle_E(specie_t *s, block_t *b)
{
  particle_t *p;
  float px;
  float deltax;
  int j0;
  int j1;
  float w1;
  float w0;
  static const char __MERCURIUM_PRETTY_FUNCTION__[44L] = "int block_particle_E(specie_t *, block_t *)";
  float *E = (*(*b).E).data;
  int size = (*(*b).E).size;
  float dx = (*s).dx;
  float x0 = (*b).x;
  float x1 = (*b).x + dx * size;
  float xhalf = (x0 + x1) / 2.00000000000000000000000000000000000000000000000000000e+00;
  size = (*(*s).E).size;
  ;
  for (p = (*b).particles; p; p = (*p).next)
    {
      px = (*p).x;
      deltax = px - x0;
      if (!(x0 <= px && px <= x1))
        {
          ;
          exit(1);
        }
      j0 = (int)floor(deltax / (*s).dx);
      j1 = j0 + 1;
      j0 >= 0 ? (void)0 : __assert_fail("j0 >= 0", "block.c", 243, __MERCURIUM_PRETTY_FUNCTION__);
      w1 = deltax / (*s).dx;
      w0 = 1.00000000000000000000000000000000000000000000000000000e+00 - w1;
      (*p).E = w0 * E[j0];
      if (px >= x1 - dx)
        {
          (*p).E += w1 * (*b).rE;
        }
      else
        {
          (*p).E += w1 * E[j1];
        }
      ;
    }
  return 0;
}
extern double fabs(double __x) __attribute__((__nothrow__)) __attribute__((__leaf__)) __attribute__((__const__));
typedef struct _IO_FILE FILE;
extern int fprintf(FILE *__restrict __stream, const char *__restrict __format, ...);
extern struct _IO_FILE *stderr;
int block_particle_x(specie_t *s, block_t *b)
{
  particle_t *p;
  float delta_u;
  float delta_x;
  float coef =  -(*s).dt / (*s).e0;
  float dt = (*s).dt;
  for (p = (*b).particles; p; p = (*p).next)
    {
      delta_u = dt * (*s).q * (*p).E / (*s).m;
      ;
      (*p).u += delta_u;
      delta_x = dt * (*p).u;
      if (fabs(delta_x) > (*s).dx)
        {
          fprintf(stderr, "Particle %d at x=%.3e has exceeded dx with delta_x=%.3e\n", (*p).i, (*p).x, delta_x);
          ;
          fprintf(stderr, "Please, reduce dt=%.3e or increase dx=%.3e\n", (*s).dt, (*s).dx);
          ;
          exit(1);
        }
      (*p).x += delta_x;
    }
  return 0;
}
int block_comm_particles(specie_t *s, block_t *left, block_t *b, block_t *right)
{
  particle_t *p;
  particle_t *tmp;
  float px;
  static const char __MERCURIUM_PRETTY_FUNCTION__[70L] = "int block_comm_particles(specie_t *, block_t *, block_t *, block_t *)";
  float x0 = (*b).x;
  float x1 = (*b).x + (*s).dx * (*s).blocksize;
  float max_x = (*s).dx * (*s).blocksize * (*s).nblocks;
  ;
  for (p = (*b).particles; p && ((tmp = (*p).next, 1)); p = tmp)
    {
      px = (*p).x;
      if (px < x0)
        {
          ;
          do
            {
              (*b).particles != (void *)0 ? (void)0 : __assert_fail("(b->particles) != ((void *)0)", "block.c", 327, __MERCURIUM_PRETTY_FUNCTION__);
              (*p).prev != (void *)0 ? (void)0 : __assert_fail("(p)->prev != ((void *)0)", "block.c", 327, __MERCURIUM_PRETTY_FUNCTION__);
              if ((*p).prev == p)
                {
                  (*b).particles = (void *)0;
                }
              else
                {
                  if (p == (*b).particles)
                    {
                      (*(*p).next).prev = (*p).prev;
                      (*b).particles = (*p).next;
                    }
                  else
                    {
                      (*(*p).prev).next = (*p).next;
                      if ((*p).next)
                        {
                          (*(*p).next).prev = (*p).prev;
                        }
                      else
                        {
                          (*(*b).particles).prev = (*p).prev;
                        }
                    }
                }
            }
          while (0);
          do
            {
              if ((*left).particles)
                {
                  (*p).prev = (*(*left).particles).prev;
                  (*(*(*left).particles).prev).next = p;
                  (*(*left).particles).prev = p;
                  (*p).next = (void *)0;
                }
              else
                {
                  (*left).particles = p;
                  (*(*left).particles).prev = (*left).particles;
                  (*(*left).particles).next = (void *)0;
                }
            }
          while (0);
        }
      else
        {
          if (x1 <= px)
            {
              ;
              do
                {
                  (*b).particles != (void *)0 ? (void)0 : __assert_fail("(b->particles) != ((void *)0)", "block.c", 334, __MERCURIUM_PRETTY_FUNCTION__);
                  (*p).prev != (void *)0 ? (void)0 : __assert_fail("(p)->prev != ((void *)0)", "block.c", 334, __MERCURIUM_PRETTY_FUNCTION__);
                  if ((*p).prev == p)
                    {
                      (*b).particles = (void *)0;
                    }
                  else
                    {
                      if (p == (*b).particles)
                        {
                          (*(*p).next).prev = (*p).prev;
                          (*b).particles = (*p).next;
                        }
                      else
                        {
                          (*(*p).prev).next = (*p).next;
                          if ((*p).next)
                            {
                              (*(*p).next).prev = (*p).prev;
                            }
                          else
                            {
                              (*(*b).particles).prev = (*p).prev;
                            }
                        }
                    }
                }
              while (0);
              do
                {
                  if ((*right).particles)
                    {
                      (*p).prev = (*(*right).particles).prev;
                      (*(*(*right).particles).prev).next = p;
                      (*(*right).particles).prev = p;
                      (*p).next = (void *)0;
                    }
                  else
                    {
                      (*right).particles = p;
                      (*(*right).particles).prev = (*right).particles;
                      (*(*right).particles).next = (void *)0;
                    }
                }
              while (0);
            }
          else
            {
              ;
            }
        }
      if ((*p).x >= max_x)
        {
          ;
          (*p).x -= max_x;
        }
      else
        {
          if ((*p).x < 0.00000000000000000000000000000000000000000000000000000e+00)
            {
              ;
              (*p).x += max_x;
            }
        }
      if ((*p).x < 0.00000000000000000000000000000000000000000000000000000e+00 || (*p).x > max_x)
        {
          fprintf(stderr, "Particle %d is at x=%.3e with max_x=%10.3e\n", (*p).i, (*p).x, max_x);
          ;
          exit(1);
        }
    }
  return 0;
}
int block_print_particles(specie_t *s, block_t *b)
{
  particle_t *p;
  for (p = (*b).particles; p; p = (*p).next)
    {
      printf("%d %10.3e %10.3e\n", (*p).i, (*p).x, (*p).u);
    }
  return 0;
}
