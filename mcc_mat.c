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
struct mcc_struct_anon_6;
typedef struct mcc_struct_anon_6 mat_t;
typedef unsigned long int size_t;
extern void *malloc(size_t __size) __attribute__((__nothrow__)) __attribute__((__leaf__)) __attribute__((__malloc__));
struct  mcc_struct_anon_6
{
  float *data;
  int *shape;
  int dim;
  int size;
};
mat_t *mat_alloc(int dim, int *shape)
{
  mat_t *m;
  int size;
  int i;
  m = malloc(sizeof(mat_t));
  (*m).dim = dim;
  size = 1;
  for (i = 0; i < dim; i++)
    {
      size *= shape[i];
    }
  (*m).size = size;
  (*m).data = malloc(sizeof(float) * size);
  return m;
}
mat_t *mat_init(int dim, int *shape, float v)
{
  mat_t *m;
  int i;
  m = mat_alloc(dim, shape);
  for (i = 0; i < (*m).size; i++)
    {
      (*m).data[i] = v;
    }
  return m;
}
void mat_set1d(mat_t *m, int pos, float v)
{
  ((float *)m)[pos] = v;
}
mat_t *vec_init(int size, float v)
{
  mat_t *m;
  int i;
  m = mat_alloc(1, &size);
  for (i = 0; i < (*m).size; i++)
    {
      (*m).data[i] = v;
    }
  return m;
}
