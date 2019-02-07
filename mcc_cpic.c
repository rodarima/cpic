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
extern char **__environ;
extern char *optarg;
extern int optind;
extern int opterr;
extern int optopt;
struct  mcc_struct_anon_4
{
  const char *invocation_source;
};
typedef __attribute__((aligned(64))) struct mcc_struct_anon_4 nanos6_task_invocation_info_t;
static nanos6_task_invocation_info_t task_invocation_info_field_E_0 = {.invocation_source = "cpic.c:26:9"};
struct mcc_struct_anon_1;
typedef struct mcc_struct_anon_1 nanos6_address_translation_entry_t;
struct mcc_struct_anon_0;
typedef struct mcc_struct_anon_0 nanos6_task_constraints_t;
struct  mcc_struct_anon_2
{
  int device_type_id;
  void (*run)(void *, void *, nanos6_address_translation_entry_t *);
  void (*get_constraints)(void *, nanos6_task_constraints_t *);
  const char *task_label;
  const char *declaration_source;
  void (*run_wrapper)(void *, void *, nanos6_address_translation_entry_t *);
};
typedef struct mcc_struct_anon_2 nanos6_task_implementation_info_t;
enum mcc_enum_anon_1
{
  nanos6_host_device = 0,
  nanos6_cuda_device = 1,
  nanos6_opencl_device = 2,
  nanos6_device_type_num = 3
};
struct block;
typedef struct block block_t;
struct specie;
typedef struct specie specie_t;
struct  nanos_task_args_field_E_0
{
  block_t *b;
  specie_t *s;
};
static void nanos6_ol_task_region_field_E_0(struct nanos_task_args_field_E_0 *const arg, void *device_env, nanos6_address_translation_entry_t *address_translation_table);
static nanos6_task_implementation_info_t implementations_var_field_E_0[1L] = {{.device_type_id = nanos6_host_device, .run = (void (*)(void *, void *, nanos6_address_translation_entry_t *))nanos6_ol_task_region_field_E_0, .get_constraints = 0, .task_label = "block_field_E", .declaration_source = "cpic.c:26:9", .run_wrapper = 0}};
typedef long int nanos6_priority_t;
typedef unsigned long int size_t;
struct  mcc_struct_anon_3
{
  int num_symbols;
  void (*register_depinfo)(void *, void *);
  void (*get_priority)(void *, nanos6_priority_t *);
  const char *type_identifier;
  int implementation_count;
  nanos6_task_implementation_info_t *implementations;
  void (*destroy_args_block)(void *);
  void (*duplicate_args_block)(const void *, void *);
  void (**reduction_initializers)(void *, void *, size_t);
  void (**reduction_combiners)(void *, void *, size_t);
};
typedef __attribute__((aligned(64))) struct mcc_struct_anon_3 nanos6_task_info_t;
static void nanos6_ol_deps_field_E_0(struct nanos_task_args_field_E_0 *const arg, void *handler);
static nanos6_task_info_t task_info_var_field_E_0 = {.num_symbols =  -1, .register_depinfo = (void (*)(void *, void *))nanos6_ol_deps_field_E_0, .get_priority = 0, .type_identifier = 0, .implementation_count = 1, .implementations = implementations_var_field_E_0, .destroy_args_block = 0, .duplicate_args_block = 0, .reduction_initializers = 0, .reduction_combiners = 0};
static nanos6_task_invocation_info_t task_invocation_info_field_E_1 = {.invocation_source = "cpic.c:41:9"};
struct  nanos_task_args_field_E_1
{
  block_t *rb;
  block_t *b;
};
static void nanos6_ol_task_region_field_E_1(struct nanos_task_args_field_E_1 *const arg, void *device_env, nanos6_address_translation_entry_t *address_translation_table);
static nanos6_task_implementation_info_t implementations_var_field_E_1[1L] = {{.device_type_id = nanos6_host_device, .run = (void (*)(void *, void *, nanos6_address_translation_entry_t *))nanos6_ol_task_region_field_E_1, .get_constraints = 0, .task_label = "block_comm_field_E", .declaration_source = "cpic.c:41:9", .run_wrapper = 0}};
static void nanos6_ol_deps_field_E_1(struct nanos_task_args_field_E_1 *const arg, void *handler);
static nanos6_task_info_t task_info_var_field_E_1 = {.num_symbols =  -1, .register_depinfo = (void (*)(void *, void *))nanos6_ol_deps_field_E_1, .get_priority = 0, .type_identifier = 0, .implementation_count = 1, .implementations = implementations_var_field_E_1, .destroy_args_block = 0, .duplicate_args_block = 0, .reduction_initializers = 0, .reduction_combiners = 0};
struct mcc_struct_anon_6;
typedef struct mcc_struct_anon_6 mat_t;
struct particle;
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
typedef struct particle particle_t;
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
int nanos6_in_final(void);
typedef unsigned int __useconds_t;
extern int usleep(__useconds_t __useconds);
int block_field_E(specie_t *s, block_t *b);
void nanos6_create_task(nanos6_task_info_t *task_info, nanos6_task_invocation_info_t *task_invocation_info, size_t args_block_size, void **args_block_pointer, void **task_pointer, size_t flags);
void nanos6_submit_task(void *task);
int block_comm_field_E(block_t *from, block_t *to);
int field_E(specie_t *s)
{
  int i;
  block_t *b;
  int ri;
  block_t *rb;
  for (i = 0; i < (*s).nblocks; i++)
    {
      b = &(*s).blocks[i];
      if (nanos6_in_final() != 0LU)
        {
          {
            usleep(100);
            block_field_E(s, b);
          }
        }
      else
        {
          unsigned long int task_flags_0;
          struct nanos_task_args_field_E_0 *nanos_data_env_0;
          void *nanos_task_ptr_0;
          task_flags_0 = (((0 != 0) << 0U | (0 != 0) << 1U) | (0 != 0) << 2U) | (0 != 0) << 3U;
          nanos6_create_task(&task_info_var_field_E_0, &task_invocation_info_field_E_0, 16LU + 0, (void **)&nanos_data_env_0, &nanos_task_ptr_0, task_flags_0);
          (*nanos_data_env_0).b = b;
          (*nanos_data_env_0).s = s;
          nanos6_submit_task(nanos_task_ptr_0);
        }
    }
  for (i = 0; i < (*s).nblocks; i++)
    {
      ri = (i + 1) % (*s).nblocks;
      b = &(*s).blocks[i];
      rb = &(*s).blocks[ri];
      if (nanos6_in_final() != 0LU)
        {
          {
            usleep(100);
            block_comm_field_E(b, rb);
          }
        }
      else
        {
          unsigned long int task_flags_1;
          struct nanos_task_args_field_E_1 *nanos_data_env_1;
          void *nanos_task_ptr_1;
          task_flags_1 = (((0 != 0) << 0U | (0 != 0) << 1U) | (0 != 0) << 2U) | (0 != 0) << 3U;
          nanos6_create_task(&task_info_var_field_E_1, &task_invocation_info_field_E_1, 16LU + 0, (void **)&nanos_data_env_1, &nanos_task_ptr_1, task_flags_1);
          (*nanos_data_env_1).rb = rb;
          (*nanos_data_env_1).b = b;
          nanos6_submit_task(nanos_task_ptr_1);
        }
    }
  return 0;
}
static void nanos6_unpacked_deps_field_E_1(block_t **const rb, block_t **const b, void *handler);
static void nanos6_ol_deps_field_E_1(struct nanos_task_args_field_E_1 *const arg, void *handler)
{
  nanos6_unpacked_deps_field_E_1(&((*arg).rb), &((*arg).b), handler);
}
void nanos6_register_region_read_depinfo1(void *handler, int symbol_index, const char *region_text, void *base_address, long int dim1size, long int dim1start, long int dim1end);
void nanos6_register_region_readwrite_depinfo1(void *handler, int symbol_index, const char *region_text, void *base_address, long int dim1size, long int dim1start, long int dim1end);
static void nanos6_unpacked_deps_field_E_1(block_t **const rb, block_t **const b, void *handler)
{
  nanos6_register_region_read_depinfo1(handler,  -1, "*rb", (void *)(*rb), sizeof(block_t), 0, sizeof(block_t));
  nanos6_register_region_readwrite_depinfo1(handler,  -1, "*b", (void *)(*b), sizeof(block_t), 0, sizeof(block_t));
  ;
}
static void nanos6_unpacked_task_region_field_E_1(block_t **const rb, block_t **const b, void *device_env, nanos6_address_translation_entry_t *address_translation_table);
static void nanos6_ol_task_region_field_E_1(struct nanos_task_args_field_E_1 *const arg, void *device_env, nanos6_address_translation_entry_t *address_translation_table)
{
  nanos6_unpacked_task_region_field_E_1(&((*arg).rb), &((*arg).b), device_env, address_translation_table);
}
static void nanos6_unpacked_deps_field_E_0(block_t **const b, specie_t **const s, void *handler);
static void nanos6_ol_deps_field_E_0(struct nanos_task_args_field_E_0 *const arg, void *handler)
{
  nanos6_unpacked_deps_field_E_0(&((*arg).b), &((*arg).s), handler);
}
static void nanos6_unpacked_deps_field_E_0(block_t **const b, specie_t **const s, void *handler)
{
  nanos6_register_region_readwrite_depinfo1(handler,  -1, "*b", (void *)(*b), sizeof(block_t), 0, sizeof(block_t));
  ;
}
static void nanos6_unpacked_task_region_field_E_0(block_t **const b, specie_t **const s, void *device_env, nanos6_address_translation_entry_t *address_translation_table);
static void nanos6_ol_task_region_field_E_0(struct nanos_task_args_field_E_0 *const arg, void *device_env, nanos6_address_translation_entry_t *address_translation_table)
{
  nanos6_unpacked_task_region_field_E_0(&((*arg).b), &((*arg).s), device_env, address_translation_table);
}
static nanos6_task_invocation_info_t task_invocation_info_field_J_2 = {.invocation_source = "cpic.c:64:9"};
struct  nanos_task_args_field_J_2
{
  block_t *b;
  specie_t *s;
};
static void nanos6_ol_task_region_field_J_2(struct nanos_task_args_field_J_2 *const arg, void *device_env, nanos6_address_translation_entry_t *address_translation_table);
static nanos6_task_implementation_info_t implementations_var_field_J_2[1L] = {{.device_type_id = nanos6_host_device, .run = (void (*)(void *, void *, nanos6_address_translation_entry_t *))nanos6_ol_task_region_field_J_2, .get_constraints = 0, .task_label = "block_field_J", .declaration_source = "cpic.c:64:9", .run_wrapper = 0}};
static void nanos6_ol_deps_field_J_2(struct nanos_task_args_field_J_2 *const arg, void *handler);
static nanos6_task_info_t task_info_var_field_J_2 = {.num_symbols =  -1, .register_depinfo = (void (*)(void *, void *))nanos6_ol_deps_field_J_2, .get_priority = 0, .type_identifier = 0, .implementation_count = 1, .implementations = implementations_var_field_J_2, .destroy_args_block = 0, .duplicate_args_block = 0, .reduction_initializers = 0, .reduction_combiners = 0};
static nanos6_task_invocation_info_t task_invocation_info_field_J_3 = {.invocation_source = "cpic.c:76:9"};
struct  nanos_task_args_field_J_3
{
  block_t *b;
  block_t *lb;
};
static void nanos6_ol_task_region_field_J_3(struct nanos_task_args_field_J_3 *const arg, void *device_env, nanos6_address_translation_entry_t *address_translation_table);
static nanos6_task_implementation_info_t implementations_var_field_J_3[1L] = {{.device_type_id = nanos6_host_device, .run = (void (*)(void *, void *, nanos6_address_translation_entry_t *))nanos6_ol_task_region_field_J_3, .get_constraints = 0, .task_label = "block_comm_field_J", .declaration_source = "cpic.c:76:9", .run_wrapper = 0}};
static void nanos6_ol_deps_field_J_3(struct nanos_task_args_field_J_3 *const arg, void *handler);
static nanos6_task_info_t task_info_var_field_J_3 = {.num_symbols =  -1, .register_depinfo = (void (*)(void *, void *))nanos6_ol_deps_field_J_3, .get_priority = 0, .type_identifier = 0, .implementation_count = 1, .implementations = implementations_var_field_J_3, .destroy_args_block = 0, .duplicate_args_block = 0, .reduction_initializers = 0, .reduction_combiners = 0};
int block_field_J(specie_t *s, block_t *b);
int block_comm_field_J(block_t *from, block_t *to);
int field_J(specie_t *s)
{
  int i;
  block_t *b;
  int li;
  block_t *lb;
  for (i = 0; i < (*s).nblocks; i++)
    {
      b = &(*s).blocks[i];
      if (nanos6_in_final() != 0LU)
        {
          block_field_J(s, b);
        }
      else
        {
          unsigned long int task_flags_2;
          struct nanos_task_args_field_J_2 *nanos_data_env_2;
          void *nanos_task_ptr_2;
          task_flags_2 = (((0 != 0) << 0U | (0 != 0) << 1U) | (0 != 0) << 2U) | (0 != 0) << 3U;
          nanos6_create_task(&task_info_var_field_J_2, &task_invocation_info_field_J_2, 16LU + 0, (void **)&nanos_data_env_2, &nanos_task_ptr_2, task_flags_2);
          (*nanos_data_env_2).b = b;
          (*nanos_data_env_2).s = s;
          nanos6_submit_task(nanos_task_ptr_2);
        }
    }
  for (i = 0; i < (*s).nblocks; i++)
    {
      li = ((*s).nblocks + i - 1) % (*s).nblocks;
      b = &(*s).blocks[i];
      lb = &(*s).blocks[li];
      if (nanos6_in_final() != 0LU)
        {
          block_comm_field_J(b, lb);
        }
      else
        {
          unsigned long int task_flags_3;
          struct nanos_task_args_field_J_3 *nanos_data_env_3;
          void *nanos_task_ptr_3;
          task_flags_3 = (((0 != 0) << 0U | (0 != 0) << 1U) | (0 != 0) << 2U) | (0 != 0) << 3U;
          nanos6_create_task(&task_info_var_field_J_3, &task_invocation_info_field_J_3, 16LU + 0, (void **)&nanos_data_env_3, &nanos_task_ptr_3, task_flags_3);
          (*nanos_data_env_3).b = b;
          (*nanos_data_env_3).lb = lb;
          nanos6_submit_task(nanos_task_ptr_3);
        }
    }
  return 0;
}
static void nanos6_unpacked_deps_field_J_3(block_t **const b, block_t **const lb, void *handler);
static void nanos6_ol_deps_field_J_3(struct nanos_task_args_field_J_3 *const arg, void *handler)
{
  nanos6_unpacked_deps_field_J_3(&((*arg).b), &((*arg).lb), handler);
}
static void nanos6_unpacked_deps_field_J_3(block_t **const b, block_t **const lb, void *handler)
{
  nanos6_register_region_readwrite_depinfo1(handler,  -1, "*b", (void *)(*b), sizeof(block_t), 0, sizeof(block_t));
  nanos6_register_region_readwrite_depinfo1(handler,  -1, "*lb", (void *)(*lb), sizeof(block_t), 0, sizeof(block_t));
  ;
}
static void nanos6_unpacked_task_region_field_J_3(block_t **const b, block_t **const lb, void *device_env, nanos6_address_translation_entry_t *address_translation_table);
static void nanos6_ol_task_region_field_J_3(struct nanos_task_args_field_J_3 *const arg, void *device_env, nanos6_address_translation_entry_t *address_translation_table)
{
  nanos6_unpacked_task_region_field_J_3(&((*arg).b), &((*arg).lb), device_env, address_translation_table);
}
static void nanos6_unpacked_deps_field_J_2(block_t **const b, specie_t **const s, void *handler);
static void nanos6_ol_deps_field_J_2(struct nanos_task_args_field_J_2 *const arg, void *handler)
{
  nanos6_unpacked_deps_field_J_2(&((*arg).b), &((*arg).s), handler);
}
static void nanos6_unpacked_deps_field_J_2(block_t **const b, specie_t **const s, void *handler)
{
  nanos6_register_region_readwrite_depinfo1(handler,  -1, "*b", (void *)(*b), sizeof(block_t), 0, sizeof(block_t));
  ;
}
static void nanos6_unpacked_task_region_field_J_2(block_t **const b, specie_t **const s, void *device_env, nanos6_address_translation_entry_t *address_translation_table);
static void nanos6_ol_task_region_field_J_2(struct nanos_task_args_field_J_2 *const arg, void *device_env, nanos6_address_translation_entry_t *address_translation_table)
{
  nanos6_unpacked_task_region_field_J_2(&((*arg).b), &((*arg).s), device_env, address_translation_table);
}
static nanos6_task_invocation_info_t task_invocation_info_particle_E_4 = {.invocation_source = "cpic.c:94:9"};
struct  nanos_task_args_particle_E_4
{
  block_t *b;
  specie_t *s;
};
static void nanos6_ol_task_region_particle_E_4(struct nanos_task_args_particle_E_4 *const arg, void *device_env, nanos6_address_translation_entry_t *address_translation_table);
static nanos6_task_implementation_info_t implementations_var_particle_E_4[1L] = {{.device_type_id = nanos6_host_device, .run = (void (*)(void *, void *, nanos6_address_translation_entry_t *))nanos6_ol_task_region_particle_E_4, .get_constraints = 0, .task_label = "block_particle_E", .declaration_source = "cpic.c:94:9", .run_wrapper = 0}};
static void nanos6_ol_deps_particle_E_4(struct nanos_task_args_particle_E_4 *const arg, void *handler);
static nanos6_task_info_t task_info_var_particle_E_4 = {.num_symbols =  -1, .register_depinfo = (void (*)(void *, void *))nanos6_ol_deps_particle_E_4, .get_priority = 0, .type_identifier = 0, .implementation_count = 1, .implementations = implementations_var_particle_E_4, .destroy_args_block = 0, .duplicate_args_block = 0, .reduction_initializers = 0, .reduction_combiners = 0};
int block_particle_E(specie_t *s, block_t *b);
int particle_E(specie_t *s)
{
  int i;
  block_t *b;
  for (i = 0; i < (*s).nblocks; i++)
    {
      b = &(*s).blocks[i];
      if (nanos6_in_final() != 0LU)
        {
          block_particle_E(s, b);
        }
      else
        {
          unsigned long int task_flags_4;
          struct nanos_task_args_particle_E_4 *nanos_data_env_4;
          void *nanos_task_ptr_4;
          task_flags_4 = (((0 != 0) << 0U | (0 != 0) << 1U) | (0 != 0) << 2U) | (0 != 0) << 3U;
          nanos6_create_task(&task_info_var_particle_E_4, &task_invocation_info_particle_E_4, 16LU + 0, (void **)&nanos_data_env_4, &nanos_task_ptr_4, task_flags_4);
          (*nanos_data_env_4).b = b;
          (*nanos_data_env_4).s = s;
          nanos6_submit_task(nanos_task_ptr_4);
        }
    }
  return 0;
}
static void nanos6_unpacked_deps_particle_E_4(block_t **const b, specie_t **const s, void *handler);
static void nanos6_ol_deps_particle_E_4(struct nanos_task_args_particle_E_4 *const arg, void *handler)
{
  nanos6_unpacked_deps_particle_E_4(&((*arg).b), &((*arg).s), handler);
}
static void nanos6_unpacked_deps_particle_E_4(block_t **const b, specie_t **const s, void *handler)
{
  nanos6_register_region_readwrite_depinfo1(handler,  -1, "*b", (void *)(*b), sizeof(block_t), 0, sizeof(block_t));
  ;
}
static void nanos6_unpacked_task_region_particle_E_4(block_t **const b, specie_t **const s, void *device_env, nanos6_address_translation_entry_t *address_translation_table);
static void nanos6_ol_task_region_particle_E_4(struct nanos_task_args_particle_E_4 *const arg, void *device_env, nanos6_address_translation_entry_t *address_translation_table)
{
  nanos6_unpacked_task_region_particle_E_4(&((*arg).b), &((*arg).s), device_env, address_translation_table);
}
static nanos6_task_invocation_info_t task_invocation_info_particle_x_5 = {.invocation_source = "cpic.c:115:9"};
struct  nanos_task_args_particle_x_5
{
  block_t *b;
  specie_t *s;
};
static void nanos6_ol_task_region_particle_x_5(struct nanos_task_args_particle_x_5 *const arg, void *device_env, nanos6_address_translation_entry_t *address_translation_table);
static nanos6_task_implementation_info_t implementations_var_particle_x_5[1L] = {{.device_type_id = nanos6_host_device, .run = (void (*)(void *, void *, nanos6_address_translation_entry_t *))nanos6_ol_task_region_particle_x_5, .get_constraints = 0, .task_label = "block_particle_x", .declaration_source = "cpic.c:115:9", .run_wrapper = 0}};
static void nanos6_ol_deps_particle_x_5(struct nanos_task_args_particle_x_5 *const arg, void *handler);
static nanos6_task_info_t task_info_var_particle_x_5 = {.num_symbols =  -1, .register_depinfo = (void (*)(void *, void *))nanos6_ol_deps_particle_x_5, .get_priority = 0, .type_identifier = 0, .implementation_count = 1, .implementations = implementations_var_particle_x_5, .destroy_args_block = 0, .duplicate_args_block = 0, .reduction_initializers = 0, .reduction_combiners = 0};
static nanos6_task_invocation_info_t task_invocation_info_particle_x_6 = {.invocation_source = "cpic.c:132:9"};
struct  nanos_task_args_particle_x_6
{
  block_t *b;
  block_t *lb;
  block_t *rb;
  specie_t *s;
};
static void nanos6_ol_task_region_particle_x_6(struct nanos_task_args_particle_x_6 *const arg, void *device_env, nanos6_address_translation_entry_t *address_translation_table);
static nanos6_task_implementation_info_t implementations_var_particle_x_6[1L] = {{.device_type_id = nanos6_host_device, .run = (void (*)(void *, void *, nanos6_address_translation_entry_t *))nanos6_ol_task_region_particle_x_6, .get_constraints = 0, .task_label = "block_comm_particle_x", .declaration_source = "cpic.c:132:9", .run_wrapper = 0}};
static void nanos6_ol_deps_particle_x_6(struct nanos_task_args_particle_x_6 *const arg, void *handler);
static nanos6_task_info_t task_info_var_particle_x_6 = {.num_symbols =  -1, .register_depinfo = (void (*)(void *, void *))nanos6_ol_deps_particle_x_6, .get_priority = 0, .type_identifier = 0, .implementation_count = 1, .implementations = implementations_var_particle_x_6, .destroy_args_block = 0, .duplicate_args_block = 0, .reduction_initializers = 0, .reduction_combiners = 0};
int block_particle_x(specie_t *s, block_t *b);
int block_comm_particles(specie_t *s, block_t *left, block_t *b, block_t *right);
int particle_x(specie_t *s)
{
  int i;
  block_t *b;
  int li;
  int ri;
  block_t *lb;
  block_t *rb;
  for (i = 0; i < (*s).nblocks; i++)
    {
      b = &(*s).blocks[i];
      if (nanos6_in_final() != 0LU)
        {
          block_particle_x(s, b);
        }
      else
        {
          unsigned long int task_flags_5;
          struct nanos_task_args_particle_x_5 *nanos_data_env_5;
          void *nanos_task_ptr_5;
          task_flags_5 = (((0 != 0) << 0U | (0 != 0) << 1U) | (0 != 0) << 2U) | (0 != 0) << 3U;
          nanos6_create_task(&task_info_var_particle_x_5, &task_invocation_info_particle_x_5, 16LU + 0, (void **)&nanos_data_env_5, &nanos_task_ptr_5, task_flags_5);
          (*nanos_data_env_5).b = b;
          (*nanos_data_env_5).s = s;
          nanos6_submit_task(nanos_task_ptr_5);
        }
    }
  for (i = 0; i < (*s).nblocks; i++)
    {
      li = ((*s).nblocks + i - 1) % (*s).nblocks;
      ri = (i + 1) % (*s).nblocks;
      b = &(*s).blocks[i];
      lb = &(*s).blocks[li];
      rb = &(*s).blocks[ri];
      if (nanos6_in_final() != 0LU)
        {
          block_comm_particles(s, lb, b, rb);
        }
      else
        {
          unsigned long int task_flags_6;
          struct nanos_task_args_particle_x_6 *nanos_data_env_6;
          void *nanos_task_ptr_6;
          task_flags_6 = (((0 != 0) << 0U | (0 != 0) << 1U) | (0 != 0) << 2U) | (0 != 0) << 3U;
          nanos6_create_task(&task_info_var_particle_x_6, &task_invocation_info_particle_x_6, 32LU + 0, (void **)&nanos_data_env_6, &nanos_task_ptr_6, task_flags_6);
          (*nanos_data_env_6).b = b;
          (*nanos_data_env_6).lb = lb;
          (*nanos_data_env_6).rb = rb;
          (*nanos_data_env_6).s = s;
          nanos6_submit_task(nanos_task_ptr_6);
        }
    }
  return 0;
}
static void nanos6_unpacked_deps_particle_x_6(block_t **const b, block_t **const lb, block_t **const rb, specie_t **const s, void *handler);
static void nanos6_ol_deps_particle_x_6(struct nanos_task_args_particle_x_6 *const arg, void *handler)
{
  nanos6_unpacked_deps_particle_x_6(&((*arg).b), &((*arg).lb), &((*arg).rb), &((*arg).s), handler);
}
static void nanos6_unpacked_deps_particle_x_6(block_t **const b, block_t **const lb, block_t **const rb, specie_t **const s, void *handler)
{
  nanos6_register_region_readwrite_depinfo1(handler,  -1, "*b", (void *)(*b), sizeof(block_t), 0, sizeof(block_t));
  nanos6_register_region_readwrite_depinfo1(handler,  -1, "*lb", (void *)(*lb), sizeof(block_t), 0, sizeof(block_t));
  nanos6_register_region_readwrite_depinfo1(handler,  -1, "*rb", (void *)(*rb), sizeof(block_t), 0, sizeof(block_t));
  ;
}
static void nanos6_unpacked_task_region_particle_x_6(block_t **const b, block_t **const lb, block_t **const rb, specie_t **const s, void *device_env, nanos6_address_translation_entry_t *address_translation_table);
static void nanos6_ol_task_region_particle_x_6(struct nanos_task_args_particle_x_6 *const arg, void *device_env, nanos6_address_translation_entry_t *address_translation_table)
{
  nanos6_unpacked_task_region_particle_x_6(&((*arg).b), &((*arg).lb), &((*arg).rb), &((*arg).s), device_env, address_translation_table);
}
static void nanos6_unpacked_deps_particle_x_5(block_t **const b, specie_t **const s, void *handler);
static void nanos6_ol_deps_particle_x_5(struct nanos_task_args_particle_x_5 *const arg, void *handler)
{
  nanos6_unpacked_deps_particle_x_5(&((*arg).b), &((*arg).s), handler);
}
static void nanos6_unpacked_deps_particle_x_5(block_t **const b, specie_t **const s, void *handler)
{
  nanos6_register_region_readwrite_depinfo1(handler,  -1, "*b", (void *)(*b), sizeof(block_t), 0, sizeof(block_t));
  ;
}
static void nanos6_unpacked_task_region_particle_x_5(block_t **const b, specie_t **const s, void *device_env, nanos6_address_translation_entry_t *address_translation_table);
static void nanos6_ol_task_region_particle_x_5(struct nanos_task_args_particle_x_5 *const arg, void *device_env, nanos6_address_translation_entry_t *address_translation_table)
{
  nanos6_unpacked_task_region_particle_x_5(&((*arg).b), &((*arg).s), device_env, address_translation_table);
}
static nanos6_task_invocation_info_t task_invocation_info_particle_J_7 = {.invocation_source = "cpic.c:151:9"};
struct  nanos_task_args_particle_J_7
{
  block_t *b;
  specie_t *s;
};
static void nanos6_ol_task_region_particle_J_7(struct nanos_task_args_particle_J_7 *const arg, void *device_env, nanos6_address_translation_entry_t *address_translation_table);
static nanos6_task_implementation_info_t implementations_var_particle_J_7[1L] = {{.device_type_id = nanos6_host_device, .run = (void (*)(void *, void *, nanos6_address_translation_entry_t *))nanos6_ol_task_region_particle_J_7, .get_constraints = 0, .task_label = 0, .declaration_source = "cpic.c:151:9", .run_wrapper = 0}};
static void nanos6_ol_deps_particle_J_7(struct nanos_task_args_particle_J_7 *const arg, void *handler);
static nanos6_task_info_t task_info_var_particle_J_7 = {.num_symbols =  -1, .register_depinfo = (void (*)(void *, void *))nanos6_ol_deps_particle_J_7, .get_priority = 0, .type_identifier = 0, .implementation_count = 1, .implementations = implementations_var_particle_J_7, .destroy_args_block = 0, .duplicate_args_block = 0, .reduction_initializers = 0, .reduction_combiners = 0};
int block_particle_J(specie_t *s, block_t *b);
int particle_J(specie_t *s)
{
  int i;
  block_t *b;
  for (i = 0; i < (*s).nblocks; i++)
    {
      b = &(*s).blocks[i];
      if (nanos6_in_final() != 0LU)
        {
          block_particle_J(s, b);
        }
      else
        {
          unsigned long int task_flags_7;
          struct nanos_task_args_particle_J_7 *nanos_data_env_7;
          void *nanos_task_ptr_7;
          task_flags_7 = (((0 != 0) << 0U | (0 != 0) << 1U) | (0 != 0) << 2U) | (0 != 0) << 3U;
          nanos6_create_task(&task_info_var_particle_J_7, &task_invocation_info_particle_J_7, 16LU + 0, (void **)&nanos_data_env_7, &nanos_task_ptr_7, task_flags_7);
          (*nanos_data_env_7).b = b;
          (*nanos_data_env_7).s = s;
          nanos6_submit_task(nanos_task_ptr_7);
        }
    }
  return 0;
}
static void nanos6_unpacked_deps_particle_J_7(block_t **const b, specie_t **const s, void *handler);
static void nanos6_ol_deps_particle_J_7(struct nanos_task_args_particle_J_7 *const arg, void *handler)
{
  nanos6_unpacked_deps_particle_J_7(&((*arg).b), &((*arg).s), handler);
}
static void nanos6_unpacked_deps_particle_J_7(block_t **const b, specie_t **const s, void *handler)
{
  nanos6_register_region_readwrite_depinfo1(handler,  -1, "*b", (void *)(*b), sizeof(block_t), 0, sizeof(block_t));
  ;
}
static void nanos6_unpacked_task_region_particle_J_7(block_t **const b, specie_t **const s, void *device_env, nanos6_address_translation_entry_t *address_translation_table);
static void nanos6_ol_task_region_particle_J_7(struct nanos_task_args_particle_J_7 *const arg, void *device_env, nanos6_address_translation_entry_t *address_translation_table)
{
  nanos6_unpacked_task_region_particle_J_7(&((*arg).b), &((*arg).s), device_env, address_translation_table);
}
void nanos6_taskwait(const char *invocation_source);
int block_print_particles(specie_t *s, block_t *b);
int print_particles(specie_t *s)
{
  int i;
  block_t *b;
  nanos6_taskwait("cpic.c:163:9");
  for (i = 0; i < (*s).nblocks; i++)
    {
      b = &(*s).blocks[i];
      block_print_particles(s, b);
    }
  return 0;
}
struct specie *specie_init();
void specie_step(struct specie *s);
int main()
{
  specie_t *s;
  int i;
  int max_it = 10;
  s = specie_init();
  particle_J(s);
  field_J(s);
  for (i = 0; i < max_it; i++)
    {
      ;
      field_E(s);
      particle_E(s);
      particle_x(s);
      particle_J(s);
      field_J(s);
      specie_step(s);
    }
  nanos6_taskwait("cpic.c:229:9");
}
static void nanos6_unpacked_task_region_field_E_0(block_t **const b, specie_t **const s, void *device_env, nanos6_address_translation_entry_t *address_translation_table)
{
  {
    usleep(100);
    block_field_E((*s), (*b));
  }
}
static void nanos6_unpacked_task_region_field_E_1(block_t **const rb, block_t **const b, void *device_env, nanos6_address_translation_entry_t *address_translation_table)
{
  {
    usleep(100);
    block_comm_field_E((*b), (*rb));
  }
}
static void nanos6_unpacked_task_region_field_J_2(block_t **const b, specie_t **const s, void *device_env, nanos6_address_translation_entry_t *address_translation_table)
{
  block_field_J((*s), (*b));
}
static void nanos6_unpacked_task_region_field_J_3(block_t **const b, block_t **const lb, void *device_env, nanos6_address_translation_entry_t *address_translation_table)
{
  block_comm_field_J((*b), (*lb));
}
static void nanos6_unpacked_task_region_particle_E_4(block_t **const b, specie_t **const s, void *device_env, nanos6_address_translation_entry_t *address_translation_table)
{
  block_particle_E((*s), (*b));
}
static void nanos6_unpacked_task_region_particle_x_5(block_t **const b, specie_t **const s, void *device_env, nanos6_address_translation_entry_t *address_translation_table)
{
  block_particle_x((*s), (*b));
}
static void nanos6_unpacked_task_region_particle_x_6(block_t **const b, block_t **const lb, block_t **const rb, specie_t **const s, void *device_env, nanos6_address_translation_entry_t *address_translation_table)
{
  block_comm_particles((*s), (*lb), (*b), (*rb));
}
static void nanos6_unpacked_task_region_particle_J_7(block_t **const b, specie_t **const s, void *device_env, nanos6_address_translation_entry_t *address_translation_table)
{
  block_particle_J((*s), (*b));
}
