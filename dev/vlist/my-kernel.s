	.file	"kernel.c"
# GNU C17 (SUSE Linux) version 8.1.0 (x86_64-suse-linux)
#	compiled by GNU C version 8.1.0, GMP version 6.1.0, MPFR version 3.1.4, MPC version 1.0.3, isl version isl-0.18-GMP

# GGC heuristics: --param ggc-min-expand=100 --param ggc-min-heapsize=131072
# options passed: 
# -iprefix /gpfs/apps/MN4/GCC/8.1.0/bin/../lib/gcc/x86_64-suse-linux/8.1.0/
# kernel.c -march=skylake-avx512 -mtune=skylake-avx512
# -mprefer-vector-width=512 -auxbase-strip gcc-kernel.s -Ofast -Wall
# -fopenmp-simd -fverbose-asm -fdump-tree-optimized -ftree-vectorize
# options enabled:  -faggressive-loop-optimizations -falign-labels
# -fassociative-math -fasynchronous-unwind-tables -fauto-inc-dec
# -fbranch-count-reg -fcaller-saves -fchkp-check-incomplete-type
# -fchkp-check-read -fchkp-check-write -fchkp-instrument-calls
# -fchkp-narrow-bounds -fchkp-optimize -fchkp-store-bounds
# -fchkp-use-static-bounds -fchkp-use-static-const-bounds
# -fchkp-use-wrappers -fcode-hoisting -fcombine-stack-adjustments -fcommon
# -fcompare-elim -fcprop-registers -fcrossjumping -fcse-follow-jumps
# -fcx-limited-range -fdefer-pop -fdelete-null-pointer-checks
# -fdevirtualize -fdevirtualize-speculatively -fdwarf2-cfi-asm
# -fearly-inlining -feliminate-unused-debug-types -fexpensive-optimizations
# -ffinite-math-only -fforward-propagate -ffp-int-builtin-inexact
# -ffunction-cse -fgcse -fgcse-after-reload -fgcse-lm -fgnu-runtime
# -fgnu-unique -fguess-branch-probability -fhoist-adjacent-loads -fident
# -fif-conversion -fif-conversion2 -findirect-inlining -finline
# -finline-atomics -finline-functions -finline-functions-called-once
# -finline-small-functions -fipa-bit-cp -fipa-cp -fipa-cp-clone -fipa-icf
# -fipa-icf-functions -fipa-icf-variables -fipa-profile -fipa-pure-const
# -fipa-ra -fipa-reference -fipa-sra -fipa-vrp -fira-hoist-pressure
# -fira-share-save-slots -fira-share-spill-slots
# -fisolate-erroneous-paths-dereference -fivopts -fkeep-static-consts
# -fleading-underscore -flifetime-dse -floop-interchange
# -floop-unroll-and-jam -flra-remat -flto-odr-type-merging
# -fmerge-constants -fmerge-debug-strings -fmove-loop-invariants
# -fomit-frame-pointer -foptimize-sibling-calls -foptimize-strlen
# -fpartial-inlining -fpeel-loops -fpeephole -fpeephole2 -fplt
# -fpredictive-commoning -fprefetch-loop-arrays -freciprocal-math -free
# -freg-struct-return -freorder-blocks -freorder-blocks-and-partition
# -freorder-functions -frerun-cse-after-loop
# -fsched-critical-path-heuristic -fsched-dep-count-heuristic
# -fsched-group-heuristic -fsched-interblock -fsched-last-insn-heuristic
# -fsched-rank-heuristic -fsched-spec -fsched-spec-insn-heuristic
# -fsched-stalled-insns-dep -fschedule-fusion -fsemantic-interposition
# -fshow-column -fshrink-wrap -fshrink-wrap-separate
# -fsplit-ivs-in-unroller -fsplit-loops -fsplit-paths -fsplit-wide-types
# -fssa-backprop -fssa-phiopt -fstdarg-opt -fstore-merging
# -fstrict-aliasing -fstrict-volatile-bitfields -fsync-libcalls
# -fthread-jumps -ftoplevel-reorder -ftree-bit-ccp -ftree-builtin-call-dce
# -ftree-ccp -ftree-ch -ftree-coalesce-vars -ftree-copy-prop -ftree-cselim
# -ftree-dce -ftree-dominator-opts -ftree-dse -ftree-forwprop -ftree-fre
# -ftree-loop-distribute-patterns -ftree-loop-distribution
# -ftree-loop-if-convert -ftree-loop-im -ftree-loop-ivcanon
# -ftree-loop-optimize -ftree-loop-vectorize -ftree-parallelize-loops=
# -ftree-partial-pre -ftree-phiprop -ftree-pre -ftree-pta -ftree-reassoc
# -ftree-scev-cprop -ftree-sink -ftree-slp-vectorize -ftree-slsr -ftree-sra
# -ftree-switch-conversion -ftree-tail-merge -ftree-ter -ftree-vrp
# -funit-at-a-time -funsafe-math-optimizations -funswitch-loops
# -funwind-tables -fverbose-asm -fzero-initialized-in-bss
# -m128bit-long-double -m64 -m80387 -madx -maes -malign-stringops -mavx
# -mavx2 -mavx512bw -mavx512cd -mavx512dq -mavx512f -mavx512vl -mbmi -mbmi2
# -mclflushopt -mclwb -mcx16 -mf16c -mfancy-math-387 -mfma -mfp-ret-in-387
# -mfsgsbase -mfxsr -mglibc -mhle -mlong-double-80 -mlzcnt -mmmx -mmovbe
# -mpclmul -mpku -mpopcnt -mprfchw -mpush-args -mrdrnd -mrdseed -mred-zone
# -msahf -msgx -msse -msse2 -msse3 -msse4 -msse4.1 -msse4.2 -mssse3 -mstv
# -mtls-direct-seg-refs -mvzeroupper -mxsave -mxsavec -mxsaveopt -mxsaves

	.text
	.p2align 4,,15
	.globl	i_boris_rotation
	.type	i_boris_rotation, @function
i_boris_rotation:
.LFB5468:
	.cfi_startproc
# kernel.c:40: 	pB = p->B[X];
	movq	80(%rsi), %rcx	# *p_36(D), pB
# kernel.c:41: 	pE = p->E[X];
	movq	56(%rsi), %r8	# *p_36(D).E, pE
# kernel.c:42: 	pu = p->u[X];
	movq	32(%rsi), %rdx	# *p_36(D).u, pu
# /gpfs/apps/MN4/GCC/8.1.0/lib/gcc/x86_64-suse-linux/8.1.0/include/avx512fintrin.h:192:   return (__m512d) __builtin_ia32_broadcastsd512 (__extension__
	vbroadcastsd	%xmm0, %zmm0	# dtqm2, tmp174
	vbroadcastsd	.LC0(%rip), %zmm4	#, tmp179
	vbroadcastsd	.LC1(%rip), %zmm6	#, tmp183
# kernel.c:51: 		B = _mm512_load_pd(&pB[ip]);
	movslq	%edi, %rdi	# ip, ip
	leaq	0(,%rdi,8), %rax	#, _2
# /gpfs/apps/MN4/GCC/8.1.0/lib/gcc/x86_64-suse-linux/8.1.0/include/avx512fintrin.h:329:   return *(__m512d *) __P;
	vmovapd	(%r8,%rdi,8), %zmm11	# MEM[(__m512d * {ref-all})_174], _175
# kernel.c:53: 		u = _mm512_load_pd(&pu[ip]);
	leaq	(%rdx,%rax), %r10	#, _176
# kernel.c:55: 		t[d] = B * dtqm2v;
	vmulpd	(%rcx,%rdi,8), %zmm0, %zmm5	# MEM[(__m512d * {ref-all})_172], tmp174, _178
# kernel.c:59: 		v_minus[d] = u + dtqm2v * E;
	vmovapd	%zmm0, %zmm12	# tmp174, _183
	vfmadd213pd	(%r10), %zmm11, %zmm12	# MEM[(__m512d * {ref-all})_176], _175, _183
# kernel.c:61: 		s[d] = k * t[d] / s_denom[d];
	vmulpd	%zmm5, %zmm4, %zmm2	# _178, tmp179, tmp188
# kernel.c:56: 		s_denom[d] += t[d] * t[d];
	vmovapd	%zmm5, %zmm1	# _178, tmp189
	vfmadd132pd	%zmm5, %zmm6, %zmm1	# _178, tmp183, tmp189
# kernel.c:61: 		s[d] = k * t[d] / s_denom[d];
	vdivpd	%zmm1, %zmm2, %zmm9	# tmp189, tmp188, _186
# kernel.c:52: 		E = _mm512_load_pd(&pE[ip]);
	leaq	8388608(%r8,%rax), %r9	#, _201
# kernel.c:55: 		t[d] = B * dtqm2v;
	vmulpd	8388608(%rcx,%rax), %zmm0, %zmm1	# MEM[(__m512d * {ref-all})_199], tmp174, _205
# kernel.c:59: 		v_minus[d] = u + dtqm2v * E;
	vmovapd	(%r9), %zmm10	# MEM[(__m512d * {ref-all})_201], _210
	vfmadd213pd	8388608(%rdx,%rax), %zmm0, %zmm10	# MEM[(__m512d * {ref-all})_203], tmp174, _210
# kernel.c:61: 		s[d] = k * t[d] / s_denom[d];
	vmulpd	%zmm1, %zmm4, %zmm8	# _205, tmp179, tmp191
# kernel.c:56: 		s_denom[d] += t[d] * t[d];
	vmovapd	%zmm1, %zmm2	# _205, tmp192
	vfmadd132pd	%zmm1, %zmm6, %zmm2	# _205, tmp183, tmp192
# kernel.c:61: 		s[d] = k * t[d] / s_denom[d];
	vdivpd	%zmm2, %zmm8, %zmm8	# tmp192, tmp191, _213
# kernel.c:52: 		E = _mm512_load_pd(&pE[ip]);
	leaq	16777216(%r8,%rax), %rdi	#, _4
# kernel.c:55: 		t[d] = B * dtqm2v;
	vmulpd	16777216(%rcx,%rax), %zmm0, %zmm2	# MEM[(__m512d * {ref-all})_3], tmp174, _6
# kernel.c:59: 		v_minus[d] = u + dtqm2v * E;
	vmovapd	(%rdi), %zmm7	# MEM[(__m512d * {ref-all})_4], _10
	vfmadd213pd	16777216(%rdx,%rax), %zmm0, %zmm7	# MEM[(__m512d * {ref-all})_5], tmp174, _10
# kernel.c:61: 		s[d] = k * t[d] / s_denom[d];
	vmulpd	%zmm4, %zmm2, %zmm4	# tmp179, _6, tmp194
# kernel.c:56: 		s_denom[d] += t[d] * t[d];
	vfmadd231pd	%zmm2, %zmm2, %zmm6	# _6, _6, tmp195
# kernel.c:61: 		s[d] = k * t[d] / s_denom[d];
	vdivpd	%zmm6, %zmm4, %zmm4	# tmp195, tmp194, _12
# kernel.c:17: 	r[X] = a[Y]*b[Z] - a[Z]*b[Y];
	vmulpd	%zmm1, %zmm7, %zmm3	# _205, _10, tmp196
# kernel.c:17: 	r[X] = a[Y]*b[Z] - a[Z]*b[Y];
	vfmsub231pd	%zmm10, %zmm2, %zmm3	# _210, _6, tmp197
# kernel.c:71: 		v_prime[d] += v_minus[d];
	vaddpd	%zmm12, %zmm3, %zmm3	# _183, tmp197, _16
# kernel.c:18: 	r[Y] = a[Z]*b[X] - a[X]*b[Z];
	vmulpd	%zmm12, %zmm2, %zmm2	# _183, _6, tmp198
# kernel.c:18: 	r[Y] = a[Z]*b[X] - a[X]*b[Z];
	vfmsub231pd	%zmm5, %zmm7, %zmm2	# _178, _10, tmp199
# kernel.c:71: 		v_prime[d] += v_minus[d];
	vaddpd	%zmm10, %zmm2, %zmm2	# _210, tmp199, _157
# kernel.c:19: 	r[Z] = a[X]*b[Y] - a[Y]*b[X];
	vmulpd	%zmm10, %zmm5, %zmm5	# _210, _178, tmp200
# kernel.c:19: 	r[Z] = a[X]*b[Y] - a[Y]*b[X];
	vfmsub132pd	%zmm12, %zmm5, %zmm1	# _183, tmp200, tmp201
# kernel.c:71: 		v_prime[d] += v_minus[d];
	vaddpd	%zmm7, %zmm1, %zmm1	# _10, tmp201, _15
# kernel.c:17: 	r[X] = a[Y]*b[Z] - a[Z]*b[Y];
	vmulpd	%zmm8, %zmm1, %zmm5	# _213, _15, tmp202
# kernel.c:17: 	r[X] = a[Y]*b[Z] - a[Z]*b[Y];
	vfmsub231pd	%zmm2, %zmm4, %zmm5	# _157, _12, tmp203
# kernel.c:82: 		v_plus[d] += v_minus[d];
	vaddpd	%zmm12, %zmm5, %zmm5	# _183, tmp203, tmp204
# kernel.c:85: 		u = v_plus[d] + dtqm2v * E;
	vfmadd231pd	%zmm11, %zmm0, %zmm5	# _175, tmp174, u
# /gpfs/apps/MN4/GCC/8.1.0/lib/gcc/x86_64-suse-linux/8.1.0/include/avx512fintrin.h:8499:   __builtin_ia32_movntpd512 (__P, (__v8df) __A);
	vmovntpd	%zmm5, (%r10)	# u,* _176
# kernel.c:88: 		_mm512_stream_pd(&p->u[d][ip], u);
	movq	40(%rsi), %rdx	# MEM[(struct particle_header *)p_36(D) + 40B], tmp207
	addq	%rax, %rdx	# _2, tmp207
# kernel.c:18: 	r[Y] = a[Z]*b[X] - a[X]*b[Z];
	vmulpd	%zmm3, %zmm4, %zmm4	# _16, _12, tmp209
# kernel.c:18: 	r[Y] = a[Z]*b[X] - a[X]*b[Z];
	vfmsub132pd	%zmm9, %zmm4, %zmm1	# _186, tmp209, tmp210
# kernel.c:82: 		v_plus[d] += v_minus[d];
	vaddpd	%zmm10, %zmm1, %zmm1	# _210, tmp210, tmp211
# kernel.c:85: 		u = v_plus[d] + dtqm2v * E;
	vfmadd231pd	(%r9), %zmm0, %zmm1	# MEM[(__m512d * {ref-all})_201], tmp174, u
# /gpfs/apps/MN4/GCC/8.1.0/lib/gcc/x86_64-suse-linux/8.1.0/include/avx512fintrin.h:8499:   __builtin_ia32_movntpd512 (__P, (__v8df) __A);
	vmovntpd	%zmm1, (%rdx)	# u,
# kernel.c:88: 		_mm512_stream_pd(&p->u[d][ip], u);
	addq	48(%rsi), %rax	# MEM[(struct particle_header *)p_36(D) + 48B], tmp214
# kernel.c:19: 	r[Z] = a[X]*b[Y] - a[Y]*b[X];
	vmulpd	%zmm9, %zmm2, %zmm2	# _186, _157, tmp216
# kernel.c:19: 	r[Z] = a[X]*b[Y] - a[Y]*b[X];
	vfmsub231pd	%zmm8, %zmm3, %zmm2	# _213, _16, tmp217
# kernel.c:82: 		v_plus[d] += v_minus[d];
	vaddpd	%zmm7, %zmm2, %zmm2	# _10, tmp217, tmp218
# kernel.c:85: 		u = v_plus[d] + dtqm2v * E;
	vfmadd132pd	(%rdi), %zmm2, %zmm0	# MEM[(__m512d * {ref-all})_4], tmp218, u
# /gpfs/apps/MN4/GCC/8.1.0/lib/gcc/x86_64-suse-linux/8.1.0/include/avx512fintrin.h:8499:   __builtin_ia32_movntpd512 (__P, (__v8df) __A);
	vmovntpd	%zmm0, (%rax)	# u,
	vzeroupper
# kernel.c:93: }
	ret	
	.cfi_endproc
.LFE5468:
	.size	i_boris_rotation, .-i_boris_rotation
	.section	.rodata.cst16,"aM",@progbits,16
	.align 16
.LC0:
	.long	0
	.long	1073741824
	.long	0
	.long	0
	.align 16
.LC1:
	.long	0
	.long	1072693248
	.long	0
	.long	0
	.ident	"GCC: (SUSE Linux) 8.1.0"
	.section	.note.GNU-stack,"",@progbits
