	.text
	.file	"main.c"
	.section	.rodata.cst8,"aM",@progbits,8
	.p2align	3               # -- Begin function main
.LCPI0_0:
	.quad	4617474832798419256     # double 5.1415000000000006
.LCPI0_1:
	.quad	4614256447914709615     # double 3.1415000000000002
	.text
	.globl	main
	.p2align	4, 0x90
	.type	main,@function
main:                                   # @main
	.cfi_startproc
# %bb.0:
	pushq	%rax
	.cfi_def_cfa_offset 16
	callq	rand@PLT
	cmpl	$3, %eax
	je	.LBB0_1
# %bb.2:
	vbroadcastsd	.LCPI0_1(%rip), %ymm0 # ymm0 = [3.1415000000000002E+0,3.1415000000000002E+0,3.1415000000000002E+0,3.1415000000000002E+0]
	vextractf128	$1, %ymm0, %xmm0
	vcvttsd2si	%xmm0, %eax
	popq	%rcx
	.cfi_def_cfa_offset 8
	vzeroupper
	retq
.LBB0_1:
	.cfi_def_cfa_offset 16
	vbroadcastsd	.LCPI0_0(%rip), %ymm0 # ymm0 = [5.1415000000000006E+0,5.1415000000000006E+0,5.1415000000000006E+0,5.1415000000000006E+0]
	vextractf128	$1, %ymm0, %xmm0
	vcvttsd2si	%xmm0, %eax
	popq	%rcx
	.cfi_def_cfa_offset 8
	vzeroupper
	retq
.Lfunc_end0:
	.size	main, .Lfunc_end0-main
	.cfi_endproc
                                        # -- End function

	.ident	"clang version 9.0.1 "
	.section	".note.GNU-stack","",@progbits
	.addrsig
