	.text
	.file	"main.c"
	.globl	main                    # -- Begin function main
	.p2align	4, 0x90
	.type	main,@function
main:                                   # @main
	.cfi_startproc
# %bb.0:
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	subq	$16, %rsp
	movl	$0, -4(%rbp)
	movl	$15, %edi
	callq	vmsk_set
	callq	vmsk_get
	cmpq	$15, %rax
	jne	.LBB0_2
# %bb.1:
	jmp	.LBB0_3
.LBB0_2:
	leaq	.L.str(%rip), %rdi
	leaq	.L.str.1(%rip), %rsi
	movl	$9, %edx
	leaq	.L__PRETTY_FUNCTION__.main(%rip), %rcx
	vzeroupper
	callq	__assert_fail@PLT
.LBB0_3:
	movl	$7, %edi
	callq	vmsk_set
	callq	vmsk_get
	cmpq	$7, %rax
	jne	.LBB0_5
# %bb.4:
	jmp	.LBB0_6
.LBB0_5:
	leaq	.L.str.2(%rip), %rdi
	leaq	.L.str.1(%rip), %rsi
	movl	$10, %edx
	leaq	.L__PRETTY_FUNCTION__.main(%rip), %rcx
	vzeroupper
	callq	__assert_fail@PLT
.LBB0_6:
	movl	$3, %edi
	callq	vmsk_set
	callq	vmsk_get
	cmpq	$3, %rax
	jne	.LBB0_8
# %bb.7:
	jmp	.LBB0_9
.LBB0_8:
	leaq	.L.str.3(%rip), %rdi
	leaq	.L.str.1(%rip), %rsi
	movl	$11, %edx
	leaq	.L__PRETTY_FUNCTION__.main(%rip), %rcx
	vzeroupper
	callq	__assert_fail@PLT
.LBB0_9:
	movl	$1, %edi
	callq	vmsk_set
	callq	vmsk_get
	cmpq	$1, %rax
	jne	.LBB0_11
# %bb.10:
	jmp	.LBB0_12
.LBB0_11:
	leaq	.L.str.4(%rip), %rdi
	leaq	.L.str.1(%rip), %rsi
	movl	$12, %edx
	leaq	.L__PRETTY_FUNCTION__.main(%rip), %rcx
	vzeroupper
	callq	__assert_fail@PLT
.LBB0_12:
	xorl	%eax, %eax
	movl	%eax, %edi
	callq	vmsk_set
	callq	vmsk_get
	cmpq	$0, %rax
	jne	.LBB0_14
# %bb.13:
	jmp	.LBB0_15
.LBB0_14:
	leaq	.L.str.5(%rip), %rdi
	leaq	.L.str.1(%rip), %rsi
	movl	$13, %edx
	leaq	.L__PRETTY_FUNCTION__.main(%rip), %rcx
	vzeroupper
	callq	__assert_fail@PLT
.LBB0_15:
	xorl	%eax, %eax
	addq	$16, %rsp
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	vzeroupper
	retq
.Lfunc_end0:
	.size	main, .Lfunc_end0-main
	.cfi_endproc
                                        # -- End function
	.p2align	4, 0x90         # -- Begin function vmsk_get
	.type	vmsk_get,@function
vmsk_get:                               # @vmsk_get
	.cfi_startproc
# %bb.0:
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	andq	$-32, %rsp
	subq	$96, %rsp
	vmovdqa	%ymm0, (%rsp)
	vmovdqa	(%rsp), %ymm0
	vmovapd	%ymm0, 32(%rsp)
	vmovapd	32(%rsp), %ymm0
	vmovmskpd	%ymm0, %eax
	cltq
	movq	%rbp, %rsp
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	vzeroupper
	retq
.Lfunc_end1:
	.size	vmsk_get, .Lfunc_end1-vmsk_get
	.cfi_endproc
                                        # -- End function
	.p2align	4, 0x90         # -- Begin function vmsk_set
	.type	vmsk_set,@function
vmsk_set:                               # @vmsk_set
	.cfi_startproc
# %bb.0:
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	andq	$-32, %rsp
	subq	$192, %rsp
	movq	%fs:40, %rax
	movq	%rax, 168(%rsp)
	movq	%rdi, 56(%rsp)
	vmovdqa	.L__const.vmsk_set.t(%rip), %xmm0
	vmovdqa	%xmm0, 144(%rsp)
	movq	56(%rsp), %rax
                                        # kill: def $eax killed $eax killed $rax
	movl	%eax, %ecx
	andl	$8, %ecx
	movl	%ecx, %edx
	movq	144(%rsp,%rdx), %rdx
	movl	%eax, %ecx
	andl	$4, %ecx
	movl	%ecx, %esi
	movq	144(%rsp,%rsi,2), %rsi
	movl	%eax, %ecx
	andl	$2, %ecx
	movl	%ecx, %edi
	movq	144(%rsp,%rdi,4), %rdi
	andl	$1, %eax
	movl	%eax, %r8d
	movq	144(%rsp,%r8,8), %r8
	movq	%rdx, 136(%rsp)
	movq	%rsi, 128(%rsp)
	movq	%rdi, 120(%rsp)
	movq	%r8, 112(%rsp)
	vmovq	120(%rsp), %xmm0        # xmm0 = mem[0],zero
	vmovq	112(%rsp), %xmm1        # xmm1 = mem[0],zero
	vpunpcklqdq	%xmm0, %xmm1, %xmm0 # xmm0 = xmm1[0],xmm0[0]
	vmovq	136(%rsp), %xmm1        # xmm1 = mem[0],zero
	vmovq	128(%rsp), %xmm2        # xmm2 = mem[0],zero
	vpunpcklqdq	%xmm1, %xmm2, %xmm1 # xmm1 = xmm2[0],xmm1[0]
	vmovdqa	%xmm1, 80(%rsp)
	vmovdqa	%xmm0, 64(%rsp)
	vmovaps	64(%rsp), %ymm0
	movq	%fs:40, %rdx
	movq	168(%rsp), %rsi
	cmpq	%rsi, %rdx
	vmovaps	%ymm0, (%rsp)           # 32-byte Spill
	jne	.LBB2_2
# %bb.1:
	vmovaps	(%rsp), %ymm0           # 32-byte Reload
	movq	%rbp, %rsp
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	retq
.LBB2_2:
	.cfi_def_cfa %rbp, 16
	vzeroupper
	callq	__stack_chk_fail@PLT
.Lfunc_end2:
	.size	vmsk_set, .Lfunc_end2-vmsk_set
	.cfi_endproc
                                        # -- End function
	.type	.L.str,@object          # @.str
	.section	.rodata.str1.1,"aMS",@progbits,1
.L.str:
	.asciz	"vmsk_get(vmsk_set(0x0f)) == 0x0f"
	.size	.L.str, 33

	.type	.L.str.1,@object        # @.str.1
.L.str.1:
	.asciz	"main.c"
	.size	.L.str.1, 7

	.type	.L__PRETTY_FUNCTION__.main,@object # @__PRETTY_FUNCTION__.main
.L__PRETTY_FUNCTION__.main:
	.asciz	"int main()"
	.size	.L__PRETTY_FUNCTION__.main, 11

	.type	.L.str.2,@object        # @.str.2
.L.str.2:
	.asciz	"vmsk_get(vmsk_set(0x07)) == 0x07"
	.size	.L.str.2, 33

	.type	.L.str.3,@object        # @.str.3
.L.str.3:
	.asciz	"vmsk_get(vmsk_set(0x03)) == 0x03"
	.size	.L.str.3, 33

	.type	.L.str.4,@object        # @.str.4
.L.str.4:
	.asciz	"vmsk_get(vmsk_set(0x01)) == 0x01"
	.size	.L.str.4, 33

	.type	.L.str.5,@object        # @.str.5
.L.str.5:
	.asciz	"vmsk_get(vmsk_set(0x00)) == 0x00"
	.size	.L.str.5, 33

	.type	.L__const.vmsk_set.t,@object # @__const.vmsk_set.t
	.section	.rodata.cst16,"aM",@progbits,16
	.p2align	4
.L__const.vmsk_set.t:
	.quad	0                       # 0x0
	.quad	-1                      # 0xffffffffffffffff
	.size	.L__const.vmsk_set.t, 16


	.ident	"clang version 9.0.1 "
	.section	".note.GNU-stack","",@progbits
	.addrsig
	.addrsig_sym vmsk_get
	.addrsig_sym vmsk_set
	.addrsig_sym __assert_fail
	.addrsig_sym __stack_chk_fail
