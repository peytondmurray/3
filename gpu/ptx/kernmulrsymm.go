package ptx

const KERNMULRSYMM = `
//
// Generated by NVIDIA NVVM Compiler
// Compiler built on Wed Aug  1 02:51:19 2012 (1343782279)
// Cuda compilation tools, release 5.0, V0.2.1221
//

.version 3.1
.target sm_30
.address_size 64

	.file	1 "/tmp/tmpxft_00000389_00000000-9_kernmulrsymm.cpp3.i"
	.file	2 "/home/arne/src/nimble-cube/gpu/ptx/kernmulrsymm.cu"

.visible .entry kernmulRSymm(
	.param .u64 kernmulRSymm_param_0,
	.param .u64 kernmulRSymm_param_1,
	.param .u64 kernmulRSymm_param_2,
	.param .u64 kernmulRSymm_param_3,
	.param .u64 kernmulRSymm_param_4,
	.param .u64 kernmulRSymm_param_5,
	.param .u64 kernmulRSymm_param_6,
	.param .u64 kernmulRSymm_param_7,
	.param .u64 kernmulRSymm_param_8,
	.param .u32 kernmulRSymm_param_9
)
{
	.reg .pred 	%p<2>;
	.reg .s32 	%r<29>;
	.reg .f32 	%f<31>;
	.reg .s64 	%rd<34>;


	ld.param.u64 	%rd9, [kernmulRSymm_param_0];
	ld.param.u64 	%rd10, [kernmulRSymm_param_1];
	ld.param.u64 	%rd11, [kernmulRSymm_param_2];
	ld.param.u64 	%rd12, [kernmulRSymm_param_3];
	ld.param.u64 	%rd13, [kernmulRSymm_param_4];
	ld.param.u64 	%rd14, [kernmulRSymm_param_5];
	ld.param.u64 	%rd15, [kernmulRSymm_param_6];
	ld.param.u64 	%rd16, [kernmulRSymm_param_7];
	ld.param.u64 	%rd17, [kernmulRSymm_param_8];
	ld.param.u32 	%r3, [kernmulRSymm_param_9];
	cvta.to.global.u64 	%rd1, %rd17;
	cvta.to.global.u64 	%rd2, %rd16;
	cvta.to.global.u64 	%rd3, %rd15;
	cvta.to.global.u64 	%rd4, %rd14;
	cvta.to.global.u64 	%rd5, %rd13;
	cvta.to.global.u64 	%rd6, %rd12;
	cvta.to.global.u64 	%rd7, %rd11;
	cvta.to.global.u64 	%rd8, %rd10;
	.loc 2 12 1
	mov.u32 	%r4, %nctaid.x;
	mov.u32 	%r5, %ctaid.y;
	mov.u32 	%r6, %ctaid.x;
	mad.lo.s32 	%r7, %r4, %r5, %r6;
	mov.u32 	%r8, %ntid.x;
	mov.u32 	%r9, %tid.x;
	mad.lo.s32 	%r1, %r7, %r8, %r9;
	.loc 2 13 1
	shl.b32 	%r2, %r1, 1;
	.loc 2 15 1
	setp.ge.s32 	%p1, %r1, %r3;
	@%p1 bra 	BB0_2;

	cvta.to.global.u64 	%rd18, %rd9;
	.loc 2 16 1
	mul.wide.s32 	%rd19, %r2, 4;
	add.s64 	%rd20, %rd18, %rd19;
	.loc 2 17 1
	add.s32 	%r10, %r2, 1;
	mul.wide.s32 	%rd21, %r10, 4;
	add.s64 	%rd22, %rd18, %rd21;
	ld.global.f32 	%f1, [%rd22];
	.loc 2 19 1
	add.s64 	%rd23, %rd8, %rd19;
	.loc 2 20 1
	add.s64 	%rd24, %rd8, %rd21;
	ld.global.f32 	%f2, [%rd24];
	.loc 2 22 1
	add.s64 	%rd25, %rd7, %rd19;
	.loc 2 23 1
	add.s64 	%rd26, %rd7, %rd21;
	ld.global.f32 	%f3, [%rd26];
	.loc 2 25 1
	mul.wide.s32 	%rd27, %r1, 4;
	add.s64 	%rd28, %rd6, %rd27;
	.loc 2 26 1
	add.s64 	%rd29, %rd5, %rd27;
	ld.global.f32 	%f4, [%rd29];
	.loc 2 27 1
	add.s64 	%rd30, %rd4, %rd27;
	ld.global.f32 	%f5, [%rd30];
	.loc 2 29 1
	add.s64 	%rd31, %rd3, %rd27;
	ld.global.f32 	%f6, [%rd31];
	.loc 2 30 1
	add.s64 	%rd32, %rd2, %rd27;
	.loc 2 31 1
	add.s64 	%rd33, %rd1, %rd27;
	.loc 2 25 1
	ld.global.f32 	%f7, [%rd28];
	.loc 2 16 1
	ld.global.f32 	%f8, [%rd20];
	.loc 2 31 1
	ld.global.f32 	%f9, [%rd33];
	.loc 2 19 1
	ld.global.f32 	%f10, [%rd23];
	.loc 2 33 1
	mul.ftz.f32 	%f11, %f10, %f9;
	fma.rn.ftz.f32 	%f12, %f8, %f7, %f11;
	.loc 2 30 1
	ld.global.f32 	%f13, [%rd32];
	.loc 2 22 1
	ld.global.f32 	%f14, [%rd25];
	.loc 2 33 1
	fma.rn.ftz.f32 	%f15, %f14, %f13, %f12;
	st.global.f32 	[%rd20], %f15;
	.loc 2 34 1
	mul.ftz.f32 	%f16, %f2, %f9;
	fma.rn.ftz.f32 	%f17, %f1, %f7, %f16;
	fma.rn.ftz.f32 	%f18, %f3, %f13, %f17;
	st.global.f32 	[%rd22], %f18;
	.loc 2 36 1
	mul.ftz.f32 	%f19, %f10, %f4;
	fma.rn.ftz.f32 	%f20, %f8, %f9, %f19;
	fma.rn.ftz.f32 	%f21, %f14, %f6, %f20;
	st.global.f32 	[%rd23], %f21;
	.loc 2 37 1
	mul.ftz.f32 	%f22, %f2, %f4;
	fma.rn.ftz.f32 	%f23, %f1, %f9, %f22;
	fma.rn.ftz.f32 	%f24, %f3, %f6, %f23;
	st.global.f32 	[%rd24], %f24;
	.loc 2 39 1
	mul.ftz.f32 	%f25, %f10, %f6;
	fma.rn.ftz.f32 	%f26, %f8, %f13, %f25;
	fma.rn.ftz.f32 	%f27, %f14, %f5, %f26;
	st.global.f32 	[%rd25], %f27;
	.loc 2 40 1
	mul.ftz.f32 	%f28, %f2, %f6;
	fma.rn.ftz.f32 	%f29, %f1, %f13, %f28;
	fma.rn.ftz.f32 	%f30, %f3, %f5, %f29;
	st.global.f32 	[%rd26], %f30;

BB0_2:
	.loc 2 42 2
	ret;
}


`
