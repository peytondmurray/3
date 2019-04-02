package cuda

/*
 THIS FILE IS AUTO-GENERATED BY CUDA2GO.
 EDITING IS FUTILE.
*/

import (
	"github.com/mumax/3/cuda/cu"
	"github.com/mumax/3/timer"
	"sync"
	"unsafe"
)

// CUDA handle for settheta kernel
var settheta_code cu.Function

// Stores the arguments for settheta kernel invocation
type settheta_args_t struct {
	arg_s  unsafe.Pointer
	arg_mx unsafe.Pointer
	arg_my unsafe.Pointer
	arg_mz unsafe.Pointer
	arg_Nx int
	arg_Ny int
	arg_Nz int
	argptr [7]unsafe.Pointer
	sync.Mutex
}

// Stores the arguments for settheta kernel invocation
var settheta_args settheta_args_t

func init() {
	// CUDA driver kernel call wants pointers to arguments, set them up once.
	settheta_args.argptr[0] = unsafe.Pointer(&settheta_args.arg_s)
	settheta_args.argptr[1] = unsafe.Pointer(&settheta_args.arg_mx)
	settheta_args.argptr[2] = unsafe.Pointer(&settheta_args.arg_my)
	settheta_args.argptr[3] = unsafe.Pointer(&settheta_args.arg_mz)
	settheta_args.argptr[4] = unsafe.Pointer(&settheta_args.arg_Nx)
	settheta_args.argptr[5] = unsafe.Pointer(&settheta_args.arg_Ny)
	settheta_args.argptr[6] = unsafe.Pointer(&settheta_args.arg_Nz)
}

// Wrapper for settheta CUDA kernel, asynchronous.
func k_settheta_async(s unsafe.Pointer, mx unsafe.Pointer, my unsafe.Pointer, mz unsafe.Pointer, Nx int, Ny int, Nz int, cfg *config) {
	if Synchronous { // debug
		Sync()
		timer.Start("settheta")
	}

	settheta_args.Lock()
	defer settheta_args.Unlock()

	if settheta_code == 0 {
		settheta_code = fatbinLoad(settheta_map, "settheta")
	}

	settheta_args.arg_s = s
	settheta_args.arg_mx = mx
	settheta_args.arg_my = my
	settheta_args.arg_mz = mz
	settheta_args.arg_Nx = Nx
	settheta_args.arg_Ny = Ny
	settheta_args.arg_Nz = Nz

	args := settheta_args.argptr[:]
	cu.LaunchKernel(settheta_code, cfg.Grid.X, cfg.Grid.Y, cfg.Grid.Z, cfg.Block.X, cfg.Block.Y, cfg.Block.Z, 0, stream0, args)

	if Synchronous { // debug
		Sync()
		timer.Stop("settheta")
	}
}

// maps compute capability on PTX code for settheta kernel.
var settheta_map = map[int]string{0: "",
	30: settheta_ptx_30,
	35: settheta_ptx_35,
	37: settheta_ptx_37,
	50: settheta_ptx_50,
	52: settheta_ptx_52,
	53: settheta_ptx_53,
	60: settheta_ptx_60,
	61: settheta_ptx_61,
	70: settheta_ptx_70,
	75: settheta_ptx_75}

// settheta PTX code for various compute capabilities.
const (
	settheta_ptx_30 = `
.version 6.4
.target sm_30
.address_size 64

	// .globl	settheta

.visible .entry settheta(
	.param .u64 settheta_param_0,
	.param .u64 settheta_param_1,
	.param .u64 settheta_param_2,
	.param .u64 settheta_param_3,
	.param .u32 settheta_param_4,
	.param .u32 settheta_param_5,
	.param .u32 settheta_param_6
)
{
	.reg .pred 	%p<8>;
	.reg .f32 	%f<27>;
	.reg .b32 	%r<18>;
	.reg .b64 	%rd<8>;


	ld.param.u64 	%rd1, [settheta_param_0];
	ld.param.u64 	%rd2, [settheta_param_3];
	ld.param.u32 	%r4, [settheta_param_4];
	ld.param.u32 	%r5, [settheta_param_5];
	ld.param.u32 	%r6, [settheta_param_6];
	mov.u32 	%r7, %ctaid.x;
	mov.u32 	%r8, %ntid.x;
	mov.u32 	%r9, %tid.x;
	mad.lo.s32 	%r1, %r8, %r7, %r9;
	mov.u32 	%r10, %ntid.y;
	mov.u32 	%r11, %ctaid.y;
	mov.u32 	%r12, %tid.y;
	mad.lo.s32 	%r2, %r10, %r11, %r12;
	mov.u32 	%r13, %ntid.z;
	mov.u32 	%r14, %ctaid.z;
	mov.u32 	%r15, %tid.z;
	mad.lo.s32 	%r3, %r13, %r14, %r15;
	setp.ge.s32	%p1, %r2, %r5;
	setp.ge.s32	%p2, %r1, %r4;
	or.pred  	%p3, %p1, %p2;
	setp.ge.s32	%p4, %r3, %r6;
	or.pred  	%p5, %p3, %p4;
	@%p5 bra 	BB0_2;

	cvta.to.global.u64 	%rd3, %rd2;
	mad.lo.s32 	%r16, %r3, %r5, %r2;
	mad.lo.s32 	%r17, %r16, %r4, %r1;
	mul.wide.s32 	%rd4, %r17, 4;
	add.s64 	%rd5, %rd3, %rd4;
	ld.global.f32 	%f1, [%rd5];
	abs.f32 	%f2, %f1;
	mov.f32 	%f3, 0f3F800000;
	sub.f32 	%f4, %f3, %f2;
	mul.f32 	%f5, %f4, 0f3F000000;
	sqrt.rn.f32 	%f6, %f5;
	setp.gt.f32	%p6, %f2, 0f3F11EB85;
	selp.f32	%f7, %f6, %f2, %p6;
	mul.f32 	%f8, %f7, %f7;
	mov.f32 	%f9, 0f3C94D2E9;
	mov.f32 	%f10, 0f3D53F941;
	fma.rn.f32 	%f11, %f10, %f8, %f9;
	mov.f32 	%f12, 0f3D3F841F;
	fma.rn.f32 	%f13, %f11, %f8, %f12;
	mov.f32 	%f14, 0f3D994929;
	fma.rn.f32 	%f15, %f13, %f8, %f14;
	mov.f32 	%f16, 0f3E2AAB94;
	fma.rn.f32 	%f17, %f15, %f8, %f16;
	mul.f32 	%f18, %f8, %f17;
	fma.rn.f32 	%f19, %f18, %f7, %f7;
	add.f32 	%f20, %f19, %f19;
	mov.f32 	%f21, 0f3FC90FDB;
	sub.f32 	%f22, %f21, %f19;
	selp.f32	%f23, %f20, %f22, %p6;
	setp.lt.f32	%p7, %f1, 0f00000000;
	mov.f32 	%f24, 0f40490FDB;
	sub.f32 	%f25, %f24, %f23;
	selp.f32	%f26, %f25, %f23, %p7;
	cvta.to.global.u64 	%rd6, %rd1;
	add.s64 	%rd7, %rd6, %rd4;
	st.global.f32 	[%rd7], %f26;

BB0_2:
	ret;
}


`
	settheta_ptx_35 = `
.version 6.4
.target sm_35
.address_size 64

	// .globl	settheta

.visible .entry settheta(
	.param .u64 settheta_param_0,
	.param .u64 settheta_param_1,
	.param .u64 settheta_param_2,
	.param .u64 settheta_param_3,
	.param .u32 settheta_param_4,
	.param .u32 settheta_param_5,
	.param .u32 settheta_param_6
)
{
	.reg .pred 	%p<8>;
	.reg .f32 	%f<27>;
	.reg .b32 	%r<18>;
	.reg .b64 	%rd<8>;


	ld.param.u64 	%rd1, [settheta_param_0];
	ld.param.u64 	%rd2, [settheta_param_3];
	ld.param.u32 	%r4, [settheta_param_4];
	ld.param.u32 	%r5, [settheta_param_5];
	ld.param.u32 	%r6, [settheta_param_6];
	mov.u32 	%r7, %ctaid.x;
	mov.u32 	%r8, %ntid.x;
	mov.u32 	%r9, %tid.x;
	mad.lo.s32 	%r1, %r8, %r7, %r9;
	mov.u32 	%r10, %ntid.y;
	mov.u32 	%r11, %ctaid.y;
	mov.u32 	%r12, %tid.y;
	mad.lo.s32 	%r2, %r10, %r11, %r12;
	mov.u32 	%r13, %ntid.z;
	mov.u32 	%r14, %ctaid.z;
	mov.u32 	%r15, %tid.z;
	mad.lo.s32 	%r3, %r13, %r14, %r15;
	setp.ge.s32	%p1, %r2, %r5;
	setp.ge.s32	%p2, %r1, %r4;
	or.pred  	%p3, %p1, %p2;
	setp.ge.s32	%p4, %r3, %r6;
	or.pred  	%p5, %p3, %p4;
	@%p5 bra 	BB0_2;

	cvta.to.global.u64 	%rd3, %rd2;
	mad.lo.s32 	%r16, %r3, %r5, %r2;
	mad.lo.s32 	%r17, %r16, %r4, %r1;
	mul.wide.s32 	%rd4, %r17, 4;
	add.s64 	%rd5, %rd3, %rd4;
	ld.global.nc.f32 	%f1, [%rd5];
	abs.f32 	%f2, %f1;
	mov.f32 	%f3, 0f3F800000;
	sub.f32 	%f4, %f3, %f2;
	mul.f32 	%f5, %f4, 0f3F000000;
	sqrt.rn.f32 	%f6, %f5;
	setp.gt.f32	%p6, %f2, 0f3F11EB85;
	selp.f32	%f7, %f6, %f2, %p6;
	mul.f32 	%f8, %f7, %f7;
	mov.f32 	%f9, 0f3C94D2E9;
	mov.f32 	%f10, 0f3D53F941;
	fma.rn.f32 	%f11, %f10, %f8, %f9;
	mov.f32 	%f12, 0f3D3F841F;
	fma.rn.f32 	%f13, %f11, %f8, %f12;
	mov.f32 	%f14, 0f3D994929;
	fma.rn.f32 	%f15, %f13, %f8, %f14;
	mov.f32 	%f16, 0f3E2AAB94;
	fma.rn.f32 	%f17, %f15, %f8, %f16;
	mul.f32 	%f18, %f8, %f17;
	fma.rn.f32 	%f19, %f18, %f7, %f7;
	add.f32 	%f20, %f19, %f19;
	mov.f32 	%f21, 0f3FC90FDB;
	sub.f32 	%f22, %f21, %f19;
	selp.f32	%f23, %f20, %f22, %p6;
	setp.lt.f32	%p7, %f1, 0f00000000;
	mov.f32 	%f24, 0f40490FDB;
	sub.f32 	%f25, %f24, %f23;
	selp.f32	%f26, %f25, %f23, %p7;
	cvta.to.global.u64 	%rd6, %rd1;
	add.s64 	%rd7, %rd6, %rd4;
	st.global.f32 	[%rd7], %f26;

BB0_2:
	ret;
}


`
	settheta_ptx_37 = `
.version 6.4
.target sm_37
.address_size 64

	// .globl	settheta

.visible .entry settheta(
	.param .u64 settheta_param_0,
	.param .u64 settheta_param_1,
	.param .u64 settheta_param_2,
	.param .u64 settheta_param_3,
	.param .u32 settheta_param_4,
	.param .u32 settheta_param_5,
	.param .u32 settheta_param_6
)
{
	.reg .pred 	%p<8>;
	.reg .f32 	%f<27>;
	.reg .b32 	%r<18>;
	.reg .b64 	%rd<8>;


	ld.param.u64 	%rd1, [settheta_param_0];
	ld.param.u64 	%rd2, [settheta_param_3];
	ld.param.u32 	%r4, [settheta_param_4];
	ld.param.u32 	%r5, [settheta_param_5];
	ld.param.u32 	%r6, [settheta_param_6];
	mov.u32 	%r7, %ctaid.x;
	mov.u32 	%r8, %ntid.x;
	mov.u32 	%r9, %tid.x;
	mad.lo.s32 	%r1, %r8, %r7, %r9;
	mov.u32 	%r10, %ntid.y;
	mov.u32 	%r11, %ctaid.y;
	mov.u32 	%r12, %tid.y;
	mad.lo.s32 	%r2, %r10, %r11, %r12;
	mov.u32 	%r13, %ntid.z;
	mov.u32 	%r14, %ctaid.z;
	mov.u32 	%r15, %tid.z;
	mad.lo.s32 	%r3, %r13, %r14, %r15;
	setp.ge.s32	%p1, %r2, %r5;
	setp.ge.s32	%p2, %r1, %r4;
	or.pred  	%p3, %p1, %p2;
	setp.ge.s32	%p4, %r3, %r6;
	or.pred  	%p5, %p3, %p4;
	@%p5 bra 	BB0_2;

	cvta.to.global.u64 	%rd3, %rd2;
	mad.lo.s32 	%r16, %r3, %r5, %r2;
	mad.lo.s32 	%r17, %r16, %r4, %r1;
	mul.wide.s32 	%rd4, %r17, 4;
	add.s64 	%rd5, %rd3, %rd4;
	ld.global.nc.f32 	%f1, [%rd5];
	abs.f32 	%f2, %f1;
	mov.f32 	%f3, 0f3F800000;
	sub.f32 	%f4, %f3, %f2;
	mul.f32 	%f5, %f4, 0f3F000000;
	sqrt.rn.f32 	%f6, %f5;
	setp.gt.f32	%p6, %f2, 0f3F11EB85;
	selp.f32	%f7, %f6, %f2, %p6;
	mul.f32 	%f8, %f7, %f7;
	mov.f32 	%f9, 0f3C94D2E9;
	mov.f32 	%f10, 0f3D53F941;
	fma.rn.f32 	%f11, %f10, %f8, %f9;
	mov.f32 	%f12, 0f3D3F841F;
	fma.rn.f32 	%f13, %f11, %f8, %f12;
	mov.f32 	%f14, 0f3D994929;
	fma.rn.f32 	%f15, %f13, %f8, %f14;
	mov.f32 	%f16, 0f3E2AAB94;
	fma.rn.f32 	%f17, %f15, %f8, %f16;
	mul.f32 	%f18, %f8, %f17;
	fma.rn.f32 	%f19, %f18, %f7, %f7;
	add.f32 	%f20, %f19, %f19;
	mov.f32 	%f21, 0f3FC90FDB;
	sub.f32 	%f22, %f21, %f19;
	selp.f32	%f23, %f20, %f22, %p6;
	setp.lt.f32	%p7, %f1, 0f00000000;
	mov.f32 	%f24, 0f40490FDB;
	sub.f32 	%f25, %f24, %f23;
	selp.f32	%f26, %f25, %f23, %p7;
	cvta.to.global.u64 	%rd6, %rd1;
	add.s64 	%rd7, %rd6, %rd4;
	st.global.f32 	[%rd7], %f26;

BB0_2:
	ret;
}


`
	settheta_ptx_50 = `
.version 6.4
.target sm_50
.address_size 64

	// .globl	settheta

.visible .entry settheta(
	.param .u64 settheta_param_0,
	.param .u64 settheta_param_1,
	.param .u64 settheta_param_2,
	.param .u64 settheta_param_3,
	.param .u32 settheta_param_4,
	.param .u32 settheta_param_5,
	.param .u32 settheta_param_6
)
{
	.reg .pred 	%p<8>;
	.reg .f32 	%f<27>;
	.reg .b32 	%r<18>;
	.reg .b64 	%rd<8>;


	ld.param.u64 	%rd1, [settheta_param_0];
	ld.param.u64 	%rd2, [settheta_param_3];
	ld.param.u32 	%r4, [settheta_param_4];
	ld.param.u32 	%r5, [settheta_param_5];
	ld.param.u32 	%r6, [settheta_param_6];
	mov.u32 	%r7, %ctaid.x;
	mov.u32 	%r8, %ntid.x;
	mov.u32 	%r9, %tid.x;
	mad.lo.s32 	%r1, %r8, %r7, %r9;
	mov.u32 	%r10, %ntid.y;
	mov.u32 	%r11, %ctaid.y;
	mov.u32 	%r12, %tid.y;
	mad.lo.s32 	%r2, %r10, %r11, %r12;
	mov.u32 	%r13, %ntid.z;
	mov.u32 	%r14, %ctaid.z;
	mov.u32 	%r15, %tid.z;
	mad.lo.s32 	%r3, %r13, %r14, %r15;
	setp.ge.s32	%p1, %r2, %r5;
	setp.ge.s32	%p2, %r1, %r4;
	or.pred  	%p3, %p1, %p2;
	setp.ge.s32	%p4, %r3, %r6;
	or.pred  	%p5, %p3, %p4;
	@%p5 bra 	BB0_2;

	cvta.to.global.u64 	%rd3, %rd2;
	mad.lo.s32 	%r16, %r3, %r5, %r2;
	mad.lo.s32 	%r17, %r16, %r4, %r1;
	mul.wide.s32 	%rd4, %r17, 4;
	add.s64 	%rd5, %rd3, %rd4;
	ld.global.nc.f32 	%f1, [%rd5];
	abs.f32 	%f2, %f1;
	mov.f32 	%f3, 0f3F800000;
	sub.f32 	%f4, %f3, %f2;
	mul.f32 	%f5, %f4, 0f3F000000;
	sqrt.rn.f32 	%f6, %f5;
	setp.gt.f32	%p6, %f2, 0f3F11EB85;
	selp.f32	%f7, %f6, %f2, %p6;
	mul.f32 	%f8, %f7, %f7;
	mov.f32 	%f9, 0f3C94D2E9;
	mov.f32 	%f10, 0f3D53F941;
	fma.rn.f32 	%f11, %f10, %f8, %f9;
	mov.f32 	%f12, 0f3D3F841F;
	fma.rn.f32 	%f13, %f11, %f8, %f12;
	mov.f32 	%f14, 0f3D994929;
	fma.rn.f32 	%f15, %f13, %f8, %f14;
	mov.f32 	%f16, 0f3E2AAB94;
	fma.rn.f32 	%f17, %f15, %f8, %f16;
	mul.f32 	%f18, %f8, %f17;
	fma.rn.f32 	%f19, %f18, %f7, %f7;
	add.f32 	%f20, %f19, %f19;
	mov.f32 	%f21, 0f3FC90FDB;
	sub.f32 	%f22, %f21, %f19;
	selp.f32	%f23, %f20, %f22, %p6;
	setp.lt.f32	%p7, %f1, 0f00000000;
	mov.f32 	%f24, 0f40490FDB;
	sub.f32 	%f25, %f24, %f23;
	selp.f32	%f26, %f25, %f23, %p7;
	cvta.to.global.u64 	%rd6, %rd1;
	add.s64 	%rd7, %rd6, %rd4;
	st.global.f32 	[%rd7], %f26;

BB0_2:
	ret;
}


`
	settheta_ptx_52 = `
.version 6.4
.target sm_52
.address_size 64

	// .globl	settheta

.visible .entry settheta(
	.param .u64 settheta_param_0,
	.param .u64 settheta_param_1,
	.param .u64 settheta_param_2,
	.param .u64 settheta_param_3,
	.param .u32 settheta_param_4,
	.param .u32 settheta_param_5,
	.param .u32 settheta_param_6
)
{
	.reg .pred 	%p<8>;
	.reg .f32 	%f<27>;
	.reg .b32 	%r<18>;
	.reg .b64 	%rd<8>;


	ld.param.u64 	%rd1, [settheta_param_0];
	ld.param.u64 	%rd2, [settheta_param_3];
	ld.param.u32 	%r4, [settheta_param_4];
	ld.param.u32 	%r5, [settheta_param_5];
	ld.param.u32 	%r6, [settheta_param_6];
	mov.u32 	%r7, %ctaid.x;
	mov.u32 	%r8, %ntid.x;
	mov.u32 	%r9, %tid.x;
	mad.lo.s32 	%r1, %r8, %r7, %r9;
	mov.u32 	%r10, %ntid.y;
	mov.u32 	%r11, %ctaid.y;
	mov.u32 	%r12, %tid.y;
	mad.lo.s32 	%r2, %r10, %r11, %r12;
	mov.u32 	%r13, %ntid.z;
	mov.u32 	%r14, %ctaid.z;
	mov.u32 	%r15, %tid.z;
	mad.lo.s32 	%r3, %r13, %r14, %r15;
	setp.ge.s32	%p1, %r2, %r5;
	setp.ge.s32	%p2, %r1, %r4;
	or.pred  	%p3, %p1, %p2;
	setp.ge.s32	%p4, %r3, %r6;
	or.pred  	%p5, %p3, %p4;
	@%p5 bra 	BB0_2;

	cvta.to.global.u64 	%rd3, %rd2;
	mad.lo.s32 	%r16, %r3, %r5, %r2;
	mad.lo.s32 	%r17, %r16, %r4, %r1;
	mul.wide.s32 	%rd4, %r17, 4;
	add.s64 	%rd5, %rd3, %rd4;
	ld.global.nc.f32 	%f1, [%rd5];
	abs.f32 	%f2, %f1;
	mov.f32 	%f3, 0f3F800000;
	sub.f32 	%f4, %f3, %f2;
	mul.f32 	%f5, %f4, 0f3F000000;
	sqrt.rn.f32 	%f6, %f5;
	setp.gt.f32	%p6, %f2, 0f3F11EB85;
	selp.f32	%f7, %f6, %f2, %p6;
	mul.f32 	%f8, %f7, %f7;
	mov.f32 	%f9, 0f3C94D2E9;
	mov.f32 	%f10, 0f3D53F941;
	fma.rn.f32 	%f11, %f10, %f8, %f9;
	mov.f32 	%f12, 0f3D3F841F;
	fma.rn.f32 	%f13, %f11, %f8, %f12;
	mov.f32 	%f14, 0f3D994929;
	fma.rn.f32 	%f15, %f13, %f8, %f14;
	mov.f32 	%f16, 0f3E2AAB94;
	fma.rn.f32 	%f17, %f15, %f8, %f16;
	mul.f32 	%f18, %f8, %f17;
	fma.rn.f32 	%f19, %f18, %f7, %f7;
	add.f32 	%f20, %f19, %f19;
	mov.f32 	%f21, 0f3FC90FDB;
	sub.f32 	%f22, %f21, %f19;
	selp.f32	%f23, %f20, %f22, %p6;
	setp.lt.f32	%p7, %f1, 0f00000000;
	mov.f32 	%f24, 0f40490FDB;
	sub.f32 	%f25, %f24, %f23;
	selp.f32	%f26, %f25, %f23, %p7;
	cvta.to.global.u64 	%rd6, %rd1;
	add.s64 	%rd7, %rd6, %rd4;
	st.global.f32 	[%rd7], %f26;

BB0_2:
	ret;
}


`
	settheta_ptx_53 = `
.version 6.4
.target sm_53
.address_size 64

	// .globl	settheta

.visible .entry settheta(
	.param .u64 settheta_param_0,
	.param .u64 settheta_param_1,
	.param .u64 settheta_param_2,
	.param .u64 settheta_param_3,
	.param .u32 settheta_param_4,
	.param .u32 settheta_param_5,
	.param .u32 settheta_param_6
)
{
	.reg .pred 	%p<8>;
	.reg .f32 	%f<27>;
	.reg .b32 	%r<18>;
	.reg .b64 	%rd<8>;


	ld.param.u64 	%rd1, [settheta_param_0];
	ld.param.u64 	%rd2, [settheta_param_3];
	ld.param.u32 	%r4, [settheta_param_4];
	ld.param.u32 	%r5, [settheta_param_5];
	ld.param.u32 	%r6, [settheta_param_6];
	mov.u32 	%r7, %ctaid.x;
	mov.u32 	%r8, %ntid.x;
	mov.u32 	%r9, %tid.x;
	mad.lo.s32 	%r1, %r8, %r7, %r9;
	mov.u32 	%r10, %ntid.y;
	mov.u32 	%r11, %ctaid.y;
	mov.u32 	%r12, %tid.y;
	mad.lo.s32 	%r2, %r10, %r11, %r12;
	mov.u32 	%r13, %ntid.z;
	mov.u32 	%r14, %ctaid.z;
	mov.u32 	%r15, %tid.z;
	mad.lo.s32 	%r3, %r13, %r14, %r15;
	setp.ge.s32	%p1, %r2, %r5;
	setp.ge.s32	%p2, %r1, %r4;
	or.pred  	%p3, %p1, %p2;
	setp.ge.s32	%p4, %r3, %r6;
	or.pred  	%p5, %p3, %p4;
	@%p5 bra 	BB0_2;

	cvta.to.global.u64 	%rd3, %rd2;
	mad.lo.s32 	%r16, %r3, %r5, %r2;
	mad.lo.s32 	%r17, %r16, %r4, %r1;
	mul.wide.s32 	%rd4, %r17, 4;
	add.s64 	%rd5, %rd3, %rd4;
	ld.global.nc.f32 	%f1, [%rd5];
	abs.f32 	%f2, %f1;
	mov.f32 	%f3, 0f3F800000;
	sub.f32 	%f4, %f3, %f2;
	mul.f32 	%f5, %f4, 0f3F000000;
	sqrt.rn.f32 	%f6, %f5;
	setp.gt.f32	%p6, %f2, 0f3F11EB85;
	selp.f32	%f7, %f6, %f2, %p6;
	mul.f32 	%f8, %f7, %f7;
	mov.f32 	%f9, 0f3C94D2E9;
	mov.f32 	%f10, 0f3D53F941;
	fma.rn.f32 	%f11, %f10, %f8, %f9;
	mov.f32 	%f12, 0f3D3F841F;
	fma.rn.f32 	%f13, %f11, %f8, %f12;
	mov.f32 	%f14, 0f3D994929;
	fma.rn.f32 	%f15, %f13, %f8, %f14;
	mov.f32 	%f16, 0f3E2AAB94;
	fma.rn.f32 	%f17, %f15, %f8, %f16;
	mul.f32 	%f18, %f8, %f17;
	fma.rn.f32 	%f19, %f18, %f7, %f7;
	add.f32 	%f20, %f19, %f19;
	mov.f32 	%f21, 0f3FC90FDB;
	sub.f32 	%f22, %f21, %f19;
	selp.f32	%f23, %f20, %f22, %p6;
	setp.lt.f32	%p7, %f1, 0f00000000;
	mov.f32 	%f24, 0f40490FDB;
	sub.f32 	%f25, %f24, %f23;
	selp.f32	%f26, %f25, %f23, %p7;
	cvta.to.global.u64 	%rd6, %rd1;
	add.s64 	%rd7, %rd6, %rd4;
	st.global.f32 	[%rd7], %f26;

BB0_2:
	ret;
}


`
	settheta_ptx_60 = `
.version 6.4
.target sm_60
.address_size 64

	// .globl	settheta

.visible .entry settheta(
	.param .u64 settheta_param_0,
	.param .u64 settheta_param_1,
	.param .u64 settheta_param_2,
	.param .u64 settheta_param_3,
	.param .u32 settheta_param_4,
	.param .u32 settheta_param_5,
	.param .u32 settheta_param_6
)
{
	.reg .pred 	%p<8>;
	.reg .f32 	%f<27>;
	.reg .b32 	%r<18>;
	.reg .b64 	%rd<8>;


	ld.param.u64 	%rd1, [settheta_param_0];
	ld.param.u64 	%rd2, [settheta_param_3];
	ld.param.u32 	%r4, [settheta_param_4];
	ld.param.u32 	%r5, [settheta_param_5];
	ld.param.u32 	%r6, [settheta_param_6];
	mov.u32 	%r7, %ctaid.x;
	mov.u32 	%r8, %ntid.x;
	mov.u32 	%r9, %tid.x;
	mad.lo.s32 	%r1, %r8, %r7, %r9;
	mov.u32 	%r10, %ntid.y;
	mov.u32 	%r11, %ctaid.y;
	mov.u32 	%r12, %tid.y;
	mad.lo.s32 	%r2, %r10, %r11, %r12;
	mov.u32 	%r13, %ntid.z;
	mov.u32 	%r14, %ctaid.z;
	mov.u32 	%r15, %tid.z;
	mad.lo.s32 	%r3, %r13, %r14, %r15;
	setp.ge.s32	%p1, %r2, %r5;
	setp.ge.s32	%p2, %r1, %r4;
	or.pred  	%p3, %p1, %p2;
	setp.ge.s32	%p4, %r3, %r6;
	or.pred  	%p5, %p3, %p4;
	@%p5 bra 	BB0_2;

	cvta.to.global.u64 	%rd3, %rd2;
	mad.lo.s32 	%r16, %r3, %r5, %r2;
	mad.lo.s32 	%r17, %r16, %r4, %r1;
	mul.wide.s32 	%rd4, %r17, 4;
	add.s64 	%rd5, %rd3, %rd4;
	ld.global.nc.f32 	%f1, [%rd5];
	abs.f32 	%f2, %f1;
	mov.f32 	%f3, 0f3F800000;
	sub.f32 	%f4, %f3, %f2;
	mul.f32 	%f5, %f4, 0f3F000000;
	sqrt.rn.f32 	%f6, %f5;
	setp.gt.f32	%p6, %f2, 0f3F11EB85;
	selp.f32	%f7, %f6, %f2, %p6;
	mul.f32 	%f8, %f7, %f7;
	mov.f32 	%f9, 0f3C94D2E9;
	mov.f32 	%f10, 0f3D53F941;
	fma.rn.f32 	%f11, %f10, %f8, %f9;
	mov.f32 	%f12, 0f3D3F841F;
	fma.rn.f32 	%f13, %f11, %f8, %f12;
	mov.f32 	%f14, 0f3D994929;
	fma.rn.f32 	%f15, %f13, %f8, %f14;
	mov.f32 	%f16, 0f3E2AAB94;
	fma.rn.f32 	%f17, %f15, %f8, %f16;
	mul.f32 	%f18, %f8, %f17;
	fma.rn.f32 	%f19, %f18, %f7, %f7;
	add.f32 	%f20, %f19, %f19;
	mov.f32 	%f21, 0f3FC90FDB;
	sub.f32 	%f22, %f21, %f19;
	selp.f32	%f23, %f20, %f22, %p6;
	setp.lt.f32	%p7, %f1, 0f00000000;
	mov.f32 	%f24, 0f40490FDB;
	sub.f32 	%f25, %f24, %f23;
	selp.f32	%f26, %f25, %f23, %p7;
	cvta.to.global.u64 	%rd6, %rd1;
	add.s64 	%rd7, %rd6, %rd4;
	st.global.f32 	[%rd7], %f26;

BB0_2:
	ret;
}


`
	settheta_ptx_61 = `
.version 6.4
.target sm_61
.address_size 64

	// .globl	settheta

.visible .entry settheta(
	.param .u64 settheta_param_0,
	.param .u64 settheta_param_1,
	.param .u64 settheta_param_2,
	.param .u64 settheta_param_3,
	.param .u32 settheta_param_4,
	.param .u32 settheta_param_5,
	.param .u32 settheta_param_6
)
{
	.reg .pred 	%p<8>;
	.reg .f32 	%f<27>;
	.reg .b32 	%r<18>;
	.reg .b64 	%rd<8>;


	ld.param.u64 	%rd1, [settheta_param_0];
	ld.param.u64 	%rd2, [settheta_param_3];
	ld.param.u32 	%r4, [settheta_param_4];
	ld.param.u32 	%r5, [settheta_param_5];
	ld.param.u32 	%r6, [settheta_param_6];
	mov.u32 	%r7, %ctaid.x;
	mov.u32 	%r8, %ntid.x;
	mov.u32 	%r9, %tid.x;
	mad.lo.s32 	%r1, %r8, %r7, %r9;
	mov.u32 	%r10, %ntid.y;
	mov.u32 	%r11, %ctaid.y;
	mov.u32 	%r12, %tid.y;
	mad.lo.s32 	%r2, %r10, %r11, %r12;
	mov.u32 	%r13, %ntid.z;
	mov.u32 	%r14, %ctaid.z;
	mov.u32 	%r15, %tid.z;
	mad.lo.s32 	%r3, %r13, %r14, %r15;
	setp.ge.s32	%p1, %r2, %r5;
	setp.ge.s32	%p2, %r1, %r4;
	or.pred  	%p3, %p1, %p2;
	setp.ge.s32	%p4, %r3, %r6;
	or.pred  	%p5, %p3, %p4;
	@%p5 bra 	BB0_2;

	cvta.to.global.u64 	%rd3, %rd2;
	mad.lo.s32 	%r16, %r3, %r5, %r2;
	mad.lo.s32 	%r17, %r16, %r4, %r1;
	mul.wide.s32 	%rd4, %r17, 4;
	add.s64 	%rd5, %rd3, %rd4;
	ld.global.nc.f32 	%f1, [%rd5];
	abs.f32 	%f2, %f1;
	mov.f32 	%f3, 0f3F800000;
	sub.f32 	%f4, %f3, %f2;
	mul.f32 	%f5, %f4, 0f3F000000;
	sqrt.rn.f32 	%f6, %f5;
	setp.gt.f32	%p6, %f2, 0f3F11EB85;
	selp.f32	%f7, %f6, %f2, %p6;
	mul.f32 	%f8, %f7, %f7;
	mov.f32 	%f9, 0f3C94D2E9;
	mov.f32 	%f10, 0f3D53F941;
	fma.rn.f32 	%f11, %f10, %f8, %f9;
	mov.f32 	%f12, 0f3D3F841F;
	fma.rn.f32 	%f13, %f11, %f8, %f12;
	mov.f32 	%f14, 0f3D994929;
	fma.rn.f32 	%f15, %f13, %f8, %f14;
	mov.f32 	%f16, 0f3E2AAB94;
	fma.rn.f32 	%f17, %f15, %f8, %f16;
	mul.f32 	%f18, %f8, %f17;
	fma.rn.f32 	%f19, %f18, %f7, %f7;
	add.f32 	%f20, %f19, %f19;
	mov.f32 	%f21, 0f3FC90FDB;
	sub.f32 	%f22, %f21, %f19;
	selp.f32	%f23, %f20, %f22, %p6;
	setp.lt.f32	%p7, %f1, 0f00000000;
	mov.f32 	%f24, 0f40490FDB;
	sub.f32 	%f25, %f24, %f23;
	selp.f32	%f26, %f25, %f23, %p7;
	cvta.to.global.u64 	%rd6, %rd1;
	add.s64 	%rd7, %rd6, %rd4;
	st.global.f32 	[%rd7], %f26;

BB0_2:
	ret;
}


`
	settheta_ptx_70 = `
.version 6.4
.target sm_70
.address_size 64

	// .globl	settheta

.visible .entry settheta(
	.param .u64 settheta_param_0,
	.param .u64 settheta_param_1,
	.param .u64 settheta_param_2,
	.param .u64 settheta_param_3,
	.param .u32 settheta_param_4,
	.param .u32 settheta_param_5,
	.param .u32 settheta_param_6
)
{
	.reg .pred 	%p<8>;
	.reg .f32 	%f<27>;
	.reg .b32 	%r<18>;
	.reg .b64 	%rd<8>;


	ld.param.u64 	%rd1, [settheta_param_0];
	ld.param.u64 	%rd2, [settheta_param_3];
	ld.param.u32 	%r4, [settheta_param_4];
	ld.param.u32 	%r5, [settheta_param_5];
	ld.param.u32 	%r6, [settheta_param_6];
	mov.u32 	%r7, %ctaid.x;
	mov.u32 	%r8, %ntid.x;
	mov.u32 	%r9, %tid.x;
	mad.lo.s32 	%r1, %r8, %r7, %r9;
	mov.u32 	%r10, %ntid.y;
	mov.u32 	%r11, %ctaid.y;
	mov.u32 	%r12, %tid.y;
	mad.lo.s32 	%r2, %r10, %r11, %r12;
	mov.u32 	%r13, %ntid.z;
	mov.u32 	%r14, %ctaid.z;
	mov.u32 	%r15, %tid.z;
	mad.lo.s32 	%r3, %r13, %r14, %r15;
	setp.ge.s32	%p1, %r2, %r5;
	setp.ge.s32	%p2, %r1, %r4;
	or.pred  	%p3, %p1, %p2;
	setp.ge.s32	%p4, %r3, %r6;
	or.pred  	%p5, %p3, %p4;
	@%p5 bra 	BB0_2;

	cvta.to.global.u64 	%rd3, %rd2;
	mad.lo.s32 	%r16, %r3, %r5, %r2;
	mad.lo.s32 	%r17, %r16, %r4, %r1;
	mul.wide.s32 	%rd4, %r17, 4;
	add.s64 	%rd5, %rd3, %rd4;
	ld.global.nc.f32 	%f1, [%rd5];
	abs.f32 	%f2, %f1;
	mov.f32 	%f3, 0f3F800000;
	sub.f32 	%f4, %f3, %f2;
	mul.f32 	%f5, %f4, 0f3F000000;
	sqrt.rn.f32 	%f6, %f5;
	setp.gt.f32	%p6, %f2, 0f3F11EB85;
	selp.f32	%f7, %f6, %f2, %p6;
	mul.f32 	%f8, %f7, %f7;
	mov.f32 	%f9, 0f3C94D2E9;
	mov.f32 	%f10, 0f3D53F941;
	fma.rn.f32 	%f11, %f10, %f8, %f9;
	mov.f32 	%f12, 0f3D3F841F;
	fma.rn.f32 	%f13, %f11, %f8, %f12;
	mov.f32 	%f14, 0f3D994929;
	fma.rn.f32 	%f15, %f13, %f8, %f14;
	mov.f32 	%f16, 0f3E2AAB94;
	fma.rn.f32 	%f17, %f15, %f8, %f16;
	mul.f32 	%f18, %f8, %f17;
	fma.rn.f32 	%f19, %f18, %f7, %f7;
	add.f32 	%f20, %f19, %f19;
	mov.f32 	%f21, 0f3FC90FDB;
	sub.f32 	%f22, %f21, %f19;
	selp.f32	%f23, %f20, %f22, %p6;
	setp.lt.f32	%p7, %f1, 0f00000000;
	mov.f32 	%f24, 0f40490FDB;
	sub.f32 	%f25, %f24, %f23;
	selp.f32	%f26, %f25, %f23, %p7;
	cvta.to.global.u64 	%rd6, %rd1;
	add.s64 	%rd7, %rd6, %rd4;
	st.global.f32 	[%rd7], %f26;

BB0_2:
	ret;
}


`
	settheta_ptx_75 = `
.version 6.4
.target sm_75
.address_size 64

	// .globl	settheta

.visible .entry settheta(
	.param .u64 settheta_param_0,
	.param .u64 settheta_param_1,
	.param .u64 settheta_param_2,
	.param .u64 settheta_param_3,
	.param .u32 settheta_param_4,
	.param .u32 settheta_param_5,
	.param .u32 settheta_param_6
)
{
	.reg .pred 	%p<8>;
	.reg .f32 	%f<27>;
	.reg .b32 	%r<18>;
	.reg .b64 	%rd<8>;


	ld.param.u64 	%rd1, [settheta_param_0];
	ld.param.u64 	%rd2, [settheta_param_3];
	ld.param.u32 	%r4, [settheta_param_4];
	ld.param.u32 	%r5, [settheta_param_5];
	ld.param.u32 	%r6, [settheta_param_6];
	mov.u32 	%r7, %ctaid.x;
	mov.u32 	%r8, %ntid.x;
	mov.u32 	%r9, %tid.x;
	mad.lo.s32 	%r1, %r8, %r7, %r9;
	mov.u32 	%r10, %ntid.y;
	mov.u32 	%r11, %ctaid.y;
	mov.u32 	%r12, %tid.y;
	mad.lo.s32 	%r2, %r10, %r11, %r12;
	mov.u32 	%r13, %ntid.z;
	mov.u32 	%r14, %ctaid.z;
	mov.u32 	%r15, %tid.z;
	mad.lo.s32 	%r3, %r13, %r14, %r15;
	setp.ge.s32	%p1, %r2, %r5;
	setp.ge.s32	%p2, %r1, %r4;
	or.pred  	%p3, %p1, %p2;
	setp.ge.s32	%p4, %r3, %r6;
	or.pred  	%p5, %p3, %p4;
	@%p5 bra 	BB0_2;

	cvta.to.global.u64 	%rd3, %rd2;
	mad.lo.s32 	%r16, %r3, %r5, %r2;
	mad.lo.s32 	%r17, %r16, %r4, %r1;
	mul.wide.s32 	%rd4, %r17, 4;
	add.s64 	%rd5, %rd3, %rd4;
	ld.global.nc.f32 	%f1, [%rd5];
	abs.f32 	%f2, %f1;
	mov.f32 	%f3, 0f3F800000;
	sub.f32 	%f4, %f3, %f2;
	mul.f32 	%f5, %f4, 0f3F000000;
	sqrt.rn.f32 	%f6, %f5;
	setp.gt.f32	%p6, %f2, 0f3F11EB85;
	selp.f32	%f7, %f6, %f2, %p6;
	mul.f32 	%f8, %f7, %f7;
	mov.f32 	%f9, 0f3C94D2E9;
	mov.f32 	%f10, 0f3D53F941;
	fma.rn.f32 	%f11, %f10, %f8, %f9;
	mov.f32 	%f12, 0f3D3F841F;
	fma.rn.f32 	%f13, %f11, %f8, %f12;
	mov.f32 	%f14, 0f3D994929;
	fma.rn.f32 	%f15, %f13, %f8, %f14;
	mov.f32 	%f16, 0f3E2AAB94;
	fma.rn.f32 	%f17, %f15, %f8, %f16;
	mul.f32 	%f18, %f8, %f17;
	fma.rn.f32 	%f19, %f18, %f7, %f7;
	add.f32 	%f20, %f19, %f19;
	mov.f32 	%f21, 0f3FC90FDB;
	sub.f32 	%f22, %f21, %f19;
	selp.f32	%f23, %f20, %f22, %p6;
	setp.lt.f32	%p7, %f1, 0f00000000;
	mov.f32 	%f24, 0f40490FDB;
	sub.f32 	%f25, %f24, %f23;
	selp.f32	%f26, %f25, %f23, %p7;
	cvta.to.global.u64 	%rd6, %rd1;
	add.s64 	%rd7, %rd6, %rd4;
	st.global.f32 	[%rd7], %f26;

BB0_2:
	ret;
}


`
)
