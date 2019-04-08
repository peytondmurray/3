package cuda

/*
 THIS FILE IS AUTO-GENERATED BY CUDA2GO.
 EDITING IS FUTILE.
*/

import(
	"unsafe"
	"github.com/mumax/3/cuda/cu"
	"github.com/mumax/3/timer"
	"sync"
)

// CUDA handle for regiondecode kernel
var regiondecode_code cu.Function

// Stores the arguments for regiondecode kernel invocation
type regiondecode_args_t struct{
	 arg_dst unsafe.Pointer
	 arg_LUT unsafe.Pointer
	 arg_regions unsafe.Pointer
	 arg_N int
	 argptr [4]unsafe.Pointer
	sync.Mutex
}

// Stores the arguments for regiondecode kernel invocation
var regiondecode_args regiondecode_args_t

func init(){
	// CUDA driver kernel call wants pointers to arguments, set them up once.
	 regiondecode_args.argptr[0] = unsafe.Pointer(&regiondecode_args.arg_dst)
	 regiondecode_args.argptr[1] = unsafe.Pointer(&regiondecode_args.arg_LUT)
	 regiondecode_args.argptr[2] = unsafe.Pointer(&regiondecode_args.arg_regions)
	 regiondecode_args.argptr[3] = unsafe.Pointer(&regiondecode_args.arg_N)
	 }

// Wrapper for regiondecode CUDA kernel, asynchronous.
func k_regiondecode_async ( dst unsafe.Pointer, LUT unsafe.Pointer, regions unsafe.Pointer, N int,  cfg *config) {
	if Synchronous{ // debug
		Sync()
		timer.Start("regiondecode")
	}

	regiondecode_args.Lock()
	defer regiondecode_args.Unlock()

	if regiondecode_code == 0{
		regiondecode_code = fatbinLoad(regiondecode_map, "regiondecode")
	}

	 regiondecode_args.arg_dst = dst
	 regiondecode_args.arg_LUT = LUT
	 regiondecode_args.arg_regions = regions
	 regiondecode_args.arg_N = N
	

	args := regiondecode_args.argptr[:]
	cu.LaunchKernel(regiondecode_code, cfg.Grid.X, cfg.Grid.Y, cfg.Grid.Z, cfg.Block.X, cfg.Block.Y, cfg.Block.Z, 0, stream0, args)

	if Synchronous{ // debug
		Sync()
		timer.Stop("regiondecode")
	}
}

// maps compute capability on PTX code for regiondecode kernel.
var regiondecode_map = map[int]string{ 0: "" ,
30: regiondecode_ptx_30 ,
35: regiondecode_ptx_35 ,
37: regiondecode_ptx_37 ,
50: regiondecode_ptx_50 ,
52: regiondecode_ptx_52 ,
53: regiondecode_ptx_53 ,
60: regiondecode_ptx_60 ,
61: regiondecode_ptx_61 ,
70: regiondecode_ptx_70 ,
75: regiondecode_ptx_75  }

// regiondecode PTX code for various compute capabilities.
const(
  regiondecode_ptx_30 = `
.version 6.4
.target sm_30
.address_size 64

	// .globl	regiondecode

.visible .entry regiondecode(
	.param .u64 regiondecode_param_0,
	.param .u64 regiondecode_param_1,
	.param .u64 regiondecode_param_2,
	.param .u32 regiondecode_param_3
)
{
	.reg .pred 	%p<2>;
	.reg .f32 	%f<2>;
	.reg .b32 	%r<10>;
	.reg .b64 	%rd<13>;


	ld.param.u64 	%rd1, [regiondecode_param_0];
	ld.param.u64 	%rd2, [regiondecode_param_1];
	ld.param.u64 	%rd3, [regiondecode_param_2];
	ld.param.u32 	%r2, [regiondecode_param_3];
	mov.u32 	%r3, %ctaid.y;
	mov.u32 	%r4, %nctaid.x;
	mov.u32 	%r5, %ctaid.x;
	mad.lo.s32 	%r6, %r4, %r3, %r5;
	mov.u32 	%r7, %ntid.x;
	mov.u32 	%r8, %tid.x;
	mad.lo.s32 	%r1, %r6, %r7, %r8;
	setp.ge.s32	%p1, %r1, %r2;
	@%p1 bra 	BB0_2;

	cvta.to.global.u64 	%rd4, %rd3;
	cvt.s64.s32	%rd5, %r1;
	add.s64 	%rd6, %rd4, %rd5;
	cvta.to.global.u64 	%rd7, %rd2;
	ld.global.u8 	%r9, [%rd6];
	mul.wide.u32 	%rd8, %r9, 4;
	add.s64 	%rd9, %rd7, %rd8;
	ld.global.f32 	%f1, [%rd9];
	cvta.to.global.u64 	%rd10, %rd1;
	mul.wide.s32 	%rd11, %r1, 4;
	add.s64 	%rd12, %rd10, %rd11;
	st.global.f32 	[%rd12], %f1;

BB0_2:
	ret;
}


`
   regiondecode_ptx_35 = `
.version 6.4
.target sm_35
.address_size 64

	// .globl	regiondecode

.visible .entry regiondecode(
	.param .u64 regiondecode_param_0,
	.param .u64 regiondecode_param_1,
	.param .u64 regiondecode_param_2,
	.param .u32 regiondecode_param_3
)
{
	.reg .pred 	%p<2>;
	.reg .b16 	%rs<2>;
	.reg .f32 	%f<2>;
	.reg .b32 	%r<11>;
	.reg .b64 	%rd<13>;


	ld.param.u64 	%rd1, [regiondecode_param_0];
	ld.param.u64 	%rd2, [regiondecode_param_1];
	ld.param.u64 	%rd3, [regiondecode_param_2];
	ld.param.u32 	%r2, [regiondecode_param_3];
	mov.u32 	%r3, %ctaid.y;
	mov.u32 	%r4, %nctaid.x;
	mov.u32 	%r5, %ctaid.x;
	mad.lo.s32 	%r6, %r4, %r3, %r5;
	mov.u32 	%r7, %ntid.x;
	mov.u32 	%r8, %tid.x;
	mad.lo.s32 	%r1, %r6, %r7, %r8;
	setp.ge.s32	%p1, %r1, %r2;
	@%p1 bra 	BB0_2;

	cvta.to.global.u64 	%rd4, %rd3;
	cvt.s64.s32	%rd5, %r1;
	add.s64 	%rd6, %rd4, %rd5;
	ld.global.nc.u8 	%rs1, [%rd6];
	cvta.to.global.u64 	%rd7, %rd2;
	cvt.u32.u16	%r9, %rs1;
	and.b32  	%r10, %r9, 255;
	mul.wide.u32 	%rd8, %r10, 4;
	add.s64 	%rd9, %rd7, %rd8;
	ld.global.nc.f32 	%f1, [%rd9];
	cvta.to.global.u64 	%rd10, %rd1;
	mul.wide.s32 	%rd11, %r1, 4;
	add.s64 	%rd12, %rd10, %rd11;
	st.global.f32 	[%rd12], %f1;

BB0_2:
	ret;
}


`
   regiondecode_ptx_37 = `
.version 6.4
.target sm_37
.address_size 64

	// .globl	regiondecode

.visible .entry regiondecode(
	.param .u64 regiondecode_param_0,
	.param .u64 regiondecode_param_1,
	.param .u64 regiondecode_param_2,
	.param .u32 regiondecode_param_3
)
{
	.reg .pred 	%p<2>;
	.reg .b16 	%rs<2>;
	.reg .f32 	%f<2>;
	.reg .b32 	%r<11>;
	.reg .b64 	%rd<13>;


	ld.param.u64 	%rd1, [regiondecode_param_0];
	ld.param.u64 	%rd2, [regiondecode_param_1];
	ld.param.u64 	%rd3, [regiondecode_param_2];
	ld.param.u32 	%r2, [regiondecode_param_3];
	mov.u32 	%r3, %ctaid.y;
	mov.u32 	%r4, %nctaid.x;
	mov.u32 	%r5, %ctaid.x;
	mad.lo.s32 	%r6, %r4, %r3, %r5;
	mov.u32 	%r7, %ntid.x;
	mov.u32 	%r8, %tid.x;
	mad.lo.s32 	%r1, %r6, %r7, %r8;
	setp.ge.s32	%p1, %r1, %r2;
	@%p1 bra 	BB0_2;

	cvta.to.global.u64 	%rd4, %rd3;
	cvt.s64.s32	%rd5, %r1;
	add.s64 	%rd6, %rd4, %rd5;
	ld.global.nc.u8 	%rs1, [%rd6];
	cvta.to.global.u64 	%rd7, %rd2;
	cvt.u32.u16	%r9, %rs1;
	and.b32  	%r10, %r9, 255;
	mul.wide.u32 	%rd8, %r10, 4;
	add.s64 	%rd9, %rd7, %rd8;
	ld.global.nc.f32 	%f1, [%rd9];
	cvta.to.global.u64 	%rd10, %rd1;
	mul.wide.s32 	%rd11, %r1, 4;
	add.s64 	%rd12, %rd10, %rd11;
	st.global.f32 	[%rd12], %f1;

BB0_2:
	ret;
}


`
   regiondecode_ptx_50 = `
.version 6.4
.target sm_50
.address_size 64

	// .globl	regiondecode

.visible .entry regiondecode(
	.param .u64 regiondecode_param_0,
	.param .u64 regiondecode_param_1,
	.param .u64 regiondecode_param_2,
	.param .u32 regiondecode_param_3
)
{
	.reg .pred 	%p<2>;
	.reg .b16 	%rs<2>;
	.reg .f32 	%f<2>;
	.reg .b32 	%r<11>;
	.reg .b64 	%rd<13>;


	ld.param.u64 	%rd1, [regiondecode_param_0];
	ld.param.u64 	%rd2, [regiondecode_param_1];
	ld.param.u64 	%rd3, [regiondecode_param_2];
	ld.param.u32 	%r2, [regiondecode_param_3];
	mov.u32 	%r3, %ctaid.y;
	mov.u32 	%r4, %nctaid.x;
	mov.u32 	%r5, %ctaid.x;
	mad.lo.s32 	%r6, %r4, %r3, %r5;
	mov.u32 	%r7, %ntid.x;
	mov.u32 	%r8, %tid.x;
	mad.lo.s32 	%r1, %r6, %r7, %r8;
	setp.ge.s32	%p1, %r1, %r2;
	@%p1 bra 	BB0_2;

	cvta.to.global.u64 	%rd4, %rd3;
	cvt.s64.s32	%rd5, %r1;
	add.s64 	%rd6, %rd4, %rd5;
	ld.global.nc.u8 	%rs1, [%rd6];
	cvta.to.global.u64 	%rd7, %rd2;
	cvt.u32.u16	%r9, %rs1;
	and.b32  	%r10, %r9, 255;
	mul.wide.u32 	%rd8, %r10, 4;
	add.s64 	%rd9, %rd7, %rd8;
	ld.global.nc.f32 	%f1, [%rd9];
	cvta.to.global.u64 	%rd10, %rd1;
	mul.wide.s32 	%rd11, %r1, 4;
	add.s64 	%rd12, %rd10, %rd11;
	st.global.f32 	[%rd12], %f1;

BB0_2:
	ret;
}


`
   regiondecode_ptx_52 = `
.version 6.4
.target sm_52
.address_size 64

	// .globl	regiondecode

.visible .entry regiondecode(
	.param .u64 regiondecode_param_0,
	.param .u64 regiondecode_param_1,
	.param .u64 regiondecode_param_2,
	.param .u32 regiondecode_param_3
)
{
	.reg .pred 	%p<2>;
	.reg .b16 	%rs<2>;
	.reg .f32 	%f<2>;
	.reg .b32 	%r<11>;
	.reg .b64 	%rd<13>;


	ld.param.u64 	%rd1, [regiondecode_param_0];
	ld.param.u64 	%rd2, [regiondecode_param_1];
	ld.param.u64 	%rd3, [regiondecode_param_2];
	ld.param.u32 	%r2, [regiondecode_param_3];
	mov.u32 	%r3, %ctaid.y;
	mov.u32 	%r4, %nctaid.x;
	mov.u32 	%r5, %ctaid.x;
	mad.lo.s32 	%r6, %r4, %r3, %r5;
	mov.u32 	%r7, %ntid.x;
	mov.u32 	%r8, %tid.x;
	mad.lo.s32 	%r1, %r6, %r7, %r8;
	setp.ge.s32	%p1, %r1, %r2;
	@%p1 bra 	BB0_2;

	cvta.to.global.u64 	%rd4, %rd3;
	cvt.s64.s32	%rd5, %r1;
	add.s64 	%rd6, %rd4, %rd5;
	ld.global.nc.u8 	%rs1, [%rd6];
	cvta.to.global.u64 	%rd7, %rd2;
	cvt.u32.u16	%r9, %rs1;
	and.b32  	%r10, %r9, 255;
	mul.wide.u32 	%rd8, %r10, 4;
	add.s64 	%rd9, %rd7, %rd8;
	ld.global.nc.f32 	%f1, [%rd9];
	cvta.to.global.u64 	%rd10, %rd1;
	mul.wide.s32 	%rd11, %r1, 4;
	add.s64 	%rd12, %rd10, %rd11;
	st.global.f32 	[%rd12], %f1;

BB0_2:
	ret;
}


`
   regiondecode_ptx_53 = `
.version 6.4
.target sm_53
.address_size 64

	// .globl	regiondecode

.visible .entry regiondecode(
	.param .u64 regiondecode_param_0,
	.param .u64 regiondecode_param_1,
	.param .u64 regiondecode_param_2,
	.param .u32 regiondecode_param_3
)
{
	.reg .pred 	%p<2>;
	.reg .b16 	%rs<2>;
	.reg .f32 	%f<2>;
	.reg .b32 	%r<11>;
	.reg .b64 	%rd<13>;


	ld.param.u64 	%rd1, [regiondecode_param_0];
	ld.param.u64 	%rd2, [regiondecode_param_1];
	ld.param.u64 	%rd3, [regiondecode_param_2];
	ld.param.u32 	%r2, [regiondecode_param_3];
	mov.u32 	%r3, %ctaid.y;
	mov.u32 	%r4, %nctaid.x;
	mov.u32 	%r5, %ctaid.x;
	mad.lo.s32 	%r6, %r4, %r3, %r5;
	mov.u32 	%r7, %ntid.x;
	mov.u32 	%r8, %tid.x;
	mad.lo.s32 	%r1, %r6, %r7, %r8;
	setp.ge.s32	%p1, %r1, %r2;
	@%p1 bra 	BB0_2;

	cvta.to.global.u64 	%rd4, %rd3;
	cvt.s64.s32	%rd5, %r1;
	add.s64 	%rd6, %rd4, %rd5;
	ld.global.nc.u8 	%rs1, [%rd6];
	cvta.to.global.u64 	%rd7, %rd2;
	cvt.u32.u16	%r9, %rs1;
	and.b32  	%r10, %r9, 255;
	mul.wide.u32 	%rd8, %r10, 4;
	add.s64 	%rd9, %rd7, %rd8;
	ld.global.nc.f32 	%f1, [%rd9];
	cvta.to.global.u64 	%rd10, %rd1;
	mul.wide.s32 	%rd11, %r1, 4;
	add.s64 	%rd12, %rd10, %rd11;
	st.global.f32 	[%rd12], %f1;

BB0_2:
	ret;
}


`
   regiondecode_ptx_60 = `
.version 6.4
.target sm_60
.address_size 64

	// .globl	regiondecode

.visible .entry regiondecode(
	.param .u64 regiondecode_param_0,
	.param .u64 regiondecode_param_1,
	.param .u64 regiondecode_param_2,
	.param .u32 regiondecode_param_3
)
{
	.reg .pred 	%p<2>;
	.reg .b16 	%rs<2>;
	.reg .f32 	%f<2>;
	.reg .b32 	%r<11>;
	.reg .b64 	%rd<13>;


	ld.param.u64 	%rd1, [regiondecode_param_0];
	ld.param.u64 	%rd2, [regiondecode_param_1];
	ld.param.u64 	%rd3, [regiondecode_param_2];
	ld.param.u32 	%r2, [regiondecode_param_3];
	mov.u32 	%r3, %ctaid.y;
	mov.u32 	%r4, %nctaid.x;
	mov.u32 	%r5, %ctaid.x;
	mad.lo.s32 	%r6, %r4, %r3, %r5;
	mov.u32 	%r7, %ntid.x;
	mov.u32 	%r8, %tid.x;
	mad.lo.s32 	%r1, %r6, %r7, %r8;
	setp.ge.s32	%p1, %r1, %r2;
	@%p1 bra 	BB0_2;

	cvta.to.global.u64 	%rd4, %rd3;
	cvt.s64.s32	%rd5, %r1;
	add.s64 	%rd6, %rd4, %rd5;
	ld.global.nc.u8 	%rs1, [%rd6];
	cvta.to.global.u64 	%rd7, %rd2;
	cvt.u32.u16	%r9, %rs1;
	and.b32  	%r10, %r9, 255;
	mul.wide.u32 	%rd8, %r10, 4;
	add.s64 	%rd9, %rd7, %rd8;
	ld.global.nc.f32 	%f1, [%rd9];
	cvta.to.global.u64 	%rd10, %rd1;
	mul.wide.s32 	%rd11, %r1, 4;
	add.s64 	%rd12, %rd10, %rd11;
	st.global.f32 	[%rd12], %f1;

BB0_2:
	ret;
}


`
   regiondecode_ptx_61 = `
.version 6.4
.target sm_61
.address_size 64

	// .globl	regiondecode

.visible .entry regiondecode(
	.param .u64 regiondecode_param_0,
	.param .u64 regiondecode_param_1,
	.param .u64 regiondecode_param_2,
	.param .u32 regiondecode_param_3
)
{
	.reg .pred 	%p<2>;
	.reg .b16 	%rs<2>;
	.reg .f32 	%f<2>;
	.reg .b32 	%r<11>;
	.reg .b64 	%rd<13>;


	ld.param.u64 	%rd1, [regiondecode_param_0];
	ld.param.u64 	%rd2, [regiondecode_param_1];
	ld.param.u64 	%rd3, [regiondecode_param_2];
	ld.param.u32 	%r2, [regiondecode_param_3];
	mov.u32 	%r3, %ctaid.y;
	mov.u32 	%r4, %nctaid.x;
	mov.u32 	%r5, %ctaid.x;
	mad.lo.s32 	%r6, %r4, %r3, %r5;
	mov.u32 	%r7, %ntid.x;
	mov.u32 	%r8, %tid.x;
	mad.lo.s32 	%r1, %r6, %r7, %r8;
	setp.ge.s32	%p1, %r1, %r2;
	@%p1 bra 	BB0_2;

	cvta.to.global.u64 	%rd4, %rd3;
	cvt.s64.s32	%rd5, %r1;
	add.s64 	%rd6, %rd4, %rd5;
	ld.global.nc.u8 	%rs1, [%rd6];
	cvta.to.global.u64 	%rd7, %rd2;
	cvt.u32.u16	%r9, %rs1;
	and.b32  	%r10, %r9, 255;
	mul.wide.u32 	%rd8, %r10, 4;
	add.s64 	%rd9, %rd7, %rd8;
	ld.global.nc.f32 	%f1, [%rd9];
	cvta.to.global.u64 	%rd10, %rd1;
	mul.wide.s32 	%rd11, %r1, 4;
	add.s64 	%rd12, %rd10, %rd11;
	st.global.f32 	[%rd12], %f1;

BB0_2:
	ret;
}


`
   regiondecode_ptx_70 = `
.version 6.4
.target sm_70
.address_size 64

	// .globl	regiondecode

.visible .entry regiondecode(
	.param .u64 regiondecode_param_0,
	.param .u64 regiondecode_param_1,
	.param .u64 regiondecode_param_2,
	.param .u32 regiondecode_param_3
)
{
	.reg .pred 	%p<2>;
	.reg .b16 	%rs<2>;
	.reg .f32 	%f<2>;
	.reg .b32 	%r<11>;
	.reg .b64 	%rd<13>;


	ld.param.u64 	%rd1, [regiondecode_param_0];
	ld.param.u64 	%rd2, [regiondecode_param_1];
	ld.param.u64 	%rd3, [regiondecode_param_2];
	ld.param.u32 	%r2, [regiondecode_param_3];
	mov.u32 	%r3, %ctaid.y;
	mov.u32 	%r4, %nctaid.x;
	mov.u32 	%r5, %ctaid.x;
	mad.lo.s32 	%r6, %r4, %r3, %r5;
	mov.u32 	%r7, %ntid.x;
	mov.u32 	%r8, %tid.x;
	mad.lo.s32 	%r1, %r6, %r7, %r8;
	setp.ge.s32	%p1, %r1, %r2;
	@%p1 bra 	BB0_2;

	cvta.to.global.u64 	%rd4, %rd3;
	cvt.s64.s32	%rd5, %r1;
	add.s64 	%rd6, %rd4, %rd5;
	ld.global.nc.u8 	%rs1, [%rd6];
	cvta.to.global.u64 	%rd7, %rd2;
	cvt.u32.u16	%r9, %rs1;
	and.b32  	%r10, %r9, 255;
	mul.wide.u32 	%rd8, %r10, 4;
	add.s64 	%rd9, %rd7, %rd8;
	ld.global.nc.f32 	%f1, [%rd9];
	cvta.to.global.u64 	%rd10, %rd1;
	mul.wide.s32 	%rd11, %r1, 4;
	add.s64 	%rd12, %rd10, %rd11;
	st.global.f32 	[%rd12], %f1;

BB0_2:
	ret;
}


`
   regiondecode_ptx_75 = `
.version 6.4
.target sm_75
.address_size 64

	// .globl	regiondecode

.visible .entry regiondecode(
	.param .u64 regiondecode_param_0,
	.param .u64 regiondecode_param_1,
	.param .u64 regiondecode_param_2,
	.param .u32 regiondecode_param_3
)
{
	.reg .pred 	%p<2>;
	.reg .b16 	%rs<2>;
	.reg .f32 	%f<2>;
	.reg .b32 	%r<11>;
	.reg .b64 	%rd<13>;


	ld.param.u64 	%rd1, [regiondecode_param_0];
	ld.param.u64 	%rd2, [regiondecode_param_1];
	ld.param.u64 	%rd3, [regiondecode_param_2];
	ld.param.u32 	%r2, [regiondecode_param_3];
	mov.u32 	%r3, %ctaid.y;
	mov.u32 	%r4, %nctaid.x;
	mov.u32 	%r5, %ctaid.x;
	mad.lo.s32 	%r6, %r4, %r3, %r5;
	mov.u32 	%r7, %ntid.x;
	mov.u32 	%r8, %tid.x;
	mad.lo.s32 	%r1, %r6, %r7, %r8;
	setp.ge.s32	%p1, %r1, %r2;
	@%p1 bra 	BB0_2;

	cvta.to.global.u64 	%rd4, %rd3;
	cvt.s64.s32	%rd5, %r1;
	add.s64 	%rd6, %rd4, %rd5;
	ld.global.nc.u8 	%rs1, [%rd6];
	cvta.to.global.u64 	%rd7, %rd2;
	cvt.u32.u16	%r9, %rs1;
	and.b32  	%r10, %r9, 255;
	mul.wide.u32 	%rd8, %r10, 4;
	add.s64 	%rd9, %rd7, %rd8;
	ld.global.nc.f32 	%f1, [%rd9];
	cvta.to.global.u64 	%rd10, %rd1;
	mul.wide.s32 	%rd11, %r1, 4;
	add.s64 	%rd12, %rd10, %rd11;
	st.global.f32 	[%rd12], %f1;

BB0_2:
	ret;
}


`
 )
