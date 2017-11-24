#ifndef SIMD_H
#define SIMD_H

#if defined(__clang__)
	// TODO
#elif defined(__GNUC__)
	#include <emmintrin.h>
	// typedef float float4 __attribute__ ((vector_size(16)));
	typedef __m128 float4;
#elif defined(_MSC_VER)
	// TODO
#endif

#endif
