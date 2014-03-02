/*
	Copyright (c) 2012 Christopher A. Taylor.  All rights reserved.

	Redistribution and use in source and binary forms, with or without
	modification, are permitted provided that the following conditions are met:

	* Redistributions of source code must retain the above copyright notice,
	  this list of conditions and the following disclaimer.
	* Redistributions in binary form must reproduce the above copyright notice,
	  this list of conditions and the following disclaimer in the documentation
	  and/or other materials provided with the distribution.
	* Neither the name of LibCat nor the names of its contributors may be used
	  to endorse or promote products derived from this software without
	  specific prior written permission.

	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
	AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
	IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
	ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
	LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
	CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
	SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
	INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
	CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
	ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
	POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef CAT_BIG_MATH_HPP
#define CAT_BIG_MATH_HPP

#include "Platform.hpp"

namespace cat {


/*
	128-bit Arithmetic

	This wraps the best ways on each compiler and architecture to emulate 128-bit
	arithmetic, assuming that a 64-bit type is provided by the compiler.
*/

#ifndef CAT_HAS_U128
	struct u128 {
		u64 hi, lo;
	};
#endif

// r += x
CAT_INLINE void u128_add(u128 &r, const u128 x);
CAT_INLINE void u128_add(u128 &r, const u64 x);

// r -= x
CAT_INLINE void u128_sub(u128 &r, const u128 x);
CAT_INLINE void u128_sub(u128 &r, const u64 x);

// r <<= shift
CAT_INLINE void u128_lshift(u128 &r, int shift);

// x * y
CAT_INLINE u128 u128_mul(const u64 x, const u64 y);

// x * y + z
CAT_INLINE u128 u128_mul_add(const u64 x, const u64 y, const u64 z);

// x * y assuming MSB(x) = MSB(y) = 0
CAT_INLINE u128 u128_mul_63(const u64 x, const u64 y);

// get high/low halves
CAT_INLINE u64 u128_high(const u128 x);
CAT_INLINE u64 u128_low(const u128 x);


/*
 * There are three main platforms to support:
 *
 * + 64-bit GCC/Clang with 128-bit emulated types
 * + 64-bit ICC/MSVC without 128-bit emulated types
 * + 32-bit with 64-bit emulated types
 */

#if defined(CAT_HAS_U128)
# include "BigMath128.inh"
#elif defined(CAT_WORD_64)
# include "BigMath64.inh"
#else
# include "BigMath32.inh"
#endif


} // namespace cat

#endif // CAT_BIG_MATH_HPP

