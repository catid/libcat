/*
	Copyright (c) 2013-2014 Christopher A. Taylor.  All rights reserved.

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

#ifndef CAT_GALOIS_256_HPP
#define CAT_GALOIS_256_HPP

#include "Platform.hpp"

/*
	GF(256) Multiply and Divide functions

	Branchless multiply and divide construction from
	"Fast Software Implementations of Finite Field Operations (Extended Abstract)"
	by Cheng Huang and Lihao Xu

	Small corrections made to paper (Q = 255):
		+ The EXP_TABLE needs to have 512*2+1 elements to handle 0*0 = 0 case.
		+ Element 255*2 should be set to 1.

	After these corrections it works properly and reduces the execution time
	to 58% of the usual version that uses branches to handle zero input.

	These tables were generated using polynomial 0x15F.  Maybe it's voodoo but
	random GF(256) matrices with this polynomial tended to be more invertible.
	There are 16 generator polynomials for GF(256), and 0x1F5 was a close second
	in terms of rate of invertibility.

	GF256_INV_TABLE[x] was also generated to accelerate GF256Divide(1, x).

	This is to my knowledge the fastest possible approach, which is further
	supported by benchmarks in the independent implementation GF-Complete:

		http://web.eecs.utk.edu/~plank/plank/papers/CS-13-716.html
*/

namespace cat {


// Precomputed tables
extern const u16 GF256_LOG_TABLE[256];
extern const u8 GF256_EXP_TABLE[512*2+1];
extern const u8 GF256_INV_TABLE[256];

// Generated tables
extern u8 * CAT_RESTRICT GF256_MUL_TABLE;
extern u8 * CAT_RESTRICT GF256_DIV_TABLE;

// Call this function to initialize the tables
void GF256Init();

// return x * y in GF(256)
// For repeated multiplication by a constant, it is faster to put the constant in y.
static CAT_INLINE u8 GF256Multiply(u8 x, u8 y) {
	return GF256_MUL_TABLE[((u32)y << 8) + x];
}

// return x / y in GF(256)
// Memory-access optimized for constant divisors in y.
static CAT_INLINE u8 GF256Divide(u8 x, u8 y) {
	return GF256_DIV_TABLE[((u32)y << 8) + x];
}

// Performs "dest[] += src[] * x" operation in GF(256)
extern void GF256MemMulAdd(void * CAT_RESTRICT vdest, u8 x, const void * CAT_RESTRICT vsrc, int bytes);

// Performs "dest[] /= x" operation in GF(256)
extern void GF256MemDivide(void * CAT_RESTRICT vdest, u8 x, int bytes);


} // namespace cat

#endif // CAT_GALOIS_256_HPP

