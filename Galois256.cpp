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

#include "Galois256.hpp"
#include "MemXOR.hpp"
using namespace cat;

#include <stdlib.h> // malloc


//// Bootstrap tables

const u16 cat::GF256_LOG_TABLE[256] = {
	512, 255, 1, 122, 2, 244, 123, 181, 3, 48, 245, 224, 124, 84, 182, 111,
	4, 233, 49, 19, 246, 107, 225, 206, 125, 56, 85, 170, 183, 91, 112, 250,
	5, 117, 234, 10, 50, 156, 20, 213, 247, 203, 108, 178, 226, 37, 207, 210,
	126, 150, 57, 100, 86, 141, 171, 40, 184, 73, 92, 164, 113, 146, 251, 229,
	6, 96, 118, 15, 235, 193, 11, 13, 51, 68, 157, 195, 21, 31, 214, 237,
	248, 168, 204, 17, 109, 222, 179, 120, 227, 162, 38, 98, 208, 176, 211, 8,
	127, 188, 151, 239, 58, 132, 101, 216, 87, 80, 142, 33, 172, 27, 41, 23,
	185, 77, 74, 197, 93, 65, 165, 159, 114, 200, 147, 70, 252, 45, 230, 53,
	7, 175, 97, 161, 119, 221, 16, 167, 236, 30, 194, 67, 12, 192, 14, 95,
	52, 44, 69, 199, 158, 64, 196, 76, 22, 26, 32, 79, 215, 131, 238, 187,
	249, 90, 169, 55, 205, 106, 18, 232, 110, 83, 223, 47, 180, 243, 121, 254,
	228, 145, 163, 72, 39, 140, 99, 149, 209, 36, 177, 202, 212, 155, 9, 116,
	128, 61, 189, 218, 152, 137, 240, 103, 59, 135, 133, 134, 102, 136, 217, 60,
	88, 104, 81, 241, 143, 138, 34, 153, 173, 219, 28, 190, 42, 62, 24, 129,
	186, 130, 78, 25, 75, 63, 198, 43, 94, 191, 66, 29, 166, 220, 160, 174,
	115, 154, 201, 35, 148, 139, 71, 144, 253, 242, 46, 82, 231, 105, 54, 89
};

const u8 cat::GF256_EXP_TABLE[512*2+1] = {
	1, 2, 4, 8, 16, 32, 64, 128, 95, 190, 35, 70, 140, 71, 142, 67,
	134, 83, 166, 19, 38, 76, 152, 111, 222, 227, 153, 109, 218, 235, 137, 77,
	154, 107, 214, 243, 185, 45, 90, 180, 55, 110, 220, 231, 145, 125, 250, 171,
	9, 18, 36, 72, 144, 127, 254, 163, 25, 50, 100, 200, 207, 193, 221, 229,
	149, 117, 234, 139, 73, 146, 123, 246, 179, 57, 114, 228, 151, 113, 226, 155,
	105, 210, 251, 169, 13, 26, 52, 104, 208, 255, 161, 29, 58, 116, 232, 143,
	65, 130, 91, 182, 51, 102, 204, 199, 209, 253, 165, 21, 42, 84, 168, 15,
	30, 60, 120, 240, 191, 33, 66, 132, 87, 174, 3, 6, 12, 24, 48, 96,
	192, 223, 225, 157, 101, 202, 203, 201, 205, 197, 213, 245, 181, 53, 106, 212,
	247, 177, 61, 122, 244, 183, 49, 98, 196, 215, 241, 189, 37, 74, 148, 119,
	238, 131, 89, 178, 59, 118, 236, 135, 81, 162, 27, 54, 108, 216, 239, 129,
	93, 186, 43, 86, 172, 7, 14, 28, 56, 112, 224, 159, 97, 194, 219, 233,
	141, 69, 138, 75, 150, 115, 230, 147, 121, 242, 187, 41, 82, 164, 23, 46,
	92, 184, 47, 94, 188, 39, 78, 156, 103, 206, 195, 217, 237, 133, 85, 170,
	11, 22, 44, 88, 176, 63, 126, 252, 167, 17, 34, 68, 136, 79, 158, 99,
	198, 211, 249, 173, 5, 10, 20, 40, 80, 160, 31, 62, 124, 248, 175, 1,
	2, 4, 8, 16, 32, 64, 128, 95, 190, 35, 70, 140, 71, 142, 67, 134,
	83, 166, 19, 38, 76, 152, 111, 222, 227, 153, 109, 218, 235, 137, 77, 154,
	107, 214, 243, 185, 45, 90, 180, 55, 110, 220, 231, 145, 125, 250, 171, 9,
	18, 36, 72, 144, 127, 254, 163, 25, 50, 100, 200, 207, 193, 221, 229, 149,
	117, 234, 139, 73, 146, 123, 246, 179, 57, 114, 228, 151, 113, 226, 155, 105,
	210, 251, 169, 13, 26, 52, 104, 208, 255, 161, 29, 58, 116, 232, 143, 65,
	130, 91, 182, 51, 102, 204, 199, 209, 253, 165, 21, 42, 84, 168, 15, 30,
	60, 120, 240, 191, 33, 66, 132, 87, 174, 3, 6, 12, 24, 48, 96, 192,
	223, 225, 157, 101, 202, 203, 201, 205, 197, 213, 245, 181, 53, 106, 212, 247,
	177, 61, 122, 244, 183, 49, 98, 196, 215, 241, 189, 37, 74, 148, 119, 238,
	131, 89, 178, 59, 118, 236, 135, 81, 162, 27, 54, 108, 216, 239, 129, 93,
	186, 43, 86, 172, 7, 14, 28, 56, 112, 224, 159, 97, 194, 219, 233, 141,
	69, 138, 75, 150, 115, 230, 147, 121, 242, 187, 41, 82, 164, 23, 46, 92,
	184, 47, 94, 188, 39, 78, 156, 103, 206, 195, 217, 237, 133, 85, 170, 11,
	22, 44, 88, 176, 63, 126, 252, 167, 17, 34, 68, 136, 79, 158, 99, 198,
	211, 249, 173, 5, 10, 20, 40, 80, 160, 31, 62, 124, 248, 175, 1, 0,
};

const u8 cat::GF256_INV_TABLE[256] = {
	0, 1, 175, 202, 248, 70, 101, 114, 124, 46, 35, 77, 157, 54, 57, 247,
	62, 152, 23, 136, 190, 244, 137, 18, 225, 147, 27, 26, 179, 59, 212, 32,
	31, 213, 76, 10, 164, 182, 68, 220, 95, 144, 122, 113, 235, 195, 9, 125,
	223, 253, 230, 189, 162, 120, 13, 156, 246, 14, 178, 29, 106, 84, 16, 153,
	160, 119, 197, 198, 38, 221, 5, 249, 82, 159, 91, 207, 34, 11, 110, 166,
	128, 104, 72, 158, 61, 107, 151, 201, 218, 116, 206, 74, 171, 155, 145, 40,
	192, 139, 209, 134, 115, 6, 241, 180, 81, 129, 60, 85, 169, 176, 78, 167,
	123, 43, 7, 100, 89, 219, 161, 65, 53, 163, 42, 112, 8, 47, 227, 187,
	80, 105, 148, 232, 205, 214, 99, 208, 19, 22, 193, 97, 173, 229, 211, 238,
	41, 94, 224, 25, 130, 233, 200, 86, 17, 63, 170, 93, 55, 12, 83, 73,
	64, 118, 52, 121, 36, 183, 79, 111, 177, 108, 154, 92, 228, 140, 203, 2,
	109, 168, 58, 28, 103, 240, 37, 165, 250, 217, 226, 127, 231, 51, 20, 245,
	96, 138, 234, 45, 199, 66, 67, 196, 150, 87, 3, 174, 215, 132, 90, 75,
	135, 98, 239, 142, 30, 33, 133, 204, 251, 185, 88, 117, 39, 69, 252, 48,
	146, 24, 186, 126, 172, 141, 50, 188, 131, 149, 194, 44, 255, 243, 143, 210,
	181, 102, 254, 237, 21, 191, 56, 15, 4, 71, 184, 216, 222, 49, 242, 236
};


//// Mul/Div tables

u8 * CAT_RESTRICT cat::GF256_MUL_TABLE = 0;
u8 * CAT_RESTRICT cat::GF256_DIV_TABLE = 0;


// Unpack 256x256 multiplication tables
void cat::GF256Init() {
	// If initialized already,
	if (GF256_MUL_TABLE) {
		return;
	}

	// Allocate table memory 65KB x 2
	GF256_MUL_TABLE = (u8 *)malloc(256 * 256 * 2);
	GF256_DIV_TABLE = GF256_MUL_TABLE + 256 * 256;
	u8 *m = GF256_MUL_TABLE, *d = GF256_DIV_TABLE;

	// Unroll y = 0 subtable
	for (int x = 0; x < 256; ++x) {
		m[x] = d[x] = 0;
	}

	// For each other y value,
	for (int y = 1; y < 256; ++y) {
		// Calculate log(y) for mult and 255 - log(y) for div
		const u8 log_y = GF256_LOG_TABLE[y];
		const u8 log_yn = 255 - log_y;

		// Next subtable
		m += 256;
		d += 256;

		// Unroll x = 0
		m[0] = 0;
		d[0] = 0;

		// Calculate x * y, x / y
		for (int x = 1; x < 256; ++x) {
			int log_x = GF256_LOG_TABLE[x];

			m[x] = GF256_EXP_TABLE[log_x + log_y];
			d[x] = GF256_EXP_TABLE[log_x + log_yn];
		}
	}
}

// Performs "dest[] += src[] * x" operation in GF(256)
void cat::GF256MemMulAdd(void * CAT_RESTRICT vdest, u8 x, const void * CAT_RESTRICT vsrc, int bytes) {
	if (x == 0) {
		return;
	}

	if (x == 1) {
		memxor(vdest, vsrc, bytes);
		return;
	}

	u8 * CAT_RESTRICT dest = reinterpret_cast<u8*>( vdest );
	const u8 * CAT_RESTRICT src = reinterpret_cast<const u8*>( vsrc );
	const u8 * CAT_RESTRICT table = GF256_MUL_TABLE + ((u32)x << 8);

	// For each block of 8 bytes,
	while (bytes >= 8) {
#ifdef CAT_ENDIAN_LITTLE
		// This optimization works because it reduces the number of memory
		// accesses on desktops by almost half.
		u64 x = table[src[0]];
		x |= (u64)table[src[1]] << 8;
		x |= (u64)table[src[2]] << 16;
		x |= (u64)table[src[3]] << 24;
		x |= (u64)table[src[4]] << 32;
		x |= (u64)table[src[5]] << 40;
		x |= (u64)table[src[6]] << 48;
		x |= (u64)table[src[7]] << 56;
		*(u64*)dest ^= x;
#else
		dest[0] ^= table[src[0]];
		dest[1] ^= table[src[1]];
		dest[2] ^= table[src[2]];
		dest[3] ^= table[src[3]];
		dest[4] ^= table[src[4]];
		dest[5] ^= table[src[5]];
		dest[6] ^= table[src[6]];
		dest[7] ^= table[src[7]];
#endif // CAT_ENDIAN_LITTLE

		src += 8;
		dest += 8;
		bytes -= 8;
	}

	// For each byte,
	while (bytes-- > 0) {
		// Multiply source byte by x and add it to destination byte
		*dest++ ^= table[*src++];
	}
}

// Performs "dest[] /= x" operation in GF(256)
void cat::GF256MemDivide(void * CAT_RESTRICT vdest, u8 x, int bytes) {
	u8 * CAT_RESTRICT dest = reinterpret_cast<u8*>( vdest );
	const u8 * CAT_RESTRICT table = GF256_DIV_TABLE + ((u32)x << 8);

	// For each block of 8 bytes,
	while (bytes >= 8) {
#ifdef CAT_ENDIAN_LITTLE
		// This optimization works because it reduces the number of memory
		// accesses on desktops by almost half.
		u64 x = table[dest[0]];
		x |= (u64)table[dest[1]] << 8;
		x |= (u64)table[dest[2]] << 16;
		x |= (u64)table[dest[3]] << 24;
		x |= (u64)table[dest[4]] << 32;
		x |= (u64)table[dest[5]] << 40;
		x |= (u64)table[dest[6]] << 48;
		x |= (u64)table[dest[7]] << 56;
		*(u64*)dest = x;
#else
		dest[0] = table[dest[0]];
		dest[1] = table[dest[1]];
		dest[2] = table[dest[2]];
		dest[3] = table[dest[3]];
		dest[4] = table[dest[4]];
		dest[5] = table[dest[5]];
		dest[6] = table[dest[6]];
		dest[7] = table[dest[7]];
#endif // CAT_ENDIAN_LITTLE

		dest += 8;
		bytes -= 8;
	}

	// For each byte,
	while (bytes-- > 0) {
		// Multiply source byte by x and add it to destination byte
		*dest = table[*dest];
		++dest;
	}
}

