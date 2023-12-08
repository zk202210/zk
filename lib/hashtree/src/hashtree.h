/*
MIT License

Copyright (c) 2021 Prysmatic Labs

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/
#ifndef HASHTREE_H
#define HASHTREE_H
#include <stdint.h>

#ifdef __aarch64__
void sha256_armv8_neon_x1(unsigned char* output, const unsigned char* input, uint64_t count);
void sha256_armv8_neon_x4(unsigned char* output, const unsigned char* input, uint64_t count);
#endif 

#ifdef __x86_64__
void sha256_1_sse(unsigned char* output, const unsigned char* input, uint64_t count);
void sha256_1_avx(unsigned char* output, const unsigned char* input, uint64_t count);
void sha256_4_avx(unsigned char* output, const unsigned char* input, uint64_t count);
void sha256_8_avx2(unsigned char* output, const unsigned char* input, uint64_t count);
void sha256_16_avx512(unsigned char* output, const unsigned char* input, uint64_t count);
void sha256_shani(unsigned char* output, const unsigned char* input, uint64_t count);

void sha256_multiblk_shani(unsigned char* output, const unsigned char* input, uint64_t count, uint64_t leaf_size);
void sha256_multiblk_16_avx512(unsigned char* output, const unsigned char* input, uint64_t count, uint64_t leaf_size);

void sha256_hash_layer(unsigned char* output, const unsigned char* input, uint64_t count);
void sha256_hash_leaves(unsigned char* output, const unsigned char* input, uint64_t count, uint64_t leaf_size);
#endif
#endif
