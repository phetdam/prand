/*******************************************************************************
* mt19937_jump.h: this file is part of the randms library.
* 
* randms: C library for generating random numbers with multiple streams.
*
* Github repository:
*       https://github.com/cheng-zhao/randms
*
* Copyright (c) 2020 Cheng Zhao <zhaocheng03@gmail.com>
* 
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
* 
* The above copyright notice and this permission notice shall be included in all
* copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
* SOFTWARE.
*
*******************************************************************************/

#ifndef __MT19927_JUMP_H__
#define __MT19927_JUMP_H__

#include <stdint.h>
#include <string.h>

/* Definitions for the representation of polynomials with 32-bit words. */
#define MT19937_N       624
#define WORD_SIZE       32
#define MUL_NBIT(x)     ((x) << 5)
#define DIV_NBIT(x)     ((x) >> 5)
#define MOD_NBIT(x)     ((x) & 0x1fU)
#define COEF(x,i)       (((x[DIV_NBIT(i)]) >> MOD_NBIT(i)) & 1)

/* Pre-computed polynomials for jumping ahead.
 * ref: https://doi.org/10.1007/978-3-540-85912-3_26
 *
 * The polynomials for skip length (g * 8^i) are computed:
 *   p = t^(g * 8^i) mod phi
 * with g from 0 to (MT19937_MAX_STEP - 1), and i from 1 to 7.
 * The idea of base-8 skip length is taken from
 * https://doi.org/10.1016/B978-0-12-384988-5.00016-4
 */
#define MT19937_MAX_STEP_B8     21      /* maximum skip length (in log8) */
#define MT19937_MAX_STEP        0x7fffffffffffffffULL   /* max skip length */
extern const uint32_t mt19937_poly[MT19937_MAX_STEP_B8][7][MT19937_N];


/* Low-level functions for polynomial multiplication and modular reduction,
 * see mt19937_poly.c for details. */

void poly_mul(uint32_t *, const uint32_t *, const uint32_t *,
    const unsigned int, uint32_t *);

void poly_mul_ub(uint32_t *, const uint32_t *, const uint32_t *,
    const unsigned int, uint32_t *);

/* This modular reduction function is only valid with the divisor being
 * the minimal polynomial of MT19937, so the divisor is not needed here. */
void poly_mod_phi(uint32_t *, uint32_t *);

#endif

