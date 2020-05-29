/*******************************************************************************
* mrg32k3a_jump.h: this file is part of the randms library.
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

#ifndef __MRG32K3A_JUMP_H__
#define __MRG32K3A_JUMP_H__

#include <stdint.h>

/* Pre-computed matrices for jumping ahead.
 * ref: https://doi.org/10.1016/B978-0-12-384988-5.00016-4
 *
 * The matrices for skip length (g * 8^i) are computed:
 *   A_k^(g * 8^i) mod m_k
 * with k in {1,2}, g from 0 to (MRG32K3A_MAX_STEP - 1), and i from 1 to 7.
 */
#define MRG32K3A_MAX_STEP_B8    21      /* maximum skip length (in log8) */
#define MRG32K3A_MAX_STEP       0x7fffffffffffffffULL   /* max skip length */

extern const uint64_t A1gj8i[MRG32K3A_MAX_STEP_B8][7][9];
extern const uint64_t A2gj8i[MRG32K3A_MAX_STEP_B8][7][9];

#endif

