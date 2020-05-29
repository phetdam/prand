/*******************************************************************************
* randms.h: this file is part of the randms library.
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

#ifndef __RANDMS_H__
#define __RANDMS_H__

#include <stdint.h>
#include <stdlib.h>

/******************************************************************************
  A list of the implemented random number generators.
******************************************************************************/
#define RANDMS_MRG32K3A         1
#define RANDMS_MT19937          2

/******************************************************************************
  Error codes.
******************************************************************************/
#define RANDMS_ERR_MEMORY               (-1)
#define RANDMS_ERR_MEMORY_JUMP          (-2)
#define RANDMS_ERR_STEP                 (-3)
#define RANDMS_ERR_UNDEF_RNG            (-4)
#define RANDMS_WARN_SEED                1

#define RANDMS_IS_ERROR(err)            ((err) < 0)
#define RANDMS_IS_WARN(err)             ((err) > 0)

/******************************************************************************
  LCG for initialising the states for the MRG32k3a generator.
******************************************************************************/
#define RANDMS_LCG(n) ((69069 * n + 1) & 0xffffffffUL)

/******************************************************************************
  The universal interface of the random number generators.
******************************************************************************/
typedef struct randms_struct {
  void *state;                  /* the state for single stream */
  void **state_stream;          /* states for multiple streams */
  int nstream;                  /* number of random streams */
  int type;                     /* type of the random number generator */
  int64_t min;                  /* minimum value of the random integer */
  int64_t max;                  /* maximum value of the random integer */
  /* function pointers for sampling numbers */
  uint64_t (*get) (void *);
  double (*get_double) (void *);
  /* function pointers for reseting states with seed and skipping steps */
  void (*reset) (void *, const uint64_t, const uint64_t, int *);
  void (*reset_all) (struct randms_struct *, const uint64_t, const uint64_t,
      int *);
  /* function pointers for jumping ahead */
  void (*jump) (void *, const uint64_t, int *);
  void (*jump_all) (struct randms_struct *, const uint64_t, int *);
} randms_t;

/******************************************************************************
  Random number generator initialisers.
******************************************************************************/
randms_t *mrg32k3a_init(const uint64_t, const unsigned int, const uint64_t,
    int *);
randms_t *mt19937_init(const uint64_t, const unsigned int, const uint64_t,
    int *);

randms_t *randms_init(const int, const uint64_t, const unsigned int,
    const uint64_t, int *);

/******************************************************************************
  Error message and clean-up functions.
******************************************************************************/
char *randms_errmsg(const int);

void randms_destroy(randms_t *);

#endif

