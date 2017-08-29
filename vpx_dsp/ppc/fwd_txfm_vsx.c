/*
 *  Copyright (c) 2017 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <assert.h>
#include <stdio.h>
#include "./vpx_dsp_rtcd.h"
#include "vpx_dsp/fwd_txfm.h"
#include "vpx_dsp/ppc/types_vsx.h"
#include "vpx_dsp/ppc/bitdepth_conversion_vsx.h"

#define MATRIX_H4_PRINT(a0)\
    printf(\
    "0[ %04hx %04hx %04hx %04hx ]\n"\
    "1[ %04hx %04hx %04hx %04hx ]\n"\
    "2[ %04hx %04hx %04hx %04hx ]\n"\
    "3[ %04hx %04hx %04hx %04hx ]\n",\
     a0[0*4 +0], a0[0*4 +1], a0[0*4 +2], a0[0*4 +3],\
     a0[1*4 +0], a0[1*4 +1], a0[1*4 +2], a0[1*4 +3],\
     a0[2*4 +0], a0[2*4 +1], a0[2*4 +2], a0[2*4 +3],\
     a0[3*4 +0], a0[3*4 +1], a0[3*4 +2], a0[3*4 +3])

#define MATRIX_SS4_PRINT(a0,a1,a2,a3)\
    printf(\
    "0[ %04hx %04hx %04hx %04hx ]\n"\
    "1[ %04hx %04hx %04hx %04hx ]\n"\
    "2[ %04hx %04hx %04hx %04hx ]\n"\
    "3[ %04hx %04hx %04hx %04hx ]\n",\
     a0[0], a0[1], a0[2], a0[3],\
     a1[0], a1[1], a1[2], a1[3],\
     a2[0], a2[1], a2[2], a2[3],\
     a3[0], a3[1], a3[2], a3[3])

#define MATRIX_SI4_PRINT(a0,a1)\
    printf(\
    "0[ %04hx %04hx %04hx %04hx ]\n"\
    "1[ %04hx %04hx %04hx %04hx ]\n"\
    "2[ %04hx %04hx %04hx %04hx ]\n"\
    "3[ %04hx %04hx %04hx %04hx ]\n",\
     a0[0], a0[1], a0[2], a0[3],\
     a0[4], a0[5], a0[6], a0[7],\
     a1[0], a1[1], a1[2], a1[3],\
     a1[4], a1[5], a1[6], a1[7])

#define MATRIX_I4_PRINT(a0)\
    printf(\
    "0[ %d %d %d %d ]\n"\
    "1[ %d %d %d %d ]\n"\
    "2[ %d %d %d %d ]\n"\
    "3[ %d %d %d %d ]\n",\
     a0[0*4 +0], a0[0*4 +1], a0[0*4 +2], a0[0*4 +3],\
     a0[1*4 +0], a0[1*4 +1], a0[1*4 +2], a0[1*4 +3],\
     a0[2*4 +0], a0[2*4 +1], a0[2*4 +2], a0[2*4 +3],\
     a0[3*4 +0], a0[3*4 +1], a0[3*4 +2], a0[3*4 +3])

static INLINE int32x4_t fdct_vector_round_shift(int32x4_t input) {
  int32x4_t one = vec_splat_s32(1);
  uint32x4_t dct_const_0 = vec_splat_u32(DCT_CONST_BITS - 1);
  uint32x4_t dct_const_1 = vec_splat_u32(DCT_CONST_BITS);
  int32x4_t tmp = vec_sl(one, dct_const_0);
  int32x4_t sum = vec_add(input, tmp);
  return vec_sra(sum, dct_const_1);
}

void vpx_fdct4x4_vsx(const int16_t *input, tran_low_t *output, int stride) {
  // The 2D transform is done with two passes which are actually pretty
  // similar. In the first one, we transform the columns and transpose
  // the results. In the second one, we transform the rows. To achieve that,
  // as the first pass results are transposed, we transpose the columns (that
  // is the transposed rows) and transpose the results (so that it goes back
  // in normal/row positions).
  // We need an intermediate buffer between passes.

  tran_low_t intermediate[16];
  printf("---------------------------------------- TRANSFORM INPUT\n");
  MATRIX_H4_PRINT(input);
  // Do the two transform/transpose passes
  {
    int16x8_t in1 = vec_vsx_ld(0, input);
    int16x8_t in2 = vec_vsx_ld(0, input + (2 * stride));

    if (in1[0]) {
      ++in1[0];
    }

    int16x8_t step1 = vec_add(in1, in2);
    int16x8_t step2 = vec_sub(in1, in2);

    uint8x16_t perm1 = {0x0, 0x1, 0x18, 0x19, 0x0, 0x1, 0x18, 0x19, 0x2, 0x3, 0x1A, 0x1B, 0x2, 0x3, 0x1A, 0x1B};
    int16x8_t x_0_0 = vec_perm(step1, step2, perm1);
    int16x8_t cospi_0_0 = {cospi_16_64, cospi_24_64, cospi_16_64, -cospi_8_64, cospi_16_64, cospi_24_64, cospi_16_64, -cospi_8_64};

    int32x4_t e_0_0 = vec_mule(x_0_0, cospi_0_0);
    int32x4_t o_0_0 = vec_mulo(x_0_0, cospi_0_0);
    int32x4_t h_0_0 = vec_mergeh(e_0_0, o_0_0);
    int32x4_t l_0_0 = vec_mergel(e_0_0, o_0_0);

    uint8x16_t perm2 = {0x8, 0x9, 0x10, 0x11, 0x8, 0x9, 0x10, 0x11, 0xA, 0xB, 0x12, 0x13, 0xA, 0xB, 0x12, 0x13};
    int16x8_t x_0_1 = vec_perm(step1, step2, perm2);
    int16x8_t cospi_0_1 = {cospi_16_64, cospi_8_64, -cospi_16_64, cospi_24_64, cospi_16_64, cospi_8_64, -cospi_16_64, cospi_24_64};

    int32x4_t e_0_1 = vec_mule(x_0_1, cospi_0_1);
    int32x4_t o_0_1 = vec_mulo(x_0_1, cospi_0_1);
    int32x4_t h_0_1 = vec_mergeh(e_0_1, o_0_1);
    int32x4_t l_0_1 = vec_mergel(e_0_1, o_0_1);

    uint8x16_t perm3 = {0x4, 0x5, 0x1C, 0x1D, 0x4, 0x5, 0x1C, 0x1D, 0x6, 0x7, 0x1E, 0x1F, 0x6, 0x7, 0x1E, 0x1F};
    int16x8_t x_0_2 = vec_perm(step1, step2, perm3);
    int16x8_t cospi_0_2 = {cospi_16_64, cospi_24_64, cospi_16_64, -cospi_8_64, cospi_16_64, cospi_24_64, cospi_16_64, -cospi_8_64};

    int32x4_t e_0_2 = vec_mule(x_0_2, cospi_0_2);
    int32x4_t o_0_2 = vec_mulo(x_0_2, cospi_0_2);
    int32x4_t h_0_2 = vec_mergeh(e_0_2, o_0_2);
    int32x4_t l_0_2 = vec_mergel(e_0_2, o_0_2);

    uint8x16_t perm4 = {0xC, 0xD, 0x14, 0x15, 0xC, 0xD, 0x14, 0x15, 0xE, 0xF, 0x16, 0x17, 0xE, 0xF, 0x16, 0x17};
    int16x8_t x_0_3 = vec_perm(step1, step2, perm4);
    int16x8_t cospi_0_3 = {cospi_16_64, cospi_8_64, -cospi_16_64, cospi_24_64, cospi_16_64, cospi_8_64, -cospi_16_64, cospi_24_64};

    int32x4_t e_0_3 = vec_mule(x_0_3, cospi_0_3);
    int32x4_t o_0_3 = vec_mulo(x_0_3, cospi_0_3);
    int32x4_t h_0_3 = vec_mergeh(e_0_3, o_0_3);
    int32x4_t l_0_3 = vec_mergel(e_0_3, o_0_3);

    int32x4_t v[4];
    v[0] = vec_add(h_0_0, h_0_1);
    v[1] = vec_add(l_0_0, l_0_1);
    v[2] = vec_add(h_0_2, h_0_3);
    v[3] = vec_add(l_0_2, l_0_3);

    vpx_transpose_s32_4x4(v[4]);

#ifdef WORDS_BIGENDIAN
    int16x8_t out1 = vec_pack(temp1, temp0);
    int16x8_t out2 = vec_pack(temp3, temp2);
#else
    int16x8_t out1 = vec_pack(temp0, temp1);
    int16x8_t out2 = vec_pack(temp2, temp3);
#endif // WORDS_BIGENDIAN

    store_tran_low(out1, 0, intermediate);
    store_tran_low(out2, 0, intermediate + (2* stride));

  // Do next column (which is a transposed row in second/horizontal pass)
  }
  // Transform.

  printf("---------------------------------------- TRANSFORM INTERMEDIATE\n");
  MATRIX_H4_PRINT(intermediate);

  int16x8_t in1 = vec_vsx_ld(0, intermediate);
  int16x8_t in0 = vec_vsx_ld(0, intermediate + (2 * stride));
  uint8x16_t perm0 = {0x8, 0x9, 0xA, 0xB, 0xC, 0xD, 0xE, 0xF, 0x0, 0x1, 0x2, 0x3, 0x4, 0x5, 0x6, 0x7};

  int16x8_t in2 = vec_perm(in0, in0, perm0);

  int16x8_t step1 = vec_add(in1, in2);
  int16x8_t step2 = vec_sub(in1, in2);

  uint8x16_t perm1 = {0x0, 0x1, 0x18, 0x19, 0x0, 0x1, 0x18, 0x19, 0x2, 0x3, 0x1A, 0x1B, 0x2, 0x3, 0x1A, 0x1B};
  int16x8_t x_0_0 = vec_perm(step1, step2, perm1);
  int16x8_t cospi_0_0 = {cospi_16_64, cospi_24_64, cospi_16_64, -cospi_8_64, cospi_16_64, cospi_24_64, cospi_16_64, -cospi_8_64};

  int32x4_t e_0_0 = vec_mule(x_0_0, cospi_0_0);
  int32x4_t o_0_0 = vec_mulo(x_0_0, cospi_0_0);
  int32x4_t h_0_0 = vec_mergeh(e_0_0, o_0_0);
  int32x4_t l_0_0 = vec_mergel(e_0_0, o_0_0);

  uint8x16_t perm2 = {0x8, 0x9, 0x10, 0x11, 0x8, 0x9, 0x10, 0x11, 0xA, 0xB, 0x12, 0x13, 0xA, 0xB, 0x12, 0x13};
  int16x8_t x_0_1 = vec_perm(step1, step2, perm2);
  int16x8_t cospi_0_1 = {cospi_16_64, cospi_8_64, -cospi_16_64, cospi_24_64, cospi_16_64, cospi_8_64, -cospi_16_64, cospi_24_64};

  int32x4_t e_0_1 = vec_mule(x_0_1, cospi_0_1);
  int32x4_t o_0_1 = vec_mulo(x_0_1, cospi_0_1);
  int32x4_t h_0_1 = vec_mergeh(e_0_1, o_0_1);
  int32x4_t l_0_1 = vec_mergel(e_0_1, o_0_1);

  uint8x16_t perm3 = {0x4, 0x5, 0x1C, 0x1D, 0x4, 0x5, 0x1C, 0x1D, 0x6, 0x7, 0x1E, 0x1F, 0x6, 0x7, 0x1E, 0x1F};
  int16x8_t x_0_2 = vec_perm(step1, step2, perm3);
  int16x8_t cospi_0_2 = {cospi_16_64, cospi_24_64, cospi_16_64, -cospi_8_64, cospi_16_64, cospi_24_64, cospi_16_64, -cospi_8_64};

  int32x4_t e_0_2 = vec_mule(x_0_2, cospi_0_2);
  int32x4_t o_0_2 = vec_mulo(x_0_2, cospi_0_2);
  int32x4_t h_0_2 = vec_mergeh(e_0_2, o_0_2);
  int32x4_t l_0_2 = vec_mergel(e_0_2, o_0_2);

  uint8x16_t perm4 = {0xC, 0xD, 0x14, 0x15, 0xC, 0xD, 0x14, 0x15, 0xE, 0xF, 0x16, 0x17, 0xE, 0xF, 0x16, 0x17};
  int16x8_t x_0_3 = vec_perm(step1, step2, perm4);
  int16x8_t cospi_0_3 = {cospi_16_64, cospi_8_64, -cospi_16_64, cospi_24_64, cospi_16_64, cospi_8_64, -cospi_16_64, cospi_24_64};

  int32x4_t e_0_3 = vec_mule(x_0_3, cospi_0_3);
  int32x4_t o_0_3 = vec_mulo(x_0_3, cospi_0_3);
  int32x4_t h_0_3 = vec_mergeh(e_0_3, o_0_3);
  int32x4_t l_0_3 = vec_mergel(e_0_3, o_0_3);

  int32x4_t tmp0 = vec_add(h_0_0, h_0_1);
  int32x4_t tmp1 = vec_add(l_0_0, l_0_1);
  int32x4_t tmp2 = vec_add(h_0_2, h_0_3);
  int32x4_t tmp3 = vec_add(l_0_2, l_0_3);

  int32x4_t temp0 = fdct_vector_round_shift(tmp0);
  int32x4_t temp1 = fdct_vector_round_shift(tmp1);
  int32x4_t temp2 = fdct_vector_round_shift(tmp2);
  int32x4_t temp3 = fdct_vector_round_shift(tmp3);

  int32x4_t one = vec_splat_s32(1);
  uint32x4_t two = vec_splat_u32(2);

  int32x4_t temp4 = vec_add(temp0, one);
  int32x4_t temp5 = vec_add(temp1, one);
  int32x4_t temp6 = vec_add(temp2, one);
  int32x4_t temp7 = vec_add(temp3, one);

  int32x4_t temp8 = vec_sra(temp4, two);
  int32x4_t temp9 = vec_sra(temp5, two);
  int32x4_t tempA = vec_sra(temp6, two);
  int32x4_t tempB = vec_sra(temp7, two);

#ifdef WORDS_BIGENDIAN
  int16x8_t out1 = vec_pack(temp9, temp8);
  int16x8_t out2 = vec_pack(tempB, tempA);
#else
  int16x8_t out1 = vec_pack(temp8, temp9);
  int16x8_t out2 = vec_pack(tempA, tempB);
#endif // WORDS_BIGENDIAN
  printf("---------------------------------------- TRANSFORM OUTPUT\n");
  MATRIX_H4_PRINT(output);

  store_tran_low(out1, 0, output);
  store_tran_low(out2, 0, output + (2* stride));
}
