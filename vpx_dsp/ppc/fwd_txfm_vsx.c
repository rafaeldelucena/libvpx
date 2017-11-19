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
#include "vpx_dsp/ppc/transpose_vsx.h"
#include "vpx_dsp/ppc/bitdepth_conversion_vsx.h"

static INLINE int32x4_t fdct_vector_round_shift(int32x4_t input) {
  uint32x4_t dct_const_0, dct_const_1;
  int32x4_t tmp, sum;
  int32x4_t one = vec_splat_s32(1);
  dct_const_0 = vec_splat_u32(DCT_CONST_BITS - 1);
  dct_const_1 = vec_splat_u32(DCT_CONST_BITS);
  tmp = vec_sl(one, dct_const_0);
  sum = vec_add(input, tmp);
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

  int16_t intermediate[16];

  // Do the two transform/transpose passes
  uint8x16_t perm0 = {0x8, 0x9, 0xA, 0xB, 0xC, 0xD, 0xE, 0xF, 0x0, 0x1, 0x2, 0x3, 0x4, 0x5, 0x6, 0x7};
  int16x8_t cospi_0 = {cospi_16_64, cospi_16_64, cospi_16_64, cospi_16_64, cospi_16_64, cospi_16_64, cospi_16_64, cospi_16_64};
  int16x8_t cospi_1 = {cospi_16_64, cospi_16_64, cospi_16_64, cospi_16_64, -cospi_16_64, -cospi_16_64, -cospi_16_64, -cospi_16_64};
  int16x8_t cospi_2 = {cospi_24_64, cospi_24_64, cospi_24_64, cospi_24_64, cospi_8_64, cospi_8_64, cospi_8_64, cospi_8_64};
  int16x8_t cospi_3 = {-cospi_8_64, -cospi_8_64, -cospi_8_64, -cospi_8_64, cospi_24_64, cospi_24_64, cospi_24_64, cospi_24_64};

  uint16x8_t four = vec_splat_u16(4);

  int16x8_t loaded1, loaded0, loaded2;
  loaded1 = vec_vsx_ld(0, input);
  loaded0 = vec_vsx_ld(0, input + (2 * stride));
  loaded2 = vec_perm(loaded0, loaded0, perm0);

  int16x8_t in_high0, in_high1;
  in_high0 = vec_sl(loaded1, four);
  in_high1 = vec_sl(loaded2, four);

  if (in_high0[0]) {
    ++in_high0[0];
  }

  int16x8_t step[3];
  step[0] = vec_perm(step[2], step[2], perm0);
  step[1] = vec_add(in_high0, in_high1);
  step[2] = vec_sub(in_high0, in_high1);

  int32x4_t e_0, o_0, h_0, l_0;
  e_0 = vec_mule(step[1], cospi_0);
  o_0 = vec_mulo(step[1], cospi_0);
  h_0 = vec_mergeh(e_0, o_0);
  l_0 = vec_mergel(e_0, o_0);

  int32x4_t e_1, o_1, h_1, l_1;
  e_1 = vec_mule(step[1], cospi_1);
  o_1 = vec_mulo(step[1], cospi_1);
  h_1 = vec_mergeh(e_1, o_1);
  l_1 = vec_mergel(e_1, o_1);

  int32x4_t e_2, o_2, h_2, l_2;
  e_2 = vec_mule(step[0], cospi_2);
  o_2 = vec_mulo(step[0], cospi_2);
  h_2 = vec_mergeh(e_2, o_2);
  l_2 = vec_mergel(e_2, o_2);

  int32x4_t e_3, o_3, h_3, l_3;
  e_3 = vec_mule(step[0], cospi_3);
  o_3 = vec_mulo(step[0], cospi_3);
  h_3 = vec_mergeh(e_3, o_3);
  l_3 = vec_mergel(e_3, o_3);

  int32x4_t v[4];
  v[0] = vec_add(h_0, l_0);
  v[2] = vec_add(h_1, l_1);
  v[1] = vec_add(h_2, l_2);
  v[3] = vec_add(h_3, l_3);

  vpx_transpose_s32_4x4(v);

  int32x4_t tmp[4];
  tmp[0] = fdct_vector_round_shift(v[0]);
  tmp[1] = fdct_vector_round_shift(v[1]);
  tmp[2] = fdct_vector_round_shift(v[2]);
  tmp[3] = fdct_vector_round_shift(v[3]);

  int16x8_t out[2];
#ifdef WORDS_BIGENDIAN
  out[0] = vec_pack(tmp[1], tmp[0]);
  out[1] = vec_pack(tmp[3], tmp[2]);
#else
  out[0] = vec_pack(tmp[0], tmp[1]);
  out[1] = vec_pack(tmp[2], tmp[3]);
#endif // WORDS_BIGENDIAN

  vec_vsx_st(out[0], 0, intermediate);
  vec_vsx_st(out[1], 0, intermediate + (2* stride));

  // Do next column (which is a transposed row in second/horizontal pass)
  // Transform.

  uint8x16_t perm1 = {0x0, 0x1, 0x18, 0x19, 0x0, 0x1, 0x18, 0x19, 0x2, 0x3, 0x1A, 0x1B, 0x2, 0x3, 0x1A, 0x1B};
  uint8x16_t perm2 = {0x8, 0x9, 0x10, 0x11, 0x8, 0x9, 0x10, 0x11, 0xA, 0xB, 0x12, 0x13, 0xA, 0xB, 0x12, 0x13};
  uint8x16_t perm3 = {0x4, 0x5, 0x1C, 0x1D, 0x4, 0x5, 0x1C, 0x1D, 0x6, 0x7, 0x1E, 0x1F, 0x6, 0x7, 0x1E, 0x1F};
  uint8x16_t perm4 = {0xC, 0xD, 0x14, 0x15, 0xC, 0xD, 0x14, 0x15, 0xE, 0xF, 0x16, 0x17, 0xE, 0xF, 0x16, 0x17};

  out[1] = vec_perm(out[1], out[1], perm0);

  int16x8_t step1, step2;
  step1 = vec_add(out[0], out[1]);
  step2 = vec_sub(out[0], out[1]);

  int16x8_t cospi_0_0 = {cospi_16_64, cospi_24_64, cospi_16_64, -cospi_8_64, cospi_16_64, cospi_24_64, cospi_16_64, -cospi_8_64};
  int16x8_t cospi_0_1 = {cospi_16_64, cospi_8_64, -cospi_16_64, cospi_24_64, cospi_16_64, cospi_8_64, -cospi_16_64, cospi_24_64};
  int16x8_t cospi_0_2 = {cospi_16_64, cospi_24_64, cospi_16_64, -cospi_8_64, cospi_16_64, cospi_24_64, cospi_16_64, -cospi_8_64};
  int16x8_t cospi_0_3 = {cospi_16_64, cospi_8_64, -cospi_16_64, cospi_24_64, cospi_16_64, cospi_8_64, -cospi_16_64, cospi_24_64};

  int16x8_t x_0_0, x_0_1, x_0_2, x_0_3;
  x_0_0 = vec_perm(step1, step2, perm1);
  x_0_1 = vec_perm(step1, step2, perm2);
  x_0_2 = vec_perm(step1, step2, perm3);
  x_0_3 = vec_perm(step1, step2, perm4);

  int32x4_t e_0_0, o_0_0, h_0_0, l_0_0;
  e_0_0 = vec_mule(x_0_0, cospi_0_0);
  o_0_0 = vec_mulo(x_0_0, cospi_0_0);
  h_0_0 = vec_mergeh(e_0_0, o_0_0);
  l_0_0 = vec_mergel(e_0_0, o_0_0);

  int32x4_t e_0_1, o_0_1, h_0_1, l_0_1;
  e_0_1 = vec_mule(x_0_1, cospi_0_1);
  o_0_1 = vec_mulo(x_0_1, cospi_0_1);
  h_0_1 = vec_mergeh(e_0_1, o_0_1);
  l_0_1 = vec_mergel(e_0_1, o_0_1);

  int32x4_t e_0_2, o_0_2, h_0_2, l_0_2;
  e_0_2 = vec_mule(x_0_2, cospi_0_2);
  o_0_2 = vec_mulo(x_0_2, cospi_0_2);
  h_0_2 = vec_mergeh(e_0_2, o_0_2);
  l_0_2 = vec_mergel(e_0_2, o_0_2);

  int32x4_t e_0_3, o_0_3, h_0_3, l_0_3;
  e_0_3 = vec_mule(x_0_3, cospi_0_3);
  o_0_3 = vec_mulo(x_0_3, cospi_0_3);
  h_0_3 = vec_mergeh(e_0_3, o_0_3);
  l_0_3 = vec_mergel(e_0_3, o_0_3);

  int32x4_t tmp0, tmp1, tmp2, tmp3;
  tmp0 = vec_add(h_0_0, h_0_1);
  tmp1 = vec_add(l_0_0, l_0_1);
  tmp2 = vec_add(h_0_2, h_0_3);
  tmp3 = vec_add(l_0_2, l_0_3);

  int32x4_t one = vec_splat_s32(1);
  uint32x4_t two = vec_splat_u32(2);

  int32x4_t temp0, temp1, temp2, temp3;
  temp0 = fdct_vector_round_shift(tmp0);
  temp1 = fdct_vector_round_shift(tmp1);
  temp2 = fdct_vector_round_shift(tmp2);
  temp3 = fdct_vector_round_shift(tmp3);

  int32x4_t temp4, temp5, temp6, temp7;
  temp4 = vec_add(temp0, one);
  temp5 = vec_add(temp1, one);
  temp6 = vec_add(temp2, one);
  temp7 = vec_add(temp3, one);

  int32x4_t temp8, temp9, tempA, tempB;
  temp8 = vec_sra(temp4, two);
  temp9 = vec_sra(temp5, two);
  tempA = vec_sra(temp6, two);
  tempB = vec_sra(temp7, two);

  int16x8_t out1, out2;
#ifdef WORDS_BIGENDIAN
  out1 = vec_pack(temp9, temp8);
  out2 = vec_pack(tempB, tempA);
#else
  out1 = vec_pack(temp8, temp9);
  out2 = vec_pack(tempA, tempB);
#endif // WORDS_BIGENDIAN

  vec_vsx_st(out1, 0, output);
  vec_vsx_st(out2, 0, output + (2* stride));
}
