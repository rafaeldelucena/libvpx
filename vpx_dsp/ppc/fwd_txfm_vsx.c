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

  // Do the two transform/transpose passes
  int16x8_t loaded[2], in_high[2], step[2], out[2], x[4];

  int32x4_t e[4], o[4], h[4], l[4], v[4], tmp[4];

  uint16x8_t four = vec_splat_u16(4);
  int32x4_t one = vec_splat_s32(1);
  uint32x4_t two = vec_splat_u32(2);

  uint8x16_t perm0 = {0x8, 0x9, 0xA, 0xB, 0xC, 0xD, 0xE, 0xF, 0x0, 0x1, 0x2, 0x3, 0x4, 0x5, 0x6, 0x7};
  uint8x16_t perm1 = {0x0, 0x1, 0x18, 0x19, 0x0, 0x1, 0x18, 0x19, 0x2, 0x3, 0x1A, 0x1B, 0x2, 0x3, 0x1A, 0x1B};
  uint8x16_t perm2 = {0x8, 0x9, 0x10, 0x11, 0x8, 0x9, 0x10, 0x11, 0xA, 0xB, 0x12, 0x13, 0xA, 0xB, 0x12, 0x13};
  uint8x16_t perm3 = {0x4, 0x5, 0x1C, 0x1D, 0x4, 0x5, 0x1C, 0x1D, 0x6, 0x7, 0x1E, 0x1F, 0x6, 0x7, 0x1E, 0x1F};
  uint8x16_t perm4 = {0xC, 0xD, 0x14, 0x15, 0xC, 0xD, 0x14, 0x15, 0xE, 0xF, 0x16, 0x17, 0xE, 0xF, 0x16, 0x17};

  int16x8_t cospi_0 = {cospi_16_64, cospi_16_64, cospi_16_64, cospi_16_64, cospi_16_64, cospi_16_64, cospi_16_64, cospi_16_64};
  int16x8_t cospi_1 = {cospi_16_64, cospi_16_64, cospi_16_64, cospi_16_64, -cospi_16_64, -cospi_16_64, -cospi_16_64, -cospi_16_64};
  int16x8_t cospi_2 = {cospi_24_64, cospi_24_64, cospi_24_64, cospi_24_64, cospi_8_64, cospi_8_64, cospi_8_64, cospi_8_64};
  int16x8_t cospi_3 = {-cospi_8_64, -cospi_8_64, -cospi_8_64, -cospi_8_64, cospi_24_64, cospi_24_64, cospi_24_64, cospi_24_64};

  int16x8_t cospi_0_0 = {cospi_16_64, cospi_24_64, cospi_16_64, -cospi_8_64, cospi_16_64, cospi_24_64, cospi_16_64, -cospi_8_64};
  int16x8_t cospi_0_1 = {cospi_16_64, cospi_8_64, -cospi_16_64, cospi_24_64, cospi_16_64, cospi_8_64, -cospi_16_64, cospi_24_64};
  int16x8_t cospi_0_2 = {cospi_16_64, cospi_24_64, cospi_16_64, -cospi_8_64, cospi_16_64, cospi_24_64, cospi_16_64, -cospi_8_64};
  int16x8_t cospi_0_3 = {cospi_16_64, cospi_8_64, -cospi_16_64, cospi_24_64, cospi_16_64, cospi_8_64, -cospi_16_64, cospi_24_64};

  loaded[0] = vec_vsx_ld(0, input);
  loaded[1] = vec_vsx_ld(0, input + (2 * stride));
  loaded[1] = vec_perm(loaded[1], loaded[1], perm0);

  in_high[0] = vec_sl(loaded[0], four);
  in_high[1] = vec_sl(loaded[1], four);

  if (in_high[0][0]) {
    ++in_high[0][0];
  }

  step[0] = vec_add(in_high[0], in_high[1]);
  step[1] = vec_sub(in_high[0], in_high[1]);
  step[1] = vec_perm(step[1], step[1], perm0);

  e[0] = vec_mule(step[0], cospi_0);
  e[1] = vec_mule(step[0], cospi_1);
  e[2] = vec_mule(step[1], cospi_2);
  e[3] = vec_mule(step[1], cospi_3);

  o[0] = vec_mulo(step[0], cospi_0);
  o[1] = vec_mulo(step[0], cospi_1);
  o[2] = vec_mulo(step[1], cospi_2);
  o[3] = vec_mulo(step[1], cospi_3);

  h[0] = vec_mergeh(e[0], o[0]);
  h[1] = vec_mergeh(e[1], o[1]);
  h[2] = vec_mergeh(e[2], o[2]);
  h[3] = vec_mergeh(e[3], o[3]);

  l[0] = vec_mergel(e[0], o[0]);
  l[1] = vec_mergel(e[1], o[1]);
  l[2] = vec_mergel(e[2], o[2]);
  l[3] = vec_mergel(e[3], o[3]);

  v[0] = vec_add(h[0], l[0]);
  v[2] = vec_add(h[1], l[1]);
  v[1] = vec_add(h[2], l[2]);
  v[3] = vec_add(h[3], l[3]);

  vpx_transpose_s32_4x4(v);

  v[0] = fdct_vector_round_shift(v[0]);
  v[1] = fdct_vector_round_shift(v[1]);
  v[2] = fdct_vector_round_shift(v[2]);
  v[3] = fdct_vector_round_shift(v[3]);

#ifdef WORDS_BIGENDIAN
  out[0] = vec_pack(v[1], v[0]);
  out[1] = vec_pack(v[3], v[2]);
#else
  out[0] = vec_pack(v[0], v[1]);
  out[1] = vec_pack(v[2], v[3]);
#endif // WORDS_BIGENDIAN

  // Do next column (which is a transposed row in second/horizontal pass)
  // Transform.

  out[1] = vec_perm(out[1], out[1], perm0);

  step[0] = vec_add(out[0], out[1]);
  step[1] = vec_sub(out[0], out[1]);

  x[0] = vec_perm(step[0], step[1], perm1);
  x[1] = vec_perm(step[0], step[1], perm2);
  x[2] = vec_perm(step[0], step[1], perm3);
  x[3] = vec_perm(step[0], step[1], perm4);

  e[0] = vec_mule(x[0], cospi_0_0);
  e[1] = vec_mule(x[1], cospi_0_1);
  e[2] = vec_mule(x[2], cospi_0_2);
  e[3] = vec_mule(x[3], cospi_0_3);

  o[0] = vec_mulo(x[0], cospi_0_0);
  o[1] = vec_mulo(x[1], cospi_0_1);
  o[2] = vec_mulo(x[2], cospi_0_2);
  o[3] = vec_mulo(x[3], cospi_0_3);

  h[0] = vec_mergeh(e[0], o[0]);
  h[1] = vec_mergeh(e[1], o[1]);
  h[2] = vec_mergeh(e[2], o[2]);
  h[3] = vec_mergeh(e[3], o[3]);

  l[0] = vec_mergel(e[0], o[0]);
  l[1] = vec_mergel(e[1], o[1]);
  l[2] = vec_mergel(e[2], o[2]);
  l[3] = vec_mergel(e[3], o[3]);

  v[0] = vec_add(h[0], h[1]);
  v[1] = vec_add(l[0], l[1]);
  v[2] = vec_add(h[2], h[3]);
  v[3] = vec_add(l[2], l[3]);

  v[0] = fdct_vector_round_shift(v[0]);
  v[1] = fdct_vector_round_shift(v[1]);
  v[2] = fdct_vector_round_shift(v[2]);
  v[3] = fdct_vector_round_shift(v[3]);

  tmp[0] = vec_add(v[0], one);
  tmp[1] = vec_add(v[1], one);
  tmp[2] = vec_add(v[2], one);
  tmp[3] = vec_add(v[3], one);

  tmp[0] = vec_sra(tmp[0], two);
  tmp[1] = vec_sra(tmp[1], two);
  tmp[2] = vec_sra(tmp[2], two);
  tmp[3] = vec_sra(tmp[3], two);

#ifdef WORDS_BIGENDIAN
  out[0] = vec_pack(tmp[1], tmp[0]);
  out[1] = vec_pack(tmp[3], tmp[2]);
#else
  out[0] = vec_pack(tmp[0], tmp[1]);
  out[1] = vec_pack(tmp[2], tmp[3]);
#endif // WORDS_BIGENDIAN

  vec_vsx_st(out[0], 0, output);
  vec_vsx_st(out[1], 0, output + (2* stride));
}
