/*
 *  Copyright (c) 2018 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "./vpx_dsp_rtcd.h"
#include "vpx_dsp/fwd_txfm.h"
#include "vpx_dsp/ppc/types_vsx.h"
#include "vpx_dsp/ppc/transpose_vsx.h"

static uint8x16_t perm_1_v = {0x00, 0x01, 0x02, 0x03, 0x08, 0x09, 0x0A, 0x0B, 0x10, 0x11, 0x12, 0x13, 0x18, 0x19, 0x1A, 0x1B};
static uint8x16_t perm_2_v = {0x06, 0x07, 0x04, 0x05, 0x0E, 0x0F, 0x0C, 0x0D, 0x16, 0x17, 0x14, 0x15, 0x1E, 0x1F, 0x1C, 0x1D};
static uint8x16_t perm_3_v = {0x02, 0x03, 0x10, 0x11, 0x02, 0x03, 0x10, 0x11, 0x06, 0x07, 0x14, 0x15, 0x06, 0x07, 0x14, 0x15};
static uint8x16_t perm_4_v = {0x0A, 0x0B, 0x18, 0x19, 0x0A, 0x0B, 0x18, 0x19, 0x0E, 0x0F, 0x1C, 0x1D, 0x0E, 0x0F, 0x1C, 0x1D};
static uint8x16_t perm_5_v = {0x00, 0x01, 0x12, 0x13, 0x00, 0x01, 0x12, 0x13, 0x04, 0x05, 0x16, 0x17, 0x04, 0x05, 0x16, 0x17};
static uint8x16_t perm_6_v = {0x08, 0x09, 0x1A, 0x1B, 0x08, 0x09, 0x1A, 0x1B, 0x0C, 0x0D, 0x1E, 0x1F, 0x0C, 0x0D, 0x1E, 0x1F};

static int16x8_t cospi_1_v = { 11585, 15137, -11585, 6270, 11585, 15137, -11585, 6270 };
static int16x8_t cospi_2_v = { 11585, 15137, -11585, 6270, 11585, 15137, -11585, 6270 };
static int16x8_t cospi_3_v = { 11585, 6270, 11585, -15137, 11585, 6270, 11585, -15137 };
static int16x8_t cospi_4_v = { 11585, 6270, 11585, -15137, 11585, 6270, 11585, -15137 };

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

static INLINE void vpx_fdct4x4_one_pass(const int16x8_t input[2], int32x4_t output[4]) {
  int16x8_t v[2];
  int16x8_t s[2];
  int16x8_t a[2];
  int16x8_t b[2];

  int32x4_t b_e[2];
  int32x4_t b_o[2];
  int32x4_t a_e[2];
  int32x4_t a_o[2];
  int32x4_t b_h[2];
  int32x4_t b_l[2];
  int32x4_t a_h[2];
  int32x4_t a_l[2];
  int32x4_t tmp[4];

  v[0] = vec_perm(input[0], input[1], perm_1_v);
  v[1] = vec_perm(input[0], input[1], perm_2_v);

  s[0] = vec_add(v[0], v[1]);
  s[1] = vec_sub(v[0], v[1]);

  b[0] = vec_perm(s[0], s[1], perm_3_v);
  b[1] = vec_perm(s[0], s[1], perm_4_v);

  a[0] = vec_perm(s[0], s[1], perm_5_v);
  a[1] = vec_perm(s[0], s[1], perm_6_v);

  b_e[0] = vec_mule(b[0], cospi_1_v);
  b_o[0] = vec_mulo(b[0], cospi_1_v);
  b_h[0] = vec_mergeh(b_e[0], b_o[0]);
  b_l[0] = vec_mergel(b_e[0], b_o[0]);

  b_e[1] = vec_mule(b[1], cospi_2_v);
  b_o[1] = vec_mulo(b[1], cospi_2_v);
  b_h[1] = vec_mergeh(b_e[1], b_o[1]);
  b_l[1] = vec_mergel(b_e[1], b_o[1]);

  a_e[0] = vec_mule(a[0], cospi_3_v);
  a_o[0] = vec_mulo(a[0], cospi_3_v);
  a_h[0] = vec_mergeh(a_e[0], a_o[0]);
  a_l[0] = vec_mergel(a_e[0], a_o[0]);

  a_e[1] = vec_mule(a[1], cospi_4_v);
  a_o[1] = vec_mulo(a[1], cospi_4_v);
  a_h[1] = vec_mergeh(a_e[1], a_o[1]);
  a_l[1] = vec_mergel(a_e[1], a_o[1]);

  tmp[0] = vec_add(a_h[0], b_h[0]);
  tmp[1] = vec_add(a_l[0], b_l[0]);
  tmp[2] = vec_add(a_h[1], b_h[1]);
  tmp[3] = vec_add(a_l[1], b_l[1]);

  output[0] = fdct_vector_round_shift(tmp[0]);
  output[1] = fdct_vector_round_shift(tmp[1]);
  output[2] = fdct_vector_round_shift(tmp[2]);
  output[3] = fdct_vector_round_shift(tmp[3]);
}

void vpx_fdct4x4_vsx(const int16_t *input, tran_low_t *output, int stride) {
  // The 2D transform is done with two passes which are actually pretty
  // similar. In the first one, we transform the columns and transpose
  // the results. In the second one, we transform the rows. To achieve that,
  // as the first pass results are transposed, we transpose the columns (that
  // is the transposed rows) and transpose the results (so that it goes back
  // in normal/row positions).
  // We need an intermediate buffer between passes.

  int16x8_t v[2];
  int16x8_t v_out[2];
  int32x4_t v_intermediate[4];
  int32x4_t last[4];
  int32x4_t x[4];

  uint16x8_t c = vec_splat_u16(4);
  int32x4_t one = vec_splat_s32(1);
  uint32x4_t two = vec_splat_u32(2);

  // Do the two transform/transpose passes

  x[0] = unpack_u16_to_s32_h(vec_vsx_ld(0, input));
  x[1] = unpack_u16_to_s32_h(vec_vsx_ld(0, input + stride));
  x[2] = unpack_u16_to_s32_h(vec_vsx_ld(0, input + stride * 2));
  x[3] = unpack_u16_to_s32_h(vec_vsx_ld(0, input + stride * 3));

  vpx_transpose_s32_4x4(x);

#ifdef WORDS_BIGENDIAN
  v[0] = vec_pack(x[1], x[0]);
  v[1] = vec_pack(x[3], x[2]);
#else
  v[0] = vec_pack(x[0], x[1]);
  v[1] = vec_pack(x[2], x[3]);
#endif // WORDS_BIGENDIAN

  v[0] = vec_sl(v[0], c);
  v[1] = vec_sl(v[1], c);

  if (v[0][0]) {
    ++v[0][0];
  }

  // Transform.
  vpx_fdct4x4_one_pass(v, v_intermediate);


  // Do next column (which is a transposed row in second/horizontal pass)
  // Transpose

  vpx_transpose_s32_4x4(v_intermediate);

#ifdef WORDS_BIGENDIAN
  v_out[0] = vec_pack(v_intermediate[1], v_intermediate[0]);
  v_out[1] = vec_pack(v_intermediate[3], v_intermediate[2]);
#else
  v_out[0] = vec_pack(v_intermediate[0], v_intermediate[1]);
  v_out[1] = vec_pack(v_intermediate[2], v_intermediate[3]);
#endif // WORDS_BIGENDIAN

  // Transform.
  vpx_fdct4x4_one_pass(v_out, last);

  last[0] = vec_add(last[0], one);
  last[1] = vec_add(last[1], one);
  last[2] = vec_add(last[2], one);
  last[3] = vec_add(last[3], one);

  last[0] = vec_sra(last[0], two);
  last[1] = vec_sra(last[1], two);
  last[2] = vec_sra(last[2], two);
  last[3] = vec_sra(last[3], two);

#ifdef WORDS_BIGENDIAN
  v_out[0] = vec_pack(last[1], last[0]);
  v_out[1] = vec_pack(last[3], last[2]);
#else
  v_out[0] = vec_pack(last[0], last[1]);
  v_out[1] = vec_pack(last[2], last[3]);
#endif // WORDS_BIGENDIAN

  vec_vsx_st(v_out[0], 0, output);
  vec_vsx_st(v_out[1], 0, output + 8);
}
