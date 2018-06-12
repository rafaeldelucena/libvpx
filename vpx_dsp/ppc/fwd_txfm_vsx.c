/*
 *  Copyright (c) 2018 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <assert.h>
#include "./vpx_dsp_rtcd.h"
#include "vpx_dsp/fwd_txfm.h"

void vpx_fdct4x4_one_pass(const tran_high_t *input, tran_low_t *output) {
  tran_high_t step1[8]; // canbe16
  tran_high_t step2[8]; // canbe16
  tran_high_t temp[16]; // needs32
  tran_high_t a[16]; // needs32
  tran_high_t b[16]; // needs32

  step1[0] = input[0] + input[3];
  step1[1] = input[1] + input[2];
  step1[2] = input[4] + input[7];
  step1[3] = input[5] + input[6];
  step1[4] = input[8] + input[11];
  step1[5] = input[9] + input[10];
  step1[6] = input[12] + input[15];
  step1[7] = input[13] + input[14];

  step2[0] = input[0] - input[3];
  step2[1] = input[1] - input[2];
  step2[2] = input[4] - input[7];
  step2[3] = input[5] - input[6];
  step2[4] = input[9] - input[10];
  step2[5] = input[8] - input[11];
  step2[6] = input[13] - input[14];
  step2[7] = input[12] - input[15];

  b[0] = step1[1] * cospi_16_64;
  b[1] = step2[0] * cospi_8_64;
  b[2] = step1[1] * -cospi_16_64;
  b[3] = step2[0] * cospi_24_64;
  b[4] = step1[3] * cospi_16_64;
  b[5] = step2[2] * cospi_8_64;
  b[6] = step1[3] * -cospi_16_64;
  b[7] = step2[2] * cospi_24_64;
  b[8] = step1[5] * cospi_16_64;
  b[9] = step2[5] * cospi_8_64;
  b[10] = step1[5] * -cospi_16_64;
  b[11] = step2[5] * cospi_24_64;
  b[12] = step1[7] * cospi_16_64;
  b[13] = step2[7] * cospi_8_64;
  b[14] = step1[7] * -cospi_16_64;
  b[15] = step2[7] * cospi_24_64;

  a[0] = step1[0] * cospi_16_64;
  a[1] = step2[1] * cospi_24_64;
  a[2] = step1[0] * cospi_16_64;
  a[3] = step2[1] * -cospi_8_64;
  a[4] = step1[2] * cospi_16_64;
  a[5] = step2[3] * cospi_24_64;
  a[6] = step1[2] * cospi_16_64;
  a[7] = step2[3] * -cospi_8_64;
  a[8] = step1[4] * cospi_16_64;
  a[9] = step2[4] * cospi_24_64;
  a[10] = step1[4] * cospi_16_64;
  a[11] = step2[4] * -cospi_8_64;
  a[12] = step1[6] * cospi_16_64;
  a[13] = step2[6] * cospi_24_64;
  a[14] = step1[6] * cospi_16_64;
  a[15] = step2[6] * -cospi_8_64;

  temp[0] = a[0] + b[0];
  temp[1] = a[1] + b[1];
  temp[2] = a[2] + b[2];
  temp[3] = a[3] + b[3];
  temp[4] = a[4] + b[4];
  temp[5] = a[5] + b[5];
  temp[6] = a[6] + b[6];
  temp[7] = a[7] + b[7];
  temp[8] = a[8] + b[8];
  temp[9] = a[9] + b[9];
  temp[10] = a[10] + b[10];
  temp[11] = a[11] + b[11];
  temp[12] = a[12] + b[12];
  temp[13] = a[13] + b[13];
  temp[14] = a[14] + b[14];
  temp[15] = a[15] + b[15];

  output[0] = (tran_low_t)fdct_round_shift(temp[0]);
  output[1] = (tran_low_t)fdct_round_shift(temp[1]);
  output[2] = (tran_low_t)fdct_round_shift(temp[2]);
  output[3] = (tran_low_t)fdct_round_shift(temp[3]);
  output[4] = (tran_low_t)fdct_round_shift(temp[4]);
  output[5] = (tran_low_t)fdct_round_shift(temp[5]);
  output[6] = (tran_low_t)fdct_round_shift(temp[6]);
  output[7] = (tran_low_t)fdct_round_shift(temp[7]);
  output[8] = (tran_low_t)fdct_round_shift(temp[8]);
  output[9] = (tran_low_t)fdct_round_shift(temp[9]);
  output[10] = (tran_low_t)fdct_round_shift(temp[10]);
  output[11] = (tran_low_t)fdct_round_shift(temp[11]);
  output[12] = (tran_low_t)fdct_round_shift(temp[12]);
  output[13] = (tran_low_t)fdct_round_shift(temp[13]);
  output[14] = (tran_low_t)fdct_round_shift(temp[14]);
  output[15] = (tran_low_t)fdct_round_shift(temp[15]);
}

void vpx_fdct4x4_vsx(const int16_t *input, tran_low_t *output, int stride) {
  // The 2D transform is done with two passes which are actually pretty
  // similar. In the first one, we transform the columns and transpose
  // the results. In the second one, we transform the rows. To achieve that,
  // as the first pass results are transposed, we transpose the columns (that
  // is the transposed rows) and transpose the results (so that it goes back
  // in normal/row positions).
  // We need an intermediate buffer between passes.
  tran_low_t intermediate[4 * 4];

  tran_high_t in[4 * 4]; // canbe16
  tran_high_t transposed[4 * 4]; // canbe16

  // Do the two transform/transpose passes

  // Load inputs.
  in[0] = input[0 * stride];
  in[1] = input[1 * stride];
  in[2] = input[2 * stride];
  in[3] = input[3 * stride];
  in[4] = input[0 * stride + 1];
  in[5] = input[1 * stride + 1];
  in[6] = input[2 * stride + 1];
  in[7] = input[3 * stride + 1];
  in[8] = input[0 * stride + 2];
  in[9] = input[1 * stride + 2];
  in[10] = input[2 * stride + 2];
  in[11] = input[3 * stride + 2];
  in[12] = input[0 * stride + 3];
  in[13] = input[1 * stride + 3];
  in[14] = input[2 * stride + 3];
  in[15] = input[3 * stride + 3];

  in[0] = in[0] * 16;
  in[1] = in[1] * 16;
  in[2] = in[2] * 16;
  in[3] = in[3] * 16;
  in[4] = in[4] * 16;
  in[5] = in[5] * 16;
  in[6] = in[6] * 16;
  in[7] = in[7] * 16;
  in[8] = in[8] * 16;
  in[9] = in[9] * 16;
  in[10] = in[10] * 16;
  in[11] = in[11] * 16;
  in[12] = in[12] * 16;
  in[13] = in[13] * 16;
  in[14] = in[14] * 16;
  in[15] = in[15] * 16;

  if (in[0]) {
    ++in[0];
  }

  // Transform.
  vpx_fdct4x4_one_pass(in, intermediate);

  // Do next column (which is a transposed row in second/horizontal pass)
  // Transpose
  transposed[0] = intermediate[0];
  transposed[1] = intermediate[4];
  transposed[2] = intermediate[8];
  transposed[3] = intermediate[12];
  transposed[4] = intermediate[1];
  transposed[5] = intermediate[5];
  transposed[6] = intermediate[9];
  transposed[7] = intermediate[13];
  transposed[8] = intermediate[2];
  transposed[9] = intermediate[6];
  transposed[10] = intermediate[10];
  transposed[11] = intermediate[14];
  transposed[12] = intermediate[3];
  transposed[13] = intermediate[7];
  transposed[14] = intermediate[11];
  transposed[15] = intermediate[15];

  // Transform.
  vpx_fdct4x4_one_pass(transposed, output);

  output[0] = (output[0] + 1) >> 2;
  output[1] = (output[1] + 1) >> 2;
  output[2] = (output[2] + 1) >> 2;
  output[3] = (output[3] + 1) >> 2;
  output[4] = (output[4] + 1) >> 2;
  output[5] = (output[5] + 1) >> 2;
  output[6] = (output[6] + 1) >> 2;
  output[7] = (output[7] + 1) >> 2;
  output[8] = (output[8] + 1) >> 2;
  output[9] = (output[9] + 1) >> 2;
  output[10] = (output[10] + 1) >> 2;
  output[11] = (output[11] + 1) >> 2;
  output[12] = (output[12] + 1) >> 2;
  output[13] = (output[13] + 1) >> 2;
  output[14] = (output[14] + 1) >> 2;
  output[15] = (output[15] + 1) >> 2;
}
