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
  tran_high_t step[16]; // canbe16
  tran_high_t temp[16]; // needs32
  tran_high_t a[16]; // needs32
  tran_high_t b[16]; // needs32

  step[0] = input[0] + input[3];
  step[1] = input[1] + input[2];
  step[2] = input[1] - input[2];
  step[3] = input[0] - input[3];
  step[4] = input[4] + input[7];
  step[5] = input[5] + input[6];
  step[6] = input[5] - input[6];
  step[7] = input[4] - input[7];
  step[8] = input[8] + input[11];
  step[9] = input[9] + input[10];
  step[10] = input[9] - input[10];
  step[11] = input[8] - input[11];
  step[12] = input[12] + input[15];
  step[13] = input[13] + input[14];
  step[14] = input[13] - input[14];
  step[15] = input[12] - input[15];

  a[0] = step[0] * cospi_16_64;
  a[1] = step[0] * cospi_16_64;
  a[2] = step[2] * cospi_24_64;
  a[3] = step[2] * -cospi_8_64;
  a[4] = step[4] * cospi_16_64;
  a[5] = step[4] * cospi_16_64;
  a[6] = step[6] * cospi_24_64;
  a[7] = step[6] * -cospi_8_64;
  a[8] = step[8] * cospi_16_64;
  a[9] = step[8] * cospi_16_64;
  a[10] = step[10] * cospi_24_64;
  a[11] = step[10] * -cospi_8_64;
  a[12] = step[12] * cospi_16_64;
  a[13] = step[12] * cospi_16_64;
  a[14] = step[14] * cospi_24_64;
  a[15] = step[14] * -cospi_8_64;

  b[0] = step[1] * cospi_16_64;
  b[1] = step[1] * -cospi_16_64;
  b[2] = step[3] * cospi_8_64;
  b[3] = step[3] * cospi_24_64;
  b[4] = step[5] * cospi_16_64;
  b[5] = step[5] * -cospi_16_64;
  b[6] = step[7] * cospi_8_64;
  b[7] = step[7] * cospi_24_64;
  b[8] = step[9] * cospi_16_64;
  b[9] = step[9] * -cospi_16_64;
  b[10] = step[11] * cospi_8_64;
  b[11] = step[11] * cospi_24_64;
  b[12] = step[13] * cospi_16_64;
  b[13] = step[13] * -cospi_16_64;
  b[14] = step[15] * cospi_8_64;
  b[15] = step[15] * cospi_24_64;

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
  output[2] = (tran_low_t)fdct_round_shift(temp[1]);
  output[1] = (tran_low_t)fdct_round_shift(temp[2]);
  output[3] = (tran_low_t)fdct_round_shift(temp[3]);
  output[4] = (tran_low_t)fdct_round_shift(temp[4]);
  output[6] = (tran_low_t)fdct_round_shift(temp[5]);
  output[5] = (tran_low_t)fdct_round_shift(temp[6]);
  output[7] = (tran_low_t)fdct_round_shift(temp[7]);
  output[8] = (tran_low_t)fdct_round_shift(temp[8]);
  output[10] = (tran_low_t)fdct_round_shift(temp[9]);
  output[9] = (tran_low_t)fdct_round_shift(temp[10]);
  output[11] = (tran_low_t)fdct_round_shift(temp[11]);
  output[12] = (tran_low_t)fdct_round_shift(temp[12]);
  output[14] = (tran_low_t)fdct_round_shift(temp[13]);
  output[13] = (tran_low_t)fdct_round_shift(temp[14]);
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
  // Do the two transform/transpose passes
  {
    tran_high_t in[4 * 4];    // canbe16
    tran_high_t *in_high = in;
    // Load inputs.
    in_high[0] = input[0 * stride] * 16;
    if (in_high[0]) {
      ++in_high[0];
    }
    in_high[1] = input[1 * stride] * 16;
    in_high[2] = input[2 * stride] * 16;
    in_high[3] = input[3 * stride] * 16;

    in_high += 4;

    in_high[0] = input[0 * stride + 1] * 16;
    in_high[1] = input[1 * stride + 1] * 16;
    in_high[2] = input[2 * stride + 1] * 16;
    in_high[3] = input[3 * stride + 1] * 16;

    in_high += 4;

    in_high[0] = input[0 * stride + 2] * 16;
    in_high[1] = input[1 * stride + 2] * 16;
    in_high[2] = input[2 * stride + 2] * 16;
    in_high[3] = input[3 * stride + 2] * 16;

    in_high += 4;

    in_high[0] = input[0 * stride + 3] * 16;
    in_high[1] = input[1 * stride + 3] * 16;
    in_high[2] = input[2 * stride + 3] * 16;
    in_high[3] = input[3 * stride + 3] * 16;

    // Transform.
    vpx_fdct4x4_one_pass(in, intermediate);
  }
  // Do next column (which is a transposed row in second/horizontal pass)
  {
    tran_high_t in[4 * 4];    // canbe16
    tran_high_t *in_high = in;
    // Load inputs.
    in_high[0] = intermediate[0 * 4];
    in_high[1] = intermediate[1 * 4];
    in_high[2] = intermediate[2 * 4];
    in_high[3] = intermediate[3 * 4];

    in_high +=4;

    in_high[0] = intermediate[0 * 4 + 1];
    in_high[1] = intermediate[1 * 4 + 1];
    in_high[2] = intermediate[2 * 4 + 1];
    in_high[3] = intermediate[3 * 4 + 1];

    in_high +=4;

    in_high[0] = intermediate[0 * 4 + 2];
    in_high[1] = intermediate[1 * 4 + 2];
    in_high[2] = intermediate[2 * 4 + 2];
    in_high[3] = intermediate[3 * 4 + 2];

    in_high +=4;

    in_high[0] = intermediate[0 * 4 + 3];
    in_high[1] = intermediate[1 * 4 + 3];
    in_high[2] = intermediate[2 * 4 + 3];
    in_high[3] = intermediate[3 * 4 + 3];

    // Transform.
    vpx_fdct4x4_one_pass(in, output);
  }

  {
    output[0 + 0 * 4] = (output[0 + 0 * 4] + 1) >> 2;
    output[1 + 0 * 4] = (output[1 + 0 * 4] + 1) >> 2;
    output[2 + 0 * 4] = (output[2 + 0 * 4] + 1) >> 2;
    output[3 + 0 * 4] = (output[3 + 0 * 4] + 1) >> 2;

    output[0 + 1 * 4] = (output[0 + 1 * 4] + 1) >> 2;
    output[1 + 1 * 4] = (output[1 + 1 * 4] + 1) >> 2;
    output[2 + 1 * 4] = (output[2 + 1 * 4] + 1) >> 2;
    output[3 + 1 * 4] = (output[3 + 1 * 4] + 1) >> 2;

    output[0 + 2 * 4] = (output[0 + 2 * 4] + 1) >> 2;
    output[1 + 2 * 4] = (output[1 + 2 * 4] + 1) >> 2;
    output[2 + 2 * 4] = (output[2 + 2 * 4] + 1) >> 2;
    output[3 + 2 * 4] = (output[3 + 2 * 4] + 1) >> 2;

    output[0 + 3 * 4] = (output[0 + 3 * 4] + 1) >> 2;
    output[1 + 3 * 4] = (output[1 + 3 * 4] + 1) >> 2;
    output[2 + 3 * 4] = (output[2 + 3 * 4] + 1) >> 2;
    output[3 + 3 * 4] = (output[3 + 3 * 4] + 1) >> 2;
  }
}
