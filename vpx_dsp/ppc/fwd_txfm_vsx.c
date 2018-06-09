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
  tran_high_t step_0[4];       // canbe16
  tran_high_t temp_0[4];       // needs32
  tran_high_t step_1[4];       // canbe16
  tran_high_t temp_1[4];       // needs32
  tran_high_t step_2[4];       // canbe16
  tran_high_t temp_2[4];       // needs32
  tran_high_t step_3[4];       // canbe16
  tran_high_t temp_3[4];       // needs32

  // Transform.
  step_0[0] = input[0 * 4 + 0] + input[0 * 4 + 3];
  step_0[1] = input[0 * 4 + 1] + input[0 * 4 + 2];
  step_0[2] = input[0 * 4 + 1] - input[0 * 4 + 2];
  step_0[3] = input[0 * 4 + 0] - input[0 * 4 + 3];
  step_1[0] = input[1 * 4 + 0] + input[1 * 4 + 3];
  step_1[1] = input[1 * 4 + 1] + input[1 * 4 + 2];
  step_1[2] = input[1 * 4 + 1] - input[1 * 4 + 2];
  step_1[3] = input[1 * 4 + 0] - input[1 * 4 + 3];
  step_2[0] = input[2 * 4 + 0] + input[2 * 4 + 3];
  step_2[1] = input[2 * 4 + 1] + input[2 * 4 + 2];
  step_2[2] = input[2 * 4 + 1] - input[2 * 4 + 2];
  step_2[3] = input[2 * 4 + 0] - input[2 * 4 + 3];
  step_3[0] = input[3 * 4 + 0] + input[3 * 4 + 3];
  step_3[1] = input[3 * 4 + 1] + input[3 * 4 + 2];
  step_3[2] = input[3 * 4 + 1] - input[3 * 4 + 2];
  step_3[3] = input[3 * 4 + 0] - input[3 * 4 + 3];

  temp_0[0] = (step_0[0] + step_0[1]) * cospi_16_64;
  temp_0[1] = (step_0[0] - step_0[1]) * cospi_16_64;
  output[0 * 4 + 0] = (tran_low_t)fdct_round_shift(temp_0[0]);
  output[0 * 4 + 2] = (tran_low_t)fdct_round_shift(temp_0[1]);
  temp_0[2] = step_0[2] * cospi_24_64 + step_0[3] * cospi_8_64;
  temp_0[3] = -step_0[2] * cospi_8_64 + step_0[3] * cospi_24_64;
  output[0 * 4 + 1] = (tran_low_t)fdct_round_shift(temp_0[2]);
  output[0 * 4 + 3] = (tran_low_t)fdct_round_shift(temp_0[3]);
  // Transform.
  temp_1[0] = (step_1[0] + step_1[1]) * cospi_16_64;
  temp_1[1] = (step_1[0] - step_1[1]) * cospi_16_64;
  output[1 * 4 + 0] = (tran_low_t)fdct_round_shift(temp_1[0]);
  output[1 * 4 + 2] = (tran_low_t)fdct_round_shift(temp_1[1]);
  temp_1[2] = step_1[2] * cospi_24_64 + step_1[3] * cospi_8_64;
  temp_1[3] = -step_1[2] * cospi_8_64 + step_1[3] * cospi_24_64;
  output[1 * 4 + 1] = (tran_low_t)fdct_round_shift(temp_1[2]);
  output[1 * 4 + 3] = (tran_low_t)fdct_round_shift(temp_1[3]);
  // Transform.
  temp_2[0] = (step_2[0] + step_2[1]) * cospi_16_64;
  temp_2[1] = (step_2[0] - step_2[1]) * cospi_16_64;
  output[2 * 4 + 0] = (tran_low_t)fdct_round_shift(temp_2[0]);
  output[2 * 4 + 2] = (tran_low_t)fdct_round_shift(temp_2[1]);
  temp_2[2] = step_2[2] * cospi_24_64 + step_2[3] * cospi_8_64;
  temp_2[3] = -step_2[2] * cospi_8_64 + step_2[3] * cospi_24_64;
  output[2 * 4 + 1] = (tran_low_t)fdct_round_shift(temp_2[2]);
  output[2 * 4 + 3] = (tran_low_t)fdct_round_shift(temp_2[3]);
  // Transform.
  temp_3[0] = (step_3[0] + step_3[1]) * cospi_16_64;
  temp_3[1] = (step_3[0] - step_3[1]) * cospi_16_64;
  output[3 * 4 + 0] = (tran_low_t)fdct_round_shift(temp_3[0]);
  output[3 * 4 + 2] = (tran_low_t)fdct_round_shift(temp_3[1]);
  temp_3[2] = step_3[2] * cospi_24_64 + step_3[3] * cospi_8_64;
  temp_3[3] = -step_3[2] * cospi_8_64 + step_3[3] * cospi_24_64;
  output[3 * 4 + 1] = (tran_low_t)fdct_round_shift(temp_3[2]);
  output[3 * 4 + 3] = (tran_low_t)fdct_round_shift(temp_3[3]);
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
