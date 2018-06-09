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
  int i;
  for (i = 0; i < 4; ++i) {
    // Transform.
    tran_high_t step[4];       // canbe16
    tran_high_t temp1, temp2;  // needs32
    step[0] = input[i * 4 + 0] + input[i * 4 + 3];
    step[1] = input[i * 4 + 1] + input[i * 4 + 2];
    step[2] = input[i * 4 + 1] - input[i * 4 + 2];
    step[3] = input[i * 4 + 0] - input[i * 4 + 3];
    temp1 = (step[0] + step[1]) * cospi_16_64;
    temp2 = (step[0] - step[1]) * cospi_16_64;
    output[i * 4 + 0] = (tran_low_t)fdct_round_shift(temp1);
    output[i * 4 + 2] = (tran_low_t)fdct_round_shift(temp2);
    temp1 = step[2] * cospi_24_64 + step[3] * cospi_8_64;
    temp2 = -step[2] * cospi_8_64 + step[3] * cospi_24_64;
    output[i * 4 + 1] = (tran_low_t)fdct_round_shift(temp1);
    output[i * 4 + 3] = (tran_low_t)fdct_round_shift(temp2);
  }
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
