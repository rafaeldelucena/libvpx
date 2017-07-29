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
#include "./vpx_dsp_rtcd.h"
#include "vpx_dsp/fwd_txfm.h"
//#include "vpx_dsp/types_vsx.h"

void vpx_fdct4x4_vsx(const int16_t *input, tran_low_t *output, int stride) {
  // The 2D transform is done with two passes which are actually pretty
  // similar. In the first one, we transform the columns and transpose
  // the results. In the second one, we transform the rows. To achieve that,
  // as the first pass results are transposed, we transpose the columns (that
  // is the transposed rows) and transpose the results (so that it goes back
  // in normal/row positions).
  // We need an intermediate buffer between passes.
  tran_low_t intermediate[16];
  // Do the two transform/transpose passes
  {
    tran_high_t in_high[16];    // canbe16
    tran_high_t step[16];       // canbe16
    tran_high_t temp[16];  // needs32

    //int16x8_t step1, step2, step3, step4;
    //int16x8_t in_high1, in_high2;
    {
      // Load inputs.
      {
        in_high[0] = input[0 * stride] * 16;
        in_high[1] = input[1 * stride] * 16;
        in_high[2] = input[2 * stride] * 16;
        in_high[3] = input[3 * stride] * 16;
        in_high[4] = (input[0 * stride + 1]) * 16;
        in_high[5] = (input[1 * stride + 1]) * 16;
        in_high[6] = (input[2 * stride + 1]) * 16;
        in_high[7] = (input[3 * stride + 1]) * 16;

        in_high[8] = (input[0 * stride + 2]) * 16;
        in_high[9] = (input[1 * stride + 2]) * 16;
        in_high[10] = (input[2 * stride + 2]) * 16;
        in_high[11] = (input[3 * stride + 2]) * 16;

        in_high[12] = (input[0 * stride + 3]) * 16;
        in_high[13] = (input[1 * stride + 3]) * 16;
        in_high[14] = (input[2 * stride + 3]) * 16;
        in_high[15] = (input[3 * stride + 3]) * 16;

        if (in_high[0]) {
          ++in_high[0];
        }
      }
      // Transform.
      //step1 = vec_add(in_high1, in_high2);
      step[0] = in_high[0] + in_high[3];
      step[1] = in_high[1] + in_high[2];
      step[4] = in_high[4] + in_high[7];
      step[5] = in_high[5] + in_high[6];
      step[8] = in_high[8] + in_high[11];
      step[9] = in_high[9] + in_high[10];
      step[12] = in_high[12] + in_high[15];
      step[13] = in_high[13] + in_high[14];

      //step2 = vec_sub(in_high1, in_high2);
      step[3] = in_high[0] - in_high[3];
      step[2] = in_high[1] - in_high[2];
      step[7] = in_high[4] - in_high[7];
      step[6] = in_high[5] - in_high[6];
      step[11] = in_high[8] - in_high[11];
      step[10] = in_high[9] - in_high[10];
      step[15] = in_high[12] - in_high[15];
      step[14] = in_high[13] - in_high[14];


      temp[0] = (step[0] + step[1]) * cospi_16_64;
      temp[1] = (step[0] - step[1]) * cospi_16_64;
      temp[2] = step[2] * cospi_24_64 + step[3] * cospi_8_64;
      temp[3] = -step[2] * cospi_8_64 + step[3] * cospi_24_64;

      temp[4] = (step[4] + step[5]) * cospi_16_64;
      temp[5] = (step[4] - step[5]) * cospi_16_64;
      temp[6] = step[6] * cospi_24_64 + step[7] * cospi_8_64;
      temp[7] = -step[6] * cospi_8_64 + step[7] * cospi_24_64;

      temp[8] = (step[8] + step[9]) * cospi_16_64;
      temp[9] = (step[8] - step[9]) * cospi_16_64;
      temp[10] = step[10] * cospi_24_64 + step[11] * cospi_8_64;
      temp[11] = -step[10] * cospi_8_64 + step[11] * cospi_24_64;

      temp[12] = (step[12] + step[13]) * cospi_16_64;
      temp[13] = (step[12] - step[13]) * cospi_16_64;
      temp[14] = step[14] * cospi_24_64 + step[15] * cospi_8_64;
      temp[15] = -step[14] * cospi_8_64 + step[15] * cospi_24_64;

      intermediate[0] = (tran_low_t)fdct_round_shift(temp[0]);
      intermediate[1] = (tran_low_t)fdct_round_shift(temp[2]);
      intermediate[2] = (tran_low_t)fdct_round_shift(temp[1]);
      intermediate[3] = (tran_low_t)fdct_round_shift(temp[3]);

      intermediate[4] = (tran_low_t)fdct_round_shift(temp[4]);
      intermediate[5] = (tran_low_t)fdct_round_shift(temp[6]);
      intermediate[6] = (tran_low_t)fdct_round_shift(temp[5]);
      intermediate[7] = (tran_low_t)fdct_round_shift(temp[7]);

      intermediate[8] = (tran_low_t)fdct_round_shift(temp[8]);
      intermediate[9] = (tran_low_t)fdct_round_shift(temp[10]);
      intermediate[10] = (tran_low_t)fdct_round_shift(temp[9]);
      intermediate[11] = (tran_low_t)fdct_round_shift(temp[11]);

      intermediate[12] = (tran_low_t)fdct_round_shift(temp[12]);
      intermediate[13] = (tran_low_t)fdct_round_shift(temp[14]);
      intermediate[14] = (tran_low_t)fdct_round_shift(temp[13]);
      intermediate[15] = (tran_low_t)fdct_round_shift(temp[15]);

      // Do next column (which is a transposed row in second/horizontal pass)
    }
  }
  {
    tran_high_t step[16];       // canbe16
    tran_high_t temp[16];  // needs32

    //int16x8_t step1, step2, step3, step4;
    //int16x8_t in_high1, in_high2;
    {
      // Transform.
      //step1 = vec_add(in_high1, in_high2);
      step[0] = intermediate[0] + intermediate[12];
      step[1] = intermediate[4] + intermediate[8];
      step[4] = intermediate[1] + intermediate[13];
      step[5] = intermediate[5] + intermediate[9];
      step[8] = intermediate[2] + intermediate[14];
      step[9] = intermediate[6] + intermediate[10];
      step[12] = intermediate[3] + intermediate[15];
      step[13] = intermediate[7] + intermediate[11];

      //step2 = vec_sub(in_high1, in_high2);
      step[3] = intermediate[0] - intermediate[12];
      step[2] = intermediate[4] - intermediate[8];
      step[7] = intermediate[1] - intermediate[13];
      step[6] = intermediate[5] - intermediate[9];
      step[11] = intermediate[2] - intermediate[14];
      step[10] = intermediate[6] - intermediate[10];
      step[15] = intermediate[3] - intermediate[15];
      step[14] = intermediate[7] - intermediate[11];


      tran_high_t x_0_0 = step[0] * cospi_16_64;
      tran_high_t x_0_1 = step[1] * cospi_16_64;
      temp[0] = x_0_0 + x_0_1;

      tran_high_t x_1_0 = step[2] * cospi_24_64;
      tran_high_t x_1_1 = step[3] * cospi_8_64;
      temp[1] = x_1_0 + x_1_1;

      tran_high_t x_2_0 = step[0] * cospi_16_64;
      tran_high_t x_2_1 = step[1] * cospi_16_64;
      temp[2] = x_2_0 - x_2_1;

      tran_high_t x_3_0 = -step[2] * cospi_8_64;
      tran_high_t x_3_1 = step[3] * cospi_24_64;
      temp[3] = x_3_0 + x_3_1;

      tran_high_t x_4_0 = step[4] * cospi_16_64;
      tran_high_t x_4_1 = step[5] * cospi_16_64;
      temp[4] = x_4_0 + x_4_1;

      tran_high_t x_5_0 = step[6] * cospi_24_64;
      tran_high_t x_5_1 = step[7] * cospi_8_64;
      temp[5] = step[6] * cospi_24_64 + step[7] * cospi_8_64;

      tran_high_t x_6_0 = step[4] * cospi_16_64;
      tran_high_t x_6_1 = step[5] * cospi_16_64;
      temp[6] = x_6_0 - x_6_1;

      tran_high_t x_7_0 = -step[6] * cospi_8_64;
      tran_high_t x_7_1 = step[7] * cospi_24_64;
      temp[7] = x_7_0 + x_7_1;

      tran_high_t x_8_0 = step[8] * cospi_16_64;
      tran_high_t x_8_1 = step[9] * cospi_16_64;
      temp[8] = x_8_0 + x_8_1;

      tran_high_t x_9_0 = step[10] * cospi_24_64;
      tran_high_t x_9_1 = step[11] * cospi_8_64;
      temp[9] = x_9_0 + x_9_1;

      tran_high_t x_10_0 = step[8] * cospi_16_64;
      tran_high_t x_10_1 = step[9] * cospi_16_64;
      temp[10] = x_10_0 - x_10_1;

      tran_high_t x_11_0 = -step[10] * cospi_8_64;
      tran_high_t x_11_1 = step[11] * cospi_24_64;
      temp[11] = x_11_0 + x_11_1;

      tran_high_t x_12_0 = step[12] * cospi_16_64;
      tran_high_t x_12_1 = step[13] * cospi_16_64;
      temp[12] = x_12_0 + x_12_1;

      tran_high_t x_13_0 = step[14] * cospi_24_64;
      tran_high_t x_13_1 = step[15] * cospi_8_64;
      temp[13] = x_13_0 + x_13_1;

      tran_high_t x_14_0 = step[12] * cospi_16_64;
      tran_high_t x_14_1 = step[13] * cospi_16_64;
      temp[14] = x_14_0 - x_14_1;

      tran_high_t x_15_0 = -step[14] * cospi_8_64;
      tran_high_t x_15_1 = step[15] * cospi_24_64;
      temp[15] = x_15_0 + x_15_1;

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

      // Do next column (which is a transposed row in second/horizontal pass)
    }
  }

  {
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
}
