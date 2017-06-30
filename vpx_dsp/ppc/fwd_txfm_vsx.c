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
  int pass;
  // We need an intermediate buffer between passes.
  tran_low_t intermediate[16];
  const tran_low_t *in_low = NULL;
  tran_low_t *out = intermediate;
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
        if (in_high[0]) {
          ++in_high[0];
        }
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

      out[0] = (tran_low_t)fdct_round_shift(temp[0]);
      out[1] = (tran_low_t)fdct_round_shift(temp[2]);
      out[2] = (tran_low_t)fdct_round_shift(temp[1]);
      out[3] = (tran_low_t)fdct_round_shift(temp[3]);

      out[4] = (tran_low_t)fdct_round_shift(temp[4]);
      out[5] = (tran_low_t)fdct_round_shift(temp[6]);
      out[6] = (tran_low_t)fdct_round_shift(temp[5]);
      out[7] = (tran_low_t)fdct_round_shift(temp[7]);

      out[8] = (tran_low_t)fdct_round_shift(temp[8]);
      out[9] = (tran_low_t)fdct_round_shift(temp[10]);
      out[10] = (tran_low_t)fdct_round_shift(temp[9]);
      out[11] = (tran_low_t)fdct_round_shift(temp[11]);

      out[12] = (tran_low_t)fdct_round_shift(temp[12]);
      out[13] = (tran_low_t)fdct_round_shift(temp[14]);
      out[14] = (tran_low_t)fdct_round_shift(temp[13]);
      out[15] = (tran_low_t)fdct_round_shift(temp[15]);

      // Do next column (which is a transposed row in second/horizontal pass)
    }
    // Setup in/out for next pass.
    in_low = intermediate;
  }
  {
    tran_high_t in_high[16];    // canbe16
    tran_high_t step[16];       // canbe16
    tran_high_t temp[16];  // needs32

    //int16x8_t step1, step2, step3, step4;
    //int16x8_t in_high1, in_high2;
    {
      // Load inputs.
      {
        // in_high1
        in_high[0] = in_low[0];
        in_high[4] = in_low[1];
        in_high[8] = in_low[2];
        in_high[12]= in_low[3];
        in_high[1] = in_low[4];
        in_high[5] = in_low[5];
        in_high[9] = in_low[6];
        in_high[13]= in_low[7];

        // in_high2
        in_high[2] = in_low[8];
        in_high[6] = in_low[9];
        in_high[10]= in_low[10];
        in_high[14]= in_low[11];
        in_high[3] = in_low[12];
        in_high[7] = in_low[13];
        in_high[11]= in_low[14];
        in_high[15]= in_low[15];
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

      out[0] = (tran_low_t)fdct_round_shift(temp[0]);
      out[1] = (tran_low_t)fdct_round_shift(temp[2]);
      out[2] = (tran_low_t)fdct_round_shift(temp[1]);
      out[3] = (tran_low_t)fdct_round_shift(temp[3]);

      out[4] = (tran_low_t)fdct_round_shift(temp[4]);
      out[5] = (tran_low_t)fdct_round_shift(temp[6]);
      out[6] = (tran_low_t)fdct_round_shift(temp[5]);
      out[7] = (tran_low_t)fdct_round_shift(temp[7]);

      out[8] = (tran_low_t)fdct_round_shift(temp[8]);
      out[9] = (tran_low_t)fdct_round_shift(temp[10]);
      out[10] = (tran_low_t)fdct_round_shift(temp[9]);
      out[11] = (tran_low_t)fdct_round_shift(temp[11]);

      out[12] = (tran_low_t)fdct_round_shift(temp[12]);
      out[13] = (tran_low_t)fdct_round_shift(temp[14]);
      out[14] = (tran_low_t)fdct_round_shift(temp[13]);
      out[15] = (tran_low_t)fdct_round_shift(temp[15]);

      // Do next column (which is a transposed row in second/horizontal pass)
    }
    // Setup in/out for next pass.
    out = output;
  }

  {
    output[0] = (output[9] + 1) >> 2;
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
