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

///* Shift down with rounding */
//#define ROUND_POWER_OF_TWO(value, n) (((value) + (1 << ((n)-1))) >> (n))
//
//static INLINE tran_high_t fdct_round_shift(tran_high_t input) {
//  tran_high_t rv = ROUND_POWER_OF_TWO(input, DCT_CONST_BITS);
//  // TODO(debargha, peter.derivaz): Find new bounds for this assert
//  // and make the bounds consts.
//  // assert(INT16_MIN <= rv && rv <= INT16_MAX);
//  return rv;
//}

static INLINE int16x8_t fdct_vector_round_shift(int16x8_t input) {
  int16x8_t one = vec_splat_s16(1);
  uint16x8_t dct_const_0 = vec_splat_u16(DCT_CONST_BITS - 1);
  uint16x8_t dct_const_1 = vec_splat_u16(DCT_CONST_BITS);
  int16x8_t tmp = vec_sl(one, dct_const_0);
  int16x8_t sum = vec_add(input, tmp);
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
  printf("---------------------------------------- INPUTS\n");
  MATRIX_H4_PRINT(input);
  tran_low_t intermediate[16];
  // Do the two transform/transpose passes
  {
    tran_high_t in_high[16];    // canbe16
    tran_high_t temp[16];  // needs32
    tran_high_t step[16];  // needs32

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
      printf("---------------------------------------- TRANSFORM INTERMEDIATE\n");
      MATRIX_H4_PRINT(intermediate);
    }
  }
  // Transform.

  //step1 = vec_add(in_high1, in_high2);
  //tran_high_t step_0 = intermediate[0] + intermediate[12];
  //tran_high_t step_1 = intermediate[1] + intermediate[13];
  //tran_high_t step_2 = intermediate[2] + intermediate[14];
  //tran_high_t step_3 = intermediate[3] + intermediate[15];
  //tran_high_t step_4 = intermediate[4] + intermediate[8];
  //tran_high_t step_5 = intermediate[5] + intermediate[9];
  //tran_high_t step_6 = intermediate[6] + intermediate[10];
  //tran_high_t step_7 = intermediate[7] + intermediate[11];

  //step2 = vec_sub(in_high1, in_high2);
  //tran_high_t step_8 = intermediate[0] - intermediate[12];
  //tran_high_t step_9 = intermediate[1] - intermediate[13];
  //tran_high_t step_10 = intermediate[2] - intermediate[14];
  //tran_high_t step_11 = intermediate[3] - intermediate[15];
  //tran_high_t step_12 = intermediate[4] - intermediate[8];
  //tran_high_t step_13 = intermediate[5] - intermediate[9];
  //tran_high_t step_14 = intermediate[6] - intermediate[10];
  //tran_high_t step_15 = intermediate[7] - intermediate[11];

  int16x8_t in1 = vec_vsx_ld(0, intermediate);
  int16x8_t in0 = vec_vsx_ld(0, intermediate + (2 * stride));
  printf("---------------------------------------- AFTER LOAD\n");
  MATRIX_SI4_PRINT(in1, in0);
  uint8x16_t perm0 = {0x8, 0x9, 0xA, 0xB, 0xC, 0xD, 0xE, 0xF, 0x0, 0x1, 0x2, 0x3, 0x4, 0x5, 0x6, 0x7};

  int16x8_t in2 = vec_perm(in0, in0, perm0);

  int16x8_t step1 = vec_add(in1, in2);
  int16x8_t step2 = vec_sub(in1, in2);

  //tran_high_t x_0_0 = step_0 * cospi_16_64;
  //tran_high_t x_1_0 = step_12 * cospi_24_64;
  //tran_high_t x_2_0 = step_0 * cospi_16_64;
  //tran_high_t x_3_0 = -step_12 * cospi_8_64;
  //tran_high_t x_4_0 = step_1 * cospi_16_64;
  //tran_high_t x_5_0 = step_13 * cospi_24_64;
  //tran_high_t x_6_0 = step_1 * cospi_16_64;
  //tran_high_t x_7_0 = -step_13 * cospi_8_64;

  uint8x16_t perm1 = {0x0, 0x1, 0x18, 0x19, 0x0, 0x1, 0x18, 0x19, 0x2, 0x3, 0x1A, 0x1B, 0x2, 0x3, 0x1A, 0x1B};
  int16x8_t x_0_0 = vec_perm(step1, step2, perm1);
  int16x8_t cospi_0_0 = {cospi_16_64, cospi_24_64, cospi_16_64, -cospi_8_64, cospi_16_64, cospi_24_64, cospi_16_64, -cospi_8_64};
  int32x4_t e_0_0 = vec_mule(x_0_0, cospi_0_0);
  int32x4_t o_0_0 = vec_mulo(x_0_0, cospi_0_0);
  int32x4_t h_0_0 = vec_mergeh(e_0_0, o_0_0);
  int32x4_t l_0_0 = vec_mergel(e_0_0, o_0_0);

  //tran_high_t x_0_1 = step_4 * cospi_16_64;
  //tran_high_t x_1_1 = step_8 * cospi_8_64;
  //tran_high_t x_2_1 = step_4 * cospi_16_64;
  //tran_high_t x_3_1 = step_8 * cospi_24_64;
  //tran_high_t x_4_1 = step_5 * cospi_16_64;
  //tran_high_t x_5_1 = step_9 * cospi_8_64;
  //tran_high_t x_6_1 = step_5 * cospi_16_64;
  //tran_high_t x_7_1 = step_9 * cospi_24_64;

  uint8x16_t perm2 = {0x8, 0x9, 0x10, 0x11, 0x8, 0x9, 0x10, 0x11, 0xA, 0xB, 0x12, 0x13, 0xA, 0xB, 0x12, 0x13};
  int16x8_t x_0_1 = vec_perm(step1, step2, perm2);

  int16x8_t cospi_0_1 = {cospi_16_64, cospi_8_64, -cospi_16_64, cospi_24_64, cospi_16_64, cospi_8_64, -cospi_16_64, cospi_24_64};
  int32x4_t e_0_1 = vec_mule(x_0_1, cospi_0_1);
  int32x4_t o_0_1 = vec_mulo(x_0_1, cospi_0_1);
  int32x4_t h_0_1 = vec_mergeh(e_0_1, o_0_1);
  int32x4_t l_0_1 = vec_mergel(e_0_1, o_0_1);

  //tran_high_t x_8_0 = step_2 * cospi_16_64;
  //tran_high_t x_9_0 = step_14 * cospi_24_64;
  //tran_high_t x_10_0 = step_2 * cospi_16_64;
  //tran_high_t x_11_0 = -step_14 * cospi_8_64;
  //tran_high_t x_12_0 = step_3 * cospi_16_64;
  //tran_high_t x_13_0 = step_15 * cospi_24_64;
  //tran_high_t x_14_0 = step_3 * cospi_16_64;
  //tran_high_t x_15_0 = -step_15 * cospi_8_64;

  uint8x16_t perm3 = {0x4, 0x5, 0x1C, 0x1D, 0x4, 0x5, 0x1C, 0x1D, 0x6, 0x7, 0x1E, 0x1F, 0x6, 0x7, 0x1E, 0x1F};
  int16x8_t x_0_2 = vec_perm(step1, step2, perm3);
  int16x8_t cospi_0_2 = {cospi_16_64, cospi_24_64, cospi_16_64, -cospi_8_64, cospi_16_64, cospi_24_64, cospi_16_64, -cospi_8_64};
  int32x4_t e_0_2 = vec_mule(x_0_2, cospi_0_2);
  int32x4_t o_0_2 = vec_mulo(x_0_2, cospi_0_2);
  int32x4_t h_0_2 = vec_mergeh(e_0_2, o_0_2);
  int32x4_t l_0_2 = vec_mergel(e_0_2, o_0_2);

  //tran_high_t x_8_1 = step_6 * cospi_16_64;
  //tran_high_t x_9_1 = step_10 * cospi_8_64;
  //tran_high_t x_10_1 = step_6 * cospi_16_64;
  //tran_high_t x_11_1 = step_10 * cospi_24_64;
  //tran_high_t x_12_1 = step_7 * cospi_16_64;
  //tran_high_t x_13_1 = step_11 * cospi_8_64;
  //tran_high_t x_14_1 = step_7 * cospi_16_64;
  //tran_high_t x_15_1 = step_11 * cospi_24_64;

  uint8x16_t perm4 = {0x14, 0x15, 0xC, 0xD, 0xC, 0xD, 0x14, 0x15, 0xE, 0xF, 0x16, 0x17, 0xE, 0xF, 0x16, 0x17};
  int16x8_t x_0_3 = vec_perm(step1, step2, perm4);
  int16x8_t cospi_0_3 = {cospi_16_64, cospi_8_64, -cospi_16_64, cospi_24_64, cospi_16_64, cospi_8_64, -cospi_16_64, cospi_24_64};
  int32x4_t e_0_3 = vec_mule(x_0_3, cospi_0_3);
  int32x4_t o_0_3 = vec_mulo(x_0_3, cospi_0_3);
  int32x4_t h_0_3 = vec_mergeh(e_0_3, o_0_3);
  int32x4_t l_0_3 = vec_mergel(e_0_3, o_0_3);

#ifdef WORDS_BIGENDIAN
  int16x8_t a_0 = vec_pack(l_0_0, h_0_0);
  int16x8_t a_1 = vec_pack(l_0_1, h_0_1);
  int16x8_t a_2 = vec_pack(l_0_2, h_0_2);
  int16x8_t a_3 = vec_pack(l_0_3, h_0_3);
#else
  int16x8_t a_0 = vec_pack(h_0_0, l_0_0);
  int16x8_t a_1 = vec_pack(h_0_1, l_0_1);
  int16x8_t a_2 = vec_pack(h_0_2, l_0_2);
  int16x8_t a_3 = vec_pack(h_0_3, l_0_3);
#endif // WORDS_BIGENDIAN

  //temp[0] = x_0_0 + x_0_1;
  //temp[1] = x_1_0 + x_1_1;
  //temp[2] = x_2_0 - x_2_1;
  //temp[3] = x_3_0 + x_3_1;
  //temp[4] = x_4_0 + x_4_1;
  //temp[5] = x_5_0 + x_5_1;
  //temp[6] = x_6_0 - x_6_1;
  //temp[7] = x_7_0 + x_7_1;

  //temp[8] = x_8_0 + x_8_1;
  //temp[9] = x_9_0 + x_9_1;
  //temp[10] = x_10_0 - x_10_1;
  //temp[11] = x_11_0 + x_11_1;
  //temp[12] = x_12_0 + x_12_1;
  //temp[13] = x_13_0 + x_13_1;
  //temp[14] = x_14_0 - x_14_1;
  //temp[15] = x_15_0 + x_15_1;

  int16x8_t temp0 = vec_add(a_0, a_1);
  int16x8_t temp1 = vec_add(a_2, a_3);

  //output[0] = (tran_low_t)fdct_round_shift(temp[0]);
  //output[1] = (tran_low_t)fdct_round_shift(temp[1]);
  //output[2] = (tran_low_t)fdct_round_shift(temp[2]);
  //output[3] = (tran_low_t)fdct_round_shift(temp[3]);

  //output[4] = (tran_low_t)fdct_round_shift(temp[4]);
  //output[5] = (tran_low_t)fdct_round_shift(temp[5]);
  //output[6] = (tran_low_t)fdct_round_shift(temp[6]);
  //output[7] = (tran_low_t)fdct_round_shift(temp[7]);

  //output[8] = (tran_low_t)fdct_round_shift(temp[8]);
  //output[9] = (tran_low_t)fdct_round_shift(temp[9]);
  //output[10] = (tran_low_t)fdct_round_shift(temp[10]);
  //output[11] = (tran_low_t)fdct_round_shift(temp[11]);

  //output[12] = (tran_low_t)fdct_round_shift(temp[12]);
  //output[13] = (tran_low_t)fdct_round_shift(temp[13]);
  //output[14] = (tran_low_t)fdct_round_shift(temp[14]);
  //output[15] = (tran_low_t)fdct_round_shift(temp[15]);
  //
  printf("---------------------------------------- BEFORE ROUDING\n");
  MATRIX_SI4_PRINT(temp0, temp1);

  int16x8_t temp2 = fdct_vector_round_shift(temp0);
  int16x8_t temp3 = fdct_vector_round_shift(temp1);

  //output[0] = (output[0] + 1) >> 2;
  //output[1] = (output[1] + 1) >> 2;
  //output[2] = (output[2] + 1) >> 2;
  //output[3] = (output[3] + 1) >> 2;

  //output[4] = (output[4] + 1) >> 2;
  //output[5] = (output[5] + 1) >> 2;
  //output[6] = (output[6] + 1) >> 2;
  //output[7] = (output[7] + 1) >> 2;

  //output[8] = (output[8] + 1) >> 2;
  //output[9] = (output[9] + 1) >> 2;
  //output[10] = (output[10] + 1) >> 2;
  //output[11] = (output[11] + 1) >> 2;

  //output[12] = (output[12] + 1) >> 2;
  //output[13] = (output[13] + 1) >> 2;
  //output[14] = (output[14] + 1) >> 2;
  //output[15] = (output[15] + 1) >> 2;

  int16x8_t one = vec_splat_s16(1);
  uint16x8_t two = vec_splat_u16(2);

  int16x8_t temp4 = vec_add(temp2, one);
  int16x8_t temp5 = vec_add(temp3, one);

  int16x8_t out1 = vec_sra(temp4, two);
  int16x8_t out2 = vec_sra(temp5, two);

  printf("---------------------------------------- BEFORE STORE\n");
  MATRIX_SI4_PRINT(out1, out2);

  store_tran_low(out1, 0, output);
  store_tran_low(out2, 0, output + (2* stride));

  printf("---------------------------------------- TRANSFORM OUTPUT\n");
  MATRIX_H4_PRINT(output);
}
