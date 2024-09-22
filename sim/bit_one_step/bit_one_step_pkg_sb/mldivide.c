/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: mldivide.c
 *
 * MATLAB Coder version            : 5.6
 * C/C++ source code generated on  : 22-Sep-2024 07:49:29
 */

/* Include Files */
#include "mldivide.h"
#include "rt_nonfinite.h"
#include <math.h>
#include <string.h>

/* Function Definitions */
/*
 * Arguments    : const float A[81]
 *                float B[9]
 * Return Type  : void
 */
void mldivide(const float A[81], float B[9])
{
  float b_A[81];
  float smax;
  int A_tmp;
  int a;
  int i;
  int j;
  int jA;
  int jp1j;
  int k;
  signed char ipiv[9];
  memcpy(&b_A[0], &A[0], 81U * sizeof(float));
  for (i = 0; i < 9; i++) {
    ipiv[i] = (signed char)(i + 1);
  }
  for (j = 0; j < 8; j++) {
    int b_tmp;
    int mmj_tmp;
    signed char i1;
    mmj_tmp = 7 - j;
    b_tmp = j * 10;
    jp1j = b_tmp + 2;
    jA = 9 - j;
    a = 0;
    smax = fabsf(b_A[b_tmp]);
    for (k = 2; k <= jA; k++) {
      float s;
      s = fabsf(b_A[(b_tmp + k) - 1]);
      if (s > smax) {
        a = k - 1;
        smax = s;
      }
    }
    if (b_A[b_tmp + a] != 0.0F) {
      if (a != 0) {
        jA = j + a;
        ipiv[j] = (signed char)(jA + 1);
        for (k = 0; k < 9; k++) {
          a = j + k * 9;
          smax = b_A[a];
          A_tmp = jA + k * 9;
          b_A[a] = b_A[A_tmp];
          b_A[A_tmp] = smax;
        }
      }
      i = (b_tmp - j) + 9;
      for (a = jp1j; a <= i; a++) {
        b_A[a - 1] /= b_A[b_tmp];
      }
    }
    jA = b_tmp;
    for (A_tmp = 0; A_tmp <= mmj_tmp; A_tmp++) {
      smax = b_A[(b_tmp + A_tmp * 9) + 9];
      if (smax != 0.0F) {
        i = jA + 11;
        a = (jA - j) + 18;
        for (jp1j = i; jp1j <= a; jp1j++) {
          b_A[jp1j - 1] += b_A[((b_tmp + jp1j) - jA) - 10] * -smax;
        }
      }
      jA += 9;
    }
    i1 = ipiv[j];
    if (i1 != j + 1) {
      smax = B[j];
      B[j] = B[i1 - 1];
      B[i1 - 1] = smax;
    }
  }
  for (k = 0; k < 9; k++) {
    jA = 9 * k;
    if (B[k] != 0.0F) {
      i = k + 2;
      for (a = i; a < 10; a++) {
        B[a - 1] -= B[k] * b_A[(a + jA) - 1];
      }
    }
  }
  for (k = 8; k >= 0; k--) {
    jA = 9 * k;
    smax = B[k];
    if (smax != 0.0F) {
      smax /= b_A[k + jA];
      B[k] = smax;
      for (a = 0; a < k; a++) {
        B[a] -= B[k] * b_A[a + jA];
      }
    }
  }
}

/*
 * File trailer for mldivide.c
 *
 * [EOF]
 */
