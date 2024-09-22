/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: axis2rot.c
 *
 * MATLAB Coder version            : 5.6
 * C/C++ source code generated on  : 22-Sep-2024 07:49:29
 */

/* Include Files */
#include "axis2rot.h"
#include "rt_nonfinite.h"
#include <math.h>
#include <string.h>

/* Function Definitions */
/*
 * This function gives the rotation matric applied to other rotation
 *  matricies, not the vector (it is transpose of the rot mat applied to the
 *  vector.
 *
 * Arguments    : const double v[3]
 *                float phi
 *                double rot[9]
 * Return Type  : void
 */
void axis2rot(const double v[3], float phi, double rot[9])
{
  double b_sign;
  float cosa;
  float sina;
  int j;
  int k;
  cosa = cosf(phi);
  sina = sinf(phi);
  b_sign = 1.0;
  memset(&rot[0], 0, 9U * sizeof(double));
  for (k = 0; k < 3; k++) {
    int i;
    i = 2 - k;
    for (j = 0; j <= i; j++) {
      float mij;
      int b_j;
      b_j = k + j;
      mij = (1.0F - cosa) * (float)v[k] * (float)v[b_j];
      if (k == b_j) {
        mij += cosa;
        rot[k + 3 * b_j] = mij;
      } else {
        float rot_tmp;
        /* index is 3 - j - k for 0 indexed programming languages */
        rot_tmp = (float)b_sign * sina * (float)v[3 - (k + b_j)];
        rot[k + 3 * b_j] = mij + rot_tmp;
        rot[b_j + 3 * k] = mij - rot_tmp;
        b_sign = -b_sign;
      }
    }
  }
}

/*
 * File trailer for axis2rot.c
 *
 * [EOF]
 */
