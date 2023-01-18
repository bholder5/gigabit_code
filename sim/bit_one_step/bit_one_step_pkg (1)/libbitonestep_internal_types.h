/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: libbitonestep_internal_types.h
 *
 * MATLAB Coder version            : 5.3
 * C/C++ source code generated on  : 16-Aug-2022 11:56:01
 */

#ifndef LIBBITONESTEP_INTERNAL_TYPES_H
#define LIBBITONESTEP_INTERNAL_TYPES_H

/* Include Files */
#include "libbitonestep_types.h"
#include "rtwtypes.h"

/* Type Definitions */
#ifndef typedef_rtRunTimeErrorInfo
#define typedef_rtRunTimeErrorInfo
typedef struct {
  int32_T lineNo;
  int32_T colNo;
  const char_T *fName;
  const char_T *pName;
} rtRunTimeErrorInfo;
#endif /* typedef_rtRunTimeErrorInfo */

#endif
/*
 * File trailer for libbitonestep_internal_types.h
 *
 * [EOF]
 */