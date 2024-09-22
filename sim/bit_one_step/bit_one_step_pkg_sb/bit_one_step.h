/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: bit_one_step.h
 *
 * MATLAB Coder version            : 5.6
 * C/C++ source code generated on  : 22-Sep-2024 07:49:29
 */

#ifndef BIT_ONE_STEP_H
#define BIT_ONE_STEP_H

/* Include Files */
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
extern void bit_one_step(const float x0[21], float tau_applied[9],
                         const float unlock[9], float w_piv, boolean_T piv_flag,
                         float dt, unsigned short num_steps, float tau_max_piv,
                         float thet_pit_nom, const float x_flex0[104],
                         const float tau_flex[5], boolean_T flexure_flag,
                         boolean_T sb_flag, float y_true[21],
                         float y_flex[104]);

#ifdef __cplusplus
}
#endif

#endif
/*
 * File trailer for bit_one_step.h
 *
 * [EOF]
 */
