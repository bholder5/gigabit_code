/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: libbitonestep.h
 *
 * MATLAB Coder version            : 5.3
 * C/C++ source code generated on  : 16-Aug-2022 11:56:01
 */

#ifndef LIBBITONESTEP_H
#define LIBBITONESTEP_H

/* Include Files */
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
extern void bit_one_step(const real_T x0[21], const real_T tau_applied[9],
                         const real_T unlock[9], real_T w_piv,
                         boolean_T piv_flag, real_T dt, uint16_T num_steps,
                         real_T tau_max_piv, real_T thet_pit_nom,
                         real_T y_true[21]);

extern void compute_angular_velocity_C(const real_T x[18], real_T z_n[9][3],
                                       real_T omega[3]);

extern void compute_rotation_mat_C(real_T z_n[9][3], const real_T theta[9],
                                   real_T C[3][3]);

extern void libbitonestep_initialize(void);

extern void libbitonestep_terminate(void);

extern void rot2axis_C(real_T C[3][3], real_T v[3], real_T *phi);

#ifdef __cplusplus
}
#endif

#endif
/*
 * File trailer for libbitonestep.h
 *
 * [EOF]
 */
