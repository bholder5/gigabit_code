/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: libbitonestep.h
 *
 * MATLAB Coder version            : 5.6
 * C/C++ source code generated on  : 14-Mar-2024 11:31:39
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
extern void bit_one_step(const real32_T x0[21], real32_T tau_applied[9],
                         const real32_T unlock[9], real32_T w_piv,
                         boolean_T piv_flag, real32_T dt, uint16_T num_steps,
                         real32_T tau_max_piv, real32_T thet_pit_nom,
                         const real32_T x_flex0[104],
                         const real32_T tau_flex[5], boolean_T flexure_flag,
                         real32_T y_true[21], real32_T y_flex[104]);

extern void compute_angular_velocity_C(const real_T x[18], real_T z_n[9][3],
                                       real_T omega[3]);

extern void compute_angular_velocity_roll_C(const real_T x[18],
                                            real_T z_n[9][3], real_T omega[3]);

extern void compute_angular_velocity_yaw_C(const real_T x[18], real_T z_n[9][3],
                                           real_T omega[3]);

extern void compute_rotation_mat_C(real_T z_n[9][3], const real_T theta[9],
                                   real_T C[3][3]);

extern void compute_rotation_mat_roll_C(real_T z_n[9][3], const real_T theta[9],
                                        real_T C[3][3]);

extern void compute_rotation_mat_yaw_C(real_T z_n[9][3], const real_T theta[9],
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
