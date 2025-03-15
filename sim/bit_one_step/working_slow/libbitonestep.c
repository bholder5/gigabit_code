/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: libbitonestep.c
 *
 * MATLAB Coder version            : 5.6
 * C/C++ source code generated on  : 23-Sep-2024 11:45:57
 */

/* Include Files */
#include "libbitonestep.h"
#include "libbitonestep_internal_types.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

/* Type Definitions */
#ifndef typedef_rtBoundsCheckInfo
#define typedef_rtBoundsCheckInfo
typedef struct {
  int32_T iFirst;
  int32_T iLast;
  int32_T lineNo;
  int32_T colNo;
  const char_T *aName;
  const char_T *fName;
  const char_T *pName;
  int32_T checkKind;
} rtBoundsCheckInfo;
#endif /* typedef_rtBoundsCheckInfo */

/* Function Declarations */
static void axis2rot(const real_T v[3], real_T phi, real_T rot[3][3]);

static void bit_one_step_anonFcn1(real_T c_n[9][3], real_T z_n[9][3],
                                  const real_T m_n[9], real_T b_r_n1_n[9][3],
                                  real_T m_w_n[9][6][6], real_T p_n[9][6],
                                  const real_T k_d[9], const real_T b_d[9],
                                  const real_T g0[3], const real_T unlock[9],
                                  const real_T hs_rw_max[3], real_T w_piv,
                                  boolean_T piv_flag, real_T tau_max_piv,
                                  real_T thet_pit_nom, const real_T y_true[21],
                                  const real_T tau_applied[9], real_T dw_piv,
                                  real_T varargout_1[24]);

static void c_rtErrorWithMessageID(const char_T *r, const char_T *aFcnName,
                                   int32_T aLineNum);

static void check_forloop_overflow_error(void);

static void d_rtErrorWithMessageID(const char_T *r, const char_T *aFcnName,
                                   int32_T aLineNum);

static int32_T div_nde_s32_floor(int32_T numerator);

static void mldivide(real_T c_A[9][9], real_T B[9]);

static void rtDynamicBoundsError(int32_T aIndexValue, int32_T aLoBound,
                                 int32_T aHiBound,
                                 const rtBoundsCheckInfo *aInfo);

static void rtErrorWithMessageID(const char_T *aFcnName, int32_T aLineNum);

static boolean_T rtIsNullOrEmptyString(const char_T *aString);

static void rtReportErrorLocation(const char_T *aFcnName, int32_T aLineNo);

/* Function Definitions */
/*
 * function rot = axis2rot( v, phi)
 *
 * This function gives the rotation matric applied to other rotation
 *  matricies, not the vector (it is transpose of the rot mat applied to the
 *  vector.
 *
 * Arguments    : const real_T v[3]
 *                real_T phi
 *                real_T rot[3][3]
 * Return Type  : void
 */
static void axis2rot(const real_T v[3], real_T phi, real_T rot[3][3])
{
  static rtBoundsCheckInfo b_emlrtBCI = {
      1,                                                       /* iFirst */
      3,                                                       /* iLast */
      19,                                                      /* lineNo */
      47,                                                      /* colNo */
      "v",                                                     /* aName */
      "axis2rot",                                              /* fName */
      "/home/bholder/bit-matlab-sim/Miscellaneous/axis2rot.m", /* pName */
      0                                                        /* checkKind */
  };
  static rtBoundsCheckInfo emlrtBCI = {
      1,                                                       /* iFirst */
      3,                                                       /* iLast */
      13,                                                      /* lineNo */
      35,                                                      /* colNo */
      "v",                                                     /* aName */
      "axis2rot",                                              /* fName */
      "/home/bholder/bit-matlab-sim/Miscellaneous/axis2rot.m", /* pName */
      0                                                        /* checkKind */
  };
  real_T b_sign;
  real_T cosa;
  real_T sina;
  int32_T b_i;
  int32_T j;
  int32_T k;
  /* 'axis2rot:5' cosa = cos(phi); */
  cosa = cos(phi);
  /* 'axis2rot:6' sina = sin(phi); */
  sina = sin(phi);
  /* 'axis2rot:8' sign = 1; */
  b_sign = 1.0;
  /* 'axis2rot:9' rot = (zeros(3,3)); */
  for (b_i = 0; b_i < 3; b_i++) {
    rot[b_i][0] = 0.0;
    rot[b_i][1] = 0.0;
    rot[b_i][2] = 0.0;
  }
  /* 'axis2rot:11' for k = 1:3 */
  for (k = 0; k < 3; k++) {
    /* 'axis2rot:12' for j = k:3 */
    b_i = 2 - k;
    for (j = 0; j <= b_i; j++) {
      real_T b_mij;
      int32_T b_j;
      b_j = k + j;
      /* 'axis2rot:13' mij = (1-cosa)*v(k)*v(j); */
      if ((b_j + 1) > 3) {
        rtDynamicBoundsError(b_j + 1, 1, 3, &emlrtBCI);
      }
      b_mij = ((1.0 - cosa) * v[k]) * v[b_j];
      /* 'axis2rot:14' if (k == j) */
      if (k == b_j) {
        /* 'axis2rot:15' mij = mij + cosa; */
        rot[b_j][k] = b_mij + cosa;
        /* 'axis2rot:16' rot(k,j) = mij; */
      } else {
        real_T rot_tmp;
        int32_T i1;
        /* 'axis2rot:17' else */
        /* index is 3 - j - k for 0 indexed programming languages */
        /* 'axis2rot:19' rot(k,j) = mij + (sign*sina*v((5-k-j)+1)); */
        i1 = 4 - (k + b_j);
        if ((i1 < 1) || (i1 > 3)) {
          rtDynamicBoundsError(i1, 1, 3, &b_emlrtBCI);
        }
        rot_tmp = (b_sign * sina) * v[i1 - 1];
        rot[b_j][k] = b_mij + rot_tmp;
        /* 'axis2rot:20' rot(j,k) = mij - (sign*sina*v((5-k-j)+1)); */
        rot[k][b_j] = b_mij - rot_tmp;
        /* 'axis2rot:21' sign = sign*-1; */
        b_sign = -b_sign;
      }
    }
  }
}

/*
 * @(y_true, tau_applied, dw_piv)
 *
 * Arguments    : real_T c_n[9][3]
 *                real_T z_n[9][3]
 *                const real_T m_n[9]
 *                real_T b_r_n1_n[9][3]
 *                real_T m_w_n[9][6][6]
 *                real_T p_n[9][6]
 *                const real_T k_d[9]
 *                const real_T b_d[9]
 *                const real_T g0[3]
 *                const real_T unlock[9]
 *                const real_T hs_rw_max[3]
 *                real_T w_piv
 *                boolean_T piv_flag
 *                real_T tau_max_piv
 *                real_T thet_pit_nom
 *                const real_T y_true[21]
 *                const real_T tau_applied[9]
 *                real_T dw_piv
 *                real_T varargout_1[24]
 * Return Type  : void
 */
static void bit_one_step_anonFcn1(real_T c_n[9][3], real_T z_n[9][3],
                                  const real_T m_n[9], real_T b_r_n1_n[9][3],
                                  real_T m_w_n[9][6][6], real_T p_n[9][6],
                                  const real_T k_d[9], const real_T b_d[9],
                                  const real_T g0[3], const real_T unlock[9],
                                  const real_T hs_rw_max[3], real_T w_piv,
                                  boolean_T piv_flag, real_T tau_max_piv,
                                  real_T thet_pit_nom, const real_T y_true[21],
                                  const real_T tau_applied[9], real_T dw_piv,
                                  real_T varargout_1[24])
{
  static rtBoundsCheckInfo
      b_emlrtBCI =
          {
              1,                     /* iFirst */
              9,                     /* iLast */
              37,                    /* lineNo */
              32,                    /* colNo */
              "T_n",                 /* aName */
              "compute_mass_matrix", /* fName */
              "/home/bholder/bit-matlab-sim/Plant_functions/"
              "compute_mass_matrix.m", /* pName */
              0                        /* checkKind */
          };
  static rtBoundsCheckInfo
      c_emlrtBCI =
          {
              1,                     /* iFirst */
              9,                     /* iLast */
              45,                    /* lineNo */
              32,                    /* colNo */
              "T_n",                 /* aName */
              "compute_mass_matrix", /* fName */
              "/home/bholder/bit-matlab-sim/Plant_functions/"
              "compute_mass_matrix.m", /* pName */
              0                        /* checkKind */
          };
  static rtBoundsCheckInfo
      d_emlrtBCI =
          {
              1,                     /* iFirst */
              9,                     /* iLast */
              52,                    /* lineNo */
              35,                    /* colNo */
              "p_n",                 /* aName */
              "compute_mass_matrix", /* fName */
              "/home/bholder/bit-matlab-sim/Plant_functions/"
              "compute_mass_matrix.m", /* pName */
              0                        /* checkKind */
          };
  static rtBoundsCheckInfo
      e_emlrtBCI =
          {
              1,                     /* iFirst */
              9,                     /* iLast */
              54,                    /* lineNo */
              38,                    /* colNo */
              "m_w_n",               /* aName */
              "compute_mass_matrix", /* fName */
              "/home/bholder/bit-matlab-sim/Plant_functions/"
              "compute_mass_matrix.m", /* pName */
              0                        /* checkKind */
          };
  static rtBoundsCheckInfo emlrtBCI = {
      1,                               /* iFirst */
      9,                               /* iLast */
      30,                              /* lineNo */
      29,                              /* colNo */
      "m_n",                           /* aName */
      "compute_potential_energy_term", /* fName */
      "/home/bholder/bit-matlab-sim/Plant_functions/"
      "compute_potential_energy_term.m", /* pName */
      0                                  /* checkKind */
  };
  static rtBoundsCheckInfo
      f_emlrtBCI =
          {
              1,                     /* iFirst */
              9,                     /* iLast */
              57,                    /* lineNo */
              20,                    /* colNo */
              "mass_mat",            /* aName */
              "compute_mass_matrix", /* fName */
              "/home/bholder/bit-matlab-sim/Plant_functions/"
              "compute_mass_matrix.m", /* pName */
              3                        /* checkKind */
          };
  static rtRunTimeErrorInfo emlrtRTEI = {
      109,   /* lineNo */
      "chol" /* fName */
  };
  real_T T_n[9][6][6];
  real_T C_n_rate[9][3][3];
  real_T c_C_n[9][3][3];
  real_T d_r_tmp[9][9];
  real_T T_nj[6][6];
  real_T b_T_n[6][6];
  real_T b_T_ni[6][6];
  real_T b_r_tmp[3][9];
  real_T c_r_tmp[3][9];
  real_T r_tmp[3][9];
  real_T s7[9][3];
  real_T Pot[9];
  real_T b_tau_applied[9];
  real_T dC_10[3][3];
  real_T dVdtheta_i[9];
  real_T d_C_n[3][3];
  real_T dtheta[9];
  real_T t1[9];
  real_T e_C_n[3];
  real_T M_ij;
  real_T d;
  real_T d1;
  real_T d2;
  real_T d3;
  real_T d_hs_idx_0;
  real_T d_hs_idx_1;
  real_T d_hs_idx_2;
  real_T int_err;
  real_T m_i;
  int32_T b_i;
  int32_T b_j;
  int32_T c_i;
  int32_T i1;
  int32_T i2;
  int32_T iac;
  int32_T idxA1j;
  int32_T idxAjj;
  int32_T idxAjjp1;
  int32_T info;
  int32_T j;
  int32_T jmax;
  int32_T k;
  int32_T n;
  boolean_T x[3];
  boolean_T exitg1;
  boolean_T y;
  /* 'bit_one_step:28' @(y_true, tau_applied, dw_piv) bit_propagator(y_true,
   * c_n, z_n, m_n, r_n1_n, m_w_n, p_n, ...  */
  /* 'bit_one_step:29'     k_d, b_d, g0, unlock, hs_rw_max, tau_applied, w_piv,
   * piv_flag, dw_piv, tau_max_piv, thet_pit_nom) */
  /* split the state */
  /* 'bit_propagator:6' theta = X(10:18); */
  /* 'bit_propagator:7' dtheta = X(1:9); */
  (void)memcpy(&b_tau_applied[0], &tau_applied[0], 9U * (sizeof(real_T)));
  (void)memcpy(&dtheta[0], &y_true[0], 9U * (sizeof(real_T)));
  /* 'bit_propagator:8' hs = X(19:21); */
  /* extract RW torque. */
  /* 'bit_propagator:11' tau_rw = tau_applied(7); */
  /* 'bit_propagator:12' tau_applied(7) = 0; */
  b_tau_applied[6] = 0.0;
  /*  set pivot speed in dtheta... */
  /* 'bit_propagator:15' if piv_flag == true */
  if (piv_flag) {
    /* 'bit_propagator:16' dtheta(6) = w_piv; */
    dtheta[5] = w_piv;
  }
  /*     %% */
  /* 'bit_propagator:20' Pot = compute_potential_energy_term(theta, c_n, z_n,
   * m_n, r_n1_n, g0); */
  /*  This function is implementing the potential energy term of the */
  /*  lagrangian. See section 3.1.1 of Romualdez Thesis, specifically */
  /*  Eq3.25-3.26. The rate is because of the partial derivitive of V wrt */
  /*  theta.  */
  /* Determine relative rotation matrices and rates for each joint */
  /* 'compute_potential_energy_term:8' ndof = length(theta); */
  /* rotation matricies */
  /* 'compute_potential_energy_term:11' C_n = zeros(3,3,ndof); */
  /* 'compute_potential_energy_term:12' C_n_rate = C_n; */
  /*  Potential energy */
  /* 'compute_potential_energy_term:14' Pot = zeros(ndof,1); */
  (void)memset(&Pot[0], 0, 9U * (sizeof(real_T)));
  /* 'compute_potential_energy_term:16' for n = 1:ndof */
  d_C_n[0][0] = 0.0;
  d_C_n[1][1] = 0.0;
  d_C_n[2][2] = 0.0;
  for (n = 0; n < 9; n++) {
    /* 'compute_potential_energy_term:17' C_n(:,:,n) = axis2rot(z_n(:,n),
     * theta(n)); */
    axis2rot(*((real_T(*)[3])(&z_n[n][0])), y_true[n + 9],
             *((real_T(*)[3][3])(&c_C_n[n][0][0])));
    /* Partial of Rot_n with respect to theta_n */
    /* 'compute_potential_energy_term:19' C_n_rate(:,:,n) = -1 * xmat(z_n(:,n))
     * * C_n(:,:,n); */
    /* 'xmat:2' mat = [        0, -vec(3),  vec(2); */
    /* 'xmat:3'             vec(3),         0, -vec(1); */
    /* 'xmat:4'            -vec(2),  vec(1),         0]; */
    d_C_n[1][0] = z_n[n][2];
    d_C_n[2][0] = -z_n[n][1];
    d_C_n[0][1] = -z_n[n][2];
    d_C_n[2][1] = z_n[n][0];
    d_C_n[0][2] = z_n[n][1];
    d_C_n[1][2] = -z_n[n][0];
    for (b_i = 0; b_i < 3; b_i++) {
      d = d_C_n[0][b_i];
      d1 = d_C_n[1][b_i];
      d2 = d_C_n[2][b_i];
      for (i1 = 0; i1 < 3; i1++) {
        C_n_rate[n][i1][b_i] =
            ((d * c_C_n[n][i1][0]) + (d1 * c_C_n[n][i1][1])) +
            (d2 * c_C_n[n][i1][2]);
      }
    }
  }
  /*  Compute potential energy term. First loop is cycling through each */
  /*  koints contribution. */
  /* 'compute_potential_energy_term:24' for n = 1:ndof */
  for (n = 0; n < 9; n++) {
    /* 'compute_potential_energy_term:25' m_i = 0; */
    m_i = 0.0;
    /* 'compute_potential_energy_term:26' dVdtheta_i = zeros(ndof,1); */
    (void)memset(&dVdtheta_i[0], 0, 9U * (sizeof(real_T)));
    /*  mass of remaining link (ex. 7 + 8 + 9) is totla mass at OF joint */
    /* 'compute_potential_energy_term:29' for q = n:ndof */
    b_i = 8 - n;
    for (jmax = 0; jmax <= b_i; jmax++) {
      /* 'compute_potential_energy_term:30' m_i = m_i + m_n(q); */
      i1 = (n + jmax) + 1;
      if (i1 > 9) {
        rtDynamicBoundsError(i1, 1, 9, &emlrtBCI);
      }
      m_i += m_n[i1 - 1];
    }
    /* 'compute_potential_energy_term:33' if (m_i ~= 0) */
    if (m_i != 0.0) {
      /* 'compute_potential_energy_term:34' t1 = zeros(ndof,1); */
      (void)memset(&t1[0], 0, 9U * (sizeof(real_T)));
      /* 'compute_potential_energy_term:35' t2 = t1; */
      /* Cycling through all joints up to joint n, the jint we are */
      /* currently calculating the contribution from */
      /*  */
      /* terms up to n-1 links.  */
      /* 'compute_potential_energy_term:41' for j = 1:n-1 */
      for (j = 0; j < n; j++) {
        /* 'compute_potential_energy_term:42' dC_10 = eye(3); */
        for (b_i = 0; b_i < 3; b_i++) {
          dC_10[b_i][0] = 0.0;
          dC_10[b_i][1] = 0.0;
          dC_10[b_i][2] = 0.0;
        }
        dC_10[0][0] = 1.0;
        dC_10[1][1] = 1.0;
        dC_10[2][2] = 1.0;
        /* 'compute_potential_energy_term:44' for k = 1:n-1 */
        for (k = 0; k < n; k++) {
          /* This multiplies the rotation matricies successively */
          /*  until the link j where the rate is inserted instead */
          /* 'compute_potential_energy_term:47' if (k == j) */
          if (k == j) {
            /* 'compute_potential_energy_term:48' dC_10 = C_n_rate(:,:,k) *
             * dC_10; */
            for (b_i = 0; b_i < 3; b_i++) {
              d = C_n_rate[k][0][b_i];
              d1 = C_n_rate[k][1][b_i];
              d2 = C_n_rate[k][2][b_i];
              for (i1 = 0; i1 < 3; i1++) {
                d_C_n[i1][b_i] = ((d * dC_10[i1][0]) + (d1 * dC_10[i1][1])) +
                                 (d2 * dC_10[i1][2]);
              }
            }
            for (b_i = 0; b_i < 3; b_i++) {
              dC_10[b_i][0] = d_C_n[b_i][0];
              dC_10[b_i][1] = d_C_n[b_i][1];
              dC_10[b_i][2] = d_C_n[b_i][2];
            }
          } else {
            /* 'compute_potential_energy_term:49' else */
            /* 'compute_potential_energy_term:50' dC_10 = C_n(:,:,k) * dC_10; */
            for (b_i = 0; b_i < 3; b_i++) {
              d = c_C_n[k][0][b_i];
              d1 = c_C_n[k][1][b_i];
              d2 = c_C_n[k][2][b_i];
              for (i1 = 0; i1 < 3; i1++) {
                d_C_n[i1][b_i] = ((d * dC_10[i1][0]) + (d1 * dC_10[i1][1])) +
                                 (d2 * dC_10[i1][2]);
              }
            }
            for (b_i = 0; b_i < 3; b_i++) {
              dC_10[b_i][0] = d_C_n[b_i][0];
              dC_10[b_i][1] = d_C_n[b_i][1];
              dC_10[b_i][2] = d_C_n[b_i][2];
            }
          }
        }
        /* 'compute_potential_energy_term:53' t1(j) = r_n1_n(:,n)' * dC_10 * g0;
         */
        M_ij = 0.0;
        d = b_r_n1_n[n][0];
        d1 = b_r_n1_n[n][1];
        d2 = b_r_n1_n[n][2];
        for (b_i = 0; b_i < 3; b_i++) {
          M_ij += (((d * dC_10[b_i][0]) + (d1 * dC_10[b_i][1])) +
                   (d2 * dC_10[b_i][2])) *
                  g0[b_i];
        }
        t1[j] = M_ij;
      }
      /* %% PE terms that go from 1:n */
      /* 'compute_potential_energy_term:56' for j = 1:n */
      d = c_n[n][0];
      d1 = c_n[n][1];
      d2 = c_n[n][2];
      for (j = 0; j <= n; j++) {
        /* 'compute_potential_energy_term:57' dC_10 = eye(3); */
        for (b_i = 0; b_i < 3; b_i++) {
          dC_10[b_i][0] = 0.0;
          dC_10[b_i][1] = 0.0;
          dC_10[b_i][2] = 0.0;
        }
        dC_10[0][0] = 1.0;
        dC_10[1][1] = 1.0;
        dC_10[2][2] = 1.0;
        /* 'compute_potential_energy_term:59' for k = 1:n */
        for (k = 0; k <= n; k++) {
          /* This multiplies the rotation matricies successively */
          /*  until the link j is reached in which case the rate */
          /*  is multiplied by the overall rotation */
          /* 'compute_potential_energy_term:63' if (k == j) */
          if (k == j) {
            /* 'compute_potential_energy_term:64' dC_10 = C_n_rate(:,:,k) *
             * dC_10; */
            for (b_i = 0; b_i < 3; b_i++) {
              d3 = C_n_rate[k][0][b_i];
              M_ij = C_n_rate[k][1][b_i];
              int_err = C_n_rate[k][2][b_i];
              for (i1 = 0; i1 < 3; i1++) {
                d_C_n[i1][b_i] = ((d3 * dC_10[i1][0]) + (M_ij * dC_10[i1][1])) +
                                 (int_err * dC_10[i1][2]);
              }
            }
            for (b_i = 0; b_i < 3; b_i++) {
              dC_10[b_i][0] = d_C_n[b_i][0];
              dC_10[b_i][1] = d_C_n[b_i][1];
              dC_10[b_i][2] = d_C_n[b_i][2];
            }
          } else {
            /* 'compute_potential_energy_term:65' else */
            /* 'compute_potential_energy_term:66' dC_10 = C_n(:,:,k) * dC_10; */
            for (b_i = 0; b_i < 3; b_i++) {
              d3 = c_C_n[k][0][b_i];
              M_ij = c_C_n[k][1][b_i];
              int_err = c_C_n[k][2][b_i];
              for (i1 = 0; i1 < 3; i1++) {
                d_C_n[i1][b_i] = ((d3 * dC_10[i1][0]) + (M_ij * dC_10[i1][1])) +
                                 (int_err * dC_10[i1][2]);
              }
            }
            for (b_i = 0; b_i < 3; b_i++) {
              dC_10[b_i][0] = d_C_n[b_i][0];
              dC_10[b_i][1] = d_C_n[b_i][1];
              dC_10[b_i][2] = d_C_n[b_i][2];
            }
          }
        }
        /* dot product */
        /* %%% made change here c_n(n) -> c_n(:,n) */
        /* 'compute_potential_energy_term:71' t2 = c_n(:,n)' * dC_10 * g0; */
        /* 'compute_potential_energy_term:73' dVdtheta_i(j) = -m_i*t1(j)-t2; */
        M_ij = 0.0;
        for (b_i = 0; b_i < 3; b_i++) {
          M_ij += (((d * dC_10[b_i][0]) + (d1 * dC_10[b_i][1])) +
                   (d2 * dC_10[b_i][2])) *
                  g0[b_i];
        }
        dVdtheta_i[j] = ((-m_i) * t1[j]) - M_ij;
      }
    }
    /* 'compute_potential_energy_term:76' Pot = Pot + dVdtheta_i; */
    for (b_i = 0; b_i < 9; b_i++) {
      Pot[b_i] += dVdtheta_i[b_i];
    }
  }
  /* 'bit_propagator:22' theta_spring = theta; */
  (void)memcpy(&t1[0], &y_true[9], 9U * (sizeof(real_T)));
  /* 'bit_propagator:23' theta_spring(9) = theta(9) - thet_pit_nom; */
  t1[8] = y_true[17] - thet_pit_nom;
  /* 'bit_propagator:25' spring = k_d.*theta_spring; */
  /* 'bit_propagator:26' damp = b_d.*dtheta; */
  /* place holder */
  /* 'bit_propagator:29' [R,r, d_hs] = RW_terms(theta, dtheta, z_n, hs, tau_rw,
   * hs_rw_max); */
  /* UNTITLED Summary of this function goes here */
  /*    Detailed explanation goes here */
  /* calculate the mapping matrix from dtheta to omega */
  /* 'RW_terms:7' s7 = zeros(3,9); */
  for (b_i = 0; b_i < 9; b_i++) {
    s7[b_i][0] = 0.0;
    s7[b_i][1] = 0.0;
    s7[b_i][2] = 0.0;
  }
  /* 'RW_terms:8' for i = 1:7 */
  for (c_i = 0; c_i < 7; c_i++) {
    real_T b_dVdtheta_i[9][3];
    /* 'RW_terms:9' Cn = axis2rot(z_n(:,i), theta(i)); */
    axis2rot(*((real_T(*)[3])(&z_n[c_i][0])), y_true[c_i + 9],
             *((real_T(*)[3][3])(&dVdtheta_i[0])));
    /* 'RW_terms:10' s7(:,i) = z_n(:,i); */
    s7[c_i][0] = z_n[c_i][0];
    s7[c_i][1] = z_n[c_i][1];
    s7[c_i][2] = z_n[c_i][2];
    /* 'RW_terms:11' s7 = Cn*s7; */
    for (b_i = 0; b_i < 3; b_i++) {
      d = dVdtheta_i[b_i];
      d1 = dVdtheta_i[b_i + 3];
      d2 = dVdtheta_i[b_i + 6];
      for (i1 = 0; i1 < 9; i1++) {
        b_dVdtheta_i[i1][b_i] =
            ((d * s7[i1][0]) + (d1 * s7[i1][1])) + (d2 * s7[i1][2]);
      }
    }
    for (b_i = 0; b_i < 9; b_i++) {
      s7[b_i][0] = b_dVdtheta_i[b_i][0];
      s7[b_i][1] = b_dVdtheta_i[b_i][1];
      s7[b_i][2] = b_dVdtheta_i[b_i][2];
    }
  }
  /* 'RW_terms:15' d_hs = tau_rw*z_n(:,7); */
  /* 'RW_terms:17' if hs(3) >= hs_rw_max */
  d_hs_idx_0 = tau_applied[6] * z_n[6][0];
  x[0] = (y_true[20] >= hs_rw_max[0]);
  d_hs_idx_1 = tau_applied[6] * z_n[6][1];
  x[1] = (y_true[20] >= hs_rw_max[1]);
  d_hs_idx_2 = tau_applied[6] * z_n[6][2];
  x[2] = (y_true[20] >= hs_rw_max[2]);
  y = true;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 3)) {
    if (!x[k]) {
      y = false;
      exitg1 = true;
    } else {
      k++;
    }
  }
  if (y) {
    /* 'RW_terms:18' if d_hs(3) > 0 */
    if (d_hs_idx_2 > 0.0) {
      /* 'RW_terms:19' d_hs(3) = 0; */
      d_hs_idx_2 = 0.0;
    }
  } else {
    x[0] = (y_true[20] <= (-hs_rw_max[0]));
    x[1] = (y_true[20] <= (-hs_rw_max[1]));
    x[2] = (y_true[20] <= (-hs_rw_max[2]));
    y = true;
    k = 0;
    exitg1 = false;
    while ((!exitg1) && (k < 3)) {
      if (!x[k]) {
        y = false;
        exitg1 = true;
      } else {
        k++;
      }
    }
    if (y && (d_hs_idx_2 < 0.0)) {
      /* 'RW_terms:21' elseif hs(3) <= -hs_rw_max */
      /* 'RW_terms:22' if d_hs(3) < 0 */
      /* 'RW_terms:23' d_hs(3) = 0; */
      d_hs_idx_2 = 0.0;
    }
  }
  /* 'RW_terms:27' r = (s7')*d_hs; */
  /* 'RW_terms:29' R = -(s7') * xmat(hs) * s7 * dtheta; */
  /* 'xmat:2' mat = [        0, -vec(3),  vec(2); */
  /* 'xmat:3'             vec(3),         0, -vec(1); */
  /* 'xmat:4'            -vec(2),  vec(1),         0]; */
  /*  rw_g2 = 2.353962297635081; */
  /*  rw_g1 = 0.034; */
  /*  w_piv = rw_g1*((-hs(3)/i_rw(3,3))-w_rw_nom) - (rw_g2*tau_rw); */
  /* calculate joint torques from gravity elasticity and damnping according */
  /* to eq 3.37 */
  /* 'bit_propagator:33' torques = tau_applied - (Pot + spring + damp + R + r);
   */
  for (b_i = 0; b_i < 3; b_i++) {
    for (i1 = 0; i1 < 9; i1++) {
      d = s7[i1][b_i];
      r_tmp[b_i][i1] = d;
      b_r_tmp[b_i][i1] = -d;
    }
  }
  d_C_n[0][0] = 0.0;
  d_C_n[1][0] = -y_true[20];
  d_C_n[2][0] = y_true[19];
  d_C_n[0][1] = y_true[20];
  d_C_n[1][1] = 0.0;
  d_C_n[2][1] = -y_true[18];
  d_C_n[0][2] = -y_true[19];
  d_C_n[1][2] = y_true[18];
  d_C_n[2][2] = 0.0;
  for (b_i = 0; b_i < 9; b_i++) {
    d = b_r_tmp[0][b_i];
    d1 = b_r_tmp[1][b_i];
    d2 = b_r_tmp[2][b_i];
    for (i1 = 0; i1 < 3; i1++) {
      c_r_tmp[i1][b_i] =
          ((d * d_C_n[i1][0]) + (d1 * d_C_n[i1][1])) + (d2 * d_C_n[i1][2]);
    }
    d = 0.0;
    d1 = c_r_tmp[0][b_i];
    d2 = c_r_tmp[1][b_i];
    d3 = c_r_tmp[2][b_i];
    for (i1 = 0; i1 < 9; i1++) {
      d += (((d1 * s7[i1][0]) + (d2 * s7[i1][1])) + (d3 * s7[i1][2])) *
           dtheta[i1];
    }
    b_tau_applied[b_i] -=
        (((Pot[b_i] + (k_d[b_i] * t1[b_i])) + (b_d[b_i] * dtheta[b_i])) + d) +
        (((r_tmp[0][b_i] * d_hs_idx_0) + (r_tmp[1][b_i] * d_hs_idx_1)) +
         (r_tmp[2][b_i] * d_hs_idx_2));
  }
  for (b_i = 0; b_i < 3; b_i++) {
    dC_10[b_i][0] = b_tau_applied[3 * b_i];
    dC_10[b_i][1] = b_tau_applied[(3 * b_i) + 1];
    dC_10[b_i][2] = b_tau_applied[(3 * b_i) + 2];
  }
  /* 'bit_propagator:35' M = compute_mass_matrix(theta, z_n, r_n1_n, m_w_n,
   * p_n); */
  /* Compute_Mass_Matrix computes the mass matrix of the 9 state system using */
  /* the angles theta to calculate the interbody transformation matrices and */
  /* the axis of rotations. */
  /*  Initialize mass matrix */
  /* 'compute_mass_matrix:7' mass_mat = (zeros(9,9)); */
  for (b_i = 0; b_i < 9; b_i++) {
    for (i1 = 0; i1 < 3; i1++) {
      c_C_n[b_i][i1][0] = 0.0;
      c_C_n[b_i][i1][1] = 0.0;
      c_C_n[b_i][i1][2] = 0.0;
    }
  }
  /*  memory to store the interbody transformations */
  /* 'compute_mass_matrix:10' T_n = (zeros(6,6,9)); */
  /* 'compute_mass_matrix:12' T_ni = (zeros(6,6)); */
  /* 'compute_mass_matrix:13' T_nj = T_ni; */
  /* 'compute_mass_matrix:15' J_ni = (zeros(6,1)); */
  /* 'compute_mass_matrix:16' J_nj = J_ni; */
  /*  equation 3.8: Interbody transformations for each frame */
  /* 'compute_mass_matrix:19' for i = 1:9 */
  d_C_n[0][0] = 0.0;
  d_C_n[1][1] = 0.0;
  d_C_n[2][2] = 0.0;
  for (c_i = 0; c_i < 9; c_i++) {
    real_T dv[3][3];
    /* 'compute_mass_matrix:20' C_n = axis2rot(z_n(:,i), theta(i)); */
    axis2rot(*((real_T(*)[3])(&z_n[c_i][0])), y_true[c_i + 9],
             *((real_T(*)[3][3])(&dVdtheta_i[0])));
    /* 'compute_mass_matrix:21' off_term = xmat(r_n1_n(:,i)); */
    /* 'xmat:2' mat = [        0, -vec(3),  vec(2); */
    /* 'xmat:3'             vec(3),         0, -vec(1); */
    /* 'xmat:4'            -vec(2),  vec(1),         0]; */
    /* 'compute_mass_matrix:22' off_term = -1*C_n*off_term; */
    /* 'compute_mass_matrix:23' T_n(:,:,i) = ([C_n, off_term; zeros(3,3), C_n]);
     */
    d = b_r_n1_n[c_i][2];
    d_C_n[1][0] = -d;
    d1 = b_r_n1_n[c_i][1];
    d_C_n[2][0] = d1;
    d_C_n[0][1] = d;
    d = b_r_n1_n[c_i][0];
    d_C_n[2][1] = -d;
    d_C_n[0][2] = -d1;
    d_C_n[1][2] = d;
    for (b_i = 0; b_i < 3; b_i++) {
      for (i1 = 0; i1 < 3; i1++) {
        dv[i1][b_i] = (((-dVdtheta_i[b_i]) * d_C_n[i1][0]) +
                       ((-dVdtheta_i[b_i + 3]) * d_C_n[i1][1])) +
                      ((-dVdtheta_i[b_i + 6]) * d_C_n[i1][2]);
        T_n[c_i][b_i][i1] = dVdtheta_i[i1 + (3 * b_i)];
      }
    }
    for (b_i = 0; b_i < 3; b_i++) {
      T_n[c_i][b_i + 3][0] = dv[b_i][0];
      T_n[c_i][b_i][3] = 0.0;
      T_n[c_i][b_i + 3][3] = dVdtheta_i[3 * b_i];
      T_n[c_i][b_i + 3][1] = dv[b_i][1];
      T_n[c_i][b_i][4] = 0.0;
      T_n[c_i][b_i + 3][4] = dVdtheta_i[(3 * b_i) + 1];
      T_n[c_i][b_i + 3][2] = dv[b_i][2];
      T_n[c_i][b_i][5] = 0.0;
      T_n[c_i][b_i + 3][5] = dVdtheta_i[(3 * b_i) + 2];
    }
  }
  /*  Generate the mass matrix */
  /*  Eq 3.12-3.13 */
  /* 'compute_mass_matrix:28' for i = 1:9 */
  for (c_i = 0; c_i < 9; c_i++) {
    /* 'compute_mass_matrix:29' for j = i:9 */
    b_i = 8 - c_i;
    for (j = 0; j <= b_i; j++) {
      b_j = (c_i + j) + 1;
      /* 'compute_mass_matrix:30' M_ij = 0; */
      M_ij = 0.0;
      /* 'compute_mass_matrix:31' for n = j:9 */
      i1 = 9 - b_j;
      if (i1 >= 0) {
        i2 = b_j - c_i;
      }
      for (n = 0; n <= i1; n++) {
        real_T T_ni[6];
        idxA1j = (b_j + n) - 1;
        /* 'compute_mass_matrix:32' T_ni = (eye(6)); */
        (void)memset(&b_T_ni[0][0], 0, 36U * (sizeof(real_T)));
        for (k = 0; k < 6; k++) {
          b_T_ni[k][k] = 1.0;
        }
        /* 'compute_mass_matrix:33' T_nj = T_ni; */
        (void)memcpy(&T_nj[0][0], &b_T_ni[0][0], 36U * (sizeof(real_T)));
        /* 'compute_mass_matrix:35' for k = i+1:j */
        for (k = 0; k <= (i2 - 2); k++) {
          /*                  tic */
          /* 'compute_mass_matrix:37' T_ni = T_n(:,:,k)*T_ni; */
          jmax = c_i + k;
          if ((jmax + 2) > 9) {
            rtDynamicBoundsError(jmax + 2, 1, 9, &b_emlrtBCI);
          }
          for (idxAjj = 0; idxAjj < 6; idxAjj++) {
            for (idxAjjp1 = 0; idxAjjp1 < 6; idxAjjp1++) {
              d = 0.0;
              for (info = 0; info < 6; info++) {
                d += T_n[jmax + 1][info][idxAjj] * b_T_ni[idxAjjp1][info];
              }
              b_T_n[idxAjjp1][idxAjj] = d;
            }
          }
          (void)memcpy(&b_T_ni[0][0], &b_T_n[0][0], 36U * (sizeof(real_T)));
          /*                  toc */
          /*                  tic */
          /*                  T_ni = mtimes(T_n(:,:,k),T_ni); */
          /*                  toc */
        }
        /* 'compute_mass_matrix:44' for k = j+1:n */
        idxAjj = idxA1j - b_j;
        for (k = 0; k <= idxAjj; k++) {
          /* 'compute_mass_matrix:45' T_nj = T_n(:,:,k) * T_nj; */
          jmax = b_j + k;
          if ((jmax + 1) > 9) {
            rtDynamicBoundsError(jmax + 1, 1, 9, &c_emlrtBCI);
          }
          for (idxAjjp1 = 0; idxAjjp1 < 6; idxAjjp1++) {
            for (info = 0; info < 6; info++) {
              d = 0.0;
              for (iac = 0; iac < 6; iac++) {
                d += T_n[jmax][iac][idxAjjp1] * T_nj[info][iac];
              }
              b_T_n[info][idxAjjp1] = d;
            }
          }
          (void)memcpy(&T_nj[0][0], &b_T_n[0][0], 36U * (sizeof(real_T)));
          /*                  T_nj = mtimes(T_n(:,:,k),T_nj); */
        }
        /* 'compute_mass_matrix:48' T_ni = T_nj*T_ni; */
        for (idxAjj = 0; idxAjj < 6; idxAjj++) {
          for (idxAjjp1 = 0; idxAjjp1 < 6; idxAjjp1++) {
            d = 0.0;
            for (info = 0; info < 6; info++) {
              d += T_nj[info][idxAjj] * b_T_ni[idxAjjp1][info];
            }
            b_T_n[idxAjjp1][idxAjj] = d;
          }
        }
        (void)memcpy(&b_T_ni[0][0], &b_T_n[0][0], 36U * (sizeof(real_T)));
        /*              T_ni = mtimes(T_nj,T_ni); */
        /* 'compute_mass_matrix:51' J_ni = (T_ni) * p_n(:,i); */
        /* 'compute_mass_matrix:52' J_nj = (T_nj) * p_n(:,j); */
        if (b_j > 9) {
          rtDynamicBoundsError(b_j, 1, 9, &d_emlrtBCI);
        }
        /* 'compute_mass_matrix:54' add = (J_ni' * m_w_n(:,:,n) * J_nj); */
        if ((idxA1j + 1) > 9) {
          rtDynamicBoundsError(idxA1j + 1, 1, 9, &e_emlrtBCI);
        }
        /* 'compute_mass_matrix:55' M_ij = M_ij + add; */
        for (idxAjj = 0; idxAjj < 6; idxAjj++) {
          d = 0.0;
          for (idxAjjp1 = 0; idxAjjp1 < 6; idxAjjp1++) {
            d += b_T_ni[idxAjjp1][idxAjj] * p_n[c_i][idxAjjp1];
          }
          T_ni[idxAjj] = d;
        }
        int_err = 0.0;
        for (idxAjj = 0; idxAjj < 6; idxAjj++) {
          d = 0.0;
          d1 = 0.0;
          for (idxAjjp1 = 0; idxAjjp1 < 6; idxAjjp1++) {
            d += T_ni[idxAjjp1] * m_w_n[idxA1j][idxAjj][idxAjjp1];
            d1 += T_nj[idxAjjp1][idxAjj] * p_n[b_j - 1][idxAjjp1];
          }
          int_err += d * d1;
        }
        M_ij += int_err;
      }
      /* 'compute_mass_matrix:57' mass_mat(i,j) = M_ij; */
      if (b_j > 9) {
        rtDynamicBoundsError(b_j, 1, 9, &f_emlrtBCI);
      }
      (*((real_T(*)[9][9])(&c_C_n[0][0][0])))[b_j - 1][c_i] = M_ij;
      /* 'compute_mass_matrix:58' if (i ~= j) */
      if ((c_i + 1) != b_j) {
        /* 'compute_mass_matrix:59' mass_mat(j,i) = M_ij; */
        (*((real_T(*)[9][9])(&c_C_n[0][0][0])))[c_i][b_j - 1] = M_ij;
      }
    }
  }
  /*  M = mass_mat_func(theta); */
  /*  M = mass_mat_func_gb(theta); */
  /* 'bit_propagator:40' M_decomp = chol(M); */
  for (b_i = 0; b_i < 9; b_i++) {
    for (i1 = 0; i1 < 3; i1++) {
      C_n_rate[b_i][i1][0] = c_C_n[b_i][i1][0];
      C_n_rate[b_i][i1][1] = c_C_n[b_i][i1][1];
      C_n_rate[b_i][i1][2] = c_C_n[b_i][i1][2];
    }
  }
  info = 0;
  j = 0;
  exitg1 = false;
  while ((!exitg1) && (j < 9)) {
    idxA1j = j * 9;
    idxAjj = idxA1j + j;
    M_ij = 0.0;
    if (j >= 1) {
      for (k = 0; k < j; k++) {
        int_err = (&C_n_rate[0][0][0])[idxA1j + k];
        M_ij += int_err * int_err;
      }
    }
    M_ij = (&C_n_rate[0][0][0])[idxAjj] - M_ij;
    if (M_ij > 0.0) {
      M_ij = sqrt(M_ij);
      (&C_n_rate[0][0][0])[idxAjj] = M_ij;
      if ((j + 1) < 9) {
        jmax = idxA1j + 10;
        idxAjjp1 = idxAjj + 10;
        if (j != 0) {
          b_i = (idxA1j + (9 * (7 - j))) + 10;
          for (iac = jmax; iac <= b_i; iac += 9) {
            int_err = 0.0;
            i1 = (iac + j) - 1;
            for (b_j = iac; b_j <= i1; b_j++) {
              int_err += (&C_n_rate[0][0][0])[b_j - 1] *
                         (&C_n_rate[0][0][0])[(idxA1j + b_j) - iac];
            }
            i1 = (idxAjj + (div_nde_s32_floor((iac - idxA1j) - 10) * 9)) + 9;
            (&C_n_rate[0][0][0])[i1] -= int_err;
          }
        }
        M_ij = 1.0 / M_ij;
        b_i = (idxAjj + (9 * (7 - j))) + 10;
        for (k = idxAjjp1; k <= b_i; k += 9) {
          (&C_n_rate[0][0][0])[k - 1] *= M_ij;
        }
      }
      j++;
    } else {
      (&C_n_rate[0][0][0])[idxAjj] = M_ij;
      info = j + 1;
      exitg1 = true;
    }
  }
  if (info == 0) {
    jmax = 7;
  } else {
    jmax = info - 3;
  }
  for (j = 0; j <= jmax; j++) {
    b_i = j + 2;
    for (c_i = b_i; c_i <= (jmax + 2); c_i++) {
      (*((real_T(*)[9][9])(&C_n_rate[0][0][0])))[j][c_i - 1] = 0.0;
    }
  }
  if (info != 0) {
    rtErrorWithMessageID(emlrtRTEI.fName, emlrtRTEI.lineNo);
  }
  /* 'bit_propagator:42' ddtheta = M_decomp\((M_decomp')\torques); */
  for (b_i = 0; b_i < 9; b_i++) {
    dVdtheta_i[b_i] = (&dC_10[0][0])[b_i];
    for (i1 = 0; i1 < 9; i1++) {
      d_r_tmp[b_i][i1] = (*((real_T(*)[9][9])(&C_n_rate[0][0][0])))[i1][b_i];
    }
  }
  mldivide(d_r_tmp, dVdtheta_i);
  mldivide(*((real_T(*)[9][9])(&C_n_rate[0][0][0])), dVdtheta_i);
  /* 'bit_propagator:43' ddtheta = ddtheta.*unlock; */
  for (b_i = 0; b_i < 9; b_i++) {
    dVdtheta_i[b_i] *= unlock[b_i];
  }
  /* 'bit_propagator:45' if piv_flag == true */
  if (piv_flag) {
    /* 'bit_propagator:46' prop_err = 10; */
    /* 'bit_propagator:47' int_err = 0; */
    /* 'bit_propagator:48' kp = 1; */
    /* 'bit_propagator:49' ki = 0.5; */
    /* 'bit_propagator:50' prop_err = dw_piv - ddtheta(6); */
    M_ij = dw_piv - dVdtheta_i[5];
    /* 'bit_propagator:51' int_err = int_err + prop_err; */
    int_err = M_ij;
    /* 'bit_propagator:52' tau_piv = torques(6); */
    m_i = (&dC_10[0][0])[5];
    /* 'bit_propagator:54' while abs(prop_err) > 1e-9 */
    exitg1 = false;
    while ((!exitg1) && (fabs(M_ij) > 1.0E-9)) {
      /* 'bit_propagator:56' tau_piv = tau_piv + ((kp*prop_err) + (ki*int_err));
       */
      m_i += M_ij + (0.5 * int_err);
      /* 'bit_propagator:57' if abs(tau_piv) > tau_max_piv */
      if (fabs(m_i) > tau_max_piv) {
        /* 'bit_propagator:58' tau_piv = sign(tau_piv) * tau_max_piv; */
        if (m_i < 0.0) {
          b_i = -1;
        } else {
          b_i = (m_i > 0.0) ? ((int32_T)1) : ((int32_T)0);
        }
        (&dC_10[0][0])[5] = ((real_T)b_i) * tau_max_piv;
        /* 'bit_propagator:59' torques(6) = tau_piv; */
        /* 'bit_propagator:61' ddtheta = M_decomp\((M_decomp')\torques); */
        for (b_i = 0; b_i < 9; b_i++) {
          for (i1 = 0; i1 < 9; i1++) {
            d_r_tmp[b_i][i1] =
                (*((real_T(*)[9][9])(&C_n_rate[0][0][0])))[i1][b_i];
          }
        }
        mldivide(d_r_tmp, &dC_10[0][0]);
        for (b_i = 0; b_i < 9; b_i++) {
          dVdtheta_i[b_i] = (&dC_10[0][0])[b_i];
        }
        mldivide(*((real_T(*)[9][9])(&C_n_rate[0][0][0])), dVdtheta_i);
        exitg1 = true;
      } else {
        /* 'bit_propagator:64' torques(6) = tau_piv; */
        (&dC_10[0][0])[5] = m_i;
        /* 'bit_propagator:66' ddtheta = M_decomp\((M_decomp')\torques); */
        for (b_i = 0; b_i < 9; b_i++) {
          dVdtheta_i[b_i] = (&dC_10[0][0])[b_i];
          for (i1 = 0; i1 < 9; i1++) {
            d_r_tmp[b_i][i1] =
                (*((real_T(*)[9][9])(&C_n_rate[0][0][0])))[i1][b_i];
          }
        }
        mldivide(d_r_tmp, dVdtheta_i);
        mldivide(*((real_T(*)[9][9])(&C_n_rate[0][0][0])), dVdtheta_i);
        /* 'bit_propagator:67' prop_err = dw_piv - ddtheta(6); */
        M_ij = dw_piv - dVdtheta_i[5];
        /* 'bit_propagator:68' int_err = int_err + prop_err; */
        int_err += M_ij;
      }
    }
  }
  /* 'bit_propagator:72' tau_gond = M(7:9,7:9) * ddtheta(7:9); */
  /*  tau_gond(1) = tau_rw */
  /* 'bit_propagator:75' Xdot = [ddtheta; dtheta; d_hs; tau_gond]; */
  d = dVdtheta_i[6];
  d1 = dVdtheta_i[7];
  d2 = dVdtheta_i[8];
  for (b_i = 0; b_i < 3; b_i++) {
    e_C_n[b_i] = (((*((real_T(*)[9][9])(&c_C_n[0][0][0])))[6][b_i + 6] * d) +
                  ((*((real_T(*)[9][9])(&c_C_n[0][0][0])))[7][b_i + 6] * d1)) +
                 ((*((real_T(*)[9][9])(&c_C_n[0][0][0])))[8][b_i + 6] * d2);
  }
  for (b_i = 0; b_i < 9; b_i++) {
    varargout_1[b_i] = dVdtheta_i[b_i];
    varargout_1[b_i + 9] = dtheta[b_i];
  }
  varargout_1[18] = d_hs_idx_0;
  varargout_1[21] = e_C_n[0];
  varargout_1[19] = d_hs_idx_1;
  varargout_1[22] = e_C_n[1];
  varargout_1[20] = d_hs_idx_2;
  varargout_1[23] = e_C_n[2];
}

/*
 * Arguments    : const char_T *r
 *                const char_T *aFcnName
 *                int32_T aLineNum
 * Return Type  : void
 */
static void c_rtErrorWithMessageID(const char_T *r, const char_T *aFcnName,
                                   int32_T aLineNum)
{
  (void)fprintf(
      stderr,
      "The loop variable of class %.*s might overflow on the last iteration of "
      "the for loop. This could lead to an infinite loop.",
      5, r);
  (void)fprintf(stderr, "\n");
  (void)fprintf(stderr, "Error in %s (line %d)", aFcnName, aLineNum);
  (void)fprintf(stderr, "\n");
  (void)fflush(stderr);
  abort();
}

/*
 * Arguments    : void
 * Return Type  : void
 */
static void check_forloop_overflow_error(void)
{
  static rtRunTimeErrorInfo emlrtRTEI = {
      87,                            /* lineNo */
      "check_forloop_overflow_error" /* fName */
  };
  c_rtErrorWithMessageID("int32", emlrtRTEI.fName, emlrtRTEI.lineNo);
}

/*
 * Arguments    : const char_T *r
 *                const char_T *aFcnName
 *                int32_T aLineNum
 * Return Type  : void
 */
static void d_rtErrorWithMessageID(const char_T *r, const char_T *aFcnName,
                                   int32_T aLineNum)
{
  (void)fprintf(stderr,
                "Domain error. To compute complex results from real x, use "
                "\'%.*s(complex(x))\'.",
                4, r);
  (void)fprintf(stderr, "\n");
  (void)fprintf(stderr, "Error in %s (line %d)", aFcnName, aLineNum);
  (void)fprintf(stderr, "\n");
  (void)fflush(stderr);
  abort();
}

/*
 * Arguments    : int32_T numerator
 * Return Type  : int32_T
 */
static int32_T div_nde_s32_floor(int32_T numerator)
{
  int32_T b_i;
  if ((numerator < 0) && ((numerator % 9) != 0)) {
    b_i = -1;
  } else {
    b_i = 0;
  }
  return (numerator / 9) + b_i;
}

/*
 * Arguments    : real_T c_A[9][9]
 *                real_T B[9]
 * Return Type  : void
 */
static void mldivide(real_T c_A[9][9], real_T B[9])
{
  real_T d_A[9][9];
  real_T smax;
  int32_T ijA;
  int32_T j;
  int32_T jA;
  int32_T jp1j;
  int32_T k;
  int32_T temp_tmp;
  int8_T ipiv[9];
  for (ijA = 0; ijA < 9; ijA++) {
    for (jA = 0; jA < 9; jA++) {
      d_A[ijA][jA] = c_A[ijA][jA];
    }
    ipiv[ijA] = (int8_T)(ijA + 1);
  }
  for (j = 0; j < 8; j++) {
    int32_T b_tmp;
    int32_T n;
    b_tmp = j * 10;
    jp1j = b_tmp + 2;
    n = 9 - j;
    jA = 0;
    smax = fabs((&d_A[0][0])[b_tmp]);
    for (k = 2; k <= n; k++) {
      real_T s;
      s = fabs((&d_A[0][0])[(b_tmp + k) - 1]);
      if (s > smax) {
        jA = k - 1;
        smax = s;
      }
    }
    if ((&d_A[0][0])[b_tmp + jA] != 0.0) {
      if (jA != 0) {
        jA += j;
        ipiv[j] = (int8_T)(jA + 1);
        for (k = 0; k < 9; k++) {
          temp_tmp = j + (k * 9);
          smax = (&d_A[0][0])[temp_tmp];
          ijA = jA + (k * 9);
          (&d_A[0][0])[temp_tmp] = (&d_A[0][0])[ijA];
          (&d_A[0][0])[ijA] = smax;
        }
      }
      ijA = (b_tmp - j) + 9;
      for (temp_tmp = jp1j; temp_tmp <= ijA; temp_tmp++) {
        (&d_A[0][0])[temp_tmp - 1] /= (&d_A[0][0])[b_tmp];
      }
    }
    n = 7 - j;
    jA = b_tmp + 11;
    for (jp1j = 0; jp1j <= n; jp1j++) {
      smax = (&d_A[0][0])[(b_tmp + (jp1j * 9)) + 9];
      if (smax != 0.0) {
        temp_tmp = (jA - j) + 7;
        if ((jA <= temp_tmp) && (temp_tmp > 2147483646)) {
          check_forloop_overflow_error();
        }
        for (ijA = jA; ijA <= temp_tmp; ijA++) {
          (&d_A[0][0])[ijA - 1] +=
              (&d_A[0][0])[((b_tmp + ijA) - jA) + 1] * (-smax);
        }
      }
      jA += 9;
    }
  }
  for (temp_tmp = 0; temp_tmp < 8; temp_tmp++) {
    int8_T b_i;
    b_i = ipiv[temp_tmp];
    if (((int32_T)b_i) != (temp_tmp + 1)) {
      smax = B[temp_tmp];
      B[temp_tmp] = B[b_i - 1];
      B[b_i - 1] = smax;
    }
  }
  for (k = 0; k < 9; k++) {
    jA = 9 * k;
    if (B[k] != 0.0) {
      ijA = k + 2;
      for (temp_tmp = ijA; temp_tmp < 10; temp_tmp++) {
        B[temp_tmp - 1] -= B[k] * (&d_A[0][0])[(temp_tmp + jA) - 1];
      }
    }
  }
  for (k = 8; k >= 0; k--) {
    jA = 9 * k;
    smax = B[k];
    if (smax != 0.0) {
      smax /= (&d_A[0][0])[k + jA];
      B[k] = smax;
      for (temp_tmp = 0; temp_tmp < k; temp_tmp++) {
        B[temp_tmp] -= B[k] * (&d_A[0][0])[temp_tmp + jA];
      }
    }
  }
}

/*
 * Arguments    : int32_T aIndexValue
 *                int32_T aLoBound
 *                int32_T aHiBound
 *                const rtBoundsCheckInfo *aInfo
 * Return Type  : void
 */
static void rtDynamicBoundsError(int32_T aIndexValue, int32_T aLoBound,
                                 int32_T aHiBound,
                                 const rtBoundsCheckInfo *aInfo)
{
  if (aLoBound == 0) {
    aIndexValue++;
    aLoBound = 1;
    aHiBound++;
  }
  if (rtIsNullOrEmptyString(aInfo->aName)) {
    (void)fprintf(stderr,
                  "Index exceeds array dimensions. Index value %d exceeds "
                  "valid range [%d-%d].",
                  aIndexValue, aLoBound, aHiBound);
    (void)fprintf(stderr, "\n");
    rtReportErrorLocation(aInfo->fName, aInfo->lineNo);
    (void)fflush(stderr);
    abort();
  } else {
    (void)fprintf(stderr,
                  "Index exceeds array dimensions. Index value %d exceeds "
                  "valid range [%d-%d] for array \'%s\'.",
                  aIndexValue, aLoBound, aHiBound, aInfo->aName);
    (void)fprintf(stderr, "\n");
    rtReportErrorLocation(aInfo->fName, aInfo->lineNo);
    (void)fflush(stderr);
    abort();
  }
}

/*
 * Arguments    : const char_T *aFcnName
 *                int32_T aLineNum
 * Return Type  : void
 */
static void rtErrorWithMessageID(const char_T *aFcnName, int32_T aLineNum)
{
  (void)fprintf(stderr, "Matrix must be positive definite.");
  (void)fprintf(stderr, "\n");
  (void)fprintf(stderr, "Error in %s (line %d)", aFcnName, aLineNum);
  (void)fprintf(stderr, "\n");
  (void)fflush(stderr);
  abort();
}

/*
 * Arguments    : const char_T *aString
 * Return Type  : boolean_T
 */
static boolean_T rtIsNullOrEmptyString(const char_T *aString)
{
  return (aString == NULL) || (((int8_T)(*aString)) == '\x00');
}

/*
 * Arguments    : const char_T *aFcnName
 *                int32_T aLineNo
 * Return Type  : void
 */
static void rtReportErrorLocation(const char_T *aFcnName, int32_T aLineNo)
{
  (void)fprintf(stderr, "Error in %s (line %d)", aFcnName, aLineNo);
  (void)fprintf(stderr, "\n");
}

/*
 * function [y_true, y_flex] = bit_one_step(x0, tau_applied, unlock, w_piv,
 * piv_flag,... dt, num_steps, tau_max_piv, thet_pit_nom, x_flex0, tau_flex,
 * flexure_flag, sb_flag)
 *
 * Run initialization script
 *
 * Arguments    : const real_T x0[21]
 *                real_T tau_applied[9]
 *                const real_T unlock[9]
 *                real_T w_piv
 *                boolean_T piv_flag
 *                real_T dt
 *                uint16_T num_steps
 *                real_T tau_max_piv
 *                real_T thet_pit_nom
 *                const real_T x_flex0[104]
 *                const real_T tau_flex[5]
 *                boolean_T flexure_flag
 *                boolean_T sb_flag
 *                real_T y_true[21]
 *                real_T y_flex[104]
 * Return Type  : void
 */
void bit_one_step(const real_T x0[21], real_T tau_applied[9],
                  const real_T unlock[9], real_T w_piv, boolean_T piv_flag,
                  real_T dt, uint16_T num_steps, real_T tau_max_piv,
                  real_T thet_pit_nom, const real_T x_flex0[104],
                  const real_T tau_flex[5], boolean_T flexure_flag,
                  boolean_T sb_flag, real_T y_true[21], real_T y_flex[104])
{
  static const real_T a_df[104][104] = {
      {0.0, -6850.0880077557886,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0},
      {1.0, -0.1655305169176462,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0},
      {0.0, 0.0, 0.0, -7962.6706577762307,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 1.0, -0.178467595465129,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, -12647.866496740349,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 1.0, -0.22492546762641491,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -13102.13472955605,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -0.22892911330414967,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -17096.533871483229,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -0.2615074291218758,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, -19067.27580068196,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 1.0, -0.27616861371764873,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -25085.498408823711,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -0.31676804389852031,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -26771.233428688109,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -0.32723834389440443,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -28195.442436642312,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -0.33582997148344168,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -34562.538682549493,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -0.37182005692296638,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -36233.2578458545,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -0.380700711036134,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -48525.3521198063,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -0.44056941391706389,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, -52667.483029779527,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0},
      {0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       1.0, -0.45898794332653031,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -64059.558804126144,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -0.5061997977246776,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -64572.011562112908,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -0.50822047011946658,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -68941.640846773153,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -0.525134804966394,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, -76967.907407482664,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0},
      {0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       1.0, -0.55486181129172218,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -81003.609059009119,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -0.56922265962981167,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, -90222.62662533854,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0},
      {0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       1.0, -0.60074163040474737,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -92586.065715088145,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -0.60855916956394029,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -95916.09088249068,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -0.619406460678255,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -101665.5083955781,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -0.6377005830186393,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, -108912.3115872177,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0},
      {0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       1.0, -0.66003730678566264,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -129428.24167601109,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -0.7195227353628546,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -152960.2015770768,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -0.78220253535021678,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, -154039.89905106369,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 1.0, -0.7849583404259457,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -163213.47678699941,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -0.80799375439912768,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -166657.26090653459,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -0.816473541289697,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, -172526.03224322811,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0},
      {0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       1.0, -0.83072506220344178,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -176768.18600326759,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -0.8408761763857211,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, -182183.82551471441,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0},
      {0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       1.0, -0.85365994521170874,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -211920.51399000589,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -0.92069650589106911,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -214126.80521207751,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -0.92547675327277135,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, -223222.93762258449,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 1.0, -0.94492949498379941,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -229895.57465282839,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -0.95894853804117852,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -240003.81612398339,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -0.979803686712769,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, -244064.29661276529,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0},
      {0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       1.0, -0.98805727893228,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, -253884.3362742863,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 1.0, -1.0077387285884889,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, -259570.1238903646,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 1.0, -1.018960497547112,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -274165.328172819,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.047215981873499,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, -283336.09120653989,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0},
      {0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       1.0, -1.0645864759737269,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -288225.92719689809,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0737335371439194,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, -297350.78483045712,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0},
      {0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       1.0, -1.09059760650839,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -301055.71221588022,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.097370880269529,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -316221.39474900393,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.1246713204292249,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, -317060.74671177007,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 1.0, -1.1261629486211491,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, -334528.25056612171,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0},
      {0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       1.0, -1.156768344252421,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0,
       0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -347548.48344749858,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.179064855633478,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -372964.2843963927,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.221416037877991,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -378028.03081502032,
       0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.229679683193994,
       0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, -388815.66217016062,
       0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 1.0, -1.247101699413742,
       0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -394431.31787162268},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.256075344669455}};
  static const real_T b_b_df[5][104] = {
      {0.0, 0.008320756923076926,   0.0, 0.0065853999999999956,
       0.0, -0.0045507713846153841, 0.0, 2.3593846153845192E-5,
       0.0, 2.7445538461513241E-5,  0.0, -0.00249863076923077,
       0.0, 1.609753846154618E-5,   0.0, -0.0022254400000000008,
       0.0, -0.002463941538461537,  0.0, 0.00223066153846153,
       0.0, -2.8643076923081969E-5, 0.0, 8.7846153846164549E-6,
       0.0, 1.9177230769230859E-5,  0.0, 0.003132307692307692,
       0.0, 0.00063538461538461828, 0.0, -3.5815384615377248E-5,
       0.0, 2.066923076923129E-6,   0.0, -1.4461538461705661E-7,
       0.0, -9.523076923055062E-6,  0.0, 5.5323076923078142E-5,
       0.0, -1.191384615384621E-5,  0.0, -3.2095384615384731E-5,
       0.0, -1.5276923076928691E-6, 0.0, -3.313846153849669E-7,
       0.0, -0.0002440030769230776, 0.0, 3.4061538461520358E-6,
       0.0, 3.0769230769294318E-7,  0.0, 5.15999999999943E-6,
       0.0, 4.5215384615379618E-5,  0.0, -8.59200000000264E-7,
       0.0, 6.1806153846194636E-6,  0.0, 9.969230769234908E-8,
       0.0, -2.6215384615387891E-5, 0.0, -5.0276923076923239E-5,
       0.0, 1.3415384615457819E-6,  0.0, 0.00027900000000000082,
       0.0, -6.5280000000000527E-5, 0.0, 0.0004809135384615377,
       0.0, -6.4572307692307713E-5, 0.0, 1.0578461538461589E-5,
       0.0, 3.8012000000001492E-5,  0.0, 5.8218461538461282E-5,
       0.0, -2.194153846153877E-5,  0.0, -0.00079818153846153721,
       0.0, 0.00013226769230769271, 0.0, 3.2464615384613912E-5,
       0.0, 8.8439384615396387E-5,  0.0, 0.0087015753846153947,
       0.0, 0.0001038030769230861,  0.0, -0.00037063692307692488,
       0.0, -0.0069612153846153836, 0.0, -0.00084236923076923},
      {0.0, -0.02028408936820085,    0.0, -0.01917573822362954,
       0.0, -0.050332117293214122,   0.0, 0.00017401236588194131,
       0.0, 7.5046009170591643E-6,   0.0, -0.0331079744840978,
       0.0, -3.5741380815567552E-5,  0.0, -0.047396984974355107,
       0.0, 0.07450339366742105,     0.0, 0.001964805512433017,
       0.0, -0.00019968240384938941, 0.0, 5.6516417671537786E-6,
       0.0, -0.00196204276914072,    0.0, -0.0063280641260071487,
       0.0, -0.0013335468115891491,  0.0, 0.00014313218770316169,
       0.0, -0.00023459760633672989, 0.0, -9.6483638659862521E-5,
       0.0, 0.0016858162049253791,   0.0, -0.01697703398445943,
       0.0, 0.0002157549203026301,   0.0, 0.0077435122981707972,
       0.0, -0.00017203518714793611, 0.0, 7.7685000194821223E-5,
       0.0, 0.01470508166875594,     0.0, 9.0665602784118777E-7,
       0.0, 0.00020593681053143761,  0.0, 0.0078081469720444052,
       0.0, 4.29262961931611E-5,     0.0, -3.6730613043539878E-5,
       0.0, -5.0244673051062043E-5,  0.0, 0.00048017905966448549,
       0.0, -9.4008861002419316E-5,  0.0, 0.00257742824112801,
       0.0, 8.6306138097192017E-6,   0.0, 0.004202317993055143,
       0.0, -0.001048493933533806,   0.0, -0.00177394200299101,
       0.0, 0.0011312292340101479,   0.0, -0.00020575356572559691,
       0.0, 0.00098064365273438261,  0.0, 0.001004444670946637,
       0.0, -0.00029764326832194073, 0.0, -0.01193821943767185,
       0.0, -0.0049886864834683478,  0.0, 0.00058980575887233346,
       0.0, -0.00021359913917956371, 0.0, -0.00065294045594938466,
       0.0, -0.00079220034843404286, 0.0, 0.0085644382159712467,
       0.0, -0.015748081391969641,   0.0, -0.00138099444723674},
      {0.0, 0.0096540544864419688,   0.0, -0.01163202554143886,
       0.0, -0.0112080844724667,     0.0, 3.8991103418184389E-5,
       0.0, 0.00018744858245097379,  0.0, 0.020152548174203429,
       0.0, 9.906953137522622E-5,    0.0, 0.02051594005585481,
       0.0, -0.070543142469148973,   0.0, -0.044373987083714719,
       0.0, 0.00057701933035178815,  0.0, 2.7359288201526761E-5,
       0.0, 0.0070067278074248491,   0.0, -0.0021995904773451308,
       0.0, -0.00045426829755807992, 0.0, 1.7266105640486431E-5,
       0.0, 0.00010334697399411461,  0.0, -2.3651989944297551E-6,
       0.0, -0.0002408698479152249,  0.0, 0.0017120551601996349,
       0.0, -1.2687209902627771E-5,  0.0, 0.0015495234145154039,
       0.0, 8.7508133047740959E-6,   0.0, 9.8497412856638663E-5,
       0.0, 0.03946244961935412,     0.0, 0.00040842752430981963,
       0.0, -9.9073891653510661E-5,  0.0, -0.0046964377497579783,
       0.0, 0.000108729033624576,    0.0, 5.6014471013789877E-5,
       0.0, -9.2546369420510632E-6,  0.0, 0.000384300455837554,
       0.0, 0.00036850405269307421,  0.0, 0.00319850945996064,
       0.0, 0.00072718931194515584,  0.0, -0.016490930816761191,
       0.0, 0.0031574213879833959,   0.0, 0.0048372685357203184,
       0.0, 0.0039254228428389224,   0.0, 0.00016628967527601281,
       0.0, -0.000407242887271996,   0.0, -0.00022437019258169381,
       0.0, -0.00046825830990990171, 0.0, -0.020657074438624759,
       0.0, -0.0069020276691464568,  0.0, 0.00077598782626265206,
       0.0, 0.0001014908067770024,   0.0, -0.001735771988292981,
       0.0, -0.000296634453801802,   0.0, 0.0058872119588950587,
       0.0, -0.0014413788705655861,  0.0, -0.00177066463838367},
      {0.0, -0.0040428003580101251,  0.0, 0.0074533822429511528,
       0.0, 0.00675336986255086,     0.0, -0.01214358152176909,
       0.0, -0.0184522418228579,     0.0, -0.0053156546824436861,
       0.0, -0.005784926297461539,   0.0, 0.01007161819628742,
       0.0, 0.0099811730068321558,   0.0, 0.016725323182585392,
       0.0, 0.058961561513053637,    0.0, 0.099212696332003261,
       0.0, -0.02273034167112653,    0.0, 0.0071072743908696511,
       0.0, -0.024456802788191359,   0.0, -0.0225282404897517,
       0.0, 0.060418491782594348,    0.0, 0.03091285972196119,
       0.0, -0.0020977979095101128,  0.0, 0.06998252935527037,
       0.0, -0.00062013179746458048, 0.0, 0.01957157968799891,
       0.0, 0.032313440511189263,    0.0, -0.0071153149194514509,
       0.0, 0.089621056139394709,    0.0, 0.0073898092281106222,
       0.0, -0.016483521469336911,   0.0, -0.048854396482386361,
       0.0, 0.0050727113573692329,   0.0, 0.027054124644752951,
       0.0, -0.0078139623198636311,  0.0, -0.0363151838560942,
       0.0, 0.01355254005566817,     0.0, -0.043992907011276658,
       0.0, 0.035938381876238032,    0.0, -0.15004597200136369,
       0.0, -0.01279734311026298,    0.0, 0.0057557493447683268,
       0.0, 0.0198531541439508,      0.0, 0.0022292579060057892,
       0.0, -0.0013548360989590251,  0.0, 0.00093487139263362558,
       0.0, 0.0084864464580716673,   0.0, -0.091530654463505251,
       0.0, -0.065284615397894624,   0.0, -0.0081798580013035145,
       0.0, 0.0263781227908255,      0.0, 0.00048200608246749322,
       0.0, 0.03940911522087126,     0.0, -0.17136697661704661,
       0.0, 0.02564680202616346,     0.0, -0.035870834732484209},
      {0.0, 0.0039324732961645709,  0.0, -0.007593121890897103,
       0.0, -0.0075416614088555586, 0.0, -0.012050109567904191,
       0.0, -0.018209846047884119,  0.0, 0.0056567649115213583,
       0.0, -0.0057411447961357581, 0.0, -0.010198305432373879,
       0.0, -0.0099500589090496978, 0.0, -0.015912022138171521,
       0.0, 0.059086734067436031,   0.0, 0.0992841302200189,
       0.0, 0.023166051028904549,   0.0, 0.0022107087942586289,
       0.0, -0.02455766460841945,   0.0, -0.021401089625040681,
       0.0, 0.0589128374195349,     0.0, 0.030095076570868708,
       0.0, 0.01097410851754513,    0.0, -0.070591404052476714,
       0.0, 0.0022084232892898648,  0.0, -0.019187685233566971,
       0.0, 0.03111064053180964,    0.0, -0.0079105713432141879,
       0.0, -0.09025196572903621,   0.0, 0.0055835836709738958,
       0.0, -0.014852297057966059,  0.0, 0.048945935752260028,
       0.0, 0.0048189169315320706,  0.0, 0.025597760264376671,
       0.0, -0.0076758314089119912, 0.0, 0.035289635636920808,
       0.0, 0.010525182605796379,   0.0, 0.041272275986769927,
       0.0, 0.026174867780662081,   0.0, 0.13343504659095989,
       0.0, -0.073044822987175323,  0.0, -0.0059794632994338674,
       0.0, -0.0220711779889558,    0.0, -0.0043628552308314743,
       0.0, 0.0199471501891688,     0.0, 0.02347865569282874,
       0.0, 0.01489152586798825,    0.0, 0.091990534816266345,
       0.0, 0.06344316825951378,    0.0, -0.016873668952046759,
       0.0, 0.02219990512445736,    0.0, -0.00076477779292602422,
       0.0, -0.00020938302908622,   0.0, 0.166599021410397,
       0.0, -0.0144830576027064,    0.0, -0.054790779983943542}};
  static const real_T c_i_n[9][9] = {
      {0.0, 0.0, 8.448E+6, 0.0, 0.0, 10.1, 2655.0, 14.0, 34.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 7.5448E+6, 0.0, 0.0, 8000.0, 3787.0, 21.0, 40.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 7.5448E+6, 0.0, 0.0, 8000.0, 3787.0, 30.0, 25.0}};
  static const real_T b_c_n[9][3] = {
      {0.0, 0.0, 0.0},   {0.0, 0.0, 0.0}, {0.0, 0.0, 35.0},
      {0.0, 0.0, -30.5}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0},
      {0.0, 0.0, -1.4},  {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
  static const real_T b_r_n1_n[9][3] = {
      {0.0, 0.0, 0.0},   {0.0, 0.0, 0.0},  {0.0, 0.0, 0.0},
      {0.0, 0.0, -61.0}, {0.0, 0.0, 0.0},  {0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0},   {0.0, 0.0, -1.4}, {0.0, 0.0, 0.0}};
  static const real_T c_n[9][3] = {
      {0.0, 0.0, 0.0},  {0.0, 0.0, 0.0}, {0.0, 0.0, -30.5},
      {0.0, 0.0, 0.0},  {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0},
      {0.0, 0.0, -1.4}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
  static const real_T b_k_d[9] = {0.0,
                                  0.0,
                                  0.0,
                                  0.0,
                                  0.0,
                                  34906.585039886588,
                                  0.0,
                                  303.09234741555724,
                                  555.66930359518835};
  static const real_T k_d[9] = {0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0017453292519943296,
                                0.0,
                                62.607671231740191,
                                62.607671231740191};
  static const real_T g0[3] = {0.0, 0.0, -9.72};
  static const real_T hs_rw_max[3] = {0.0, 0.0, 56.548667764616276};
  static const int32_T b_i_n[9][9] = {{0, 0, 5448000, 0, 0, 1, 246, 151, 213},
                                      {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                      {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                      {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                      {0, 0, 5448000, 0, 0, 1, 455, 405, 134},
                                      {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                      {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                      {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                      {0, 0, 5448000, 0, 0, 1, 408, 339, 244}};
  static const int32_T iv[9] = {0, 0, 100000, 0, 0, 1, 350, 73, 150};
  static const int32_T m_n[9] = {0, 0, 100000, 0, 0, 1, 350, 73, 150};
  static const int16_T b_m_n[9] = {0, 0, 10000, 0, 0, 1, 1850, 60, 200};
  static const int16_T iv1[9] = {0, 0, 10000, 0, 0, 1, 1850, 60, 200};
  static const int8_T p_n[9][6] = {
      {0, 0, 0, 0, 0, 1}, {0, 0, 0, 1, 0, 0}, {0, 0, 0, 0, 1, 0},
      {0, 0, 0, 1, 0, 0}, {0, 0, 0, 0, 1, 0}, {0, 0, 0, 0, 0, 1},
      {0, 0, 0, 0, 0, 1}, {0, 0, 0, 1, 0, 0}, {0, 0, 0, 0, 1, 0}};
  static const int8_T z_n[9][3] = {{0, 0, 1}, {1, 0, 0}, {0, 1, 0},
                                   {1, 0, 0}, {0, 1, 0}, {0, 0, 1},
                                   {0, 0, 1}, {1, 0, 0}, {0, 1, 0}};
  static const int8_T a[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  real_T m_w_n[9][6][6];
  real_T sys_workspace_p_n[9][6];
  real_T sys_workspace_c_n[9][3];
  real_T sys_workspace_z_n[9][3];
  real_T i_n[9];
  real_T sys_workspace_b_d[9];
  real_T sys_workspace_k_d[9];
  real_T b_tau[5];
  real_T c_tau[5];
  real_T d_tau[5];
  real_T tau[5];
  real_T tau_app_flex_idx_0;
  real_T tau_app_flex_idx_1;
  real_T tau_app_flex_idx_2;
  int32_T b_i;
  int32_T i1;
  int32_T k;
  int32_T step;
  /* 'bit_one_step:4' if sb_flag */
  if (sb_flag) {
    real_T offterm[3][3];
    /* 'bit_one_step:5' [ndof, g0, r_n1_n, z_n, p_n, m_n, c_n, ... */
    /* 'bit_one_step:6'         i_n, m_w_n,  i_rw, bear_k_cst, bear_c_cst, k_d,
     * b_d, ... */
    /* 'bit_one_step:7'         w_rw_max, w_rw_nom, hs_rw, hs_rw_max, a_flex,
     * b_flex, a_df, b_df] = init_func_sb(); */
    /* %% This file initializes are parameters for simulating BIT (9 DOF) */
    /* 'init_func_sb:6' ndof = 9; */
    /* 'init_func_sb:7' g0 = [0; 0; -9.72]; */
    /* 'init_func_sb:9' tel_offset = [0; 40.0*pi/180; 0]; */
    /* each column is vector (F4 is -61m from F3) */
    /* 'init_func_sb:12' r_n1_n = [0,0,0,  0,0,0,0,   0,0; */
    /* 'init_func_sb:13'           0,0,0,  0,0,0,0,   0,0; */
    /* 'init_func_sb:14'           0,0,0,-61,0,0,0,-1.4,0]; */
    /*  r_n1_n = [0,0,0,  0,0,0,0,   0,0; */
    /*            0,0,0,  0,0,0,0,   0,0; */
    /*            0,0,0,-61,0,0,-1.4,0,0]; */
    /* each column is vector */
    /* 'init_func_sb:19' z_n = [0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0; */
    /* 'init_func_sb:20'        0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0; */
    /* 'init_func_sb:21'        1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0]; */
    /* 'init_func_sb:23' p_n = [zeros(3,9); z_n]; */
    /*  was 350 */
    /*  m_n = [0.0,0.0,100000.0,0.0,0.0,1.0,350.0,73.0,150.0]; */
    /* 'init_func_sb:26' m_n = [0.0,0.0,10000.0,0.0,0.0,1.0,1850.0,60.0,200.0];
     */
    /* each column is vector (COM of B3 (flight train) is 30.5m along z */
    /* 'init_func_sb:29' c_n = [0.,0.,  0.0,0.,0., 0.0,0.,0.,0.; */
    /* 'init_func_sb:30'        0.,0.,  0.0,0.,0., 0.0,0.,0.,0.; */
    /* 'init_func_sb:31'        0.,0.,35.0,-30.5*1,0.,0.,-1.4,0,0.]; */
    /* each row is row wise matric % Updated gondola on april 13 according to */
    /* michaels model */
    /* 'init_func_sb:34' yaw_sc = 5; */
    /*  i_n = [ [  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0], */
    /*          [  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0], */
    /*          [5448000.0,  0.0,  0.0,  0.0, 5448000.0,  0.0,  0.0,  0.0,
     * 5448000.0], */
    /*          [  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0], */
    /*          [  1.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,   1.0], */
    /*          [  246*yaw_sc,    0,    0,    0,  455*yaw_sc,    0,    0,    0,
     * 408*yaw_sc], */
    /*          [ 151,  0.0,  0.0,  0.0, 405,  0.0,  0.0,  0.0,  339]/1.1, */
    /*          [   213,  0.0,  0.0,  0.0, 134.0,  0.0, 0, 0, 244]/1.1]; */
    /*  ball_i_z = 3744800.0; */
    /*  ball_i_x = 0.01*(.158/.1)*ball_i_z */
    /*   */
    /*  i_n = [ [  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0], */
    /*          [  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0], */
    /*          [ball_i_z,  0.0,  0.0,  0.0, ball_i_x,  0.0,  0.0,  0.0,
     * ball_i_x], */
    /*          [  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0], */
    /*          [  1.0,  0.0,  0.0,  0.0,  50000.0,  0.0,  0.0,  0.0, 50000.0],
     */
    /*          [  3778,    0,    0,    0,  1.0*3787,    0,    0,
     * 0,   1.0*2655], */
    /*          [ 14,  0.0,  0.0,  0.0, 21,  0.0,  0.0,  0.0,  30], */
    /*          [   34,  0.0,  0.0,  0.0, 40,  0.0, 0, 0, 25]]; */
    /* 'init_func_sb:57' i_n = [ [  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     * 0.0,   0.0], */
    /* 'init_func_sb:58'         [  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     * 0.0,   0.0], */
    /* 'init_func_sb:59'         [8448000.0,  0.0,  0.0,  0.0, 7544800.0,  0.0,
     * 0.0,  0.0, 7544800.0], */
    /* 'init_func_sb:60'         [  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     * 0.0,   0.0], */
    /* 'init_func_sb:61'         [  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     * 0.0,   0.0], */
    /* 'init_func_sb:62'         [  10.1,  0.0,  0.0,  0.0,  8000.0,  0.0,  0.0,
     * 0.0,   8000.0], */
    /* 'init_func_sb:63'         [  2655,    0,    0,    0,  3787,    0,    0,
     * 0,   3787], */
    /* 'init_func_sb:64'         [ 14,  0.0,  0.0,  0.0, 21,  0.0,  0.0,  0.0,
     * 30], */
    /* 'init_func_sb:65'         [   34,  0.0,  0.0,  0.0, 40,  0.0, 0, 0, 25]];
     */
    /* 'init_func_sb:67' m_w_n = zeros(6,6,9); */
    /* 'init_func_sb:68' for k = 1:9 */
    /* 'init_func_sb:77' i_rw = reshape([2.5,0.,0.,0.,2.5,0.,0.,0.,4.5], 3, 3);
     */
    /*  bear_k_cst = 2.0*0.9486*(0.0254*4.44822162)*180.0/pi; */
    /* 'init_func_sb:79' bear_k_cst = 9.67 * 0.113 * 180 / pi; */
    /* SB spring constant */
    /* 'init_func_sb:82' bear_k_cst_r = 6*7.8023 * 0.113 * 180 / pi; */
    /* 'init_func_sb:83' bear_k_cst_p = 11*7.8023 * 0.113 * 180 / pi; */
    /* 'init_func_sb:84' bear_c_cst = 0.; */
    /* 'init_func_sb:86' k_d =
     * [0.,0.,0.,0.,0.,1*2000000*pi/180.0,0.,bear_k_cst_r,bear_k_cst_p]'; */
    /* 'init_func_sb:87' b_d = [0.,0.,0.,0.,0.,0.,0.,0.,0.]'; */
    /* 'init_func_sb:89' w_rw_max = 4.0*pi; */
    /* 'init_func_sb:90' w_rw_nom = 2*pi; */
    /* 'init_func_sb:91' hs_rw = i_rw * w_rw_nom * z_n(:,7); */
    /* 'init_func_sb:92' hs_rw_max = i_rw * w_rw_max * z_n(:,7); */
    /*  theta_0 =
     * [0.,0.4*pi/180.0,0.4*pi/180.0,0.1*pi/180.0,0.1*pi/180.0,0.1*pi/180.0,0.,0.1,-40*pi/180]';
     */
    /* 'init_func_sb:96' theta_0 = [0,0,0,0,0,0,0,0,0,0]'; */
    /* setting IC from past sim */
    /*  y0 = [0.017552353814854, -0.002156992032555, -0.002273627285241, ... */
    /*      -0.004091940730352,  -0.002796089196615,   0.019674817779806,... */
    /*      -0.017606183923045,                   0,                   0, ... */
    /*       0.207860712172010,  -0.003878840466313,  -0.004340266988222, ... */
    /*      -0.001098037684871,  -0.001085183886166,  -0.001924742862772, ... */
    /*       2.937417436471931,                   0,                   0, ... */
    /*                       0,                   0, 28.274274172758336]'; */
    /* 'init_func_sb:106' theta_des =
     * [0.,0.4*pi/180.0,0.4*pi/180.0,0.1*pi/180.0,0.1*pi/180.0,0.1*pi/180.0,0.,0.1,40*pi/180]';
     */
    /* 'init_func_sb:107' d_theta_dt_0 = [0.,0.,0.,0.,0.,0.,0.,0.,0.]'; */
    /* 'init_func_sb:109' unlock = [1.0,1.0,1.0,1.0,1.0,1.0,1,1,1]'; */
    /* 'init_func_sb:111' a_flex = a_f_func; */
    /* 'init_func_sb:112' b_flex = b_f_func(); */
    /* 'init_func_sb:113' a_df = a_mf_func(); */
    /* 'init_func_sb:114' b_df = b_mf_func(); */
    offterm[0][0] = 0.0;
    offterm[1][1] = 0.0;
    offterm[2][2] = 0.0;
    for (k = 0; k < 9; k++) {
      int16_T i2;
      /* 'init_func_sb:69' mass = eye(3) * m_n(k); */
      /* 'init_func_sb:70' i_k = reshape(i_n(k,:), 3, 3); */
      /* 'init_func_sb:71' offterm = xmat(c_n(:,k)); */
      tau_app_flex_idx_0 = b_c_n[k][0];
      tau_app_flex_idx_1 = b_c_n[k][1];
      tau_app_flex_idx_2 = b_c_n[k][2];
      /* 'xmat:2' mat = [        0, -vec(3),  vec(2); */
      /* 'xmat:3'             vec(3),         0, -vec(1); */
      /* 'xmat:4'            -vec(2),  vec(1),         0]; */
      offterm[1][0] = -tau_app_flex_idx_2;
      offterm[2][0] = tau_app_flex_idx_1;
      offterm[0][1] = tau_app_flex_idx_2;
      offterm[2][1] = -tau_app_flex_idx_0;
      offterm[0][2] = -tau_app_flex_idx_1;
      offterm[1][2] = tau_app_flex_idx_0;
      /* 'init_func_sb:73' m_w_n_i = [mass, -offterm; offterm, i_k]; */
      for (b_i = 0; b_i < 9; b_i++) {
        i_n[b_i] = c_i_n[b_i][k];
      }
      /* 'init_func_sb:74' m_w_n(:,:,k) = m_w_n_i; */
      i2 = iv1[k];
      for (b_i = 0; b_i < 3; b_i++) {
        m_w_n[k][b_i][0] =
            (real_T)((int32_T)(((int32_T)a[b_i][0]) * ((int32_T)i2)));
        tau_app_flex_idx_0 = offterm[b_i][0];
        m_w_n[k][b_i + 3][0] = -tau_app_flex_idx_0;
        m_w_n[k][b_i][3] = tau_app_flex_idx_0;
        m_w_n[k][b_i + 3][3] = i_n[3 * b_i];
        m_w_n[k][b_i][1] =
            (real_T)((int32_T)(((int32_T)a[b_i][1]) * ((int32_T)i2)));
        tau_app_flex_idx_0 = offterm[b_i][1];
        m_w_n[k][b_i + 3][1] = -tau_app_flex_idx_0;
        m_w_n[k][b_i][4] = tau_app_flex_idx_0;
        m_w_n[k][b_i + 3][4] = i_n[(3 * b_i) + 1];
        m_w_n[k][b_i][2] =
            (real_T)((int32_T)(((int32_T)a[b_i][2]) * ((int32_T)i2)));
        tau_app_flex_idx_0 = offterm[b_i][2];
        m_w_n[k][b_i + 3][2] = -tau_app_flex_idx_0;
        m_w_n[k][b_i][5] = tau_app_flex_idx_0;
        m_w_n[k][b_i + 3][5] = i_n[(3 * b_i) + 2];
        sys_workspace_z_n[k][b_i] = (real_T)z_n[k][b_i];
      }
      for (b_i = 0; b_i < 6; b_i++) {
        sys_workspace_p_n[k][b_i] = (real_T)p_n[k][b_i];
      }
    }
    for (k = 0; k < 9; k++) {
      i_n[k] = (real_T)b_m_n[k];
      sys_workspace_c_n[k][0] = b_c_n[k][0];
      sys_workspace_c_n[k][1] = b_c_n[k][1];
      sys_workspace_c_n[k][2] = b_c_n[k][2];
      sys_workspace_k_d[k] = b_k_d[k];
      sys_workspace_b_d[k] = 0.0;
    }
  } else {
    real_T offterm[3][3];
    /* 'bit_one_step:8' else */
    /* 'bit_one_step:9' [ndof, g0, r_n1_n, z_n, p_n, m_n, c_n, ... */
    /* 'bit_one_step:10'         i_n, m_w_n,  i_rw, bear_k_cst, bear_c_cst, k_d,
     * b_d, ... */
    /* 'bit_one_step:11'         w_rw_max, w_rw_nom, hs_rw, hs_rw_max, a_flex,
     * b_flex, a_df, b_df] = init_func(); */
    /* %% This file initializes are parameters for simulating BIT (9 DOF) */
    /* 'init_func:6' ndof = 9; */
    /* 'init_func:7' g0 = [0; 0; -9.72]; */
    /* 'init_func:9' tel_offset = [0; 40.0*pi/180; 0]; */
    /* each column is vector (F4 is -61m from F3) */
    /* 'init_func:12' r_n1_n = [0,0,0,  0,0,0,0,   0,0; */
    /* 'init_func:13'           0,0,0,  0,0,0,0,   0,0; */
    /* 'init_func:14'           0,0,0,-61,0,0,0,-1.4,0]; */
    /* each column is vector */
    /* 'init_func:16' z_n = [0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0; */
    /* 'init_func:17'        0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0; */
    /* 'init_func:18'        1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0]; */
    /* 'init_func:20' p_n = [zeros(3,9); z_n]; */
    /*  was 350 */
    /* 'init_func:22' m_n = [0.0,0.0,100000.0,0.0,0.0,1.0,350.0,73.0,150.0]; */
    /* each column is vector (COM of B3 (flight train) is 30.5m along z */
    /* 'init_func:26' c_n = [0.,0.,  0.0,0.,0., 0.0,0.,0.,0.; */
    /* 'init_func:27'        0.,0.,  0.0,0.,0., 0.0,0.,0.,0.; */
    /* 'init_func:28'       0,0.,-30.5,0.,0.,0.,-1.4,0.,0.]; */
    /* each row is row wise matric % Updated gondola on april 13 according to */
    /* michaels model */
    /* 'init_func:31' yaw_sc = 1; */
    /* 'init_func:32' i_n = [ [  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     * 0.0], */
    /* 'init_func:33'         [  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     * 0.0], */
    /* 'init_func:34'         [5448000.0,  0.0,  0.0,  0.0, 5448000.0,  0.0,
     * 0.0,  0.0, 5448000.0], */
    /* 'init_func:35'         [  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     * 0.0], */
    /* 'init_func:36'         [  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     * 0.0], */
    /* 'init_func:37'         [  1.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,
     * 0.0,   1.0], */
    /* 'init_func:38'         [  246*yaw_sc,    0,    0,    0,  455*yaw_sc, 0,
     * 0,    0,   408*yaw_sc], */
    /* 'init_func:39'         [ 151,  0.0,  0.0,  0.0, 405,  0.0,  0.0,  0.0,
     * 339], */
    /* 'init_func:40'         [   213,  0.0,  0.0,  0.0, 134.0,  0.0, 0, 0,
     * 244]]; */
    /* 'init_func:44' m_w_n = zeros(6,6,9); */
    /* 'init_func:45' for k = 1:9 */
    /* 'init_func:54' i_rw = reshape([2.5,0.,0.,0.,2.5,0.,0.,0.,4.5], 3, 3); */
    /*  bear_k_cst = 2.0*0.9486*(0.0254*4.44822162)*180.0/pi; */
    /* 'init_func:56' bear_k_cst = 9.67 * 0.113 * 180 / pi; */
    /* 'init_func:58' bear_c_cst = 0.; */
    /* 'init_func:60' k_d =
     * [0.,0.,0.,0.,0.,0.1*pi/180.0,0.,bear_k_cst,bear_k_cst]'; */
    /* 'init_func:61' b_d = [0.,0.,0.,0.,0.,0.,0.,0.,0.]'; */
    /* 'init_func:63' w_rw_max = 4.0*pi; */
    /* 'init_func:64' w_rw_nom = 2*pi; */
    /* 'init_func:65' hs_rw = i_rw * w_rw_nom * z_n(:,7); */
    /* 'init_func:66' hs_rw_max = i_rw * w_rw_max * z_n(:,7); */
    /*  theta_0 =
     * [0.,0.4*pi/180.0,0.4*pi/180.0,0.1*pi/180.0,0.1*pi/180.0,0.1*pi/180.0,0.,0.1,-40*pi/180]';
     */
    /* 'init_func:70' theta_0 = [0,0,0,0,0,0,0,0,0,0]'; */
    /* setting IC from past sim */
    /*  y0 = [0.017552353814854, -0.002156992032555, -0.002273627285241, ... */
    /*      -0.004091940730352,  -0.002796089196615,   0.019674817779806,... */
    /*      -0.017606183923045,                   0,                   0, ... */
    /*       0.207860712172010,  -0.003878840466313,  -0.004340266988222, ... */
    /*      -0.001098037684871,  -0.001085183886166,  -0.001924742862772, ... */
    /*       2.937417436471931,                   0,                   0, ... */
    /*                       0,                   0, 28.274274172758336]'; */
    /* 'init_func:80' theta_des =
     * [0.,0.4*pi/180.0,0.4*pi/180.0,0.1*pi/180.0,0.1*pi/180.0,0.1*pi/180.0,0.,0.1,40*pi/180]';
     */
    /* 'init_func:81' d_theta_dt_0 = [0.,0.,0.,0.,0.,0.,0.,0.,0.]'; */
    /* 'init_func:83' unlock = [1.0,1.0,1.0,1.0,1.0,1.0,1,1,1]'; */
    /* 'init_func:85' a_flex = a_f_func; */
    /* 'init_func:86' b_flex = b_f_func(); */
    /* 'init_func:87' a_df = a_mf_func(); */
    /* 'init_func:88' b_df = b_mf_func(); */
    offterm[0][0] = 0.0;
    offterm[1][1] = 0.0;
    offterm[2][2] = 0.0;
    for (k = 0; k < 9; k++) {
      /* 'init_func:46' mass = eye(3) * m_n(k); */
      /* 'init_func:47' i_k = reshape(i_n(k,:), 3, 3); */
      /* 'init_func:48' offterm = xmat(c_n(:,k)); */
      tau_app_flex_idx_0 = c_n[k][0];
      tau_app_flex_idx_1 = c_n[k][1];
      tau_app_flex_idx_2 = c_n[k][2];
      /* 'xmat:2' mat = [        0, -vec(3),  vec(2); */
      /* 'xmat:3'             vec(3),         0, -vec(1); */
      /* 'xmat:4'            -vec(2),  vec(1),         0]; */
      offterm[1][0] = -tau_app_flex_idx_2;
      offterm[2][0] = tau_app_flex_idx_1;
      offterm[0][1] = tau_app_flex_idx_2;
      offterm[2][1] = -tau_app_flex_idx_0;
      offterm[0][2] = -tau_app_flex_idx_1;
      offterm[1][2] = tau_app_flex_idx_0;
      /* 'init_func:50' m_w_n_i = [mass, -offterm; offterm, i_k]; */
      for (b_i = 0; b_i < 9; b_i++) {
        i_n[b_i] = (real_T)b_i_n[b_i][k];
      }
      /* 'init_func:51' m_w_n(:,:,k) = m_w_n_i; */
      b_i = iv[k];
      for (i1 = 0; i1 < 3; i1++) {
        m_w_n[k][i1][0] = (real_T)((int32_T)(((int32_T)a[i1][0]) * b_i));
        tau_app_flex_idx_0 = offterm[i1][0];
        m_w_n[k][i1 + 3][0] = -tau_app_flex_idx_0;
        m_w_n[k][i1][3] = tau_app_flex_idx_0;
        m_w_n[k][i1 + 3][3] = i_n[3 * i1];
        m_w_n[k][i1][1] = (real_T)((int32_T)(((int32_T)a[i1][1]) * b_i));
        tau_app_flex_idx_0 = offterm[i1][1];
        m_w_n[k][i1 + 3][1] = -tau_app_flex_idx_0;
        m_w_n[k][i1][4] = tau_app_flex_idx_0;
        m_w_n[k][i1 + 3][4] = i_n[(3 * i1) + 1];
        m_w_n[k][i1][2] = (real_T)((int32_T)(((int32_T)a[i1][2]) * b_i));
        tau_app_flex_idx_0 = offterm[i1][2];
        m_w_n[k][i1 + 3][2] = -tau_app_flex_idx_0;
        m_w_n[k][i1][5] = tau_app_flex_idx_0;
        m_w_n[k][i1 + 3][5] = i_n[(3 * i1) + 2];
        sys_workspace_z_n[k][i1] = (real_T)z_n[k][i1];
      }
      for (b_i = 0; b_i < 6; b_i++) {
        sys_workspace_p_n[k][b_i] = (real_T)p_n[k][b_i];
      }
    }
    for (k = 0; k < 9; k++) {
      i_n[k] = (real_T)m_n[k];
      sys_workspace_c_n[k][0] = c_n[k][0];
      sys_workspace_c_n[k][1] = c_n[k][1];
      sys_workspace_c_n[k][2] = c_n[k][2];
      sys_workspace_k_d[k] = k_d[k];
      sys_workspace_b_d[k] = 0.0;
    }
  }
  /* 'bit_one_step:14' if ~flexure_flag */
  if (!flexure_flag) {
    /* 'bit_one_step:15' k_d(8) = 0; */
    sys_workspace_k_d[7] = 0.0;
    /* 'bit_one_step:16' k_d(9) = 0; */
    sys_workspace_k_d[8] = 0.0;
  }
  /*     %% Setup Simulation */
  /*  initial conditions, state is dtheta; theta */
  /* 'bit_one_step:21' y_true = x0; */
  (void)memcpy(&y_true[0], &x0[0], 21U * (sizeof(real_T)));
  /* 'bit_one_step:22' y_flex = x_flex0; */
  (void)memcpy(&y_flex[0], &x_flex0[0], 104U * (sizeof(real_T)));
  /*  Sim Parameters */
  /*  y_all1 = zeros(18, tf/(dt)); */
  /* 'bit_one_step:26' step = 0; */
  /* 'bit_one_step:28' sys = @(y_true, tau_applied, dw_piv)
   * bit_propagator(y_true, c_n, z_n, m_n, r_n1_n, m_w_n, p_n, ...  */
  /* 'bit_one_step:29'     k_d, b_d, g0, unlock, hs_rw_max, tau_applied, w_piv,
   * piv_flag, dw_piv, tau_max_piv, thet_pit_nom); */
  /* 'bit_one_step:31' tau_app_flex = tau_applied(7:9); */
  /* 'bit_one_step:33' tau_applied(7) = tau_applied(7) + tau_flex(1); */
  tau_applied[6] += tau_flex[0];
  /* 'bit_one_step:34' tau_applied(8) = tau_applied(8) + tau_flex(2) +
   * tau_flex(3); */
  tau_applied[7] = (tau_flex[1] + tau_applied[7]) + tau_flex[2];
  /* 'bit_one_step:35' tau_applied(9) = tau_applied(9) + tau_flex(4) +
   * tau_flex(5); */
  tau_applied[8] = (tau_flex[3] + tau_applied[8]) + tau_flex[4];
  /* 'bit_one_step:37' sys_flex = @(y_flex, tau_app_flex, tau_flex)
   * flex_propogate(a_df, b_df, tau_app_flex, tau_flex, y_flex); */
  /*  sim */
  /* 'bit_one_step:41' for step = 1:num_steps */
  b_i = (int32_T)num_steps;
  for (step = 0; step < b_i; step++) {
    real_T b_df[104];
    real_T kf1[104];
    real_T kf2[104];
    real_T kf3[104];
    real_T k1[24];
    real_T k2[24];
    real_T k3[24];
    real_T varargout_1[24];
    real_T b_y_true[21];
    real_T b_tau_tmp;
    real_T c_tau_tmp;
    real_T d_tau_tmp;
    real_T tau_tmp;
    boolean_T th_over[9];
    boolean_T th_under[9];
    /*         %% Propagate the system */
    /* RK4 solver */
    /* 'bit_one_step:44' dw_piv = (w_piv - y_true(6))/dt; */
    tau_app_flex_idx_0 = (w_piv - y_true[5]) / dt;
    /* 'bit_one_step:46' [k1] = sys(y_true, tau_applied, dw_piv) * dt; */
    bit_one_step_anonFcn1(sys_workspace_c_n, sys_workspace_z_n, i_n,
                          *((real_T(*)[9][3])(&b_r_n1_n[0][0])), m_w_n,
                          sys_workspace_p_n, sys_workspace_k_d,
                          sys_workspace_b_d, g0, unlock, hs_rw_max, w_piv,
                          piv_flag, tau_max_piv, thet_pit_nom, y_true,
                          tau_applied, tau_app_flex_idx_0, k1);
    for (i1 = 0; i1 < 24; i1++) {
      k1[i1] *= dt;
    }
    /* 'bit_one_step:47' [k2] = sys(y_true + (k1(1:21)/2), tau_applied, dw_piv)
     * * dt; */
    for (i1 = 0; i1 < 21; i1++) {
      b_y_true[i1] = y_true[i1] + (k1[i1] / 2.0);
    }
    bit_one_step_anonFcn1(sys_workspace_c_n, sys_workspace_z_n, i_n,
                          *((real_T(*)[9][3])(&b_r_n1_n[0][0])), m_w_n,
                          sys_workspace_p_n, sys_workspace_k_d,
                          sys_workspace_b_d, g0, unlock, hs_rw_max, w_piv,
                          piv_flag, tau_max_piv, thet_pit_nom, b_y_true,
                          tau_applied, tau_app_flex_idx_0, k2);
    for (i1 = 0; i1 < 24; i1++) {
      k2[i1] *= dt;
    }
    /* 'bit_one_step:48' [k3] = sys(y_true + (k2(1:21)/2), tau_applied, dw_piv)
     * * dt; */
    for (i1 = 0; i1 < 21; i1++) {
      b_y_true[i1] = y_true[i1] + (k2[i1] / 2.0);
    }
    bit_one_step_anonFcn1(sys_workspace_c_n, sys_workspace_z_n, i_n,
                          *((real_T(*)[9][3])(&b_r_n1_n[0][0])), m_w_n,
                          sys_workspace_p_n, sys_workspace_k_d,
                          sys_workspace_b_d, g0, unlock, hs_rw_max, w_piv,
                          piv_flag, tau_max_piv, thet_pit_nom, b_y_true,
                          tau_applied, tau_app_flex_idx_0, k3);
    for (i1 = 0; i1 < 24; i1++) {
      k3[i1] *= dt;
    }
    /* 'bit_one_step:49' [k4] = sys(y_true + k3(1:21), tau_applied, dw_piv) *
     * dt; */
    for (i1 = 0; i1 < 21; i1++) {
      b_y_true[i1] = y_true[i1] + k3[i1];
    }
    bit_one_step_anonFcn1(sys_workspace_c_n, sys_workspace_z_n, i_n,
                          *((real_T(*)[9][3])(&b_r_n1_n[0][0])), m_w_n,
                          sys_workspace_p_n, sys_workspace_k_d,
                          sys_workspace_b_d, g0, unlock, hs_rw_max, w_piv,
                          piv_flag, tau_max_piv, thet_pit_nom, b_y_true,
                          tau_applied, tau_app_flex_idx_0, varargout_1);
    /* 'bit_one_step:51' temp = ((k1+(2*k2)+(2*k3)+k4)/6); */
    for (i1 = 0; i1 < 24; i1++) {
      k1[i1] = (((k1[i1] + (2.0 * k2[i1])) + (2.0 * k3[i1])) +
                (varargout_1[i1] * dt)) /
               6.0;
    }
    /* 'bit_one_step:52' tdd = temp(1:21); */
    /* 'bit_one_step:53' tau_app_flex = temp(22:24)/dt */
    /* 'bit_one_step:54' y_true = y_true + tdd; */
    for (i1 = 0; i1 < 21; i1++) {
      y_true[i1] += k1[i1];
    }
    /* 'bit_one_step:56' th_over = y_true(10:18) > pi; */
    /* 'bit_one_step:57' th_under = y_true(10:18) < -pi; */
    for (k = 0; k < 9; k++) {
      tau_app_flex_idx_0 = y_true[k + 9];
      th_over[k] = (tau_app_flex_idx_0 > 3.1415926535897931);
      th_under[k] = (tau_app_flex_idx_0 < -3.1415926535897931);
    }
    /* 'bit_one_step:58' y_true(10:14) = y_true(10:14) -(2*pi*th_over(1:5)) +
     * (2*pi*th_under(1:5)); */
    for (i1 = 0; i1 < 5; i1++) {
      y_true[i1 + 9] =
          (y_true[i1 + 9] -
           (6.2831853071795862 * ((real_T)(th_over[i1] ? 1.0 : 0.0)))) +
          (6.2831853071795862 * ((real_T)(th_under[i1] ? 1.0 : 0.0)));
    }
    /* 'bit_one_step:59' y_true(16:18) = y_true(16:18) -(2*pi*th_over(7:9)) +
     * (2*pi*th_under(7:9)); */
    y_true[15] = (y_true[15] -
                  (6.2831853071795862 * ((real_T)(th_over[6] ? 1.0 : 0.0)))) +
                 (6.2831853071795862 * ((real_T)(th_under[6] ? 1.0 : 0.0)));
    y_true[16] = (y_true[16] -
                  (6.2831853071795862 * ((real_T)(th_over[7] ? 1.0 : 0.0)))) +
                 (6.2831853071795862 * ((real_T)(th_under[7] ? 1.0 : 0.0)));
    y_true[17] = (y_true[17] -
                  (6.2831853071795862 * ((real_T)(th_over[8] ? 1.0 : 0.0)))) +
                 (6.2831853071795862 * ((real_T)(th_under[8] ? 1.0 : 0.0)));
    /*          fprintf('current state:  %0.15f \n  %0.15f \n %0.15f \n %0.15f
     * \n %0.15f \n %0.15f \n %0.15f \n %0.15f \n %0.15f \n %0.15f \n %0.15f \n
     * %0.15f \n %0.15f \n %0.15f \n %0.15f \n %0.15f \n %0.15f \n %0.15f \n
     * %0.15f \n %0.15f \n %0.15f \n \n', ... */
    /*              y_true(1), y_true(2), y_true(3), y_true(4), y_true(5),
     * y_true(6),... */
    /*               y_true(7), y_true(8), y_true(9), y_true(10), y_true(11),
     * y_true(12),... */
    /*                y_true(13), y_true(14), y_true(15), y_true(16),
     * y_true(17), y_true(18),... */
    /*                 y_true(19), y_true(20), y_true(21));       */
    /*         %% Propogate flexible system */
    /* 'bit_one_step:67' kf1 = sys_flex(y_flex, tau_app_flex, tau_flex) * dt; */
    /* 'bit_one_step:37' @(y_flex, tau_app_flex, tau_flex) flex_propogate(a_df,
     * b_df, tau_app_flex, tau_flex, y_flex) */
    /* UNTITLED Summary of this function goes here */
    /*    Detailed explanation goes here */
    /* 'flex_propogate:4' tau = tau_flex; */
    /* 'flex_propogate:6' tau_yaw = tau_applied(1) - tau(1); */
    /* 'flex_propogate:7' tau_roll = tau_applied(2) - (tau(2) + tau(3)); */
    tau_app_flex_idx_0 = (k1[22] / dt) - (tau_flex[1] + tau_flex[2]);
    /* 'flex_propogate:8' tau_pitch = tau_applied(3) - (tau(4) + tau(5)); */
    tau_app_flex_idx_1 = (k1[23] / dt) - (tau_flex[3] + tau_flex[4]);
    /*   */
    /* 'flex_propogate:11' tau(1) = -tau(1) + tau_yaw; */
    tau_tmp = (-tau_flex[0]) + ((k1[21] / dt) - tau_flex[0]);
    tau[0] = tau_tmp;
    /* 'flex_propogate:12' tau(2) = tau(2) + (tau_roll/2); */
    b_tau_tmp = tau_flex[1] + (tau_app_flex_idx_0 / 2.0);
    tau[1] = b_tau_tmp;
    /* 'flex_propogate:13' tau(3) = tau(3) + (tau_roll/2); */
    c_tau_tmp = tau_flex[2] + (tau_app_flex_idx_0 / 2.0);
    tau[2] = c_tau_tmp;
    /* 'flex_propogate:14' tau(4) = tau(4) + (tau_pitch/2); */
    d_tau_tmp = tau_flex[3] + (tau_app_flex_idx_1 / 2.0);
    tau[3] = d_tau_tmp;
    /* 'flex_propogate:15' tau(5) = tau(5) + (tau_pitch/2); */
    tau_app_flex_idx_1 = tau_flex[4] + (tau_app_flex_idx_1 / 2.0);
    tau[4] = tau_app_flex_idx_1;
    /* 'flex_propogate:18' eta_dot = (a_flex * x0_flex) + (b_flex * tau); */
    /* 'bit_one_step:68' kf2 = sys_flex(y_flex + (kf1/2), tau_app_flex,
     * tau_flex) * dt; */
    /* 'bit_one_step:37' @(y_flex, tau_app_flex, tau_flex) flex_propogate(a_df,
     * b_df, tau_app_flex, tau_flex, y_flex) */
    /* UNTITLED Summary of this function goes here */
    /*    Detailed explanation goes here */
    /* 'flex_propogate:4' tau = tau_flex; */
    /* 'flex_propogate:6' tau_yaw = tau_applied(1) - tau(1); */
    /* 'flex_propogate:7' tau_roll = tau_applied(2) - (tau(2) + tau(3)); */
    /* 'flex_propogate:8' tau_pitch = tau_applied(3) - (tau(4) + tau(5)); */
    /*   */
    /* 'flex_propogate:11' tau(1) = -tau(1) + tau_yaw; */
    b_tau[0] = tau_tmp;
    /* 'flex_propogate:12' tau(2) = tau(2) + (tau_roll/2); */
    b_tau[1] = b_tau_tmp;
    /* 'flex_propogate:13' tau(3) = tau(3) + (tau_roll/2); */
    b_tau[2] = c_tau_tmp;
    /* 'flex_propogate:14' tau(4) = tau(4) + (tau_pitch/2); */
    b_tau[3] = d_tau_tmp;
    /* 'flex_propogate:15' tau(5) = tau(5) + (tau_pitch/2); */
    b_tau[4] = tau_app_flex_idx_1;
    /* 'flex_propogate:18' eta_dot = (a_flex * x0_flex) + (b_flex * tau); */
    for (i1 = 0; i1 < 104; i1++) {
      tau_app_flex_idx_0 = 0.0;
      for (k = 0; k < 104; k++) {
        tau_app_flex_idx_0 += a_df[k][i1] * y_flex[k];
      }
      tau_app_flex_idx_2 = 0.0;
      for (k = 0; k < 5; k++) {
        tau_app_flex_idx_2 += b_b_df[k][i1] * tau[k];
      }
      tau_app_flex_idx_0 = (tau_app_flex_idx_0 + tau_app_flex_idx_2) * dt;
      kf1[i1] = tau_app_flex_idx_0;
      b_df[i1] = y_flex[i1] + (tau_app_flex_idx_0 / 2.0);
    }
    for (i1 = 0; i1 < 104; i1++) {
      tau_app_flex_idx_0 = 0.0;
      for (k = 0; k < 104; k++) {
        tau_app_flex_idx_0 += a_df[k][i1] * b_df[k];
      }
      kf2[i1] = tau_app_flex_idx_0;
    }
    /* 'bit_one_step:69' kf3 = sys_flex(y_flex + (kf2/2), tau_app_flex,
     * tau_flex) * dt; */
    /* 'bit_one_step:37' @(y_flex, tau_app_flex, tau_flex) flex_propogate(a_df,
     * b_df, tau_app_flex, tau_flex, y_flex) */
    /* UNTITLED Summary of this function goes here */
    /*    Detailed explanation goes here */
    /* 'flex_propogate:4' tau = tau_flex; */
    /* 'flex_propogate:6' tau_yaw = tau_applied(1) - tau(1); */
    /* 'flex_propogate:7' tau_roll = tau_applied(2) - (tau(2) + tau(3)); */
    /* 'flex_propogate:8' tau_pitch = tau_applied(3) - (tau(4) + tau(5)); */
    /*   */
    /* 'flex_propogate:11' tau(1) = -tau(1) + tau_yaw; */
    c_tau[0] = tau_tmp;
    /* 'flex_propogate:12' tau(2) = tau(2) + (tau_roll/2); */
    c_tau[1] = b_tau_tmp;
    /* 'flex_propogate:13' tau(3) = tau(3) + (tau_roll/2); */
    c_tau[2] = c_tau_tmp;
    /* 'flex_propogate:14' tau(4) = tau(4) + (tau_pitch/2); */
    c_tau[3] = d_tau_tmp;
    /* 'flex_propogate:15' tau(5) = tau(5) + (tau_pitch/2); */
    c_tau[4] = tau_app_flex_idx_1;
    /* 'flex_propogate:18' eta_dot = (a_flex * x0_flex) + (b_flex * tau); */
    for (i1 = 0; i1 < 104; i1++) {
      tau_app_flex_idx_0 = 0.0;
      for (k = 0; k < 5; k++) {
        tau_app_flex_idx_0 += b_b_df[k][i1] * b_tau[k];
      }
      tau_app_flex_idx_0 = (kf2[i1] + tau_app_flex_idx_0) * dt;
      kf2[i1] = tau_app_flex_idx_0;
      b_df[i1] = y_flex[i1] + (tau_app_flex_idx_0 / 2.0);
    }
    for (i1 = 0; i1 < 104; i1++) {
      tau_app_flex_idx_0 = 0.0;
      for (k = 0; k < 104; k++) {
        tau_app_flex_idx_0 += a_df[k][i1] * b_df[k];
      }
      kf3[i1] = tau_app_flex_idx_0;
    }
    /* 'bit_one_step:70' kf4 = sys_flex(y_flex + kf3, tau_app_flex, tau_flex) *
     * dt; */
    /* 'bit_one_step:37' @(y_flex, tau_app_flex, tau_flex) flex_propogate(a_df,
     * b_df, tau_app_flex, tau_flex, y_flex) */
    /* UNTITLED Summary of this function goes here */
    /*    Detailed explanation goes here */
    /* 'flex_propogate:4' tau = tau_flex; */
    /* 'flex_propogate:6' tau_yaw = tau_applied(1) - tau(1); */
    /* 'flex_propogate:7' tau_roll = tau_applied(2) - (tau(2) + tau(3)); */
    /* 'flex_propogate:8' tau_pitch = tau_applied(3) - (tau(4) + tau(5)); */
    /*   */
    /* 'flex_propogate:11' tau(1) = -tau(1) + tau_yaw; */
    d_tau[0] = tau_tmp;
    /* 'flex_propogate:12' tau(2) = tau(2) + (tau_roll/2); */
    d_tau[1] = b_tau_tmp;
    /* 'flex_propogate:13' tau(3) = tau(3) + (tau_roll/2); */
    d_tau[2] = c_tau_tmp;
    /* 'flex_propogate:14' tau(4) = tau(4) + (tau_pitch/2); */
    d_tau[3] = d_tau_tmp;
    /* 'flex_propogate:15' tau(5) = tau(5) + (tau_pitch/2); */
    d_tau[4] = tau_app_flex_idx_1;
    /* 'flex_propogate:18' eta_dot = (a_flex * x0_flex) + (b_flex * tau); */
    /* 'bit_one_step:72' eta_dd = ((kf1+(2*kf2)+(2*kf3)+kf4)/6); */
    /* 'bit_one_step:73' y_flex = y_flex + eta_dd; */
    for (i1 = 0; i1 < 104; i1++) {
      tau_app_flex_idx_0 = 0.0;
      for (k = 0; k < 5; k++) {
        tau_app_flex_idx_0 += b_b_df[k][i1] * c_tau[k];
      }
      tau_app_flex_idx_0 = (kf3[i1] + tau_app_flex_idx_0) * dt;
      kf3[i1] = tau_app_flex_idx_0;
      b_df[i1] = y_flex[i1] + tau_app_flex_idx_0;
    }
    for (i1 = 0; i1 < 104; i1++) {
      tau_app_flex_idx_0 = 0.0;
      for (k = 0; k < 104; k++) {
        tau_app_flex_idx_0 += a_df[k][i1] * b_df[k];
      }
      tau_app_flex_idx_2 = 0.0;
      for (k = 0; k < 5; k++) {
        tau_app_flex_idx_2 += b_b_df[k][i1] * d_tau[k];
      }
      y_flex[i1] += (((kf1[i1] + (2.0 * kf2[i1])) + (2.0 * kf3[i1])) +
                     ((tau_app_flex_idx_0 + tau_app_flex_idx_2) * dt)) /
                    6.0;
    }
  }
}

/*
 * function [omega] = compute_angular_velocity_C(x, z_n)
 *
 * UNTITLED2 Summary of this function goes here
 *    Detailed explanation goes here
 *
 * Arguments    : const real_T x[18]
 *                real_T z_n[9][3]
 *                real_T omega[3]
 * Return Type  : void
 */
void compute_angular_velocity_C(const real_T x[18], real_T z_n[9][3],
                                real_T omega[3])
{
  real_T s9[9][3];
  real_T d;
  int32_T b_i;
  int32_T c_i;
  int32_T i1;
  /* 'compute_angular_velocity_C:4' theta = x(10:18); */
  /* 'compute_angular_velocity_C:5' dtheta = x(1:9); */
  /* 'compute_angular_velocity_C:7' s9 = zeros(3,9); */
  for (b_i = 0; b_i < 9; b_i++) {
    s9[b_i][0] = 0.0;
    s9[b_i][1] = 0.0;
    s9[b_i][2] = 0.0;
  }
  /* 'compute_angular_velocity_C:8' for i = 1:9 */
  for (c_i = 0; c_i < 9; c_i++) {
    real_T g_Cn[9][3];
    real_T f_Cn[3][3];
    /* 'compute_angular_velocity_C:9' Cn = axis2rot(z_n(:,i), theta(i)); */
    axis2rot(*((real_T(*)[3])(&z_n[c_i][0])), x[c_i + 9], f_Cn);
    /* 'compute_angular_velocity_C:10' s9(:,i) = z_n(:,i); */
    s9[c_i][0] = z_n[c_i][0];
    s9[c_i][1] = z_n[c_i][1];
    s9[c_i][2] = z_n[c_i][2];
    /* 'compute_angular_velocity_C:11' s9 = Cn*s9; */
    for (b_i = 0; b_i < 3; b_i++) {
      real_T d1;
      real_T d2;
      d = f_Cn[0][b_i];
      d1 = f_Cn[1][b_i];
      d2 = f_Cn[2][b_i];
      for (i1 = 0; i1 < 9; i1++) {
        g_Cn[i1][b_i] = ((d * s9[i1][0]) + (d1 * s9[i1][1])) + (d2 * s9[i1][2]);
      }
    }
    for (b_i = 0; b_i < 9; b_i++) {
      s9[b_i][0] = g_Cn[b_i][0];
      s9[b_i][1] = g_Cn[b_i][1];
      s9[b_i][2] = g_Cn[b_i][2];
    }
  }
  /* 'compute_angular_velocity_C:14' omega = s9 * dtheta; */
  for (b_i = 0; b_i < 3; b_i++) {
    d = 0.0;
    for (i1 = 0; i1 < 9; i1++) {
      d += s9[i1][b_i] * x[i1];
    }
    omega[b_i] = d;
  }
}

/*
 * function [omega] = compute_angular_velocity_roll_C(x, z_n)
 *
 * UNTITLED2 Summary of this function goes here
 *    Detailed explanation goes here
 *
 * Arguments    : const real_T x[18]
 *                real_T z_n[9][3]
 *                real_T omega[3]
 * Return Type  : void
 */
void compute_angular_velocity_roll_C(const real_T x[18], real_T z_n[9][3],
                                     real_T omega[3])
{
  real_T s8[8][3];
  real_T d;
  int32_T b_i;
  int32_T c_i;
  int32_T i1;
  /* 'compute_angular_velocity_roll_C:4' theta = x(10:18); */
  /* 'compute_angular_velocity_roll_C:5' dtheta = x(1:8); */
  /* 'compute_angular_velocity_roll_C:7' s8 = zeros(3,8); */
  for (b_i = 0; b_i < 8; b_i++) {
    s8[b_i][0] = 0.0;
    s8[b_i][1] = 0.0;
    s8[b_i][2] = 0.0;
  }
  /* 'compute_angular_velocity_roll_C:8' for i = 1:8 */
  for (c_i = 0; c_i < 8; c_i++) {
    real_T g_Cn[8][3];
    real_T f_Cn[3][3];
    /* 'compute_angular_velocity_roll_C:9' Cn = axis2rot(z_n(:,i), theta(i)); */
    axis2rot(*((real_T(*)[3])(&z_n[c_i][0])), x[c_i + 9], f_Cn);
    /* 'compute_angular_velocity_roll_C:10' s8(:,i) = z_n(:,i); */
    s8[c_i][0] = z_n[c_i][0];
    s8[c_i][1] = z_n[c_i][1];
    s8[c_i][2] = z_n[c_i][2];
    /* 'compute_angular_velocity_roll_C:11' s8 = Cn*s8; */
    for (b_i = 0; b_i < 3; b_i++) {
      real_T d1;
      real_T d2;
      d = f_Cn[0][b_i];
      d1 = f_Cn[1][b_i];
      d2 = f_Cn[2][b_i];
      for (i1 = 0; i1 < 8; i1++) {
        g_Cn[i1][b_i] = ((d * s8[i1][0]) + (d1 * s8[i1][1])) + (d2 * s8[i1][2]);
      }
    }
    for (b_i = 0; b_i < 8; b_i++) {
      s8[b_i][0] = g_Cn[b_i][0];
      s8[b_i][1] = g_Cn[b_i][1];
      s8[b_i][2] = g_Cn[b_i][2];
    }
  }
  /* 'compute_angular_velocity_roll_C:14' omega = s8 * dtheta; */
  for (b_i = 0; b_i < 3; b_i++) {
    d = 0.0;
    for (i1 = 0; i1 < 8; i1++) {
      d += s8[i1][b_i] * x[i1];
    }
    omega[b_i] = d;
  }
}

/*
 * function [omega] = compute_angular_velocity_yaw_C(x, z_n)
 *
 * UNTITLED2 Summary of this function goes here
 *    Detailed explanation goes here
 *
 * Arguments    : const real_T x[18]
 *                real_T z_n[9][3]
 *                real_T omega[3]
 * Return Type  : void
 */
void compute_angular_velocity_yaw_C(const real_T x[18], real_T z_n[9][3],
                                    real_T omega[3])
{
  real_T s7[7][3];
  real_T d;
  int32_T b_i;
  int32_T c_i;
  int32_T i1;
  /* 'compute_angular_velocity_yaw_C:4' theta = x(10:18); */
  /* 'compute_angular_velocity_yaw_C:5' dtheta = x(1:7); */
  /* 'compute_angular_velocity_yaw_C:7' s7 = zeros(3,7); */
  for (b_i = 0; b_i < 7; b_i++) {
    s7[b_i][0] = 0.0;
    s7[b_i][1] = 0.0;
    s7[b_i][2] = 0.0;
  }
  /* 'compute_angular_velocity_yaw_C:8' for i = 1:7 */
  for (c_i = 0; c_i < 7; c_i++) {
    real_T g_Cn[7][3];
    real_T f_Cn[3][3];
    /* 'compute_angular_velocity_yaw_C:9' Cn = axis2rot(z_n(:,i), theta(i)); */
    axis2rot(*((real_T(*)[3])(&z_n[c_i][0])), x[c_i + 9], f_Cn);
    /* 'compute_angular_velocity_yaw_C:10' s7(:,i) = z_n(:,i); */
    s7[c_i][0] = z_n[c_i][0];
    s7[c_i][1] = z_n[c_i][1];
    s7[c_i][2] = z_n[c_i][2];
    /* 'compute_angular_velocity_yaw_C:11' s7 = Cn*s7; */
    for (b_i = 0; b_i < 3; b_i++) {
      real_T d1;
      real_T d2;
      d = f_Cn[0][b_i];
      d1 = f_Cn[1][b_i];
      d2 = f_Cn[2][b_i];
      for (i1 = 0; i1 < 7; i1++) {
        g_Cn[i1][b_i] = ((d * s7[i1][0]) + (d1 * s7[i1][1])) + (d2 * s7[i1][2]);
      }
    }
    for (b_i = 0; b_i < 7; b_i++) {
      s7[b_i][0] = g_Cn[b_i][0];
      s7[b_i][1] = g_Cn[b_i][1];
      s7[b_i][2] = g_Cn[b_i][2];
    }
  }
  /* 'compute_angular_velocity_yaw_C:14' omega = s7 * dtheta; */
  for (b_i = 0; b_i < 3; b_i++) {
    d = 0.0;
    for (i1 = 0; i1 < 7; i1++) {
      d += s7[i1][b_i] * x[i1];
    }
    omega[b_i] = d;
  }
}

/*
 * function [C] = compute_rotation_mat_C(z_n, theta)
 *
 * UNTITLED3 Summary of this function goes here
 *    Detailed explanation goes here
 *
 * Arguments    : real_T z_n[9][3]
 *                const real_T theta[9]
 *                real_T C[3][3]
 * Return Type  : void
 */
void compute_rotation_mat_C(real_T z_n[9][3], const real_T theta[9],
                            real_T C[3][3])
{
  real_T b_a[3][3];
  int32_T b_i;
  int32_T c_i;
  int32_T i1;
  /* 'compute_rotation_mat_C:4' C = (eye(3)); */
  for (b_i = 0; b_i < 3; b_i++) {
    C[b_i][0] = 0.0;
    C[b_i][1] = 0.0;
    C[b_i][2] = 0.0;
  }
  C[0][0] = 1.0;
  C[1][1] = 1.0;
  C[2][2] = 1.0;
  /* 'compute_rotation_mat_C:5' for i = 1:9 */
  for (c_i = 0; c_i < 9; c_i++) {
    real_T a[3][3];
    /* 'compute_rotation_mat_C:6' C = axis2rot(z_n(:,i), theta(i)) * C; */
    axis2rot(*((real_T(*)[3])(&z_n[c_i][0])), theta[c_i], b_a);
    for (b_i = 0; b_i < 3; b_i++) {
      real_T d;
      real_T d1;
      real_T d2;
      d = b_a[0][b_i];
      d1 = b_a[1][b_i];
      d2 = b_a[2][b_i];
      for (i1 = 0; i1 < 3; i1++) {
        a[i1][b_i] = ((d * C[i1][0]) + (d1 * C[i1][1])) + (d2 * C[i1][2]);
      }
    }
    for (b_i = 0; b_i < 3; b_i++) {
      C[b_i][0] = a[b_i][0];
      C[b_i][1] = a[b_i][1];
      C[b_i][2] = a[b_i][2];
    }
  }
  /* 'compute_rotation_mat_C:8' C = C'; */
  for (b_i = 0; b_i < 3; b_i++) {
    b_a[b_i][0] = C[0][b_i];
    b_a[b_i][1] = C[1][b_i];
    b_a[b_i][2] = C[2][b_i];
  }
  for (b_i = 0; b_i < 3; b_i++) {
    C[b_i][0] = b_a[b_i][0];
    C[b_i][1] = b_a[b_i][1];
    C[b_i][2] = b_a[b_i][2];
  }
}

/*
 * function [C] = compute_rotation_mat_roll_C(z_n, theta)
 *
 * UNTITLED3 Summary of this function goes here
 *    Detailed explanation goes here
 *
 * Arguments    : real_T z_n[9][3]
 *                const real_T theta[9]
 *                real_T C[3][3]
 * Return Type  : void
 */
void compute_rotation_mat_roll_C(real_T z_n[9][3], const real_T theta[9],
                                 real_T C[3][3])
{
  real_T b_a[3][3];
  int32_T b_i;
  int32_T c_i;
  int32_T i1;
  /* 'compute_rotation_mat_roll_C:4' C = (eye(3)); */
  for (b_i = 0; b_i < 3; b_i++) {
    C[b_i][0] = 0.0;
    C[b_i][1] = 0.0;
    C[b_i][2] = 0.0;
  }
  C[0][0] = 1.0;
  C[1][1] = 1.0;
  C[2][2] = 1.0;
  /* 'compute_rotation_mat_roll_C:5' for i = 1:8 */
  for (c_i = 0; c_i < 8; c_i++) {
    real_T a[3][3];
    /* 'compute_rotation_mat_roll_C:6' C = axis2rot(z_n(:,i), theta(i)) * C; */
    axis2rot(*((real_T(*)[3])(&z_n[c_i][0])), theta[c_i], b_a);
    for (b_i = 0; b_i < 3; b_i++) {
      real_T d;
      real_T d1;
      real_T d2;
      d = b_a[0][b_i];
      d1 = b_a[1][b_i];
      d2 = b_a[2][b_i];
      for (i1 = 0; i1 < 3; i1++) {
        a[i1][b_i] = ((d * C[i1][0]) + (d1 * C[i1][1])) + (d2 * C[i1][2]);
      }
    }
    for (b_i = 0; b_i < 3; b_i++) {
      C[b_i][0] = a[b_i][0];
      C[b_i][1] = a[b_i][1];
      C[b_i][2] = a[b_i][2];
    }
  }
  /* 'compute_rotation_mat_roll_C:8' C = C'; */
  for (b_i = 0; b_i < 3; b_i++) {
    b_a[b_i][0] = C[0][b_i];
    b_a[b_i][1] = C[1][b_i];
    b_a[b_i][2] = C[2][b_i];
  }
  for (b_i = 0; b_i < 3; b_i++) {
    C[b_i][0] = b_a[b_i][0];
    C[b_i][1] = b_a[b_i][1];
    C[b_i][2] = b_a[b_i][2];
  }
}

/*
 * function [C] = compute_rotation_mat_yaw_C(z_n, theta)
 *
 * UNTITLED3 Summary of this function goes here
 *    Detailed explanation goes here
 *
 * Arguments    : real_T z_n[9][3]
 *                const real_T theta[9]
 *                real_T C[3][3]
 * Return Type  : void
 */
void compute_rotation_mat_yaw_C(real_T z_n[9][3], const real_T theta[9],
                                real_T C[3][3])
{
  real_T b_a[3][3];
  int32_T b_i;
  int32_T c_i;
  int32_T i1;
  /* 'compute_rotation_mat_yaw_C:4' C = (eye(3)); */
  for (b_i = 0; b_i < 3; b_i++) {
    C[b_i][0] = 0.0;
    C[b_i][1] = 0.0;
    C[b_i][2] = 0.0;
  }
  C[0][0] = 1.0;
  C[1][1] = 1.0;
  C[2][2] = 1.0;
  /* 'compute_rotation_mat_yaw_C:5' for i = 1:7 */
  for (c_i = 0; c_i < 7; c_i++) {
    real_T a[3][3];
    /* 'compute_rotation_mat_yaw_C:6' C = axis2rot(z_n(:,i), theta(i)) * C; */
    axis2rot(*((real_T(*)[3])(&z_n[c_i][0])), theta[c_i], b_a);
    for (b_i = 0; b_i < 3; b_i++) {
      real_T d;
      real_T d1;
      real_T d2;
      d = b_a[0][b_i];
      d1 = b_a[1][b_i];
      d2 = b_a[2][b_i];
      for (i1 = 0; i1 < 3; i1++) {
        a[i1][b_i] = ((d * C[i1][0]) + (d1 * C[i1][1])) + (d2 * C[i1][2]);
      }
    }
    for (b_i = 0; b_i < 3; b_i++) {
      C[b_i][0] = a[b_i][0];
      C[b_i][1] = a[b_i][1];
      C[b_i][2] = a[b_i][2];
    }
  }
  /* 'compute_rotation_mat_yaw_C:8' C = C'; */
  for (b_i = 0; b_i < 3; b_i++) {
    b_a[b_i][0] = C[0][b_i];
    b_a[b_i][1] = C[1][b_i];
    b_a[b_i][2] = C[2][b_i];
  }
  for (b_i = 0; b_i < 3; b_i++) {
    C[b_i][0] = b_a[b_i][0];
    C[b_i][1] = b_a[b_i][1];
    C[b_i][2] = b_a[b_i][2];
  }
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void libbitonestep_initialize(void)
{
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void libbitonestep_terminate(void)
{
}

/*
 * function [v, phi] = rot2axisC(C)
 *
 * Arguments    : real_T C[3][3]
 *                real_T v[3]
 *                real_T *phi
 * Return Type  : void
 */
void rot2axis_C(real_T C[3][3], real_T v[3], real_T *phi)
{
  static rtRunTimeErrorInfo b_emlrtRTEI = {
      13,    /* lineNo */
      "sqrt" /* fName */
  };
  static rtRunTimeErrorInfo emlrtRTEI = {
      14,    /* lineNo */
      "acos" /* fName */
  };
  real_T b_v;
  /* 'rot2axis_C:2' phi = acos((C(1,1) + C(2,2) + C(3,3) - 1)/2); */
  *phi = (((C[0][0] + C[1][1]) + C[2][2]) - 1.0) / 2.0;
  if (((*phi) < -1.0) || ((*phi) > 1.0)) {
    d_rtErrorWithMessageID("acos", emlrtRTEI.fName, emlrtRTEI.lineNo);
  }
  *phi = acos(*phi);
  /* 'rot2axis_C:4' sinphi = sin(phi); */
  /* 'rot2axis_C:5' v = [C(3,2) - C(2,3); C(1,3) - C(3,1); C(2,1) -
   * C(1,2)]/(2*sinphi); */
  b_v = 2.0 * sin(*phi);
  v[0] = (C[1][2] - C[2][1]) / b_v;
  v[1] = (C[2][0] - C[0][2]) / b_v;
  v[2] = (C[0][1] - C[1][0]) / b_v;
  /* 'rot2axis_C:6' v = v/sqrt(v'*v); */
  b_v = ((v[0] * v[0]) + (v[1] * v[1])) + (v[2] * v[2]);
  if (b_v < 0.0) {
    d_rtErrorWithMessageID("sqrt", b_emlrtRTEI.fName, b_emlrtRTEI.lineNo);
  }
  b_v = sqrt(b_v);
  v[0] /= b_v;
  v[1] /= b_v;
  v[2] /= b_v;
}

/*
 * File trailer for libbitonestep.c
 *
 * [EOF]
 */
