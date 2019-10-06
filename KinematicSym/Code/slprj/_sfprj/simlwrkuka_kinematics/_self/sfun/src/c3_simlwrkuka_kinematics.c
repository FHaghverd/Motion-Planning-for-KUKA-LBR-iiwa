/* Include files */

#include <stddef.h>
#include "blas.h"
#include "simlwrkuka_kinematics_sfun.h"
#include "c3_simlwrkuka_kinematics.h"
#include "mwmathutil.h"
#define CHARTINSTANCE_CHARTNUMBER      (chartInstance->chartNumber)
#define CHARTINSTANCE_INSTANCENUMBER   (chartInstance->instanceNumber)
#include "simlwrkuka_kinematics_sfun_debug_macros.h"
#define _SF_MEX_LISTEN_FOR_CTRL_C(S)   sf_mex_listen_for_ctrl_c(sfGlobalDebugInstanceStruct,S);

/* Type Definitions */

/* Named Constants */
#define CALL_EVENT                     (-1)

/* Variable Declarations */

/* Variable Definitions */
static real_T _sfTime_;
static const char * c3_debug_family_names[19] = { "q1", "q2", "q3", "q4", "q5",
  "q6", "q7", "Te", "p", "R7", "eta", "e", "Q", "nargin", "nargout", "q", "x",
  "xp", "xe" };

static const char * c3_b_debug_family_names[49] = { "d3", "d5", "d0", "d7", "f",
  "a", "d", "t", "A0", "A1", "A2", "A3", "A4", "A5", "A6", "A7", "Ae", "T1",
  "T2", "T3", "T4", "T5", "T6", "z0", "z1", "z2", "z3", "z4", "z5", "z6", "p0",
  "p1", "p2", "p3", "p4", "p5", "p6", "pe", "nargin", "nargout", "q1", "q2",
  "q3", "q4", "q5", "q6", "q7", "je", "Te" };

/* Function Declarations */
static void initialize_c3_simlwrkuka_kinematics
  (SFc3_simlwrkuka_kinematicsInstanceStruct *chartInstance);
static void initialize_params_c3_simlwrkuka_kinematics
  (SFc3_simlwrkuka_kinematicsInstanceStruct *chartInstance);
static void enable_c3_simlwrkuka_kinematics
  (SFc3_simlwrkuka_kinematicsInstanceStruct *chartInstance);
static void disable_c3_simlwrkuka_kinematics
  (SFc3_simlwrkuka_kinematicsInstanceStruct *chartInstance);
static void c3_update_debugger_state_c3_simlwrkuka_kinematics
  (SFc3_simlwrkuka_kinematicsInstanceStruct *chartInstance);
static const mxArray *get_sim_state_c3_simlwrkuka_kinematics
  (SFc3_simlwrkuka_kinematicsInstanceStruct *chartInstance);
static void set_sim_state_c3_simlwrkuka_kinematics
  (SFc3_simlwrkuka_kinematicsInstanceStruct *chartInstance, const mxArray *c3_st);
static void finalize_c3_simlwrkuka_kinematics
  (SFc3_simlwrkuka_kinematicsInstanceStruct *chartInstance);
static void sf_gateway_c3_simlwrkuka_kinematics
  (SFc3_simlwrkuka_kinematicsInstanceStruct *chartInstance);
static void initSimStructsc3_simlwrkuka_kinematics
  (SFc3_simlwrkuka_kinematicsInstanceStruct *chartInstance);
static void c3_jacobin(SFc3_simlwrkuka_kinematicsInstanceStruct *chartInstance,
  real_T c3_q1, real_T c3_q2, real_T c3_q3, real_T c3_q4, real_T c3_q5, real_T
  c3_q6, real_T c3_q7, real_T c3_je[42], real_T c3_Te[16]);
static void init_script_number_translation(uint32_T c3_machineNumber, uint32_T
  c3_chartNumber, uint32_T c3_instanceNumber);
static const mxArray *c3_sf_marshallOut(void *chartInstanceVoid, void *c3_inData);
static void c3_emlrt_marshallIn(SFc3_simlwrkuka_kinematicsInstanceStruct
  *chartInstance, const mxArray *c3_xe, const char_T *c3_identifier, real_T
  c3_y[3]);
static void c3_b_emlrt_marshallIn(SFc3_simlwrkuka_kinematicsInstanceStruct
  *chartInstance, const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId,
  real_T c3_y[3]);
static void c3_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData);
static const mxArray *c3_b_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData);
static void c3_c_emlrt_marshallIn(SFc3_simlwrkuka_kinematicsInstanceStruct
  *chartInstance, const mxArray *c3_x, const char_T *c3_identifier, real_T c3_y
  [7]);
static void c3_d_emlrt_marshallIn(SFc3_simlwrkuka_kinematicsInstanceStruct
  *chartInstance, const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId,
  real_T c3_y[7]);
static void c3_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData);
static const mxArray *c3_c_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData);
static real_T c3_e_emlrt_marshallIn(SFc3_simlwrkuka_kinematicsInstanceStruct
  *chartInstance, const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId);
static void c3_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData);
static const mxArray *c3_d_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData);
static void c3_f_emlrt_marshallIn(SFc3_simlwrkuka_kinematicsInstanceStruct
  *chartInstance, const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId,
  real_T c3_y[4]);
static void c3_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData);
static const mxArray *c3_e_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData);
static void c3_g_emlrt_marshallIn(SFc3_simlwrkuka_kinematicsInstanceStruct
  *chartInstance, const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId,
  real_T c3_y[9]);
static void c3_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData);
static const mxArray *c3_f_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData);
static void c3_h_emlrt_marshallIn(SFc3_simlwrkuka_kinematicsInstanceStruct
  *chartInstance, const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId,
  real_T c3_y[16]);
static void c3_f_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData);
static const mxArray *c3_g_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData);
static void c3_i_emlrt_marshallIn(SFc3_simlwrkuka_kinematicsInstanceStruct
  *chartInstance, const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId,
  real_T c3_y[42]);
static void c3_g_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData);
static void c3_info_helper(const mxArray **c3_info);
static const mxArray *c3_emlrt_marshallOut(const char * c3_u);
static const mxArray *c3_b_emlrt_marshallOut(const uint32_T c3_u);
static void c3_eml_scalar_eg(SFc3_simlwrkuka_kinematicsInstanceStruct
  *chartInstance);
static void c3_eml_xgemm(SFc3_simlwrkuka_kinematicsInstanceStruct *chartInstance,
  real_T c3_A[16], real_T c3_B[16], real_T c3_C[16], real_T c3_b_C[16]);
static void c3_eml_error(SFc3_simlwrkuka_kinematicsInstanceStruct *chartInstance);
static const mxArray *c3_h_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData);
static int32_T c3_j_emlrt_marshallIn(SFc3_simlwrkuka_kinematicsInstanceStruct
  *chartInstance, const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId);
static void c3_h_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData);
static uint8_T c3_k_emlrt_marshallIn(SFc3_simlwrkuka_kinematicsInstanceStruct
  *chartInstance, const mxArray *c3_b_is_active_c3_simlwrkuka_kinematics, const
  char_T *c3_identifier);
static uint8_T c3_l_emlrt_marshallIn(SFc3_simlwrkuka_kinematicsInstanceStruct
  *chartInstance, const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId);
static void c3_b_eml_xgemm(SFc3_simlwrkuka_kinematicsInstanceStruct
  *chartInstance, real_T c3_A[16], real_T c3_B[16], real_T c3_C[16]);
static void init_dsm_address_info(SFc3_simlwrkuka_kinematicsInstanceStruct
  *chartInstance);

/* Function Definitions */
static void initialize_c3_simlwrkuka_kinematics
  (SFc3_simlwrkuka_kinematicsInstanceStruct *chartInstance)
{
  chartInstance->c3_sfEvent = CALL_EVENT;
  _sfTime_ = sf_get_time(chartInstance->S);
  chartInstance->c3_is_active_c3_simlwrkuka_kinematics = 0U;
}

static void initialize_params_c3_simlwrkuka_kinematics
  (SFc3_simlwrkuka_kinematicsInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void enable_c3_simlwrkuka_kinematics
  (SFc3_simlwrkuka_kinematicsInstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void disable_c3_simlwrkuka_kinematics
  (SFc3_simlwrkuka_kinematicsInstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void c3_update_debugger_state_c3_simlwrkuka_kinematics
  (SFc3_simlwrkuka_kinematicsInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static const mxArray *get_sim_state_c3_simlwrkuka_kinematics
  (SFc3_simlwrkuka_kinematicsInstanceStruct *chartInstance)
{
  const mxArray *c3_st;
  const mxArray *c3_y = NULL;
  int32_T c3_i0;
  real_T c3_u[7];
  const mxArray *c3_b_y = NULL;
  int32_T c3_i1;
  real_T c3_b_u[3];
  const mxArray *c3_c_y = NULL;
  int32_T c3_i2;
  real_T c3_c_u[3];
  const mxArray *c3_d_y = NULL;
  uint8_T c3_hoistedGlobal;
  uint8_T c3_d_u;
  const mxArray *c3_e_y = NULL;
  real_T (*c3_xp)[3];
  real_T (*c3_xe)[3];
  real_T (*c3_x)[7];
  c3_xe = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S, 3);
  c3_xp = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S, 2);
  c3_x = (real_T (*)[7])ssGetOutputPortSignal(chartInstance->S, 1);
  c3_st = NULL;
  c3_st = NULL;
  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_createcellmatrix(4, 1), false);
  for (c3_i0 = 0; c3_i0 < 7; c3_i0++) {
    c3_u[c3_i0] = (*c3_x)[c3_i0];
  }

  c3_b_y = NULL;
  sf_mex_assign(&c3_b_y, sf_mex_create("y", c3_u, 0, 0U, 1U, 0U, 1, 7), false);
  sf_mex_setcell(c3_y, 0, c3_b_y);
  for (c3_i1 = 0; c3_i1 < 3; c3_i1++) {
    c3_b_u[c3_i1] = (*c3_xe)[c3_i1];
  }

  c3_c_y = NULL;
  sf_mex_assign(&c3_c_y, sf_mex_create("y", c3_b_u, 0, 0U, 1U, 0U, 1, 3), false);
  sf_mex_setcell(c3_y, 1, c3_c_y);
  for (c3_i2 = 0; c3_i2 < 3; c3_i2++) {
    c3_c_u[c3_i2] = (*c3_xp)[c3_i2];
  }

  c3_d_y = NULL;
  sf_mex_assign(&c3_d_y, sf_mex_create("y", c3_c_u, 0, 0U, 1U, 0U, 1, 3), false);
  sf_mex_setcell(c3_y, 2, c3_d_y);
  c3_hoistedGlobal = chartInstance->c3_is_active_c3_simlwrkuka_kinematics;
  c3_d_u = c3_hoistedGlobal;
  c3_e_y = NULL;
  sf_mex_assign(&c3_e_y, sf_mex_create("y", &c3_d_u, 3, 0U, 0U, 0U, 0), false);
  sf_mex_setcell(c3_y, 3, c3_e_y);
  sf_mex_assign(&c3_st, c3_y, false);
  return c3_st;
}

static void set_sim_state_c3_simlwrkuka_kinematics
  (SFc3_simlwrkuka_kinematicsInstanceStruct *chartInstance, const mxArray *c3_st)
{
  const mxArray *c3_u;
  real_T c3_dv0[7];
  int32_T c3_i3;
  real_T c3_dv1[3];
  int32_T c3_i4;
  real_T c3_dv2[3];
  int32_T c3_i5;
  real_T (*c3_x)[7];
  real_T (*c3_xe)[3];
  real_T (*c3_xp)[3];
  c3_xe = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S, 3);
  c3_xp = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S, 2);
  c3_x = (real_T (*)[7])ssGetOutputPortSignal(chartInstance->S, 1);
  chartInstance->c3_doneDoubleBufferReInit = true;
  c3_u = sf_mex_dup(c3_st);
  c3_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c3_u, 0)), "x",
                        c3_dv0);
  for (c3_i3 = 0; c3_i3 < 7; c3_i3++) {
    (*c3_x)[c3_i3] = c3_dv0[c3_i3];
  }

  c3_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c3_u, 1)), "xe",
                      c3_dv1);
  for (c3_i4 = 0; c3_i4 < 3; c3_i4++) {
    (*c3_xe)[c3_i4] = c3_dv1[c3_i4];
  }

  c3_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c3_u, 2)), "xp",
                      c3_dv2);
  for (c3_i5 = 0; c3_i5 < 3; c3_i5++) {
    (*c3_xp)[c3_i5] = c3_dv2[c3_i5];
  }

  chartInstance->c3_is_active_c3_simlwrkuka_kinematics = c3_k_emlrt_marshallIn
    (chartInstance, sf_mex_dup(sf_mex_getcell(c3_u, 3)),
     "is_active_c3_simlwrkuka_kinematics");
  sf_mex_destroy(&c3_u);
  c3_update_debugger_state_c3_simlwrkuka_kinematics(chartInstance);
  sf_mex_destroy(&c3_st);
}

static void finalize_c3_simlwrkuka_kinematics
  (SFc3_simlwrkuka_kinematicsInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void sf_gateway_c3_simlwrkuka_kinematics
  (SFc3_simlwrkuka_kinematicsInstanceStruct *chartInstance)
{
  int32_T c3_i6;
  int32_T c3_i7;
  real_T c3_q[7];
  uint32_T c3_debug_family_var_map[19];
  real_T c3_q1;
  real_T c3_q2;
  real_T c3_q3;
  real_T c3_q4;
  real_T c3_q5;
  real_T c3_q6;
  real_T c3_q7;
  real_T c3_Te[16];
  real_T c3_p[3];
  real_T c3_R7[9];
  real_T c3_eta;
  real_T c3_e[3];
  real_T c3_Q[4];
  real_T c3_nargin = 1.0;
  real_T c3_nargout = 3.0;
  real_T c3_x[7];
  real_T c3_xp[3];
  real_T c3_xe[3];
  real_T c3_b_Te[16];
  real_T c3_unusedU0[42];
  int32_T c3_i8;
  int32_T c3_i9;
  int32_T c3_i10;
  int32_T c3_i11;
  int32_T c3_i12;
  int32_T c3_i13;
  real_T c3_b_x;
  real_T c3_c_x;
  real_T c3_d_x;
  real_T c3_e_x;
  real_T c3_f_x;
  real_T c3_g_x;
  real_T c3_h_x;
  real_T c3_i_x;
  real_T c3_j_x;
  real_T c3_k_x;
  real_T c3_l_x;
  real_T c3_m_x;
  real_T c3_n_x;
  real_T c3_o_x;
  real_T c3_b[3];
  int32_T c3_i14;
  int32_T c3_i15;
  int32_T c3_i16;
  int32_T c3_i17;
  int32_T c3_i18;
  int32_T c3_i19;
  int32_T c3_i20;
  int32_T c3_i21;
  int32_T c3_i22;
  int32_T c3_i23;
  int32_T c3_i24;
  int32_T c3_i25;
  real_T (*c3_p_x)[7];
  real_T (*c3_b_xp)[3];
  real_T (*c3_b_xe)[3];
  real_T (*c3_b_q)[7];
  c3_b_xe = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S, 3);
  c3_b_xp = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S, 2);
  c3_p_x = (real_T (*)[7])ssGetOutputPortSignal(chartInstance->S, 1);
  c3_b_q = (real_T (*)[7])ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = sf_get_time(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 2U, chartInstance->c3_sfEvent);
  for (c3_i6 = 0; c3_i6 < 7; c3_i6++) {
    _SFD_DATA_RANGE_CHECK((*c3_b_q)[c3_i6], 0U);
  }

  chartInstance->c3_sfEvent = CALL_EVENT;
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 2U, chartInstance->c3_sfEvent);
  for (c3_i7 = 0; c3_i7 < 7; c3_i7++) {
    c3_q[c3_i7] = (*c3_b_q)[c3_i7];
  }

  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 19U, 19U, c3_debug_family_names,
    c3_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_q1, 0U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_q2, 1U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_q3, 2U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_q4, 3U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_q5, 4U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_q6, 5U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_q7, 6U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_Te, 7U, c3_f_sf_marshallOut,
    c3_f_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_p, 8U, c3_sf_marshallOut,
    c3_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_R7, 9U, c3_e_sf_marshallOut,
    c3_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_eta, 10U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_e, 11U, c3_sf_marshallOut,
    c3_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_Q, 12U, c3_d_sf_marshallOut,
    c3_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_nargin, 13U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_nargout, 14U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c3_q, 15U, c3_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_x, 16U, c3_b_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_xp, 17U, c3_sf_marshallOut,
    c3_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_xe, 18U, c3_sf_marshallOut,
    c3_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 2);
  c3_q1 = c3_q[0];
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 2);
  c3_q2 = c3_q[1];
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 2);
  c3_q3 = c3_q[2];
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 2);
  c3_q4 = c3_q[3];
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 2);
  c3_q5 = c3_q[4];
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 2);
  c3_q6 = c3_q[5];
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 2);
  c3_q7 = c3_q[6];
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 3);
  c3_jacobin(chartInstance, c3_q1, c3_q2, c3_q3, c3_q4, c3_q5, c3_q6, c3_q7,
             c3_unusedU0, c3_b_Te);
  for (c3_i8 = 0; c3_i8 < 16; c3_i8++) {
    c3_Te[c3_i8] = c3_b_Te[c3_i8];
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 3);
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 4);
  for (c3_i9 = 0; c3_i9 < 3; c3_i9++) {
    c3_p[c3_i9] = c3_Te[c3_i9 + 12];
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 5);
  c3_i10 = 0;
  c3_i11 = 0;
  for (c3_i12 = 0; c3_i12 < 3; c3_i12++) {
    for (c3_i13 = 0; c3_i13 < 3; c3_i13++) {
      c3_R7[c3_i13 + c3_i10] = c3_Te[c3_i13 + c3_i11];
    }

    c3_i10 += 3;
    c3_i11 += 4;
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 15);
  c3_b_x = ((c3_R7[0] + c3_R7[4]) + c3_R7[8]) + 1.0;
  c3_c_x = c3_b_x;
  if (c3_c_x < 0.0) {
    c3_eml_error(chartInstance);
  }

  c3_c_x = muDoubleScalarSqrt(c3_c_x);
  c3_eta = 0.5 * c3_c_x;
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 16);
  c3_d_x = c3_R7[5] - c3_R7[7];
  c3_e_x = c3_d_x;
  c3_e_x = muDoubleScalarSign(c3_e_x);
  c3_f_x = ((c3_R7[0] - c3_R7[4]) - c3_R7[8]) + 1.0;
  c3_g_x = c3_f_x;
  if (c3_g_x < 0.0) {
    c3_eml_error(chartInstance);
  }

  c3_g_x = muDoubleScalarSqrt(c3_g_x);
  c3_h_x = c3_R7[6] - c3_R7[2];
  c3_i_x = c3_h_x;
  c3_i_x = muDoubleScalarSign(c3_i_x);
  c3_j_x = ((c3_R7[4] - c3_R7[8]) - c3_R7[0]) + 1.0;
  c3_k_x = c3_j_x;
  if (c3_k_x < 0.0) {
    c3_eml_error(chartInstance);
  }

  c3_k_x = muDoubleScalarSqrt(c3_k_x);
  c3_l_x = c3_R7[1] - c3_R7[3];
  c3_m_x = c3_l_x;
  c3_m_x = muDoubleScalarSign(c3_m_x);
  c3_n_x = ((c3_R7[8] - c3_R7[0]) - c3_R7[4]) + 1.0;
  c3_o_x = c3_n_x;
  if (c3_o_x < 0.0) {
    c3_eml_error(chartInstance);
  }

  c3_o_x = muDoubleScalarSqrt(c3_o_x);
  c3_b[0] = c3_e_x * c3_g_x;
  c3_b[1] = c3_i_x * c3_k_x;
  c3_b[2] = c3_m_x * c3_o_x;
  for (c3_i14 = 0; c3_i14 < 3; c3_i14++) {
    c3_e[c3_i14] = 0.5 * c3_b[c3_i14];
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 19);
  c3_Q[0] = c3_eta;
  for (c3_i15 = 0; c3_i15 < 3; c3_i15++) {
    c3_Q[c3_i15 + 1] = c3_e[c3_i15];
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 20);
  for (c3_i16 = 0; c3_i16 < 3; c3_i16++) {
    c3_x[c3_i16] = c3_p[c3_i16];
  }

  for (c3_i17 = 0; c3_i17 < 4; c3_i17++) {
    c3_x[c3_i17 + 3] = c3_Q[c3_i17];
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 21);
  for (c3_i18 = 0; c3_i18 < 3; c3_i18++) {
    c3_xp[c3_i18] = c3_x[c3_i18];
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 22);
  for (c3_i19 = 0; c3_i19 < 3; c3_i19++) {
    c3_xe[c3_i19] = c3_x[c3_i19 + 3];
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, -22);
  _SFD_SYMBOL_SCOPE_POP();
  for (c3_i20 = 0; c3_i20 < 7; c3_i20++) {
    (*c3_p_x)[c3_i20] = c3_x[c3_i20];
  }

  for (c3_i21 = 0; c3_i21 < 3; c3_i21++) {
    (*c3_b_xp)[c3_i21] = c3_xp[c3_i21];
  }

  for (c3_i22 = 0; c3_i22 < 3; c3_i22++) {
    (*c3_b_xe)[c3_i22] = c3_xe[c3_i22];
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 2U, chartInstance->c3_sfEvent);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_simlwrkuka_kinematicsMachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
  for (c3_i23 = 0; c3_i23 < 7; c3_i23++) {
    _SFD_DATA_RANGE_CHECK((*c3_p_x)[c3_i23], 1U);
  }

  for (c3_i24 = 0; c3_i24 < 3; c3_i24++) {
    _SFD_DATA_RANGE_CHECK((*c3_b_xp)[c3_i24], 2U);
  }

  for (c3_i25 = 0; c3_i25 < 3; c3_i25++) {
    _SFD_DATA_RANGE_CHECK((*c3_b_xe)[c3_i25], 3U);
  }
}

static void initSimStructsc3_simlwrkuka_kinematics
  (SFc3_simlwrkuka_kinematicsInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c3_jacobin(SFc3_simlwrkuka_kinematicsInstanceStruct *chartInstance,
  real_T c3_q1, real_T c3_q2, real_T c3_q3, real_T c3_q4, real_T c3_q5, real_T
  c3_q6, real_T c3_q7, real_T c3_je[42], real_T c3_Te[16])
{
  uint32_T c3_debug_family_var_map[49];
  real_T c3_d3;
  real_T c3_d5;
  real_T c3_d0;
  real_T c3_d7;
  real_T c3_f[7];
  real_T c3_a[7];
  real_T c3_d[7];
  real_T c3_t[7];
  real_T c3_A0[16];
  real_T c3_A1[16];
  real_T c3_A2[16];
  real_T c3_A3[16];
  real_T c3_A4[16];
  real_T c3_A5[16];
  real_T c3_A6[16];
  real_T c3_A7[16];
  real_T c3_Ae[16];
  real_T c3_T1[16];
  real_T c3_T2[16];
  real_T c3_T3[16];
  real_T c3_T4[16];
  real_T c3_T5[16];
  real_T c3_T6[16];
  real_T c3_z0[3];
  real_T c3_z1[3];
  real_T c3_z2[3];
  real_T c3_z3[3];
  real_T c3_z4[3];
  real_T c3_z5[3];
  real_T c3_z6[3];
  real_T c3_p0[3];
  real_T c3_p1[3];
  real_T c3_p2[3];
  real_T c3_p3[3];
  real_T c3_p4[3];
  real_T c3_p5[3];
  real_T c3_p6[3];
  real_T c3_pe[3];
  real_T c3_nargin = 7.0;
  real_T c3_nargout = 2.0;
  int32_T c3_i26;
  static real_T c3_dv3[7] = { 1.5707963267948966, -1.5707963267948966,
    -1.5707963267948966, 1.5707963267948966, 1.5707963267948966,
    -1.5707963267948966, 0.0 };

  int32_T c3_i27;
  int32_T c3_i28;
  static real_T c3_b_a[16] = { 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.31, 1.0 };

  real_T c3_x;
  real_T c3_b_x;
  real_T c3_c_x;
  real_T c3_d_x;
  real_T c3_e_x;
  real_T c3_f_x;
  real_T c3_g_x;
  real_T c3_h_x;
  real_T c3_i_x;
  real_T c3_j_x;
  real_T c3_k_x;
  real_T c3_l_x;
  real_T c3_m_x;
  real_T c3_n_x;
  real_T c3_o_x;
  real_T c3_p_x;
  int32_T c3_i29;
  int32_T c3_i30;
  static real_T c3_dv4[4] = { 0.0, 0.0, 0.0, 1.0 };

  real_T c3_q_x;
  real_T c3_r_x;
  real_T c3_s_x;
  real_T c3_t_x;
  real_T c3_u_x;
  real_T c3_v_x;
  real_T c3_w_x;
  real_T c3_x_x;
  real_T c3_y_x;
  real_T c3_ab_x;
  real_T c3_bb_x;
  real_T c3_cb_x;
  real_T c3_db_x;
  real_T c3_eb_x;
  real_T c3_fb_x;
  real_T c3_gb_x;
  int32_T c3_i31;
  int32_T c3_i32;
  real_T c3_hb_x;
  real_T c3_ib_x;
  real_T c3_jb_x;
  real_T c3_kb_x;
  real_T c3_lb_x;
  real_T c3_mb_x;
  real_T c3_nb_x;
  real_T c3_ob_x;
  real_T c3_pb_x;
  real_T c3_qb_x;
  real_T c3_rb_x;
  real_T c3_sb_x;
  real_T c3_tb_x;
  real_T c3_ub_x;
  real_T c3_vb_x;
  real_T c3_wb_x;
  int32_T c3_i33;
  int32_T c3_i34;
  real_T c3_xb_x;
  real_T c3_yb_x;
  real_T c3_ac_x;
  real_T c3_bc_x;
  real_T c3_cc_x;
  real_T c3_dc_x;
  real_T c3_ec_x;
  real_T c3_fc_x;
  real_T c3_gc_x;
  real_T c3_hc_x;
  real_T c3_ic_x;
  real_T c3_jc_x;
  real_T c3_kc_x;
  real_T c3_lc_x;
  real_T c3_mc_x;
  real_T c3_nc_x;
  int32_T c3_i35;
  int32_T c3_i36;
  real_T c3_oc_x;
  real_T c3_pc_x;
  real_T c3_qc_x;
  real_T c3_rc_x;
  real_T c3_sc_x;
  real_T c3_tc_x;
  real_T c3_uc_x;
  real_T c3_vc_x;
  real_T c3_wc_x;
  real_T c3_xc_x;
  real_T c3_yc_x;
  real_T c3_ad_x;
  real_T c3_bd_x;
  real_T c3_cd_x;
  real_T c3_dd_x;
  real_T c3_ed_x;
  int32_T c3_i37;
  int32_T c3_i38;
  real_T c3_fd_x;
  real_T c3_gd_x;
  real_T c3_hd_x;
  real_T c3_id_x;
  real_T c3_jd_x;
  real_T c3_kd_x;
  real_T c3_ld_x;
  real_T c3_md_x;
  real_T c3_nd_x;
  real_T c3_od_x;
  real_T c3_pd_x;
  real_T c3_qd_x;
  real_T c3_rd_x;
  real_T c3_sd_x;
  real_T c3_td_x;
  real_T c3_ud_x;
  int32_T c3_i39;
  int32_T c3_i40;
  real_T c3_vd_x;
  real_T c3_wd_x;
  real_T c3_xd_x;
  real_T c3_yd_x;
  real_T c3_ae_x;
  real_T c3_be_x;
  real_T c3_ce_x;
  real_T c3_de_x;
  real_T c3_ee_x;
  real_T c3_fe_x;
  real_T c3_ge_x;
  real_T c3_he_x;
  real_T c3_ie_x;
  real_T c3_je_x;
  real_T c3_ke_x;
  real_T c3_le_x;
  int32_T c3_i41;
  int32_T c3_i42;
  int32_T c3_i43;
  static real_T c3_b[16] = { 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.078, 1.0 };

  int32_T c3_i44;
  real_T c3_b_b[16];
  int32_T c3_i45;
  int32_T c3_i46;
  int32_T c3_i47;
  real_T c3_dv5[16];
  int32_T c3_i48;
  real_T c3_dv6[16];
  int32_T c3_i49;
  real_T c3_dv7[16];
  int32_T c3_i50;
  real_T c3_dv8[16];
  int32_T c3_i51;
  int32_T c3_i52;
  real_T c3_y[16];
  int32_T c3_i53;
  real_T c3_c_a[16];
  int32_T c3_i54;
  real_T c3_c_b[16];
  int32_T c3_i55;
  int32_T c3_i56;
  int32_T c3_i57;
  int32_T c3_i58;
  real_T c3_dv9[16];
  int32_T c3_i59;
  real_T c3_dv10[16];
  int32_T c3_i60;
  real_T c3_dv11[16];
  int32_T c3_i61;
  real_T c3_dv12[16];
  int32_T c3_i62;
  int32_T c3_i63;
  int32_T c3_i64;
  real_T c3_d_a[16];
  int32_T c3_i65;
  real_T c3_d_b[16];
  int32_T c3_i66;
  int32_T c3_i67;
  real_T c3_b_y[16];
  int32_T c3_i68;
  real_T c3_c_y[16];
  int32_T c3_i69;
  real_T c3_e_b[16];
  int32_T c3_i70;
  int32_T c3_i71;
  int32_T c3_i72;
  int32_T c3_i73;
  real_T c3_dv13[16];
  int32_T c3_i74;
  real_T c3_dv14[16];
  int32_T c3_i75;
  real_T c3_dv15[16];
  int32_T c3_i76;
  real_T c3_dv16[16];
  int32_T c3_i77;
  int32_T c3_i78;
  int32_T c3_i79;
  real_T c3_e_a[16];
  int32_T c3_i80;
  real_T c3_f_b[16];
  int32_T c3_i81;
  int32_T c3_i82;
  int32_T c3_i83;
  real_T c3_d_y[16];
  int32_T c3_i84;
  real_T c3_g_b[16];
  int32_T c3_i85;
  int32_T c3_i86;
  int32_T c3_i87;
  real_T c3_e_y[16];
  int32_T c3_i88;
  real_T c3_h_b[16];
  int32_T c3_i89;
  int32_T c3_i90;
  int32_T c3_i91;
  int32_T c3_i92;
  real_T c3_dv17[16];
  int32_T c3_i93;
  real_T c3_dv18[16];
  int32_T c3_i94;
  real_T c3_dv19[16];
  int32_T c3_i95;
  real_T c3_dv20[16];
  int32_T c3_i96;
  int32_T c3_i97;
  int32_T c3_i98;
  real_T c3_f_a[16];
  int32_T c3_i99;
  real_T c3_i_b[16];
  int32_T c3_i100;
  int32_T c3_i101;
  int32_T c3_i102;
  real_T c3_f_y[16];
  int32_T c3_i103;
  real_T c3_j_b[16];
  int32_T c3_i104;
  int32_T c3_i105;
  int32_T c3_i106;
  real_T c3_g_y[16];
  int32_T c3_i107;
  real_T c3_k_b[16];
  int32_T c3_i108;
  int32_T c3_i109;
  int32_T c3_i110;
  real_T c3_h_y[16];
  int32_T c3_i111;
  real_T c3_l_b[16];
  int32_T c3_i112;
  int32_T c3_i113;
  int32_T c3_i114;
  int32_T c3_i115;
  real_T c3_dv21[16];
  int32_T c3_i116;
  real_T c3_dv22[16];
  int32_T c3_i117;
  real_T c3_dv23[16];
  int32_T c3_i118;
  real_T c3_dv24[16];
  int32_T c3_i119;
  int32_T c3_i120;
  int32_T c3_i121;
  real_T c3_g_a[16];
  int32_T c3_i122;
  real_T c3_m_b[16];
  int32_T c3_i123;
  int32_T c3_i124;
  int32_T c3_i125;
  real_T c3_i_y[16];
  int32_T c3_i126;
  real_T c3_n_b[16];
  int32_T c3_i127;
  int32_T c3_i128;
  int32_T c3_i129;
  real_T c3_j_y[16];
  int32_T c3_i130;
  real_T c3_o_b[16];
  int32_T c3_i131;
  int32_T c3_i132;
  int32_T c3_i133;
  real_T c3_k_y[16];
  int32_T c3_i134;
  real_T c3_p_b[16];
  int32_T c3_i135;
  int32_T c3_i136;
  int32_T c3_i137;
  real_T c3_l_y[16];
  int32_T c3_i138;
  real_T c3_q_b[16];
  int32_T c3_i139;
  int32_T c3_i140;
  int32_T c3_i141;
  int32_T c3_i142;
  real_T c3_dv25[16];
  int32_T c3_i143;
  real_T c3_dv26[16];
  int32_T c3_i144;
  real_T c3_dv27[16];
  int32_T c3_i145;
  real_T c3_dv28[16];
  int32_T c3_i146;
  int32_T c3_i147;
  int32_T c3_i148;
  real_T c3_h_a[16];
  int32_T c3_i149;
  real_T c3_r_b[16];
  int32_T c3_i150;
  int32_T c3_i151;
  int32_T c3_i152;
  real_T c3_m_y[16];
  int32_T c3_i153;
  real_T c3_s_b[16];
  int32_T c3_i154;
  int32_T c3_i155;
  int32_T c3_i156;
  real_T c3_n_y[16];
  int32_T c3_i157;
  real_T c3_t_b[16];
  int32_T c3_i158;
  int32_T c3_i159;
  int32_T c3_i160;
  real_T c3_o_y[16];
  int32_T c3_i161;
  real_T c3_u_b[16];
  int32_T c3_i162;
  int32_T c3_i163;
  int32_T c3_i164;
  real_T c3_p_y[16];
  int32_T c3_i165;
  real_T c3_v_b[16];
  int32_T c3_i166;
  int32_T c3_i167;
  int32_T c3_i168;
  real_T c3_q_y[16];
  int32_T c3_i169;
  real_T c3_w_b[16];
  int32_T c3_i170;
  int32_T c3_i171;
  int32_T c3_i172;
  real_T c3_r_y[16];
  int32_T c3_i173;
  real_T c3_x_b[16];
  int32_T c3_i174;
  int32_T c3_i175;
  int32_T c3_i176;
  real_T c3_dv29[16];
  int32_T c3_i177;
  real_T c3_dv30[16];
  int32_T c3_i178;
  real_T c3_dv31[16];
  int32_T c3_i179;
  real_T c3_dv32[16];
  int32_T c3_i180;
  static real_T c3_i_a[3] = { 0.0, 0.0, 1.0 };

  int32_T c3_i181;
  int32_T c3_i182;
  int32_T c3_i183;
  int32_T c3_i184;
  int32_T c3_i185;
  int32_T c3_i186;
  int32_T c3_i187;
  int32_T c3_i188;
  int32_T c3_i189;
  int32_T c3_i190;
  int32_T c3_i191;
  int32_T c3_i192;
  int32_T c3_i193;
  int32_T c3_i194;
  int32_T c3_i195;
  real_T c3_y_b[3];
  real_T c3_c1;
  real_T c3_c2;
  real_T c3_c3;
  real_T c3_dv33[3];
  int32_T c3_i196;
  real_T c3_j_a[3];
  int32_T c3_i197;
  real_T c3_b_c1;
  real_T c3_b_c2;
  real_T c3_b_c3;
  real_T c3_dv34[3];
  int32_T c3_i198;
  int32_T c3_i199;
  real_T c3_c_c1;
  real_T c3_c_c2;
  real_T c3_c_c3;
  real_T c3_dv35[3];
  int32_T c3_i200;
  int32_T c3_i201;
  real_T c3_d_c1;
  real_T c3_d_c2;
  real_T c3_d_c3;
  real_T c3_dv36[3];
  int32_T c3_i202;
  int32_T c3_i203;
  real_T c3_e_c1;
  real_T c3_e_c2;
  real_T c3_e_c3;
  real_T c3_dv37[3];
  int32_T c3_i204;
  int32_T c3_i205;
  real_T c3_f_c1;
  real_T c3_f_c2;
  real_T c3_f_c3;
  real_T c3_dv38[3];
  int32_T c3_i206;
  int32_T c3_i207;
  real_T c3_g_c1;
  real_T c3_g_c2;
  real_T c3_g_c3;
  int32_T c3_i208;
  int32_T c3_i209;
  int32_T c3_i210;
  int32_T c3_i211;
  int32_T c3_i212;
  int32_T c3_i213;
  int32_T c3_i214;
  int32_T c3_i215;
  int32_T c3_i216;
  int32_T c3_i217;
  int32_T c3_i218;
  int32_T c3_i219;
  int32_T c3_i220;
  int32_T c3_i221;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 49U, 49U, c3_b_debug_family_names,
    c3_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_d3, 0U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_d5, 1U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c3_d0, 2U, c3_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c3_d7, 3U, c3_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c3_f, 4U, c3_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c3_a, 5U, c3_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_d, 6U, c3_b_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_t, 7U, c3_b_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c3_A0, 8U, c3_f_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_A1, 9U, c3_f_sf_marshallOut,
    c3_f_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_A2, 10U, c3_f_sf_marshallOut,
    c3_f_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_A3, 11U, c3_f_sf_marshallOut,
    c3_f_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_A4, 12U, c3_f_sf_marshallOut,
    c3_f_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_A5, 13U, c3_f_sf_marshallOut,
    c3_f_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_A6, 14U, c3_f_sf_marshallOut,
    c3_f_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_A7, 15U, c3_f_sf_marshallOut,
    c3_f_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c3_Ae, 16U, c3_f_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_T1, 17U, c3_f_sf_marshallOut,
    c3_f_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_T2, 18U, c3_f_sf_marshallOut,
    c3_f_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_T3, 19U, c3_f_sf_marshallOut,
    c3_f_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_T4, 20U, c3_f_sf_marshallOut,
    c3_f_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_T5, 21U, c3_f_sf_marshallOut,
    c3_f_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_T6, 22U, c3_f_sf_marshallOut,
    c3_f_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c3_z0, 23U, c3_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_z1, 24U, c3_sf_marshallOut,
    c3_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_z2, 25U, c3_sf_marshallOut,
    c3_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_z3, 26U, c3_sf_marshallOut,
    c3_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_z4, 27U, c3_sf_marshallOut,
    c3_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_z5, 28U, c3_sf_marshallOut,
    c3_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_z6, 29U, c3_sf_marshallOut,
    c3_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_p0, 30U, c3_sf_marshallOut,
    c3_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_p1, 31U, c3_sf_marshallOut,
    c3_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_p2, 32U, c3_sf_marshallOut,
    c3_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_p3, 33U, c3_sf_marshallOut,
    c3_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_p4, 34U, c3_sf_marshallOut,
    c3_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_p5, 35U, c3_sf_marshallOut,
    c3_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_p6, 36U, c3_sf_marshallOut,
    c3_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_pe, 37U, c3_sf_marshallOut,
    c3_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_nargin, 38U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_nargout, 39U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_q1, 40U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_q2, 41U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_q3, 42U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_q4, 43U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_q5, 44U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_q6, 45U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_q7, 46U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_je, 47U, c3_g_sf_marshallOut,
    c3_g_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_Te, 48U, c3_f_sf_marshallOut,
    c3_f_sf_marshallIn);
  CV_SCRIPT_FCN(0, 0);
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 2);
  c3_d3 = 0.4;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 3);
  c3_d5 = 0.39;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 4);
  c3_d0 = 0.31;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 5);
  c3_d7 = 0.078;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 6);
  for (c3_i26 = 0; c3_i26 < 7; c3_i26++) {
    c3_f[c3_i26] = c3_dv3[c3_i26];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 7);
  for (c3_i27 = 0; c3_i27 < 7; c3_i27++) {
    c3_a[c3_i27] = 0.0;
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 8);
  c3_d[0] = 0.0;
  c3_d[1] = 0.0;
  c3_d[2] = c3_d3;
  c3_d[3] = 0.0;
  c3_d[4] = c3_d5;
  c3_d[5] = 0.0;
  c3_d[6] = 0.0;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 9);
  c3_t[0] = c3_q1;
  c3_t[1] = c3_q2;
  c3_t[2] = c3_q3;
  c3_t[3] = c3_q4;
  c3_t[4] = c3_q5;
  c3_t[5] = c3_q6;
  c3_t[6] = c3_q7;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 11);
  for (c3_i28 = 0; c3_i28 < 16; c3_i28++) {
    c3_A0[c3_i28] = c3_b_a[c3_i28];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 12);
  c3_x = c3_t[0];
  c3_b_x = c3_x;
  c3_b_x = muDoubleScalarCos(c3_b_x);
  c3_c_x = c3_t[0];
  c3_d_x = c3_c_x;
  c3_d_x = muDoubleScalarSin(c3_d_x);
  c3_e_x = c3_t[0];
  c3_f_x = c3_e_x;
  c3_f_x = muDoubleScalarSin(c3_f_x);
  c3_g_x = c3_t[0];
  c3_h_x = c3_g_x;
  c3_h_x = muDoubleScalarCos(c3_h_x);
  c3_i_x = c3_t[0];
  c3_j_x = c3_i_x;
  c3_j_x = muDoubleScalarSin(c3_j_x);
  c3_k_x = c3_t[0];
  c3_l_x = c3_k_x;
  c3_l_x = muDoubleScalarCos(c3_l_x);
  c3_m_x = c3_t[0];
  c3_n_x = c3_m_x;
  c3_n_x = muDoubleScalarCos(c3_n_x);
  c3_o_x = c3_t[0];
  c3_p_x = c3_o_x;
  c3_p_x = muDoubleScalarSin(c3_p_x);
  c3_A1[0] = c3_b_x;
  c3_A1[4] = -c3_d_x * 6.123233995736766E-17;
  c3_A1[8] = c3_f_x;
  c3_A1[12] = 0.0 * c3_h_x;
  c3_A1[1] = c3_j_x;
  c3_A1[5] = c3_l_x * 6.123233995736766E-17;
  c3_A1[9] = -c3_n_x;
  c3_A1[13] = 0.0 * c3_p_x;
  c3_A1[2] = 0.0;
  c3_A1[6] = 1.0;
  c3_A1[10] = 6.123233995736766E-17;
  c3_A1[14] = c3_d[0];
  c3_i29 = 0;
  for (c3_i30 = 0; c3_i30 < 4; c3_i30++) {
    c3_A1[c3_i29 + 3] = c3_dv4[c3_i30];
    c3_i29 += 4;
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 13);
  c3_q_x = c3_t[1];
  c3_r_x = c3_q_x;
  c3_r_x = muDoubleScalarCos(c3_r_x);
  c3_s_x = c3_t[1];
  c3_t_x = c3_s_x;
  c3_t_x = muDoubleScalarSin(c3_t_x);
  c3_u_x = c3_t[1];
  c3_v_x = c3_u_x;
  c3_v_x = muDoubleScalarSin(c3_v_x);
  c3_w_x = c3_t[1];
  c3_x_x = c3_w_x;
  c3_x_x = muDoubleScalarCos(c3_x_x);
  c3_y_x = c3_t[1];
  c3_ab_x = c3_y_x;
  c3_ab_x = muDoubleScalarSin(c3_ab_x);
  c3_bb_x = c3_t[1];
  c3_cb_x = c3_bb_x;
  c3_cb_x = muDoubleScalarCos(c3_cb_x);
  c3_db_x = c3_t[1];
  c3_eb_x = c3_db_x;
  c3_eb_x = muDoubleScalarCos(c3_eb_x);
  c3_fb_x = c3_t[1];
  c3_gb_x = c3_fb_x;
  c3_gb_x = muDoubleScalarSin(c3_gb_x);
  c3_A2[0] = c3_r_x;
  c3_A2[4] = -c3_t_x * 6.123233995736766E-17;
  c3_A2[8] = -c3_v_x;
  c3_A2[12] = 0.0 * c3_x_x;
  c3_A2[1] = c3_ab_x;
  c3_A2[5] = c3_cb_x * 6.123233995736766E-17;
  c3_A2[9] = -(-c3_eb_x);
  c3_A2[13] = 0.0 * c3_gb_x;
  c3_A2[2] = 0.0;
  c3_A2[6] = -1.0;
  c3_A2[10] = 6.123233995736766E-17;
  c3_A2[14] = c3_d[1];
  c3_i31 = 0;
  for (c3_i32 = 0; c3_i32 < 4; c3_i32++) {
    c3_A2[c3_i31 + 3] = c3_dv4[c3_i32];
    c3_i31 += 4;
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 14);
  c3_hb_x = c3_t[2];
  c3_ib_x = c3_hb_x;
  c3_ib_x = muDoubleScalarCos(c3_ib_x);
  c3_jb_x = c3_t[2];
  c3_kb_x = c3_jb_x;
  c3_kb_x = muDoubleScalarSin(c3_kb_x);
  c3_lb_x = c3_t[2];
  c3_mb_x = c3_lb_x;
  c3_mb_x = muDoubleScalarSin(c3_mb_x);
  c3_nb_x = c3_t[2];
  c3_ob_x = c3_nb_x;
  c3_ob_x = muDoubleScalarCos(c3_ob_x);
  c3_pb_x = c3_t[2];
  c3_qb_x = c3_pb_x;
  c3_qb_x = muDoubleScalarSin(c3_qb_x);
  c3_rb_x = c3_t[2];
  c3_sb_x = c3_rb_x;
  c3_sb_x = muDoubleScalarCos(c3_sb_x);
  c3_tb_x = c3_t[2];
  c3_ub_x = c3_tb_x;
  c3_ub_x = muDoubleScalarCos(c3_ub_x);
  c3_vb_x = c3_t[2];
  c3_wb_x = c3_vb_x;
  c3_wb_x = muDoubleScalarSin(c3_wb_x);
  c3_A3[0] = c3_ib_x;
  c3_A3[4] = -c3_kb_x * 6.123233995736766E-17;
  c3_A3[8] = -c3_mb_x;
  c3_A3[12] = 0.0 * c3_ob_x;
  c3_A3[1] = c3_qb_x;
  c3_A3[5] = c3_sb_x * 6.123233995736766E-17;
  c3_A3[9] = -(-c3_ub_x);
  c3_A3[13] = 0.0 * c3_wb_x;
  c3_A3[2] = 0.0;
  c3_A3[6] = -1.0;
  c3_A3[10] = 6.123233995736766E-17;
  c3_A3[14] = c3_d[2];
  c3_i33 = 0;
  for (c3_i34 = 0; c3_i34 < 4; c3_i34++) {
    c3_A3[c3_i33 + 3] = c3_dv4[c3_i34];
    c3_i33 += 4;
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 15);
  c3_xb_x = c3_t[3];
  c3_yb_x = c3_xb_x;
  c3_yb_x = muDoubleScalarCos(c3_yb_x);
  c3_ac_x = c3_t[3];
  c3_bc_x = c3_ac_x;
  c3_bc_x = muDoubleScalarSin(c3_bc_x);
  c3_cc_x = c3_t[3];
  c3_dc_x = c3_cc_x;
  c3_dc_x = muDoubleScalarSin(c3_dc_x);
  c3_ec_x = c3_t[3];
  c3_fc_x = c3_ec_x;
  c3_fc_x = muDoubleScalarCos(c3_fc_x);
  c3_gc_x = c3_t[3];
  c3_hc_x = c3_gc_x;
  c3_hc_x = muDoubleScalarSin(c3_hc_x);
  c3_ic_x = c3_t[3];
  c3_jc_x = c3_ic_x;
  c3_jc_x = muDoubleScalarCos(c3_jc_x);
  c3_kc_x = c3_t[3];
  c3_lc_x = c3_kc_x;
  c3_lc_x = muDoubleScalarCos(c3_lc_x);
  c3_mc_x = c3_t[3];
  c3_nc_x = c3_mc_x;
  c3_nc_x = muDoubleScalarSin(c3_nc_x);
  c3_A4[0] = c3_yb_x;
  c3_A4[4] = -c3_bc_x * 6.123233995736766E-17;
  c3_A4[8] = c3_dc_x;
  c3_A4[12] = 0.0 * c3_fc_x;
  c3_A4[1] = c3_hc_x;
  c3_A4[5] = c3_jc_x * 6.123233995736766E-17;
  c3_A4[9] = -c3_lc_x;
  c3_A4[13] = 0.0 * c3_nc_x;
  c3_A4[2] = 0.0;
  c3_A4[6] = 1.0;
  c3_A4[10] = 6.123233995736766E-17;
  c3_A4[14] = c3_d[3];
  c3_i35 = 0;
  for (c3_i36 = 0; c3_i36 < 4; c3_i36++) {
    c3_A4[c3_i35 + 3] = c3_dv4[c3_i36];
    c3_i35 += 4;
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 16);
  c3_oc_x = c3_t[4];
  c3_pc_x = c3_oc_x;
  c3_pc_x = muDoubleScalarCos(c3_pc_x);
  c3_qc_x = c3_t[4];
  c3_rc_x = c3_qc_x;
  c3_rc_x = muDoubleScalarSin(c3_rc_x);
  c3_sc_x = c3_t[4];
  c3_tc_x = c3_sc_x;
  c3_tc_x = muDoubleScalarSin(c3_tc_x);
  c3_uc_x = c3_t[4];
  c3_vc_x = c3_uc_x;
  c3_vc_x = muDoubleScalarCos(c3_vc_x);
  c3_wc_x = c3_t[4];
  c3_xc_x = c3_wc_x;
  c3_xc_x = muDoubleScalarSin(c3_xc_x);
  c3_yc_x = c3_t[4];
  c3_ad_x = c3_yc_x;
  c3_ad_x = muDoubleScalarCos(c3_ad_x);
  c3_bd_x = c3_t[4];
  c3_cd_x = c3_bd_x;
  c3_cd_x = muDoubleScalarCos(c3_cd_x);
  c3_dd_x = c3_t[4];
  c3_ed_x = c3_dd_x;
  c3_ed_x = muDoubleScalarSin(c3_ed_x);
  c3_A5[0] = c3_pc_x;
  c3_A5[4] = -c3_rc_x * 6.123233995736766E-17;
  c3_A5[8] = c3_tc_x;
  c3_A5[12] = 0.0 * c3_vc_x;
  c3_A5[1] = c3_xc_x;
  c3_A5[5] = c3_ad_x * 6.123233995736766E-17;
  c3_A5[9] = -c3_cd_x;
  c3_A5[13] = 0.0 * c3_ed_x;
  c3_A5[2] = 0.0;
  c3_A5[6] = 1.0;
  c3_A5[10] = 6.123233995736766E-17;
  c3_A5[14] = c3_d[4];
  c3_i37 = 0;
  for (c3_i38 = 0; c3_i38 < 4; c3_i38++) {
    c3_A5[c3_i37 + 3] = c3_dv4[c3_i38];
    c3_i37 += 4;
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 17);
  c3_fd_x = c3_t[5];
  c3_gd_x = c3_fd_x;
  c3_gd_x = muDoubleScalarCos(c3_gd_x);
  c3_hd_x = c3_t[5];
  c3_id_x = c3_hd_x;
  c3_id_x = muDoubleScalarSin(c3_id_x);
  c3_jd_x = c3_t[5];
  c3_kd_x = c3_jd_x;
  c3_kd_x = muDoubleScalarSin(c3_kd_x);
  c3_ld_x = c3_t[5];
  c3_md_x = c3_ld_x;
  c3_md_x = muDoubleScalarCos(c3_md_x);
  c3_nd_x = c3_t[5];
  c3_od_x = c3_nd_x;
  c3_od_x = muDoubleScalarSin(c3_od_x);
  c3_pd_x = c3_t[5];
  c3_qd_x = c3_pd_x;
  c3_qd_x = muDoubleScalarCos(c3_qd_x);
  c3_rd_x = c3_t[5];
  c3_sd_x = c3_rd_x;
  c3_sd_x = muDoubleScalarCos(c3_sd_x);
  c3_td_x = c3_t[5];
  c3_ud_x = c3_td_x;
  c3_ud_x = muDoubleScalarSin(c3_ud_x);
  c3_A6[0] = c3_gd_x;
  c3_A6[4] = -c3_id_x * 6.123233995736766E-17;
  c3_A6[8] = -c3_kd_x;
  c3_A6[12] = 0.0 * c3_md_x;
  c3_A6[1] = c3_od_x;
  c3_A6[5] = c3_qd_x * 6.123233995736766E-17;
  c3_A6[9] = -(-c3_sd_x);
  c3_A6[13] = 0.0 * c3_ud_x;
  c3_A6[2] = 0.0;
  c3_A6[6] = -1.0;
  c3_A6[10] = 6.123233995736766E-17;
  c3_A6[14] = c3_d[5];
  c3_i39 = 0;
  for (c3_i40 = 0; c3_i40 < 4; c3_i40++) {
    c3_A6[c3_i39 + 3] = c3_dv4[c3_i40];
    c3_i39 += 4;
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 18);
  c3_vd_x = c3_t[6];
  c3_wd_x = c3_vd_x;
  c3_wd_x = muDoubleScalarCos(c3_wd_x);
  c3_xd_x = c3_t[6];
  c3_yd_x = c3_xd_x;
  c3_yd_x = muDoubleScalarSin(c3_yd_x);
  c3_ae_x = c3_t[6];
  c3_be_x = c3_ae_x;
  c3_be_x = muDoubleScalarSin(c3_be_x);
  c3_ce_x = c3_t[6];
  c3_de_x = c3_ce_x;
  c3_de_x = muDoubleScalarCos(c3_de_x);
  c3_ee_x = c3_t[6];
  c3_fe_x = c3_ee_x;
  c3_fe_x = muDoubleScalarSin(c3_fe_x);
  c3_ge_x = c3_t[6];
  c3_he_x = c3_ge_x;
  c3_he_x = muDoubleScalarCos(c3_he_x);
  c3_ie_x = c3_t[6];
  c3_je_x = c3_ie_x;
  c3_je_x = muDoubleScalarCos(c3_je_x);
  c3_ke_x = c3_t[6];
  c3_le_x = c3_ke_x;
  c3_le_x = muDoubleScalarSin(c3_le_x);
  c3_A7[0] = c3_wd_x;
  c3_A7[4] = -c3_yd_x;
  c3_A7[8] = c3_be_x * 0.0;
  c3_A7[12] = 0.0 * c3_de_x;
  c3_A7[1] = c3_fe_x;
  c3_A7[5] = c3_he_x;
  c3_A7[9] = -c3_je_x * 0.0;
  c3_A7[13] = 0.0 * c3_le_x;
  c3_A7[2] = 0.0;
  c3_A7[6] = 0.0;
  c3_A7[10] = 1.0;
  c3_A7[14] = c3_d[6];
  c3_i41 = 0;
  for (c3_i42 = 0; c3_i42 < 4; c3_i42++) {
    c3_A7[c3_i41 + 3] = c3_dv4[c3_i42];
    c3_i41 += 4;
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 19);
  for (c3_i43 = 0; c3_i43 < 16; c3_i43++) {
    c3_Ae[c3_i43] = c3_b[c3_i43];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 21);
  for (c3_i44 = 0; c3_i44 < 16; c3_i44++) {
    c3_b_b[c3_i44] = c3_A1[c3_i44];
  }

  c3_eml_scalar_eg(chartInstance);
  c3_eml_scalar_eg(chartInstance);
  for (c3_i45 = 0; c3_i45 < 16; c3_i45++) {
    c3_T1[c3_i45] = 0.0;
  }

  for (c3_i46 = 0; c3_i46 < 16; c3_i46++) {
    c3_T1[c3_i46] = 0.0;
  }

  for (c3_i47 = 0; c3_i47 < 16; c3_i47++) {
    c3_dv5[c3_i47] = c3_b_a[c3_i47];
  }

  for (c3_i48 = 0; c3_i48 < 16; c3_i48++) {
    c3_dv6[c3_i48] = c3_b_b[c3_i48];
  }

  for (c3_i49 = 0; c3_i49 < 16; c3_i49++) {
    c3_dv7[c3_i49] = c3_dv5[c3_i49];
  }

  for (c3_i50 = 0; c3_i50 < 16; c3_i50++) {
    c3_dv8[c3_i50] = c3_dv6[c3_i50];
  }

  c3_b_eml_xgemm(chartInstance, c3_dv7, c3_dv8, c3_T1);
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 22);
  for (c3_i51 = 0; c3_i51 < 16; c3_i51++) {
    c3_b_b[c3_i51] = c3_A1[c3_i51];
  }

  c3_eml_scalar_eg(chartInstance);
  c3_eml_scalar_eg(chartInstance);
  for (c3_i52 = 0; c3_i52 < 16; c3_i52++) {
    c3_y[c3_i52] = 0.0;
  }

  for (c3_i53 = 0; c3_i53 < 16; c3_i53++) {
    c3_c_a[c3_i53] = c3_b_a[c3_i53];
  }

  for (c3_i54 = 0; c3_i54 < 16; c3_i54++) {
    c3_c_b[c3_i54] = c3_b_b[c3_i54];
  }

  c3_b_eml_xgemm(chartInstance, c3_c_a, c3_c_b, c3_y);
  for (c3_i55 = 0; c3_i55 < 16; c3_i55++) {
    c3_b_b[c3_i55] = c3_A2[c3_i55];
  }

  c3_eml_scalar_eg(chartInstance);
  c3_eml_scalar_eg(chartInstance);
  for (c3_i56 = 0; c3_i56 < 16; c3_i56++) {
    c3_T2[c3_i56] = 0.0;
  }

  for (c3_i57 = 0; c3_i57 < 16; c3_i57++) {
    c3_T2[c3_i57] = 0.0;
  }

  for (c3_i58 = 0; c3_i58 < 16; c3_i58++) {
    c3_dv9[c3_i58] = c3_y[c3_i58];
  }

  for (c3_i59 = 0; c3_i59 < 16; c3_i59++) {
    c3_dv10[c3_i59] = c3_b_b[c3_i59];
  }

  for (c3_i60 = 0; c3_i60 < 16; c3_i60++) {
    c3_dv11[c3_i60] = c3_dv9[c3_i60];
  }

  for (c3_i61 = 0; c3_i61 < 16; c3_i61++) {
    c3_dv12[c3_i61] = c3_dv10[c3_i61];
  }

  c3_b_eml_xgemm(chartInstance, c3_dv11, c3_dv12, c3_T2);
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 23);
  for (c3_i62 = 0; c3_i62 < 16; c3_i62++) {
    c3_b_b[c3_i62] = c3_A1[c3_i62];
  }

  c3_eml_scalar_eg(chartInstance);
  c3_eml_scalar_eg(chartInstance);
  for (c3_i63 = 0; c3_i63 < 16; c3_i63++) {
    c3_y[c3_i63] = 0.0;
  }

  for (c3_i64 = 0; c3_i64 < 16; c3_i64++) {
    c3_d_a[c3_i64] = c3_b_a[c3_i64];
  }

  for (c3_i65 = 0; c3_i65 < 16; c3_i65++) {
    c3_d_b[c3_i65] = c3_b_b[c3_i65];
  }

  c3_b_eml_xgemm(chartInstance, c3_d_a, c3_d_b, c3_y);
  for (c3_i66 = 0; c3_i66 < 16; c3_i66++) {
    c3_b_b[c3_i66] = c3_A2[c3_i66];
  }

  c3_eml_scalar_eg(chartInstance);
  c3_eml_scalar_eg(chartInstance);
  for (c3_i67 = 0; c3_i67 < 16; c3_i67++) {
    c3_b_y[c3_i67] = 0.0;
  }

  for (c3_i68 = 0; c3_i68 < 16; c3_i68++) {
    c3_c_y[c3_i68] = c3_y[c3_i68];
  }

  for (c3_i69 = 0; c3_i69 < 16; c3_i69++) {
    c3_e_b[c3_i69] = c3_b_b[c3_i69];
  }

  c3_b_eml_xgemm(chartInstance, c3_c_y, c3_e_b, c3_b_y);
  for (c3_i70 = 0; c3_i70 < 16; c3_i70++) {
    c3_b_b[c3_i70] = c3_A3[c3_i70];
  }

  c3_eml_scalar_eg(chartInstance);
  c3_eml_scalar_eg(chartInstance);
  for (c3_i71 = 0; c3_i71 < 16; c3_i71++) {
    c3_T3[c3_i71] = 0.0;
  }

  for (c3_i72 = 0; c3_i72 < 16; c3_i72++) {
    c3_T3[c3_i72] = 0.0;
  }

  for (c3_i73 = 0; c3_i73 < 16; c3_i73++) {
    c3_dv13[c3_i73] = c3_b_y[c3_i73];
  }

  for (c3_i74 = 0; c3_i74 < 16; c3_i74++) {
    c3_dv14[c3_i74] = c3_b_b[c3_i74];
  }

  for (c3_i75 = 0; c3_i75 < 16; c3_i75++) {
    c3_dv15[c3_i75] = c3_dv13[c3_i75];
  }

  for (c3_i76 = 0; c3_i76 < 16; c3_i76++) {
    c3_dv16[c3_i76] = c3_dv14[c3_i76];
  }

  c3_b_eml_xgemm(chartInstance, c3_dv15, c3_dv16, c3_T3);
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 24);
  for (c3_i77 = 0; c3_i77 < 16; c3_i77++) {
    c3_b_b[c3_i77] = c3_A1[c3_i77];
  }

  c3_eml_scalar_eg(chartInstance);
  c3_eml_scalar_eg(chartInstance);
  for (c3_i78 = 0; c3_i78 < 16; c3_i78++) {
    c3_y[c3_i78] = 0.0;
  }

  for (c3_i79 = 0; c3_i79 < 16; c3_i79++) {
    c3_e_a[c3_i79] = c3_b_a[c3_i79];
  }

  for (c3_i80 = 0; c3_i80 < 16; c3_i80++) {
    c3_f_b[c3_i80] = c3_b_b[c3_i80];
  }

  c3_b_eml_xgemm(chartInstance, c3_e_a, c3_f_b, c3_y);
  for (c3_i81 = 0; c3_i81 < 16; c3_i81++) {
    c3_b_b[c3_i81] = c3_A2[c3_i81];
  }

  c3_eml_scalar_eg(chartInstance);
  c3_eml_scalar_eg(chartInstance);
  for (c3_i82 = 0; c3_i82 < 16; c3_i82++) {
    c3_b_y[c3_i82] = 0.0;
  }

  for (c3_i83 = 0; c3_i83 < 16; c3_i83++) {
    c3_d_y[c3_i83] = c3_y[c3_i83];
  }

  for (c3_i84 = 0; c3_i84 < 16; c3_i84++) {
    c3_g_b[c3_i84] = c3_b_b[c3_i84];
  }

  c3_b_eml_xgemm(chartInstance, c3_d_y, c3_g_b, c3_b_y);
  for (c3_i85 = 0; c3_i85 < 16; c3_i85++) {
    c3_b_b[c3_i85] = c3_A3[c3_i85];
  }

  c3_eml_scalar_eg(chartInstance);
  c3_eml_scalar_eg(chartInstance);
  for (c3_i86 = 0; c3_i86 < 16; c3_i86++) {
    c3_y[c3_i86] = 0.0;
  }

  for (c3_i87 = 0; c3_i87 < 16; c3_i87++) {
    c3_e_y[c3_i87] = c3_b_y[c3_i87];
  }

  for (c3_i88 = 0; c3_i88 < 16; c3_i88++) {
    c3_h_b[c3_i88] = c3_b_b[c3_i88];
  }

  c3_b_eml_xgemm(chartInstance, c3_e_y, c3_h_b, c3_y);
  for (c3_i89 = 0; c3_i89 < 16; c3_i89++) {
    c3_b_b[c3_i89] = c3_A4[c3_i89];
  }

  c3_eml_scalar_eg(chartInstance);
  c3_eml_scalar_eg(chartInstance);
  for (c3_i90 = 0; c3_i90 < 16; c3_i90++) {
    c3_T4[c3_i90] = 0.0;
  }

  for (c3_i91 = 0; c3_i91 < 16; c3_i91++) {
    c3_T4[c3_i91] = 0.0;
  }

  for (c3_i92 = 0; c3_i92 < 16; c3_i92++) {
    c3_dv17[c3_i92] = c3_y[c3_i92];
  }

  for (c3_i93 = 0; c3_i93 < 16; c3_i93++) {
    c3_dv18[c3_i93] = c3_b_b[c3_i93];
  }

  for (c3_i94 = 0; c3_i94 < 16; c3_i94++) {
    c3_dv19[c3_i94] = c3_dv17[c3_i94];
  }

  for (c3_i95 = 0; c3_i95 < 16; c3_i95++) {
    c3_dv20[c3_i95] = c3_dv18[c3_i95];
  }

  c3_b_eml_xgemm(chartInstance, c3_dv19, c3_dv20, c3_T4);
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 25);
  for (c3_i96 = 0; c3_i96 < 16; c3_i96++) {
    c3_b_b[c3_i96] = c3_A1[c3_i96];
  }

  c3_eml_scalar_eg(chartInstance);
  c3_eml_scalar_eg(chartInstance);
  for (c3_i97 = 0; c3_i97 < 16; c3_i97++) {
    c3_y[c3_i97] = 0.0;
  }

  for (c3_i98 = 0; c3_i98 < 16; c3_i98++) {
    c3_f_a[c3_i98] = c3_b_a[c3_i98];
  }

  for (c3_i99 = 0; c3_i99 < 16; c3_i99++) {
    c3_i_b[c3_i99] = c3_b_b[c3_i99];
  }

  c3_b_eml_xgemm(chartInstance, c3_f_a, c3_i_b, c3_y);
  for (c3_i100 = 0; c3_i100 < 16; c3_i100++) {
    c3_b_b[c3_i100] = c3_A2[c3_i100];
  }

  c3_eml_scalar_eg(chartInstance);
  c3_eml_scalar_eg(chartInstance);
  for (c3_i101 = 0; c3_i101 < 16; c3_i101++) {
    c3_b_y[c3_i101] = 0.0;
  }

  for (c3_i102 = 0; c3_i102 < 16; c3_i102++) {
    c3_f_y[c3_i102] = c3_y[c3_i102];
  }

  for (c3_i103 = 0; c3_i103 < 16; c3_i103++) {
    c3_j_b[c3_i103] = c3_b_b[c3_i103];
  }

  c3_b_eml_xgemm(chartInstance, c3_f_y, c3_j_b, c3_b_y);
  for (c3_i104 = 0; c3_i104 < 16; c3_i104++) {
    c3_b_b[c3_i104] = c3_A3[c3_i104];
  }

  c3_eml_scalar_eg(chartInstance);
  c3_eml_scalar_eg(chartInstance);
  for (c3_i105 = 0; c3_i105 < 16; c3_i105++) {
    c3_y[c3_i105] = 0.0;
  }

  for (c3_i106 = 0; c3_i106 < 16; c3_i106++) {
    c3_g_y[c3_i106] = c3_b_y[c3_i106];
  }

  for (c3_i107 = 0; c3_i107 < 16; c3_i107++) {
    c3_k_b[c3_i107] = c3_b_b[c3_i107];
  }

  c3_b_eml_xgemm(chartInstance, c3_g_y, c3_k_b, c3_y);
  for (c3_i108 = 0; c3_i108 < 16; c3_i108++) {
    c3_b_b[c3_i108] = c3_A4[c3_i108];
  }

  c3_eml_scalar_eg(chartInstance);
  c3_eml_scalar_eg(chartInstance);
  for (c3_i109 = 0; c3_i109 < 16; c3_i109++) {
    c3_b_y[c3_i109] = 0.0;
  }

  for (c3_i110 = 0; c3_i110 < 16; c3_i110++) {
    c3_h_y[c3_i110] = c3_y[c3_i110];
  }

  for (c3_i111 = 0; c3_i111 < 16; c3_i111++) {
    c3_l_b[c3_i111] = c3_b_b[c3_i111];
  }

  c3_b_eml_xgemm(chartInstance, c3_h_y, c3_l_b, c3_b_y);
  for (c3_i112 = 0; c3_i112 < 16; c3_i112++) {
    c3_b_b[c3_i112] = c3_A5[c3_i112];
  }

  c3_eml_scalar_eg(chartInstance);
  c3_eml_scalar_eg(chartInstance);
  for (c3_i113 = 0; c3_i113 < 16; c3_i113++) {
    c3_T5[c3_i113] = 0.0;
  }

  for (c3_i114 = 0; c3_i114 < 16; c3_i114++) {
    c3_T5[c3_i114] = 0.0;
  }

  for (c3_i115 = 0; c3_i115 < 16; c3_i115++) {
    c3_dv21[c3_i115] = c3_b_y[c3_i115];
  }

  for (c3_i116 = 0; c3_i116 < 16; c3_i116++) {
    c3_dv22[c3_i116] = c3_b_b[c3_i116];
  }

  for (c3_i117 = 0; c3_i117 < 16; c3_i117++) {
    c3_dv23[c3_i117] = c3_dv21[c3_i117];
  }

  for (c3_i118 = 0; c3_i118 < 16; c3_i118++) {
    c3_dv24[c3_i118] = c3_dv22[c3_i118];
  }

  c3_b_eml_xgemm(chartInstance, c3_dv23, c3_dv24, c3_T5);
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 26);
  for (c3_i119 = 0; c3_i119 < 16; c3_i119++) {
    c3_b_b[c3_i119] = c3_A1[c3_i119];
  }

  c3_eml_scalar_eg(chartInstance);
  c3_eml_scalar_eg(chartInstance);
  for (c3_i120 = 0; c3_i120 < 16; c3_i120++) {
    c3_y[c3_i120] = 0.0;
  }

  for (c3_i121 = 0; c3_i121 < 16; c3_i121++) {
    c3_g_a[c3_i121] = c3_b_a[c3_i121];
  }

  for (c3_i122 = 0; c3_i122 < 16; c3_i122++) {
    c3_m_b[c3_i122] = c3_b_b[c3_i122];
  }

  c3_b_eml_xgemm(chartInstance, c3_g_a, c3_m_b, c3_y);
  for (c3_i123 = 0; c3_i123 < 16; c3_i123++) {
    c3_b_b[c3_i123] = c3_A2[c3_i123];
  }

  c3_eml_scalar_eg(chartInstance);
  c3_eml_scalar_eg(chartInstance);
  for (c3_i124 = 0; c3_i124 < 16; c3_i124++) {
    c3_b_y[c3_i124] = 0.0;
  }

  for (c3_i125 = 0; c3_i125 < 16; c3_i125++) {
    c3_i_y[c3_i125] = c3_y[c3_i125];
  }

  for (c3_i126 = 0; c3_i126 < 16; c3_i126++) {
    c3_n_b[c3_i126] = c3_b_b[c3_i126];
  }

  c3_b_eml_xgemm(chartInstance, c3_i_y, c3_n_b, c3_b_y);
  for (c3_i127 = 0; c3_i127 < 16; c3_i127++) {
    c3_b_b[c3_i127] = c3_A3[c3_i127];
  }

  c3_eml_scalar_eg(chartInstance);
  c3_eml_scalar_eg(chartInstance);
  for (c3_i128 = 0; c3_i128 < 16; c3_i128++) {
    c3_y[c3_i128] = 0.0;
  }

  for (c3_i129 = 0; c3_i129 < 16; c3_i129++) {
    c3_j_y[c3_i129] = c3_b_y[c3_i129];
  }

  for (c3_i130 = 0; c3_i130 < 16; c3_i130++) {
    c3_o_b[c3_i130] = c3_b_b[c3_i130];
  }

  c3_b_eml_xgemm(chartInstance, c3_j_y, c3_o_b, c3_y);
  for (c3_i131 = 0; c3_i131 < 16; c3_i131++) {
    c3_b_b[c3_i131] = c3_A4[c3_i131];
  }

  c3_eml_scalar_eg(chartInstance);
  c3_eml_scalar_eg(chartInstance);
  for (c3_i132 = 0; c3_i132 < 16; c3_i132++) {
    c3_b_y[c3_i132] = 0.0;
  }

  for (c3_i133 = 0; c3_i133 < 16; c3_i133++) {
    c3_k_y[c3_i133] = c3_y[c3_i133];
  }

  for (c3_i134 = 0; c3_i134 < 16; c3_i134++) {
    c3_p_b[c3_i134] = c3_b_b[c3_i134];
  }

  c3_b_eml_xgemm(chartInstance, c3_k_y, c3_p_b, c3_b_y);
  for (c3_i135 = 0; c3_i135 < 16; c3_i135++) {
    c3_b_b[c3_i135] = c3_A5[c3_i135];
  }

  c3_eml_scalar_eg(chartInstance);
  c3_eml_scalar_eg(chartInstance);
  for (c3_i136 = 0; c3_i136 < 16; c3_i136++) {
    c3_y[c3_i136] = 0.0;
  }

  for (c3_i137 = 0; c3_i137 < 16; c3_i137++) {
    c3_l_y[c3_i137] = c3_b_y[c3_i137];
  }

  for (c3_i138 = 0; c3_i138 < 16; c3_i138++) {
    c3_q_b[c3_i138] = c3_b_b[c3_i138];
  }

  c3_b_eml_xgemm(chartInstance, c3_l_y, c3_q_b, c3_y);
  for (c3_i139 = 0; c3_i139 < 16; c3_i139++) {
    c3_b_b[c3_i139] = c3_A6[c3_i139];
  }

  c3_eml_scalar_eg(chartInstance);
  c3_eml_scalar_eg(chartInstance);
  for (c3_i140 = 0; c3_i140 < 16; c3_i140++) {
    c3_T6[c3_i140] = 0.0;
  }

  for (c3_i141 = 0; c3_i141 < 16; c3_i141++) {
    c3_T6[c3_i141] = 0.0;
  }

  for (c3_i142 = 0; c3_i142 < 16; c3_i142++) {
    c3_dv25[c3_i142] = c3_y[c3_i142];
  }

  for (c3_i143 = 0; c3_i143 < 16; c3_i143++) {
    c3_dv26[c3_i143] = c3_b_b[c3_i143];
  }

  for (c3_i144 = 0; c3_i144 < 16; c3_i144++) {
    c3_dv27[c3_i144] = c3_dv25[c3_i144];
  }

  for (c3_i145 = 0; c3_i145 < 16; c3_i145++) {
    c3_dv28[c3_i145] = c3_dv26[c3_i145];
  }

  c3_b_eml_xgemm(chartInstance, c3_dv27, c3_dv28, c3_T6);
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 27);
  for (c3_i146 = 0; c3_i146 < 16; c3_i146++) {
    c3_b_b[c3_i146] = c3_A1[c3_i146];
  }

  c3_eml_scalar_eg(chartInstance);
  c3_eml_scalar_eg(chartInstance);
  for (c3_i147 = 0; c3_i147 < 16; c3_i147++) {
    c3_y[c3_i147] = 0.0;
  }

  for (c3_i148 = 0; c3_i148 < 16; c3_i148++) {
    c3_h_a[c3_i148] = c3_b_a[c3_i148];
  }

  for (c3_i149 = 0; c3_i149 < 16; c3_i149++) {
    c3_r_b[c3_i149] = c3_b_b[c3_i149];
  }

  c3_b_eml_xgemm(chartInstance, c3_h_a, c3_r_b, c3_y);
  for (c3_i150 = 0; c3_i150 < 16; c3_i150++) {
    c3_b_b[c3_i150] = c3_A2[c3_i150];
  }

  c3_eml_scalar_eg(chartInstance);
  c3_eml_scalar_eg(chartInstance);
  for (c3_i151 = 0; c3_i151 < 16; c3_i151++) {
    c3_b_y[c3_i151] = 0.0;
  }

  for (c3_i152 = 0; c3_i152 < 16; c3_i152++) {
    c3_m_y[c3_i152] = c3_y[c3_i152];
  }

  for (c3_i153 = 0; c3_i153 < 16; c3_i153++) {
    c3_s_b[c3_i153] = c3_b_b[c3_i153];
  }

  c3_b_eml_xgemm(chartInstance, c3_m_y, c3_s_b, c3_b_y);
  for (c3_i154 = 0; c3_i154 < 16; c3_i154++) {
    c3_b_b[c3_i154] = c3_A3[c3_i154];
  }

  c3_eml_scalar_eg(chartInstance);
  c3_eml_scalar_eg(chartInstance);
  for (c3_i155 = 0; c3_i155 < 16; c3_i155++) {
    c3_y[c3_i155] = 0.0;
  }

  for (c3_i156 = 0; c3_i156 < 16; c3_i156++) {
    c3_n_y[c3_i156] = c3_b_y[c3_i156];
  }

  for (c3_i157 = 0; c3_i157 < 16; c3_i157++) {
    c3_t_b[c3_i157] = c3_b_b[c3_i157];
  }

  c3_b_eml_xgemm(chartInstance, c3_n_y, c3_t_b, c3_y);
  for (c3_i158 = 0; c3_i158 < 16; c3_i158++) {
    c3_b_b[c3_i158] = c3_A4[c3_i158];
  }

  c3_eml_scalar_eg(chartInstance);
  c3_eml_scalar_eg(chartInstance);
  for (c3_i159 = 0; c3_i159 < 16; c3_i159++) {
    c3_b_y[c3_i159] = 0.0;
  }

  for (c3_i160 = 0; c3_i160 < 16; c3_i160++) {
    c3_o_y[c3_i160] = c3_y[c3_i160];
  }

  for (c3_i161 = 0; c3_i161 < 16; c3_i161++) {
    c3_u_b[c3_i161] = c3_b_b[c3_i161];
  }

  c3_b_eml_xgemm(chartInstance, c3_o_y, c3_u_b, c3_b_y);
  for (c3_i162 = 0; c3_i162 < 16; c3_i162++) {
    c3_b_b[c3_i162] = c3_A5[c3_i162];
  }

  c3_eml_scalar_eg(chartInstance);
  c3_eml_scalar_eg(chartInstance);
  for (c3_i163 = 0; c3_i163 < 16; c3_i163++) {
    c3_y[c3_i163] = 0.0;
  }

  for (c3_i164 = 0; c3_i164 < 16; c3_i164++) {
    c3_p_y[c3_i164] = c3_b_y[c3_i164];
  }

  for (c3_i165 = 0; c3_i165 < 16; c3_i165++) {
    c3_v_b[c3_i165] = c3_b_b[c3_i165];
  }

  c3_b_eml_xgemm(chartInstance, c3_p_y, c3_v_b, c3_y);
  for (c3_i166 = 0; c3_i166 < 16; c3_i166++) {
    c3_b_b[c3_i166] = c3_A6[c3_i166];
  }

  c3_eml_scalar_eg(chartInstance);
  c3_eml_scalar_eg(chartInstance);
  for (c3_i167 = 0; c3_i167 < 16; c3_i167++) {
    c3_b_y[c3_i167] = 0.0;
  }

  for (c3_i168 = 0; c3_i168 < 16; c3_i168++) {
    c3_q_y[c3_i168] = c3_y[c3_i168];
  }

  for (c3_i169 = 0; c3_i169 < 16; c3_i169++) {
    c3_w_b[c3_i169] = c3_b_b[c3_i169];
  }

  c3_b_eml_xgemm(chartInstance, c3_q_y, c3_w_b, c3_b_y);
  for (c3_i170 = 0; c3_i170 < 16; c3_i170++) {
    c3_b_b[c3_i170] = c3_A7[c3_i170];
  }

  c3_eml_scalar_eg(chartInstance);
  c3_eml_scalar_eg(chartInstance);
  for (c3_i171 = 0; c3_i171 < 16; c3_i171++) {
    c3_y[c3_i171] = 0.0;
  }

  for (c3_i172 = 0; c3_i172 < 16; c3_i172++) {
    c3_r_y[c3_i172] = c3_b_y[c3_i172];
  }

  for (c3_i173 = 0; c3_i173 < 16; c3_i173++) {
    c3_x_b[c3_i173] = c3_b_b[c3_i173];
  }

  c3_b_eml_xgemm(chartInstance, c3_r_y, c3_x_b, c3_y);
  c3_eml_scalar_eg(chartInstance);
  c3_eml_scalar_eg(chartInstance);
  for (c3_i174 = 0; c3_i174 < 16; c3_i174++) {
    c3_Te[c3_i174] = 0.0;
  }

  for (c3_i175 = 0; c3_i175 < 16; c3_i175++) {
    c3_Te[c3_i175] = 0.0;
  }

  for (c3_i176 = 0; c3_i176 < 16; c3_i176++) {
    c3_dv29[c3_i176] = c3_y[c3_i176];
  }

  for (c3_i177 = 0; c3_i177 < 16; c3_i177++) {
    c3_dv30[c3_i177] = c3_b[c3_i177];
  }

  for (c3_i178 = 0; c3_i178 < 16; c3_i178++) {
    c3_dv31[c3_i178] = c3_dv29[c3_i178];
  }

  for (c3_i179 = 0; c3_i179 < 16; c3_i179++) {
    c3_dv32[c3_i179] = c3_dv30[c3_i179];
  }

  c3_b_eml_xgemm(chartInstance, c3_dv31, c3_dv32, c3_Te);
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 38);
  for (c3_i180 = 0; c3_i180 < 3; c3_i180++) {
    c3_z0[c3_i180] = c3_i_a[c3_i180];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 39);
  for (c3_i181 = 0; c3_i181 < 3; c3_i181++) {
    c3_z1[c3_i181] = c3_T1[c3_i181 + 8];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 40);
  for (c3_i182 = 0; c3_i182 < 3; c3_i182++) {
    c3_z2[c3_i182] = c3_T2[c3_i182 + 8];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 41);
  for (c3_i183 = 0; c3_i183 < 3; c3_i183++) {
    c3_z3[c3_i183] = c3_T3[c3_i183 + 8];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 42);
  for (c3_i184 = 0; c3_i184 < 3; c3_i184++) {
    c3_z4[c3_i184] = c3_T4[c3_i184 + 8];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 43);
  for (c3_i185 = 0; c3_i185 < 3; c3_i185++) {
    c3_z5[c3_i185] = c3_T5[c3_i185 + 8];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 44);
  for (c3_i186 = 0; c3_i186 < 3; c3_i186++) {
    c3_z6[c3_i186] = c3_T6[c3_i186 + 8];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 46);
  for (c3_i187 = 0; c3_i187 < 3; c3_i187++) {
    c3_p0[c3_i187] = 0.0;
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 47);
  for (c3_i188 = 0; c3_i188 < 3; c3_i188++) {
    c3_p1[c3_i188] = c3_T1[c3_i188 + 12];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 48);
  for (c3_i189 = 0; c3_i189 < 3; c3_i189++) {
    c3_p2[c3_i189] = c3_T2[c3_i189 + 12];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 49);
  for (c3_i190 = 0; c3_i190 < 3; c3_i190++) {
    c3_p3[c3_i190] = c3_T3[c3_i190 + 12];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 50);
  for (c3_i191 = 0; c3_i191 < 3; c3_i191++) {
    c3_p4[c3_i191] = c3_T4[c3_i191 + 12];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 51);
  for (c3_i192 = 0; c3_i192 < 3; c3_i192++) {
    c3_p5[c3_i192] = c3_T5[c3_i192 + 12];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 52);
  for (c3_i193 = 0; c3_i193 < 3; c3_i193++) {
    c3_p6[c3_i193] = c3_T6[c3_i193 + 12];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 53);
  for (c3_i194 = 0; c3_i194 < 3; c3_i194++) {
    c3_pe[c3_i194] = c3_Te[c3_i194 + 12];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 54);
  for (c3_i195 = 0; c3_i195 < 3; c3_i195++) {
    c3_y_b[c3_i195] = c3_pe[c3_i195] - c3_p0[c3_i195];
  }

  c3_c1 = 0.0 * c3_y_b[2] - c3_y_b[1];
  c3_c2 = c3_y_b[0] - 0.0 * c3_y_b[2];
  c3_c3 = 0.0 * c3_y_b[1] - 0.0 * c3_y_b[0];
  c3_dv33[0] = c3_c1;
  c3_dv33[1] = c3_c2;
  c3_dv33[2] = c3_c3;
  for (c3_i196 = 0; c3_i196 < 3; c3_i196++) {
    c3_j_a[c3_i196] = c3_z1[c3_i196];
  }

  for (c3_i197 = 0; c3_i197 < 3; c3_i197++) {
    c3_y_b[c3_i197] = c3_pe[c3_i197] - c3_p1[c3_i197];
  }

  c3_b_c1 = c3_j_a[1] * c3_y_b[2] - c3_j_a[2] * c3_y_b[1];
  c3_b_c2 = c3_j_a[2] * c3_y_b[0] - c3_j_a[0] * c3_y_b[2];
  c3_b_c3 = c3_j_a[0] * c3_y_b[1] - c3_j_a[1] * c3_y_b[0];
  c3_dv34[0] = c3_b_c1;
  c3_dv34[1] = c3_b_c2;
  c3_dv34[2] = c3_b_c3;
  for (c3_i198 = 0; c3_i198 < 3; c3_i198++) {
    c3_j_a[c3_i198] = c3_z2[c3_i198];
  }

  for (c3_i199 = 0; c3_i199 < 3; c3_i199++) {
    c3_y_b[c3_i199] = c3_pe[c3_i199] - c3_p2[c3_i199];
  }

  c3_c_c1 = c3_j_a[1] * c3_y_b[2] - c3_j_a[2] * c3_y_b[1];
  c3_c_c2 = c3_j_a[2] * c3_y_b[0] - c3_j_a[0] * c3_y_b[2];
  c3_c_c3 = c3_j_a[0] * c3_y_b[1] - c3_j_a[1] * c3_y_b[0];
  c3_dv35[0] = c3_c_c1;
  c3_dv35[1] = c3_c_c2;
  c3_dv35[2] = c3_c_c3;
  for (c3_i200 = 0; c3_i200 < 3; c3_i200++) {
    c3_j_a[c3_i200] = c3_z3[c3_i200];
  }

  for (c3_i201 = 0; c3_i201 < 3; c3_i201++) {
    c3_y_b[c3_i201] = c3_pe[c3_i201] - c3_p3[c3_i201];
  }

  c3_d_c1 = c3_j_a[1] * c3_y_b[2] - c3_j_a[2] * c3_y_b[1];
  c3_d_c2 = c3_j_a[2] * c3_y_b[0] - c3_j_a[0] * c3_y_b[2];
  c3_d_c3 = c3_j_a[0] * c3_y_b[1] - c3_j_a[1] * c3_y_b[0];
  c3_dv36[0] = c3_d_c1;
  c3_dv36[1] = c3_d_c2;
  c3_dv36[2] = c3_d_c3;
  for (c3_i202 = 0; c3_i202 < 3; c3_i202++) {
    c3_j_a[c3_i202] = c3_z4[c3_i202];
  }

  for (c3_i203 = 0; c3_i203 < 3; c3_i203++) {
    c3_y_b[c3_i203] = c3_pe[c3_i203] - c3_p4[c3_i203];
  }

  c3_e_c1 = c3_j_a[1] * c3_y_b[2] - c3_j_a[2] * c3_y_b[1];
  c3_e_c2 = c3_j_a[2] * c3_y_b[0] - c3_j_a[0] * c3_y_b[2];
  c3_e_c3 = c3_j_a[0] * c3_y_b[1] - c3_j_a[1] * c3_y_b[0];
  c3_dv37[0] = c3_e_c1;
  c3_dv37[1] = c3_e_c2;
  c3_dv37[2] = c3_e_c3;
  for (c3_i204 = 0; c3_i204 < 3; c3_i204++) {
    c3_j_a[c3_i204] = c3_z5[c3_i204];
  }

  for (c3_i205 = 0; c3_i205 < 3; c3_i205++) {
    c3_y_b[c3_i205] = c3_pe[c3_i205] - c3_p5[c3_i205];
  }

  c3_f_c1 = c3_j_a[1] * c3_y_b[2] - c3_j_a[2] * c3_y_b[1];
  c3_f_c2 = c3_j_a[2] * c3_y_b[0] - c3_j_a[0] * c3_y_b[2];
  c3_f_c3 = c3_j_a[0] * c3_y_b[1] - c3_j_a[1] * c3_y_b[0];
  c3_dv38[0] = c3_f_c1;
  c3_dv38[1] = c3_f_c2;
  c3_dv38[2] = c3_f_c3;
  for (c3_i206 = 0; c3_i206 < 3; c3_i206++) {
    c3_j_a[c3_i206] = c3_z6[c3_i206];
  }

  for (c3_i207 = 0; c3_i207 < 3; c3_i207++) {
    c3_y_b[c3_i207] = c3_pe[c3_i207] - c3_p6[c3_i207];
  }

  c3_g_c1 = c3_j_a[1] * c3_y_b[2] - c3_j_a[2] * c3_y_b[1];
  c3_g_c2 = c3_j_a[2] * c3_y_b[0] - c3_j_a[0] * c3_y_b[2];
  c3_g_c3 = c3_j_a[0] * c3_y_b[1] - c3_j_a[1] * c3_y_b[0];
  c3_y_b[0] = c3_g_c1;
  c3_y_b[1] = c3_g_c2;
  c3_y_b[2] = c3_g_c3;
  for (c3_i208 = 0; c3_i208 < 3; c3_i208++) {
    c3_je[c3_i208] = c3_dv33[c3_i208];
  }

  for (c3_i209 = 0; c3_i209 < 3; c3_i209++) {
    c3_je[c3_i209 + 6] = c3_dv34[c3_i209];
  }

  for (c3_i210 = 0; c3_i210 < 3; c3_i210++) {
    c3_je[c3_i210 + 12] = c3_dv35[c3_i210];
  }

  for (c3_i211 = 0; c3_i211 < 3; c3_i211++) {
    c3_je[c3_i211 + 18] = c3_dv36[c3_i211];
  }

  for (c3_i212 = 0; c3_i212 < 3; c3_i212++) {
    c3_je[c3_i212 + 24] = c3_dv37[c3_i212];
  }

  for (c3_i213 = 0; c3_i213 < 3; c3_i213++) {
    c3_je[c3_i213 + 30] = c3_dv38[c3_i213];
  }

  for (c3_i214 = 0; c3_i214 < 3; c3_i214++) {
    c3_je[c3_i214 + 36] = c3_y_b[c3_i214];
  }

  for (c3_i215 = 0; c3_i215 < 3; c3_i215++) {
    c3_je[c3_i215 + 3] = c3_z0[c3_i215];
  }

  for (c3_i216 = 0; c3_i216 < 3; c3_i216++) {
    c3_je[c3_i216 + 9] = c3_z1[c3_i216];
  }

  for (c3_i217 = 0; c3_i217 < 3; c3_i217++) {
    c3_je[c3_i217 + 15] = c3_z2[c3_i217];
  }

  for (c3_i218 = 0; c3_i218 < 3; c3_i218++) {
    c3_je[c3_i218 + 21] = c3_z3[c3_i218];
  }

  for (c3_i219 = 0; c3_i219 < 3; c3_i219++) {
    c3_je[c3_i219 + 27] = c3_z4[c3_i219];
  }

  for (c3_i220 = 0; c3_i220 < 3; c3_i220++) {
    c3_je[c3_i220 + 33] = c3_z5[c3_i220];
  }

  for (c3_i221 = 0; c3_i221 < 3; c3_i221++) {
    c3_je[c3_i221 + 39] = c3_z6[c3_i221];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, -54);
  _SFD_SYMBOL_SCOPE_POP();
}

static void init_script_number_translation(uint32_T c3_machineNumber, uint32_T
  c3_chartNumber, uint32_T c3_instanceNumber)
{
  (void)c3_machineNumber;
  _SFD_SCRIPT_TRANSLATION(c3_chartNumber, c3_instanceNumber, 0U,
    sf_debug_get_script_id(
    "C:\\Users\\faezeh\\Desktop\\GitHub\\Motion Planning for KUKA LBR iiwa\\KinematicSym\\jacobin.m"));
}

static const mxArray *c3_sf_marshallOut(void *chartInstanceVoid, void *c3_inData)
{
  const mxArray *c3_mxArrayOutData = NULL;
  int32_T c3_i222;
  real_T c3_b_inData[3];
  int32_T c3_i223;
  real_T c3_u[3];
  const mxArray *c3_y = NULL;
  SFc3_simlwrkuka_kinematicsInstanceStruct *chartInstance;
  chartInstance = (SFc3_simlwrkuka_kinematicsInstanceStruct *)chartInstanceVoid;
  c3_mxArrayOutData = NULL;
  for (c3_i222 = 0; c3_i222 < 3; c3_i222++) {
    c3_b_inData[c3_i222] = (*(real_T (*)[3])c3_inData)[c3_i222];
  }

  for (c3_i223 = 0; c3_i223 < 3; c3_i223++) {
    c3_u[c3_i223] = c3_b_inData[c3_i223];
  }

  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", c3_u, 0, 0U, 1U, 0U, 1, 3), false);
  sf_mex_assign(&c3_mxArrayOutData, c3_y, false);
  return c3_mxArrayOutData;
}

static void c3_emlrt_marshallIn(SFc3_simlwrkuka_kinematicsInstanceStruct
  *chartInstance, const mxArray *c3_xe, const char_T *c3_identifier, real_T
  c3_y[3])
{
  emlrtMsgIdentifier c3_thisId;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_xe), &c3_thisId, c3_y);
  sf_mex_destroy(&c3_xe);
}

static void c3_b_emlrt_marshallIn(SFc3_simlwrkuka_kinematicsInstanceStruct
  *chartInstance, const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId,
  real_T c3_y[3])
{
  real_T c3_dv39[3];
  int32_T c3_i224;
  (void)chartInstance;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), c3_dv39, 1, 0, 0U, 1, 0U, 1, 3);
  for (c3_i224 = 0; c3_i224 < 3; c3_i224++) {
    c3_y[c3_i224] = c3_dv39[c3_i224];
  }

  sf_mex_destroy(&c3_u);
}

static void c3_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData)
{
  const mxArray *c3_xe;
  const char_T *c3_identifier;
  emlrtMsgIdentifier c3_thisId;
  real_T c3_y[3];
  int32_T c3_i225;
  SFc3_simlwrkuka_kinematicsInstanceStruct *chartInstance;
  chartInstance = (SFc3_simlwrkuka_kinematicsInstanceStruct *)chartInstanceVoid;
  c3_xe = sf_mex_dup(c3_mxArrayInData);
  c3_identifier = c3_varName;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_xe), &c3_thisId, c3_y);
  sf_mex_destroy(&c3_xe);
  for (c3_i225 = 0; c3_i225 < 3; c3_i225++) {
    (*(real_T (*)[3])c3_outData)[c3_i225] = c3_y[c3_i225];
  }

  sf_mex_destroy(&c3_mxArrayInData);
}

static const mxArray *c3_b_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData)
{
  const mxArray *c3_mxArrayOutData = NULL;
  int32_T c3_i226;
  real_T c3_b_inData[7];
  int32_T c3_i227;
  real_T c3_u[7];
  const mxArray *c3_y = NULL;
  SFc3_simlwrkuka_kinematicsInstanceStruct *chartInstance;
  chartInstance = (SFc3_simlwrkuka_kinematicsInstanceStruct *)chartInstanceVoid;
  c3_mxArrayOutData = NULL;
  for (c3_i226 = 0; c3_i226 < 7; c3_i226++) {
    c3_b_inData[c3_i226] = (*(real_T (*)[7])c3_inData)[c3_i226];
  }

  for (c3_i227 = 0; c3_i227 < 7; c3_i227++) {
    c3_u[c3_i227] = c3_b_inData[c3_i227];
  }

  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", c3_u, 0, 0U, 1U, 0U, 1, 7), false);
  sf_mex_assign(&c3_mxArrayOutData, c3_y, false);
  return c3_mxArrayOutData;
}

static void c3_c_emlrt_marshallIn(SFc3_simlwrkuka_kinematicsInstanceStruct
  *chartInstance, const mxArray *c3_x, const char_T *c3_identifier, real_T c3_y
  [7])
{
  emlrtMsgIdentifier c3_thisId;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_x), &c3_thisId, c3_y);
  sf_mex_destroy(&c3_x);
}

static void c3_d_emlrt_marshallIn(SFc3_simlwrkuka_kinematicsInstanceStruct
  *chartInstance, const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId,
  real_T c3_y[7])
{
  real_T c3_dv40[7];
  int32_T c3_i228;
  (void)chartInstance;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), c3_dv40, 1, 0, 0U, 1, 0U, 1, 7);
  for (c3_i228 = 0; c3_i228 < 7; c3_i228++) {
    c3_y[c3_i228] = c3_dv40[c3_i228];
  }

  sf_mex_destroy(&c3_u);
}

static void c3_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData)
{
  const mxArray *c3_x;
  const char_T *c3_identifier;
  emlrtMsgIdentifier c3_thisId;
  real_T c3_y[7];
  int32_T c3_i229;
  SFc3_simlwrkuka_kinematicsInstanceStruct *chartInstance;
  chartInstance = (SFc3_simlwrkuka_kinematicsInstanceStruct *)chartInstanceVoid;
  c3_x = sf_mex_dup(c3_mxArrayInData);
  c3_identifier = c3_varName;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_x), &c3_thisId, c3_y);
  sf_mex_destroy(&c3_x);
  for (c3_i229 = 0; c3_i229 < 7; c3_i229++) {
    (*(real_T (*)[7])c3_outData)[c3_i229] = c3_y[c3_i229];
  }

  sf_mex_destroy(&c3_mxArrayInData);
}

static const mxArray *c3_c_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData)
{
  const mxArray *c3_mxArrayOutData = NULL;
  real_T c3_u;
  const mxArray *c3_y = NULL;
  SFc3_simlwrkuka_kinematicsInstanceStruct *chartInstance;
  chartInstance = (SFc3_simlwrkuka_kinematicsInstanceStruct *)chartInstanceVoid;
  c3_mxArrayOutData = NULL;
  c3_u = *(real_T *)c3_inData;
  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", &c3_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c3_mxArrayOutData, c3_y, false);
  return c3_mxArrayOutData;
}

static real_T c3_e_emlrt_marshallIn(SFc3_simlwrkuka_kinematicsInstanceStruct
  *chartInstance, const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId)
{
  real_T c3_y;
  real_T c3_d0;
  (void)chartInstance;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), &c3_d0, 1, 0, 0U, 0, 0U, 0);
  c3_y = c3_d0;
  sf_mex_destroy(&c3_u);
  return c3_y;
}

static void c3_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData)
{
  const mxArray *c3_nargout;
  const char_T *c3_identifier;
  emlrtMsgIdentifier c3_thisId;
  real_T c3_y;
  SFc3_simlwrkuka_kinematicsInstanceStruct *chartInstance;
  chartInstance = (SFc3_simlwrkuka_kinematicsInstanceStruct *)chartInstanceVoid;
  c3_nargout = sf_mex_dup(c3_mxArrayInData);
  c3_identifier = c3_varName;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_y = c3_e_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_nargout), &c3_thisId);
  sf_mex_destroy(&c3_nargout);
  *(real_T *)c3_outData = c3_y;
  sf_mex_destroy(&c3_mxArrayInData);
}

static const mxArray *c3_d_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData)
{
  const mxArray *c3_mxArrayOutData = NULL;
  int32_T c3_i230;
  real_T c3_b_inData[4];
  int32_T c3_i231;
  real_T c3_u[4];
  const mxArray *c3_y = NULL;
  SFc3_simlwrkuka_kinematicsInstanceStruct *chartInstance;
  chartInstance = (SFc3_simlwrkuka_kinematicsInstanceStruct *)chartInstanceVoid;
  c3_mxArrayOutData = NULL;
  for (c3_i230 = 0; c3_i230 < 4; c3_i230++) {
    c3_b_inData[c3_i230] = (*(real_T (*)[4])c3_inData)[c3_i230];
  }

  for (c3_i231 = 0; c3_i231 < 4; c3_i231++) {
    c3_u[c3_i231] = c3_b_inData[c3_i231];
  }

  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", c3_u, 0, 0U, 1U, 0U, 1, 4), false);
  sf_mex_assign(&c3_mxArrayOutData, c3_y, false);
  return c3_mxArrayOutData;
}

static void c3_f_emlrt_marshallIn(SFc3_simlwrkuka_kinematicsInstanceStruct
  *chartInstance, const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId,
  real_T c3_y[4])
{
  real_T c3_dv41[4];
  int32_T c3_i232;
  (void)chartInstance;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), c3_dv41, 1, 0, 0U, 1, 0U, 1, 4);
  for (c3_i232 = 0; c3_i232 < 4; c3_i232++) {
    c3_y[c3_i232] = c3_dv41[c3_i232];
  }

  sf_mex_destroy(&c3_u);
}

static void c3_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData)
{
  const mxArray *c3_Q;
  const char_T *c3_identifier;
  emlrtMsgIdentifier c3_thisId;
  real_T c3_y[4];
  int32_T c3_i233;
  SFc3_simlwrkuka_kinematicsInstanceStruct *chartInstance;
  chartInstance = (SFc3_simlwrkuka_kinematicsInstanceStruct *)chartInstanceVoid;
  c3_Q = sf_mex_dup(c3_mxArrayInData);
  c3_identifier = c3_varName;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_Q), &c3_thisId, c3_y);
  sf_mex_destroy(&c3_Q);
  for (c3_i233 = 0; c3_i233 < 4; c3_i233++) {
    (*(real_T (*)[4])c3_outData)[c3_i233] = c3_y[c3_i233];
  }

  sf_mex_destroy(&c3_mxArrayInData);
}

static const mxArray *c3_e_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData)
{
  const mxArray *c3_mxArrayOutData = NULL;
  int32_T c3_i234;
  int32_T c3_i235;
  int32_T c3_i236;
  real_T c3_b_inData[9];
  int32_T c3_i237;
  int32_T c3_i238;
  int32_T c3_i239;
  real_T c3_u[9];
  const mxArray *c3_y = NULL;
  SFc3_simlwrkuka_kinematicsInstanceStruct *chartInstance;
  chartInstance = (SFc3_simlwrkuka_kinematicsInstanceStruct *)chartInstanceVoid;
  c3_mxArrayOutData = NULL;
  c3_i234 = 0;
  for (c3_i235 = 0; c3_i235 < 3; c3_i235++) {
    for (c3_i236 = 0; c3_i236 < 3; c3_i236++) {
      c3_b_inData[c3_i236 + c3_i234] = (*(real_T (*)[9])c3_inData)[c3_i236 +
        c3_i234];
    }

    c3_i234 += 3;
  }

  c3_i237 = 0;
  for (c3_i238 = 0; c3_i238 < 3; c3_i238++) {
    for (c3_i239 = 0; c3_i239 < 3; c3_i239++) {
      c3_u[c3_i239 + c3_i237] = c3_b_inData[c3_i239 + c3_i237];
    }

    c3_i237 += 3;
  }

  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", c3_u, 0, 0U, 1U, 0U, 2, 3, 3), false);
  sf_mex_assign(&c3_mxArrayOutData, c3_y, false);
  return c3_mxArrayOutData;
}

static void c3_g_emlrt_marshallIn(SFc3_simlwrkuka_kinematicsInstanceStruct
  *chartInstance, const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId,
  real_T c3_y[9])
{
  real_T c3_dv42[9];
  int32_T c3_i240;
  (void)chartInstance;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), c3_dv42, 1, 0, 0U, 1, 0U, 2, 3, 3);
  for (c3_i240 = 0; c3_i240 < 9; c3_i240++) {
    c3_y[c3_i240] = c3_dv42[c3_i240];
  }

  sf_mex_destroy(&c3_u);
}

static void c3_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData)
{
  const mxArray *c3_R7;
  const char_T *c3_identifier;
  emlrtMsgIdentifier c3_thisId;
  real_T c3_y[9];
  int32_T c3_i241;
  int32_T c3_i242;
  int32_T c3_i243;
  SFc3_simlwrkuka_kinematicsInstanceStruct *chartInstance;
  chartInstance = (SFc3_simlwrkuka_kinematicsInstanceStruct *)chartInstanceVoid;
  c3_R7 = sf_mex_dup(c3_mxArrayInData);
  c3_identifier = c3_varName;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_g_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_R7), &c3_thisId, c3_y);
  sf_mex_destroy(&c3_R7);
  c3_i241 = 0;
  for (c3_i242 = 0; c3_i242 < 3; c3_i242++) {
    for (c3_i243 = 0; c3_i243 < 3; c3_i243++) {
      (*(real_T (*)[9])c3_outData)[c3_i243 + c3_i241] = c3_y[c3_i243 + c3_i241];
    }

    c3_i241 += 3;
  }

  sf_mex_destroy(&c3_mxArrayInData);
}

static const mxArray *c3_f_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData)
{
  const mxArray *c3_mxArrayOutData = NULL;
  int32_T c3_i244;
  int32_T c3_i245;
  int32_T c3_i246;
  real_T c3_b_inData[16];
  int32_T c3_i247;
  int32_T c3_i248;
  int32_T c3_i249;
  real_T c3_u[16];
  const mxArray *c3_y = NULL;
  SFc3_simlwrkuka_kinematicsInstanceStruct *chartInstance;
  chartInstance = (SFc3_simlwrkuka_kinematicsInstanceStruct *)chartInstanceVoid;
  c3_mxArrayOutData = NULL;
  c3_i244 = 0;
  for (c3_i245 = 0; c3_i245 < 4; c3_i245++) {
    for (c3_i246 = 0; c3_i246 < 4; c3_i246++) {
      c3_b_inData[c3_i246 + c3_i244] = (*(real_T (*)[16])c3_inData)[c3_i246 +
        c3_i244];
    }

    c3_i244 += 4;
  }

  c3_i247 = 0;
  for (c3_i248 = 0; c3_i248 < 4; c3_i248++) {
    for (c3_i249 = 0; c3_i249 < 4; c3_i249++) {
      c3_u[c3_i249 + c3_i247] = c3_b_inData[c3_i249 + c3_i247];
    }

    c3_i247 += 4;
  }

  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", c3_u, 0, 0U, 1U, 0U, 2, 4, 4), false);
  sf_mex_assign(&c3_mxArrayOutData, c3_y, false);
  return c3_mxArrayOutData;
}

static void c3_h_emlrt_marshallIn(SFc3_simlwrkuka_kinematicsInstanceStruct
  *chartInstance, const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId,
  real_T c3_y[16])
{
  real_T c3_dv43[16];
  int32_T c3_i250;
  (void)chartInstance;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), c3_dv43, 1, 0, 0U, 1, 0U, 2, 4, 4);
  for (c3_i250 = 0; c3_i250 < 16; c3_i250++) {
    c3_y[c3_i250] = c3_dv43[c3_i250];
  }

  sf_mex_destroy(&c3_u);
}

static void c3_f_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData)
{
  const mxArray *c3_Te;
  const char_T *c3_identifier;
  emlrtMsgIdentifier c3_thisId;
  real_T c3_y[16];
  int32_T c3_i251;
  int32_T c3_i252;
  int32_T c3_i253;
  SFc3_simlwrkuka_kinematicsInstanceStruct *chartInstance;
  chartInstance = (SFc3_simlwrkuka_kinematicsInstanceStruct *)chartInstanceVoid;
  c3_Te = sf_mex_dup(c3_mxArrayInData);
  c3_identifier = c3_varName;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_h_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_Te), &c3_thisId, c3_y);
  sf_mex_destroy(&c3_Te);
  c3_i251 = 0;
  for (c3_i252 = 0; c3_i252 < 4; c3_i252++) {
    for (c3_i253 = 0; c3_i253 < 4; c3_i253++) {
      (*(real_T (*)[16])c3_outData)[c3_i253 + c3_i251] = c3_y[c3_i253 + c3_i251];
    }

    c3_i251 += 4;
  }

  sf_mex_destroy(&c3_mxArrayInData);
}

static const mxArray *c3_g_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData)
{
  const mxArray *c3_mxArrayOutData = NULL;
  int32_T c3_i254;
  int32_T c3_i255;
  int32_T c3_i256;
  real_T c3_b_inData[42];
  int32_T c3_i257;
  int32_T c3_i258;
  int32_T c3_i259;
  real_T c3_u[42];
  const mxArray *c3_y = NULL;
  SFc3_simlwrkuka_kinematicsInstanceStruct *chartInstance;
  chartInstance = (SFc3_simlwrkuka_kinematicsInstanceStruct *)chartInstanceVoid;
  c3_mxArrayOutData = NULL;
  c3_i254 = 0;
  for (c3_i255 = 0; c3_i255 < 7; c3_i255++) {
    for (c3_i256 = 0; c3_i256 < 6; c3_i256++) {
      c3_b_inData[c3_i256 + c3_i254] = (*(real_T (*)[42])c3_inData)[c3_i256 +
        c3_i254];
    }

    c3_i254 += 6;
  }

  c3_i257 = 0;
  for (c3_i258 = 0; c3_i258 < 7; c3_i258++) {
    for (c3_i259 = 0; c3_i259 < 6; c3_i259++) {
      c3_u[c3_i259 + c3_i257] = c3_b_inData[c3_i259 + c3_i257];
    }

    c3_i257 += 6;
  }

  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", c3_u, 0, 0U, 1U, 0U, 2, 6, 7), false);
  sf_mex_assign(&c3_mxArrayOutData, c3_y, false);
  return c3_mxArrayOutData;
}

static void c3_i_emlrt_marshallIn(SFc3_simlwrkuka_kinematicsInstanceStruct
  *chartInstance, const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId,
  real_T c3_y[42])
{
  real_T c3_dv44[42];
  int32_T c3_i260;
  (void)chartInstance;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), c3_dv44, 1, 0, 0U, 1, 0U, 2, 6, 7);
  for (c3_i260 = 0; c3_i260 < 42; c3_i260++) {
    c3_y[c3_i260] = c3_dv44[c3_i260];
  }

  sf_mex_destroy(&c3_u);
}

static void c3_g_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData)
{
  const mxArray *c3_je;
  const char_T *c3_identifier;
  emlrtMsgIdentifier c3_thisId;
  real_T c3_y[42];
  int32_T c3_i261;
  int32_T c3_i262;
  int32_T c3_i263;
  SFc3_simlwrkuka_kinematicsInstanceStruct *chartInstance;
  chartInstance = (SFc3_simlwrkuka_kinematicsInstanceStruct *)chartInstanceVoid;
  c3_je = sf_mex_dup(c3_mxArrayInData);
  c3_identifier = c3_varName;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_i_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_je), &c3_thisId, c3_y);
  sf_mex_destroy(&c3_je);
  c3_i261 = 0;
  for (c3_i262 = 0; c3_i262 < 7; c3_i262++) {
    for (c3_i263 = 0; c3_i263 < 6; c3_i263++) {
      (*(real_T (*)[42])c3_outData)[c3_i263 + c3_i261] = c3_y[c3_i263 + c3_i261];
    }

    c3_i261 += 6;
  }

  sf_mex_destroy(&c3_mxArrayInData);
}

const mxArray *sf_c3_simlwrkuka_kinematics_get_eml_resolved_functions_info(void)
{
  const mxArray *c3_nameCaptureInfo = NULL;
  c3_nameCaptureInfo = NULL;
  sf_mex_assign(&c3_nameCaptureInfo, sf_mex_createstruct("structure", 2, 33, 1),
                false);
  c3_info_helper(&c3_nameCaptureInfo);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c3_nameCaptureInfo);
  return c3_nameCaptureInfo;
}

static void c3_info_helper(const mxArray **c3_info)
{
  const mxArray *c3_rhs0 = NULL;
  const mxArray *c3_lhs0 = NULL;
  const mxArray *c3_rhs1 = NULL;
  const mxArray *c3_lhs1 = NULL;
  const mxArray *c3_rhs2 = NULL;
  const mxArray *c3_lhs2 = NULL;
  const mxArray *c3_rhs3 = NULL;
  const mxArray *c3_lhs3 = NULL;
  const mxArray *c3_rhs4 = NULL;
  const mxArray *c3_lhs4 = NULL;
  const mxArray *c3_rhs5 = NULL;
  const mxArray *c3_lhs5 = NULL;
  const mxArray *c3_rhs6 = NULL;
  const mxArray *c3_lhs6 = NULL;
  const mxArray *c3_rhs7 = NULL;
  const mxArray *c3_lhs7 = NULL;
  const mxArray *c3_rhs8 = NULL;
  const mxArray *c3_lhs8 = NULL;
  const mxArray *c3_rhs9 = NULL;
  const mxArray *c3_lhs9 = NULL;
  const mxArray *c3_rhs10 = NULL;
  const mxArray *c3_lhs10 = NULL;
  const mxArray *c3_rhs11 = NULL;
  const mxArray *c3_lhs11 = NULL;
  const mxArray *c3_rhs12 = NULL;
  const mxArray *c3_lhs12 = NULL;
  const mxArray *c3_rhs13 = NULL;
  const mxArray *c3_lhs13 = NULL;
  const mxArray *c3_rhs14 = NULL;
  const mxArray *c3_lhs14 = NULL;
  const mxArray *c3_rhs15 = NULL;
  const mxArray *c3_lhs15 = NULL;
  const mxArray *c3_rhs16 = NULL;
  const mxArray *c3_lhs16 = NULL;
  const mxArray *c3_rhs17 = NULL;
  const mxArray *c3_lhs17 = NULL;
  const mxArray *c3_rhs18 = NULL;
  const mxArray *c3_lhs18 = NULL;
  const mxArray *c3_rhs19 = NULL;
  const mxArray *c3_lhs19 = NULL;
  const mxArray *c3_rhs20 = NULL;
  const mxArray *c3_lhs20 = NULL;
  const mxArray *c3_rhs21 = NULL;
  const mxArray *c3_lhs21 = NULL;
  const mxArray *c3_rhs22 = NULL;
  const mxArray *c3_lhs22 = NULL;
  const mxArray *c3_rhs23 = NULL;
  const mxArray *c3_lhs23 = NULL;
  const mxArray *c3_rhs24 = NULL;
  const mxArray *c3_lhs24 = NULL;
  const mxArray *c3_rhs25 = NULL;
  const mxArray *c3_lhs25 = NULL;
  const mxArray *c3_rhs26 = NULL;
  const mxArray *c3_lhs26 = NULL;
  const mxArray *c3_rhs27 = NULL;
  const mxArray *c3_lhs27 = NULL;
  const mxArray *c3_rhs28 = NULL;
  const mxArray *c3_lhs28 = NULL;
  const mxArray *c3_rhs29 = NULL;
  const mxArray *c3_lhs29 = NULL;
  const mxArray *c3_rhs30 = NULL;
  const mxArray *c3_lhs30 = NULL;
  const mxArray *c3_rhs31 = NULL;
  const mxArray *c3_lhs31 = NULL;
  const mxArray *c3_rhs32 = NULL;
  const mxArray *c3_lhs32 = NULL;
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "context", "context", 0);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("jacobin"), "name", "name", 0);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 0);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[E]C:/Users/faezeh/Desktop/GitHub/Motion Planning for KUKA LBR iiwa/KinematicSym/jacobin.m"),
                  "resolved", "resolved", 0);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1500209896U), "fileTimeLo",
                  "fileTimeLo", 0);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 0);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 0);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 0);
  sf_mex_assign(&c3_rhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs0), "rhs", "rhs", 0);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs0), "lhs", "lhs", 0);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[E]C:/Users/faezeh/Desktop/GitHub/Motion Planning for KUKA LBR iiwa/KinematicSym/jacobin.m"),
                  "context", "context", 1);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("mrdivide"), "name", "name", 1);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 1);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "resolved",
                  "resolved", 1);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1388463696U), "fileTimeLo",
                  "fileTimeLo", 1);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 1);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1370017086U), "mFileTimeLo",
                  "mFileTimeLo", 1);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 1);
  sf_mex_assign(&c3_rhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs1), "rhs", "rhs", 1);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs1), "lhs", "lhs", 1);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 2);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.assert"),
                  "name", "name", 2);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 2);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/assert.m"),
                  "resolved", "resolved", 2);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 2);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 2);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 2);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 2);
  sf_mex_assign(&c3_rhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs2), "rhs", "rhs", 2);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs2), "lhs", "lhs", 2);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 3);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("rdivide"), "name", "name", 3);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 3);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "resolved",
                  "resolved", 3);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363717480U), "fileTimeLo",
                  "fileTimeLo", 3);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 3);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 3);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 3);
  sf_mex_assign(&c3_rhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs3), "rhs", "rhs", 3);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs3), "lhs", "lhs", 3);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 4);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 4);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 4);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 4);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 4);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 4);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 4);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 4);
  sf_mex_assign(&c3_rhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs4), "rhs", "rhs", 4);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs4), "lhs", "lhs", 4);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 5);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_scalexp_compatible"),
                  "name", "name", 5);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 5);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_compatible.m"),
                  "resolved", "resolved", 5);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1286825996U), "fileTimeLo",
                  "fileTimeLo", 5);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 5);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 5);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 5);
  sf_mex_assign(&c3_rhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs5), "rhs", "rhs", 5);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs5), "lhs", "lhs", 5);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 6);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_div"), "name", "name", 6);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 6);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 6);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 6);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 6);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 6);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 6);
  sf_mex_assign(&c3_rhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs6), "rhs", "rhs", 6);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs6), "lhs", "lhs", 6);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "context",
                  "context", 7);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.div"), "name",
                  "name", 7);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 7);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/div.p"), "resolved",
                  "resolved", 7);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 7);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 7);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 7);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 7);
  sf_mex_assign(&c3_rhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs7), "rhs", "rhs", 7);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs7), "lhs", "lhs", 7);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[E]C:/Users/faezeh/Desktop/GitHub/Motion Planning for KUKA LBR iiwa/KinematicSym/jacobin.m"),
                  "context", "context", 8);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("cos"), "name", "name", 8);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 8);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "resolved",
                  "resolved", 8);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1343837572U), "fileTimeLo",
                  "fileTimeLo", 8);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 8);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 8);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 8);
  sf_mex_assign(&c3_rhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs8), "rhs", "rhs", 8);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs8), "lhs", "lhs", 8);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "context",
                  "context", 9);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_scalar_cos"), "name",
                  "name", 9);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 9);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_cos.m"),
                  "resolved", "resolved", 9);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1286825922U), "fileTimeLo",
                  "fileTimeLo", 9);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 9);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 9);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 9);
  sf_mex_assign(&c3_rhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs9), "rhs", "rhs", 9);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs9), "lhs", "lhs", 9);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[E]C:/Users/faezeh/Desktop/GitHub/Motion Planning for KUKA LBR iiwa/KinematicSym/jacobin.m"),
                  "context", "context", 10);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("sin"), "name", "name", 10);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 10);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "resolved",
                  "resolved", 10);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1343837586U), "fileTimeLo",
                  "fileTimeLo", 10);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 10);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 10);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 10);
  sf_mex_assign(&c3_rhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs10), "rhs", "rhs",
                  10);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs10), "lhs", "lhs",
                  10);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "context",
                  "context", 11);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_scalar_sin"), "name",
                  "name", 11);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 11);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sin.m"),
                  "resolved", "resolved", 11);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1286825936U), "fileTimeLo",
                  "fileTimeLo", 11);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 11);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 11);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 11);
  sf_mex_assign(&c3_rhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs11), "rhs", "rhs",
                  11);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs11), "lhs", "lhs",
                  11);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[E]C:/Users/faezeh/Desktop/GitHub/Motion Planning for KUKA LBR iiwa/KinematicSym/jacobin.m"),
                  "context", "context", 12);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 12);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 12);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 12);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1383880894U), "fileTimeLo",
                  "fileTimeLo", 12);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 12);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 12);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 12);
  sf_mex_assign(&c3_rhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs12), "rhs", "rhs",
                  12);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs12), "lhs", "lhs",
                  12);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m!common_checks"),
                  "context", "context", 13);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 13);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 13);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 13);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 13);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 13);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 13);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 13);
  sf_mex_assign(&c3_rhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs13), "rhs", "rhs",
                  13);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs13), "lhs", "lhs",
                  13);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 14);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 14);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 14);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 14);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 14);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 14);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 14);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 14);
  sf_mex_assign(&c3_rhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs14), "rhs", "rhs",
                  14);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs14), "lhs", "lhs",
                  14);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 15);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 15);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 15);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 15);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 15);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 15);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 15);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 15);
  sf_mex_assign(&c3_rhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs15), "rhs", "rhs",
                  15);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs15), "lhs", "lhs",
                  15);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "context",
                  "context", 16);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 16);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 16);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 16);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 16);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 16);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 16);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 16);
  sf_mex_assign(&c3_rhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs16), "rhs", "rhs",
                  16);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs16), "lhs", "lhs",
                  16);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 17);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_xgemm"), "name", "name",
                  17);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 17);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"),
                  "resolved", "resolved", 17);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1375987890U), "fileTimeLo",
                  "fileTimeLo", 17);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 17);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 17);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 17);
  sf_mex_assign(&c3_rhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs17), "rhs", "rhs",
                  17);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs17), "lhs", "lhs",
                  17);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 18);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 18);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 18);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 18);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 18);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 18);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 18);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 18);
  sf_mex_assign(&c3_rhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs18), "rhs", "rhs",
                  18);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs18), "lhs", "lhs",
                  18);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 19);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.blas.xgemm"),
                  "name", "name", 19);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 19);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "resolved", "resolved", 19);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 19);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 19);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 19);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 19);
  sf_mex_assign(&c3_rhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs19), "rhs", "rhs",
                  19);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs19), "lhs", "lhs",
                  19);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 20);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 20);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 20);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 20);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 20);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 20);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 20);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 20);
  sf_mex_assign(&c3_rhs20, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs20, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs20), "rhs", "rhs",
                  20);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs20), "lhs", "lhs",
                  20);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p!below_threshold"),
                  "context", "context", 21);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.blas.threshold"),
                  "name", "name", 21);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 21);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 21);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 21);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 21);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 21);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 21);
  sf_mex_assign(&c3_rhs21, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs21, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs21), "rhs", "rhs",
                  21);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs21), "lhs", "lhs",
                  21);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "context", "context", 22);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 22);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 22);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 22);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1381857500U), "fileTimeLo",
                  "fileTimeLo", 22);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 22);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 22);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 22);
  sf_mex_assign(&c3_rhs22, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs22, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs22), "rhs", "rhs",
                  22);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs22), "lhs", "lhs",
                  22);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 23);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 23);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 23);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 23);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 23);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 23);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 23);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 23);
  sf_mex_assign(&c3_rhs23, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs23, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs23), "rhs", "rhs",
                  23);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs23), "lhs", "lhs",
                  23);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 24);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.refblas.xgemm"),
                  "name", "name", 24);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 24);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgemm.p"),
                  "resolved", "resolved", 24);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 24);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 24);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 24);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 24);
  sf_mex_assign(&c3_rhs24, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs24, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs24), "rhs", "rhs",
                  24);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs24), "lhs", "lhs",
                  24);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[E]C:/Users/faezeh/Desktop/GitHub/Motion Planning for KUKA LBR iiwa/KinematicSym/jacobin.m"),
                  "context", "context", 25);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("cross"), "name", "name", 25);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 25);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/specfun/cross.m"), "resolved",
                  "resolved", 25);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1286826042U), "fileTimeLo",
                  "fileTimeLo", 25);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 25);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 25);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 25);
  sf_mex_assign(&c3_rhs25, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs25, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs25), "rhs", "rhs",
                  25);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs25), "lhs", "lhs",
                  25);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "context", "context", 26);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("sqrt"), "name", "name", 26);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 26);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "resolved",
                  "resolved", 26);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1343837586U), "fileTimeLo",
                  "fileTimeLo", 26);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 26);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 26);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 26);
  sf_mex_assign(&c3_rhs26, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs26, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs26), "rhs", "rhs",
                  26);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs26), "lhs", "lhs",
                  26);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "context",
                  "context", 27);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_error"), "name", "name",
                  27);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 27);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_error.m"), "resolved",
                  "resolved", 27);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1343837558U), "fileTimeLo",
                  "fileTimeLo", 27);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 27);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 27);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 27);
  sf_mex_assign(&c3_rhs27, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs27, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs27), "rhs", "rhs",
                  27);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs27), "lhs", "lhs",
                  27);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "context",
                  "context", 28);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_scalar_sqrt"), "name",
                  "name", 28);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 28);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sqrt.m"),
                  "resolved", "resolved", 28);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1286825938U), "fileTimeLo",
                  "fileTimeLo", 28);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 28);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 28);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 28);
  sf_mex_assign(&c3_rhs28, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs28, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs28), "rhs", "rhs",
                  28);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs28), "lhs", "lhs",
                  28);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "context", "context", 29);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("sign"), "name", "name", 29);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 29);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sign.m"), "resolved",
                  "resolved", 29);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363717456U), "fileTimeLo",
                  "fileTimeLo", 29);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 29);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 29);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 29);
  sf_mex_assign(&c3_rhs29, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs29, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs29), "rhs", "rhs",
                  29);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs29), "lhs", "lhs",
                  29);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sign.m"), "context",
                  "context", 30);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 30);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 30);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 30);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 30);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 30);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 30);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 30);
  sf_mex_assign(&c3_rhs30, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs30, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs30), "rhs", "rhs",
                  30);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs30), "lhs", "lhs",
                  30);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sign.m"), "context",
                  "context", 31);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_scalar_sign"), "name",
                  "name", 31);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 31);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sign.m"),
                  "resolved", "resolved", 31);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1356545094U), "fileTimeLo",
                  "fileTimeLo", 31);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 31);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 31);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 31);
  sf_mex_assign(&c3_rhs31, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs31, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs31), "rhs", "rhs",
                  31);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs31), "lhs", "lhs",
                  31);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "context", "context", 32);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 32);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 32);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 32);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1383880894U), "fileTimeLo",
                  "fileTimeLo", 32);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 32);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 32);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 32);
  sf_mex_assign(&c3_rhs32, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs32, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs32), "rhs", "rhs",
                  32);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs32), "lhs", "lhs",
                  32);
  sf_mex_destroy(&c3_rhs0);
  sf_mex_destroy(&c3_lhs0);
  sf_mex_destroy(&c3_rhs1);
  sf_mex_destroy(&c3_lhs1);
  sf_mex_destroy(&c3_rhs2);
  sf_mex_destroy(&c3_lhs2);
  sf_mex_destroy(&c3_rhs3);
  sf_mex_destroy(&c3_lhs3);
  sf_mex_destroy(&c3_rhs4);
  sf_mex_destroy(&c3_lhs4);
  sf_mex_destroy(&c3_rhs5);
  sf_mex_destroy(&c3_lhs5);
  sf_mex_destroy(&c3_rhs6);
  sf_mex_destroy(&c3_lhs6);
  sf_mex_destroy(&c3_rhs7);
  sf_mex_destroy(&c3_lhs7);
  sf_mex_destroy(&c3_rhs8);
  sf_mex_destroy(&c3_lhs8);
  sf_mex_destroy(&c3_rhs9);
  sf_mex_destroy(&c3_lhs9);
  sf_mex_destroy(&c3_rhs10);
  sf_mex_destroy(&c3_lhs10);
  sf_mex_destroy(&c3_rhs11);
  sf_mex_destroy(&c3_lhs11);
  sf_mex_destroy(&c3_rhs12);
  sf_mex_destroy(&c3_lhs12);
  sf_mex_destroy(&c3_rhs13);
  sf_mex_destroy(&c3_lhs13);
  sf_mex_destroy(&c3_rhs14);
  sf_mex_destroy(&c3_lhs14);
  sf_mex_destroy(&c3_rhs15);
  sf_mex_destroy(&c3_lhs15);
  sf_mex_destroy(&c3_rhs16);
  sf_mex_destroy(&c3_lhs16);
  sf_mex_destroy(&c3_rhs17);
  sf_mex_destroy(&c3_lhs17);
  sf_mex_destroy(&c3_rhs18);
  sf_mex_destroy(&c3_lhs18);
  sf_mex_destroy(&c3_rhs19);
  sf_mex_destroy(&c3_lhs19);
  sf_mex_destroy(&c3_rhs20);
  sf_mex_destroy(&c3_lhs20);
  sf_mex_destroy(&c3_rhs21);
  sf_mex_destroy(&c3_lhs21);
  sf_mex_destroy(&c3_rhs22);
  sf_mex_destroy(&c3_lhs22);
  sf_mex_destroy(&c3_rhs23);
  sf_mex_destroy(&c3_lhs23);
  sf_mex_destroy(&c3_rhs24);
  sf_mex_destroy(&c3_lhs24);
  sf_mex_destroy(&c3_rhs25);
  sf_mex_destroy(&c3_lhs25);
  sf_mex_destroy(&c3_rhs26);
  sf_mex_destroy(&c3_lhs26);
  sf_mex_destroy(&c3_rhs27);
  sf_mex_destroy(&c3_lhs27);
  sf_mex_destroy(&c3_rhs28);
  sf_mex_destroy(&c3_lhs28);
  sf_mex_destroy(&c3_rhs29);
  sf_mex_destroy(&c3_lhs29);
  sf_mex_destroy(&c3_rhs30);
  sf_mex_destroy(&c3_lhs30);
  sf_mex_destroy(&c3_rhs31);
  sf_mex_destroy(&c3_lhs31);
  sf_mex_destroy(&c3_rhs32);
  sf_mex_destroy(&c3_lhs32);
}

static const mxArray *c3_emlrt_marshallOut(const char * c3_u)
{
  const mxArray *c3_y = NULL;
  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", c3_u, 15, 0U, 0U, 0U, 2, 1, strlen
    (c3_u)), false);
  return c3_y;
}

static const mxArray *c3_b_emlrt_marshallOut(const uint32_T c3_u)
{
  const mxArray *c3_y = NULL;
  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", &c3_u, 7, 0U, 0U, 0U, 0), false);
  return c3_y;
}

static void c3_eml_scalar_eg(SFc3_simlwrkuka_kinematicsInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void c3_eml_xgemm(SFc3_simlwrkuka_kinematicsInstanceStruct *chartInstance,
  real_T c3_A[16], real_T c3_B[16], real_T c3_C[16], real_T c3_b_C[16])
{
  int32_T c3_i264;
  int32_T c3_i265;
  real_T c3_b_A[16];
  int32_T c3_i266;
  real_T c3_b_B[16];
  for (c3_i264 = 0; c3_i264 < 16; c3_i264++) {
    c3_b_C[c3_i264] = c3_C[c3_i264];
  }

  for (c3_i265 = 0; c3_i265 < 16; c3_i265++) {
    c3_b_A[c3_i265] = c3_A[c3_i265];
  }

  for (c3_i266 = 0; c3_i266 < 16; c3_i266++) {
    c3_b_B[c3_i266] = c3_B[c3_i266];
  }

  c3_b_eml_xgemm(chartInstance, c3_b_A, c3_b_B, c3_b_C);
}

static void c3_eml_error(SFc3_simlwrkuka_kinematicsInstanceStruct *chartInstance)
{
  int32_T c3_i267;
  static char_T c3_cv0[30] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o', 'l',
    'b', 'o', 'x', ':', 'E', 'l', 'F', 'u', 'n', 'D', 'o', 'm', 'a', 'i', 'n',
    'E', 'r', 'r', 'o', 'r' };

  char_T c3_u[30];
  const mxArray *c3_y = NULL;
  int32_T c3_i268;
  static char_T c3_cv1[4] = { 's', 'q', 'r', 't' };

  char_T c3_b_u[4];
  const mxArray *c3_b_y = NULL;
  (void)chartInstance;
  for (c3_i267 = 0; c3_i267 < 30; c3_i267++) {
    c3_u[c3_i267] = c3_cv0[c3_i267];
  }

  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", c3_u, 10, 0U, 1U, 0U, 2, 1, 30), false);
  for (c3_i268 = 0; c3_i268 < 4; c3_i268++) {
    c3_b_u[c3_i268] = c3_cv1[c3_i268];
  }

  c3_b_y = NULL;
  sf_mex_assign(&c3_b_y, sf_mex_create("y", c3_b_u, 10, 0U, 1U, 0U, 2, 1, 4),
                false);
  sf_mex_call_debug(sfGlobalDebugInstanceStruct, "error", 0U, 1U, 14,
                    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message", 1U,
    2U, 14, c3_y, 14, c3_b_y));
}

static const mxArray *c3_h_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData)
{
  const mxArray *c3_mxArrayOutData = NULL;
  int32_T c3_u;
  const mxArray *c3_y = NULL;
  SFc3_simlwrkuka_kinematicsInstanceStruct *chartInstance;
  chartInstance = (SFc3_simlwrkuka_kinematicsInstanceStruct *)chartInstanceVoid;
  c3_mxArrayOutData = NULL;
  c3_u = *(int32_T *)c3_inData;
  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", &c3_u, 6, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c3_mxArrayOutData, c3_y, false);
  return c3_mxArrayOutData;
}

static int32_T c3_j_emlrt_marshallIn(SFc3_simlwrkuka_kinematicsInstanceStruct
  *chartInstance, const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId)
{
  int32_T c3_y;
  int32_T c3_i269;
  (void)chartInstance;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), &c3_i269, 1, 6, 0U, 0, 0U, 0);
  c3_y = c3_i269;
  sf_mex_destroy(&c3_u);
  return c3_y;
}

static void c3_h_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData)
{
  const mxArray *c3_b_sfEvent;
  const char_T *c3_identifier;
  emlrtMsgIdentifier c3_thisId;
  int32_T c3_y;
  SFc3_simlwrkuka_kinematicsInstanceStruct *chartInstance;
  chartInstance = (SFc3_simlwrkuka_kinematicsInstanceStruct *)chartInstanceVoid;
  c3_b_sfEvent = sf_mex_dup(c3_mxArrayInData);
  c3_identifier = c3_varName;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_y = c3_j_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_b_sfEvent),
    &c3_thisId);
  sf_mex_destroy(&c3_b_sfEvent);
  *(int32_T *)c3_outData = c3_y;
  sf_mex_destroy(&c3_mxArrayInData);
}

static uint8_T c3_k_emlrt_marshallIn(SFc3_simlwrkuka_kinematicsInstanceStruct
  *chartInstance, const mxArray *c3_b_is_active_c3_simlwrkuka_kinematics, const
  char_T *c3_identifier)
{
  uint8_T c3_y;
  emlrtMsgIdentifier c3_thisId;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_y = c3_l_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c3_b_is_active_c3_simlwrkuka_kinematics), &c3_thisId);
  sf_mex_destroy(&c3_b_is_active_c3_simlwrkuka_kinematics);
  return c3_y;
}

static uint8_T c3_l_emlrt_marshallIn(SFc3_simlwrkuka_kinematicsInstanceStruct
  *chartInstance, const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId)
{
  uint8_T c3_y;
  uint8_T c3_u0;
  (void)chartInstance;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), &c3_u0, 1, 3, 0U, 0, 0U, 0);
  c3_y = c3_u0;
  sf_mex_destroy(&c3_u);
  return c3_y;
}

static void c3_b_eml_xgemm(SFc3_simlwrkuka_kinematicsInstanceStruct
  *chartInstance, real_T c3_A[16], real_T c3_B[16], real_T c3_C[16])
{
  int32_T c3_i270;
  int32_T c3_i271;
  int32_T c3_i272;
  int32_T c3_i273;
  int32_T c3_i274;
  (void)chartInstance;
  for (c3_i270 = 0; c3_i270 < 4; c3_i270++) {
    c3_i271 = 0;
    for (c3_i272 = 0; c3_i272 < 4; c3_i272++) {
      c3_C[c3_i271 + c3_i270] = 0.0;
      c3_i273 = 0;
      for (c3_i274 = 0; c3_i274 < 4; c3_i274++) {
        c3_C[c3_i271 + c3_i270] += c3_A[c3_i273 + c3_i270] * c3_B[c3_i274 +
          c3_i271];
        c3_i273 += 4;
      }

      c3_i271 += 4;
    }
  }
}

static void init_dsm_address_info(SFc3_simlwrkuka_kinematicsInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

/* SFunction Glue Code */
#ifdef utFree
#undef utFree
#endif

#ifdef utMalloc
#undef utMalloc
#endif

#ifdef __cplusplus

extern "C" void *utMalloc(size_t size);
extern "C" void utFree(void*);

#else

extern void *utMalloc(size_t size);
extern void utFree(void*);

#endif

void sf_c3_simlwrkuka_kinematics_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(512527961U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(1795202584U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(405859728U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(2565562944U);
}

mxArray *sf_c3_simlwrkuka_kinematics_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1,1,5,
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("hTpfPrOMuMHSkShfpYu0WB");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,1,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(7);
      pr[1] = (double)(1);
      mxSetField(mxData,0,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,0,"type",mxType);
    }

    mxSetField(mxData,0,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"inputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"parameters",mxCreateDoubleMatrix(0,0,
                mxREAL));
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,3,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(7);
      pr[1] = (double)(1);
      mxSetField(mxData,0,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,0,"type",mxType);
    }

    mxSetField(mxData,0,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(3);
      pr[1] = (double)(1);
      mxSetField(mxData,1,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,1,"type",mxType);
    }

    mxSetField(mxData,1,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(3);
      pr[1] = (double)(1);
      mxSetField(mxData,2,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,2,"type",mxType);
    }

    mxSetField(mxData,2,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"outputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"locals",mxCreateDoubleMatrix(0,0,mxREAL));
  }

  return(mxAutoinheritanceInfo);
}

mxArray *sf_c3_simlwrkuka_kinematics_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c3_simlwrkuka_kinematics_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

static const mxArray *sf_get_sim_state_info_c3_simlwrkuka_kinematics(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x4'type','srcId','name','auxInfo'{{M[1],M[5],T\"x\",},{M[1],M[9],T\"xe\",},{M[1],M[8],T\"xp\",},{M[8],M[0],T\"is_active_c3_simlwrkuka_kinematics\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 4, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c3_simlwrkuka_kinematics_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc3_simlwrkuka_kinematicsInstanceStruct *chartInstance;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
    chartInstance = (SFc3_simlwrkuka_kinematicsInstanceStruct *)
      chartInfo->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _simlwrkuka_kinematicsMachineNumber_,
           3,
           1,
           1,
           0,
           4,
           0,
           0,
           0,
           0,
           1,
           &(chartInstance->chartNumber),
           &(chartInstance->instanceNumber),
           (void *)S);

        /* Each instance must initialize ist own list of scripts */
        init_script_number_translation(_simlwrkuka_kinematicsMachineNumber_,
          chartInstance->chartNumber,chartInstance->instanceNumber);
        if (chartAlreadyPresent==0) {
          /* this is the first instance */
          sf_debug_set_chart_disable_implicit_casting
            (sfGlobalDebugInstanceStruct,_simlwrkuka_kinematicsMachineNumber_,
             chartInstance->chartNumber,1);
          sf_debug_set_chart_event_thresholds(sfGlobalDebugInstanceStruct,
            _simlwrkuka_kinematicsMachineNumber_,
            chartInstance->chartNumber,
            0,
            0,
            0);
          _SFD_SET_DATA_PROPS(0,1,1,0,"q");
          _SFD_SET_DATA_PROPS(1,2,0,1,"x");
          _SFD_SET_DATA_PROPS(2,2,0,1,"xp");
          _SFD_SET_DATA_PROPS(3,2,0,1,"xe");
          _SFD_STATE_INFO(0,0,2);
          _SFD_CH_SUBSTATE_COUNT(0);
          _SFD_CH_SUBSTATE_DECOMP(0);
        }

        _SFD_CV_INIT_CHART(0,0,0,0);

        {
          _SFD_CV_INIT_STATE(0,0,0,0,0,0,NULL,NULL);
        }

        _SFD_CV_INIT_TRANS(0,0,NULL,NULL,0,NULL);

        /* Initialization of MATLAB Function Model Coverage */
        _SFD_CV_INIT_EML(0,1,1,0,0,0,0,0,0,0,0);
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,722);
        _SFD_CV_INIT_SCRIPT(0,1,0,0,0,0,0,0,0,0);
        _SFD_CV_INIT_SCRIPT_FCN(0,0,"jacobin",0,-1,2531);

        {
          unsigned int dimVector[1];
          dimVector[0]= 7;
          _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c3_b_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 7;
          _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c3_b_sf_marshallOut,(MexInFcnForType)
            c3_b_sf_marshallIn);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c3_sf_marshallOut,(MexInFcnForType)
            c3_sf_marshallIn);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c3_sf_marshallOut,(MexInFcnForType)
            c3_sf_marshallIn);
        }

        {
          real_T (*c3_q)[7];
          real_T (*c3_x)[7];
          real_T (*c3_xp)[3];
          real_T (*c3_xe)[3];
          c3_xe = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S, 3);
          c3_xp = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S, 2);
          c3_x = (real_T (*)[7])ssGetOutputPortSignal(chartInstance->S, 1);
          c3_q = (real_T (*)[7])ssGetInputPortSignal(chartInstance->S, 0);
          _SFD_SET_DATA_VALUE_PTR(0U, *c3_q);
          _SFD_SET_DATA_VALUE_PTR(1U, *c3_x);
          _SFD_SET_DATA_VALUE_PTR(2U, *c3_xp);
          _SFD_SET_DATA_VALUE_PTR(3U, *c3_xe);
        }
      }
    } else {
      sf_debug_reset_current_state_configuration(sfGlobalDebugInstanceStruct,
        _simlwrkuka_kinematicsMachineNumber_,chartInstance->chartNumber,
        chartInstance->instanceNumber);
    }
  }
}

static const char* sf_get_instance_specialization(void)
{
  return "5V0DCKpOfiNTHguuu1xyoC";
}

static void sf_opaque_initialize_c3_simlwrkuka_kinematics(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc3_simlwrkuka_kinematicsInstanceStruct*)
    chartInstanceVar)->S,0);
  initialize_params_c3_simlwrkuka_kinematics
    ((SFc3_simlwrkuka_kinematicsInstanceStruct*) chartInstanceVar);
  initialize_c3_simlwrkuka_kinematics((SFc3_simlwrkuka_kinematicsInstanceStruct*)
    chartInstanceVar);
}

static void sf_opaque_enable_c3_simlwrkuka_kinematics(void *chartInstanceVar)
{
  enable_c3_simlwrkuka_kinematics((SFc3_simlwrkuka_kinematicsInstanceStruct*)
    chartInstanceVar);
}

static void sf_opaque_disable_c3_simlwrkuka_kinematics(void *chartInstanceVar)
{
  disable_c3_simlwrkuka_kinematics((SFc3_simlwrkuka_kinematicsInstanceStruct*)
    chartInstanceVar);
}

static void sf_opaque_gateway_c3_simlwrkuka_kinematics(void *chartInstanceVar)
{
  sf_gateway_c3_simlwrkuka_kinematics((SFc3_simlwrkuka_kinematicsInstanceStruct*)
    chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c3_simlwrkuka_kinematics
  (SimStruct* S)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c3_simlwrkuka_kinematics
    ((SFc3_simlwrkuka_kinematicsInstanceStruct*)chartInfo->chartInstance);/* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c3_simlwrkuka_kinematics();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 4, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  mxDestroyArray(prhs[3]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_raw2high'.\n");
  }

  return plhs[0];
}

extern void sf_internal_set_sim_state_c3_simlwrkuka_kinematics(SimStruct* S,
  const mxArray *st)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[3];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxDuplicateArray(st);      /* high level simctx */
  prhs[2] = (mxArray*) sf_get_sim_state_info_c3_simlwrkuka_kinematics();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 3, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c3_simlwrkuka_kinematics
    ((SFc3_simlwrkuka_kinematicsInstanceStruct*)chartInfo->chartInstance,
     mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c3_simlwrkuka_kinematics(SimStruct*
  S)
{
  return sf_internal_get_sim_state_c3_simlwrkuka_kinematics(S);
}

static void sf_opaque_set_sim_state_c3_simlwrkuka_kinematics(SimStruct* S, const
  mxArray *st)
{
  sf_internal_set_sim_state_c3_simlwrkuka_kinematics(S, st);
}

static void sf_opaque_terminate_c3_simlwrkuka_kinematics(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc3_simlwrkuka_kinematicsInstanceStruct*) chartInstanceVar)
      ->S;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_simlwrkuka_kinematics_optimization_info();
    }

    finalize_c3_simlwrkuka_kinematics((SFc3_simlwrkuka_kinematicsInstanceStruct*)
      chartInstanceVar);
    utFree((void *)chartInstanceVar);
    if (crtInfo != NULL) {
      utFree((void *)crtInfo);
    }

    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc3_simlwrkuka_kinematics
    ((SFc3_simlwrkuka_kinematicsInstanceStruct*) chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c3_simlwrkuka_kinematics(SimStruct *S)
{
  int i;
  for (i=0;i<ssGetNumRunTimeParams(S);i++) {
    if (ssGetSFcnParamTunable(S,i)) {
      ssUpdateDlgParamAsRunTimeParam(S,i);
    }
  }

  if (sf_machine_global_initializer_called()) {
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
    initialize_params_c3_simlwrkuka_kinematics
      ((SFc3_simlwrkuka_kinematicsInstanceStruct*)(chartInfo->chartInstance));
  }
}

static void mdlSetWorkWidths_c3_simlwrkuka_kinematics(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_simlwrkuka_kinematics_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(sf_get_instance_specialization(),infoStruct,3);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(sf_get_instance_specialization(),
                infoStruct,3,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop
      (sf_get_instance_specialization(),infoStruct,3,
       "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(sf_get_instance_specialization(),infoStruct,3);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,3,1);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,3,3);
    }

    {
      unsigned int outPortIdx;
      for (outPortIdx=1; outPortIdx<=3; ++outPortIdx) {
        ssSetOutputPortOptimizeInIR(S, outPortIdx, 1U);
      }
    }

    {
      unsigned int inPortIdx;
      for (inPortIdx=0; inPortIdx < 1; ++inPortIdx) {
        ssSetInputPortOptimizeInIR(S, inPortIdx, 1U);
      }
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,3);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(1233942600U));
  ssSetChecksum1(S,(646405626U));
  ssSetChecksum2(S,(2364274580U));
  ssSetChecksum3(S,(454570854U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c3_simlwrkuka_kinematics(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c3_simlwrkuka_kinematics(SimStruct *S)
{
  SFc3_simlwrkuka_kinematicsInstanceStruct *chartInstance;
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)utMalloc(sizeof
    (ChartRunTimeInfo));
  chartInstance = (SFc3_simlwrkuka_kinematicsInstanceStruct *)utMalloc(sizeof
    (SFc3_simlwrkuka_kinematicsInstanceStruct));
  memset(chartInstance, 0, sizeof(SFc3_simlwrkuka_kinematicsInstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway =
    sf_opaque_gateway_c3_simlwrkuka_kinematics;
  chartInstance->chartInfo.initializeChart =
    sf_opaque_initialize_c3_simlwrkuka_kinematics;
  chartInstance->chartInfo.terminateChart =
    sf_opaque_terminate_c3_simlwrkuka_kinematics;
  chartInstance->chartInfo.enableChart =
    sf_opaque_enable_c3_simlwrkuka_kinematics;
  chartInstance->chartInfo.disableChart =
    sf_opaque_disable_c3_simlwrkuka_kinematics;
  chartInstance->chartInfo.getSimState =
    sf_opaque_get_sim_state_c3_simlwrkuka_kinematics;
  chartInstance->chartInfo.setSimState =
    sf_opaque_set_sim_state_c3_simlwrkuka_kinematics;
  chartInstance->chartInfo.getSimStateInfo =
    sf_get_sim_state_info_c3_simlwrkuka_kinematics;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c3_simlwrkuka_kinematics;
  chartInstance->chartInfo.mdlStart = mdlStart_c3_simlwrkuka_kinematics;
  chartInstance->chartInfo.mdlSetWorkWidths =
    mdlSetWorkWidths_c3_simlwrkuka_kinematics;
  chartInstance->chartInfo.extModeExec = NULL;
  chartInstance->chartInfo.restoreLastMajorStepConfiguration = NULL;
  chartInstance->chartInfo.restoreBeforeLastMajorStepConfiguration = NULL;
  chartInstance->chartInfo.storeCurrentConfiguration = NULL;
  chartInstance->chartInfo.debugInstance = sfGlobalDebugInstanceStruct;
  chartInstance->S = S;
  crtInfo->instanceInfo = (&(chartInstance->chartInfo));
  crtInfo->isJITEnabled = false;
  ssSetUserData(S,(void *)(crtInfo));  /* register the chart instance with simstruct */
  init_dsm_address_info(chartInstance);
  if (!sim_mode_is_rtw_gen(S)) {
  }

  sf_opaque_init_subchart_simstructs(chartInstance->chartInfo.chartInstance);
  chart_debug_initialization(S,1);
}

void c3_simlwrkuka_kinematics_method_dispatcher(SimStruct *S, int_T method, void
  *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c3_simlwrkuka_kinematics(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c3_simlwrkuka_kinematics(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c3_simlwrkuka_kinematics(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c3_simlwrkuka_kinematics_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
