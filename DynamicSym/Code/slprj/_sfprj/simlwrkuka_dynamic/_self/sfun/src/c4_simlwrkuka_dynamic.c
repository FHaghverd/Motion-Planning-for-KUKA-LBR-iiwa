/* Include files */

#include <stddef.h>
#include "blas.h"
#include "simlwrkuka_dynamic_sfun.h"
#include "c4_simlwrkuka_dynamic.h"
#include "mwmathutil.h"
#define CHARTINSTANCE_CHARTNUMBER      (chartInstance->chartNumber)
#define CHARTINSTANCE_INSTANCENUMBER   (chartInstance->instanceNumber)
#include "simlwrkuka_dynamic_sfun_debug_macros.h"
#define _SF_MEX_LISTEN_FOR_CTRL_C(S)   sf_mex_listen_for_ctrl_c(sfGlobalDebugInstanceStruct,S);

/* Type Definitions */

/* Named Constants */
#define CALL_EVENT                     (-1)

/* Variable Declarations */

/* Variable Definitions */
static real_T _sfTime_;
static const char * c4_debug_family_names[20] = { "q1", "q2", "q3", "q4", "q5",
  "q6", "q7", "J", "Te", "pe", "Re", "eta", "e", "Q", "nargin", "nargout", "q",
  "dq", "x", "xd" };

static const char * c4_b_debug_family_names[50] = { "q", "d3", "d5", "d0", "d7",
  "f", "a", "d", "t", "A0", "A1", "A2", "A3", "A4", "A5", "A6", "A7", "Ae", "T1",
  "T2", "T3", "T4", "T5", "T6", "z0", "z1", "z2", "z3", "z4", "z5", "z6", "p0",
  "p1", "p2", "p3", "p4", "p5", "p6", "pe", "nargin", "nargout", "q1", "q2",
  "q3", "q4", "q5", "q6", "q7", "J", "Te" };

/* Function Declarations */
static void initialize_c4_simlwrkuka_dynamic
  (SFc4_simlwrkuka_dynamicInstanceStruct *chartInstance);
static void initialize_params_c4_simlwrkuka_dynamic
  (SFc4_simlwrkuka_dynamicInstanceStruct *chartInstance);
static void enable_c4_simlwrkuka_dynamic(SFc4_simlwrkuka_dynamicInstanceStruct
  *chartInstance);
static void disable_c4_simlwrkuka_dynamic(SFc4_simlwrkuka_dynamicInstanceStruct *
  chartInstance);
static void c4_update_debugger_state_c4_simlwrkuka_dynamic
  (SFc4_simlwrkuka_dynamicInstanceStruct *chartInstance);
static const mxArray *get_sim_state_c4_simlwrkuka_dynamic
  (SFc4_simlwrkuka_dynamicInstanceStruct *chartInstance);
static void set_sim_state_c4_simlwrkuka_dynamic
  (SFc4_simlwrkuka_dynamicInstanceStruct *chartInstance, const mxArray *c4_st);
static void finalize_c4_simlwrkuka_dynamic(SFc4_simlwrkuka_dynamicInstanceStruct
  *chartInstance);
static void sf_gateway_c4_simlwrkuka_dynamic
  (SFc4_simlwrkuka_dynamicInstanceStruct *chartInstance);
static void c4_chartstep_c4_simlwrkuka_dynamic
  (SFc4_simlwrkuka_dynamicInstanceStruct *chartInstance);
static void initSimStructsc4_simlwrkuka_dynamic
  (SFc4_simlwrkuka_dynamicInstanceStruct *chartInstance);
static void c4_jacobin(SFc4_simlwrkuka_dynamicInstanceStruct *chartInstance,
  real_T c4_q1, real_T c4_q2, real_T c4_q3, real_T c4_q4, real_T c4_q5, real_T
  c4_q6, real_T c4_q7, real_T c4_J[42], real_T c4_Te[16]);
static void init_script_number_translation(uint32_T c4_machineNumber, uint32_T
  c4_chartNumber, uint32_T c4_instanceNumber);
static const mxArray *c4_sf_marshallOut(void *chartInstanceVoid, void *c4_inData);
static void c4_emlrt_marshallIn(SFc4_simlwrkuka_dynamicInstanceStruct
  *chartInstance, const mxArray *c4_xd, const char_T *c4_identifier, real_T
  c4_y[6]);
static void c4_b_emlrt_marshallIn(SFc4_simlwrkuka_dynamicInstanceStruct
  *chartInstance, const mxArray *c4_u, const emlrtMsgIdentifier *c4_parentId,
  real_T c4_y[6]);
static void c4_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c4_mxArrayInData, const char_T *c4_varName, void *c4_outData);
static const mxArray *c4_b_sf_marshallOut(void *chartInstanceVoid, void
  *c4_inData);
static void c4_c_emlrt_marshallIn(SFc4_simlwrkuka_dynamicInstanceStruct
  *chartInstance, const mxArray *c4_x, const char_T *c4_identifier, real_T c4_y
  [7]);
static void c4_d_emlrt_marshallIn(SFc4_simlwrkuka_dynamicInstanceStruct
  *chartInstance, const mxArray *c4_u, const emlrtMsgIdentifier *c4_parentId,
  real_T c4_y[7]);
static void c4_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c4_mxArrayInData, const char_T *c4_varName, void *c4_outData);
static const mxArray *c4_c_sf_marshallOut(void *chartInstanceVoid, void
  *c4_inData);
static real_T c4_e_emlrt_marshallIn(SFc4_simlwrkuka_dynamicInstanceStruct
  *chartInstance, const mxArray *c4_u, const emlrtMsgIdentifier *c4_parentId);
static void c4_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c4_mxArrayInData, const char_T *c4_varName, void *c4_outData);
static const mxArray *c4_d_sf_marshallOut(void *chartInstanceVoid, void
  *c4_inData);
static void c4_f_emlrt_marshallIn(SFc4_simlwrkuka_dynamicInstanceStruct
  *chartInstance, const mxArray *c4_u, const emlrtMsgIdentifier *c4_parentId,
  real_T c4_y[4]);
static void c4_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c4_mxArrayInData, const char_T *c4_varName, void *c4_outData);
static const mxArray *c4_e_sf_marshallOut(void *chartInstanceVoid, void
  *c4_inData);
static void c4_g_emlrt_marshallIn(SFc4_simlwrkuka_dynamicInstanceStruct
  *chartInstance, const mxArray *c4_u, const emlrtMsgIdentifier *c4_parentId,
  real_T c4_y[3]);
static void c4_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c4_mxArrayInData, const char_T *c4_varName, void *c4_outData);
static const mxArray *c4_f_sf_marshallOut(void *chartInstanceVoid, void
  *c4_inData);
static void c4_h_emlrt_marshallIn(SFc4_simlwrkuka_dynamicInstanceStruct
  *chartInstance, const mxArray *c4_u, const emlrtMsgIdentifier *c4_parentId,
  real_T c4_y[9]);
static void c4_f_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c4_mxArrayInData, const char_T *c4_varName, void *c4_outData);
static const mxArray *c4_g_sf_marshallOut(void *chartInstanceVoid, void
  *c4_inData);
static void c4_i_emlrt_marshallIn(SFc4_simlwrkuka_dynamicInstanceStruct
  *chartInstance, const mxArray *c4_u, const emlrtMsgIdentifier *c4_parentId,
  real_T c4_y[16]);
static void c4_g_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c4_mxArrayInData, const char_T *c4_varName, void *c4_outData);
static const mxArray *c4_h_sf_marshallOut(void *chartInstanceVoid, void
  *c4_inData);
static void c4_j_emlrt_marshallIn(SFc4_simlwrkuka_dynamicInstanceStruct
  *chartInstance, const mxArray *c4_u, const emlrtMsgIdentifier *c4_parentId,
  real_T c4_y[42]);
static void c4_h_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c4_mxArrayInData, const char_T *c4_varName, void *c4_outData);
static void c4_info_helper(const mxArray **c4_info);
static const mxArray *c4_emlrt_marshallOut(const char * c4_u);
static const mxArray *c4_b_emlrt_marshallOut(const uint32_T c4_u);
static void c4_eml_scalar_eg(SFc4_simlwrkuka_dynamicInstanceStruct
  *chartInstance);
static void c4_threshold(SFc4_simlwrkuka_dynamicInstanceStruct *chartInstance);
static void c4_eml_error(SFc4_simlwrkuka_dynamicInstanceStruct *chartInstance);
static void c4_b_eml_scalar_eg(SFc4_simlwrkuka_dynamicInstanceStruct
  *chartInstance);
static const mxArray *c4_i_sf_marshallOut(void *chartInstanceVoid, void
  *c4_inData);
static int32_T c4_k_emlrt_marshallIn(SFc4_simlwrkuka_dynamicInstanceStruct
  *chartInstance, const mxArray *c4_u, const emlrtMsgIdentifier *c4_parentId);
static void c4_i_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c4_mxArrayInData, const char_T *c4_varName, void *c4_outData);
static uint8_T c4_l_emlrt_marshallIn(SFc4_simlwrkuka_dynamicInstanceStruct
  *chartInstance, const mxArray *c4_b_is_active_c4_simlwrkuka_dynamic, const
  char_T *c4_identifier);
static uint8_T c4_m_emlrt_marshallIn(SFc4_simlwrkuka_dynamicInstanceStruct
  *chartInstance, const mxArray *c4_u, const emlrtMsgIdentifier *c4_parentId);
static void init_dsm_address_info(SFc4_simlwrkuka_dynamicInstanceStruct
  *chartInstance);

/* Function Definitions */
static void initialize_c4_simlwrkuka_dynamic
  (SFc4_simlwrkuka_dynamicInstanceStruct *chartInstance)
{
  chartInstance->c4_sfEvent = CALL_EVENT;
  _sfTime_ = sf_get_time(chartInstance->S);
  chartInstance->c4_is_active_c4_simlwrkuka_dynamic = 0U;
}

static void initialize_params_c4_simlwrkuka_dynamic
  (SFc4_simlwrkuka_dynamicInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void enable_c4_simlwrkuka_dynamic(SFc4_simlwrkuka_dynamicInstanceStruct
  *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void disable_c4_simlwrkuka_dynamic(SFc4_simlwrkuka_dynamicInstanceStruct *
  chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void c4_update_debugger_state_c4_simlwrkuka_dynamic
  (SFc4_simlwrkuka_dynamicInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static const mxArray *get_sim_state_c4_simlwrkuka_dynamic
  (SFc4_simlwrkuka_dynamicInstanceStruct *chartInstance)
{
  const mxArray *c4_st;
  const mxArray *c4_y = NULL;
  int32_T c4_i0;
  real_T c4_u[7];
  const mxArray *c4_b_y = NULL;
  int32_T c4_i1;
  real_T c4_b_u[6];
  const mxArray *c4_c_y = NULL;
  uint8_T c4_hoistedGlobal;
  uint8_T c4_c_u;
  const mxArray *c4_d_y = NULL;
  real_T (*c4_xd)[6];
  real_T (*c4_x)[7];
  c4_xd = (real_T (*)[6])ssGetOutputPortSignal(chartInstance->S, 2);
  c4_x = (real_T (*)[7])ssGetOutputPortSignal(chartInstance->S, 1);
  c4_st = NULL;
  c4_st = NULL;
  c4_y = NULL;
  sf_mex_assign(&c4_y, sf_mex_createcellmatrix(3, 1), false);
  for (c4_i0 = 0; c4_i0 < 7; c4_i0++) {
    c4_u[c4_i0] = (*c4_x)[c4_i0];
  }

  c4_b_y = NULL;
  sf_mex_assign(&c4_b_y, sf_mex_create("y", c4_u, 0, 0U, 1U, 0U, 1, 7), false);
  sf_mex_setcell(c4_y, 0, c4_b_y);
  for (c4_i1 = 0; c4_i1 < 6; c4_i1++) {
    c4_b_u[c4_i1] = (*c4_xd)[c4_i1];
  }

  c4_c_y = NULL;
  sf_mex_assign(&c4_c_y, sf_mex_create("y", c4_b_u, 0, 0U, 1U, 0U, 1, 6), false);
  sf_mex_setcell(c4_y, 1, c4_c_y);
  c4_hoistedGlobal = chartInstance->c4_is_active_c4_simlwrkuka_dynamic;
  c4_c_u = c4_hoistedGlobal;
  c4_d_y = NULL;
  sf_mex_assign(&c4_d_y, sf_mex_create("y", &c4_c_u, 3, 0U, 0U, 0U, 0), false);
  sf_mex_setcell(c4_y, 2, c4_d_y);
  sf_mex_assign(&c4_st, c4_y, false);
  return c4_st;
}

static void set_sim_state_c4_simlwrkuka_dynamic
  (SFc4_simlwrkuka_dynamicInstanceStruct *chartInstance, const mxArray *c4_st)
{
  const mxArray *c4_u;
  real_T c4_dv0[7];
  int32_T c4_i2;
  real_T c4_dv1[6];
  int32_T c4_i3;
  real_T (*c4_x)[7];
  real_T (*c4_xd)[6];
  c4_xd = (real_T (*)[6])ssGetOutputPortSignal(chartInstance->S, 2);
  c4_x = (real_T (*)[7])ssGetOutputPortSignal(chartInstance->S, 1);
  chartInstance->c4_doneDoubleBufferReInit = true;
  c4_u = sf_mex_dup(c4_st);
  c4_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c4_u, 0)), "x",
                        c4_dv0);
  for (c4_i2 = 0; c4_i2 < 7; c4_i2++) {
    (*c4_x)[c4_i2] = c4_dv0[c4_i2];
  }

  c4_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c4_u, 1)), "xd",
                      c4_dv1);
  for (c4_i3 = 0; c4_i3 < 6; c4_i3++) {
    (*c4_xd)[c4_i3] = c4_dv1[c4_i3];
  }

  chartInstance->c4_is_active_c4_simlwrkuka_dynamic = c4_l_emlrt_marshallIn
    (chartInstance, sf_mex_dup(sf_mex_getcell(c4_u, 2)),
     "is_active_c4_simlwrkuka_dynamic");
  sf_mex_destroy(&c4_u);
  c4_update_debugger_state_c4_simlwrkuka_dynamic(chartInstance);
  sf_mex_destroy(&c4_st);
}

static void finalize_c4_simlwrkuka_dynamic(SFc4_simlwrkuka_dynamicInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void sf_gateway_c4_simlwrkuka_dynamic
  (SFc4_simlwrkuka_dynamicInstanceStruct *chartInstance)
{
  int32_T c4_i4;
  int32_T c4_i5;
  int32_T c4_i6;
  int32_T c4_i7;
  real_T (*c4_dq)[7];
  real_T (*c4_xd)[6];
  real_T (*c4_x)[7];
  real_T (*c4_q)[7];
  c4_dq = (real_T (*)[7])ssGetInputPortSignal(chartInstance->S, 1);
  c4_xd = (real_T (*)[6])ssGetOutputPortSignal(chartInstance->S, 2);
  c4_x = (real_T (*)[7])ssGetOutputPortSignal(chartInstance->S, 1);
  c4_q = (real_T (*)[7])ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = sf_get_time(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 3U, chartInstance->c4_sfEvent);
  for (c4_i4 = 0; c4_i4 < 7; c4_i4++) {
    _SFD_DATA_RANGE_CHECK((*c4_q)[c4_i4], 0U);
  }

  chartInstance->c4_sfEvent = CALL_EVENT;
  c4_chartstep_c4_simlwrkuka_dynamic(chartInstance);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_simlwrkuka_dynamicMachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
  for (c4_i5 = 0; c4_i5 < 7; c4_i5++) {
    _SFD_DATA_RANGE_CHECK((*c4_x)[c4_i5], 1U);
  }

  for (c4_i6 = 0; c4_i6 < 6; c4_i6++) {
    _SFD_DATA_RANGE_CHECK((*c4_xd)[c4_i6], 2U);
  }

  for (c4_i7 = 0; c4_i7 < 7; c4_i7++) {
    _SFD_DATA_RANGE_CHECK((*c4_dq)[c4_i7], 3U);
  }
}

static void c4_chartstep_c4_simlwrkuka_dynamic
  (SFc4_simlwrkuka_dynamicInstanceStruct *chartInstance)
{
  int32_T c4_i8;
  real_T c4_q[7];
  int32_T c4_i9;
  real_T c4_dq[7];
  uint32_T c4_debug_family_var_map[20];
  real_T c4_q1;
  real_T c4_q2;
  real_T c4_q3;
  real_T c4_q4;
  real_T c4_q5;
  real_T c4_q6;
  real_T c4_q7;
  real_T c4_J[42];
  real_T c4_Te[16];
  real_T c4_pe[3];
  real_T c4_Re[9];
  real_T c4_eta;
  real_T c4_e[3];
  real_T c4_Q[4];
  real_T c4_nargin = 2.0;
  real_T c4_nargout = 2.0;
  real_T c4_x[7];
  real_T c4_xd[6];
  real_T c4_b_Te[16];
  real_T c4_b_J[42];
  int32_T c4_i10;
  int32_T c4_i11;
  int32_T c4_i12;
  int32_T c4_i13;
  int32_T c4_i14;
  int32_T c4_i15;
  int32_T c4_i16;
  real_T c4_b_x;
  real_T c4_c_x;
  real_T c4_d_x;
  real_T c4_e_x;
  real_T c4_f_x;
  real_T c4_g_x;
  real_T c4_h_x;
  real_T c4_i_x;
  real_T c4_j_x;
  real_T c4_k_x;
  real_T c4_l_x;
  real_T c4_m_x;
  real_T c4_n_x;
  real_T c4_o_x;
  real_T c4_b[3];
  int32_T c4_i17;
  int32_T c4_i18;
  int32_T c4_i19;
  int32_T c4_i20;
  int32_T c4_i21;
  int32_T c4_i22;
  real_T c4_b_b[7];
  int32_T c4_i23;
  int32_T c4_i24;
  int32_T c4_i25;
  real_T c4_C[6];
  int32_T c4_i26;
  int32_T c4_i27;
  int32_T c4_i28;
  int32_T c4_i29;
  int32_T c4_i30;
  int32_T c4_i31;
  int32_T c4_i32;
  int32_T c4_i33;
  real_T (*c4_p_x)[7];
  real_T (*c4_b_xd)[6];
  real_T (*c4_b_dq)[7];
  real_T (*c4_b_q)[7];
  c4_b_dq = (real_T (*)[7])ssGetInputPortSignal(chartInstance->S, 1);
  c4_b_xd = (real_T (*)[6])ssGetOutputPortSignal(chartInstance->S, 2);
  c4_p_x = (real_T (*)[7])ssGetOutputPortSignal(chartInstance->S, 1);
  c4_b_q = (real_T (*)[7])ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 3U, chartInstance->c4_sfEvent);
  for (c4_i8 = 0; c4_i8 < 7; c4_i8++) {
    c4_q[c4_i8] = (*c4_b_q)[c4_i8];
  }

  for (c4_i9 = 0; c4_i9 < 7; c4_i9++) {
    c4_dq[c4_i9] = (*c4_b_dq)[c4_i9];
  }

  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 20U, 20U, c4_debug_family_names,
    c4_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c4_q1, 0U, c4_c_sf_marshallOut,
    c4_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c4_q2, 1U, c4_c_sf_marshallOut,
    c4_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c4_q3, 2U, c4_c_sf_marshallOut,
    c4_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c4_q4, 3U, c4_c_sf_marshallOut,
    c4_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c4_q5, 4U, c4_c_sf_marshallOut,
    c4_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c4_q6, 5U, c4_c_sf_marshallOut,
    c4_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c4_q7, 6U, c4_c_sf_marshallOut,
    c4_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c4_J, 7U, c4_h_sf_marshallOut,
    c4_h_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c4_Te, 8U, c4_g_sf_marshallOut,
    c4_g_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c4_pe, 9U, c4_e_sf_marshallOut,
    c4_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c4_Re, 10U, c4_f_sf_marshallOut,
    c4_f_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c4_eta, 11U, c4_c_sf_marshallOut,
    c4_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c4_e, 12U, c4_e_sf_marshallOut,
    c4_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c4_Q, 13U, c4_d_sf_marshallOut,
    c4_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c4_nargin, 14U, c4_c_sf_marshallOut,
    c4_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c4_nargout, 15U, c4_c_sf_marshallOut,
    c4_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c4_q, 16U, c4_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c4_dq, 17U, c4_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c4_x, 18U, c4_b_sf_marshallOut,
    c4_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c4_xd, 19U, c4_sf_marshallOut,
    c4_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c4_sfEvent, 2);
  c4_q1 = c4_q[0];
  _SFD_EML_CALL(0U, chartInstance->c4_sfEvent, 2);
  c4_q2 = c4_q[1];
  _SFD_EML_CALL(0U, chartInstance->c4_sfEvent, 2);
  c4_q3 = c4_q[2];
  _SFD_EML_CALL(0U, chartInstance->c4_sfEvent, 2);
  c4_q4 = c4_q[3];
  _SFD_EML_CALL(0U, chartInstance->c4_sfEvent, 2);
  c4_q5 = c4_q[4];
  _SFD_EML_CALL(0U, chartInstance->c4_sfEvent, 2);
  c4_q6 = c4_q[5];
  _SFD_EML_CALL(0U, chartInstance->c4_sfEvent, 2);
  c4_q7 = c4_q[6];
  _SFD_EML_CALL(0U, chartInstance->c4_sfEvent, 3);
  c4_jacobin(chartInstance, c4_q1, c4_q2, c4_q3, c4_q4, c4_q5, c4_q6, c4_q7,
             c4_b_J, c4_b_Te);
  for (c4_i10 = 0; c4_i10 < 42; c4_i10++) {
    c4_J[c4_i10] = c4_b_J[c4_i10];
  }

  for (c4_i11 = 0; c4_i11 < 16; c4_i11++) {
    c4_Te[c4_i11] = c4_b_Te[c4_i11];
  }

  _SFD_EML_CALL(0U, chartInstance->c4_sfEvent, 4);
  for (c4_i12 = 0; c4_i12 < 3; c4_i12++) {
    c4_pe[c4_i12] = c4_Te[c4_i12 + 12];
  }

  _SFD_EML_CALL(0U, chartInstance->c4_sfEvent, 5);
  c4_i13 = 0;
  c4_i14 = 0;
  for (c4_i15 = 0; c4_i15 < 3; c4_i15++) {
    for (c4_i16 = 0; c4_i16 < 3; c4_i16++) {
      c4_Re[c4_i16 + c4_i13] = c4_Te[c4_i16 + c4_i14];
    }

    c4_i13 += 3;
    c4_i14 += 4;
  }

  _SFD_EML_CALL(0U, chartInstance->c4_sfEvent, 15);
  c4_b_x = ((c4_Re[0] + c4_Re[4]) + c4_Re[8]) + 1.0;
  c4_c_x = c4_b_x;
  if (c4_c_x < 0.0) {
    c4_eml_error(chartInstance);
  }

  c4_c_x = muDoubleScalarSqrt(c4_c_x);
  c4_eta = 0.5 * c4_c_x;
  _SFD_EML_CALL(0U, chartInstance->c4_sfEvent, 16);
  c4_d_x = c4_Re[5] - c4_Re[7];
  c4_e_x = c4_d_x;
  c4_e_x = muDoubleScalarSign(c4_e_x);
  c4_f_x = ((c4_Re[0] - c4_Re[4]) - c4_Re[8]) + 1.0;
  c4_g_x = c4_f_x;
  if (c4_g_x < 0.0) {
    c4_eml_error(chartInstance);
  }

  c4_g_x = muDoubleScalarSqrt(c4_g_x);
  c4_h_x = c4_Re[6] - c4_Re[2];
  c4_i_x = c4_h_x;
  c4_i_x = muDoubleScalarSign(c4_i_x);
  c4_j_x = ((c4_Re[4] - c4_Re[8]) - c4_Re[0]) + 1.0;
  c4_k_x = c4_j_x;
  if (c4_k_x < 0.0) {
    c4_eml_error(chartInstance);
  }

  c4_k_x = muDoubleScalarSqrt(c4_k_x);
  c4_l_x = c4_Re[1] - c4_Re[3];
  c4_m_x = c4_l_x;
  c4_m_x = muDoubleScalarSign(c4_m_x);
  c4_n_x = ((c4_Re[8] - c4_Re[0]) - c4_Re[4]) + 1.0;
  c4_o_x = c4_n_x;
  if (c4_o_x < 0.0) {
    c4_eml_error(chartInstance);
  }

  c4_o_x = muDoubleScalarSqrt(c4_o_x);
  c4_b[0] = c4_e_x * c4_g_x;
  c4_b[1] = c4_i_x * c4_k_x;
  c4_b[2] = c4_m_x * c4_o_x;
  for (c4_i17 = 0; c4_i17 < 3; c4_i17++) {
    c4_e[c4_i17] = 0.5 * c4_b[c4_i17];
  }

  _SFD_EML_CALL(0U, chartInstance->c4_sfEvent, 19);
  c4_Q[0] = c4_eta;
  for (c4_i18 = 0; c4_i18 < 3; c4_i18++) {
    c4_Q[c4_i18 + 1] = c4_e[c4_i18];
  }

  _SFD_EML_CALL(0U, chartInstance->c4_sfEvent, 20);
  for (c4_i19 = 0; c4_i19 < 3; c4_i19++) {
    c4_x[c4_i19] = c4_pe[c4_i19];
  }

  for (c4_i20 = 0; c4_i20 < 4; c4_i20++) {
    c4_x[c4_i20 + 3] = c4_Q[c4_i20];
  }

  _SFD_EML_CALL(0U, chartInstance->c4_sfEvent, 21);
  for (c4_i21 = 0; c4_i21 < 42; c4_i21++) {
    c4_b_J[c4_i21] = c4_J[c4_i21];
  }

  for (c4_i22 = 0; c4_i22 < 7; c4_i22++) {
    c4_b_b[c4_i22] = c4_dq[c4_i22];
  }

  c4_b_eml_scalar_eg(chartInstance);
  c4_b_eml_scalar_eg(chartInstance);
  for (c4_i23 = 0; c4_i23 < 6; c4_i23++) {
    c4_xd[c4_i23] = 0.0;
  }

  for (c4_i24 = 0; c4_i24 < 6; c4_i24++) {
    c4_xd[c4_i24] = 0.0;
  }

  for (c4_i25 = 0; c4_i25 < 6; c4_i25++) {
    c4_C[c4_i25] = c4_xd[c4_i25];
  }

  for (c4_i26 = 0; c4_i26 < 6; c4_i26++) {
    c4_xd[c4_i26] = c4_C[c4_i26];
  }

  c4_threshold(chartInstance);
  for (c4_i27 = 0; c4_i27 < 6; c4_i27++) {
    c4_C[c4_i27] = c4_xd[c4_i27];
  }

  for (c4_i28 = 0; c4_i28 < 6; c4_i28++) {
    c4_xd[c4_i28] = c4_C[c4_i28];
  }

  for (c4_i29 = 0; c4_i29 < 6; c4_i29++) {
    c4_xd[c4_i29] = 0.0;
    c4_i30 = 0;
    for (c4_i31 = 0; c4_i31 < 7; c4_i31++) {
      c4_xd[c4_i29] += c4_b_J[c4_i30 + c4_i29] * c4_b_b[c4_i31];
      c4_i30 += 6;
    }
  }

  _SFD_EML_CALL(0U, chartInstance->c4_sfEvent, -21);
  _SFD_SYMBOL_SCOPE_POP();
  for (c4_i32 = 0; c4_i32 < 7; c4_i32++) {
    (*c4_p_x)[c4_i32] = c4_x[c4_i32];
  }

  for (c4_i33 = 0; c4_i33 < 6; c4_i33++) {
    (*c4_b_xd)[c4_i33] = c4_xd[c4_i33];
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 3U, chartInstance->c4_sfEvent);
}

static void initSimStructsc4_simlwrkuka_dynamic
  (SFc4_simlwrkuka_dynamicInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c4_jacobin(SFc4_simlwrkuka_dynamicInstanceStruct *chartInstance,
  real_T c4_q1, real_T c4_q2, real_T c4_q3, real_T c4_q4, real_T c4_q5, real_T
  c4_q6, real_T c4_q7, real_T c4_J[42], real_T c4_Te[16])
{
  uint32_T c4_debug_family_var_map[50];
  real_T c4_q[7];
  real_T c4_d3;
  real_T c4_d5;
  real_T c4_d0;
  real_T c4_d7;
  real_T c4_f[7];
  real_T c4_a[7];
  real_T c4_d[7];
  real_T c4_t[7];
  real_T c4_A0[16];
  real_T c4_A1[16];
  real_T c4_A2[16];
  real_T c4_A3[16];
  real_T c4_A4[16];
  real_T c4_A5[16];
  real_T c4_A6[16];
  real_T c4_A7[16];
  real_T c4_Ae[16];
  real_T c4_T1[16];
  real_T c4_T2[16];
  real_T c4_T3[16];
  real_T c4_T4[16];
  real_T c4_T5[16];
  real_T c4_T6[16];
  real_T c4_z0[3];
  real_T c4_z1[3];
  real_T c4_z2[3];
  real_T c4_z3[3];
  real_T c4_z4[3];
  real_T c4_z5[3];
  real_T c4_z6[3];
  real_T c4_p0[3];
  real_T c4_p1[3];
  real_T c4_p2[3];
  real_T c4_p3[3];
  real_T c4_p4[3];
  real_T c4_p5[3];
  real_T c4_p6[3];
  real_T c4_pe[3];
  real_T c4_nargin = 7.0;
  real_T c4_nargout = 2.0;
  int32_T c4_i34;
  static real_T c4_dv2[7] = { 1.5707963267948966, -1.5707963267948966,
    -1.5707963267948966, 1.5707963267948966, 1.5707963267948966,
    -1.5707963267948966, 0.0 };

  int32_T c4_i35;
  int32_T c4_i36;
  static real_T c4_b_a[16] = { 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.31, 1.0 };

  real_T c4_x;
  real_T c4_b_x;
  real_T c4_c_x;
  real_T c4_d_x;
  real_T c4_e_x;
  real_T c4_f_x;
  real_T c4_g_x;
  real_T c4_h_x;
  real_T c4_i_x;
  real_T c4_j_x;
  real_T c4_k_x;
  real_T c4_l_x;
  real_T c4_m_x;
  real_T c4_n_x;
  real_T c4_o_x;
  real_T c4_p_x;
  int32_T c4_i37;
  int32_T c4_i38;
  static real_T c4_dv3[4] = { 0.0, 0.0, 0.0, 1.0 };

  real_T c4_q_x;
  real_T c4_r_x;
  real_T c4_s_x;
  real_T c4_t_x;
  real_T c4_u_x;
  real_T c4_v_x;
  real_T c4_w_x;
  real_T c4_x_x;
  real_T c4_y_x;
  real_T c4_ab_x;
  real_T c4_bb_x;
  real_T c4_cb_x;
  real_T c4_db_x;
  real_T c4_eb_x;
  real_T c4_fb_x;
  real_T c4_gb_x;
  int32_T c4_i39;
  int32_T c4_i40;
  real_T c4_hb_x;
  real_T c4_ib_x;
  real_T c4_jb_x;
  real_T c4_kb_x;
  real_T c4_lb_x;
  real_T c4_mb_x;
  real_T c4_nb_x;
  real_T c4_ob_x;
  real_T c4_pb_x;
  real_T c4_qb_x;
  real_T c4_rb_x;
  real_T c4_sb_x;
  real_T c4_tb_x;
  real_T c4_ub_x;
  real_T c4_vb_x;
  real_T c4_wb_x;
  int32_T c4_i41;
  int32_T c4_i42;
  real_T c4_xb_x;
  real_T c4_yb_x;
  real_T c4_ac_x;
  real_T c4_bc_x;
  real_T c4_cc_x;
  real_T c4_dc_x;
  real_T c4_ec_x;
  real_T c4_fc_x;
  real_T c4_gc_x;
  real_T c4_hc_x;
  real_T c4_ic_x;
  real_T c4_jc_x;
  real_T c4_kc_x;
  real_T c4_lc_x;
  real_T c4_mc_x;
  real_T c4_nc_x;
  int32_T c4_i43;
  int32_T c4_i44;
  real_T c4_oc_x;
  real_T c4_pc_x;
  real_T c4_qc_x;
  real_T c4_rc_x;
  real_T c4_sc_x;
  real_T c4_tc_x;
  real_T c4_uc_x;
  real_T c4_vc_x;
  real_T c4_wc_x;
  real_T c4_xc_x;
  real_T c4_yc_x;
  real_T c4_ad_x;
  real_T c4_bd_x;
  real_T c4_cd_x;
  real_T c4_dd_x;
  real_T c4_ed_x;
  int32_T c4_i45;
  int32_T c4_i46;
  real_T c4_fd_x;
  real_T c4_gd_x;
  real_T c4_hd_x;
  real_T c4_id_x;
  real_T c4_jd_x;
  real_T c4_kd_x;
  real_T c4_ld_x;
  real_T c4_md_x;
  real_T c4_nd_x;
  real_T c4_od_x;
  real_T c4_pd_x;
  real_T c4_qd_x;
  real_T c4_rd_x;
  real_T c4_sd_x;
  real_T c4_td_x;
  real_T c4_ud_x;
  int32_T c4_i47;
  int32_T c4_i48;
  real_T c4_vd_x;
  real_T c4_wd_x;
  real_T c4_xd_x;
  real_T c4_yd_x;
  real_T c4_ae_x;
  real_T c4_be_x;
  real_T c4_ce_x;
  real_T c4_de_x;
  real_T c4_ee_x;
  real_T c4_fe_x;
  real_T c4_ge_x;
  real_T c4_he_x;
  real_T c4_ie_x;
  real_T c4_je_x;
  real_T c4_ke_x;
  real_T c4_le_x;
  int32_T c4_i49;
  int32_T c4_i50;
  int32_T c4_i51;
  static real_T c4_b[16] = { 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.078, 1.0 };

  int32_T c4_i52;
  real_T c4_b_b[16];
  int32_T c4_i53;
  int32_T c4_i54;
  int32_T c4_i55;
  real_T c4_C[16];
  int32_T c4_i56;
  int32_T c4_i57;
  int32_T c4_i58;
  int32_T c4_i59;
  int32_T c4_i60;
  int32_T c4_i61;
  int32_T c4_i62;
  int32_T c4_i63;
  int32_T c4_i64;
  int32_T c4_i65;
  int32_T c4_i66;
  int32_T c4_i67;
  real_T c4_y[16];
  int32_T c4_i68;
  int32_T c4_i69;
  int32_T c4_i70;
  int32_T c4_i71;
  int32_T c4_i72;
  int32_T c4_i73;
  int32_T c4_i74;
  int32_T c4_i75;
  int32_T c4_i76;
  int32_T c4_i77;
  int32_T c4_i78;
  int32_T c4_i79;
  int32_T c4_i80;
  int32_T c4_i81;
  int32_T c4_i82;
  int32_T c4_i83;
  int32_T c4_i84;
  int32_T c4_i85;
  int32_T c4_i86;
  int32_T c4_i87;
  int32_T c4_i88;
  int32_T c4_i89;
  int32_T c4_i90;
  int32_T c4_i91;
  real_T c4_b_y[16];
  int32_T c4_i92;
  int32_T c4_i93;
  int32_T c4_i94;
  int32_T c4_i95;
  int32_T c4_i96;
  int32_T c4_i97;
  int32_T c4_i98;
  int32_T c4_i99;
  int32_T c4_i100;
  int32_T c4_i101;
  int32_T c4_i102;
  int32_T c4_i103;
  int32_T c4_i104;
  int32_T c4_i105;
  int32_T c4_i106;
  int32_T c4_i107;
  int32_T c4_i108;
  int32_T c4_i109;
  int32_T c4_i110;
  int32_T c4_i111;
  int32_T c4_i112;
  int32_T c4_i113;
  int32_T c4_i114;
  int32_T c4_i115;
  int32_T c4_i116;
  int32_T c4_i117;
  int32_T c4_i118;
  int32_T c4_i119;
  int32_T c4_i120;
  int32_T c4_i121;
  int32_T c4_i122;
  int32_T c4_i123;
  int32_T c4_i124;
  int32_T c4_i125;
  int32_T c4_i126;
  int32_T c4_i127;
  int32_T c4_i128;
  int32_T c4_i129;
  int32_T c4_i130;
  int32_T c4_i131;
  int32_T c4_i132;
  int32_T c4_i133;
  int32_T c4_i134;
  int32_T c4_i135;
  int32_T c4_i136;
  int32_T c4_i137;
  int32_T c4_i138;
  int32_T c4_i139;
  int32_T c4_i140;
  int32_T c4_i141;
  int32_T c4_i142;
  int32_T c4_i143;
  int32_T c4_i144;
  int32_T c4_i145;
  int32_T c4_i146;
  int32_T c4_i147;
  int32_T c4_i148;
  int32_T c4_i149;
  int32_T c4_i150;
  int32_T c4_i151;
  int32_T c4_i152;
  int32_T c4_i153;
  int32_T c4_i154;
  int32_T c4_i155;
  int32_T c4_i156;
  int32_T c4_i157;
  int32_T c4_i158;
  int32_T c4_i159;
  int32_T c4_i160;
  int32_T c4_i161;
  int32_T c4_i162;
  int32_T c4_i163;
  int32_T c4_i164;
  int32_T c4_i165;
  int32_T c4_i166;
  int32_T c4_i167;
  int32_T c4_i168;
  int32_T c4_i169;
  int32_T c4_i170;
  int32_T c4_i171;
  int32_T c4_i172;
  int32_T c4_i173;
  int32_T c4_i174;
  int32_T c4_i175;
  int32_T c4_i176;
  int32_T c4_i177;
  int32_T c4_i178;
  int32_T c4_i179;
  int32_T c4_i180;
  int32_T c4_i181;
  int32_T c4_i182;
  int32_T c4_i183;
  int32_T c4_i184;
  int32_T c4_i185;
  int32_T c4_i186;
  int32_T c4_i187;
  int32_T c4_i188;
  int32_T c4_i189;
  int32_T c4_i190;
  int32_T c4_i191;
  int32_T c4_i192;
  int32_T c4_i193;
  int32_T c4_i194;
  int32_T c4_i195;
  int32_T c4_i196;
  int32_T c4_i197;
  int32_T c4_i198;
  int32_T c4_i199;
  int32_T c4_i200;
  int32_T c4_i201;
  int32_T c4_i202;
  int32_T c4_i203;
  int32_T c4_i204;
  int32_T c4_i205;
  int32_T c4_i206;
  int32_T c4_i207;
  int32_T c4_i208;
  int32_T c4_i209;
  int32_T c4_i210;
  int32_T c4_i211;
  int32_T c4_i212;
  int32_T c4_i213;
  int32_T c4_i214;
  int32_T c4_i215;
  int32_T c4_i216;
  int32_T c4_i217;
  int32_T c4_i218;
  int32_T c4_i219;
  int32_T c4_i220;
  int32_T c4_i221;
  int32_T c4_i222;
  int32_T c4_i223;
  int32_T c4_i224;
  int32_T c4_i225;
  int32_T c4_i226;
  int32_T c4_i227;
  int32_T c4_i228;
  int32_T c4_i229;
  int32_T c4_i230;
  int32_T c4_i231;
  int32_T c4_i232;
  int32_T c4_i233;
  int32_T c4_i234;
  int32_T c4_i235;
  int32_T c4_i236;
  int32_T c4_i237;
  int32_T c4_i238;
  int32_T c4_i239;
  int32_T c4_i240;
  int32_T c4_i241;
  int32_T c4_i242;
  int32_T c4_i243;
  int32_T c4_i244;
  int32_T c4_i245;
  int32_T c4_i246;
  int32_T c4_i247;
  int32_T c4_i248;
  int32_T c4_i249;
  int32_T c4_i250;
  int32_T c4_i251;
  int32_T c4_i252;
  int32_T c4_i253;
  int32_T c4_i254;
  int32_T c4_i255;
  int32_T c4_i256;
  int32_T c4_i257;
  int32_T c4_i258;
  int32_T c4_i259;
  int32_T c4_i260;
  int32_T c4_i261;
  int32_T c4_i262;
  int32_T c4_i263;
  int32_T c4_i264;
  int32_T c4_i265;
  int32_T c4_i266;
  int32_T c4_i267;
  static real_T c4_c_a[3] = { 0.0, 0.0, 1.0 };

  int32_T c4_i268;
  int32_T c4_i269;
  int32_T c4_i270;
  int32_T c4_i271;
  int32_T c4_i272;
  int32_T c4_i273;
  int32_T c4_i274;
  int32_T c4_i275;
  int32_T c4_i276;
  int32_T c4_i277;
  int32_T c4_i278;
  int32_T c4_i279;
  int32_T c4_i280;
  int32_T c4_i281;
  int32_T c4_i282;
  real_T c4_c_b[3];
  real_T c4_c1;
  real_T c4_c2;
  real_T c4_c3;
  real_T c4_dv4[3];
  int32_T c4_i283;
  real_T c4_d_a[3];
  int32_T c4_i284;
  real_T c4_b_c1;
  real_T c4_b_c2;
  real_T c4_b_c3;
  real_T c4_dv5[3];
  int32_T c4_i285;
  int32_T c4_i286;
  real_T c4_c_c1;
  real_T c4_c_c2;
  real_T c4_c_c3;
  real_T c4_dv6[3];
  int32_T c4_i287;
  int32_T c4_i288;
  real_T c4_d_c1;
  real_T c4_d_c2;
  real_T c4_d_c3;
  real_T c4_dv7[3];
  int32_T c4_i289;
  int32_T c4_i290;
  real_T c4_e_c1;
  real_T c4_e_c2;
  real_T c4_e_c3;
  real_T c4_dv8[3];
  int32_T c4_i291;
  int32_T c4_i292;
  real_T c4_f_c1;
  real_T c4_f_c2;
  real_T c4_f_c3;
  real_T c4_dv9[3];
  int32_T c4_i293;
  int32_T c4_i294;
  real_T c4_g_c1;
  real_T c4_g_c2;
  real_T c4_g_c3;
  int32_T c4_i295;
  int32_T c4_i296;
  int32_T c4_i297;
  int32_T c4_i298;
  int32_T c4_i299;
  int32_T c4_i300;
  int32_T c4_i301;
  int32_T c4_i302;
  int32_T c4_i303;
  int32_T c4_i304;
  int32_T c4_i305;
  int32_T c4_i306;
  int32_T c4_i307;
  int32_T c4_i308;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 50U, 50U, c4_b_debug_family_names,
    c4_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c4_q, 0U, c4_b_sf_marshallOut,
    c4_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c4_d3, 1U, c4_c_sf_marshallOut,
    c4_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c4_d5, 2U, c4_c_sf_marshallOut,
    c4_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c4_d0, 3U, c4_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c4_d7, 4U, c4_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c4_f, 5U, c4_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c4_a, 6U, c4_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c4_d, 7U, c4_b_sf_marshallOut,
    c4_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c4_t, 8U, c4_b_sf_marshallOut,
    c4_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c4_A0, 9U, c4_g_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c4_A1, 10U, c4_g_sf_marshallOut,
    c4_g_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c4_A2, 11U, c4_g_sf_marshallOut,
    c4_g_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c4_A3, 12U, c4_g_sf_marshallOut,
    c4_g_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c4_A4, 13U, c4_g_sf_marshallOut,
    c4_g_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c4_A5, 14U, c4_g_sf_marshallOut,
    c4_g_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c4_A6, 15U, c4_g_sf_marshallOut,
    c4_g_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c4_A7, 16U, c4_g_sf_marshallOut,
    c4_g_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c4_Ae, 17U, c4_g_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c4_T1, 18U, c4_g_sf_marshallOut,
    c4_g_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c4_T2, 19U, c4_g_sf_marshallOut,
    c4_g_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c4_T3, 20U, c4_g_sf_marshallOut,
    c4_g_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c4_T4, 21U, c4_g_sf_marshallOut,
    c4_g_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c4_T5, 22U, c4_g_sf_marshallOut,
    c4_g_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c4_T6, 23U, c4_g_sf_marshallOut,
    c4_g_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c4_z0, 24U, c4_e_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c4_z1, 25U, c4_e_sf_marshallOut,
    c4_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c4_z2, 26U, c4_e_sf_marshallOut,
    c4_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c4_z3, 27U, c4_e_sf_marshallOut,
    c4_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c4_z4, 28U, c4_e_sf_marshallOut,
    c4_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c4_z5, 29U, c4_e_sf_marshallOut,
    c4_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c4_z6, 30U, c4_e_sf_marshallOut,
    c4_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c4_p0, 31U, c4_e_sf_marshallOut,
    c4_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c4_p1, 32U, c4_e_sf_marshallOut,
    c4_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c4_p2, 33U, c4_e_sf_marshallOut,
    c4_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c4_p3, 34U, c4_e_sf_marshallOut,
    c4_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c4_p4, 35U, c4_e_sf_marshallOut,
    c4_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c4_p5, 36U, c4_e_sf_marshallOut,
    c4_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c4_p6, 37U, c4_e_sf_marshallOut,
    c4_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c4_pe, 38U, c4_e_sf_marshallOut,
    c4_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c4_nargin, 39U, c4_c_sf_marshallOut,
    c4_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c4_nargout, 40U, c4_c_sf_marshallOut,
    c4_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c4_q1, 41U, c4_c_sf_marshallOut,
    c4_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c4_q2, 42U, c4_c_sf_marshallOut,
    c4_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c4_q3, 43U, c4_c_sf_marshallOut,
    c4_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c4_q4, 44U, c4_c_sf_marshallOut,
    c4_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c4_q5, 45U, c4_c_sf_marshallOut,
    c4_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c4_q6, 46U, c4_c_sf_marshallOut,
    c4_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c4_q7, 47U, c4_c_sf_marshallOut,
    c4_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c4_J, 48U, c4_h_sf_marshallOut,
    c4_h_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c4_Te, 49U, c4_g_sf_marshallOut,
    c4_g_sf_marshallIn);
  CV_SCRIPT_FCN(0, 0);
  _SFD_SCRIPT_CALL(0U, chartInstance->c4_sfEvent, 2);
  c4_q[0] = c4_q1;
  c4_q[1] = c4_q2;
  c4_q[2] = c4_q3;
  c4_q[3] = c4_q4;
  c4_q[4] = c4_q5;
  c4_q[5] = c4_q6;
  c4_q[6] = c4_q7;
  _SFD_SCRIPT_CALL(0U, chartInstance->c4_sfEvent, 3);
  c4_d3 = 0.4;
  _SFD_SCRIPT_CALL(0U, chartInstance->c4_sfEvent, 4);
  c4_d5 = 0.39;
  _SFD_SCRIPT_CALL(0U, chartInstance->c4_sfEvent, 5);
  c4_d0 = 0.31;
  _SFD_SCRIPT_CALL(0U, chartInstance->c4_sfEvent, 6);
  c4_d7 = 0.078;
  _SFD_SCRIPT_CALL(0U, chartInstance->c4_sfEvent, 7);
  for (c4_i34 = 0; c4_i34 < 7; c4_i34++) {
    c4_f[c4_i34] = c4_dv2[c4_i34];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c4_sfEvent, 8);
  for (c4_i35 = 0; c4_i35 < 7; c4_i35++) {
    c4_a[c4_i35] = 0.0;
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c4_sfEvent, 9);
  c4_d[0] = 0.0;
  c4_d[1] = 0.0;
  c4_d[2] = c4_d3;
  c4_d[3] = 0.0;
  c4_d[4] = c4_d5;
  c4_d[5] = 0.0;
  c4_d[6] = 0.0;
  _SFD_SCRIPT_CALL(0U, chartInstance->c4_sfEvent, 10);
  c4_t[0] = c4_q1;
  c4_t[1] = c4_q2;
  c4_t[2] = c4_q3;
  c4_t[3] = c4_q4;
  c4_t[4] = c4_q5;
  c4_t[5] = c4_q6;
  c4_t[6] = c4_q7;
  _SFD_SCRIPT_CALL(0U, chartInstance->c4_sfEvent, 12);
  for (c4_i36 = 0; c4_i36 < 16; c4_i36++) {
    c4_A0[c4_i36] = c4_b_a[c4_i36];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c4_sfEvent, 13);
  c4_x = c4_t[0];
  c4_b_x = c4_x;
  c4_b_x = muDoubleScalarCos(c4_b_x);
  c4_c_x = c4_t[0];
  c4_d_x = c4_c_x;
  c4_d_x = muDoubleScalarSin(c4_d_x);
  c4_e_x = c4_t[0];
  c4_f_x = c4_e_x;
  c4_f_x = muDoubleScalarSin(c4_f_x);
  c4_g_x = c4_t[0];
  c4_h_x = c4_g_x;
  c4_h_x = muDoubleScalarCos(c4_h_x);
  c4_i_x = c4_t[0];
  c4_j_x = c4_i_x;
  c4_j_x = muDoubleScalarSin(c4_j_x);
  c4_k_x = c4_t[0];
  c4_l_x = c4_k_x;
  c4_l_x = muDoubleScalarCos(c4_l_x);
  c4_m_x = c4_t[0];
  c4_n_x = c4_m_x;
  c4_n_x = muDoubleScalarCos(c4_n_x);
  c4_o_x = c4_t[0];
  c4_p_x = c4_o_x;
  c4_p_x = muDoubleScalarSin(c4_p_x);
  c4_A1[0] = c4_b_x;
  c4_A1[4] = -c4_d_x * 6.123233995736766E-17;
  c4_A1[8] = c4_f_x;
  c4_A1[12] = 0.0 * c4_h_x;
  c4_A1[1] = c4_j_x;
  c4_A1[5] = c4_l_x * 6.123233995736766E-17;
  c4_A1[9] = -c4_n_x;
  c4_A1[13] = 0.0 * c4_p_x;
  c4_A1[2] = 0.0;
  c4_A1[6] = 1.0;
  c4_A1[10] = 6.123233995736766E-17;
  c4_A1[14] = c4_d[0];
  c4_i37 = 0;
  for (c4_i38 = 0; c4_i38 < 4; c4_i38++) {
    c4_A1[c4_i37 + 3] = c4_dv3[c4_i38];
    c4_i37 += 4;
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c4_sfEvent, 14);
  c4_q_x = c4_t[1];
  c4_r_x = c4_q_x;
  c4_r_x = muDoubleScalarCos(c4_r_x);
  c4_s_x = c4_t[1];
  c4_t_x = c4_s_x;
  c4_t_x = muDoubleScalarSin(c4_t_x);
  c4_u_x = c4_t[1];
  c4_v_x = c4_u_x;
  c4_v_x = muDoubleScalarSin(c4_v_x);
  c4_w_x = c4_t[1];
  c4_x_x = c4_w_x;
  c4_x_x = muDoubleScalarCos(c4_x_x);
  c4_y_x = c4_t[1];
  c4_ab_x = c4_y_x;
  c4_ab_x = muDoubleScalarSin(c4_ab_x);
  c4_bb_x = c4_t[1];
  c4_cb_x = c4_bb_x;
  c4_cb_x = muDoubleScalarCos(c4_cb_x);
  c4_db_x = c4_t[1];
  c4_eb_x = c4_db_x;
  c4_eb_x = muDoubleScalarCos(c4_eb_x);
  c4_fb_x = c4_t[1];
  c4_gb_x = c4_fb_x;
  c4_gb_x = muDoubleScalarSin(c4_gb_x);
  c4_A2[0] = c4_r_x;
  c4_A2[4] = -c4_t_x * 6.123233995736766E-17;
  c4_A2[8] = -c4_v_x;
  c4_A2[12] = 0.0 * c4_x_x;
  c4_A2[1] = c4_ab_x;
  c4_A2[5] = c4_cb_x * 6.123233995736766E-17;
  c4_A2[9] = -(-c4_eb_x);
  c4_A2[13] = 0.0 * c4_gb_x;
  c4_A2[2] = 0.0;
  c4_A2[6] = -1.0;
  c4_A2[10] = 6.123233995736766E-17;
  c4_A2[14] = c4_d[1];
  c4_i39 = 0;
  for (c4_i40 = 0; c4_i40 < 4; c4_i40++) {
    c4_A2[c4_i39 + 3] = c4_dv3[c4_i40];
    c4_i39 += 4;
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c4_sfEvent, 15);
  c4_hb_x = c4_t[2];
  c4_ib_x = c4_hb_x;
  c4_ib_x = muDoubleScalarCos(c4_ib_x);
  c4_jb_x = c4_t[2];
  c4_kb_x = c4_jb_x;
  c4_kb_x = muDoubleScalarSin(c4_kb_x);
  c4_lb_x = c4_t[2];
  c4_mb_x = c4_lb_x;
  c4_mb_x = muDoubleScalarSin(c4_mb_x);
  c4_nb_x = c4_t[2];
  c4_ob_x = c4_nb_x;
  c4_ob_x = muDoubleScalarCos(c4_ob_x);
  c4_pb_x = c4_t[2];
  c4_qb_x = c4_pb_x;
  c4_qb_x = muDoubleScalarSin(c4_qb_x);
  c4_rb_x = c4_t[2];
  c4_sb_x = c4_rb_x;
  c4_sb_x = muDoubleScalarCos(c4_sb_x);
  c4_tb_x = c4_t[2];
  c4_ub_x = c4_tb_x;
  c4_ub_x = muDoubleScalarCos(c4_ub_x);
  c4_vb_x = c4_t[2];
  c4_wb_x = c4_vb_x;
  c4_wb_x = muDoubleScalarSin(c4_wb_x);
  c4_A3[0] = c4_ib_x;
  c4_A3[4] = -c4_kb_x * 6.123233995736766E-17;
  c4_A3[8] = -c4_mb_x;
  c4_A3[12] = 0.0 * c4_ob_x;
  c4_A3[1] = c4_qb_x;
  c4_A3[5] = c4_sb_x * 6.123233995736766E-17;
  c4_A3[9] = -(-c4_ub_x);
  c4_A3[13] = 0.0 * c4_wb_x;
  c4_A3[2] = 0.0;
  c4_A3[6] = -1.0;
  c4_A3[10] = 6.123233995736766E-17;
  c4_A3[14] = c4_d[2];
  c4_i41 = 0;
  for (c4_i42 = 0; c4_i42 < 4; c4_i42++) {
    c4_A3[c4_i41 + 3] = c4_dv3[c4_i42];
    c4_i41 += 4;
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c4_sfEvent, 16);
  c4_xb_x = c4_t[3];
  c4_yb_x = c4_xb_x;
  c4_yb_x = muDoubleScalarCos(c4_yb_x);
  c4_ac_x = c4_t[3];
  c4_bc_x = c4_ac_x;
  c4_bc_x = muDoubleScalarSin(c4_bc_x);
  c4_cc_x = c4_t[3];
  c4_dc_x = c4_cc_x;
  c4_dc_x = muDoubleScalarSin(c4_dc_x);
  c4_ec_x = c4_t[3];
  c4_fc_x = c4_ec_x;
  c4_fc_x = muDoubleScalarCos(c4_fc_x);
  c4_gc_x = c4_t[3];
  c4_hc_x = c4_gc_x;
  c4_hc_x = muDoubleScalarSin(c4_hc_x);
  c4_ic_x = c4_t[3];
  c4_jc_x = c4_ic_x;
  c4_jc_x = muDoubleScalarCos(c4_jc_x);
  c4_kc_x = c4_t[3];
  c4_lc_x = c4_kc_x;
  c4_lc_x = muDoubleScalarCos(c4_lc_x);
  c4_mc_x = c4_t[3];
  c4_nc_x = c4_mc_x;
  c4_nc_x = muDoubleScalarSin(c4_nc_x);
  c4_A4[0] = c4_yb_x;
  c4_A4[4] = -c4_bc_x * 6.123233995736766E-17;
  c4_A4[8] = c4_dc_x;
  c4_A4[12] = 0.0 * c4_fc_x;
  c4_A4[1] = c4_hc_x;
  c4_A4[5] = c4_jc_x * 6.123233995736766E-17;
  c4_A4[9] = -c4_lc_x;
  c4_A4[13] = 0.0 * c4_nc_x;
  c4_A4[2] = 0.0;
  c4_A4[6] = 1.0;
  c4_A4[10] = 6.123233995736766E-17;
  c4_A4[14] = c4_d[3];
  c4_i43 = 0;
  for (c4_i44 = 0; c4_i44 < 4; c4_i44++) {
    c4_A4[c4_i43 + 3] = c4_dv3[c4_i44];
    c4_i43 += 4;
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c4_sfEvent, 17);
  c4_oc_x = c4_t[4];
  c4_pc_x = c4_oc_x;
  c4_pc_x = muDoubleScalarCos(c4_pc_x);
  c4_qc_x = c4_t[4];
  c4_rc_x = c4_qc_x;
  c4_rc_x = muDoubleScalarSin(c4_rc_x);
  c4_sc_x = c4_t[4];
  c4_tc_x = c4_sc_x;
  c4_tc_x = muDoubleScalarSin(c4_tc_x);
  c4_uc_x = c4_t[4];
  c4_vc_x = c4_uc_x;
  c4_vc_x = muDoubleScalarCos(c4_vc_x);
  c4_wc_x = c4_t[4];
  c4_xc_x = c4_wc_x;
  c4_xc_x = muDoubleScalarSin(c4_xc_x);
  c4_yc_x = c4_t[4];
  c4_ad_x = c4_yc_x;
  c4_ad_x = muDoubleScalarCos(c4_ad_x);
  c4_bd_x = c4_t[4];
  c4_cd_x = c4_bd_x;
  c4_cd_x = muDoubleScalarCos(c4_cd_x);
  c4_dd_x = c4_t[4];
  c4_ed_x = c4_dd_x;
  c4_ed_x = muDoubleScalarSin(c4_ed_x);
  c4_A5[0] = c4_pc_x;
  c4_A5[4] = -c4_rc_x * 6.123233995736766E-17;
  c4_A5[8] = c4_tc_x;
  c4_A5[12] = 0.0 * c4_vc_x;
  c4_A5[1] = c4_xc_x;
  c4_A5[5] = c4_ad_x * 6.123233995736766E-17;
  c4_A5[9] = -c4_cd_x;
  c4_A5[13] = 0.0 * c4_ed_x;
  c4_A5[2] = 0.0;
  c4_A5[6] = 1.0;
  c4_A5[10] = 6.123233995736766E-17;
  c4_A5[14] = c4_d[4];
  c4_i45 = 0;
  for (c4_i46 = 0; c4_i46 < 4; c4_i46++) {
    c4_A5[c4_i45 + 3] = c4_dv3[c4_i46];
    c4_i45 += 4;
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c4_sfEvent, 18);
  c4_fd_x = c4_t[5];
  c4_gd_x = c4_fd_x;
  c4_gd_x = muDoubleScalarCos(c4_gd_x);
  c4_hd_x = c4_t[5];
  c4_id_x = c4_hd_x;
  c4_id_x = muDoubleScalarSin(c4_id_x);
  c4_jd_x = c4_t[5];
  c4_kd_x = c4_jd_x;
  c4_kd_x = muDoubleScalarSin(c4_kd_x);
  c4_ld_x = c4_t[5];
  c4_md_x = c4_ld_x;
  c4_md_x = muDoubleScalarCos(c4_md_x);
  c4_nd_x = c4_t[5];
  c4_od_x = c4_nd_x;
  c4_od_x = muDoubleScalarSin(c4_od_x);
  c4_pd_x = c4_t[5];
  c4_qd_x = c4_pd_x;
  c4_qd_x = muDoubleScalarCos(c4_qd_x);
  c4_rd_x = c4_t[5];
  c4_sd_x = c4_rd_x;
  c4_sd_x = muDoubleScalarCos(c4_sd_x);
  c4_td_x = c4_t[5];
  c4_ud_x = c4_td_x;
  c4_ud_x = muDoubleScalarSin(c4_ud_x);
  c4_A6[0] = c4_gd_x;
  c4_A6[4] = -c4_id_x * 6.123233995736766E-17;
  c4_A6[8] = -c4_kd_x;
  c4_A6[12] = 0.0 * c4_md_x;
  c4_A6[1] = c4_od_x;
  c4_A6[5] = c4_qd_x * 6.123233995736766E-17;
  c4_A6[9] = -(-c4_sd_x);
  c4_A6[13] = 0.0 * c4_ud_x;
  c4_A6[2] = 0.0;
  c4_A6[6] = -1.0;
  c4_A6[10] = 6.123233995736766E-17;
  c4_A6[14] = c4_d[5];
  c4_i47 = 0;
  for (c4_i48 = 0; c4_i48 < 4; c4_i48++) {
    c4_A6[c4_i47 + 3] = c4_dv3[c4_i48];
    c4_i47 += 4;
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c4_sfEvent, 19);
  c4_vd_x = c4_t[6];
  c4_wd_x = c4_vd_x;
  c4_wd_x = muDoubleScalarCos(c4_wd_x);
  c4_xd_x = c4_t[6];
  c4_yd_x = c4_xd_x;
  c4_yd_x = muDoubleScalarSin(c4_yd_x);
  c4_ae_x = c4_t[6];
  c4_be_x = c4_ae_x;
  c4_be_x = muDoubleScalarSin(c4_be_x);
  c4_ce_x = c4_t[6];
  c4_de_x = c4_ce_x;
  c4_de_x = muDoubleScalarCos(c4_de_x);
  c4_ee_x = c4_t[6];
  c4_fe_x = c4_ee_x;
  c4_fe_x = muDoubleScalarSin(c4_fe_x);
  c4_ge_x = c4_t[6];
  c4_he_x = c4_ge_x;
  c4_he_x = muDoubleScalarCos(c4_he_x);
  c4_ie_x = c4_t[6];
  c4_je_x = c4_ie_x;
  c4_je_x = muDoubleScalarCos(c4_je_x);
  c4_ke_x = c4_t[6];
  c4_le_x = c4_ke_x;
  c4_le_x = muDoubleScalarSin(c4_le_x);
  c4_A7[0] = c4_wd_x;
  c4_A7[4] = -c4_yd_x;
  c4_A7[8] = c4_be_x * 0.0;
  c4_A7[12] = 0.0 * c4_de_x;
  c4_A7[1] = c4_fe_x;
  c4_A7[5] = c4_he_x;
  c4_A7[9] = -c4_je_x * 0.0;
  c4_A7[13] = 0.0 * c4_le_x;
  c4_A7[2] = 0.0;
  c4_A7[6] = 0.0;
  c4_A7[10] = 1.0;
  c4_A7[14] = c4_d[6];
  c4_i49 = 0;
  for (c4_i50 = 0; c4_i50 < 4; c4_i50++) {
    c4_A7[c4_i49 + 3] = c4_dv3[c4_i50];
    c4_i49 += 4;
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c4_sfEvent, 20);
  for (c4_i51 = 0; c4_i51 < 16; c4_i51++) {
    c4_Ae[c4_i51] = c4_b[c4_i51];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c4_sfEvent, 22);
  for (c4_i52 = 0; c4_i52 < 16; c4_i52++) {
    c4_b_b[c4_i52] = c4_A1[c4_i52];
  }

  c4_eml_scalar_eg(chartInstance);
  c4_eml_scalar_eg(chartInstance);
  for (c4_i53 = 0; c4_i53 < 16; c4_i53++) {
    c4_T1[c4_i53] = 0.0;
  }

  for (c4_i54 = 0; c4_i54 < 16; c4_i54++) {
    c4_T1[c4_i54] = 0.0;
  }

  for (c4_i55 = 0; c4_i55 < 16; c4_i55++) {
    c4_C[c4_i55] = c4_T1[c4_i55];
  }

  for (c4_i56 = 0; c4_i56 < 16; c4_i56++) {
    c4_T1[c4_i56] = c4_C[c4_i56];
  }

  c4_threshold(chartInstance);
  for (c4_i57 = 0; c4_i57 < 16; c4_i57++) {
    c4_C[c4_i57] = c4_T1[c4_i57];
  }

  for (c4_i58 = 0; c4_i58 < 16; c4_i58++) {
    c4_T1[c4_i58] = c4_C[c4_i58];
  }

  for (c4_i59 = 0; c4_i59 < 4; c4_i59++) {
    c4_i60 = 0;
    for (c4_i61 = 0; c4_i61 < 4; c4_i61++) {
      c4_T1[c4_i60 + c4_i59] = 0.0;
      c4_i62 = 0;
      for (c4_i63 = 0; c4_i63 < 4; c4_i63++) {
        c4_T1[c4_i60 + c4_i59] += c4_b_a[c4_i62 + c4_i59] * c4_b_b[c4_i63 +
          c4_i60];
        c4_i62 += 4;
      }

      c4_i60 += 4;
    }
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c4_sfEvent, 23);
  for (c4_i64 = 0; c4_i64 < 16; c4_i64++) {
    c4_b_b[c4_i64] = c4_A1[c4_i64];
  }

  c4_eml_scalar_eg(chartInstance);
  c4_eml_scalar_eg(chartInstance);
  c4_threshold(chartInstance);
  for (c4_i65 = 0; c4_i65 < 4; c4_i65++) {
    c4_i66 = 0;
    for (c4_i67 = 0; c4_i67 < 4; c4_i67++) {
      c4_y[c4_i66 + c4_i65] = 0.0;
      c4_i68 = 0;
      for (c4_i69 = 0; c4_i69 < 4; c4_i69++) {
        c4_y[c4_i66 + c4_i65] += c4_b_a[c4_i68 + c4_i65] * c4_b_b[c4_i69 +
          c4_i66];
        c4_i68 += 4;
      }

      c4_i66 += 4;
    }
  }

  for (c4_i70 = 0; c4_i70 < 16; c4_i70++) {
    c4_b_b[c4_i70] = c4_A2[c4_i70];
  }

  c4_eml_scalar_eg(chartInstance);
  c4_eml_scalar_eg(chartInstance);
  for (c4_i71 = 0; c4_i71 < 16; c4_i71++) {
    c4_T2[c4_i71] = 0.0;
  }

  for (c4_i72 = 0; c4_i72 < 16; c4_i72++) {
    c4_T2[c4_i72] = 0.0;
  }

  for (c4_i73 = 0; c4_i73 < 16; c4_i73++) {
    c4_C[c4_i73] = c4_T2[c4_i73];
  }

  for (c4_i74 = 0; c4_i74 < 16; c4_i74++) {
    c4_T2[c4_i74] = c4_C[c4_i74];
  }

  c4_threshold(chartInstance);
  for (c4_i75 = 0; c4_i75 < 16; c4_i75++) {
    c4_C[c4_i75] = c4_T2[c4_i75];
  }

  for (c4_i76 = 0; c4_i76 < 16; c4_i76++) {
    c4_T2[c4_i76] = c4_C[c4_i76];
  }

  for (c4_i77 = 0; c4_i77 < 4; c4_i77++) {
    c4_i78 = 0;
    for (c4_i79 = 0; c4_i79 < 4; c4_i79++) {
      c4_T2[c4_i78 + c4_i77] = 0.0;
      c4_i80 = 0;
      for (c4_i81 = 0; c4_i81 < 4; c4_i81++) {
        c4_T2[c4_i78 + c4_i77] += c4_y[c4_i80 + c4_i77] * c4_b_b[c4_i81 + c4_i78];
        c4_i80 += 4;
      }

      c4_i78 += 4;
    }
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c4_sfEvent, 24);
  for (c4_i82 = 0; c4_i82 < 16; c4_i82++) {
    c4_b_b[c4_i82] = c4_A1[c4_i82];
  }

  c4_eml_scalar_eg(chartInstance);
  c4_eml_scalar_eg(chartInstance);
  c4_threshold(chartInstance);
  for (c4_i83 = 0; c4_i83 < 4; c4_i83++) {
    c4_i84 = 0;
    for (c4_i85 = 0; c4_i85 < 4; c4_i85++) {
      c4_y[c4_i84 + c4_i83] = 0.0;
      c4_i86 = 0;
      for (c4_i87 = 0; c4_i87 < 4; c4_i87++) {
        c4_y[c4_i84 + c4_i83] += c4_b_a[c4_i86 + c4_i83] * c4_b_b[c4_i87 +
          c4_i84];
        c4_i86 += 4;
      }

      c4_i84 += 4;
    }
  }

  for (c4_i88 = 0; c4_i88 < 16; c4_i88++) {
    c4_b_b[c4_i88] = c4_A2[c4_i88];
  }

  c4_eml_scalar_eg(chartInstance);
  c4_eml_scalar_eg(chartInstance);
  c4_threshold(chartInstance);
  for (c4_i89 = 0; c4_i89 < 4; c4_i89++) {
    c4_i90 = 0;
    for (c4_i91 = 0; c4_i91 < 4; c4_i91++) {
      c4_b_y[c4_i90 + c4_i89] = 0.0;
      c4_i92 = 0;
      for (c4_i93 = 0; c4_i93 < 4; c4_i93++) {
        c4_b_y[c4_i90 + c4_i89] += c4_y[c4_i92 + c4_i89] * c4_b_b[c4_i93 +
          c4_i90];
        c4_i92 += 4;
      }

      c4_i90 += 4;
    }
  }

  for (c4_i94 = 0; c4_i94 < 16; c4_i94++) {
    c4_b_b[c4_i94] = c4_A3[c4_i94];
  }

  c4_eml_scalar_eg(chartInstance);
  c4_eml_scalar_eg(chartInstance);
  for (c4_i95 = 0; c4_i95 < 16; c4_i95++) {
    c4_T3[c4_i95] = 0.0;
  }

  for (c4_i96 = 0; c4_i96 < 16; c4_i96++) {
    c4_T3[c4_i96] = 0.0;
  }

  for (c4_i97 = 0; c4_i97 < 16; c4_i97++) {
    c4_C[c4_i97] = c4_T3[c4_i97];
  }

  for (c4_i98 = 0; c4_i98 < 16; c4_i98++) {
    c4_T3[c4_i98] = c4_C[c4_i98];
  }

  c4_threshold(chartInstance);
  for (c4_i99 = 0; c4_i99 < 16; c4_i99++) {
    c4_C[c4_i99] = c4_T3[c4_i99];
  }

  for (c4_i100 = 0; c4_i100 < 16; c4_i100++) {
    c4_T3[c4_i100] = c4_C[c4_i100];
  }

  for (c4_i101 = 0; c4_i101 < 4; c4_i101++) {
    c4_i102 = 0;
    for (c4_i103 = 0; c4_i103 < 4; c4_i103++) {
      c4_T3[c4_i102 + c4_i101] = 0.0;
      c4_i104 = 0;
      for (c4_i105 = 0; c4_i105 < 4; c4_i105++) {
        c4_T3[c4_i102 + c4_i101] += c4_b_y[c4_i104 + c4_i101] * c4_b_b[c4_i105 +
          c4_i102];
        c4_i104 += 4;
      }

      c4_i102 += 4;
    }
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c4_sfEvent, 25);
  for (c4_i106 = 0; c4_i106 < 16; c4_i106++) {
    c4_b_b[c4_i106] = c4_A1[c4_i106];
  }

  c4_eml_scalar_eg(chartInstance);
  c4_eml_scalar_eg(chartInstance);
  c4_threshold(chartInstance);
  for (c4_i107 = 0; c4_i107 < 4; c4_i107++) {
    c4_i108 = 0;
    for (c4_i109 = 0; c4_i109 < 4; c4_i109++) {
      c4_y[c4_i108 + c4_i107] = 0.0;
      c4_i110 = 0;
      for (c4_i111 = 0; c4_i111 < 4; c4_i111++) {
        c4_y[c4_i108 + c4_i107] += c4_b_a[c4_i110 + c4_i107] * c4_b_b[c4_i111 +
          c4_i108];
        c4_i110 += 4;
      }

      c4_i108 += 4;
    }
  }

  for (c4_i112 = 0; c4_i112 < 16; c4_i112++) {
    c4_b_b[c4_i112] = c4_A2[c4_i112];
  }

  c4_eml_scalar_eg(chartInstance);
  c4_eml_scalar_eg(chartInstance);
  c4_threshold(chartInstance);
  for (c4_i113 = 0; c4_i113 < 4; c4_i113++) {
    c4_i114 = 0;
    for (c4_i115 = 0; c4_i115 < 4; c4_i115++) {
      c4_b_y[c4_i114 + c4_i113] = 0.0;
      c4_i116 = 0;
      for (c4_i117 = 0; c4_i117 < 4; c4_i117++) {
        c4_b_y[c4_i114 + c4_i113] += c4_y[c4_i116 + c4_i113] * c4_b_b[c4_i117 +
          c4_i114];
        c4_i116 += 4;
      }

      c4_i114 += 4;
    }
  }

  for (c4_i118 = 0; c4_i118 < 16; c4_i118++) {
    c4_b_b[c4_i118] = c4_A3[c4_i118];
  }

  c4_eml_scalar_eg(chartInstance);
  c4_eml_scalar_eg(chartInstance);
  c4_threshold(chartInstance);
  for (c4_i119 = 0; c4_i119 < 4; c4_i119++) {
    c4_i120 = 0;
    for (c4_i121 = 0; c4_i121 < 4; c4_i121++) {
      c4_y[c4_i120 + c4_i119] = 0.0;
      c4_i122 = 0;
      for (c4_i123 = 0; c4_i123 < 4; c4_i123++) {
        c4_y[c4_i120 + c4_i119] += c4_b_y[c4_i122 + c4_i119] * c4_b_b[c4_i123 +
          c4_i120];
        c4_i122 += 4;
      }

      c4_i120 += 4;
    }
  }

  for (c4_i124 = 0; c4_i124 < 16; c4_i124++) {
    c4_b_b[c4_i124] = c4_A4[c4_i124];
  }

  c4_eml_scalar_eg(chartInstance);
  c4_eml_scalar_eg(chartInstance);
  for (c4_i125 = 0; c4_i125 < 16; c4_i125++) {
    c4_T4[c4_i125] = 0.0;
  }

  for (c4_i126 = 0; c4_i126 < 16; c4_i126++) {
    c4_T4[c4_i126] = 0.0;
  }

  for (c4_i127 = 0; c4_i127 < 16; c4_i127++) {
    c4_C[c4_i127] = c4_T4[c4_i127];
  }

  for (c4_i128 = 0; c4_i128 < 16; c4_i128++) {
    c4_T4[c4_i128] = c4_C[c4_i128];
  }

  c4_threshold(chartInstance);
  for (c4_i129 = 0; c4_i129 < 16; c4_i129++) {
    c4_C[c4_i129] = c4_T4[c4_i129];
  }

  for (c4_i130 = 0; c4_i130 < 16; c4_i130++) {
    c4_T4[c4_i130] = c4_C[c4_i130];
  }

  for (c4_i131 = 0; c4_i131 < 4; c4_i131++) {
    c4_i132 = 0;
    for (c4_i133 = 0; c4_i133 < 4; c4_i133++) {
      c4_T4[c4_i132 + c4_i131] = 0.0;
      c4_i134 = 0;
      for (c4_i135 = 0; c4_i135 < 4; c4_i135++) {
        c4_T4[c4_i132 + c4_i131] += c4_y[c4_i134 + c4_i131] * c4_b_b[c4_i135 +
          c4_i132];
        c4_i134 += 4;
      }

      c4_i132 += 4;
    }
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c4_sfEvent, 26);
  for (c4_i136 = 0; c4_i136 < 16; c4_i136++) {
    c4_b_b[c4_i136] = c4_A1[c4_i136];
  }

  c4_eml_scalar_eg(chartInstance);
  c4_eml_scalar_eg(chartInstance);
  c4_threshold(chartInstance);
  for (c4_i137 = 0; c4_i137 < 4; c4_i137++) {
    c4_i138 = 0;
    for (c4_i139 = 0; c4_i139 < 4; c4_i139++) {
      c4_y[c4_i138 + c4_i137] = 0.0;
      c4_i140 = 0;
      for (c4_i141 = 0; c4_i141 < 4; c4_i141++) {
        c4_y[c4_i138 + c4_i137] += c4_b_a[c4_i140 + c4_i137] * c4_b_b[c4_i141 +
          c4_i138];
        c4_i140 += 4;
      }

      c4_i138 += 4;
    }
  }

  for (c4_i142 = 0; c4_i142 < 16; c4_i142++) {
    c4_b_b[c4_i142] = c4_A2[c4_i142];
  }

  c4_eml_scalar_eg(chartInstance);
  c4_eml_scalar_eg(chartInstance);
  c4_threshold(chartInstance);
  for (c4_i143 = 0; c4_i143 < 4; c4_i143++) {
    c4_i144 = 0;
    for (c4_i145 = 0; c4_i145 < 4; c4_i145++) {
      c4_b_y[c4_i144 + c4_i143] = 0.0;
      c4_i146 = 0;
      for (c4_i147 = 0; c4_i147 < 4; c4_i147++) {
        c4_b_y[c4_i144 + c4_i143] += c4_y[c4_i146 + c4_i143] * c4_b_b[c4_i147 +
          c4_i144];
        c4_i146 += 4;
      }

      c4_i144 += 4;
    }
  }

  for (c4_i148 = 0; c4_i148 < 16; c4_i148++) {
    c4_b_b[c4_i148] = c4_A3[c4_i148];
  }

  c4_eml_scalar_eg(chartInstance);
  c4_eml_scalar_eg(chartInstance);
  c4_threshold(chartInstance);
  for (c4_i149 = 0; c4_i149 < 4; c4_i149++) {
    c4_i150 = 0;
    for (c4_i151 = 0; c4_i151 < 4; c4_i151++) {
      c4_y[c4_i150 + c4_i149] = 0.0;
      c4_i152 = 0;
      for (c4_i153 = 0; c4_i153 < 4; c4_i153++) {
        c4_y[c4_i150 + c4_i149] += c4_b_y[c4_i152 + c4_i149] * c4_b_b[c4_i153 +
          c4_i150];
        c4_i152 += 4;
      }

      c4_i150 += 4;
    }
  }

  for (c4_i154 = 0; c4_i154 < 16; c4_i154++) {
    c4_b_b[c4_i154] = c4_A4[c4_i154];
  }

  c4_eml_scalar_eg(chartInstance);
  c4_eml_scalar_eg(chartInstance);
  c4_threshold(chartInstance);
  for (c4_i155 = 0; c4_i155 < 4; c4_i155++) {
    c4_i156 = 0;
    for (c4_i157 = 0; c4_i157 < 4; c4_i157++) {
      c4_b_y[c4_i156 + c4_i155] = 0.0;
      c4_i158 = 0;
      for (c4_i159 = 0; c4_i159 < 4; c4_i159++) {
        c4_b_y[c4_i156 + c4_i155] += c4_y[c4_i158 + c4_i155] * c4_b_b[c4_i159 +
          c4_i156];
        c4_i158 += 4;
      }

      c4_i156 += 4;
    }
  }

  for (c4_i160 = 0; c4_i160 < 16; c4_i160++) {
    c4_b_b[c4_i160] = c4_A5[c4_i160];
  }

  c4_eml_scalar_eg(chartInstance);
  c4_eml_scalar_eg(chartInstance);
  for (c4_i161 = 0; c4_i161 < 16; c4_i161++) {
    c4_T5[c4_i161] = 0.0;
  }

  for (c4_i162 = 0; c4_i162 < 16; c4_i162++) {
    c4_T5[c4_i162] = 0.0;
  }

  for (c4_i163 = 0; c4_i163 < 16; c4_i163++) {
    c4_C[c4_i163] = c4_T5[c4_i163];
  }

  for (c4_i164 = 0; c4_i164 < 16; c4_i164++) {
    c4_T5[c4_i164] = c4_C[c4_i164];
  }

  c4_threshold(chartInstance);
  for (c4_i165 = 0; c4_i165 < 16; c4_i165++) {
    c4_C[c4_i165] = c4_T5[c4_i165];
  }

  for (c4_i166 = 0; c4_i166 < 16; c4_i166++) {
    c4_T5[c4_i166] = c4_C[c4_i166];
  }

  for (c4_i167 = 0; c4_i167 < 4; c4_i167++) {
    c4_i168 = 0;
    for (c4_i169 = 0; c4_i169 < 4; c4_i169++) {
      c4_T5[c4_i168 + c4_i167] = 0.0;
      c4_i170 = 0;
      for (c4_i171 = 0; c4_i171 < 4; c4_i171++) {
        c4_T5[c4_i168 + c4_i167] += c4_b_y[c4_i170 + c4_i167] * c4_b_b[c4_i171 +
          c4_i168];
        c4_i170 += 4;
      }

      c4_i168 += 4;
    }
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c4_sfEvent, 27);
  for (c4_i172 = 0; c4_i172 < 16; c4_i172++) {
    c4_b_b[c4_i172] = c4_A1[c4_i172];
  }

  c4_eml_scalar_eg(chartInstance);
  c4_eml_scalar_eg(chartInstance);
  c4_threshold(chartInstance);
  for (c4_i173 = 0; c4_i173 < 4; c4_i173++) {
    c4_i174 = 0;
    for (c4_i175 = 0; c4_i175 < 4; c4_i175++) {
      c4_y[c4_i174 + c4_i173] = 0.0;
      c4_i176 = 0;
      for (c4_i177 = 0; c4_i177 < 4; c4_i177++) {
        c4_y[c4_i174 + c4_i173] += c4_b_a[c4_i176 + c4_i173] * c4_b_b[c4_i177 +
          c4_i174];
        c4_i176 += 4;
      }

      c4_i174 += 4;
    }
  }

  for (c4_i178 = 0; c4_i178 < 16; c4_i178++) {
    c4_b_b[c4_i178] = c4_A2[c4_i178];
  }

  c4_eml_scalar_eg(chartInstance);
  c4_eml_scalar_eg(chartInstance);
  c4_threshold(chartInstance);
  for (c4_i179 = 0; c4_i179 < 4; c4_i179++) {
    c4_i180 = 0;
    for (c4_i181 = 0; c4_i181 < 4; c4_i181++) {
      c4_b_y[c4_i180 + c4_i179] = 0.0;
      c4_i182 = 0;
      for (c4_i183 = 0; c4_i183 < 4; c4_i183++) {
        c4_b_y[c4_i180 + c4_i179] += c4_y[c4_i182 + c4_i179] * c4_b_b[c4_i183 +
          c4_i180];
        c4_i182 += 4;
      }

      c4_i180 += 4;
    }
  }

  for (c4_i184 = 0; c4_i184 < 16; c4_i184++) {
    c4_b_b[c4_i184] = c4_A3[c4_i184];
  }

  c4_eml_scalar_eg(chartInstance);
  c4_eml_scalar_eg(chartInstance);
  c4_threshold(chartInstance);
  for (c4_i185 = 0; c4_i185 < 4; c4_i185++) {
    c4_i186 = 0;
    for (c4_i187 = 0; c4_i187 < 4; c4_i187++) {
      c4_y[c4_i186 + c4_i185] = 0.0;
      c4_i188 = 0;
      for (c4_i189 = 0; c4_i189 < 4; c4_i189++) {
        c4_y[c4_i186 + c4_i185] += c4_b_y[c4_i188 + c4_i185] * c4_b_b[c4_i189 +
          c4_i186];
        c4_i188 += 4;
      }

      c4_i186 += 4;
    }
  }

  for (c4_i190 = 0; c4_i190 < 16; c4_i190++) {
    c4_b_b[c4_i190] = c4_A4[c4_i190];
  }

  c4_eml_scalar_eg(chartInstance);
  c4_eml_scalar_eg(chartInstance);
  c4_threshold(chartInstance);
  for (c4_i191 = 0; c4_i191 < 4; c4_i191++) {
    c4_i192 = 0;
    for (c4_i193 = 0; c4_i193 < 4; c4_i193++) {
      c4_b_y[c4_i192 + c4_i191] = 0.0;
      c4_i194 = 0;
      for (c4_i195 = 0; c4_i195 < 4; c4_i195++) {
        c4_b_y[c4_i192 + c4_i191] += c4_y[c4_i194 + c4_i191] * c4_b_b[c4_i195 +
          c4_i192];
        c4_i194 += 4;
      }

      c4_i192 += 4;
    }
  }

  for (c4_i196 = 0; c4_i196 < 16; c4_i196++) {
    c4_b_b[c4_i196] = c4_A5[c4_i196];
  }

  c4_eml_scalar_eg(chartInstance);
  c4_eml_scalar_eg(chartInstance);
  c4_threshold(chartInstance);
  for (c4_i197 = 0; c4_i197 < 4; c4_i197++) {
    c4_i198 = 0;
    for (c4_i199 = 0; c4_i199 < 4; c4_i199++) {
      c4_y[c4_i198 + c4_i197] = 0.0;
      c4_i200 = 0;
      for (c4_i201 = 0; c4_i201 < 4; c4_i201++) {
        c4_y[c4_i198 + c4_i197] += c4_b_y[c4_i200 + c4_i197] * c4_b_b[c4_i201 +
          c4_i198];
        c4_i200 += 4;
      }

      c4_i198 += 4;
    }
  }

  for (c4_i202 = 0; c4_i202 < 16; c4_i202++) {
    c4_b_b[c4_i202] = c4_A6[c4_i202];
  }

  c4_eml_scalar_eg(chartInstance);
  c4_eml_scalar_eg(chartInstance);
  for (c4_i203 = 0; c4_i203 < 16; c4_i203++) {
    c4_T6[c4_i203] = 0.0;
  }

  for (c4_i204 = 0; c4_i204 < 16; c4_i204++) {
    c4_T6[c4_i204] = 0.0;
  }

  for (c4_i205 = 0; c4_i205 < 16; c4_i205++) {
    c4_C[c4_i205] = c4_T6[c4_i205];
  }

  for (c4_i206 = 0; c4_i206 < 16; c4_i206++) {
    c4_T6[c4_i206] = c4_C[c4_i206];
  }

  c4_threshold(chartInstance);
  for (c4_i207 = 0; c4_i207 < 16; c4_i207++) {
    c4_C[c4_i207] = c4_T6[c4_i207];
  }

  for (c4_i208 = 0; c4_i208 < 16; c4_i208++) {
    c4_T6[c4_i208] = c4_C[c4_i208];
  }

  for (c4_i209 = 0; c4_i209 < 4; c4_i209++) {
    c4_i210 = 0;
    for (c4_i211 = 0; c4_i211 < 4; c4_i211++) {
      c4_T6[c4_i210 + c4_i209] = 0.0;
      c4_i212 = 0;
      for (c4_i213 = 0; c4_i213 < 4; c4_i213++) {
        c4_T6[c4_i210 + c4_i209] += c4_y[c4_i212 + c4_i209] * c4_b_b[c4_i213 +
          c4_i210];
        c4_i212 += 4;
      }

      c4_i210 += 4;
    }
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c4_sfEvent, 28);
  for (c4_i214 = 0; c4_i214 < 16; c4_i214++) {
    c4_b_b[c4_i214] = c4_A1[c4_i214];
  }

  c4_eml_scalar_eg(chartInstance);
  c4_eml_scalar_eg(chartInstance);
  c4_threshold(chartInstance);
  for (c4_i215 = 0; c4_i215 < 4; c4_i215++) {
    c4_i216 = 0;
    for (c4_i217 = 0; c4_i217 < 4; c4_i217++) {
      c4_y[c4_i216 + c4_i215] = 0.0;
      c4_i218 = 0;
      for (c4_i219 = 0; c4_i219 < 4; c4_i219++) {
        c4_y[c4_i216 + c4_i215] += c4_b_a[c4_i218 + c4_i215] * c4_b_b[c4_i219 +
          c4_i216];
        c4_i218 += 4;
      }

      c4_i216 += 4;
    }
  }

  for (c4_i220 = 0; c4_i220 < 16; c4_i220++) {
    c4_b_b[c4_i220] = c4_A2[c4_i220];
  }

  c4_eml_scalar_eg(chartInstance);
  c4_eml_scalar_eg(chartInstance);
  c4_threshold(chartInstance);
  for (c4_i221 = 0; c4_i221 < 4; c4_i221++) {
    c4_i222 = 0;
    for (c4_i223 = 0; c4_i223 < 4; c4_i223++) {
      c4_b_y[c4_i222 + c4_i221] = 0.0;
      c4_i224 = 0;
      for (c4_i225 = 0; c4_i225 < 4; c4_i225++) {
        c4_b_y[c4_i222 + c4_i221] += c4_y[c4_i224 + c4_i221] * c4_b_b[c4_i225 +
          c4_i222];
        c4_i224 += 4;
      }

      c4_i222 += 4;
    }
  }

  for (c4_i226 = 0; c4_i226 < 16; c4_i226++) {
    c4_b_b[c4_i226] = c4_A3[c4_i226];
  }

  c4_eml_scalar_eg(chartInstance);
  c4_eml_scalar_eg(chartInstance);
  c4_threshold(chartInstance);
  for (c4_i227 = 0; c4_i227 < 4; c4_i227++) {
    c4_i228 = 0;
    for (c4_i229 = 0; c4_i229 < 4; c4_i229++) {
      c4_y[c4_i228 + c4_i227] = 0.0;
      c4_i230 = 0;
      for (c4_i231 = 0; c4_i231 < 4; c4_i231++) {
        c4_y[c4_i228 + c4_i227] += c4_b_y[c4_i230 + c4_i227] * c4_b_b[c4_i231 +
          c4_i228];
        c4_i230 += 4;
      }

      c4_i228 += 4;
    }
  }

  for (c4_i232 = 0; c4_i232 < 16; c4_i232++) {
    c4_b_b[c4_i232] = c4_A4[c4_i232];
  }

  c4_eml_scalar_eg(chartInstance);
  c4_eml_scalar_eg(chartInstance);
  c4_threshold(chartInstance);
  for (c4_i233 = 0; c4_i233 < 4; c4_i233++) {
    c4_i234 = 0;
    for (c4_i235 = 0; c4_i235 < 4; c4_i235++) {
      c4_b_y[c4_i234 + c4_i233] = 0.0;
      c4_i236 = 0;
      for (c4_i237 = 0; c4_i237 < 4; c4_i237++) {
        c4_b_y[c4_i234 + c4_i233] += c4_y[c4_i236 + c4_i233] * c4_b_b[c4_i237 +
          c4_i234];
        c4_i236 += 4;
      }

      c4_i234 += 4;
    }
  }

  for (c4_i238 = 0; c4_i238 < 16; c4_i238++) {
    c4_b_b[c4_i238] = c4_A5[c4_i238];
  }

  c4_eml_scalar_eg(chartInstance);
  c4_eml_scalar_eg(chartInstance);
  c4_threshold(chartInstance);
  for (c4_i239 = 0; c4_i239 < 4; c4_i239++) {
    c4_i240 = 0;
    for (c4_i241 = 0; c4_i241 < 4; c4_i241++) {
      c4_y[c4_i240 + c4_i239] = 0.0;
      c4_i242 = 0;
      for (c4_i243 = 0; c4_i243 < 4; c4_i243++) {
        c4_y[c4_i240 + c4_i239] += c4_b_y[c4_i242 + c4_i239] * c4_b_b[c4_i243 +
          c4_i240];
        c4_i242 += 4;
      }

      c4_i240 += 4;
    }
  }

  for (c4_i244 = 0; c4_i244 < 16; c4_i244++) {
    c4_b_b[c4_i244] = c4_A6[c4_i244];
  }

  c4_eml_scalar_eg(chartInstance);
  c4_eml_scalar_eg(chartInstance);
  c4_threshold(chartInstance);
  for (c4_i245 = 0; c4_i245 < 4; c4_i245++) {
    c4_i246 = 0;
    for (c4_i247 = 0; c4_i247 < 4; c4_i247++) {
      c4_b_y[c4_i246 + c4_i245] = 0.0;
      c4_i248 = 0;
      for (c4_i249 = 0; c4_i249 < 4; c4_i249++) {
        c4_b_y[c4_i246 + c4_i245] += c4_y[c4_i248 + c4_i245] * c4_b_b[c4_i249 +
          c4_i246];
        c4_i248 += 4;
      }

      c4_i246 += 4;
    }
  }

  for (c4_i250 = 0; c4_i250 < 16; c4_i250++) {
    c4_b_b[c4_i250] = c4_A7[c4_i250];
  }

  c4_eml_scalar_eg(chartInstance);
  c4_eml_scalar_eg(chartInstance);
  c4_threshold(chartInstance);
  for (c4_i251 = 0; c4_i251 < 4; c4_i251++) {
    c4_i252 = 0;
    for (c4_i253 = 0; c4_i253 < 4; c4_i253++) {
      c4_y[c4_i252 + c4_i251] = 0.0;
      c4_i254 = 0;
      for (c4_i255 = 0; c4_i255 < 4; c4_i255++) {
        c4_y[c4_i252 + c4_i251] += c4_b_y[c4_i254 + c4_i251] * c4_b_b[c4_i255 +
          c4_i252];
        c4_i254 += 4;
      }

      c4_i252 += 4;
    }
  }

  c4_eml_scalar_eg(chartInstance);
  c4_eml_scalar_eg(chartInstance);
  for (c4_i256 = 0; c4_i256 < 16; c4_i256++) {
    c4_Te[c4_i256] = 0.0;
  }

  for (c4_i257 = 0; c4_i257 < 16; c4_i257++) {
    c4_Te[c4_i257] = 0.0;
  }

  for (c4_i258 = 0; c4_i258 < 16; c4_i258++) {
    c4_C[c4_i258] = c4_Te[c4_i258];
  }

  for (c4_i259 = 0; c4_i259 < 16; c4_i259++) {
    c4_Te[c4_i259] = c4_C[c4_i259];
  }

  c4_threshold(chartInstance);
  for (c4_i260 = 0; c4_i260 < 16; c4_i260++) {
    c4_C[c4_i260] = c4_Te[c4_i260];
  }

  for (c4_i261 = 0; c4_i261 < 16; c4_i261++) {
    c4_Te[c4_i261] = c4_C[c4_i261];
  }

  for (c4_i262 = 0; c4_i262 < 4; c4_i262++) {
    c4_i263 = 0;
    for (c4_i264 = 0; c4_i264 < 4; c4_i264++) {
      c4_Te[c4_i263 + c4_i262] = 0.0;
      c4_i265 = 0;
      for (c4_i266 = 0; c4_i266 < 4; c4_i266++) {
        c4_Te[c4_i263 + c4_i262] += c4_y[c4_i265 + c4_i262] * c4_b[c4_i266 +
          c4_i263];
        c4_i265 += 4;
      }

      c4_i263 += 4;
    }
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c4_sfEvent, 39);
  for (c4_i267 = 0; c4_i267 < 3; c4_i267++) {
    c4_z0[c4_i267] = c4_c_a[c4_i267];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c4_sfEvent, 40);
  for (c4_i268 = 0; c4_i268 < 3; c4_i268++) {
    c4_z1[c4_i268] = c4_T1[c4_i268 + 8];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c4_sfEvent, 41);
  for (c4_i269 = 0; c4_i269 < 3; c4_i269++) {
    c4_z2[c4_i269] = c4_T2[c4_i269 + 8];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c4_sfEvent, 42);
  for (c4_i270 = 0; c4_i270 < 3; c4_i270++) {
    c4_z3[c4_i270] = c4_T3[c4_i270 + 8];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c4_sfEvent, 43);
  for (c4_i271 = 0; c4_i271 < 3; c4_i271++) {
    c4_z4[c4_i271] = c4_T4[c4_i271 + 8];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c4_sfEvent, 44);
  for (c4_i272 = 0; c4_i272 < 3; c4_i272++) {
    c4_z5[c4_i272] = c4_T5[c4_i272 + 8];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c4_sfEvent, 45);
  for (c4_i273 = 0; c4_i273 < 3; c4_i273++) {
    c4_z6[c4_i273] = c4_T6[c4_i273 + 8];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c4_sfEvent, 47);
  for (c4_i274 = 0; c4_i274 < 3; c4_i274++) {
    c4_p0[c4_i274] = 0.0;
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c4_sfEvent, 48);
  for (c4_i275 = 0; c4_i275 < 3; c4_i275++) {
    c4_p1[c4_i275] = c4_T1[c4_i275 + 12];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c4_sfEvent, 49);
  for (c4_i276 = 0; c4_i276 < 3; c4_i276++) {
    c4_p2[c4_i276] = c4_T2[c4_i276 + 12];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c4_sfEvent, 50);
  for (c4_i277 = 0; c4_i277 < 3; c4_i277++) {
    c4_p3[c4_i277] = c4_T3[c4_i277 + 12];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c4_sfEvent, 51);
  for (c4_i278 = 0; c4_i278 < 3; c4_i278++) {
    c4_p4[c4_i278] = c4_T4[c4_i278 + 12];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c4_sfEvent, 52);
  for (c4_i279 = 0; c4_i279 < 3; c4_i279++) {
    c4_p5[c4_i279] = c4_T5[c4_i279 + 12];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c4_sfEvent, 53);
  for (c4_i280 = 0; c4_i280 < 3; c4_i280++) {
    c4_p6[c4_i280] = c4_T6[c4_i280 + 12];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c4_sfEvent, 54);
  for (c4_i281 = 0; c4_i281 < 3; c4_i281++) {
    c4_pe[c4_i281] = c4_Te[c4_i281 + 12];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c4_sfEvent, 55);
  for (c4_i282 = 0; c4_i282 < 3; c4_i282++) {
    c4_c_b[c4_i282] = c4_pe[c4_i282] - c4_p0[c4_i282];
  }

  c4_c1 = 0.0 * c4_c_b[2] - c4_c_b[1];
  c4_c2 = c4_c_b[0] - 0.0 * c4_c_b[2];
  c4_c3 = 0.0 * c4_c_b[1] - 0.0 * c4_c_b[0];
  c4_dv4[0] = c4_c1;
  c4_dv4[1] = c4_c2;
  c4_dv4[2] = c4_c3;
  for (c4_i283 = 0; c4_i283 < 3; c4_i283++) {
    c4_d_a[c4_i283] = c4_z1[c4_i283];
  }

  for (c4_i284 = 0; c4_i284 < 3; c4_i284++) {
    c4_c_b[c4_i284] = c4_pe[c4_i284] - c4_p1[c4_i284];
  }

  c4_b_c1 = c4_d_a[1] * c4_c_b[2] - c4_d_a[2] * c4_c_b[1];
  c4_b_c2 = c4_d_a[2] * c4_c_b[0] - c4_d_a[0] * c4_c_b[2];
  c4_b_c3 = c4_d_a[0] * c4_c_b[1] - c4_d_a[1] * c4_c_b[0];
  c4_dv5[0] = c4_b_c1;
  c4_dv5[1] = c4_b_c2;
  c4_dv5[2] = c4_b_c3;
  for (c4_i285 = 0; c4_i285 < 3; c4_i285++) {
    c4_d_a[c4_i285] = c4_z2[c4_i285];
  }

  for (c4_i286 = 0; c4_i286 < 3; c4_i286++) {
    c4_c_b[c4_i286] = c4_pe[c4_i286] - c4_p2[c4_i286];
  }

  c4_c_c1 = c4_d_a[1] * c4_c_b[2] - c4_d_a[2] * c4_c_b[1];
  c4_c_c2 = c4_d_a[2] * c4_c_b[0] - c4_d_a[0] * c4_c_b[2];
  c4_c_c3 = c4_d_a[0] * c4_c_b[1] - c4_d_a[1] * c4_c_b[0];
  c4_dv6[0] = c4_c_c1;
  c4_dv6[1] = c4_c_c2;
  c4_dv6[2] = c4_c_c3;
  for (c4_i287 = 0; c4_i287 < 3; c4_i287++) {
    c4_d_a[c4_i287] = c4_z3[c4_i287];
  }

  for (c4_i288 = 0; c4_i288 < 3; c4_i288++) {
    c4_c_b[c4_i288] = c4_pe[c4_i288] - c4_p3[c4_i288];
  }

  c4_d_c1 = c4_d_a[1] * c4_c_b[2] - c4_d_a[2] * c4_c_b[1];
  c4_d_c2 = c4_d_a[2] * c4_c_b[0] - c4_d_a[0] * c4_c_b[2];
  c4_d_c3 = c4_d_a[0] * c4_c_b[1] - c4_d_a[1] * c4_c_b[0];
  c4_dv7[0] = c4_d_c1;
  c4_dv7[1] = c4_d_c2;
  c4_dv7[2] = c4_d_c3;
  for (c4_i289 = 0; c4_i289 < 3; c4_i289++) {
    c4_d_a[c4_i289] = c4_z4[c4_i289];
  }

  for (c4_i290 = 0; c4_i290 < 3; c4_i290++) {
    c4_c_b[c4_i290] = c4_pe[c4_i290] - c4_p4[c4_i290];
  }

  c4_e_c1 = c4_d_a[1] * c4_c_b[2] - c4_d_a[2] * c4_c_b[1];
  c4_e_c2 = c4_d_a[2] * c4_c_b[0] - c4_d_a[0] * c4_c_b[2];
  c4_e_c3 = c4_d_a[0] * c4_c_b[1] - c4_d_a[1] * c4_c_b[0];
  c4_dv8[0] = c4_e_c1;
  c4_dv8[1] = c4_e_c2;
  c4_dv8[2] = c4_e_c3;
  for (c4_i291 = 0; c4_i291 < 3; c4_i291++) {
    c4_d_a[c4_i291] = c4_z5[c4_i291];
  }

  for (c4_i292 = 0; c4_i292 < 3; c4_i292++) {
    c4_c_b[c4_i292] = c4_pe[c4_i292] - c4_p5[c4_i292];
  }

  c4_f_c1 = c4_d_a[1] * c4_c_b[2] - c4_d_a[2] * c4_c_b[1];
  c4_f_c2 = c4_d_a[2] * c4_c_b[0] - c4_d_a[0] * c4_c_b[2];
  c4_f_c3 = c4_d_a[0] * c4_c_b[1] - c4_d_a[1] * c4_c_b[0];
  c4_dv9[0] = c4_f_c1;
  c4_dv9[1] = c4_f_c2;
  c4_dv9[2] = c4_f_c3;
  for (c4_i293 = 0; c4_i293 < 3; c4_i293++) {
    c4_d_a[c4_i293] = c4_z6[c4_i293];
  }

  for (c4_i294 = 0; c4_i294 < 3; c4_i294++) {
    c4_c_b[c4_i294] = c4_pe[c4_i294] - c4_p6[c4_i294];
  }

  c4_g_c1 = c4_d_a[1] * c4_c_b[2] - c4_d_a[2] * c4_c_b[1];
  c4_g_c2 = c4_d_a[2] * c4_c_b[0] - c4_d_a[0] * c4_c_b[2];
  c4_g_c3 = c4_d_a[0] * c4_c_b[1] - c4_d_a[1] * c4_c_b[0];
  c4_c_b[0] = c4_g_c1;
  c4_c_b[1] = c4_g_c2;
  c4_c_b[2] = c4_g_c3;
  for (c4_i295 = 0; c4_i295 < 3; c4_i295++) {
    c4_J[c4_i295] = c4_dv4[c4_i295];
  }

  for (c4_i296 = 0; c4_i296 < 3; c4_i296++) {
    c4_J[c4_i296 + 6] = c4_dv5[c4_i296];
  }

  for (c4_i297 = 0; c4_i297 < 3; c4_i297++) {
    c4_J[c4_i297 + 12] = c4_dv6[c4_i297];
  }

  for (c4_i298 = 0; c4_i298 < 3; c4_i298++) {
    c4_J[c4_i298 + 18] = c4_dv7[c4_i298];
  }

  for (c4_i299 = 0; c4_i299 < 3; c4_i299++) {
    c4_J[c4_i299 + 24] = c4_dv8[c4_i299];
  }

  for (c4_i300 = 0; c4_i300 < 3; c4_i300++) {
    c4_J[c4_i300 + 30] = c4_dv9[c4_i300];
  }

  for (c4_i301 = 0; c4_i301 < 3; c4_i301++) {
    c4_J[c4_i301 + 36] = c4_c_b[c4_i301];
  }

  for (c4_i302 = 0; c4_i302 < 3; c4_i302++) {
    c4_J[c4_i302 + 3] = c4_z0[c4_i302];
  }

  for (c4_i303 = 0; c4_i303 < 3; c4_i303++) {
    c4_J[c4_i303 + 9] = c4_z1[c4_i303];
  }

  for (c4_i304 = 0; c4_i304 < 3; c4_i304++) {
    c4_J[c4_i304 + 15] = c4_z2[c4_i304];
  }

  for (c4_i305 = 0; c4_i305 < 3; c4_i305++) {
    c4_J[c4_i305 + 21] = c4_z3[c4_i305];
  }

  for (c4_i306 = 0; c4_i306 < 3; c4_i306++) {
    c4_J[c4_i306 + 27] = c4_z4[c4_i306];
  }

  for (c4_i307 = 0; c4_i307 < 3; c4_i307++) {
    c4_J[c4_i307 + 33] = c4_z5[c4_i307];
  }

  for (c4_i308 = 0; c4_i308 < 3; c4_i308++) {
    c4_J[c4_i308 + 39] = c4_z6[c4_i308];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c4_sfEvent, -55);
  _SFD_SYMBOL_SCOPE_POP();
}

static void init_script_number_translation(uint32_T c4_machineNumber, uint32_T
  c4_chartNumber, uint32_T c4_instanceNumber)
{
  (void)c4_machineNumber;
  _SFD_SCRIPT_TRANSLATION(c4_chartNumber, c4_instanceNumber, 0U,
    sf_debug_get_script_id(
    "C:\\Users\\Admin\\Desktop\\kuka lwr\\dynamic\\jacobin.m"));
}

static const mxArray *c4_sf_marshallOut(void *chartInstanceVoid, void *c4_inData)
{
  const mxArray *c4_mxArrayOutData = NULL;
  int32_T c4_i309;
  real_T c4_b_inData[6];
  int32_T c4_i310;
  real_T c4_u[6];
  const mxArray *c4_y = NULL;
  SFc4_simlwrkuka_dynamicInstanceStruct *chartInstance;
  chartInstance = (SFc4_simlwrkuka_dynamicInstanceStruct *)chartInstanceVoid;
  c4_mxArrayOutData = NULL;
  for (c4_i309 = 0; c4_i309 < 6; c4_i309++) {
    c4_b_inData[c4_i309] = (*(real_T (*)[6])c4_inData)[c4_i309];
  }

  for (c4_i310 = 0; c4_i310 < 6; c4_i310++) {
    c4_u[c4_i310] = c4_b_inData[c4_i310];
  }

  c4_y = NULL;
  sf_mex_assign(&c4_y, sf_mex_create("y", c4_u, 0, 0U, 1U, 0U, 1, 6), false);
  sf_mex_assign(&c4_mxArrayOutData, c4_y, false);
  return c4_mxArrayOutData;
}

static void c4_emlrt_marshallIn(SFc4_simlwrkuka_dynamicInstanceStruct
  *chartInstance, const mxArray *c4_xd, const char_T *c4_identifier, real_T
  c4_y[6])
{
  emlrtMsgIdentifier c4_thisId;
  c4_thisId.fIdentifier = c4_identifier;
  c4_thisId.fParent = NULL;
  c4_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c4_xd), &c4_thisId, c4_y);
  sf_mex_destroy(&c4_xd);
}

static void c4_b_emlrt_marshallIn(SFc4_simlwrkuka_dynamicInstanceStruct
  *chartInstance, const mxArray *c4_u, const emlrtMsgIdentifier *c4_parentId,
  real_T c4_y[6])
{
  real_T c4_dv10[6];
  int32_T c4_i311;
  (void)chartInstance;
  sf_mex_import(c4_parentId, sf_mex_dup(c4_u), c4_dv10, 1, 0, 0U, 1, 0U, 1, 6);
  for (c4_i311 = 0; c4_i311 < 6; c4_i311++) {
    c4_y[c4_i311] = c4_dv10[c4_i311];
  }

  sf_mex_destroy(&c4_u);
}

static void c4_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c4_mxArrayInData, const char_T *c4_varName, void *c4_outData)
{
  const mxArray *c4_xd;
  const char_T *c4_identifier;
  emlrtMsgIdentifier c4_thisId;
  real_T c4_y[6];
  int32_T c4_i312;
  SFc4_simlwrkuka_dynamicInstanceStruct *chartInstance;
  chartInstance = (SFc4_simlwrkuka_dynamicInstanceStruct *)chartInstanceVoid;
  c4_xd = sf_mex_dup(c4_mxArrayInData);
  c4_identifier = c4_varName;
  c4_thisId.fIdentifier = c4_identifier;
  c4_thisId.fParent = NULL;
  c4_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c4_xd), &c4_thisId, c4_y);
  sf_mex_destroy(&c4_xd);
  for (c4_i312 = 0; c4_i312 < 6; c4_i312++) {
    (*(real_T (*)[6])c4_outData)[c4_i312] = c4_y[c4_i312];
  }

  sf_mex_destroy(&c4_mxArrayInData);
}

static const mxArray *c4_b_sf_marshallOut(void *chartInstanceVoid, void
  *c4_inData)
{
  const mxArray *c4_mxArrayOutData = NULL;
  int32_T c4_i313;
  real_T c4_b_inData[7];
  int32_T c4_i314;
  real_T c4_u[7];
  const mxArray *c4_y = NULL;
  SFc4_simlwrkuka_dynamicInstanceStruct *chartInstance;
  chartInstance = (SFc4_simlwrkuka_dynamicInstanceStruct *)chartInstanceVoid;
  c4_mxArrayOutData = NULL;
  for (c4_i313 = 0; c4_i313 < 7; c4_i313++) {
    c4_b_inData[c4_i313] = (*(real_T (*)[7])c4_inData)[c4_i313];
  }

  for (c4_i314 = 0; c4_i314 < 7; c4_i314++) {
    c4_u[c4_i314] = c4_b_inData[c4_i314];
  }

  c4_y = NULL;
  sf_mex_assign(&c4_y, sf_mex_create("y", c4_u, 0, 0U, 1U, 0U, 1, 7), false);
  sf_mex_assign(&c4_mxArrayOutData, c4_y, false);
  return c4_mxArrayOutData;
}

static void c4_c_emlrt_marshallIn(SFc4_simlwrkuka_dynamicInstanceStruct
  *chartInstance, const mxArray *c4_x, const char_T *c4_identifier, real_T c4_y
  [7])
{
  emlrtMsgIdentifier c4_thisId;
  c4_thisId.fIdentifier = c4_identifier;
  c4_thisId.fParent = NULL;
  c4_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c4_x), &c4_thisId, c4_y);
  sf_mex_destroy(&c4_x);
}

static void c4_d_emlrt_marshallIn(SFc4_simlwrkuka_dynamicInstanceStruct
  *chartInstance, const mxArray *c4_u, const emlrtMsgIdentifier *c4_parentId,
  real_T c4_y[7])
{
  real_T c4_dv11[7];
  int32_T c4_i315;
  (void)chartInstance;
  sf_mex_import(c4_parentId, sf_mex_dup(c4_u), c4_dv11, 1, 0, 0U, 1, 0U, 1, 7);
  for (c4_i315 = 0; c4_i315 < 7; c4_i315++) {
    c4_y[c4_i315] = c4_dv11[c4_i315];
  }

  sf_mex_destroy(&c4_u);
}

static void c4_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c4_mxArrayInData, const char_T *c4_varName, void *c4_outData)
{
  const mxArray *c4_x;
  const char_T *c4_identifier;
  emlrtMsgIdentifier c4_thisId;
  real_T c4_y[7];
  int32_T c4_i316;
  SFc4_simlwrkuka_dynamicInstanceStruct *chartInstance;
  chartInstance = (SFc4_simlwrkuka_dynamicInstanceStruct *)chartInstanceVoid;
  c4_x = sf_mex_dup(c4_mxArrayInData);
  c4_identifier = c4_varName;
  c4_thisId.fIdentifier = c4_identifier;
  c4_thisId.fParent = NULL;
  c4_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c4_x), &c4_thisId, c4_y);
  sf_mex_destroy(&c4_x);
  for (c4_i316 = 0; c4_i316 < 7; c4_i316++) {
    (*(real_T (*)[7])c4_outData)[c4_i316] = c4_y[c4_i316];
  }

  sf_mex_destroy(&c4_mxArrayInData);
}

static const mxArray *c4_c_sf_marshallOut(void *chartInstanceVoid, void
  *c4_inData)
{
  const mxArray *c4_mxArrayOutData = NULL;
  real_T c4_u;
  const mxArray *c4_y = NULL;
  SFc4_simlwrkuka_dynamicInstanceStruct *chartInstance;
  chartInstance = (SFc4_simlwrkuka_dynamicInstanceStruct *)chartInstanceVoid;
  c4_mxArrayOutData = NULL;
  c4_u = *(real_T *)c4_inData;
  c4_y = NULL;
  sf_mex_assign(&c4_y, sf_mex_create("y", &c4_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c4_mxArrayOutData, c4_y, false);
  return c4_mxArrayOutData;
}

static real_T c4_e_emlrt_marshallIn(SFc4_simlwrkuka_dynamicInstanceStruct
  *chartInstance, const mxArray *c4_u, const emlrtMsgIdentifier *c4_parentId)
{
  real_T c4_y;
  real_T c4_d0;
  (void)chartInstance;
  sf_mex_import(c4_parentId, sf_mex_dup(c4_u), &c4_d0, 1, 0, 0U, 0, 0U, 0);
  c4_y = c4_d0;
  sf_mex_destroy(&c4_u);
  return c4_y;
}

static void c4_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c4_mxArrayInData, const char_T *c4_varName, void *c4_outData)
{
  const mxArray *c4_nargout;
  const char_T *c4_identifier;
  emlrtMsgIdentifier c4_thisId;
  real_T c4_y;
  SFc4_simlwrkuka_dynamicInstanceStruct *chartInstance;
  chartInstance = (SFc4_simlwrkuka_dynamicInstanceStruct *)chartInstanceVoid;
  c4_nargout = sf_mex_dup(c4_mxArrayInData);
  c4_identifier = c4_varName;
  c4_thisId.fIdentifier = c4_identifier;
  c4_thisId.fParent = NULL;
  c4_y = c4_e_emlrt_marshallIn(chartInstance, sf_mex_dup(c4_nargout), &c4_thisId);
  sf_mex_destroy(&c4_nargout);
  *(real_T *)c4_outData = c4_y;
  sf_mex_destroy(&c4_mxArrayInData);
}

static const mxArray *c4_d_sf_marshallOut(void *chartInstanceVoid, void
  *c4_inData)
{
  const mxArray *c4_mxArrayOutData = NULL;
  int32_T c4_i317;
  real_T c4_b_inData[4];
  int32_T c4_i318;
  real_T c4_u[4];
  const mxArray *c4_y = NULL;
  SFc4_simlwrkuka_dynamicInstanceStruct *chartInstance;
  chartInstance = (SFc4_simlwrkuka_dynamicInstanceStruct *)chartInstanceVoid;
  c4_mxArrayOutData = NULL;
  for (c4_i317 = 0; c4_i317 < 4; c4_i317++) {
    c4_b_inData[c4_i317] = (*(real_T (*)[4])c4_inData)[c4_i317];
  }

  for (c4_i318 = 0; c4_i318 < 4; c4_i318++) {
    c4_u[c4_i318] = c4_b_inData[c4_i318];
  }

  c4_y = NULL;
  sf_mex_assign(&c4_y, sf_mex_create("y", c4_u, 0, 0U, 1U, 0U, 1, 4), false);
  sf_mex_assign(&c4_mxArrayOutData, c4_y, false);
  return c4_mxArrayOutData;
}

static void c4_f_emlrt_marshallIn(SFc4_simlwrkuka_dynamicInstanceStruct
  *chartInstance, const mxArray *c4_u, const emlrtMsgIdentifier *c4_parentId,
  real_T c4_y[4])
{
  real_T c4_dv12[4];
  int32_T c4_i319;
  (void)chartInstance;
  sf_mex_import(c4_parentId, sf_mex_dup(c4_u), c4_dv12, 1, 0, 0U, 1, 0U, 1, 4);
  for (c4_i319 = 0; c4_i319 < 4; c4_i319++) {
    c4_y[c4_i319] = c4_dv12[c4_i319];
  }

  sf_mex_destroy(&c4_u);
}

static void c4_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c4_mxArrayInData, const char_T *c4_varName, void *c4_outData)
{
  const mxArray *c4_Q;
  const char_T *c4_identifier;
  emlrtMsgIdentifier c4_thisId;
  real_T c4_y[4];
  int32_T c4_i320;
  SFc4_simlwrkuka_dynamicInstanceStruct *chartInstance;
  chartInstance = (SFc4_simlwrkuka_dynamicInstanceStruct *)chartInstanceVoid;
  c4_Q = sf_mex_dup(c4_mxArrayInData);
  c4_identifier = c4_varName;
  c4_thisId.fIdentifier = c4_identifier;
  c4_thisId.fParent = NULL;
  c4_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c4_Q), &c4_thisId, c4_y);
  sf_mex_destroy(&c4_Q);
  for (c4_i320 = 0; c4_i320 < 4; c4_i320++) {
    (*(real_T (*)[4])c4_outData)[c4_i320] = c4_y[c4_i320];
  }

  sf_mex_destroy(&c4_mxArrayInData);
}

static const mxArray *c4_e_sf_marshallOut(void *chartInstanceVoid, void
  *c4_inData)
{
  const mxArray *c4_mxArrayOutData = NULL;
  int32_T c4_i321;
  real_T c4_b_inData[3];
  int32_T c4_i322;
  real_T c4_u[3];
  const mxArray *c4_y = NULL;
  SFc4_simlwrkuka_dynamicInstanceStruct *chartInstance;
  chartInstance = (SFc4_simlwrkuka_dynamicInstanceStruct *)chartInstanceVoid;
  c4_mxArrayOutData = NULL;
  for (c4_i321 = 0; c4_i321 < 3; c4_i321++) {
    c4_b_inData[c4_i321] = (*(real_T (*)[3])c4_inData)[c4_i321];
  }

  for (c4_i322 = 0; c4_i322 < 3; c4_i322++) {
    c4_u[c4_i322] = c4_b_inData[c4_i322];
  }

  c4_y = NULL;
  sf_mex_assign(&c4_y, sf_mex_create("y", c4_u, 0, 0U, 1U, 0U, 1, 3), false);
  sf_mex_assign(&c4_mxArrayOutData, c4_y, false);
  return c4_mxArrayOutData;
}

static void c4_g_emlrt_marshallIn(SFc4_simlwrkuka_dynamicInstanceStruct
  *chartInstance, const mxArray *c4_u, const emlrtMsgIdentifier *c4_parentId,
  real_T c4_y[3])
{
  real_T c4_dv13[3];
  int32_T c4_i323;
  (void)chartInstance;
  sf_mex_import(c4_parentId, sf_mex_dup(c4_u), c4_dv13, 1, 0, 0U, 1, 0U, 1, 3);
  for (c4_i323 = 0; c4_i323 < 3; c4_i323++) {
    c4_y[c4_i323] = c4_dv13[c4_i323];
  }

  sf_mex_destroy(&c4_u);
}

static void c4_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c4_mxArrayInData, const char_T *c4_varName, void *c4_outData)
{
  const mxArray *c4_e;
  const char_T *c4_identifier;
  emlrtMsgIdentifier c4_thisId;
  real_T c4_y[3];
  int32_T c4_i324;
  SFc4_simlwrkuka_dynamicInstanceStruct *chartInstance;
  chartInstance = (SFc4_simlwrkuka_dynamicInstanceStruct *)chartInstanceVoid;
  c4_e = sf_mex_dup(c4_mxArrayInData);
  c4_identifier = c4_varName;
  c4_thisId.fIdentifier = c4_identifier;
  c4_thisId.fParent = NULL;
  c4_g_emlrt_marshallIn(chartInstance, sf_mex_dup(c4_e), &c4_thisId, c4_y);
  sf_mex_destroy(&c4_e);
  for (c4_i324 = 0; c4_i324 < 3; c4_i324++) {
    (*(real_T (*)[3])c4_outData)[c4_i324] = c4_y[c4_i324];
  }

  sf_mex_destroy(&c4_mxArrayInData);
}

static const mxArray *c4_f_sf_marshallOut(void *chartInstanceVoid, void
  *c4_inData)
{
  const mxArray *c4_mxArrayOutData = NULL;
  int32_T c4_i325;
  int32_T c4_i326;
  int32_T c4_i327;
  real_T c4_b_inData[9];
  int32_T c4_i328;
  int32_T c4_i329;
  int32_T c4_i330;
  real_T c4_u[9];
  const mxArray *c4_y = NULL;
  SFc4_simlwrkuka_dynamicInstanceStruct *chartInstance;
  chartInstance = (SFc4_simlwrkuka_dynamicInstanceStruct *)chartInstanceVoid;
  c4_mxArrayOutData = NULL;
  c4_i325 = 0;
  for (c4_i326 = 0; c4_i326 < 3; c4_i326++) {
    for (c4_i327 = 0; c4_i327 < 3; c4_i327++) {
      c4_b_inData[c4_i327 + c4_i325] = (*(real_T (*)[9])c4_inData)[c4_i327 +
        c4_i325];
    }

    c4_i325 += 3;
  }

  c4_i328 = 0;
  for (c4_i329 = 0; c4_i329 < 3; c4_i329++) {
    for (c4_i330 = 0; c4_i330 < 3; c4_i330++) {
      c4_u[c4_i330 + c4_i328] = c4_b_inData[c4_i330 + c4_i328];
    }

    c4_i328 += 3;
  }

  c4_y = NULL;
  sf_mex_assign(&c4_y, sf_mex_create("y", c4_u, 0, 0U, 1U, 0U, 2, 3, 3), false);
  sf_mex_assign(&c4_mxArrayOutData, c4_y, false);
  return c4_mxArrayOutData;
}

static void c4_h_emlrt_marshallIn(SFc4_simlwrkuka_dynamicInstanceStruct
  *chartInstance, const mxArray *c4_u, const emlrtMsgIdentifier *c4_parentId,
  real_T c4_y[9])
{
  real_T c4_dv14[9];
  int32_T c4_i331;
  (void)chartInstance;
  sf_mex_import(c4_parentId, sf_mex_dup(c4_u), c4_dv14, 1, 0, 0U, 1, 0U, 2, 3, 3);
  for (c4_i331 = 0; c4_i331 < 9; c4_i331++) {
    c4_y[c4_i331] = c4_dv14[c4_i331];
  }

  sf_mex_destroy(&c4_u);
}

static void c4_f_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c4_mxArrayInData, const char_T *c4_varName, void *c4_outData)
{
  const mxArray *c4_Re;
  const char_T *c4_identifier;
  emlrtMsgIdentifier c4_thisId;
  real_T c4_y[9];
  int32_T c4_i332;
  int32_T c4_i333;
  int32_T c4_i334;
  SFc4_simlwrkuka_dynamicInstanceStruct *chartInstance;
  chartInstance = (SFc4_simlwrkuka_dynamicInstanceStruct *)chartInstanceVoid;
  c4_Re = sf_mex_dup(c4_mxArrayInData);
  c4_identifier = c4_varName;
  c4_thisId.fIdentifier = c4_identifier;
  c4_thisId.fParent = NULL;
  c4_h_emlrt_marshallIn(chartInstance, sf_mex_dup(c4_Re), &c4_thisId, c4_y);
  sf_mex_destroy(&c4_Re);
  c4_i332 = 0;
  for (c4_i333 = 0; c4_i333 < 3; c4_i333++) {
    for (c4_i334 = 0; c4_i334 < 3; c4_i334++) {
      (*(real_T (*)[9])c4_outData)[c4_i334 + c4_i332] = c4_y[c4_i334 + c4_i332];
    }

    c4_i332 += 3;
  }

  sf_mex_destroy(&c4_mxArrayInData);
}

static const mxArray *c4_g_sf_marshallOut(void *chartInstanceVoid, void
  *c4_inData)
{
  const mxArray *c4_mxArrayOutData = NULL;
  int32_T c4_i335;
  int32_T c4_i336;
  int32_T c4_i337;
  real_T c4_b_inData[16];
  int32_T c4_i338;
  int32_T c4_i339;
  int32_T c4_i340;
  real_T c4_u[16];
  const mxArray *c4_y = NULL;
  SFc4_simlwrkuka_dynamicInstanceStruct *chartInstance;
  chartInstance = (SFc4_simlwrkuka_dynamicInstanceStruct *)chartInstanceVoid;
  c4_mxArrayOutData = NULL;
  c4_i335 = 0;
  for (c4_i336 = 0; c4_i336 < 4; c4_i336++) {
    for (c4_i337 = 0; c4_i337 < 4; c4_i337++) {
      c4_b_inData[c4_i337 + c4_i335] = (*(real_T (*)[16])c4_inData)[c4_i337 +
        c4_i335];
    }

    c4_i335 += 4;
  }

  c4_i338 = 0;
  for (c4_i339 = 0; c4_i339 < 4; c4_i339++) {
    for (c4_i340 = 0; c4_i340 < 4; c4_i340++) {
      c4_u[c4_i340 + c4_i338] = c4_b_inData[c4_i340 + c4_i338];
    }

    c4_i338 += 4;
  }

  c4_y = NULL;
  sf_mex_assign(&c4_y, sf_mex_create("y", c4_u, 0, 0U, 1U, 0U, 2, 4, 4), false);
  sf_mex_assign(&c4_mxArrayOutData, c4_y, false);
  return c4_mxArrayOutData;
}

static void c4_i_emlrt_marshallIn(SFc4_simlwrkuka_dynamicInstanceStruct
  *chartInstance, const mxArray *c4_u, const emlrtMsgIdentifier *c4_parentId,
  real_T c4_y[16])
{
  real_T c4_dv15[16];
  int32_T c4_i341;
  (void)chartInstance;
  sf_mex_import(c4_parentId, sf_mex_dup(c4_u), c4_dv15, 1, 0, 0U, 1, 0U, 2, 4, 4);
  for (c4_i341 = 0; c4_i341 < 16; c4_i341++) {
    c4_y[c4_i341] = c4_dv15[c4_i341];
  }

  sf_mex_destroy(&c4_u);
}

static void c4_g_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c4_mxArrayInData, const char_T *c4_varName, void *c4_outData)
{
  const mxArray *c4_Te;
  const char_T *c4_identifier;
  emlrtMsgIdentifier c4_thisId;
  real_T c4_y[16];
  int32_T c4_i342;
  int32_T c4_i343;
  int32_T c4_i344;
  SFc4_simlwrkuka_dynamicInstanceStruct *chartInstance;
  chartInstance = (SFc4_simlwrkuka_dynamicInstanceStruct *)chartInstanceVoid;
  c4_Te = sf_mex_dup(c4_mxArrayInData);
  c4_identifier = c4_varName;
  c4_thisId.fIdentifier = c4_identifier;
  c4_thisId.fParent = NULL;
  c4_i_emlrt_marshallIn(chartInstance, sf_mex_dup(c4_Te), &c4_thisId, c4_y);
  sf_mex_destroy(&c4_Te);
  c4_i342 = 0;
  for (c4_i343 = 0; c4_i343 < 4; c4_i343++) {
    for (c4_i344 = 0; c4_i344 < 4; c4_i344++) {
      (*(real_T (*)[16])c4_outData)[c4_i344 + c4_i342] = c4_y[c4_i344 + c4_i342];
    }

    c4_i342 += 4;
  }

  sf_mex_destroy(&c4_mxArrayInData);
}

static const mxArray *c4_h_sf_marshallOut(void *chartInstanceVoid, void
  *c4_inData)
{
  const mxArray *c4_mxArrayOutData = NULL;
  int32_T c4_i345;
  int32_T c4_i346;
  int32_T c4_i347;
  real_T c4_b_inData[42];
  int32_T c4_i348;
  int32_T c4_i349;
  int32_T c4_i350;
  real_T c4_u[42];
  const mxArray *c4_y = NULL;
  SFc4_simlwrkuka_dynamicInstanceStruct *chartInstance;
  chartInstance = (SFc4_simlwrkuka_dynamicInstanceStruct *)chartInstanceVoid;
  c4_mxArrayOutData = NULL;
  c4_i345 = 0;
  for (c4_i346 = 0; c4_i346 < 7; c4_i346++) {
    for (c4_i347 = 0; c4_i347 < 6; c4_i347++) {
      c4_b_inData[c4_i347 + c4_i345] = (*(real_T (*)[42])c4_inData)[c4_i347 +
        c4_i345];
    }

    c4_i345 += 6;
  }

  c4_i348 = 0;
  for (c4_i349 = 0; c4_i349 < 7; c4_i349++) {
    for (c4_i350 = 0; c4_i350 < 6; c4_i350++) {
      c4_u[c4_i350 + c4_i348] = c4_b_inData[c4_i350 + c4_i348];
    }

    c4_i348 += 6;
  }

  c4_y = NULL;
  sf_mex_assign(&c4_y, sf_mex_create("y", c4_u, 0, 0U, 1U, 0U, 2, 6, 7), false);
  sf_mex_assign(&c4_mxArrayOutData, c4_y, false);
  return c4_mxArrayOutData;
}

static void c4_j_emlrt_marshallIn(SFc4_simlwrkuka_dynamicInstanceStruct
  *chartInstance, const mxArray *c4_u, const emlrtMsgIdentifier *c4_parentId,
  real_T c4_y[42])
{
  real_T c4_dv16[42];
  int32_T c4_i351;
  (void)chartInstance;
  sf_mex_import(c4_parentId, sf_mex_dup(c4_u), c4_dv16, 1, 0, 0U, 1, 0U, 2, 6, 7);
  for (c4_i351 = 0; c4_i351 < 42; c4_i351++) {
    c4_y[c4_i351] = c4_dv16[c4_i351];
  }

  sf_mex_destroy(&c4_u);
}

static void c4_h_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c4_mxArrayInData, const char_T *c4_varName, void *c4_outData)
{
  const mxArray *c4_J;
  const char_T *c4_identifier;
  emlrtMsgIdentifier c4_thisId;
  real_T c4_y[42];
  int32_T c4_i352;
  int32_T c4_i353;
  int32_T c4_i354;
  SFc4_simlwrkuka_dynamicInstanceStruct *chartInstance;
  chartInstance = (SFc4_simlwrkuka_dynamicInstanceStruct *)chartInstanceVoid;
  c4_J = sf_mex_dup(c4_mxArrayInData);
  c4_identifier = c4_varName;
  c4_thisId.fIdentifier = c4_identifier;
  c4_thisId.fParent = NULL;
  c4_j_emlrt_marshallIn(chartInstance, sf_mex_dup(c4_J), &c4_thisId, c4_y);
  sf_mex_destroy(&c4_J);
  c4_i352 = 0;
  for (c4_i353 = 0; c4_i353 < 7; c4_i353++) {
    for (c4_i354 = 0; c4_i354 < 6; c4_i354++) {
      (*(real_T (*)[42])c4_outData)[c4_i354 + c4_i352] = c4_y[c4_i354 + c4_i352];
    }

    c4_i352 += 6;
  }

  sf_mex_destroy(&c4_mxArrayInData);
}

const mxArray *sf_c4_simlwrkuka_dynamic_get_eml_resolved_functions_info(void)
{
  const mxArray *c4_nameCaptureInfo = NULL;
  c4_nameCaptureInfo = NULL;
  sf_mex_assign(&c4_nameCaptureInfo, sf_mex_createstruct("structure", 2, 33, 1),
                false);
  c4_info_helper(&c4_nameCaptureInfo);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c4_nameCaptureInfo);
  return c4_nameCaptureInfo;
}

static void c4_info_helper(const mxArray **c4_info)
{
  const mxArray *c4_rhs0 = NULL;
  const mxArray *c4_lhs0 = NULL;
  const mxArray *c4_rhs1 = NULL;
  const mxArray *c4_lhs1 = NULL;
  const mxArray *c4_rhs2 = NULL;
  const mxArray *c4_lhs2 = NULL;
  const mxArray *c4_rhs3 = NULL;
  const mxArray *c4_lhs3 = NULL;
  const mxArray *c4_rhs4 = NULL;
  const mxArray *c4_lhs4 = NULL;
  const mxArray *c4_rhs5 = NULL;
  const mxArray *c4_lhs5 = NULL;
  const mxArray *c4_rhs6 = NULL;
  const mxArray *c4_lhs6 = NULL;
  const mxArray *c4_rhs7 = NULL;
  const mxArray *c4_lhs7 = NULL;
  const mxArray *c4_rhs8 = NULL;
  const mxArray *c4_lhs8 = NULL;
  const mxArray *c4_rhs9 = NULL;
  const mxArray *c4_lhs9 = NULL;
  const mxArray *c4_rhs10 = NULL;
  const mxArray *c4_lhs10 = NULL;
  const mxArray *c4_rhs11 = NULL;
  const mxArray *c4_lhs11 = NULL;
  const mxArray *c4_rhs12 = NULL;
  const mxArray *c4_lhs12 = NULL;
  const mxArray *c4_rhs13 = NULL;
  const mxArray *c4_lhs13 = NULL;
  const mxArray *c4_rhs14 = NULL;
  const mxArray *c4_lhs14 = NULL;
  const mxArray *c4_rhs15 = NULL;
  const mxArray *c4_lhs15 = NULL;
  const mxArray *c4_rhs16 = NULL;
  const mxArray *c4_lhs16 = NULL;
  const mxArray *c4_rhs17 = NULL;
  const mxArray *c4_lhs17 = NULL;
  const mxArray *c4_rhs18 = NULL;
  const mxArray *c4_lhs18 = NULL;
  const mxArray *c4_rhs19 = NULL;
  const mxArray *c4_lhs19 = NULL;
  const mxArray *c4_rhs20 = NULL;
  const mxArray *c4_lhs20 = NULL;
  const mxArray *c4_rhs21 = NULL;
  const mxArray *c4_lhs21 = NULL;
  const mxArray *c4_rhs22 = NULL;
  const mxArray *c4_lhs22 = NULL;
  const mxArray *c4_rhs23 = NULL;
  const mxArray *c4_lhs23 = NULL;
  const mxArray *c4_rhs24 = NULL;
  const mxArray *c4_lhs24 = NULL;
  const mxArray *c4_rhs25 = NULL;
  const mxArray *c4_lhs25 = NULL;
  const mxArray *c4_rhs26 = NULL;
  const mxArray *c4_lhs26 = NULL;
  const mxArray *c4_rhs27 = NULL;
  const mxArray *c4_lhs27 = NULL;
  const mxArray *c4_rhs28 = NULL;
  const mxArray *c4_lhs28 = NULL;
  const mxArray *c4_rhs29 = NULL;
  const mxArray *c4_lhs29 = NULL;
  const mxArray *c4_rhs30 = NULL;
  const mxArray *c4_lhs30 = NULL;
  const mxArray *c4_rhs31 = NULL;
  const mxArray *c4_lhs31 = NULL;
  const mxArray *c4_rhs32 = NULL;
  const mxArray *c4_lhs32 = NULL;
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "context", "context", 0);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("jacobin"), "name", "name", 0);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 0);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[E]C:/Users/Admin/Desktop/kuka lwr/dynamic/jacobin.m"), "resolved",
                  "resolved", 0);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1500456527U), "fileTimeLo",
                  "fileTimeLo", 0);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 0);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 0);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 0);
  sf_mex_assign(&c4_rhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs0), "rhs", "rhs", 0);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs0), "lhs", "lhs", 0);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[E]C:/Users/Admin/Desktop/kuka lwr/dynamic/jacobin.m"), "context",
                  "context", 1);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("mrdivide"), "name", "name", 1);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 1);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "resolved",
                  "resolved", 1);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1388463696U), "fileTimeLo",
                  "fileTimeLo", 1);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 1);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1370017086U), "mFileTimeLo",
                  "mFileTimeLo", 1);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 1);
  sf_mex_assign(&c4_rhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs1), "rhs", "rhs", 1);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs1), "lhs", "lhs", 1);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 2);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.assert"),
                  "name", "name", 2);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 2);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/assert.m"),
                  "resolved", "resolved", 2);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 2);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 2);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 2);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 2);
  sf_mex_assign(&c4_rhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs2), "rhs", "rhs", 2);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs2), "lhs", "lhs", 2);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 3);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("rdivide"), "name", "name", 3);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 3);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "resolved",
                  "resolved", 3);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363717480U), "fileTimeLo",
                  "fileTimeLo", 3);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 3);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 3);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 3);
  sf_mex_assign(&c4_rhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs3), "rhs", "rhs", 3);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs3), "lhs", "lhs", 3);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 4);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 4);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 4);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 4);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 4);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 4);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 4);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 4);
  sf_mex_assign(&c4_rhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs4), "rhs", "rhs", 4);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs4), "lhs", "lhs", 4);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 5);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_scalexp_compatible"),
                  "name", "name", 5);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 5);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_compatible.m"),
                  "resolved", "resolved", 5);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1286825996U), "fileTimeLo",
                  "fileTimeLo", 5);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 5);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 5);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 5);
  sf_mex_assign(&c4_rhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs5), "rhs", "rhs", 5);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs5), "lhs", "lhs", 5);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 6);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_div"), "name", "name", 6);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 6);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 6);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 6);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 6);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 6);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 6);
  sf_mex_assign(&c4_rhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs6), "rhs", "rhs", 6);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs6), "lhs", "lhs", 6);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "context",
                  "context", 7);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.div"), "name",
                  "name", 7);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 7);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/div.p"), "resolved",
                  "resolved", 7);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 7);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 7);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 7);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 7);
  sf_mex_assign(&c4_rhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs7), "rhs", "rhs", 7);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs7), "lhs", "lhs", 7);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[E]C:/Users/Admin/Desktop/kuka lwr/dynamic/jacobin.m"), "context",
                  "context", 8);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("cos"), "name", "name", 8);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 8);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "resolved",
                  "resolved", 8);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1343837572U), "fileTimeLo",
                  "fileTimeLo", 8);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 8);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 8);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 8);
  sf_mex_assign(&c4_rhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs8), "rhs", "rhs", 8);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs8), "lhs", "lhs", 8);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "context",
                  "context", 9);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_scalar_cos"), "name",
                  "name", 9);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 9);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_cos.m"),
                  "resolved", "resolved", 9);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1286825922U), "fileTimeLo",
                  "fileTimeLo", 9);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 9);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 9);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 9);
  sf_mex_assign(&c4_rhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs9), "rhs", "rhs", 9);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs9), "lhs", "lhs", 9);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[E]C:/Users/Admin/Desktop/kuka lwr/dynamic/jacobin.m"), "context",
                  "context", 10);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("sin"), "name", "name", 10);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 10);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "resolved",
                  "resolved", 10);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1343837586U), "fileTimeLo",
                  "fileTimeLo", 10);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 10);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 10);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 10);
  sf_mex_assign(&c4_rhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs10), "rhs", "rhs",
                  10);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs10), "lhs", "lhs",
                  10);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "context",
                  "context", 11);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_scalar_sin"), "name",
                  "name", 11);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 11);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sin.m"),
                  "resolved", "resolved", 11);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1286825936U), "fileTimeLo",
                  "fileTimeLo", 11);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 11);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 11);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 11);
  sf_mex_assign(&c4_rhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs11), "rhs", "rhs",
                  11);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs11), "lhs", "lhs",
                  11);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[E]C:/Users/Admin/Desktop/kuka lwr/dynamic/jacobin.m"), "context",
                  "context", 12);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 12);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 12);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 12);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1383880894U), "fileTimeLo",
                  "fileTimeLo", 12);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 12);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 12);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 12);
  sf_mex_assign(&c4_rhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs12), "rhs", "rhs",
                  12);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs12), "lhs", "lhs",
                  12);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m!common_checks"),
                  "context", "context", 13);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 13);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 13);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 13);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 13);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 13);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 13);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 13);
  sf_mex_assign(&c4_rhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs13), "rhs", "rhs",
                  13);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs13), "lhs", "lhs",
                  13);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 14);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 14);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 14);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 14);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 14);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 14);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 14);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 14);
  sf_mex_assign(&c4_rhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs14), "rhs", "rhs",
                  14);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs14), "lhs", "lhs",
                  14);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 15);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 15);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 15);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 15);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 15);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 15);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 15);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 15);
  sf_mex_assign(&c4_rhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs15), "rhs", "rhs",
                  15);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs15), "lhs", "lhs",
                  15);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "context",
                  "context", 16);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 16);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 16);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 16);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 16);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 16);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 16);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 16);
  sf_mex_assign(&c4_rhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs16), "rhs", "rhs",
                  16);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs16), "lhs", "lhs",
                  16);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 17);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_xgemm"), "name", "name",
                  17);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 17);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"),
                  "resolved", "resolved", 17);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1375987890U), "fileTimeLo",
                  "fileTimeLo", 17);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 17);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 17);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 17);
  sf_mex_assign(&c4_rhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs17), "rhs", "rhs",
                  17);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs17), "lhs", "lhs",
                  17);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 18);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 18);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 18);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 18);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 18);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 18);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 18);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 18);
  sf_mex_assign(&c4_rhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs18), "rhs", "rhs",
                  18);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs18), "lhs", "lhs",
                  18);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 19);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.blas.xgemm"),
                  "name", "name", 19);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 19);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "resolved", "resolved", 19);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 19);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 19);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 19);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 19);
  sf_mex_assign(&c4_rhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs19), "rhs", "rhs",
                  19);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs19), "lhs", "lhs",
                  19);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 20);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 20);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 20);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 20);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 20);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 20);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 20);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 20);
  sf_mex_assign(&c4_rhs20, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs20, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs20), "rhs", "rhs",
                  20);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs20), "lhs", "lhs",
                  20);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p!below_threshold"),
                  "context", "context", 21);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.blas.threshold"),
                  "name", "name", 21);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 21);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 21);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 21);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 21);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 21);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 21);
  sf_mex_assign(&c4_rhs21, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs21, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs21), "rhs", "rhs",
                  21);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs21), "lhs", "lhs",
                  21);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "context", "context", 22);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 22);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 22);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 22);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1381857500U), "fileTimeLo",
                  "fileTimeLo", 22);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 22);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 22);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 22);
  sf_mex_assign(&c4_rhs22, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs22, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs22), "rhs", "rhs",
                  22);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs22), "lhs", "lhs",
                  22);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 23);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 23);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 23);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 23);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 23);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 23);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 23);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 23);
  sf_mex_assign(&c4_rhs23, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs23, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs23), "rhs", "rhs",
                  23);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs23), "lhs", "lhs",
                  23);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 24);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.refblas.xgemm"),
                  "name", "name", 24);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 24);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgemm.p"),
                  "resolved", "resolved", 24);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 24);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 24);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 24);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 24);
  sf_mex_assign(&c4_rhs24, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs24, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs24), "rhs", "rhs",
                  24);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs24), "lhs", "lhs",
                  24);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[E]C:/Users/Admin/Desktop/kuka lwr/dynamic/jacobin.m"), "context",
                  "context", 25);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("cross"), "name", "name", 25);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 25);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/specfun/cross.m"), "resolved",
                  "resolved", 25);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1286826042U), "fileTimeLo",
                  "fileTimeLo", 25);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 25);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 25);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 25);
  sf_mex_assign(&c4_rhs25, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs25, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs25), "rhs", "rhs",
                  25);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs25), "lhs", "lhs",
                  25);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "context", "context", 26);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("sqrt"), "name", "name", 26);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 26);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "resolved",
                  "resolved", 26);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1343837586U), "fileTimeLo",
                  "fileTimeLo", 26);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 26);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 26);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 26);
  sf_mex_assign(&c4_rhs26, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs26, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs26), "rhs", "rhs",
                  26);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs26), "lhs", "lhs",
                  26);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "context",
                  "context", 27);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_error"), "name", "name",
                  27);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 27);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_error.m"), "resolved",
                  "resolved", 27);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1343837558U), "fileTimeLo",
                  "fileTimeLo", 27);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 27);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 27);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 27);
  sf_mex_assign(&c4_rhs27, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs27, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs27), "rhs", "rhs",
                  27);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs27), "lhs", "lhs",
                  27);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "context",
                  "context", 28);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_scalar_sqrt"), "name",
                  "name", 28);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 28);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sqrt.m"),
                  "resolved", "resolved", 28);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1286825938U), "fileTimeLo",
                  "fileTimeLo", 28);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 28);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 28);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 28);
  sf_mex_assign(&c4_rhs28, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs28, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs28), "rhs", "rhs",
                  28);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs28), "lhs", "lhs",
                  28);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "context", "context", 29);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("sign"), "name", "name", 29);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 29);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sign.m"), "resolved",
                  "resolved", 29);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363717456U), "fileTimeLo",
                  "fileTimeLo", 29);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 29);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 29);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 29);
  sf_mex_assign(&c4_rhs29, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs29, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs29), "rhs", "rhs",
                  29);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs29), "lhs", "lhs",
                  29);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sign.m"), "context",
                  "context", 30);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 30);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 30);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 30);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 30);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 30);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 30);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 30);
  sf_mex_assign(&c4_rhs30, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs30, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs30), "rhs", "rhs",
                  30);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs30), "lhs", "lhs",
                  30);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sign.m"), "context",
                  "context", 31);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_scalar_sign"), "name",
                  "name", 31);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 31);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sign.m"),
                  "resolved", "resolved", 31);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1356545094U), "fileTimeLo",
                  "fileTimeLo", 31);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 31);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 31);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 31);
  sf_mex_assign(&c4_rhs31, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs31, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs31), "rhs", "rhs",
                  31);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs31), "lhs", "lhs",
                  31);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "context", "context", 32);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 32);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 32);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 32);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1383880894U), "fileTimeLo",
                  "fileTimeLo", 32);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 32);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 32);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 32);
  sf_mex_assign(&c4_rhs32, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs32, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs32), "rhs", "rhs",
                  32);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs32), "lhs", "lhs",
                  32);
  sf_mex_destroy(&c4_rhs0);
  sf_mex_destroy(&c4_lhs0);
  sf_mex_destroy(&c4_rhs1);
  sf_mex_destroy(&c4_lhs1);
  sf_mex_destroy(&c4_rhs2);
  sf_mex_destroy(&c4_lhs2);
  sf_mex_destroy(&c4_rhs3);
  sf_mex_destroy(&c4_lhs3);
  sf_mex_destroy(&c4_rhs4);
  sf_mex_destroy(&c4_lhs4);
  sf_mex_destroy(&c4_rhs5);
  sf_mex_destroy(&c4_lhs5);
  sf_mex_destroy(&c4_rhs6);
  sf_mex_destroy(&c4_lhs6);
  sf_mex_destroy(&c4_rhs7);
  sf_mex_destroy(&c4_lhs7);
  sf_mex_destroy(&c4_rhs8);
  sf_mex_destroy(&c4_lhs8);
  sf_mex_destroy(&c4_rhs9);
  sf_mex_destroy(&c4_lhs9);
  sf_mex_destroy(&c4_rhs10);
  sf_mex_destroy(&c4_lhs10);
  sf_mex_destroy(&c4_rhs11);
  sf_mex_destroy(&c4_lhs11);
  sf_mex_destroy(&c4_rhs12);
  sf_mex_destroy(&c4_lhs12);
  sf_mex_destroy(&c4_rhs13);
  sf_mex_destroy(&c4_lhs13);
  sf_mex_destroy(&c4_rhs14);
  sf_mex_destroy(&c4_lhs14);
  sf_mex_destroy(&c4_rhs15);
  sf_mex_destroy(&c4_lhs15);
  sf_mex_destroy(&c4_rhs16);
  sf_mex_destroy(&c4_lhs16);
  sf_mex_destroy(&c4_rhs17);
  sf_mex_destroy(&c4_lhs17);
  sf_mex_destroy(&c4_rhs18);
  sf_mex_destroy(&c4_lhs18);
  sf_mex_destroy(&c4_rhs19);
  sf_mex_destroy(&c4_lhs19);
  sf_mex_destroy(&c4_rhs20);
  sf_mex_destroy(&c4_lhs20);
  sf_mex_destroy(&c4_rhs21);
  sf_mex_destroy(&c4_lhs21);
  sf_mex_destroy(&c4_rhs22);
  sf_mex_destroy(&c4_lhs22);
  sf_mex_destroy(&c4_rhs23);
  sf_mex_destroy(&c4_lhs23);
  sf_mex_destroy(&c4_rhs24);
  sf_mex_destroy(&c4_lhs24);
  sf_mex_destroy(&c4_rhs25);
  sf_mex_destroy(&c4_lhs25);
  sf_mex_destroy(&c4_rhs26);
  sf_mex_destroy(&c4_lhs26);
  sf_mex_destroy(&c4_rhs27);
  sf_mex_destroy(&c4_lhs27);
  sf_mex_destroy(&c4_rhs28);
  sf_mex_destroy(&c4_lhs28);
  sf_mex_destroy(&c4_rhs29);
  sf_mex_destroy(&c4_lhs29);
  sf_mex_destroy(&c4_rhs30);
  sf_mex_destroy(&c4_lhs30);
  sf_mex_destroy(&c4_rhs31);
  sf_mex_destroy(&c4_lhs31);
  sf_mex_destroy(&c4_rhs32);
  sf_mex_destroy(&c4_lhs32);
}

static const mxArray *c4_emlrt_marshallOut(const char * c4_u)
{
  const mxArray *c4_y = NULL;
  c4_y = NULL;
  sf_mex_assign(&c4_y, sf_mex_create("y", c4_u, 15, 0U, 0U, 0U, 2, 1, strlen
    (c4_u)), false);
  return c4_y;
}

static const mxArray *c4_b_emlrt_marshallOut(const uint32_T c4_u)
{
  const mxArray *c4_y = NULL;
  c4_y = NULL;
  sf_mex_assign(&c4_y, sf_mex_create("y", &c4_u, 7, 0U, 0U, 0U, 0), false);
  return c4_y;
}

static void c4_eml_scalar_eg(SFc4_simlwrkuka_dynamicInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void c4_threshold(SFc4_simlwrkuka_dynamicInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c4_eml_error(SFc4_simlwrkuka_dynamicInstanceStruct *chartInstance)
{
  int32_T c4_i355;
  static char_T c4_cv0[30] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o', 'l',
    'b', 'o', 'x', ':', 'E', 'l', 'F', 'u', 'n', 'D', 'o', 'm', 'a', 'i', 'n',
    'E', 'r', 'r', 'o', 'r' };

  char_T c4_u[30];
  const mxArray *c4_y = NULL;
  int32_T c4_i356;
  static char_T c4_cv1[4] = { 's', 'q', 'r', 't' };

  char_T c4_b_u[4];
  const mxArray *c4_b_y = NULL;
  (void)chartInstance;
  for (c4_i355 = 0; c4_i355 < 30; c4_i355++) {
    c4_u[c4_i355] = c4_cv0[c4_i355];
  }

  c4_y = NULL;
  sf_mex_assign(&c4_y, sf_mex_create("y", c4_u, 10, 0U, 1U, 0U, 2, 1, 30), false);
  for (c4_i356 = 0; c4_i356 < 4; c4_i356++) {
    c4_b_u[c4_i356] = c4_cv1[c4_i356];
  }

  c4_b_y = NULL;
  sf_mex_assign(&c4_b_y, sf_mex_create("y", c4_b_u, 10, 0U, 1U, 0U, 2, 1, 4),
                false);
  sf_mex_call_debug(sfGlobalDebugInstanceStruct, "error", 0U, 1U, 14,
                    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message", 1U,
    2U, 14, c4_y, 14, c4_b_y));
}

static void c4_b_eml_scalar_eg(SFc4_simlwrkuka_dynamicInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static const mxArray *c4_i_sf_marshallOut(void *chartInstanceVoid, void
  *c4_inData)
{
  const mxArray *c4_mxArrayOutData = NULL;
  int32_T c4_u;
  const mxArray *c4_y = NULL;
  SFc4_simlwrkuka_dynamicInstanceStruct *chartInstance;
  chartInstance = (SFc4_simlwrkuka_dynamicInstanceStruct *)chartInstanceVoid;
  c4_mxArrayOutData = NULL;
  c4_u = *(int32_T *)c4_inData;
  c4_y = NULL;
  sf_mex_assign(&c4_y, sf_mex_create("y", &c4_u, 6, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c4_mxArrayOutData, c4_y, false);
  return c4_mxArrayOutData;
}

static int32_T c4_k_emlrt_marshallIn(SFc4_simlwrkuka_dynamicInstanceStruct
  *chartInstance, const mxArray *c4_u, const emlrtMsgIdentifier *c4_parentId)
{
  int32_T c4_y;
  int32_T c4_i357;
  (void)chartInstance;
  sf_mex_import(c4_parentId, sf_mex_dup(c4_u), &c4_i357, 1, 6, 0U, 0, 0U, 0);
  c4_y = c4_i357;
  sf_mex_destroy(&c4_u);
  return c4_y;
}

static void c4_i_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c4_mxArrayInData, const char_T *c4_varName, void *c4_outData)
{
  const mxArray *c4_b_sfEvent;
  const char_T *c4_identifier;
  emlrtMsgIdentifier c4_thisId;
  int32_T c4_y;
  SFc4_simlwrkuka_dynamicInstanceStruct *chartInstance;
  chartInstance = (SFc4_simlwrkuka_dynamicInstanceStruct *)chartInstanceVoid;
  c4_b_sfEvent = sf_mex_dup(c4_mxArrayInData);
  c4_identifier = c4_varName;
  c4_thisId.fIdentifier = c4_identifier;
  c4_thisId.fParent = NULL;
  c4_y = c4_k_emlrt_marshallIn(chartInstance, sf_mex_dup(c4_b_sfEvent),
    &c4_thisId);
  sf_mex_destroy(&c4_b_sfEvent);
  *(int32_T *)c4_outData = c4_y;
  sf_mex_destroy(&c4_mxArrayInData);
}

static uint8_T c4_l_emlrt_marshallIn(SFc4_simlwrkuka_dynamicInstanceStruct
  *chartInstance, const mxArray *c4_b_is_active_c4_simlwrkuka_dynamic, const
  char_T *c4_identifier)
{
  uint8_T c4_y;
  emlrtMsgIdentifier c4_thisId;
  c4_thisId.fIdentifier = c4_identifier;
  c4_thisId.fParent = NULL;
  c4_y = c4_m_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c4_b_is_active_c4_simlwrkuka_dynamic), &c4_thisId);
  sf_mex_destroy(&c4_b_is_active_c4_simlwrkuka_dynamic);
  return c4_y;
}

static uint8_T c4_m_emlrt_marshallIn(SFc4_simlwrkuka_dynamicInstanceStruct
  *chartInstance, const mxArray *c4_u, const emlrtMsgIdentifier *c4_parentId)
{
  uint8_T c4_y;
  uint8_T c4_u0;
  (void)chartInstance;
  sf_mex_import(c4_parentId, sf_mex_dup(c4_u), &c4_u0, 1, 3, 0U, 0, 0U, 0);
  c4_y = c4_u0;
  sf_mex_destroy(&c4_u);
  return c4_y;
}

static void init_dsm_address_info(SFc4_simlwrkuka_dynamicInstanceStruct
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

void sf_c4_simlwrkuka_dynamic_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(2790790916U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(1196804898U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(730481837U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(4106425310U);
}

mxArray *sf_c4_simlwrkuka_dynamic_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1,1,5,
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("OsY4dXMxbVOSvQ8EqUlzfH");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,2,3,dataFields);

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
      pr[0] = (double)(7);
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
    mxSetField(mxAutoinheritanceInfo,0,"inputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"parameters",mxCreateDoubleMatrix(0,0,
                mxREAL));
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,2,3,dataFields);

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
      pr[0] = (double)(6);
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
    mxSetField(mxAutoinheritanceInfo,0,"outputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"locals",mxCreateDoubleMatrix(0,0,mxREAL));
  }

  return(mxAutoinheritanceInfo);
}

mxArray *sf_c4_simlwrkuka_dynamic_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c4_simlwrkuka_dynamic_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

static const mxArray *sf_get_sim_state_info_c4_simlwrkuka_dynamic(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x3'type','srcId','name','auxInfo'{{M[1],M[5],T\"x\",},{M[1],M[8],T\"xd\",},{M[8],M[0],T\"is_active_c4_simlwrkuka_dynamic\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 3, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c4_simlwrkuka_dynamic_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc4_simlwrkuka_dynamicInstanceStruct *chartInstance;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
    chartInstance = (SFc4_simlwrkuka_dynamicInstanceStruct *)
      chartInfo->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _simlwrkuka_dynamicMachineNumber_,
           4,
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
        init_script_number_translation(_simlwrkuka_dynamicMachineNumber_,
          chartInstance->chartNumber,chartInstance->instanceNumber);
        if (chartAlreadyPresent==0) {
          /* this is the first instance */
          sf_debug_set_chart_disable_implicit_casting
            (sfGlobalDebugInstanceStruct,_simlwrkuka_dynamicMachineNumber_,
             chartInstance->chartNumber,1);
          sf_debug_set_chart_event_thresholds(sfGlobalDebugInstanceStruct,
            _simlwrkuka_dynamicMachineNumber_,
            chartInstance->chartNumber,
            0,
            0,
            0);
          _SFD_SET_DATA_PROPS(0,1,1,0,"q");
          _SFD_SET_DATA_PROPS(1,2,0,1,"x");
          _SFD_SET_DATA_PROPS(2,2,0,1,"xd");
          _SFD_SET_DATA_PROPS(3,1,1,0,"dq");
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
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,741);
        _SFD_CV_INIT_SCRIPT(0,1,0,0,0,0,0,0,0,0);
        _SFD_CV_INIT_SCRIPT_FCN(0,0,"jacobin",0,-1,2551);

        {
          unsigned int dimVector[1];
          dimVector[0]= 7;
          _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c4_b_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 7;
          _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c4_b_sf_marshallOut,(MexInFcnForType)
            c4_b_sf_marshallIn);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 6;
          _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c4_sf_marshallOut,(MexInFcnForType)
            c4_sf_marshallIn);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 7;
          _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c4_b_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          real_T (*c4_q)[7];
          real_T (*c4_x)[7];
          real_T (*c4_xd)[6];
          real_T (*c4_dq)[7];
          c4_dq = (real_T (*)[7])ssGetInputPortSignal(chartInstance->S, 1);
          c4_xd = (real_T (*)[6])ssGetOutputPortSignal(chartInstance->S, 2);
          c4_x = (real_T (*)[7])ssGetOutputPortSignal(chartInstance->S, 1);
          c4_q = (real_T (*)[7])ssGetInputPortSignal(chartInstance->S, 0);
          _SFD_SET_DATA_VALUE_PTR(0U, *c4_q);
          _SFD_SET_DATA_VALUE_PTR(1U, *c4_x);
          _SFD_SET_DATA_VALUE_PTR(2U, *c4_xd);
          _SFD_SET_DATA_VALUE_PTR(3U, *c4_dq);
        }
      }
    } else {
      sf_debug_reset_current_state_configuration(sfGlobalDebugInstanceStruct,
        _simlwrkuka_dynamicMachineNumber_,chartInstance->chartNumber,
        chartInstance->instanceNumber);
    }
  }
}

static const char* sf_get_instance_specialization(void)
{
  return "Z9YQKMHUe0Q3vUwzQIUyz";
}

static void sf_opaque_initialize_c4_simlwrkuka_dynamic(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc4_simlwrkuka_dynamicInstanceStruct*)
    chartInstanceVar)->S,0);
  initialize_params_c4_simlwrkuka_dynamic((SFc4_simlwrkuka_dynamicInstanceStruct*)
    chartInstanceVar);
  initialize_c4_simlwrkuka_dynamic((SFc4_simlwrkuka_dynamicInstanceStruct*)
    chartInstanceVar);
}

static void sf_opaque_enable_c4_simlwrkuka_dynamic(void *chartInstanceVar)
{
  enable_c4_simlwrkuka_dynamic((SFc4_simlwrkuka_dynamicInstanceStruct*)
    chartInstanceVar);
}

static void sf_opaque_disable_c4_simlwrkuka_dynamic(void *chartInstanceVar)
{
  disable_c4_simlwrkuka_dynamic((SFc4_simlwrkuka_dynamicInstanceStruct*)
    chartInstanceVar);
}

static void sf_opaque_gateway_c4_simlwrkuka_dynamic(void *chartInstanceVar)
{
  sf_gateway_c4_simlwrkuka_dynamic((SFc4_simlwrkuka_dynamicInstanceStruct*)
    chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c4_simlwrkuka_dynamic(SimStruct*
  S)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c4_simlwrkuka_dynamic
    ((SFc4_simlwrkuka_dynamicInstanceStruct*)chartInfo->chartInstance);/* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c4_simlwrkuka_dynamic();/* state var info */
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

extern void sf_internal_set_sim_state_c4_simlwrkuka_dynamic(SimStruct* S, const
  mxArray *st)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[3];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxDuplicateArray(st);      /* high level simctx */
  prhs[2] = (mxArray*) sf_get_sim_state_info_c4_simlwrkuka_dynamic();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 3, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c4_simlwrkuka_dynamic((SFc4_simlwrkuka_dynamicInstanceStruct*)
    chartInfo->chartInstance, mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c4_simlwrkuka_dynamic(SimStruct* S)
{
  return sf_internal_get_sim_state_c4_simlwrkuka_dynamic(S);
}

static void sf_opaque_set_sim_state_c4_simlwrkuka_dynamic(SimStruct* S, const
  mxArray *st)
{
  sf_internal_set_sim_state_c4_simlwrkuka_dynamic(S, st);
}

static void sf_opaque_terminate_c4_simlwrkuka_dynamic(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc4_simlwrkuka_dynamicInstanceStruct*) chartInstanceVar)
      ->S;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_simlwrkuka_dynamic_optimization_info();
    }

    finalize_c4_simlwrkuka_dynamic((SFc4_simlwrkuka_dynamicInstanceStruct*)
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
  initSimStructsc4_simlwrkuka_dynamic((SFc4_simlwrkuka_dynamicInstanceStruct*)
    chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c4_simlwrkuka_dynamic(SimStruct *S)
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
    initialize_params_c4_simlwrkuka_dynamic
      ((SFc4_simlwrkuka_dynamicInstanceStruct*)(chartInfo->chartInstance));
  }
}

static void mdlSetWorkWidths_c4_simlwrkuka_dynamic(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_simlwrkuka_dynamic_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(sf_get_instance_specialization(),infoStruct,4);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(sf_get_instance_specialization(),
                infoStruct,4,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop
      (sf_get_instance_specialization(),infoStruct,4,
       "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(sf_get_instance_specialization(),infoStruct,4);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,4,2);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,4,2);
    }

    {
      unsigned int outPortIdx;
      for (outPortIdx=1; outPortIdx<=2; ++outPortIdx) {
        ssSetOutputPortOptimizeInIR(S, outPortIdx, 1U);
      }
    }

    {
      unsigned int inPortIdx;
      for (inPortIdx=0; inPortIdx < 2; ++inPortIdx) {
        ssSetInputPortOptimizeInIR(S, inPortIdx, 1U);
      }
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,4);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(92822046U));
  ssSetChecksum1(S,(1701338294U));
  ssSetChecksum2(S,(3484333648U));
  ssSetChecksum3(S,(1091223799U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c4_simlwrkuka_dynamic(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c4_simlwrkuka_dynamic(SimStruct *S)
{
  SFc4_simlwrkuka_dynamicInstanceStruct *chartInstance;
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)utMalloc(sizeof
    (ChartRunTimeInfo));
  chartInstance = (SFc4_simlwrkuka_dynamicInstanceStruct *)utMalloc(sizeof
    (SFc4_simlwrkuka_dynamicInstanceStruct));
  memset(chartInstance, 0, sizeof(SFc4_simlwrkuka_dynamicInstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway =
    sf_opaque_gateway_c4_simlwrkuka_dynamic;
  chartInstance->chartInfo.initializeChart =
    sf_opaque_initialize_c4_simlwrkuka_dynamic;
  chartInstance->chartInfo.terminateChart =
    sf_opaque_terminate_c4_simlwrkuka_dynamic;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c4_simlwrkuka_dynamic;
  chartInstance->chartInfo.disableChart =
    sf_opaque_disable_c4_simlwrkuka_dynamic;
  chartInstance->chartInfo.getSimState =
    sf_opaque_get_sim_state_c4_simlwrkuka_dynamic;
  chartInstance->chartInfo.setSimState =
    sf_opaque_set_sim_state_c4_simlwrkuka_dynamic;
  chartInstance->chartInfo.getSimStateInfo =
    sf_get_sim_state_info_c4_simlwrkuka_dynamic;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c4_simlwrkuka_dynamic;
  chartInstance->chartInfo.mdlStart = mdlStart_c4_simlwrkuka_dynamic;
  chartInstance->chartInfo.mdlSetWorkWidths =
    mdlSetWorkWidths_c4_simlwrkuka_dynamic;
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

void c4_simlwrkuka_dynamic_method_dispatcher(SimStruct *S, int_T method, void
  *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c4_simlwrkuka_dynamic(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c4_simlwrkuka_dynamic(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c4_simlwrkuka_dynamic(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c4_simlwrkuka_dynamic_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
