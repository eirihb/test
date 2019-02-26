/*
 * p2p4.c
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "p2p4".
 *
 * Model version              : 1.175
 * Simulink Coder version : 8.9 (R2015b) 13-Aug-2015
 * C source code generated on : Wed Feb 20 11:14:13 2019
 *
 * Target selection: quarc_win64.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "p2p4.h"
#include "p2p4_private.h"
#include "p2p4_dt.h"

/* Block signals (auto storage) */
B_p2p4_T p2p4_B;

/* Continuous states */
X_p2p4_T p2p4_X;

/* Block states (auto storage) */
DW_p2p4_T p2p4_DW;

/* Real-time model */
RT_MODEL_p2p4_T p2p4_M_;
RT_MODEL_p2p4_T *const p2p4_M = &p2p4_M_;

/*
 * Writes out MAT-file header.  Returns success or failure.
 * Returns:
 *      0 - success
 *      1 - failure
 */
int_T rt_WriteMat4FileHeader(FILE *fp, int32_T m, int32_T n, const char *name)
{
  typedef enum { ELITTLE_ENDIAN, EBIG_ENDIAN } ByteOrder;

  int16_T one = 1;
  ByteOrder byteOrder = (*((int8_T *)&one)==1) ? ELITTLE_ENDIAN : EBIG_ENDIAN;
  int32_T type = (byteOrder == ELITTLE_ENDIAN) ? 0: 1000;
  int32_T imagf = 0;
  int32_T name_len = (int32_T)strlen(name) + 1;
  if ((fwrite(&type, sizeof(int32_T), 1, fp) == 0) ||
      (fwrite(&m, sizeof(int32_T), 1, fp) == 0) ||
      (fwrite(&n, sizeof(int32_T), 1, fp) == 0) ||
      (fwrite(&imagf, sizeof(int32_T), 1, fp) == 0) ||
      (fwrite(&name_len, sizeof(int32_T), 1, fp) == 0) ||
      (fwrite(name, sizeof(char), name_len, fp) == 0)) {
    return(1);
  } else {
    return(0);
  }
}                                      /* end rt_WriteMat4FileHeader */

/*
 * This function updates continuous states using the ODE1 fixed-step
 * solver algorithm
 */
static void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
{
  time_T tnew = rtsiGetSolverStopTime(si);
  time_T h = rtsiGetStepSize(si);
  real_T *x = rtsiGetContStates(si);
  ODE1_IntgData *id = (ODE1_IntgData *)rtsiGetSolverData(si);
  real_T *f0 = id->f[0];
  int_T i;
  int_T nXc = 4;
  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);
  rtsiSetdX(si, f0);
  p2p4_derivatives();
  rtsiSetT(si, tnew);
  for (i = 0; i < nXc; ++i) {
    x[i] += h * f0[i];
  }

  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/* Model output function */
void p2p4_output(void)
{
  /* local block i/o variables */
  real_T rtb_Derivative;
  real_T rtb_HILReadEncoderTimebase_o1;
  real_T rtb_HILReadEncoderTimebase_o2;
  real_T rtb_HILReadEncoderTimebase_o3;
  real_T rtb_TmpSignalConversionAtToFile[5];
  real_T *lastU;
  real_T rtb_Backgain;
  if (rtmIsMajorTimeStep(p2p4_M)) {
    /* set solver stop time */
    if (!(p2p4_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&p2p4_M->solverInfo, ((p2p4_M->Timing.clockTickH0 +
        1) * p2p4_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&p2p4_M->solverInfo, ((p2p4_M->Timing.clockTick0 + 1)
        * p2p4_M->Timing.stepSize0 + p2p4_M->Timing.clockTickH0 *
        p2p4_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(p2p4_M)) {
    p2p4_M->Timing.t[0] = rtsiGetT(&p2p4_M->solverInfo);
  }

  if (rtmIsMajorTimeStep(p2p4_M)) {
    /* S-Function (hil_read_encoder_timebase_block): '<S4>/HIL Read Encoder Timebase' */

    /* S-Function Block: p2p4/Helicopter_interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
    {
      t_error result;
      result = hil_task_read_encoder(p2p4_DW.HILReadEncoderTimebase_Task, 1,
        &p2p4_DW.HILReadEncoderTimebase_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(p2p4_M, _rt_error_message);
      } else {
        rtb_HILReadEncoderTimebase_o1 = p2p4_DW.HILReadEncoderTimebase_Buffer[0];
        rtb_HILReadEncoderTimebase_o2 = p2p4_DW.HILReadEncoderTimebase_Buffer[1];
        rtb_HILReadEncoderTimebase_o3 = p2p4_DW.HILReadEncoderTimebase_Buffer[2];
      }
    }
  }

  /* FromWorkspace: '<Root>/From Workspace' */
  {
    real_T *pDataValues = (real_T *) p2p4_DW.FromWorkspace_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *) p2p4_DW.FromWorkspace_PWORK.TimePtr;
    int_T currTimeIndex = p2p4_DW.FromWorkspace_IWORK.PrevIndex;
    real_T t = p2p4_M->Timing.t[0];

    /* Get index */
    if (t <= pTimeValues[0]) {
      currTimeIndex = 0;
    } else if (t >= pTimeValues[140]) {
      currTimeIndex = 139;
    } else {
      if (t < pTimeValues[currTimeIndex]) {
        while (t < pTimeValues[currTimeIndex]) {
          currTimeIndex--;
        }
      } else {
        while (t >= pTimeValues[currTimeIndex + 1]) {
          currTimeIndex++;
        }
      }
    }

    p2p4_DW.FromWorkspace_IWORK.PrevIndex = currTimeIndex;

    /* Post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          rtb_Derivative = pDataValues[currTimeIndex];
        } else {
          rtb_Derivative = pDataValues[currTimeIndex + 1];
        }
      } else {
        real_T f1 = (t2 - t) / (t2 - t1);
        real_T f2 = 1.0 - f1;
        real_T d1;
        real_T d2;
        int_T TimeIndex= currTimeIndex;
        d1 = pDataValues[TimeIndex];
        d2 = pDataValues[TimeIndex + 1];
        rtb_Derivative = (real_T) rtInterpolate(d1, d2, f1, f2);
        pDataValues += 141;
      }
    }
  }

  if (rtmIsMajorTimeStep(p2p4_M)) {
    /* Gain: '<S4>/Travel: Count to rad' incorporates:
     *  Gain: '<S4>/Travel_gain'
     */
    p2p4_B.TravelCounttorad = p2p4_P.travel_gain * rtb_HILReadEncoderTimebase_o1
      * p2p4_P.TravelCounttorad_Gain;

    /* Gain: '<S11>/Gain' */
    p2p4_B.Gain = p2p4_P.Gain_Gain * p2p4_B.TravelCounttorad;

    /* Gain: '<S4>/Pitch: Count to rad' */
    p2p4_B.PitchCounttorad = p2p4_P.PitchCounttorad_Gain *
      rtb_HILReadEncoderTimebase_o2;

    /* Gain: '<S8>/Gain' */
    p2p4_B.Gain_i = p2p4_P.Gain_Gain_a * p2p4_B.PitchCounttorad;
  }

  /* Gain: '<S12>/Gain' incorporates:
   *  TransferFcn: '<S4>/Travel: Transfer Fcn'
   */
  p2p4_B.Gain_d = (p2p4_P.TravelTransferFcn_C * p2p4_X.TravelTransferFcn_CSTATE
                   + p2p4_P.TravelTransferFcn_D * p2p4_B.TravelCounttorad) *
    p2p4_P.Gain_Gain_l;

  /* Gain: '<S9>/Gain' incorporates:
   *  TransferFcn: '<S4>/Pitch: Transfer Fcn'
   */
  p2p4_B.Gain_b = (p2p4_P.PitchTransferFcn_C * p2p4_X.PitchTransferFcn_CSTATE +
                   p2p4_P.PitchTransferFcn_D * p2p4_B.PitchCounttorad) *
    p2p4_P.Gain_Gain_ae;
  if (rtmIsMajorTimeStep(p2p4_M)) {
    /* Gain: '<S4>/Elevation: Count to rad' */
    p2p4_B.ElevationCounttorad = p2p4_P.ElevationCounttorad_Gain *
      rtb_HILReadEncoderTimebase_o3;

    /* Gain: '<S6>/Gain' */
    p2p4_B.Gain_e = p2p4_P.Gain_Gain_lv * p2p4_B.ElevationCounttorad;

    /* Sum: '<Root>/Sum' incorporates:
     *  Constant: '<Root>/elavation_offset [deg]'
     */
    p2p4_B.Sum = p2p4_B.Gain_e + p2p4_P.elavation_offsetdeg_Value;
  }

  /* Gain: '<S7>/Gain' incorporates:
   *  TransferFcn: '<S4>/Elevation: Transfer Fcn'
   */
  p2p4_B.Gain_dg = (p2p4_P.ElevationTransferFcn_C *
                    p2p4_X.ElevationTransferFcn_CSTATE +
                    p2p4_P.ElevationTransferFcn_D * p2p4_B.ElevationCounttorad) *
    p2p4_P.Gain_Gain_n;

  /* Gain: '<S2>/Gain1' */
  p2p4_B.Gain1[0] = p2p4_P.Gain1_Gain * p2p4_B.Gain;
  p2p4_B.Gain1[1] = p2p4_P.Gain1_Gain * p2p4_B.Gain_d;
  p2p4_B.Gain1[2] = p2p4_P.Gain1_Gain * p2p4_B.Gain_i;
  p2p4_B.Gain1[3] = p2p4_P.Gain1_Gain * p2p4_B.Gain_b;
  p2p4_B.Gain1[4] = p2p4_P.Gain1_Gain * p2p4_B.Sum;
  p2p4_B.Gain1[5] = p2p4_P.Gain1_Gain * p2p4_B.Gain_dg;

  /* Sum: '<Root>/Sum1' incorporates:
   *  Constant: '<Root>/Vd_bias'
   *  Gain: '<S5>/K_pd'
   *  Gain: '<S5>/K_pp'
   *  Sum: '<S5>/Sum2'
   *  Sum: '<S5>/Sum3'
   */
  p2p4_B.Sum1 = ((rtb_Derivative - p2p4_B.Gain1[2]) * p2p4_P.K_pp - p2p4_P.K_pd *
                 p2p4_B.Gain1[3]) + p2p4_P.Vd_ff;
  if (rtmIsMajorTimeStep(p2p4_M)) {
    /* SignalConversion: '<Root>/TmpSignal ConversionAtTo FileInport1' */
    rtb_TmpSignalConversionAtToFile[0] = p2p4_B.Sum1;
    rtb_TmpSignalConversionAtToFile[1] = p2p4_B.Gain1[0];
    rtb_TmpSignalConversionAtToFile[2] = p2p4_B.Gain1[1];
    rtb_TmpSignalConversionAtToFile[3] = p2p4_B.Gain1[2];
    rtb_TmpSignalConversionAtToFile[4] = p2p4_B.Gain1[3];

    /* ToFile: '<Root>/To File' */
    {
      if (!(++p2p4_DW.ToFile_IWORK.Decimation % 1) &&
          (p2p4_DW.ToFile_IWORK.Count*6)+1 < 100000000 ) {
        FILE *fp = (FILE *) p2p4_DW.ToFile_PWORK.FilePtr;
        if (fp != (NULL)) {
          real_T u[6];
          p2p4_DW.ToFile_IWORK.Decimation = 0;
          u[0] = p2p4_M->Timing.t[1];
          u[1] = rtb_TmpSignalConversionAtToFile[0];
          u[2] = rtb_TmpSignalConversionAtToFile[1];
          u[3] = rtb_TmpSignalConversionAtToFile[2];
          u[4] = rtb_TmpSignalConversionAtToFile[3];
          u[5] = rtb_TmpSignalConversionAtToFile[4];
          if (fwrite(u, sizeof(real_T), 6, fp) != 6) {
            rtmSetErrorStatus(p2p4_M, "Error writing to MAT-file p2p4_1.mat");
            return;
          }

          if (((++p2p4_DW.ToFile_IWORK.Count)*6)+1 >= 100000000) {
            (void)fprintf(stdout,
                          "*** The ToFile block will stop logging data before\n"
                          "    the simulation has ended, because it has reached\n"
                          "    the maximum number of elements (100000000)\n"
                          "    allowed in MAT-file p2p4_1.mat.\n");
          }
        }
      }
    }
  }

  /* Integrator: '<S3>/Integrator' */
  /* Limited  Integrator  */
  if (p2p4_X.Integrator_CSTATE >= p2p4_P.Integrator_UpperSat) {
    p2p4_X.Integrator_CSTATE = p2p4_P.Integrator_UpperSat;
  } else {
    if (p2p4_X.Integrator_CSTATE <= p2p4_P.Integrator_LowerSat) {
      p2p4_X.Integrator_CSTATE = p2p4_P.Integrator_LowerSat;
    }
  }

  /* Sum: '<S3>/Sum' incorporates:
   *  Constant: '<Root>/elevation_ref'
   */
  rtb_Derivative = p2p4_P.elevation_ref_Value - p2p4_B.Gain1[4];

  /* Sum: '<Root>/Sum2' incorporates:
   *  Constant: '<Root>/Vs_bias'
   *  Gain: '<S3>/K_ed'
   *  Gain: '<S3>/K_ep'
   *  Integrator: '<S3>/Integrator'
   *  Sum: '<S3>/Sum1'
   */
  p2p4_B.Sum2 = ((p2p4_P.K_ep * rtb_Derivative + p2p4_X.Integrator_CSTATE) -
                 p2p4_P.K_ed * p2p4_B.Gain1[5]) + p2p4_P.Vs_ff;
  if (rtmIsMajorTimeStep(p2p4_M)) {
  }

  /* Gain: '<S3>/K_ei' */
  p2p4_B.K_ei = p2p4_P.K_ei * rtb_Derivative;
  if (rtmIsMajorTimeStep(p2p4_M)) {
  }

  /* Derivative: '<S4>/Derivative' */
  if ((p2p4_DW.TimeStampA >= p2p4_M->Timing.t[0]) && (p2p4_DW.TimeStampB >=
       p2p4_M->Timing.t[0])) {
    rtb_Derivative = 0.0;
  } else {
    rtb_Backgain = p2p4_DW.TimeStampA;
    lastU = &p2p4_DW.LastUAtTimeA;
    if (p2p4_DW.TimeStampA < p2p4_DW.TimeStampB) {
      if (p2p4_DW.TimeStampB < p2p4_M->Timing.t[0]) {
        rtb_Backgain = p2p4_DW.TimeStampB;
        lastU = &p2p4_DW.LastUAtTimeB;
      }
    } else {
      if (p2p4_DW.TimeStampA >= p2p4_M->Timing.t[0]) {
        rtb_Backgain = p2p4_DW.TimeStampB;
        lastU = &p2p4_DW.LastUAtTimeB;
      }
    }

    rtb_Derivative = (p2p4_B.PitchCounttorad - *lastU) / (p2p4_M->Timing.t[0] -
      rtb_Backgain);
  }

  /* End of Derivative: '<S4>/Derivative' */

  /* Gain: '<S10>/Gain' */
  p2p4_B.Gain_l = p2p4_P.Gain_Gain_a1 * rtb_Derivative;
  if (rtmIsMajorTimeStep(p2p4_M)) {
  }

  /* Gain: '<S1>/Back gain' incorporates:
   *  Sum: '<S1>/Subtract'
   */
  rtb_Backgain = (p2p4_B.Sum2 - p2p4_B.Sum1) * p2p4_P.Backgain_Gain;

  /* Saturate: '<S4>/Back motor: Saturation' */
  if (rtb_Backgain > p2p4_P.BackmotorSaturation_UpperSat) {
    p2p4_B.BackmotorSaturation = p2p4_P.BackmotorSaturation_UpperSat;
  } else if (rtb_Backgain < p2p4_P.BackmotorSaturation_LowerSat) {
    p2p4_B.BackmotorSaturation = p2p4_P.BackmotorSaturation_LowerSat;
  } else {
    p2p4_B.BackmotorSaturation = rtb_Backgain;
  }

  /* End of Saturate: '<S4>/Back motor: Saturation' */
  if (rtmIsMajorTimeStep(p2p4_M)) {
  }

  /* Gain: '<S1>/Front gain' incorporates:
   *  Sum: '<S1>/Add'
   */
  rtb_Backgain = (p2p4_B.Sum1 + p2p4_B.Sum2) * p2p4_P.Frontgain_Gain;

  /* Saturate: '<S4>/Front motor: Saturation' */
  if (rtb_Backgain > p2p4_P.FrontmotorSaturation_UpperSat) {
    p2p4_B.FrontmotorSaturation = p2p4_P.FrontmotorSaturation_UpperSat;
  } else if (rtb_Backgain < p2p4_P.FrontmotorSaturation_LowerSat) {
    p2p4_B.FrontmotorSaturation = p2p4_P.FrontmotorSaturation_LowerSat;
  } else {
    p2p4_B.FrontmotorSaturation = rtb_Backgain;
  }

  /* End of Saturate: '<S4>/Front motor: Saturation' */
  if (rtmIsMajorTimeStep(p2p4_M)) {
    /* S-Function (hil_write_analog_block): '<S4>/HIL Write Analog' */

    /* S-Function Block: p2p4/Helicopter_interface/HIL Write Analog (hil_write_analog_block) */
    {
      t_error result;
      p2p4_DW.HILWriteAnalog_Buffer[0] = p2p4_B.FrontmotorSaturation;
      p2p4_DW.HILWriteAnalog_Buffer[1] = p2p4_B.BackmotorSaturation;
      result = hil_write_analog(p2p4_DW.HILInitialize_Card,
        p2p4_P.HILWriteAnalog_channels, 2, &p2p4_DW.HILWriteAnalog_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(p2p4_M, _rt_error_message);
      }
    }
  }
}

/* Model update function */
void p2p4_update(void)
{
  real_T *lastU;

  /* Update for Derivative: '<S4>/Derivative' */
  if (p2p4_DW.TimeStampA == (rtInf)) {
    p2p4_DW.TimeStampA = p2p4_M->Timing.t[0];
    lastU = &p2p4_DW.LastUAtTimeA;
  } else if (p2p4_DW.TimeStampB == (rtInf)) {
    p2p4_DW.TimeStampB = p2p4_M->Timing.t[0];
    lastU = &p2p4_DW.LastUAtTimeB;
  } else if (p2p4_DW.TimeStampA < p2p4_DW.TimeStampB) {
    p2p4_DW.TimeStampA = p2p4_M->Timing.t[0];
    lastU = &p2p4_DW.LastUAtTimeA;
  } else {
    p2p4_DW.TimeStampB = p2p4_M->Timing.t[0];
    lastU = &p2p4_DW.LastUAtTimeB;
  }

  *lastU = p2p4_B.PitchCounttorad;

  /* End of Update for Derivative: '<S4>/Derivative' */
  if (rtmIsMajorTimeStep(p2p4_M)) {
    rt_ertODEUpdateContinuousStates(&p2p4_M->solverInfo);
  }

  /* Update absolute time for base rate */
  /* The "clockTick0" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick0"
   * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
   * overflow during the application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick0 and the high bits
   * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
   */
  if (!(++p2p4_M->Timing.clockTick0)) {
    ++p2p4_M->Timing.clockTickH0;
  }

  p2p4_M->Timing.t[0] = rtsiGetSolverStopTime(&p2p4_M->solverInfo);

  {
    /* Update absolute timer for sample time: [0.002s, 0.0s] */
    /* The "clockTick1" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick1"
     * and "Timing.stepSize1". Size of "clockTick1" ensures timer will not
     * overflow during the application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick1 and the high bits
     * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
     */
    if (!(++p2p4_M->Timing.clockTick1)) {
      ++p2p4_M->Timing.clockTickH1;
    }

    p2p4_M->Timing.t[1] = p2p4_M->Timing.clockTick1 * p2p4_M->Timing.stepSize1 +
      p2p4_M->Timing.clockTickH1 * p2p4_M->Timing.stepSize1 * 4294967296.0;
  }
}

/* Derivatives for root system: '<Root>' */
void p2p4_derivatives(void)
{
  boolean_T lsat;
  boolean_T usat;
  XDot_p2p4_T *_rtXdot;
  _rtXdot = ((XDot_p2p4_T *) p2p4_M->ModelData.derivs);

  /* Derivatives for TransferFcn: '<S4>/Travel: Transfer Fcn' */
  _rtXdot->TravelTransferFcn_CSTATE = 0.0;
  _rtXdot->TravelTransferFcn_CSTATE += p2p4_P.TravelTransferFcn_A *
    p2p4_X.TravelTransferFcn_CSTATE;
  _rtXdot->TravelTransferFcn_CSTATE += p2p4_B.TravelCounttorad;

  /* Derivatives for TransferFcn: '<S4>/Pitch: Transfer Fcn' */
  _rtXdot->PitchTransferFcn_CSTATE = 0.0;
  _rtXdot->PitchTransferFcn_CSTATE += p2p4_P.PitchTransferFcn_A *
    p2p4_X.PitchTransferFcn_CSTATE;
  _rtXdot->PitchTransferFcn_CSTATE += p2p4_B.PitchCounttorad;

  /* Derivatives for TransferFcn: '<S4>/Elevation: Transfer Fcn' */
  _rtXdot->ElevationTransferFcn_CSTATE = 0.0;
  _rtXdot->ElevationTransferFcn_CSTATE += p2p4_P.ElevationTransferFcn_A *
    p2p4_X.ElevationTransferFcn_CSTATE;
  _rtXdot->ElevationTransferFcn_CSTATE += p2p4_B.ElevationCounttorad;

  /* Derivatives for Integrator: '<S3>/Integrator' */
  lsat = (p2p4_X.Integrator_CSTATE <= p2p4_P.Integrator_LowerSat);
  usat = (p2p4_X.Integrator_CSTATE >= p2p4_P.Integrator_UpperSat);
  if (((!lsat) && (!usat)) || (lsat && (p2p4_B.K_ei > 0.0)) || (usat &&
       (p2p4_B.K_ei < 0.0))) {
    _rtXdot->Integrator_CSTATE = p2p4_B.K_ei;
  } else {
    /* in saturation */
    _rtXdot->Integrator_CSTATE = 0.0;
  }

  /* End of Derivatives for Integrator: '<S3>/Integrator' */
}

/* Model initialize function */
void p2p4_initialize(void)
{
  /* Start for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: p2p4/HIL Initialize (hil_initialize_block) */
  {
    t_int result;
    t_boolean is_switching;
    result = hil_open("q8_usb", "0", &p2p4_DW.HILInitialize_Card);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(p2p4_M, _rt_error_message);
      return;
    }

    is_switching = false;
    result = hil_set_card_specific_options(p2p4_DW.HILInitialize_Card,
      "update_rate=normal;decimation=1", 32);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(p2p4_M, _rt_error_message);
      return;
    }

    result = hil_watchdog_clear(p2p4_DW.HILInitialize_Card);
    if (result < 0 && result != -QERR_HIL_WATCHDOG_CLEAR) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(p2p4_M, _rt_error_message);
      return;
    }

    if ((p2p4_P.HILInitialize_set_analog_input_ && !is_switching) ||
        (p2p4_P.HILInitialize_set_analog_inpu_m && is_switching)) {
      {
        int_T i1;
        real_T *dw_AIMinimums = &p2p4_DW.HILInitialize_AIMinimums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMinimums[i1] = p2p4_P.HILInitialize_analog_input_mini;
        }
      }

      {
        int_T i1;
        real_T *dw_AIMaximums = &p2p4_DW.HILInitialize_AIMaximums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMaximums[i1] = p2p4_P.HILInitialize_analog_input_maxi;
        }
      }

      result = hil_set_analog_input_ranges(p2p4_DW.HILInitialize_Card,
        p2p4_P.HILInitialize_analog_input_chan, 8U,
        &p2p4_DW.HILInitialize_AIMinimums[0], &p2p4_DW.HILInitialize_AIMaximums
        [0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(p2p4_M, _rt_error_message);
        return;
      }
    }

    if ((p2p4_P.HILInitialize_set_analog_output && !is_switching) ||
        (p2p4_P.HILInitialize_set_analog_outp_b && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOMinimums = &p2p4_DW.HILInitialize_AOMinimums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMinimums[i1] = p2p4_P.HILInitialize_analog_output_min;
        }
      }

      {
        int_T i1;
        real_T *dw_AOMaximums = &p2p4_DW.HILInitialize_AOMaximums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMaximums[i1] = p2p4_P.HILInitialize_analog_output_max;
        }
      }

      result = hil_set_analog_output_ranges(p2p4_DW.HILInitialize_Card,
        p2p4_P.HILInitialize_analog_output_cha, 8U,
        &p2p4_DW.HILInitialize_AOMinimums[0], &p2p4_DW.HILInitialize_AOMaximums
        [0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(p2p4_M, _rt_error_message);
        return;
      }
    }

    if ((p2p4_P.HILInitialize_set_analog_outp_e && !is_switching) ||
        (p2p4_P.HILInitialize_set_analog_outp_j && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &p2p4_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = p2p4_P.HILInitialize_initial_analog_ou;
        }
      }

      result = hil_write_analog(p2p4_DW.HILInitialize_Card,
        p2p4_P.HILInitialize_analog_output_cha, 8U,
        &p2p4_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(p2p4_M, _rt_error_message);
        return;
      }
    }

    if (p2p4_P.HILInitialize_set_analog_outp_p) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &p2p4_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = p2p4_P.HILInitialize_watchdog_analog_o;
        }
      }

      result = hil_watchdog_set_analog_expiration_state
        (p2p4_DW.HILInitialize_Card, p2p4_P.HILInitialize_analog_output_cha, 8U,
         &p2p4_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(p2p4_M, _rt_error_message);
        return;
      }
    }

    if ((p2p4_P.HILInitialize_set_encoder_param && !is_switching) ||
        (p2p4_P.HILInitialize_set_encoder_par_m && is_switching)) {
      {
        int_T i1;
        int32_T *dw_QuadratureModes = &p2p4_DW.HILInitialize_QuadratureModes[0];
        for (i1=0; i1 < 8; i1++) {
          dw_QuadratureModes[i1] = p2p4_P.HILInitialize_quadrature;
        }
      }

      result = hil_set_encoder_quadrature_mode(p2p4_DW.HILInitialize_Card,
        p2p4_P.HILInitialize_encoder_channels, 8U, (t_encoder_quadrature_mode *)
        &p2p4_DW.HILInitialize_QuadratureModes[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(p2p4_M, _rt_error_message);
        return;
      }
    }

    if ((p2p4_P.HILInitialize_set_encoder_count && !is_switching) ||
        (p2p4_P.HILInitialize_set_encoder_cou_k && is_switching)) {
      {
        int_T i1;
        int32_T *dw_InitialEICounts = &p2p4_DW.HILInitialize_InitialEICounts[0];
        for (i1=0; i1 < 8; i1++) {
          dw_InitialEICounts[i1] = p2p4_P.HILInitialize_initial_encoder_c;
        }
      }

      result = hil_set_encoder_counts(p2p4_DW.HILInitialize_Card,
        p2p4_P.HILInitialize_encoder_channels, 8U,
        &p2p4_DW.HILInitialize_InitialEICounts[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(p2p4_M, _rt_error_message);
        return;
      }
    }

    if ((p2p4_P.HILInitialize_set_pwm_params_at && !is_switching) ||
        (p2p4_P.HILInitialize_set_pwm_params__f && is_switching)) {
      uint32_T num_duty_cycle_modes = 0;
      uint32_T num_frequency_modes = 0;

      {
        int_T i1;
        int32_T *dw_POModeValues = &p2p4_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] = p2p4_P.HILInitialize_pwm_modes;
        }
      }

      result = hil_set_pwm_mode(p2p4_DW.HILInitialize_Card,
        p2p4_P.HILInitialize_pwm_channels, 8U, (t_pwm_mode *)
        &p2p4_DW.HILInitialize_POModeValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(p2p4_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        const uint32_T *p_HILInitialize_pwm_channels =
          p2p4_P.HILInitialize_pwm_channels;
        int32_T *dw_POModeValues = &p2p4_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          if (dw_POModeValues[i1] == PWM_DUTY_CYCLE_MODE || dw_POModeValues[i1] ==
              PWM_ONE_SHOT_MODE || dw_POModeValues[i1] == PWM_TIME_MODE) {
            p2p4_DW.HILInitialize_POSortedChans[num_duty_cycle_modes] =
              p_HILInitialize_pwm_channels[i1];
            p2p4_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes] =
              p2p4_P.HILInitialize_pwm_frequency;
            num_duty_cycle_modes++;
          } else {
            p2p4_DW.HILInitialize_POSortedChans[7U - num_frequency_modes] =
              p_HILInitialize_pwm_channels[i1];
            p2p4_DW.HILInitialize_POSortedFreqs[7U - num_frequency_modes] =
              p2p4_P.HILInitialize_pwm_frequency;
            num_frequency_modes++;
          }
        }
      }

      if (num_duty_cycle_modes > 0) {
        result = hil_set_pwm_frequency(p2p4_DW.HILInitialize_Card,
          &p2p4_DW.HILInitialize_POSortedChans[0], num_duty_cycle_modes,
          &p2p4_DW.HILInitialize_POSortedFreqs[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(p2p4_M, _rt_error_message);
          return;
        }
      }

      if (num_frequency_modes > 0) {
        result = hil_set_pwm_duty_cycle(p2p4_DW.HILInitialize_Card,
          &p2p4_DW.HILInitialize_POSortedChans[num_duty_cycle_modes],
          num_frequency_modes,
          &p2p4_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(p2p4_M, _rt_error_message);
          return;
        }
      }

      {
        int_T i1;
        int32_T *dw_POModeValues = &p2p4_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] = p2p4_P.HILInitialize_pwm_configuration;
        }
      }

      {
        int_T i1;
        int32_T *dw_POAlignValues = &p2p4_DW.HILInitialize_POAlignValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POAlignValues[i1] = p2p4_P.HILInitialize_pwm_alignment;
        }
      }

      {
        int_T i1;
        int32_T *dw_POPolarityVals = &p2p4_DW.HILInitialize_POPolarityVals[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POPolarityVals[i1] = p2p4_P.HILInitialize_pwm_polarity;
        }
      }

      result = hil_set_pwm_configuration(p2p4_DW.HILInitialize_Card,
        p2p4_P.HILInitialize_pwm_channels, 8U,
        (t_pwm_configuration *) &p2p4_DW.HILInitialize_POModeValues[0],
        (t_pwm_alignment *) &p2p4_DW.HILInitialize_POAlignValues[0],
        (t_pwm_polarity *) &p2p4_DW.HILInitialize_POPolarityVals[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(p2p4_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        real_T *dw_POSortedFreqs = &p2p4_DW.HILInitialize_POSortedFreqs[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POSortedFreqs[i1] = p2p4_P.HILInitialize_pwm_leading_deadb;
        }
      }

      {
        int_T i1;
        real_T *dw_POValues = &p2p4_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = p2p4_P.HILInitialize_pwm_trailing_dead;
        }
      }

      result = hil_set_pwm_deadband(p2p4_DW.HILInitialize_Card,
        p2p4_P.HILInitialize_pwm_channels, 8U,
        &p2p4_DW.HILInitialize_POSortedFreqs[0],
        &p2p4_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(p2p4_M, _rt_error_message);
        return;
      }
    }

    if ((p2p4_P.HILInitialize_set_pwm_outputs_a && !is_switching) ||
        (p2p4_P.HILInitialize_set_pwm_outputs_g && is_switching)) {
      {
        int_T i1;
        real_T *dw_POValues = &p2p4_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = p2p4_P.HILInitialize_initial_pwm_outpu;
        }
      }

      result = hil_write_pwm(p2p4_DW.HILInitialize_Card,
        p2p4_P.HILInitialize_pwm_channels, 8U, &p2p4_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(p2p4_M, _rt_error_message);
        return;
      }
    }

    if (p2p4_P.HILInitialize_set_pwm_outputs_o) {
      {
        int_T i1;
        real_T *dw_POValues = &p2p4_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = p2p4_P.HILInitialize_watchdog_pwm_outp;
        }
      }

      result = hil_watchdog_set_pwm_expiration_state(p2p4_DW.HILInitialize_Card,
        p2p4_P.HILInitialize_pwm_channels, 8U, &p2p4_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(p2p4_M, _rt_error_message);
        return;
      }
    }
  }

  /* Start for S-Function (hil_read_encoder_timebase_block): '<S4>/HIL Read Encoder Timebase' */

  /* S-Function Block: p2p4/Helicopter_interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
  {
    t_error result;
    result = hil_task_create_encoder_reader(p2p4_DW.HILInitialize_Card,
      p2p4_P.HILReadEncoderTimebase_samples_,
      p2p4_P.HILReadEncoderTimebase_channels, 3,
      &p2p4_DW.HILReadEncoderTimebase_Task);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(p2p4_M, _rt_error_message);
    }
  }

  /* Start for FromWorkspace: '<Root>/From Workspace' */
  {
    static real_T pTimeValues0[] = { 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75,
      2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0,
      5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.25,
      8.5, 8.75, 9.0, 9.25, 9.5, 9.75, 10.0, 10.25, 10.5, 10.75, 11.0, 11.25,
      11.5, 11.75, 12.0, 12.25, 12.5, 12.75, 13.0, 13.25, 13.5, 13.75, 14.0,
      14.25, 14.5, 14.75, 15.0, 15.25, 15.5, 15.75, 16.0, 16.25, 16.5, 16.75,
      17.0, 17.25, 17.5, 17.75, 18.0, 18.25, 18.5, 18.75, 19.0, 19.25, 19.5,
      19.75, 20.0, 20.25, 20.5, 20.75, 21.0, 21.25, 21.5, 21.75, 22.0, 22.25,
      22.5, 22.75, 23.0, 23.25, 23.5, 23.75, 24.0, 24.25, 24.5, 24.75, 25.0,
      25.25, 25.5, 25.75, 26.0, 26.25, 26.5, 26.75, 27.0, 27.25, 27.5, 27.75,
      28.0, 28.25, 28.5, 28.75, 29.0, 29.25, 29.5, 29.75, 30.0, 30.25, 30.5,
      30.75, 31.0, 31.25, 31.5, 31.75, 32.0, 32.25, 32.5, 32.75, 33.0, 33.25,
      33.5, 33.75, 34.0, 34.25, 34.5, 34.75, 35.0 } ;

    static real_T pDataValues0[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.7606093265025207,
      2.0496073918518563, 1.4382927083499606, 0.92157629518122308,
      0.49323105099131137, 0.14610511569282503, -0.12760502713794386,
      -0.33606725191635389, -0.48754270917254883, -0.590148703482257,
      -0.65166399356075511, -0.6793787131646698, -0.67998716766644085,
      -0.65951947344542594, -0.62330684835752459, -0.575974945702791,
      -0.52145968099742368, -0.46304034562984331, -0.40338530948396945,
      -0.34460620252010288, -0.28831707788517291, -0.23569566025313327,
      -0.18754434972771383, -0.14434916978229739, -0.10633530970930938,
      -0.073518314714462463, -0.045750320079502194, -0.022761011888698279,
      -0.00419322931932764, 0.010366692955967038, 0.02136360395583381,
      0.029252897511381715, 0.034484147304728929, 0.037488567209517654,
      0.038669830558082152, 0.0383977825026004, 0.037004594391614076,
      0.03478293443028535, 0.0319857616935626, 0.028827388146451027,
      0.025485493472495846, 0.022103818414220378, 0.018795302560499327,
      0.015645470987160182, 0.012715910076626424, 0.010047705673226215,
      0.0076647461506429426, 0.0055768188314175848, 0.0037824505018039245,
      0.0022714616154070333, 0.0010272193651329475, 2.8587372612021156E-5,
      -0.00074842041735889708, -0.001329266636928228, -0.0017397778897885186,
      -0.00200526884603347, -0.0021498772596666438, -0.0021960833829286486,
      -0.0021643874864025747, -0.0020731201853727566, -0.0019383618220252252,
      -0.0017739490878645459, -0.0015915492456661795, -0.0014007846049355344,
      -0.0012093922218763264, -0.0010234060579714654, -0.00084735098279264953,
      -0.000684440005219944, -0.00053676793542756229, -0.00040549630201924813,
      -0.0002910257683661488, -0.00019315351095672284, -0.00011121404769309784,
      -4.4202847244267929E-5, 9.1172735158240314E-6, 5.0126431703656493E-5,
      8.02737157993283E-5, 0.00010101777681525981, 0.00011378059275624075,
      0.00011991292783586349, 0.00012066997537696356, 0.00011719569189794221,
      0.00011051438650855339, 0.00010152821253883726, 9.1019306492579055E-5,
      7.96554228262858E-5, 6.7998011830106009E-5, 5.6511771790224152E-5,
      4.5574764561684411E-5, 3.5488203458777789E-5, 2.6484992005859631E-5,
      1.8736004565589326E-5, 1.2352966934914278E-5, 7.3866806228339854E-6,
      3.8194316429596142E-6, 1.5512232403036563E-6, 3.8211673136778757E-7, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } ;

    p2p4_DW.FromWorkspace_PWORK.TimePtr = (void *) pTimeValues0;
    p2p4_DW.FromWorkspace_PWORK.DataPtr = (void *) pDataValues0;
    p2p4_DW.FromWorkspace_IWORK.PrevIndex = 0;
  }

  /* Start for ToFile: '<Root>/To File' */
  {
    FILE *fp = (NULL);
    char fileName[509] = "p2p4_1.mat";
    if ((fp = fopen(fileName, "wb")) == (NULL)) {
      rtmSetErrorStatus(p2p4_M, "Error creating .mat file p2p4_1.mat");
      return;
    }

    if (rt_WriteMat4FileHeader(fp,6,0,"ulambdarppdot")) {
      rtmSetErrorStatus(p2p4_M,
                        "Error writing mat file header to file p2p4_1.mat");
      return;
    }

    p2p4_DW.ToFile_IWORK.Count = 0;
    p2p4_DW.ToFile_IWORK.Decimation = -1;
    p2p4_DW.ToFile_PWORK.FilePtr = fp;
  }

  /* InitializeConditions for TransferFcn: '<S4>/Travel: Transfer Fcn' */
  p2p4_X.TravelTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S4>/Pitch: Transfer Fcn' */
  p2p4_X.PitchTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S4>/Elevation: Transfer Fcn' */
  p2p4_X.ElevationTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for Integrator: '<S3>/Integrator' */
  p2p4_X.Integrator_CSTATE = p2p4_P.Integrator_IC;

  /* InitializeConditions for Derivative: '<S4>/Derivative' */
  p2p4_DW.TimeStampA = (rtInf);
  p2p4_DW.TimeStampB = (rtInf);
}

/* Model terminate function */
void p2p4_terminate(void)
{
  /* Terminate for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: p2p4/HIL Initialize (hil_initialize_block) */
  {
    t_boolean is_switching;
    t_int result;
    t_uint32 num_final_analog_outputs = 0;
    t_uint32 num_final_pwm_outputs = 0;
    hil_task_stop_all(p2p4_DW.HILInitialize_Card);
    hil_monitor_stop_all(p2p4_DW.HILInitialize_Card);
    is_switching = false;
    if ((p2p4_P.HILInitialize_set_analog_out_ex && !is_switching) ||
        (p2p4_P.HILInitialize_set_analog_outp_c && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &p2p4_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = p2p4_P.HILInitialize_final_analog_outp;
        }
      }

      num_final_analog_outputs = 8U;
    }

    if ((p2p4_P.HILInitialize_set_pwm_output_ap && !is_switching) ||
        (p2p4_P.HILInitialize_set_pwm_outputs_p && is_switching)) {
      {
        int_T i1;
        real_T *dw_POValues = &p2p4_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = p2p4_P.HILInitialize_final_pwm_outputs;
        }
      }

      num_final_pwm_outputs = 8U;
    }

    if (0
        || num_final_analog_outputs > 0
        || num_final_pwm_outputs > 0
        ) {
      /* Attempt to write the final outputs atomically (due to firmware issue in old Q2-USB). Otherwise write channels individually */
      result = hil_write(p2p4_DW.HILInitialize_Card
                         , p2p4_P.HILInitialize_analog_output_cha,
                         num_final_analog_outputs
                         , p2p4_P.HILInitialize_pwm_channels,
                         num_final_pwm_outputs
                         , NULL, 0
                         , NULL, 0
                         , &p2p4_DW.HILInitialize_AOVoltages[0]
                         , &p2p4_DW.HILInitialize_POValues[0]
                         , (t_boolean *) NULL
                         , NULL
                         );
      if (result == -QERR_HIL_WRITE_NOT_SUPPORTED) {
        t_error local_result;
        result = 0;

        /* The hil_write operation is not supported by this card. Write final outputs for each channel type */
        if (num_final_analog_outputs > 0) {
          local_result = hil_write_analog(p2p4_DW.HILInitialize_Card,
            p2p4_P.HILInitialize_analog_output_cha, num_final_analog_outputs,
            &p2p4_DW.HILInitialize_AOVoltages[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (num_final_pwm_outputs > 0) {
          local_result = hil_write_pwm(p2p4_DW.HILInitialize_Card,
            p2p4_P.HILInitialize_pwm_channels, num_final_pwm_outputs,
            &p2p4_DW.HILInitialize_POValues[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(p2p4_M, _rt_error_message);
        }
      }
    }

    hil_task_delete_all(p2p4_DW.HILInitialize_Card);
    hil_monitor_delete_all(p2p4_DW.HILInitialize_Card);
    hil_close(p2p4_DW.HILInitialize_Card);
    p2p4_DW.HILInitialize_Card = NULL;
  }

  /* Terminate for ToFile: '<Root>/To File' */
  {
    FILE *fp = (FILE *) p2p4_DW.ToFile_PWORK.FilePtr;
    if (fp != (NULL)) {
      char fileName[509] = "p2p4_1.mat";
      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(p2p4_M, "Error closing MAT-file p2p4_1.mat");
        return;
      }

      if ((fp = fopen(fileName, "r+b")) == (NULL)) {
        rtmSetErrorStatus(p2p4_M, "Error reopening MAT-file p2p4_1.mat");
        return;
      }

      if (rt_WriteMat4FileHeader(fp, 6, p2p4_DW.ToFile_IWORK.Count,
           "ulambdarppdot")) {
        rtmSetErrorStatus(p2p4_M,
                          "Error writing header for ulambdarppdot to MAT-file p2p4_1.mat");
      }

      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(p2p4_M, "Error closing MAT-file p2p4_1.mat");
        return;
      }

      p2p4_DW.ToFile_PWORK.FilePtr = (NULL);
    }
  }
}

/*========================================================================*
 * Start of Classic call interface                                        *
 *========================================================================*/

/* Solver interface called by GRT_Main */
#ifndef USE_GENERATED_SOLVER

void rt_ODECreateIntegrationData(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

void rt_ODEDestroyIntegrationData(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

void rt_ODEUpdateContinuousStates(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

#endif

void MdlOutputs(int_T tid)
{
  p2p4_output();
  UNUSED_PARAMETER(tid);
}

void MdlUpdate(int_T tid)
{
  p2p4_update();
  UNUSED_PARAMETER(tid);
}

void MdlInitializeSizes(void)
{
}

void MdlInitializeSampleTimes(void)
{
}

void MdlInitialize(void)
{
}

void MdlStart(void)
{
  p2p4_initialize();
}

void MdlTerminate(void)
{
  p2p4_terminate();
}

/* Registration function */
RT_MODEL_p2p4_T *p2p4(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* non-finite (run-time) assignments */
  p2p4_P.Integrator_UpperSat = rtInf;
  p2p4_P.Integrator_LowerSat = rtMinusInf;

  /* initialize real-time model */
  (void) memset((void *)p2p4_M, 0,
                sizeof(RT_MODEL_p2p4_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&p2p4_M->solverInfo, &p2p4_M->Timing.simTimeStep);
    rtsiSetTPtr(&p2p4_M->solverInfo, &rtmGetTPtr(p2p4_M));
    rtsiSetStepSizePtr(&p2p4_M->solverInfo, &p2p4_M->Timing.stepSize0);
    rtsiSetdXPtr(&p2p4_M->solverInfo, &p2p4_M->ModelData.derivs);
    rtsiSetContStatesPtr(&p2p4_M->solverInfo, (real_T **)
                         &p2p4_M->ModelData.contStates);
    rtsiSetNumContStatesPtr(&p2p4_M->solverInfo, &p2p4_M->Sizes.numContStates);
    rtsiSetNumPeriodicContStatesPtr(&p2p4_M->solverInfo,
      &p2p4_M->Sizes.numPeriodicContStates);
    rtsiSetPeriodicContStateIndicesPtr(&p2p4_M->solverInfo,
      &p2p4_M->ModelData.periodicContStateIndices);
    rtsiSetPeriodicContStateRangesPtr(&p2p4_M->solverInfo,
      &p2p4_M->ModelData.periodicContStateRanges);
    rtsiSetErrorStatusPtr(&p2p4_M->solverInfo, (&rtmGetErrorStatus(p2p4_M)));
    rtsiSetRTModelPtr(&p2p4_M->solverInfo, p2p4_M);
  }

  rtsiSetSimTimeStep(&p2p4_M->solverInfo, MAJOR_TIME_STEP);
  p2p4_M->ModelData.intgData.f[0] = p2p4_M->ModelData.odeF[0];
  p2p4_M->ModelData.contStates = ((real_T *) &p2p4_X);
  rtsiSetSolverData(&p2p4_M->solverInfo, (void *)&p2p4_M->ModelData.intgData);
  rtsiSetSolverName(&p2p4_M->solverInfo,"ode1");

  /* Initialize timing info */
  {
    int_T *mdlTsMap = p2p4_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    mdlTsMap[1] = 1;
    p2p4_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    p2p4_M->Timing.sampleTimes = (&p2p4_M->Timing.sampleTimesArray[0]);
    p2p4_M->Timing.offsetTimes = (&p2p4_M->Timing.offsetTimesArray[0]);

    /* task periods */
    p2p4_M->Timing.sampleTimes[0] = (0.0);
    p2p4_M->Timing.sampleTimes[1] = (0.002);

    /* task offsets */
    p2p4_M->Timing.offsetTimes[0] = (0.0);
    p2p4_M->Timing.offsetTimes[1] = (0.0);
  }

  rtmSetTPtr(p2p4_M, &p2p4_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = p2p4_M->Timing.sampleHitArray;
    mdlSampleHits[0] = 1;
    mdlSampleHits[1] = 1;
    p2p4_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(p2p4_M, -1);
  p2p4_M->Timing.stepSize0 = 0.002;
  p2p4_M->Timing.stepSize1 = 0.002;

  /* External mode info */
  p2p4_M->Sizes.checksums[0] = (1655397102U);
  p2p4_M->Sizes.checksums[1] = (2726450126U);
  p2p4_M->Sizes.checksums[2] = (888731761U);
  p2p4_M->Sizes.checksums[3] = (172524935U);

  {
    static const sysRanDType rtAlwaysEnabled = SUBSYS_RAN_BC_ENABLE;
    static RTWExtModeInfo rt_ExtModeInfo;
    static const sysRanDType *systemRan[1];
    p2p4_M->extModeInfo = (&rt_ExtModeInfo);
    rteiSetSubSystemActiveVectorAddresses(&rt_ExtModeInfo, systemRan);
    systemRan[0] = &rtAlwaysEnabled;
    rteiSetModelMappingInfoPtr(p2p4_M->extModeInfo,
      &p2p4_M->SpecialInfo.mappingInfo);
    rteiSetChecksumsPtr(p2p4_M->extModeInfo, p2p4_M->Sizes.checksums);
    rteiSetTPtr(p2p4_M->extModeInfo, rtmGetTPtr(p2p4_M));
  }

  p2p4_M->solverInfoPtr = (&p2p4_M->solverInfo);
  p2p4_M->Timing.stepSize = (0.002);
  rtsiSetFixedStepSize(&p2p4_M->solverInfo, 0.002);
  rtsiSetSolverMode(&p2p4_M->solverInfo, SOLVER_MODE_SINGLETASKING);

  /* block I/O */
  p2p4_M->ModelData.blockIO = ((void *) &p2p4_B);

  {
    int32_T i;
    for (i = 0; i < 6; i++) {
      p2p4_B.Gain1[i] = 0.0;
    }

    p2p4_B.TravelCounttorad = 0.0;
    p2p4_B.Gain = 0.0;
    p2p4_B.Gain_d = 0.0;
    p2p4_B.PitchCounttorad = 0.0;
    p2p4_B.Gain_i = 0.0;
    p2p4_B.Gain_b = 0.0;
    p2p4_B.ElevationCounttorad = 0.0;
    p2p4_B.Gain_e = 0.0;
    p2p4_B.Sum = 0.0;
    p2p4_B.Gain_dg = 0.0;
    p2p4_B.Sum1 = 0.0;
    p2p4_B.Sum2 = 0.0;
    p2p4_B.K_ei = 0.0;
    p2p4_B.Gain_l = 0.0;
    p2p4_B.BackmotorSaturation = 0.0;
    p2p4_B.FrontmotorSaturation = 0.0;
  }

  /* parameters */
  p2p4_M->ModelData.defaultParam = ((real_T *)&p2p4_P);

  /* states (continuous) */
  {
    real_T *x = (real_T *) &p2p4_X;
    p2p4_M->ModelData.contStates = (x);
    (void) memset((void *)&p2p4_X, 0,
                  sizeof(X_p2p4_T));
  }

  /* states (dwork) */
  p2p4_M->ModelData.dwork = ((void *) &p2p4_DW);
  (void) memset((void *)&p2p4_DW, 0,
                sizeof(DW_p2p4_T));

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      p2p4_DW.HILInitialize_AIMinimums[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      p2p4_DW.HILInitialize_AIMaximums[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      p2p4_DW.HILInitialize_AOMinimums[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      p2p4_DW.HILInitialize_AOMaximums[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      p2p4_DW.HILInitialize_AOVoltages[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      p2p4_DW.HILInitialize_FilterFrequency[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      p2p4_DW.HILInitialize_POSortedFreqs[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      p2p4_DW.HILInitialize_POValues[i] = 0.0;
    }
  }

  p2p4_DW.TimeStampA = 0.0;
  p2p4_DW.LastUAtTimeA = 0.0;
  p2p4_DW.TimeStampB = 0.0;
  p2p4_DW.LastUAtTimeB = 0.0;
  p2p4_DW.HILWriteAnalog_Buffer[0] = 0.0;
  p2p4_DW.HILWriteAnalog_Buffer[1] = 0.0;

  /* data type transition information */
  {
    static DataTypeTransInfo dtInfo;
    (void) memset((char_T *) &dtInfo, 0,
                  sizeof(dtInfo));
    p2p4_M->SpecialInfo.mappingInfo = (&dtInfo);
    dtInfo.numDataTypes = 16;
    dtInfo.dataTypeSizes = &rtDataTypeSizes[0];
    dtInfo.dataTypeNames = &rtDataTypeNames[0];

    /* Block I/O transition table */
    dtInfo.B = &rtBTransTable;

    /* Parameters transition table */
    dtInfo.P = &rtPTransTable;
  }

  /* Initialize Sizes */
  p2p4_M->Sizes.numContStates = (4);   /* Number of continuous states */
  p2p4_M->Sizes.numPeriodicContStates = (0);/* Number of periodic continuous states */
  p2p4_M->Sizes.numY = (0);            /* Number of model outputs */
  p2p4_M->Sizes.numU = (0);            /* Number of model inputs */
  p2p4_M->Sizes.sysDirFeedThru = (0);  /* The model is not direct feedthrough */
  p2p4_M->Sizes.numSampTimes = (2);    /* Number of sample times */
  p2p4_M->Sizes.numBlocks = (56);      /* Number of blocks */
  p2p4_M->Sizes.numBlockIO = (17);     /* Number of block outputs */
  p2p4_M->Sizes.numBlockPrms = (142);  /* Sum of parameter "widths" */
  return p2p4_M;
}

/*========================================================================*
 * End of Classic call interface                                          *
 *========================================================================*/
