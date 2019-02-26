/*
 * p4p4.c
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "p4p4".
 *
 * Model version              : 1.182
 * Simulink Coder version : 8.9 (R2015b) 13-Aug-2015
 * C source code generated on : Wed Feb 20 12:13:33 2019
 *
 * Target selection: quarc_win64.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "p4p4.h"
#include "p4p4_private.h"
#include "p4p4_dt.h"

/* Block signals (auto storage) */
B_p4p4_T p4p4_B;

/* Continuous states */
X_p4p4_T p4p4_X;

/* Block states (auto storage) */
DW_p4p4_T p4p4_DW;

/* Real-time model */
RT_MODEL_p4p4_T p4p4_M_;
RT_MODEL_p4p4_T *const p4p4_M = &p4p4_M_;

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
  p4p4_derivatives();
  rtsiSetT(si, tnew);
  for (i = 0; i < nXc; ++i) {
    x[i] += h * f0[i];
  }

  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/* Model output function */
void p4p4_output(void)
{
  /* local block i/o variables */
  real_T rtb_HILReadEncoderTimebase_o1;
  real_T rtb_HILReadEncoderTimebase_o2;
  real_T rtb_HILReadEncoderTimebase_o3;
  real_T rtb_TmpSignalConversionAtToFile[4];
  real_T *lastU;
  real_T rtb_Backgain;
  int32_T i;
  real_T tmp[6];
  real_T tmp_0[6];
  int32_T i_0;
  if (rtmIsMajorTimeStep(p4p4_M)) {
    /* set solver stop time */
    if (!(p4p4_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&p4p4_M->solverInfo, ((p4p4_M->Timing.clockTickH0 +
        1) * p4p4_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&p4p4_M->solverInfo, ((p4p4_M->Timing.clockTick0 + 1)
        * p4p4_M->Timing.stepSize0 + p4p4_M->Timing.clockTickH0 *
        p4p4_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(p4p4_M)) {
    p4p4_M->Timing.t[0] = rtsiGetT(&p4p4_M->solverInfo);
  }

  if (rtmIsMajorTimeStep(p4p4_M)) {
    /* S-Function (hil_read_encoder_timebase_block): '<S4>/HIL Read Encoder Timebase' */

    /* S-Function Block: p4p4/Helicopter_interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
    {
      t_error result;
      result = hil_task_read_encoder(p4p4_DW.HILReadEncoderTimebase_Task, 1,
        &p4p4_DW.HILReadEncoderTimebase_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(p4p4_M, _rt_error_message);
      } else {
        rtb_HILReadEncoderTimebase_o1 = p4p4_DW.HILReadEncoderTimebase_Buffer[0];
        rtb_HILReadEncoderTimebase_o2 = p4p4_DW.HILReadEncoderTimebase_Buffer[1];
        rtb_HILReadEncoderTimebase_o3 = p4p4_DW.HILReadEncoderTimebase_Buffer[2];
      }
    }
  }

  /* FromWorkspace: '<Root>/u Workspace' */
  {
    real_T *pDataValues = (real_T *) p4p4_DW.uWorkspace_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *) p4p4_DW.uWorkspace_PWORK.TimePtr;
    int_T currTimeIndex = p4p4_DW.uWorkspace_IWORK.PrevIndex;
    real_T t = p4p4_M->Timing.t[0];

    /* Get index */
    if (t <= pTimeValues[0]) {
      currTimeIndex = 0;
    } else if (t >= pTimeValues[104]) {
      currTimeIndex = 103;
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

    p4p4_DW.uWorkspace_IWORK.PrevIndex = currTimeIndex;

    /* Post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 2; ++elIdx) {
              (&p4p4_B.uWorkspace[0])[elIdx] = pDataValues[currTimeIndex];
              pDataValues += 105;
            }
          }
        } else {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 2; ++elIdx) {
              (&p4p4_B.uWorkspace[0])[elIdx] = pDataValues[currTimeIndex + 1];
              pDataValues += 105;
            }
          }
        }
      } else {
        real_T f1 = (t2 - t) / (t2 - t1);
        real_T f2 = 1.0 - f1;
        real_T d1;
        real_T d2;
        int_T TimeIndex= currTimeIndex;

        {
          int_T elIdx;
          for (elIdx = 0; elIdx < 2; ++elIdx) {
            d1 = pDataValues[TimeIndex];
            d2 = pDataValues[TimeIndex + 1];
            (&p4p4_B.uWorkspace[0])[elIdx] = (real_T) rtInterpolate(d1, d2, f1,
              f2);
            pDataValues += 105;
          }
        }
      }
    }
  }

  /* FromWorkspace: '<Root>/x Workspace1' */
  {
    real_T *pDataValues = (real_T *) p4p4_DW.xWorkspace1_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *) p4p4_DW.xWorkspace1_PWORK.TimePtr;
    int_T currTimeIndex = p4p4_DW.xWorkspace1_IWORK.PrevIndex;
    real_T t = p4p4_M->Timing.t[0];

    /* Get index */
    if (t <= pTimeValues[0]) {
      currTimeIndex = 0;
    } else if (t >= pTimeValues[104]) {
      currTimeIndex = 103;
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

    p4p4_DW.xWorkspace1_IWORK.PrevIndex = currTimeIndex;

    /* Post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 6; ++elIdx) {
              (&p4p4_B.xWorkspace1[0])[elIdx] = pDataValues[currTimeIndex];
              pDataValues += 105;
            }
          }
        } else {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 6; ++elIdx) {
              (&p4p4_B.xWorkspace1[0])[elIdx] = pDataValues[currTimeIndex + 1];
              pDataValues += 105;
            }
          }
        }
      } else {
        real_T f1 = (t2 - t) / (t2 - t1);
        real_T f2 = 1.0 - f1;
        real_T d1;
        real_T d2;
        int_T TimeIndex= currTimeIndex;

        {
          int_T elIdx;
          for (elIdx = 0; elIdx < 6; ++elIdx) {
            d1 = pDataValues[TimeIndex];
            d2 = pDataValues[TimeIndex + 1];
            (&p4p4_B.xWorkspace1[0])[elIdx] = (real_T) rtInterpolate(d1, d2, f1,
              f2);
            pDataValues += 105;
          }
        }
      }
    }
  }

  if (rtmIsMajorTimeStep(p4p4_M)) {
    /* Gain: '<S4>/Travel: Count to rad' incorporates:
     *  Gain: '<S4>/Travel_gain'
     */
    p4p4_B.TravelCounttorad = p4p4_P.travel_gain * rtb_HILReadEncoderTimebase_o1
      * p4p4_P.TravelCounttorad_Gain;

    /* Gain: '<S12>/Gain' */
    p4p4_B.Gain = p4p4_P.Gain_Gain * p4p4_B.TravelCounttorad;

    /* Gain: '<S4>/Pitch: Count to rad' */
    p4p4_B.PitchCounttorad = p4p4_P.PitchCounttorad_Gain *
      rtb_HILReadEncoderTimebase_o2;

    /* Gain: '<S9>/Gain' */
    p4p4_B.Gain_i = p4p4_P.Gain_Gain_a * p4p4_B.PitchCounttorad;
  }

  /* Gain: '<S13>/Gain' incorporates:
   *  TransferFcn: '<S4>/Travel: Transfer Fcn'
   */
  p4p4_B.Gain_d = (p4p4_P.TravelTransferFcn_C * p4p4_X.TravelTransferFcn_CSTATE
                   + p4p4_P.TravelTransferFcn_D * p4p4_B.TravelCounttorad) *
    p4p4_P.Gain_Gain_l;

  /* Gain: '<S10>/Gain' incorporates:
   *  TransferFcn: '<S4>/Pitch: Transfer Fcn'
   */
  p4p4_B.Gain_b = (p4p4_P.PitchTransferFcn_C * p4p4_X.PitchTransferFcn_CSTATE +
                   p4p4_P.PitchTransferFcn_D * p4p4_B.PitchCounttorad) *
    p4p4_P.Gain_Gain_ae;
  if (rtmIsMajorTimeStep(p4p4_M)) {
    /* Gain: '<S4>/Elevation: Count to rad' */
    p4p4_B.ElevationCounttorad = p4p4_P.ElevationCounttorad_Gain *
      rtb_HILReadEncoderTimebase_o3;

    /* Gain: '<S7>/Gain' */
    p4p4_B.Gain_e = p4p4_P.Gain_Gain_lv * p4p4_B.ElevationCounttorad;

    /* Sum: '<Root>/Sum' incorporates:
     *  Constant: '<Root>/elavation_offset [deg]'
     */
    p4p4_B.Sum = p4p4_B.Gain_e + p4p4_P.elavation_offsetdeg_Value;
  }

  /* Gain: '<S8>/Gain' incorporates:
   *  TransferFcn: '<S4>/Elevation: Transfer Fcn'
   */
  p4p4_B.Gain_dg = (p4p4_P.ElevationTransferFcn_C *
                    p4p4_X.ElevationTransferFcn_CSTATE +
                    p4p4_P.ElevationTransferFcn_D * p4p4_B.ElevationCounttorad) *
    p4p4_P.Gain_Gain_n;

  /* Gain: '<S2>/Gain1' */
  p4p4_B.Gain1[0] = p4p4_P.Gain1_Gain * p4p4_B.Gain;
  p4p4_B.Gain1[1] = p4p4_P.Gain1_Gain * p4p4_B.Gain_d;
  p4p4_B.Gain1[2] = p4p4_P.Gain1_Gain * p4p4_B.Gain_i;
  p4p4_B.Gain1[3] = p4p4_P.Gain1_Gain * p4p4_B.Gain_b;
  p4p4_B.Gain1[4] = p4p4_P.Gain1_Gain * p4p4_B.Sum;
  p4p4_B.Gain1[5] = p4p4_P.Gain1_Gain * p4p4_B.Gain_dg;

  /* Sum: '<Root>/Sum5' incorporates:
   *  Constant: '<Root>/pi '
   */
  p4p4_B.Sum5 = p4p4_P.pi_Value + p4p4_B.Gain1[0];

  /* Sum: '<S5>/Sum4' incorporates:
   *  Gain: '<S5>/Gain'
   */
  tmp[0] = p4p4_B.Sum5;
  for (i = 0; i < 5; i++) {
    tmp[i + 1] = p4p4_B.Gain1[1 + i];
  }

  for (i = 0; i < 6; i++) {
    tmp_0[i] = tmp[i] - p4p4_B.xWorkspace1[i];
  }

  /* End of Sum: '<S5>/Sum4' */

  /* Sum: '<S5>/Sum3' incorporates:
   *  Gain: '<S5>/Gain'
   */
  for (i = 0; i < 2; i++) {
    rtb_Backgain = 0.0;
    for (i_0 = 0; i_0 < 6; i_0++) {
      rtb_Backgain += p4p4_P.K_lqr[(i_0 << 1) + i] * tmp_0[i_0];
    }

    p4p4_B.Sum3[i] = p4p4_B.uWorkspace[i] - rtb_Backgain;
  }

  /* End of Sum: '<S5>/Sum3' */
  if (rtmIsMajorTimeStep(p4p4_M)) {
  }

  /* Sum: '<Root>/Sum1' incorporates:
   *  Constant: '<Root>/Vd_bias'
   *  Gain: '<S6>/K_pd'
   *  Gain: '<S6>/K_pp'
   *  Sum: '<S6>/Sum2'
   *  Sum: '<S6>/Sum3'
   */
  p4p4_B.Sum1 = ((p4p4_B.Sum3[0] - p4p4_B.Gain1[2]) * p4p4_P.K_pp - p4p4_P.K_pd *
                 p4p4_B.Gain1[3]) + p4p4_P.Vd_ff;
  if (rtmIsMajorTimeStep(p4p4_M)) {
  }

  /* Integrator: '<S3>/Integrator' */
  /* Limited  Integrator  */
  if (p4p4_X.Integrator_CSTATE >= p4p4_P.Integrator_UpperSat) {
    p4p4_X.Integrator_CSTATE = p4p4_P.Integrator_UpperSat;
  } else {
    if (p4p4_X.Integrator_CSTATE <= p4p4_P.Integrator_LowerSat) {
      p4p4_X.Integrator_CSTATE = p4p4_P.Integrator_LowerSat;
    }
  }

  /* Sum: '<S3>/Sum' */
  rtb_Backgain = p4p4_B.Sum3[1] - p4p4_B.Gain1[4];

  /* Sum: '<Root>/Sum2' incorporates:
   *  Constant: '<Root>/Vs_bias'
   *  Gain: '<S3>/K_ed'
   *  Gain: '<S3>/K_ep'
   *  Integrator: '<S3>/Integrator'
   *  Sum: '<S3>/Sum1'
   */
  p4p4_B.Sum2 = ((p4p4_P.K_ep * rtb_Backgain + p4p4_X.Integrator_CSTATE) -
                 p4p4_P.K_ed * p4p4_B.Gain1[5]) + p4p4_P.Vs_ff;
  if (rtmIsMajorTimeStep(p4p4_M)) {
  }

  /* Gain: '<S3>/K_ei' */
  p4p4_B.K_ei = p4p4_P.K_ei * rtb_Backgain;
  if (rtmIsMajorTimeStep(p4p4_M)) {
  }

  /* Derivative: '<S4>/Derivative' */
  if ((p4p4_DW.TimeStampA >= p4p4_M->Timing.t[0]) && (p4p4_DW.TimeStampB >=
       p4p4_M->Timing.t[0])) {
    rtb_Backgain = 0.0;
  } else {
    rtb_Backgain = p4p4_DW.TimeStampA;
    lastU = &p4p4_DW.LastUAtTimeA;
    if (p4p4_DW.TimeStampA < p4p4_DW.TimeStampB) {
      if (p4p4_DW.TimeStampB < p4p4_M->Timing.t[0]) {
        rtb_Backgain = p4p4_DW.TimeStampB;
        lastU = &p4p4_DW.LastUAtTimeB;
      }
    } else {
      if (p4p4_DW.TimeStampA >= p4p4_M->Timing.t[0]) {
        rtb_Backgain = p4p4_DW.TimeStampB;
        lastU = &p4p4_DW.LastUAtTimeB;
      }
    }

    rtb_Backgain = (p4p4_B.PitchCounttorad - *lastU) / (p4p4_M->Timing.t[0] -
      rtb_Backgain);
  }

  /* End of Derivative: '<S4>/Derivative' */

  /* Gain: '<S11>/Gain' */
  p4p4_B.Gain_l = p4p4_P.Gain_Gain_a1 * rtb_Backgain;
  if (rtmIsMajorTimeStep(p4p4_M)) {
  }

  /* Gain: '<S1>/Back gain' incorporates:
   *  Sum: '<S1>/Subtract'
   */
  rtb_Backgain = (p4p4_B.Sum2 - p4p4_B.Sum1) * p4p4_P.Backgain_Gain;

  /* Saturate: '<S4>/Back motor: Saturation' */
  if (rtb_Backgain > p4p4_P.BackmotorSaturation_UpperSat) {
    p4p4_B.BackmotorSaturation = p4p4_P.BackmotorSaturation_UpperSat;
  } else if (rtb_Backgain < p4p4_P.BackmotorSaturation_LowerSat) {
    p4p4_B.BackmotorSaturation = p4p4_P.BackmotorSaturation_LowerSat;
  } else {
    p4p4_B.BackmotorSaturation = rtb_Backgain;
  }

  /* End of Saturate: '<S4>/Back motor: Saturation' */
  if (rtmIsMajorTimeStep(p4p4_M)) {
  }

  /* Gain: '<S1>/Front gain' incorporates:
   *  Sum: '<S1>/Add'
   */
  rtb_Backgain = (p4p4_B.Sum1 + p4p4_B.Sum2) * p4p4_P.Frontgain_Gain;

  /* Saturate: '<S4>/Front motor: Saturation' */
  if (rtb_Backgain > p4p4_P.FrontmotorSaturation_UpperSat) {
    p4p4_B.FrontmotorSaturation = p4p4_P.FrontmotorSaturation_UpperSat;
  } else if (rtb_Backgain < p4p4_P.FrontmotorSaturation_LowerSat) {
    p4p4_B.FrontmotorSaturation = p4p4_P.FrontmotorSaturation_LowerSat;
  } else {
    p4p4_B.FrontmotorSaturation = rtb_Backgain;
  }

  /* End of Saturate: '<S4>/Front motor: Saturation' */
  if (rtmIsMajorTimeStep(p4p4_M)) {
    /* S-Function (hil_write_analog_block): '<S4>/HIL Write Analog' */

    /* S-Function Block: p4p4/Helicopter_interface/HIL Write Analog (hil_write_analog_block) */
    {
      t_error result;
      p4p4_DW.HILWriteAnalog_Buffer[0] = p4p4_B.FrontmotorSaturation;
      p4p4_DW.HILWriteAnalog_Buffer[1] = p4p4_B.BackmotorSaturation;
      result = hil_write_analog(p4p4_DW.HILInitialize_Card,
        p4p4_P.HILWriteAnalog_channels, 2, &p4p4_DW.HILWriteAnalog_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(p4p4_M, _rt_error_message);
      }
    }

    /* SignalConversion: '<S5>/TmpSignal ConversionAtTo FileInport1' */
    for (i = 0; i < 6; i++) {
      p4p4_B.TmpSignalConversionAtToFileInpo[i] = p4p4_B.xWorkspace1[i];
    }

    p4p4_B.TmpSignalConversionAtToFileInpo[6] = p4p4_B.Sum5;
    for (i = 0; i < 5; i++) {
      p4p4_B.TmpSignalConversionAtToFileInpo[i + 7] = p4p4_B.Gain1[i + 1];
    }

    /* End of SignalConversion: '<S5>/TmpSignal ConversionAtTo FileInport1' */

    /* ToFile: '<S5>/To File' */
    {
      if (!(++p4p4_DW.ToFile_IWORK.Decimation % 1) &&
          (p4p4_DW.ToFile_IWORK.Count*13)+1 < 100000000 ) {
        FILE *fp = (FILE *) p4p4_DW.ToFile_PWORK.FilePtr;
        if (fp != (NULL)) {
          real_T u[13];
          p4p4_DW.ToFile_IWORK.Decimation = 0;
          u[0] = p4p4_M->Timing.t[1];
          u[1] = p4p4_B.TmpSignalConversionAtToFileInpo[0];
          u[2] = p4p4_B.TmpSignalConversionAtToFileInpo[1];
          u[3] = p4p4_B.TmpSignalConversionAtToFileInpo[2];
          u[4] = p4p4_B.TmpSignalConversionAtToFileInpo[3];
          u[5] = p4p4_B.TmpSignalConversionAtToFileInpo[4];
          u[6] = p4p4_B.TmpSignalConversionAtToFileInpo[5];
          u[7] = p4p4_B.TmpSignalConversionAtToFileInpo[6];
          u[8] = p4p4_B.TmpSignalConversionAtToFileInpo[7];
          u[9] = p4p4_B.TmpSignalConversionAtToFileInpo[8];
          u[10] = p4p4_B.TmpSignalConversionAtToFileInpo[9];
          u[11] = p4p4_B.TmpSignalConversionAtToFileInpo[10];
          u[12] = p4p4_B.TmpSignalConversionAtToFileInpo[11];
          if (fwrite(u, sizeof(real_T), 13, fp) != 13) {
            rtmSetErrorStatus(p4p4_M,
                              "Error writing to MAT-file p4p4withFeedback_x.mat");
            return;
          }

          if (((++p4p4_DW.ToFile_IWORK.Count)*13)+1 >= 100000000) {
            (void)fprintf(stdout,
                          "*** The ToFile block will stop logging data before\n"
                          "    the simulation has ended, because it has reached\n"
                          "    the maximum number of elements (100000000)\n"
                          "    allowed in MAT-file p4p4withFeedback_x.mat.\n");
          }
        }
      }
    }

    /* SignalConversion: '<S5>/TmpSignal ConversionAtTo File1Inport1' */
    rtb_TmpSignalConversionAtToFile[0] = p4p4_B.uWorkspace[0];
    rtb_TmpSignalConversionAtToFile[2] = p4p4_B.Sum3[0];
    rtb_TmpSignalConversionAtToFile[1] = p4p4_B.uWorkspace[1];
    rtb_TmpSignalConversionAtToFile[3] = p4p4_B.Sum3[1];

    /* ToFile: '<S5>/To File1' */
    {
      if (!(++p4p4_DW.ToFile1_IWORK.Decimation % 1) &&
          (p4p4_DW.ToFile1_IWORK.Count*5)+1 < 100000000 ) {
        FILE *fp = (FILE *) p4p4_DW.ToFile1_PWORK.FilePtr;
        if (fp != (NULL)) {
          real_T u[5];
          p4p4_DW.ToFile1_IWORK.Decimation = 0;
          u[0] = p4p4_M->Timing.t[1];
          u[1] = rtb_TmpSignalConversionAtToFile[0];
          u[2] = rtb_TmpSignalConversionAtToFile[1];
          u[3] = rtb_TmpSignalConversionAtToFile[2];
          u[4] = rtb_TmpSignalConversionAtToFile[3];
          if (fwrite(u, sizeof(real_T), 5, fp) != 5) {
            rtmSetErrorStatus(p4p4_M,
                              "Error writing to MAT-file p4p4withFeedback_u.mat");
            return;
          }

          if (((++p4p4_DW.ToFile1_IWORK.Count)*5)+1 >= 100000000) {
            (void)fprintf(stdout,
                          "*** The ToFile block will stop logging data before\n"
                          "    the simulation has ended, because it has reached\n"
                          "    the maximum number of elements (100000000)\n"
                          "    allowed in MAT-file p4p4withFeedback_u.mat.\n");
          }
        }
      }
    }
  }
}

/* Model update function */
void p4p4_update(void)
{
  real_T *lastU;

  /* Update for Derivative: '<S4>/Derivative' */
  if (p4p4_DW.TimeStampA == (rtInf)) {
    p4p4_DW.TimeStampA = p4p4_M->Timing.t[0];
    lastU = &p4p4_DW.LastUAtTimeA;
  } else if (p4p4_DW.TimeStampB == (rtInf)) {
    p4p4_DW.TimeStampB = p4p4_M->Timing.t[0];
    lastU = &p4p4_DW.LastUAtTimeB;
  } else if (p4p4_DW.TimeStampA < p4p4_DW.TimeStampB) {
    p4p4_DW.TimeStampA = p4p4_M->Timing.t[0];
    lastU = &p4p4_DW.LastUAtTimeA;
  } else {
    p4p4_DW.TimeStampB = p4p4_M->Timing.t[0];
    lastU = &p4p4_DW.LastUAtTimeB;
  }

  *lastU = p4p4_B.PitchCounttorad;

  /* End of Update for Derivative: '<S4>/Derivative' */
  if (rtmIsMajorTimeStep(p4p4_M)) {
    rt_ertODEUpdateContinuousStates(&p4p4_M->solverInfo);
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
  if (!(++p4p4_M->Timing.clockTick0)) {
    ++p4p4_M->Timing.clockTickH0;
  }

  p4p4_M->Timing.t[0] = rtsiGetSolverStopTime(&p4p4_M->solverInfo);

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
    if (!(++p4p4_M->Timing.clockTick1)) {
      ++p4p4_M->Timing.clockTickH1;
    }

    p4p4_M->Timing.t[1] = p4p4_M->Timing.clockTick1 * p4p4_M->Timing.stepSize1 +
      p4p4_M->Timing.clockTickH1 * p4p4_M->Timing.stepSize1 * 4294967296.0;
  }
}

/* Derivatives for root system: '<Root>' */
void p4p4_derivatives(void)
{
  boolean_T lsat;
  boolean_T usat;
  XDot_p4p4_T *_rtXdot;
  _rtXdot = ((XDot_p4p4_T *) p4p4_M->ModelData.derivs);

  /* Derivatives for TransferFcn: '<S4>/Travel: Transfer Fcn' */
  _rtXdot->TravelTransferFcn_CSTATE = 0.0;
  _rtXdot->TravelTransferFcn_CSTATE += p4p4_P.TravelTransferFcn_A *
    p4p4_X.TravelTransferFcn_CSTATE;
  _rtXdot->TravelTransferFcn_CSTATE += p4p4_B.TravelCounttorad;

  /* Derivatives for TransferFcn: '<S4>/Pitch: Transfer Fcn' */
  _rtXdot->PitchTransferFcn_CSTATE = 0.0;
  _rtXdot->PitchTransferFcn_CSTATE += p4p4_P.PitchTransferFcn_A *
    p4p4_X.PitchTransferFcn_CSTATE;
  _rtXdot->PitchTransferFcn_CSTATE += p4p4_B.PitchCounttorad;

  /* Derivatives for TransferFcn: '<S4>/Elevation: Transfer Fcn' */
  _rtXdot->ElevationTransferFcn_CSTATE = 0.0;
  _rtXdot->ElevationTransferFcn_CSTATE += p4p4_P.ElevationTransferFcn_A *
    p4p4_X.ElevationTransferFcn_CSTATE;
  _rtXdot->ElevationTransferFcn_CSTATE += p4p4_B.ElevationCounttorad;

  /* Derivatives for Integrator: '<S3>/Integrator' */
  lsat = (p4p4_X.Integrator_CSTATE <= p4p4_P.Integrator_LowerSat);
  usat = (p4p4_X.Integrator_CSTATE >= p4p4_P.Integrator_UpperSat);
  if (((!lsat) && (!usat)) || (lsat && (p4p4_B.K_ei > 0.0)) || (usat &&
       (p4p4_B.K_ei < 0.0))) {
    _rtXdot->Integrator_CSTATE = p4p4_B.K_ei;
  } else {
    /* in saturation */
    _rtXdot->Integrator_CSTATE = 0.0;
  }

  /* End of Derivatives for Integrator: '<S3>/Integrator' */
}

/* Model initialize function */
void p4p4_initialize(void)
{
  /* Start for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: p4p4/HIL Initialize (hil_initialize_block) */
  {
    t_int result;
    t_boolean is_switching;
    result = hil_open("q8_usb", "0", &p4p4_DW.HILInitialize_Card);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(p4p4_M, _rt_error_message);
      return;
    }

    is_switching = false;
    result = hil_set_card_specific_options(p4p4_DW.HILInitialize_Card,
      "update_rate=normal;decimation=1", 32);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(p4p4_M, _rt_error_message);
      return;
    }

    result = hil_watchdog_clear(p4p4_DW.HILInitialize_Card);
    if (result < 0 && result != -QERR_HIL_WATCHDOG_CLEAR) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(p4p4_M, _rt_error_message);
      return;
    }

    if ((p4p4_P.HILInitialize_set_analog_input_ && !is_switching) ||
        (p4p4_P.HILInitialize_set_analog_inpu_m && is_switching)) {
      {
        int_T i1;
        real_T *dw_AIMinimums = &p4p4_DW.HILInitialize_AIMinimums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMinimums[i1] = p4p4_P.HILInitialize_analog_input_mini;
        }
      }

      {
        int_T i1;
        real_T *dw_AIMaximums = &p4p4_DW.HILInitialize_AIMaximums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMaximums[i1] = p4p4_P.HILInitialize_analog_input_maxi;
        }
      }

      result = hil_set_analog_input_ranges(p4p4_DW.HILInitialize_Card,
        p4p4_P.HILInitialize_analog_input_chan, 8U,
        &p4p4_DW.HILInitialize_AIMinimums[0], &p4p4_DW.HILInitialize_AIMaximums
        [0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(p4p4_M, _rt_error_message);
        return;
      }
    }

    if ((p4p4_P.HILInitialize_set_analog_output && !is_switching) ||
        (p4p4_P.HILInitialize_set_analog_outp_b && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOMinimums = &p4p4_DW.HILInitialize_AOMinimums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMinimums[i1] = p4p4_P.HILInitialize_analog_output_min;
        }
      }

      {
        int_T i1;
        real_T *dw_AOMaximums = &p4p4_DW.HILInitialize_AOMaximums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMaximums[i1] = p4p4_P.HILInitialize_analog_output_max;
        }
      }

      result = hil_set_analog_output_ranges(p4p4_DW.HILInitialize_Card,
        p4p4_P.HILInitialize_analog_output_cha, 8U,
        &p4p4_DW.HILInitialize_AOMinimums[0], &p4p4_DW.HILInitialize_AOMaximums
        [0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(p4p4_M, _rt_error_message);
        return;
      }
    }

    if ((p4p4_P.HILInitialize_set_analog_outp_e && !is_switching) ||
        (p4p4_P.HILInitialize_set_analog_outp_j && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &p4p4_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = p4p4_P.HILInitialize_initial_analog_ou;
        }
      }

      result = hil_write_analog(p4p4_DW.HILInitialize_Card,
        p4p4_P.HILInitialize_analog_output_cha, 8U,
        &p4p4_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(p4p4_M, _rt_error_message);
        return;
      }
    }

    if (p4p4_P.HILInitialize_set_analog_outp_p) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &p4p4_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = p4p4_P.HILInitialize_watchdog_analog_o;
        }
      }

      result = hil_watchdog_set_analog_expiration_state
        (p4p4_DW.HILInitialize_Card, p4p4_P.HILInitialize_analog_output_cha, 8U,
         &p4p4_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(p4p4_M, _rt_error_message);
        return;
      }
    }

    if ((p4p4_P.HILInitialize_set_encoder_param && !is_switching) ||
        (p4p4_P.HILInitialize_set_encoder_par_m && is_switching)) {
      {
        int_T i1;
        int32_T *dw_QuadratureModes = &p4p4_DW.HILInitialize_QuadratureModes[0];
        for (i1=0; i1 < 8; i1++) {
          dw_QuadratureModes[i1] = p4p4_P.HILInitialize_quadrature;
        }
      }

      result = hil_set_encoder_quadrature_mode(p4p4_DW.HILInitialize_Card,
        p4p4_P.HILInitialize_encoder_channels, 8U, (t_encoder_quadrature_mode *)
        &p4p4_DW.HILInitialize_QuadratureModes[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(p4p4_M, _rt_error_message);
        return;
      }
    }

    if ((p4p4_P.HILInitialize_set_encoder_count && !is_switching) ||
        (p4p4_P.HILInitialize_set_encoder_cou_k && is_switching)) {
      {
        int_T i1;
        int32_T *dw_InitialEICounts = &p4p4_DW.HILInitialize_InitialEICounts[0];
        for (i1=0; i1 < 8; i1++) {
          dw_InitialEICounts[i1] = p4p4_P.HILInitialize_initial_encoder_c;
        }
      }

      result = hil_set_encoder_counts(p4p4_DW.HILInitialize_Card,
        p4p4_P.HILInitialize_encoder_channels, 8U,
        &p4p4_DW.HILInitialize_InitialEICounts[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(p4p4_M, _rt_error_message);
        return;
      }
    }

    if ((p4p4_P.HILInitialize_set_pwm_params_at && !is_switching) ||
        (p4p4_P.HILInitialize_set_pwm_params__f && is_switching)) {
      uint32_T num_duty_cycle_modes = 0;
      uint32_T num_frequency_modes = 0;

      {
        int_T i1;
        int32_T *dw_POModeValues = &p4p4_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] = p4p4_P.HILInitialize_pwm_modes;
        }
      }

      result = hil_set_pwm_mode(p4p4_DW.HILInitialize_Card,
        p4p4_P.HILInitialize_pwm_channels, 8U, (t_pwm_mode *)
        &p4p4_DW.HILInitialize_POModeValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(p4p4_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        const uint32_T *p_HILInitialize_pwm_channels =
          p4p4_P.HILInitialize_pwm_channels;
        int32_T *dw_POModeValues = &p4p4_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          if (dw_POModeValues[i1] == PWM_DUTY_CYCLE_MODE || dw_POModeValues[i1] ==
              PWM_ONE_SHOT_MODE || dw_POModeValues[i1] == PWM_TIME_MODE) {
            p4p4_DW.HILInitialize_POSortedChans[num_duty_cycle_modes] =
              p_HILInitialize_pwm_channels[i1];
            p4p4_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes] =
              p4p4_P.HILInitialize_pwm_frequency;
            num_duty_cycle_modes++;
          } else {
            p4p4_DW.HILInitialize_POSortedChans[7U - num_frequency_modes] =
              p_HILInitialize_pwm_channels[i1];
            p4p4_DW.HILInitialize_POSortedFreqs[7U - num_frequency_modes] =
              p4p4_P.HILInitialize_pwm_frequency;
            num_frequency_modes++;
          }
        }
      }

      if (num_duty_cycle_modes > 0) {
        result = hil_set_pwm_frequency(p4p4_DW.HILInitialize_Card,
          &p4p4_DW.HILInitialize_POSortedChans[0], num_duty_cycle_modes,
          &p4p4_DW.HILInitialize_POSortedFreqs[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(p4p4_M, _rt_error_message);
          return;
        }
      }

      if (num_frequency_modes > 0) {
        result = hil_set_pwm_duty_cycle(p4p4_DW.HILInitialize_Card,
          &p4p4_DW.HILInitialize_POSortedChans[num_duty_cycle_modes],
          num_frequency_modes,
          &p4p4_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(p4p4_M, _rt_error_message);
          return;
        }
      }

      {
        int_T i1;
        int32_T *dw_POModeValues = &p4p4_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] = p4p4_P.HILInitialize_pwm_configuration;
        }
      }

      {
        int_T i1;
        int32_T *dw_POAlignValues = &p4p4_DW.HILInitialize_POAlignValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POAlignValues[i1] = p4p4_P.HILInitialize_pwm_alignment;
        }
      }

      {
        int_T i1;
        int32_T *dw_POPolarityVals = &p4p4_DW.HILInitialize_POPolarityVals[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POPolarityVals[i1] = p4p4_P.HILInitialize_pwm_polarity;
        }
      }

      result = hil_set_pwm_configuration(p4p4_DW.HILInitialize_Card,
        p4p4_P.HILInitialize_pwm_channels, 8U,
        (t_pwm_configuration *) &p4p4_DW.HILInitialize_POModeValues[0],
        (t_pwm_alignment *) &p4p4_DW.HILInitialize_POAlignValues[0],
        (t_pwm_polarity *) &p4p4_DW.HILInitialize_POPolarityVals[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(p4p4_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        real_T *dw_POSortedFreqs = &p4p4_DW.HILInitialize_POSortedFreqs[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POSortedFreqs[i1] = p4p4_P.HILInitialize_pwm_leading_deadb;
        }
      }

      {
        int_T i1;
        real_T *dw_POValues = &p4p4_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = p4p4_P.HILInitialize_pwm_trailing_dead;
        }
      }

      result = hil_set_pwm_deadband(p4p4_DW.HILInitialize_Card,
        p4p4_P.HILInitialize_pwm_channels, 8U,
        &p4p4_DW.HILInitialize_POSortedFreqs[0],
        &p4p4_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(p4p4_M, _rt_error_message);
        return;
      }
    }

    if ((p4p4_P.HILInitialize_set_pwm_outputs_a && !is_switching) ||
        (p4p4_P.HILInitialize_set_pwm_outputs_g && is_switching)) {
      {
        int_T i1;
        real_T *dw_POValues = &p4p4_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = p4p4_P.HILInitialize_initial_pwm_outpu;
        }
      }

      result = hil_write_pwm(p4p4_DW.HILInitialize_Card,
        p4p4_P.HILInitialize_pwm_channels, 8U, &p4p4_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(p4p4_M, _rt_error_message);
        return;
      }
    }

    if (p4p4_P.HILInitialize_set_pwm_outputs_o) {
      {
        int_T i1;
        real_T *dw_POValues = &p4p4_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = p4p4_P.HILInitialize_watchdog_pwm_outp;
        }
      }

      result = hil_watchdog_set_pwm_expiration_state(p4p4_DW.HILInitialize_Card,
        p4p4_P.HILInitialize_pwm_channels, 8U, &p4p4_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(p4p4_M, _rt_error_message);
        return;
      }
    }
  }

  /* Start for S-Function (hil_read_encoder_timebase_block): '<S4>/HIL Read Encoder Timebase' */

  /* S-Function Block: p4p4/Helicopter_interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
  {
    t_error result;
    result = hil_task_create_encoder_reader(p4p4_DW.HILInitialize_Card,
      p4p4_P.HILReadEncoderTimebase_samples_,
      p4p4_P.HILReadEncoderTimebase_channels, 3,
      &p4p4_DW.HILReadEncoderTimebase_Task);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(p4p4_M, _rt_error_message);
    }
  }

  /* Start for FromWorkspace: '<Root>/u Workspace' */
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
      25.25, 25.5, 25.75, 26.0 } ;

    static real_T pDataValues0[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.52359877559829882,
      0.52359877559829882, 0.52359877559829882, 0.52359877559829882,
      0.52359877559829882, 0.52359877559829882, 0.52359877559829882,
      0.52359877559829882, 0.52359877559829882, 0.52359877559829882,
      0.52359877559829882, 0.37449554142187541, 0.10869007190766065,
      -0.10914750605850169, -0.274264207029001, -0.39346467279040909,
      -0.47334738386700415, -0.52020634236314156, -0.52359877559829882,
      -0.52359877559829882, -0.518857701827955, -0.48725660071661991,
      -0.44678777539728914, -0.400660531791523, -0.35157519966897555,
      -0.30176391976365552, -0.25302792105199462, -0.20679534980941483,
      -0.16415626837481248, -0.12591420891716446, -0.092625955752187461,
      -0.064626271866578516, -0.042046867911532, -0.024825647300555292,
      -0.012679925309128488, -0.0050913986861696927, -0.0012401638476062615,
      1.8481290700824468E-7, -6.9644262422972682E-7, -3.4273446076956679E-7,
      -3.4273446076956679E-7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.10546426372090498,
      0.12241411881071285, 0.14124576121399854, 0.16184430581173498,
      0.18391624781679877, 0.2069053713815204, 0.22989374578282748,
      0.2514457238938, 0.26940616294048864, 0.28063063784413683,
      0.28063099091881855, 0.26309109586236373, 0.21924217185105446,
      0.13702419450585776, -1.6765096109064615E-7, -3.0980755120112962E-7,
      9.9353215556353818E-7, 3.5512375052871475E-7, 9.4716913389672854E-7,
      -2.2069417717136092E-7, 7.67288415702648E-7, -9.6427139228740528E-7,
      -4.2142333403541169E-7, 1.2283111267108221E-6, 6.7362834610691961E-7,
      2.6852007574875143E-7, -1.0722715185112349E-6, -1.0000757331445463E-6,
      -9.3503890860910028E-7, 4.1991653807692209E-7, 4.1309481278745285E-7,
      -4.1218892913562836E-7, -1.0335005810793461E-7, 6.8280994260927113E-7,
      -8.1918448720634959E-7, 4.9491689022008053E-7, -4.9353864228560787E-7,
      -4.5649820535068751E-9, -2.3042532010655019E-8, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } ;

    p4p4_DW.uWorkspace_PWORK.TimePtr = (void *) pTimeValues0;
    p4p4_DW.uWorkspace_PWORK.DataPtr = (void *) pDataValues0;
    p4p4_DW.uWorkspace_IWORK.PrevIndex = 0;
  }

  /* Start for FromWorkspace: '<Root>/x Workspace1' */
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
      25.25, 25.5, 25.75, 26.0 } ;

    static real_T pDataValues0[] = { 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1378421413625261, 3.1262155534579983,
      3.1033093000299643, 3.0666274151911783, 3.0144539223941584,
      2.9456562771175667, 2.8595077632935446, 2.7555515879651526,
      2.633505110490284, 2.4931956060320961, 2.334518576064299,
      2.158479340589889, 1.9681193138959088, 1.768040302233385,
      1.5635669874489087, 1.3600183002449817, 1.1622065440586762,
      0.97415192454381672, 0.79884524596198447, 0.63836028079063711,
      0.494039500388067, 0.36650594994397029, 0.25574900537962392,
      0.16123810186569873, 0.082041236058621042, 0.016936904398012511,
      -0.0354857095782299, -0.07674016078092788, -0.10837000927131045,
      -0.13188924674357652, -0.14873477925983722, -0.1602297813568,
      -0.1675568160988975, -0.17173979731424627, -0.17363430409488095,
      -0.17392609575128082, -0.17313828105499685, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, -0.015002048909068423, -0.046506351618112111, -0.091625013712135384,
      -0.14672753935514371, -0.20869397118807925, -0.27519058110636674,
      -0.34459405529608839, -0.41582470131356869, -0.48818590989947458,
      -0.56123801783275173, -0.63470811987118847, -0.70415694189763944,
      -0.76144010677592144, -0.80031604665009459, -0.81789325913790445,
      -0.81419474881570753, -0.79124702474522257, -0.7522184780594382,
      -0.70122671432732864, -0.64193986068538955, -0.57728312161028061,
      -0.51013420177638691, -0.44302777825738537, -0.37804361405570069,
      -0.3167874632283107, -0.26041732664243417, -0.20969045590496965,
      -0.16501780481079192, -0.12651939396153028, -0.094076949889064357,
      -0.06738213006504272, -0.045980008387851073, -0.029308138968390133,
      -0.016731924861395086, -0.0075780271225386652, -0.0011671666255994848,
      0.0031512587851357885, 0.0059622361413975049, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.10602875205865551, 0.22266037932317656, 0.31888147181640641,
      0.38944360631144165, 0.43795507377677839, 0.46997264230390062,
      0.4905172487754707, 0.5034310014147434, 0.51142138586029329,
      0.51630439857701826, 0.51925862127063682, 0.49083774996608226,
      0.40485553158962412, 0.27476029540700003, 0.12422902471846602,
      -0.026139658410362736, -0.16218574946180334, -0.27583886206738512,
      -0.36039031112058362, -0.41901683851591348, -0.4569691378736892,
      -0.47458291965510979, -0.47428257193178086, -0.45928325361707162,
      -0.43293507890238264, -0.398402596325597, -0.35851814860321485,
      -0.31572923641144979, -0.2720920644460218, -0.22929080418188622,
      -0.18866879114497118, -0.15126202204742989, -0.11783040568290433,
      -0.088883878160209664, -0.0646962531322506, -0.045309513535961583,
      -0.030520981496192211, -0.019866914376569383, -0.012621149941259584, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.42411500823462206, 0.46652650905808424,
      0.38488436997291947, 0.28224853798014093, 0.19404586986134692,
      0.128070274108489, 0.082178425886280437, 0.051655010557090583,
      0.031961537782199824, 0.019532050866899908, 0.011816890774474512,
      -0.11368348521821825, -0.34392887350583257, -0.52038094473049634,
      -0.602125082754136, -0.601474732515315, -0.54418436420576233,
      -0.45461245042232717, -0.33820579621279412, -0.23450610958131954,
      -0.15180919743110274, -0.0704551271256825, 0.0012013908933158472,
      0.059997273258836893, 0.10539269885875602, 0.13812993030714263,
      0.15953779088952857, 0.17115564876706038, 0.17454868786171202,
      0.1712050410565423, 0.16248805214766004, 0.14962707639016523,
      0.13372646545810227, 0.11578611009077866, 0.096750500111836307,
      0.077546958385156017, 0.059154128159077508, 0.042616268478491315,
      0.028983057741239195, 0.018990228804231898, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0042185705488361992, 0.011646277630566434, 0.021583989506211067,
      0.033554537758844594, 0.04723015704284822, 0.062369561958157425,
      0.078759748456942114, 0.096156906833638645, 0.11422105843899841,
      0.13243849864263421, 0.15002536006400849, 0.16580358080562224,
      0.17803918572207234, 0.18422937321995181, 0.18082191158375815,
      0.17140824728094184, 0.15852721197718794, 0.14394227510864796,
      0.12885026239520181, 0.1140373549350211, 0.099995630654641218,
      0.087009063318156776, 0.07521728083314709, 0.064661897941860072,
      0.055320003918895827, 0.0471284023282459, 0.040000598326239374,
      0.033838739828876294, 0.02854156339585253, 0.024009724739544741,
      0.020148975533718475, 0.016872120533083762, 0.014100044377351871,
      0.011761941174987089, 0.0097950447110946564, 0.0081444489824353244,
      0.0067622700152502463, 0.0056071844930425047, 0.0046436414574065688, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.016874282195344797, 0.029710828326920936,
      0.03975084750257854, 0.0478821930105341, 0.054702477136014491,
      0.060557619661236829, 0.065560745995138781, 0.069588633506786082,
      0.0722566064214391, 0.072869760814543175, 0.07034744568549714,
      0.063112882966455008, 0.04894241966580036, 0.024760749991517898,
      -0.013629846544774617, -0.037654657211265256, -0.051524141215015584,
      -0.058339747474159966, -0.060368050853784626, -0.059251629840722804,
      -0.056166897121519462, -0.051946269345937823, -0.047167129940038736,
      -0.042221531565148053, -0.037367576091856994, -0.032766406362599691,
      -0.02851121600802611, -0.024647433989452312, -0.021188705732095064,
      -0.018127354625231153, -0.015442996823305053, -0.013107420002538851,
      -0.011088304622927564, -0.0093524128094591239, -0.007867585855569726,
      -0.0066023829146373341, -0.0055287158687403115, -0.0046203420888309681,
      -0.0038541721425437419, -0.0032096528044130461, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } ;

    p4p4_DW.xWorkspace1_PWORK.TimePtr = (void *) pTimeValues0;
    p4p4_DW.xWorkspace1_PWORK.DataPtr = (void *) pDataValues0;
    p4p4_DW.xWorkspace1_IWORK.PrevIndex = 0;
  }

  /* Start for ToFile: '<S5>/To File' */
  {
    FILE *fp = (NULL);
    char fileName[509] = "p4p4withFeedback_x.mat";
    if ((fp = fopen(fileName, "wb")) == (NULL)) {
      rtmSetErrorStatus(p4p4_M,
                        "Error creating .mat file p4p4withFeedback_x.mat");
      return;
    }

    if (rt_WriteMat4FileHeader(fp,13,0,
         "allstatestrajectoryandmeasuredwithfeedback_x")) {
      rtmSetErrorStatus(p4p4_M,
                        "Error writing mat file header to file p4p4withFeedback_x.mat");
      return;
    }

    p4p4_DW.ToFile_IWORK.Count = 0;
    p4p4_DW.ToFile_IWORK.Decimation = -1;
    p4p4_DW.ToFile_PWORK.FilePtr = fp;
  }

  /* Start for ToFile: '<S5>/To File1' */
  {
    FILE *fp = (NULL);
    char fileName[509] = "p4p4withFeedback_u.mat";
    if ((fp = fopen(fileName, "wb")) == (NULL)) {
      rtmSetErrorStatus(p4p4_M,
                        "Error creating .mat file p4p4withFeedback_u.mat");
      return;
    }

    if (rt_WriteMat4FileHeader(fp,5,0,
         "allstatestrajectoryandmeasuredwithfeedback_u")) {
      rtmSetErrorStatus(p4p4_M,
                        "Error writing mat file header to file p4p4withFeedback_u.mat");
      return;
    }

    p4p4_DW.ToFile1_IWORK.Count = 0;
    p4p4_DW.ToFile1_IWORK.Decimation = -1;
    p4p4_DW.ToFile1_PWORK.FilePtr = fp;
  }

  /* InitializeConditions for TransferFcn: '<S4>/Travel: Transfer Fcn' */
  p4p4_X.TravelTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S4>/Pitch: Transfer Fcn' */
  p4p4_X.PitchTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S4>/Elevation: Transfer Fcn' */
  p4p4_X.ElevationTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for Integrator: '<S3>/Integrator' */
  p4p4_X.Integrator_CSTATE = p4p4_P.Integrator_IC;

  /* InitializeConditions for Derivative: '<S4>/Derivative' */
  p4p4_DW.TimeStampA = (rtInf);
  p4p4_DW.TimeStampB = (rtInf);
}

/* Model terminate function */
void p4p4_terminate(void)
{
  /* Terminate for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: p4p4/HIL Initialize (hil_initialize_block) */
  {
    t_boolean is_switching;
    t_int result;
    t_uint32 num_final_analog_outputs = 0;
    t_uint32 num_final_pwm_outputs = 0;
    hil_task_stop_all(p4p4_DW.HILInitialize_Card);
    hil_monitor_stop_all(p4p4_DW.HILInitialize_Card);
    is_switching = false;
    if ((p4p4_P.HILInitialize_set_analog_out_ex && !is_switching) ||
        (p4p4_P.HILInitialize_set_analog_outp_c && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &p4p4_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = p4p4_P.HILInitialize_final_analog_outp;
        }
      }

      num_final_analog_outputs = 8U;
    }

    if ((p4p4_P.HILInitialize_set_pwm_output_ap && !is_switching) ||
        (p4p4_P.HILInitialize_set_pwm_outputs_p && is_switching)) {
      {
        int_T i1;
        real_T *dw_POValues = &p4p4_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = p4p4_P.HILInitialize_final_pwm_outputs;
        }
      }

      num_final_pwm_outputs = 8U;
    }

    if (0
        || num_final_analog_outputs > 0
        || num_final_pwm_outputs > 0
        ) {
      /* Attempt to write the final outputs atomically (due to firmware issue in old Q2-USB). Otherwise write channels individually */
      result = hil_write(p4p4_DW.HILInitialize_Card
                         , p4p4_P.HILInitialize_analog_output_cha,
                         num_final_analog_outputs
                         , p4p4_P.HILInitialize_pwm_channels,
                         num_final_pwm_outputs
                         , NULL, 0
                         , NULL, 0
                         , &p4p4_DW.HILInitialize_AOVoltages[0]
                         , &p4p4_DW.HILInitialize_POValues[0]
                         , (t_boolean *) NULL
                         , NULL
                         );
      if (result == -QERR_HIL_WRITE_NOT_SUPPORTED) {
        t_error local_result;
        result = 0;

        /* The hil_write operation is not supported by this card. Write final outputs for each channel type */
        if (num_final_analog_outputs > 0) {
          local_result = hil_write_analog(p4p4_DW.HILInitialize_Card,
            p4p4_P.HILInitialize_analog_output_cha, num_final_analog_outputs,
            &p4p4_DW.HILInitialize_AOVoltages[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (num_final_pwm_outputs > 0) {
          local_result = hil_write_pwm(p4p4_DW.HILInitialize_Card,
            p4p4_P.HILInitialize_pwm_channels, num_final_pwm_outputs,
            &p4p4_DW.HILInitialize_POValues[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(p4p4_M, _rt_error_message);
        }
      }
    }

    hil_task_delete_all(p4p4_DW.HILInitialize_Card);
    hil_monitor_delete_all(p4p4_DW.HILInitialize_Card);
    hil_close(p4p4_DW.HILInitialize_Card);
    p4p4_DW.HILInitialize_Card = NULL;
  }

  /* Terminate for ToFile: '<S5>/To File' */
  {
    FILE *fp = (FILE *) p4p4_DW.ToFile_PWORK.FilePtr;
    if (fp != (NULL)) {
      char fileName[509] = "p4p4withFeedback_x.mat";
      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(p4p4_M,
                          "Error closing MAT-file p4p4withFeedback_x.mat");
        return;
      }

      if ((fp = fopen(fileName, "r+b")) == (NULL)) {
        rtmSetErrorStatus(p4p4_M,
                          "Error reopening MAT-file p4p4withFeedback_x.mat");
        return;
      }

      if (rt_WriteMat4FileHeader(fp, 13, p4p4_DW.ToFile_IWORK.Count,
           "allstatestrajectoryandmeasuredwithfeedback_x")) {
        rtmSetErrorStatus(p4p4_M,
                          "Error writing header for allstatestrajectoryandmeasuredwithfeedback_x to MAT-file p4p4withFeedback_x.mat");
      }

      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(p4p4_M,
                          "Error closing MAT-file p4p4withFeedback_x.mat");
        return;
      }

      p4p4_DW.ToFile_PWORK.FilePtr = (NULL);
    }
  }

  /* Terminate for ToFile: '<S5>/To File1' */
  {
    FILE *fp = (FILE *) p4p4_DW.ToFile1_PWORK.FilePtr;
    if (fp != (NULL)) {
      char fileName[509] = "p4p4withFeedback_u.mat";
      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(p4p4_M,
                          "Error closing MAT-file p4p4withFeedback_u.mat");
        return;
      }

      if ((fp = fopen(fileName, "r+b")) == (NULL)) {
        rtmSetErrorStatus(p4p4_M,
                          "Error reopening MAT-file p4p4withFeedback_u.mat");
        return;
      }

      if (rt_WriteMat4FileHeader(fp, 5, p4p4_DW.ToFile1_IWORK.Count,
           "allstatestrajectoryandmeasuredwithfeedback_u")) {
        rtmSetErrorStatus(p4p4_M,
                          "Error writing header for allstatestrajectoryandmeasuredwithfeedback_u to MAT-file p4p4withFeedback_u.mat");
      }

      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(p4p4_M,
                          "Error closing MAT-file p4p4withFeedback_u.mat");
        return;
      }

      p4p4_DW.ToFile1_PWORK.FilePtr = (NULL);
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
  p4p4_output();
  UNUSED_PARAMETER(tid);
}

void MdlUpdate(int_T tid)
{
  p4p4_update();
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
  p4p4_initialize();
}

void MdlTerminate(void)
{
  p4p4_terminate();
}

/* Registration function */
RT_MODEL_p4p4_T *p4p4(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* non-finite (run-time) assignments */
  p4p4_P.Integrator_UpperSat = rtInf;
  p4p4_P.Integrator_LowerSat = rtMinusInf;

  /* initialize real-time model */
  (void) memset((void *)p4p4_M, 0,
                sizeof(RT_MODEL_p4p4_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&p4p4_M->solverInfo, &p4p4_M->Timing.simTimeStep);
    rtsiSetTPtr(&p4p4_M->solverInfo, &rtmGetTPtr(p4p4_M));
    rtsiSetStepSizePtr(&p4p4_M->solverInfo, &p4p4_M->Timing.stepSize0);
    rtsiSetdXPtr(&p4p4_M->solverInfo, &p4p4_M->ModelData.derivs);
    rtsiSetContStatesPtr(&p4p4_M->solverInfo, (real_T **)
                         &p4p4_M->ModelData.contStates);
    rtsiSetNumContStatesPtr(&p4p4_M->solverInfo, &p4p4_M->Sizes.numContStates);
    rtsiSetNumPeriodicContStatesPtr(&p4p4_M->solverInfo,
      &p4p4_M->Sizes.numPeriodicContStates);
    rtsiSetPeriodicContStateIndicesPtr(&p4p4_M->solverInfo,
      &p4p4_M->ModelData.periodicContStateIndices);
    rtsiSetPeriodicContStateRangesPtr(&p4p4_M->solverInfo,
      &p4p4_M->ModelData.periodicContStateRanges);
    rtsiSetErrorStatusPtr(&p4p4_M->solverInfo, (&rtmGetErrorStatus(p4p4_M)));
    rtsiSetRTModelPtr(&p4p4_M->solverInfo, p4p4_M);
  }

  rtsiSetSimTimeStep(&p4p4_M->solverInfo, MAJOR_TIME_STEP);
  p4p4_M->ModelData.intgData.f[0] = p4p4_M->ModelData.odeF[0];
  p4p4_M->ModelData.contStates = ((real_T *) &p4p4_X);
  rtsiSetSolverData(&p4p4_M->solverInfo, (void *)&p4p4_M->ModelData.intgData);
  rtsiSetSolverName(&p4p4_M->solverInfo,"ode1");

  /* Initialize timing info */
  {
    int_T *mdlTsMap = p4p4_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    mdlTsMap[1] = 1;
    p4p4_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    p4p4_M->Timing.sampleTimes = (&p4p4_M->Timing.sampleTimesArray[0]);
    p4p4_M->Timing.offsetTimes = (&p4p4_M->Timing.offsetTimesArray[0]);

    /* task periods */
    p4p4_M->Timing.sampleTimes[0] = (0.0);
    p4p4_M->Timing.sampleTimes[1] = (0.002);

    /* task offsets */
    p4p4_M->Timing.offsetTimes[0] = (0.0);
    p4p4_M->Timing.offsetTimes[1] = (0.0);
  }

  rtmSetTPtr(p4p4_M, &p4p4_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = p4p4_M->Timing.sampleHitArray;
    mdlSampleHits[0] = 1;
    mdlSampleHits[1] = 1;
    p4p4_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(p4p4_M, -1);
  p4p4_M->Timing.stepSize0 = 0.002;
  p4p4_M->Timing.stepSize1 = 0.002;

  /* External mode info */
  p4p4_M->Sizes.checksums[0] = (271447400U);
  p4p4_M->Sizes.checksums[1] = (3133095282U);
  p4p4_M->Sizes.checksums[2] = (2221599985U);
  p4p4_M->Sizes.checksums[3] = (1586353086U);

  {
    static const sysRanDType rtAlwaysEnabled = SUBSYS_RAN_BC_ENABLE;
    static RTWExtModeInfo rt_ExtModeInfo;
    static const sysRanDType *systemRan[1];
    p4p4_M->extModeInfo = (&rt_ExtModeInfo);
    rteiSetSubSystemActiveVectorAddresses(&rt_ExtModeInfo, systemRan);
    systemRan[0] = &rtAlwaysEnabled;
    rteiSetModelMappingInfoPtr(p4p4_M->extModeInfo,
      &p4p4_M->SpecialInfo.mappingInfo);
    rteiSetChecksumsPtr(p4p4_M->extModeInfo, p4p4_M->Sizes.checksums);
    rteiSetTPtr(p4p4_M->extModeInfo, rtmGetTPtr(p4p4_M));
  }

  p4p4_M->solverInfoPtr = (&p4p4_M->solverInfo);
  p4p4_M->Timing.stepSize = (0.002);
  rtsiSetFixedStepSize(&p4p4_M->solverInfo, 0.002);
  rtsiSetSolverMode(&p4p4_M->solverInfo, SOLVER_MODE_SINGLETASKING);

  /* block I/O */
  p4p4_M->ModelData.blockIO = ((void *) &p4p4_B);

  {
    int32_T i;
    for (i = 0; i < 6; i++) {
      p4p4_B.xWorkspace1[i] = 0.0;
    }

    for (i = 0; i < 6; i++) {
      p4p4_B.Gain1[i] = 0.0;
    }

    for (i = 0; i < 12; i++) {
      p4p4_B.TmpSignalConversionAtToFileInpo[i] = 0.0;
    }

    p4p4_B.uWorkspace[0] = 0.0;
    p4p4_B.uWorkspace[1] = 0.0;
    p4p4_B.TravelCounttorad = 0.0;
    p4p4_B.Gain = 0.0;
    p4p4_B.Gain_d = 0.0;
    p4p4_B.PitchCounttorad = 0.0;
    p4p4_B.Gain_i = 0.0;
    p4p4_B.Gain_b = 0.0;
    p4p4_B.ElevationCounttorad = 0.0;
    p4p4_B.Gain_e = 0.0;
    p4p4_B.Sum = 0.0;
    p4p4_B.Gain_dg = 0.0;
    p4p4_B.Sum5 = 0.0;
    p4p4_B.Sum3[0] = 0.0;
    p4p4_B.Sum3[1] = 0.0;
    p4p4_B.Sum1 = 0.0;
    p4p4_B.Sum2 = 0.0;
    p4p4_B.K_ei = 0.0;
    p4p4_B.Gain_l = 0.0;
    p4p4_B.BackmotorSaturation = 0.0;
    p4p4_B.FrontmotorSaturation = 0.0;
  }

  /* parameters */
  p4p4_M->ModelData.defaultParam = ((real_T *)&p4p4_P);

  /* states (continuous) */
  {
    real_T *x = (real_T *) &p4p4_X;
    p4p4_M->ModelData.contStates = (x);
    (void) memset((void *)&p4p4_X, 0,
                  sizeof(X_p4p4_T));
  }

  /* states (dwork) */
  p4p4_M->ModelData.dwork = ((void *) &p4p4_DW);
  (void) memset((void *)&p4p4_DW, 0,
                sizeof(DW_p4p4_T));

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      p4p4_DW.HILInitialize_AIMinimums[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      p4p4_DW.HILInitialize_AIMaximums[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      p4p4_DW.HILInitialize_AOMinimums[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      p4p4_DW.HILInitialize_AOMaximums[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      p4p4_DW.HILInitialize_AOVoltages[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      p4p4_DW.HILInitialize_FilterFrequency[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      p4p4_DW.HILInitialize_POSortedFreqs[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      p4p4_DW.HILInitialize_POValues[i] = 0.0;
    }
  }

  p4p4_DW.TimeStampA = 0.0;
  p4p4_DW.LastUAtTimeA = 0.0;
  p4p4_DW.TimeStampB = 0.0;
  p4p4_DW.LastUAtTimeB = 0.0;
  p4p4_DW.HILWriteAnalog_Buffer[0] = 0.0;
  p4p4_DW.HILWriteAnalog_Buffer[1] = 0.0;

  /* data type transition information */
  {
    static DataTypeTransInfo dtInfo;
    (void) memset((char_T *) &dtInfo, 0,
                  sizeof(dtInfo));
    p4p4_M->SpecialInfo.mappingInfo = (&dtInfo);
    dtInfo.numDataTypes = 16;
    dtInfo.dataTypeSizes = &rtDataTypeSizes[0];
    dtInfo.dataTypeNames = &rtDataTypeNames[0];

    /* Block I/O transition table */
    dtInfo.B = &rtBTransTable;

    /* Parameters transition table */
    dtInfo.P = &rtPTransTable;
  }

  /* Initialize Sizes */
  p4p4_M->Sizes.numContStates = (4);   /* Number of continuous states */
  p4p4_M->Sizes.numPeriodicContStates = (0);/* Number of periodic continuous states */
  p4p4_M->Sizes.numY = (0);            /* Number of model outputs */
  p4p4_M->Sizes.numU = (0);            /* Number of model inputs */
  p4p4_M->Sizes.sysDirFeedThru = (0);  /* The model is not direct feedthrough */
  p4p4_M->Sizes.numSampTimes = (2);    /* Number of sample times */
  p4p4_M->Sizes.numBlocks = (69);      /* Number of blocks */
  p4p4_M->Sizes.numBlockIO = (22);     /* Number of block outputs */
  p4p4_M->Sizes.numBlockPrms = (154);  /* Sum of parameter "widths" */
  return p4p4_M;
}

/*========================================================================*
 * End of Classic call interface                                          *
 *========================================================================*/
