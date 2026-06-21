#ifndef SERDES_SIM_H
#define SERDES_SIM_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>

/* ── Global Parameters ───────────────────────────────────────────────── */
#define OSF             16          /* oversampling factor                */
#define DATA_RATE       60e9        /* symbol rate (Hz)                   */
#define FS              (OSF * DATA_RATE)
#define N_BIT           2048        /* 2^11 PRBS length                   */
                                    /* NOTE: CTLE sweep needs              */
                                    /* CTLE_NA * CTLE_NZ * CTLE_WINDOW     */
                                    /* symbols to finish (default 17500).  */
                                    /* Set N_BIT >= 20000 for full sweep.  */
#define ADC_BITS        5           /* ADC quantization bits              */
#define NUM_LEVELS      3        /* PAM-4                              */

/* TX FFE */
#define TX_FFE_PRE      0
#define TX_FFE_POST     0
#define TX_FFE_LEN      (TX_FFE_PRE + 1 + TX_FFE_POST)

/* RX FFE */
#define RX_FFE_PRE      3
#define RX_FFE_POST     10
#define RX_FFE_LEN      (RX_FFE_PRE + 1 + RX_FFE_POST)

/* DFE */
#define N_DFE           0

/* CDR */
#define LEN_CDR         1000

/* CTLE sweep */
#define CTLE_NA         7           /* # gain steps                       */
#define CTLE_NZ         5           /* # zero-frequency steps             */
#define CTLE_WINDOW     500         /* symbols per sweep point            */

/* Channel */
#define MAX_CHANNEL_TAPS 4096

/* Random Seed (-1 for current time) */
#define RANDOM_SEED 1

/* ── CTLE Structure ──────────────────────────────────────────────────── */
/*  HP: 1st-order Butterworth highpass  (2 b-coeffs, 2 a-coeffs, 1 state)
 *  LP: 2nd-order Butterworth lowpass   (3 b-coeffs, 3 a-coeffs, 2 states)
 */
typedef struct {
    /* highpass filter coefficients */
    double bhp[2];
    double ahp[2];      /* ahp[0] = 1 always */
    double zi_hp[1];     /* 1 state for 1st order */

    /* lowpass filter coefficients */
    double blp[3];
    double alp[3];      /* alp[0] = 1 always */
    double zi_lp[2];     /* 2 states for 2nd order */

    double A;            /* CTLE boost gain */
} CTLE;

/* ── Function Declarations ───────────────────────────────────────────── */

/* Butterworth filter design (matching MATLAB's butter) */
void butter1_highpass(double Wn, double b[2], double a[2]);
void butter2_lowpass (double Wn, double b[3], double a[3]);

/* CTLE design and single-sample step */
void   ctle_design(CTLE *ctle, double Fs, double zHz, double pHz, double A);
double ctle_step  (CTLE *ctle, double x);

/* Channel loading: reads FIR taps from a text file (one tap per line) */
int load_channel_taps(const char *filename, double *h_fir, int max_taps);

/* ADC quantisation (B-bit, range [-1, +1]) */
double adc_quantize(double x, int B);

/* Find the mode (most frequent value) in an integer array */
int int_mode(const int *arr, int n);

#endif /* SERDES_SIM_H */