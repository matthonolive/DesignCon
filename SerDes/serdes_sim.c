/*
 * serdes_sim.c
 *
 * C translation of a MATLAB SerDes link simulation:
 *   CDR → TX FFE → CTLE sweep → RX FFE + DFE
 *
 * Build:
 *   gcc -O2 -o serdes_sim serdes_sim.c -lm
 *
 * Usage:
 *   ./serdes_sim channel_taps.txt
 *
 *   channel_taps.txt: one FIR tap per line (ASCII, e.g. exported from MATLAB
 *   via  writematrix(h_fir(:), 'channel_taps.txt') ).
 *
 * NOTE: The MATLAB FitChannel() does S-parameter loading, interpolation,
 *       and IFFT. Those steps must be performed offline (in MATLAB, Python,
 *       etc.) to produce the FIR tap file consumed here.
 */

#include "serdes_sim.h"
#include "s4p_channel.h"
#include <ctype.h>

/* True if 'name' ends in .s4p / .s2p / .sNp (case-insensitive). */
static int is_touchstone(const char *name)
{
    size_t n = strlen(name);
    if (n < 4) return 0;
    const char *e = name + n - 4;
    return e[0] == '.' &&
           (e[1] == 's' || e[1] == 'S') &&
           isdigit((unsigned char)e[2]) &&
           (e[3] == 'p' || e[3] == 'P');
}

/* ═══════════════════════════════════════════════════════════════════════
 *  Utility helpers
 * ═══════════════════════════════════════════════════════════════════════ */

/* ── 1st-order Butterworth highpass ────────────────────────────────────
 *  Matches MATLAB: [b,a] = butter(1, Wn, 'high')
 *  Wn = f_cutoff / (Fs/2),  range (0,1)
 */
void butter1_highpass(double Wn, double b[2], double a[2])
{
    double C  = 1.0 / tan(M_PI * Wn / 2.0);
    double D  = 1.0 + C;
    b[0] =  C / D;
    b[1] = -C / D;
    a[0] =  1.0;
    a[1] =  (1.0 - C) / D;
}

/* ── 2nd-order Butterworth lowpass ─────────────────────────────────────
 *  Matches MATLAB: [b,a] = butter(2, Wn, 'low')
 */
void butter2_lowpass(double Wn, double b[3], double a[3])
{
    double C  = 1.0 / tan(M_PI * Wn / 2.0);
    double C2 = C * C;
    double S2 = sqrt(2.0);
    double D  = C2 + S2 * C + 1.0;

    b[0] = 1.0 / D;
    b[1] = 2.0 / D;
    b[2] = 1.0 / D;
    a[0] = 1.0;
    a[1] = 2.0 * (1.0 - C2) / D;
    a[2] = (C2 - S2 * C + 1.0) / D;
}

/* ── CTLE design ──────────────────────────────────────────────────────
 *  Equivalent to MATLAB ctleDesign(Fs, zHz, pHz, A)
 */
void ctle_design(CTLE *ctle, double Fs, double zHz, double pHz, double A)
{
    double Wn_hp = zHz / (Fs / 2.0);
    double Wn_lp = pHz / (Fs / 2.0);

    /* clamp to valid range */
    if (Wn_hp <= 0.0) Wn_hp = 1e-12;
    if (Wn_hp >= 1.0) Wn_hp = 1.0 - 1e-12;
    if (Wn_lp <= 0.0) Wn_lp = 1e-12;
    if (Wn_lp >= 1.0) Wn_lp = 1.0 - 1e-12;

    butter1_highpass(Wn_hp, ctle->bhp, ctle->ahp);
    butter2_lowpass (Wn_lp, ctle->blp, ctle->alp);
    ctle->A = A;

    /* zero filter states */
    ctle->zi_hp[0] = 0.0;
    ctle->zi_lp[0] = 0.0;
    ctle->zi_lp[1] = 0.0;
}

/* ── CTLE single-sample step ──────────────────────────────────────────
 *  Equivalent to MATLAB ctleStep(ctle, x)
 *
 *  Uses Direct-Form-II Transposed (same as MATLAB's `filter` with zi).
 *  For a filter H(z) = B(z)/A(z) with states zi:
 *      y     = b[0]*x + zi[0]
 *      zi[0] = b[1]*x - a[1]*y + zi[1]       (if order >= 2)
 *      zi[1] = b[2]*x - a[2]*y               (if order == 2)
 */
double ctle_step(CTLE *ctle, double x)
{
    double hp, v, y;

    /* ── highpass branch (1st order) ── */
    hp = ctle->bhp[0] * x + ctle->zi_hp[0];
    ctle->zi_hp[0] = ctle->bhp[1] * x - ctle->ahp[1] * hp;

    /* ── blend: v = x + A * hp ── */
    v = x + ctle->A * hp;

    /* ── lowpass branch (2nd order) ── */
    y = ctle->blp[0] * v + ctle->zi_lp[0];
    ctle->zi_lp[0] = ctle->blp[1] * v - ctle->alp[1] * y + ctle->zi_lp[1];
    ctle->zi_lp[1] = ctle->blp[2] * v - ctle->alp[2] * y;

    return y;
}

/* ── Load channel FIR taps from text file ─────────────────────────────
 *  Returns number of taps loaded (L).
 */
int load_channel_taps(const char *filename, double *h_fir, int max_taps)
{
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "ERROR: cannot open '%s'\n", filename);
        return -1;
    }
    int L = 0;
    while (L < max_taps && fscanf(fp, "%lf", &h_fir[L]) == 1)
        L++;
    fclose(fp);
    printf("Loaded %d channel taps from '%s'\n", L, filename);
    return L;
}

/* ── ADC quantisation ─────────────────────────────────────────────── */
double adc_quantize(double x, int B)
{
    double clamped = x;
    if (clamped >  1.0) clamped =  1.0;
    if (clamped < -1.0) clamped = -1.0;

    int levels = 1 << B;                                  /* 2^B */
    int xq = (int)round((clamped + 1.0) * (levels - 1) / 2.0);
    return xq * 2.0 / (levels - 1) - 1.0;
}

/* ── Mode of an integer array ─────────────────────────────────────── */
int int_mode(const int *arr, int n)
{
    if (n == 0) return 0;

    /* find range */
    int mn = arr[0], mx = arr[0];
    for (int i = 1; i < n; i++) {
        if (arr[i] < mn) mn = arr[i];
        if (arr[i] > mx) mx = arr[i];
    }
    int range = mx - mn + 1;
    int *counts = (int *)calloc(range, sizeof(int));
    for (int i = 0; i < n; i++)
        counts[arr[i] - mn]++;

    int best_val = mn, best_cnt = 0;
    for (int i = 0; i < range; i++) {
        if (counts[i] > best_cnt) {
            best_cnt = counts[i];
            best_val = mn + i;
        }
    }
    free(counts);
    return best_val;
}

void csv_write_array_header(FILE *fp, const char *prefix, int n)
{
    for (int i = 0; i < n; i++)
        fprintf(fp, ",%s[%d]", prefix, i);
}

void csv_write_array_row(FILE *fp, const double *x, int n)
{
    for (int i = 0; i < n; i++)
        fprintf(fp, ",%.8e", x[i]);
}

/* ═══════════════════════════════════════════════════════════════════════
 *  Main simulation
 * ═══════════════════════════════════════════════════════════════════════ */
int main(int argc, char *argv[])
{
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <channel_taps.txt> [tune_tx_ffe=0|1]\n", argv[0]);
        return 1;
    }

    int tune_tx_ffe = 0;
    if (argc >= 3) tune_tx_ffe = atoi(argv[1 + 1]);

    srand((RANDOM_SEED == -1) ? ( (unsigned)time(NULL) ) : RANDOM_SEED);

    /* ── Load channel ───────────────────────────────────────────────── */
    double h_fir[MAX_CHANNEL_TAPS];
    int L;
    if (is_touchstone(argv[1])) {
        S4pOpts opt; s4p_default_opts(&opt);
        L = channel_from_s4p(argv[1], h_fir, MAX_CHANNEL_TAPS, FS, &opt);
    } else {
        L = load_channel_taps(argv[1], h_fir, MAX_CHANNEL_TAPS);
    }
    if (L <= 0) return 1;

    /* ── Generate input bitstream (PAM-4) ──────────────────────────── */
    double *bits     = (double *)malloc(N_BIT * sizeof(double));
    double *bits_osf = (double *)malloc(N_BIT * OSF * sizeof(double));
    int    *data_clk_osf = (int *)malloc(2 * N_BIT * OSF * sizeof(int));

    for (int i = 0; i < N_BIT; i++) {
        int sym = rand() % NUM_LEVELS;                       /* 0..3 */
        bits[i] = (2.0 * sym - (NUM_LEVELS - 1)) / (NUM_LEVELS - 1);
        for (int j = 0; j < OSF; j++)
            bits_osf[i * OSF + j] = bits[i];
    }
    /* data clock: repeating pattern [0 0 … 0  1 1 … 1] each OSF samples,
       tiled N_BIT times → total length 2*N_BIT*OSF                         */
    for (int i = 0; i < N_BIT; i++) {
        for (int j = 0; j < OSF; j++) {
            data_clk_osf[2 * i * OSF + j]       = 0;
            data_clk_osf[(2 * i + 1) * OSF + j] = 1;
        }
    }

    /* ═══════════════════════════════════════════════════════════════════
     *  CDR
     * ═══════════════════════════════════════════════════════════════════ */
    printf("\n=== CDR ===\n");

    double *channel_buffer = (double *)calloc(L, sizeof(double));
    int     total_cdr = LEN_CDR * OSF;
    double *d_edge = (double *)malloc(total_cdr * sizeof(double));
    for (int i = 0; i < total_cdr; i++) d_edge[i] = NAN;

    double post_channel = 0.0;
    int    sample_instant = 0;
    int    lag = 0;

    for (int pt = 0; pt < total_cdr; pt++) {
        /* shift in new sample */
        memmove(channel_buffer + 1, channel_buffer, (L - 1) * sizeof(double));
        channel_buffer[0] = data_clk_osf[pt];

        double post_prev = post_channel;
        post_channel = 0.0;
        for (int k = 0; k < L; k++)
            post_channel += h_fir[k] * channel_buffer[k];

        d_edge[pt] = post_channel - post_prev;
    }

    /* find zero crossings of d_edge */
    int *cross_raw = (int *)malloc(total_cdr * sizeof(int));
    int  n_cross = 0;
    for (int i = 0; i < total_cdr - 1; i++) {
        if (d_edge[i] * d_edge[i + 1] <= 0.0) {
            int zc = (fabs(d_edge[i]) <= fabs(d_edge[i + 1])) ? i : i + 1;
            cross_raw[n_cross++] = zc;
        }
    }

    /* mode of (zero_crossing mod OSF) */
    int *cross_mod = (int *)malloc(n_cross * sizeof(int));
    for (int i = 0; i < n_cross; i++)
        cross_mod[i] = cross_raw[i] % OSF;
    sample_instant = int_mode(cross_mod, n_cross);

    /* lag: first index where d_edge > 0.8 * max(d_edge) */
    double d_max = -DBL_MAX;
    for (int i = 0; i < total_cdr; i++)
        if (!isnan(d_edge[i]) && d_edge[i] > d_max) d_max = d_edge[i];

    int first_high = 0;
    for (int i = 0; i < total_cdr; i++) {
        if (!isnan(d_edge[i]) && d_edge[i] > d_max * 0.8) {
            first_high = i;
            break;
        }
    }
    lag = first_high - (int)round((double)sample_instant / 2.0);

    printf("sample_instant = %d\n", sample_instant);
    printf("lag            = %d\n", lag);

    free(cross_raw);
    free(cross_mod);
    free(d_edge);
    free(channel_buffer);

    /* ═══════════════════════════════════════════════════════════════════
     *  TX FFE
     * ═══════════════════════════════════════════════════════════════════ */
    double TX_FFE[TX_FFE_LEN];
    memset(TX_FFE, 0, sizeof(TX_FFE));
    TX_FFE[TX_FFE_PRE] = 1.0;

    if (tune_tx_ffe) {
        printf("\n=== TX FFE Training ===\n");

        channel_buffer = (double *)calloc(L, sizeof(double));
        double mu = 0.01;

        /* temporary CTLE for TX FFE tuning */
        double ctle_z_tmp = 80e9, ctle_p_tmp = 100e9, ctle_A_tmp = 1.0;
        CTLE ctle_tmp;
        ctle_design(&ctle_tmp, FS, ctle_z_tmp, ctle_p_tmp, ctle_A_tmp);

        int N_tx_samp = (N_BIT - TX_FFE_POST) * OSF;

        FILE *fp_tx = fopen("tx_ffe_train_trace.csv", "w");
        if (!fp_tx) {
            fprintf(stderr, "WARNING: could not open tx_ffe_train_trace.csv\n");
        } else {
            fprintf(fp_tx,
                "pt,symbol_idx,osf_idx,tx_out,post_channel_raw,post_ctle,post_sampler,desired,bit_error,updated");
            csv_write_array_header(fp_tx, "TX_FFE", TX_FFE_LEN);
            fprintf(fp_tx, "\n");
        }

        for (int pt = 0; pt < N_tx_samp; pt++) {
            double tx_out;
            double post_channel_raw;
            double post_ctle;
            double post_sampler = NAN;
            double desired = NAN;
            double bit_error = NAN;
            int updated = 0;

            /* apply TX FFE */
            if (pt > TX_FFE_PRE * OSF) {
                tx_out = 0.0;
                for (int k = 0; k < TX_FFE_LEN; k++) {
                    int idx = pt + (k - TX_FFE_PRE) * OSF;
                    if (idx >= 0 && idx < N_BIT * OSF)
                        tx_out += TX_FFE[k] * bits_osf[idx];
                }
            } else {
                tx_out = bits_osf[pt];
            }

            /* apply channel */
            memmove(channel_buffer + 1, channel_buffer, (L - 1) * sizeof(double));
            channel_buffer[0] = tx_out;
            post_channel = 0.0;
            for (int k = 0; k < L; k++)
                post_channel += h_fir[k] * channel_buffer[k];

            post_channel_raw = post_channel;

            /* apply CTLE */
            post_ctle = ctle_step(&ctle_tmp, post_channel_raw);
            post_channel = post_ctle;

            /* sample at decision instant */
            if (pt % OSF == sample_instant &&
                (pt - lag - TX_FFE_PRE * OSF > 0)) {

                post_sampler = adc_quantize(post_channel, ADC_BITS);

                if (pt - lag > 0 && pt - lag < N_BIT * OSF) {
                    desired   = bits_osf[pt - lag];
                    bit_error = desired - post_sampler;

                    /* LMS update */
                    for (int k = 0; k < TX_FFE_LEN; k++) {
                        int idx = pt - lag + (k - TX_FFE_PRE) * OSF;
                        if (idx >= 0 && idx < N_BIT * OSF)
                            TX_FFE[k] += mu * bit_error * bits_osf[idx];
                    }
                    updated = 1;
                }
            }

            if (fp_tx) {
                fprintf(fp_tx, "%d,%d,%d,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%d",
                        pt, pt / OSF, pt % OSF,
                        tx_out, post_channel_raw, post_ctle,
                        post_sampler, desired, bit_error, updated);
                csv_write_array_row(fp_tx, TX_FFE, TX_FFE_LEN);
                fprintf(fp_tx, "\n");
            }
        }

        if (fp_tx) fclose(fp_tx);

        free(channel_buffer);

        printf("TX FFE taps:\n");
        for (int k = 0; k < TX_FFE_LEN; k++)
            printf("  TX_FFE[%d] = %+.6f\n", k, TX_FFE[k]);
    }

    /* ═══════════════════════════════════════════════════════════════════
     *  CTLE Sweep
     * ═══════════════════════════════════════════════════════════════════ */
    printf("\n=== CTLE Sweep ===\n");

    double ctle_p = 100e9;
    double ctle_A_max = 2.0;

    double A_vec[CTLE_NA], z_vec[CTLE_NZ];
    for (int i = 0; i < CTLE_NA; i++)
        A_vec[i] = 1.0 + i * (ctle_A_max - 1.0) / (CTLE_NA - 1);
    for (int i = 0; i < CTLE_NZ; i++)
        z_vec[i] = ctle_p / 2.0 + i * (ctle_p / 2.0) / (CTLE_NZ - 1);

    double ctle_z = z_vec[0];
    double ctle_A = A_vec[0];
    CTLE   ctle;
    ctle_design(&ctle, FS, ctle_z, ctle_p, ctle_A);

    channel_buffer = (double *)calloc(L, sizeof(double));

    int    ia = 0, iz = 0;
    int    ctle_cnt = 0;
    double err_acc  = 0.0;
    int    ctle_train_done = 0;

    double J[CTLE_NA][CTLE_NZ];
    for (int a = 0; a < CTLE_NA; a++)
        for (int z = 0; z < CTLE_NZ; z++)
            J[a][z] = 1e30;

    int N_ctle_samp = (N_BIT - TX_FFE_POST) * OSF;

    FILE *fp_ctle = fopen("ctle_sweep_trace.csv", "w");
    if (!fp_ctle) {
        fprintf(stderr, "WARNING: could not open ctle_sweep_trace.csv\n");
    } else {
        fprintf(fp_ctle,
            "pt,symbol_idx,osf_idx,ia,iz,ctle_A,ctle_z,post_channel_raw,post_ctle,desired,bit_error,ctle_cnt,err_acc,ctle_train_done");
        fprintf(fp_ctle, "\n");
    }

    for (int pt = 0; pt < N_ctle_samp; pt++) {
        double tx_out;
        double post_channel_raw;
        double post_ctle;
        double desired = NAN;
        double bit_error = NAN;

        /* apply TX FFE */
        if (pt > TX_FFE_PRE * OSF) {
            tx_out = 0.0;
            for (int k = 0; k < TX_FFE_LEN; k++) {
                int idx = pt + (k - TX_FFE_PRE) * OSF;
                if (idx >= 0 && idx < N_BIT * OSF)
                    tx_out += TX_FFE[k] * bits_osf[idx];
            }
        } else {
            tx_out = bits_osf[pt];
        }

        /* apply channel */
        memmove(channel_buffer + 1, channel_buffer, (L - 1) * sizeof(double));
        channel_buffer[0] = tx_out;
        post_channel = 0.0;
        for (int k = 0; k < L; k++)
            post_channel += h_fir[k] * channel_buffer[k];

        post_channel_raw = post_channel;

        /* apply CTLE */
        post_ctle = ctle_step(&ctle, post_channel_raw);
        post_channel = post_ctle;

        /* sample at decision instant */
        if (pt % OSF == sample_instant &&
            (pt - lag - TX_FFE_PRE * OSF > 0)) {

            double post_sampler = post_channel;

            if (pt - lag > 0 && pt - lag < N_BIT * OSF) {
                desired   = bits_osf[pt - lag];
                bit_error = desired - post_sampler;

                if (!ctle_train_done) {
                    ctle_cnt++;
                    err_acc += bit_error * bit_error;

                    if (ctle_cnt == CTLE_WINDOW) {
                        J[ia][iz] = err_acc / CTLE_WINDOW;
                        ctle_cnt = 0;
                        err_acc  = 0.0;

                        ia++;
                        if (ia >= CTLE_NA) {
                            ia = 0;
                            iz++;
                        }
                        if (iz >= CTLE_NZ) {
                            double best_J = 1e30;
                            int ia_best = 0, iz_best = 0;
                            for (int a = 0; a < CTLE_NA; a++)
                                for (int z = 0; z < CTLE_NZ; z++)
                                    if (J[a][z] < best_J) {
                                        best_J  = J[a][z];
                                        ia_best = a;
                                        iz_best = z;
                                    }
                            ctle_A = A_vec[ia_best];
                            ctle_z = z_vec[iz_best];
                            ctle_design(&ctle, FS, ctle_z, ctle_p, ctle_A);
                            ctle_train_done = 1;
                        } else {
                            ctle_A = A_vec[ia];
                            ctle_z = z_vec[iz];
                            ctle_design(&ctle, FS, ctle_z, ctle_p, ctle_A);
                        }
                    }
                }
            }
        }

        if (fp_ctle) {
            fprintf(fp_ctle,
                "%d,%d,%d,%d,%d,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%d,%.8e,%d\n",
                pt, pt / OSF, pt % OSF,
                ia, iz, ctle_A, ctle_z,
                post_channel_raw, post_ctle,
                desired, bit_error,
                ctle_cnt, err_acc, ctle_train_done);
        }
    }

    if (fp_ctle) fclose(fp_ctle);
    free(channel_buffer);

    printf("ctle_z = %.4e Hz\n", ctle_z);
    printf("ctle_A = %.4f\n",    ctle_A);

    /* ═══════════════════════════════════════════════════════════════════
     *  RX FFE + DFE  —  PASS 1: TRAINING  (taps adapt, no eye capture)
     * ═══════════════════════════════════════════════════════════════════ */
    printf("\n=== RX FFE + DFE Training ===\n");

    double RX_FFE[RX_FFE_LEN];
    memset(RX_FFE, 0, sizeof(RX_FFE));
    RX_FFE[RX_FFE_PRE] = 1.0;

    double DFE[N_DFE];
    memset(DFE, 0, sizeof(DFE));
    double d_hist[N_DFE];
    memset(d_hist, 0, sizeof(d_hist));
    int en_DFE = 1;

    double rx_buffer[RX_FFE_LEN];
    memset(rx_buffer, 0, sizeof(rx_buffer));

    channel_buffer = (double *)calloc(L, sizeof(double));

    double mu_ffe = 0.01;
    double mu_dfe = 0.005;

    int N_rx_samp = (N_BIT - TX_FFE_POST) * OSF;

    /* re-init CTLE with final parameters */
    ctle_design(&ctle, FS, ctle_z, ctle_p, ctle_A);

    FILE *fp_rx = fopen("rx_ffe_dfe_train_trace.csv", "w");
    if (!fp_rx) {
        fprintf(stderr, "WARNING: could not open rx_ffe_dfe_train_trace.csv\n");
    } else {
        fprintf(fp_rx,
            "pt,symbol_idx,osf_idx,post_adc,y_ffe,y_total,d_hat,desired,bit_error,updated");
        csv_write_array_header(fp_rx, "RX_FFE", RX_FFE_LEN);
        csv_write_array_header(fp_rx, "DFE", N_DFE);
        fprintf(fp_rx, "\n");
    }

        for (int pt = 0; pt < N_rx_samp; pt++) {
        double tx_out;
        double y_ffe_all = NAN;
        double y_all = NAN;
        double d_hat = NAN;
        double desired = NAN;
        double bit_error = NAN;
        int updated = 0;

        /* apply TX FFE */
        if (pt > TX_FFE_PRE * OSF) {
            tx_out = 0.0;
            for (int k = 0; k < TX_FFE_LEN; k++) {
                int idx = pt + (k - TX_FFE_PRE) * OSF;
                if (idx >= 0 && idx < N_BIT * OSF)
                    tx_out += TX_FFE[k] * bits_osf[idx];
            }
        } else {
            tx_out = bits_osf[pt];
        }

        /* apply channel */
        memmove(channel_buffer + 1, channel_buffer, (L - 1) * sizeof(double));
        channel_buffer[0] = tx_out;
        post_channel = 0.0;
        for (int k = 0; k < L; k++)
            post_channel += h_fir[k] * channel_buffer[k];

        /* apply CTLE */
        post_channel = ctle_step(&ctle, post_channel);

        /* ADC quantise */
        post_channel = adc_quantize(post_channel, ADC_BITS);

        /* update RX FFE input buffer */
        memmove(rx_buffer + 1, rx_buffer, (RX_FFE_LEN - 1) * sizeof(double));
        rx_buffer[0] = post_channel;

        /* compute outputs every timestep for logging */
        y_ffe_all = 0.0;
        for (int k = 0; k < RX_FFE_LEN; k++)
            y_ffe_all += RX_FFE[k] * rx_buffer[k];

        y_all = y_ffe_all;
        if (en_DFE) {
            for (int k = 0; k < N_DFE; k++)
                y_all -= DFE[k] * d_hist[k];
        }

        /* sample at decision instant */
        if (pt % OSF == sample_instant && (pt - lag > 0)) {

            double y_ffe = y_ffe_all;
            double y = y_all;

            /* slicer */
            d_hat = (y < 0.0) ? -1.0 : 1.0;

            int lag_idx = pt - lag;
            if (lag_idx >= 0 && lag_idx < N_BIT * OSF) {
                desired   = bits_osf[lag_idx];
                bit_error = desired - y;

                for (int k = 0; k < RX_FFE_LEN; k++)
                    RX_FFE[k] += mu_ffe * bit_error * rx_buffer[k];

                if (en_DFE) {
                    for (int k = 0; k < N_DFE; k++)
                        DFE[k] -= mu_dfe * bit_error * d_hist[k];
                }

                updated = 1;
            }

            if (N_DFE > 0)
                memmove(d_hist + 1, d_hist, (N_DFE - 1) * sizeof(double));

            d_hist[0] = d_hat;
        }

        if (fp_rx) {
            fprintf(fp_rx, "%d,%d,%d,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%d",
                    pt, pt / OSF, pt % OSF,
                    post_channel, y_ffe_all, y_all,
                    d_hat, desired, bit_error, updated);
            csv_write_array_row(fp_rx, RX_FFE, RX_FFE_LEN);
            csv_write_array_row(fp_rx, DFE, N_DFE);
            fprintf(fp_rx, "\n");
        }
    }
    if (fp_rx) fclose(fp_rx);
    free(channel_buffer);

    /* ── Print trained tap values ──────────────────────────────────── */
    printf("\nRX FFE taps (after training):\n");
    for (int k = 0; k < RX_FFE_LEN; k++)
        printf("  RX_FFE[%2d] = %+.6f\n", k, RX_FFE[k]);
    printf("\nDFE taps (after training):\n");
    for (int k = 0; k < N_DFE; k++)
        printf("  DFE[%d] = %+.6f\n", k, DFE[k]);

    /* ═══════════════════════════════════════════════════════════════════
     *  RX FFE + DFE  —  PASS 2: CAPTURE  (taps frozen, record eye)
     * ═══════════════════════════════════════════════════════════════════ */
    printf("\n=== Eye Diagram Capture (taps frozen) ===\n");

    /* re-init all stateful blocks so the signal path is clean */
    channel_buffer = (double *)calloc(L, sizeof(double));
    memset(rx_buffer, 0, sizeof(rx_buffer));
    memset(d_hist,    0, sizeof(d_hist));
    ctle_design(&ctle, FS, ctle_z, ctle_p, ctle_A);

    /* eye trace buffers */
    double *yffe_trace = (double *)malloc(N_rx_samp * sizeof(double));
    double *ydfe_trace = (double *)malloc(N_rx_samp * sizeof(double));
    for (int i = 0; i < N_rx_samp; i++) {
        yffe_trace[i] = NAN;
        ydfe_trace[i] = NAN;
    }

    for (int pt = 0; pt < N_rx_samp; pt++) {
        /* apply TX FFE (frozen) */
        double tx_out;
        if (pt > TX_FFE_PRE * OSF) {
            tx_out = 0.0;
            for (int k = 0; k < TX_FFE_LEN; k++) {
                int idx = pt + (k - TX_FFE_PRE) * OSF;
                if (idx >= 0 && idx < N_BIT * OSF)
                    tx_out += TX_FFE[k] * bits_osf[idx];
            }
        } else {
            tx_out = bits_osf[pt];
        }

        /* apply channel */
        memmove(channel_buffer + 1, channel_buffer, (L - 1) * sizeof(double));
        channel_buffer[0] = tx_out;
        post_channel = 0.0;
        for (int k = 0; k < L; k++)
            post_channel += h_fir[k] * channel_buffer[k];

        /* apply CTLE (frozen) */
        post_channel = ctle_step(&ctle, post_channel);

        /* ADC quantise */
        post_channel = adc_quantize(post_channel, ADC_BITS);

        /* update RX FFE input buffer */
        memmove(rx_buffer + 1, rx_buffer, (RX_FFE_LEN - 1) * sizeof(double));
        rx_buffer[0] = post_channel;

        /* ── compute eye traces (frozen taps) ── */
        {
            double yf = 0.0;
            for (int k = 0; k < RX_FFE_LEN; k++)
                yf += RX_FFE[k] * rx_buffer[k];
            yffe_trace[pt] = yf;

            if (en_DFE) {
                double yd = yf;
                for (int k = 0; k < N_DFE; k++)
                    yd -= DFE[k] * d_hist[k];
                ydfe_trace[pt] = yd;
            } else {
                ydfe_trace[pt] = yf;
            }
        }

        /* sample at decision instant — slicer only, NO tap updates */
        if (pt % OSF == sample_instant && (pt - lag > 0)) {
            double y_ffe = 0.0;
            for (int k = 0; k < RX_FFE_LEN; k++)
                y_ffe += RX_FFE[k] * rx_buffer[k];

            double y = y_ffe;
            if (en_DFE) {
                for (int k = 0; k < N_DFE; k++)
                    y -= DFE[k] * d_hist[k];
            }

            double d_hat = (y < 0.0) ? -1.0 : 1.0;

            if (N_DFE > 0) 
                memmove(d_hist + 1, d_hist, (N_DFE - 1) * sizeof(double));

            d_hist[0] = d_hat;
        }
    }
    free(channel_buffer);

    /* ═══════════════════════════════════════════════════════════════════
     *  Write eye diagram CSVs
     * ═══════════════════════════════════════════════════════════════════ */
    {
        /* skip a few symbols for filter settling */
        int skip_symbols = 50;
        int start_pt     = skip_symbols * OSF;
        int n_symbols    = (N_rx_samp - start_pt) / OSF;

        const char *ffe_csv = "eye_ffe.csv";
        const char *dfe_csv = "eye_dfe.csv";

        FILE *fp_ffe = fopen(ffe_csv, "w");
        FILE *fp_dfe = fopen(dfe_csv, "w");
        if (!fp_ffe || !fp_dfe) {
            fprintf(stderr, "ERROR: cannot open eye CSV files for writing\n");
        } else {
            for (int row = 0; row < OSF; row++) {
                for (int col = 0; col < n_symbols; col++) {
                    int idx = start_pt + col * OSF + row;
                    if (col > 0) { fprintf(fp_ffe, ","); fprintf(fp_dfe, ","); }
                    fprintf(fp_ffe, "%.8e", yffe_trace[idx]);
                    fprintf(fp_dfe, "%.8e", ydfe_trace[idx]);
                }
                fprintf(fp_ffe, "\n");
                fprintf(fp_dfe, "\n");
            }
            fclose(fp_ffe);
            fclose(fp_dfe);
            printf("\nEye diagram data written to '%s' and '%s'\n", ffe_csv, dfe_csv);
            printf("  (%d symbols x %d samples/UI)\n", n_symbols, OSF);
        }
    }

    free(yffe_trace);
    free(ydfe_trace);

    /* ── Cleanup ───────────────────────────────────────────────────── */
    free(bits);
    free(bits_osf);
    free(data_clk_osf);

    printf("\nDone.\n");
    return 0;
}