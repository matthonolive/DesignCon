/*
 * s4p_channel.c  —  Touchstone .s4p -> channel FIR taps (no external deps)
 *
 * Pipeline:
 *   1. Parse Touchstone option line (units / format / Rref) and all data,
 *      tolerant of the line-wrapping used by 4-port files.
 *   2. Form differential through response Sdd21(f) by mixed-mode combination.
 *   3. Interpolate Sdd21 onto a uniform half-spectrum at df = Fs/NFFT,
 *      with DC extrapolation and a smooth raised-cosine band-edge taper.
 *   4. Enforce conjugate symmetry and IFFT (radix-2) to a real h[n] at Fs.
 *   5. Trim leading bulk delay + trailing decay, return contiguous taps.
 */

#include "s4p_channel.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ───────────────────────── small complex helpers ───────────────────────── */
typedef struct { double re, im; } cplx;

static cplx c_add(cplx a, cplx b){ return (cplx){a.re+b.re, a.im+b.im}; }
static cplx c_sub(cplx a, cplx b){ return (cplx){a.re-b.re, a.im-b.im}; }
static cplx c_scale(cplx a, double s){ return (cplx){a.re*s, a.im*s}; }
static double c_abs(cplx a){ return sqrt(a.re*a.re + a.im*a.im); }

/* ───────────────────────── radix-2 in-place FFT ────────────────────────── */
/* dir = -1 forward, +1 inverse (inverse includes the 1/n scaling).          */
static void fft_radix2(double *re, double *im, int n, int dir)
{
    /* bit-reversal permutation */
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) {
            double tr = re[i]; re[i] = re[j]; re[j] = tr;
            double ti = im[i]; im[i] = im[j]; im[j] = ti;
        }
    }
    for (int len = 2; len <= n; len <<= 1) {
        double ang = dir * 2.0 * M_PI / len;
        double wr = cos(ang), wi = sin(ang);
        for (int i = 0; i < n; i += len) {
            double cur_r = 1.0, cur_i = 0.0;
            for (int k = 0; k < len / 2; k++) {
                int a = i + k, b = i + k + len / 2;
                double ur = re[a],            ui = im[a];
                double vr = re[b]*cur_r - im[b]*cur_i;
                double vi = re[b]*cur_i + im[b]*cur_r;
                re[a] = ur + vr; im[a] = ui + vi;
                re[b] = ur - vr; im[b] = ui - vi;
                double nr = cur_r*wr - cur_i*wi;
                cur_i     = cur_r*wi + cur_i*wr;
                cur_r     = nr;
            }
        }
    }
    if (dir > 0)
        for (int i = 0; i < n; i++) { re[i] /= n; im[i] /= n; }
}

/* ───────────────────────── Touchstone parser ───────────────────────────── */
/* Reads freq[] and the full 4x4 S-matrix per point into s[16][Np] (cplx).
 * Returns Np (>0) or -1. Caller frees *freq and the 16 column arrays via the
 * single allocation pointed to by *blob.                                     */
typedef struct {
    double *freq;        /* [Np]                       */
    cplx   *S;           /* [Np*16], row-major per pt  */
    int     Np;
    double  Rref;
} TouchstoneData;

static double freq_scale(const char *unit)
{
    char u[8] = {0}; int j = 0;
    for (int i = 0; unit[i] && j < 7; i++) u[j++] = (char)tolower((unsigned char)unit[i]);
    if (!strcmp(u, "hz"))  return 1.0;
    if (!strcmp(u, "khz")) return 1e3;
    if (!strcmp(u, "mhz")) return 1e6;
    if (!strcmp(u, "ghz")) return 1e9;
    return 1.0;
}

static void free_touchstone(TouchstoneData *td)
{
    if (!td) return;
    free(td->freq);
    free(td->S);
    td->freq = NULL; td->S = NULL; td->Np = 0;
}

static int parse_touchstone(const char *filename, TouchstoneData *td)
{
    FILE *fp = fopen(filename, "r");
    if (!fp) { fprintf(stderr, "ERROR: cannot open '%s'\n", filename); return -1; }

    double fscale = 1.0;     /* default Hz                          */
    char   fmt[4] = "ma";    /* ma | db | ri                        */
    double Rref   = 50.0;
    int    have_opt = 0;

    /* token stream of all numeric data (option/comment lines skipped) */
    size_t cap = 1 << 16, ndat = 0;
    double *dat = (double *)malloc(cap * sizeof(double));
    if (!dat) { fclose(fp); return -1; }

    char line[8192];
    while (fgets(line, sizeof line, fp)) {
        /* strip trailing CR/LF */
        char *p = line;
        while (*p == ' ' || *p == '\t') p++;
        if (*p == '!' ) continue;                 /* comment */
        if (*p == '#') {                          /* option line */
            if (!have_opt) {
                /* format: # <funit> <type> <fmt> R <Rref> */
                char funit[16] = "hz", type[8] = "s", f2[8] = "ma", rtag[8] = "";
                double r = 50.0;
                int got = sscanf(p + 1, " %15s %7s %7s %7s %lf",
                                 funit, type, f2, rtag, &r);
                if (got >= 1) fscale = freq_scale(funit);
                if (got >= 3) { fmt[0]=f2[0]; fmt[1]=f2[1]; fmt[2]=0; }
                if (got >= 5 && (rtag[0]=='R'||rtag[0]=='r')) Rref = r;
                have_opt = 1;
            }
            continue;
        }
        if (*p == '\0' || *p == '\n' || *p == '\r') continue;  /* blank */

        /* tokenize numbers on this line */
        char *tok = strtok(p, " \t\r\n");
        while (tok) {
            char *end; double v = strtod(tok, &end);
            if (end != tok) {
                if (ndat == cap) {
                    cap <<= 1;
                    double *nd = (double *)realloc(dat, cap * sizeof(double));
                    if (!nd) { free(dat); fclose(fp); return -1; }
                    dat = nd;
                }
                dat[ndat++] = v;
            }
            tok = strtok(NULL, " \t\r\n");
        }
    }
    fclose(fp);

    /* 4-port: each point = 1 freq + 16 complex pairs = 33 numbers */
    const int PER = 1 + 16 * 2;
    if (ndat % PER != 0 || ndat == 0) {
        fprintf(stderr, "ERROR: %zu data tokens not a multiple of %d "
                        "(is this a 4-port file?)\n", ndat, PER);
        free(dat); return -1;
    }
    int Np = (int)(ndat / PER);

    double *freq = (double *)malloc((size_t)Np * sizeof(double));
    cplx   *S    = (cplx  *)malloc((size_t)Np * 16 * sizeof(cplx));
    if (!freq || !S) { free(dat); free(freq); free(S); return -1; }

    int ri = (fmt[0]=='r' || fmt[0]=='R');      /* real/imag         */
    int db = (fmt[0]=='d' || fmt[0]=='D');      /* dB/angle          */

    for (int n = 0; n < Np; n++) {
        const double *rec = dat + (size_t)n * PER;
        freq[n] = rec[0] * fscale;
        for (int k = 0; k < 16; k++) {
            double a = rec[1 + 2*k], b = rec[2 + 2*k];
            cplx z;
            if (ri) {
                z.re = a; z.im = b;
            } else {
                double mag = db ? pow(10.0, a / 20.0) : a;
                double ang = b * M_PI / 180.0;
                z.re = mag * cos(ang);
                z.im = mag * sin(ang);
            }
            S[(size_t)n*16 + k] = z;            /* row-major: S11,S12,...     */
        }
    }
    free(dat);

    td->freq = freq; td->S = S; td->Np = Np; td->Rref = Rref;
    return Np;
}

/* index into row-major 4x4 S using 1-based port numbers */
static inline cplx Sij(const cplx *row16, int i, int j)
{
    return row16[(size_t)(i-1)*4 + (j-1)];
}

/* ───────────────────────── main conversion ─────────────────────────────── */
void s4p_default_opts(S4pOpts *o)
{
    o->in_p = 1; o->in_n = 3; o->out_p = 2; o->out_n = 4;
    o->nfft = 16384;
    o->taper_frac = 0.10;
    o->trim_thresh = 1e-3;
    o->verbose = 1;
}

/* Verify the assumed through topology by checking that the chosen forward
 * paths dominate at the lowest measured frequency. Warn (don't fail) if not. */
static void sanity_check_ports(const TouchstoneData *td, const S4pOpts *o)
{
    const cplx *r0 = td->S;  /* first frequency point */
    double thru = 0.5 * (c_abs(Sij(r0, o->out_p, o->in_p)) +
                         c_abs(Sij(r0, o->out_n, o->in_n)));
    /* largest off-diagonal magnitude overall at f_min */
    double best = 0.0; int bi=1, bj=1;
    for (int i = 1; i <= 4; i++)
        for (int j = 1; j <= 4; j++)
            if (i != j) {
                double m = c_abs(Sij(r0, i, j));
                if (m > best) { best = m; bi = i; bj = j; }
            }
    if (thru < 0.5 * best) {
        fprintf(stderr,
            "WARNING: assumed through paths (S%d%d / S%d%d) are weak at f_min "
            "(|thru|=%.3f). Largest term is S%d%d=%.3f.\n"
            "         Check in_p/in_n/out_p/out_n in S4pOpts.\n",
            o->out_p, o->in_p, o->out_n, o->in_n, thru, bi, bj, best);
    }
}

int channel_from_s4p(const char *filename, double *h_fir, int max_taps,
                     double Fs, const S4pOpts *opts)
{
    S4pOpts def; s4p_default_opts(&def);
    const S4pOpts *o = opts ? opts : &def;
    int    NFFT  = (o->nfft > 0) ? o->nfft : 16384;
    double tfrac = (o->taper_frac > 0.0) ? o->taper_frac : 0.10;
    double tthr  = (o->trim_thresh > 0.0) ? o->trim_thresh : 1e-3;

    /* require power of two */
    if (NFFT & (NFFT - 1)) {
        fprintf(stderr, "ERROR: nfft (%d) must be a power of two\n", NFFT);
        return -1;
    }

    TouchstoneData td = {0};
    if (parse_touchstone(filename, &td) < 0) return -1;
    sanity_check_ports(&td, o);

    int Np = td.Np;
    double fmin = td.freq[0], fmax = td.freq[Np - 1];

    /* differential through response Sdd21 on the measured grid */
    cplx *sdd = (cplx *)malloc((size_t)Np * sizeof(cplx));
    if (!sdd) { free_touchstone(&td); return -1; }
    for (int n = 0; n < Np; n++) {
        const cplx *r = td.S + (size_t)n * 16;
        cplx v = c_sub(Sij(r, o->out_p, o->in_p), Sij(r, o->out_p, o->in_n));
        v      = c_sub(v, Sij(r, o->out_n, o->in_p));
        v      = c_add(v, Sij(r, o->out_n, o->in_n));
        sdd[n] = c_scale(v, 0.5);
    }

    /* build conjugate-symmetric spectrum on uniform grid df = Fs/NFFT */
    double df = Fs / NFFT;
    double *Hr = (double *)calloc(NFFT, sizeof(double));
    double *Hi = (double *)calloc(NFFT, sizeof(double));
    if (!Hr || !Hi) { free(sdd); free(Hr); free(Hi); free_touchstone(&td); return -1; }

    double f_taper0 = fmax * (1.0 - tfrac);   /* start of band-edge taper */
    double Hdc = c_abs(sdd[0]);               /* DC magnitude (real)      */
    int    half = NFFT / 2;
    int    seg  = 0;                          /* running interp segment   */

    for (int k = 0; k <= half; k++) {
        double f = k * df;
        cplx H;
        if (f <= 0.0) {
            H = (cplx){ Hdc, 0.0 };
        } else if (f < fmin) {
            /* linear ramp from (DC: Hdc real) to first measured point */
            double t = f / fmin;
            H.re = (1.0 - t) * Hdc + t * sdd[0].re;
            H.im =                   t * sdd[0].im;
        } else if (f >= fmax) {
            H = (cplx){ 0.0, 0.0 };           /* no data above fmax        */
        } else {
            while (seg < Np - 2 && td.freq[seg + 1] < f) seg++;
            double f0 = td.freq[seg], f1 = td.freq[seg + 1];
            double t  = (f - f0) / (f1 - f0);
            H.re = (1.0 - t) * sdd[seg].re + t * sdd[seg + 1].re;
            H.im = (1.0 - t) * sdd[seg].im + t * sdd[seg + 1].im;
            if (f > f_taper0) {               /* raised-cosine to 0 at fmax */
                double w = 0.5 * (1.0 + cos(M_PI * (f - f_taper0) /
                                            (fmax - f_taper0)));
                H = c_scale(H, w);
            }
        }
        Hr[k] = H.re; Hi[k] = H.im;
    }
    Hi[0] = 0.0;                              /* DC strictly real          */
    Hi[half] = 0.0;                           /* Nyquist strictly real     */
    for (int k = 1; k < half; k++) {          /* conjugate symmetry        */
        Hr[NFFT - k] =  Hr[k];
        Hi[NFFT - k] = -Hi[k];
    }
    free(sdd);

    /* IFFT -> real impulse response at Fs */
    fft_radix2(Hr, Hi, NFFT, +1);
    double *h = Hr;                           /* real part is the response */

    /* locate peak */
    int    pk = 0; double pkv = 0.0;
    for (int n = 0; n < NFFT; n++) {
        double a = fabs(h[n]);
        if (a > pkv) { pkv = a; pk = n; }
    }
    double thr = tthr * pkv;
    const int GAP = 8;                        /* tolerate short sub-thr runs */

    /* expand left from peak */
    int start = pk, run = 0;
    for (int n = pk; n >= 0; n--) {
        if (fabs(h[n]) >= thr) { start = n; run = 0; }
        else if (++run > GAP) break;
    }
    /* expand right from peak */
    int end = pk; run = 0;
    for (int n = pk; n < NFFT; n++) {
        if (fabs(h[n]) >= thr) { end = n; run = 0; }
        else if (++run > GAP) break;
    }

    int L = end - start + 1;
    if (L > max_taps) L = max_taps;
    for (int i = 0; i < L; i++) h_fir[i] = h[start + i];

    /* diagnostics */
    if (o->verbose) {
        double dcgain = 0.0;
        for (int i = 0; i < L; i++) dcgain += h_fir[i];
        printf("s4p channel: '%s'\n", filename);
        printf("  points=%d  f=[%.3g, %.3g] Hz  Rref=%.0f\n",
               Np, fmin, fmax, td.Rref);
        printf("  Fs=%.4g Hz  dt=%.4g ps  NFFT=%d  df=%.4g MHz\n",
               Fs, 1e12 / Fs, NFFT, df / 1e6);
        printf("  Sdd21(DC)~%.4f  peak tap=%.4f @ n=%d (%.3f ns)\n",
               Hdc, h[pk], pk, pk * 1e9 / Fs);
        printf("  taps kept=%d  (sum=%.4f, ~%.2f dB DC loss)\n",
               L, dcgain, 20.0 * log10(fabs(dcgain) + 1e-30));
        printf("  trimmed bulk delay = %d samples (%.3f ns)\n",
               start, start * 1e9 / Fs);
    }

    free(Hr); free(Hi);
    free_touchstone(&td);
    return L;
}
