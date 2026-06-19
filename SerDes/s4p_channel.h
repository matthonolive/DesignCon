#ifndef S4P_CHANNEL_H
#define S4P_CHANNEL_H

/*
 * s4p_channel.h
 *
 * Self-contained Touchstone (.s4p) -> channel FIR loader.
 *
 * Replaces the "offline MATLAB FitChannel()" step referenced in serdes_sim.c:
 * reads a 4-port Touchstone file, forms the differential through response
 * Sdd21(f), interpolates it onto the simulation sample grid, performs an
 * IFFT to a real impulse response, trims the bulk propagation delay, and
 * returns FIR taps sampled at Fs = OSF * DATA_RATE.
 *
 * No external dependencies (just libm).
 */

/* Options controlling the s4p -> taps conversion. Pass NULL to use defaults. */
typedef struct {
    /* Single-ended port indices (1-based) forming the differential ports.
     * Defaults match a {1<->2, 3<->4} through topology:
     *   diff-IN  = (in_p, in_n) = (1, 3)
     *   diff-OUT = (out_p, out_n) = (2, 4)
     * Sdd21 = 0.5*(S[out_p,in_p] - S[out_p,in_n] - S[out_n,in_p] + S[out_n,in_n])
     */
    int in_p, in_n, out_p, out_n;

    int    nfft;          /* IFFT length (power of 2). 0 -> 16384            */
    double taper_frac;    /* raised-cosine band-edge taper as a fraction of  */
                          /* f_max, e.g. 0.10 tapers the top 10%. 0 -> 0.10  */
    double trim_thresh;   /* tap-trim threshold, fraction of peak. 0 -> 1e-3 */
    int    verbose;       /* 1 -> print a diagnostic summary to stdout       */
} S4pOpts;

/* Fill *opts with defaults. */
void s4p_default_opts(S4pOpts *opts);

/*
 * Read 'filename', synthesize the channel impulse response at sample rate Fs,
 * and write up to max_taps FIR taps into h_fir[].
 *
 * Returns the number of taps written (>0) on success, or -1 on error.
 */
int channel_from_s4p(const char *filename, double *h_fir, int max_taps,
                     double Fs, const S4pOpts *opts);

#endif /* S4P_CHANNEL_H */