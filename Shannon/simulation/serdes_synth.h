/* serdes_synth.h — host-only synthetic link model (symbol-spaced, OSF=1).
 *
 * Physical model only: symbols -> channel FIR -> AWGN -> ADC-quantised
 * samples. No equalisation here (that is the receiver's job, on firmware).
 * A self-contained stand-in for the .s4p path so the whole flow runs with
 * no external file; swap channel_from_s4p() in for hardware-grade channels.
 */
#ifndef SERDES_SYNTH_H
#define SERDES_SYNTH_H

#define SYN_M       2          /* PAM levels (NRZ) */
#define SYN_TGT_LEN 3          /* partial-response target length (memory L-1) */

/* Symbol amplitude for an index 0..SYN_M-1 (NRZ: 0->-1, 1->+1). */
double syn_amp(int idx);

/* Channel impulse response (symbol-spaced). Returns length, fills h[]. */
int    syn_channel(double *h, int maxlen);

/* Recommended PR target the FFE shapes toward. Fills ht[0..*L-1], main first. */
void   syn_pr_target(double *ht, int *L);

/* Generate rx[n] = sum_k h[k]*amp(sym[n-k]) + noise, ADC-quantised to [-1,1).
 * sym[] are symbol indices; rx[] length nsym. */
void   syn_generate(const int *sym, int nsym, const double *h, int hlen,
                    double noise_sigma, int adc_bits, double *rx);

#endif /* SERDES_SYNTH_H */
