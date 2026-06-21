/* serdes_synth.h — host-only synthetic link model (symbol-spaced, OSF=1).
 * Physical model only: symbols -> channel FIR -> AWGN -> ADC samples.
 * PAM-M with amplitudes normalised to +/-1 full scale. */
#ifndef SERDES_SYNTH_H
#define SERDES_SYNTH_H

#define SYN_M        4         /* PAM levels (PAM-4) */
#define SYN_TGT_LEN  3         /* partial-response target length (memory L-1) */
#define SYN_TGT_SUM  0.80      /* |target| L1 norm: keep equalised eye < full scale */

/* Normalised amplitude for index 0..M-1: PAM-4 -> {-1,-1/3,+1/3,+1}. */
double syn_amp(int idx);
/* Gray code for symbol index (for BER counting). */
int    syn_gray(int idx);
/* Slice a normalised sample to the nearest PAM-M index. */
int    syn_slice(double y);

int    syn_channel(double *h, int maxlen);          /* fills h[], returns len */
void   syn_pr_target(double *ht, int *L);           /* normalised PR target   */
void   syn_generate(const int *sym, int nsym, const double *h, int hlen,
                    double noise_sigma, int adc_bits, double *rx);

#endif /* SERDES_SYNTH_H */