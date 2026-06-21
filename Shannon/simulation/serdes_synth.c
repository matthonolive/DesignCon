#include "serdes_synth.h"
#include <math.h>
#include <stdlib.h>

double syn_amp(int idx){ return 2.0 * idx - 1.0; }   /* {-1,+1} */

int syn_channel(double *h, int maxlen){
    /* lossy channel with one precursor and a postcursor tail; main at idx 1 */
    static const double c[] = { 0.12, 1.00, 0.55, 0.28, 0.12, 0.05 };
    int n = (int)(sizeof(c)/sizeof(c[0]));
    if (n > maxlen) n = maxlen;
    for (int i = 0; i < n; i++) h[i] = c[i];
    return n;
}

void syn_pr_target(double *ht, int *L){
    /* deliberate controlled ISI so MLSE beats a memoryless slicer */
    ht[0] = 1.00; ht[1] = 0.50; ht[2] = 0.20;
    *L = SYN_TGT_LEN;
}

static double randn(void){
    double u1 = (rand() + 1.0) / (RAND_MAX + 2.0);
    double u2 = (rand() + 1.0) / (RAND_MAX + 2.0);
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

static double adc_q(double x, int B){
    if (x >  0.999969) x =  0.999969;
    if (x < -1.0)      x = -1.0;
    int lv = 1 << B;
    int xq = (int)lround((x + 1.0) * (lv - 1) / 2.0);
    return xq * 2.0 / (lv - 1) - 1.0;
}

void syn_generate(const int *sym, int nsym, const double *h, int hlen,
                  double noise_sigma, int adc_bits, double *rx){
    for (int n = 0; n < nsym; n++){
        double s = 0.0;
        for (int k = 0; k < hlen; k++){
            int m = n - k;
            if (m >= 0) s += h[k] * syn_amp(sym[m]);
        }
        s += noise_sigma * randn();
        rx[n] = adc_q(s, adc_bits);
    }
}
