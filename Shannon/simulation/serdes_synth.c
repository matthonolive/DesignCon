#include "serdes_synth.h"
#include <math.h>
#include <stdlib.h>

double syn_amp(int idx){ return (2.0*idx - (SYN_M-1)) / (double)(SYN_M-1); }

int syn_gray(int idx){ return idx ^ (idx >> 1); }

int syn_slice(double y){
    int best=0; double bd=1e300;
    for (int j=0;j<SYN_M;j++){ double d=y-syn_amp(j); d*=d; if(d<bd){bd=d;best=j;} }
    return best;
}

int syn_channel(double *h, int maxlen){
    /* lossy channel, one precursor + postcursor tail; normalise to unit DC */
    static const double c[] = { 0.12, 1.00, 0.55, 0.28, 0.12, 0.05 };
    int n = (int)(sizeof(c)/sizeof(c[0])); if (n>maxlen) n=maxlen;
    double s=0; for(int i=0;i<n;i++) s+=c[i];
    for (int i=0;i<n;i++) h[i]=c[i]/s;
    return n;
}

void syn_pr_target(double *ht, int *L){
    double raw[SYN_TGT_LEN] = { 1.00, 0.50, 0.20 };
    double s=0; for(int i=0;i<SYN_TGT_LEN;i++) s+=fabs(raw[i]);
    for (int i=0;i<SYN_TGT_LEN;i++) ht[i] = raw[i] * (SYN_TGT_SUM / s);
    *L = SYN_TGT_LEN;
}

static double randn(void){
    double u1=(rand()+1.0)/(RAND_MAX+2.0), u2=(rand()+1.0)/(RAND_MAX+2.0);
    return sqrt(-2.0*log(u1))*cos(2.0*M_PI*u2);
}
static double adc_q(double x,int B){
    if (x> 0.999969) x= 0.999969; if (x<-1.0) x=-1.0;
    int lv=1<<B; int xq=(int)lround((x+1.0)*(lv-1)/2.0);
    return xq*2.0/(lv-1)-1.0;
}
void syn_generate(const int *sym, int nsym, const double *h, int hlen,
                  double noise_sigma, int adc_bits, double *rx){
    for (int n=0;n<nsym;n++){
        double s=0;
        for (int k=0;k<hlen;k++){ int m=n-k; if(m>=0) s+=h[k]*syn_amp(sym[m]); }
        s += noise_sigma*randn();
        rx[n]=adc_q(s,adc_bits);
    }
}