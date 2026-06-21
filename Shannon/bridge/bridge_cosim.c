/* bridge_cosim.c — fixed-point sign-off for the receiver datapath.
 * Runs the firmware's integer FFE+MLSE against a floating-point reference
 * using the SAME generated data, and reports:
 *   - decision agreement (integer detector vs float detector)
 *   - MLSE symbol errors vs ground truth
 *   - memoryless-slicer symbol errors vs ground truth (MLSE should beat it)
 *   - FFE-output quantisation SNR
 *
 *   gcc -O2 -Itasks -o cosim bridge/bridge_cosim.c -lm && ./cosim
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>

#include "serdes_workload.h"
#include "eq_fixed.h"
#include "mlse_fixed.h"

/* ── float reference MLSE (mirrors mlse_fixed.h, double metrics) ──────── */
typedef struct { double pm[256]; int8_t sym[MLSE_TB*MLSE_NS], pre[MLSE_TB*MLSE_NS];
                 int pw[8]; int head, filled; } mlse_f;
static void mlf_init(mlse_f*s){ for(int v=0;v<MLSE_NS;v++)s->pm[v]=0; s->head=0;s->filled=0;
    s->pw[0]=1; for(int d=1;d<MLSE_L;d++)s->pw[d]=s->pw[d-1]*MLSE_M; }
static double mlf_exp(int v,int j){ double e=(mlse_target[0]/32768.0)*mlse_levels[j];
    for(int d=1;d<MLSE_L;d++){int idx=(v/((int[]){1,MLSE_M,MLSE_M*MLSE_M,MLSE_M*MLSE_M*MLSE_M})[d-1])%MLSE_M;
        e+=(mlse_target[d]/32768.0)*mlse_levels[idx];} return e; }
static int mlf_step(mlse_f*s,double y){ int ns=MLSE_NS,m=MLSE_M;
    int8_t*bs=s->sym+s->head*ns,*bp=s->pre+s->head*ns; double nm[256];
    for(int v=0;v<ns;v++)nm[v]=1e300;
    for(int v=0;v<ns;v++){ double b=s->pm[v]; if(b>=1e299)continue;
        for(int j=0;j<m;j++){ double e=mlf_exp(v,j),d=y-e,c=b+d*d; int vp=(j+m*v)%ns;
            if(c<nm[vp]){nm[vp]=c;bs[vp]=j;bp[vp]=v;} } }
    double mn=nm[0];int best=0; for(int v=1;v<ns;v++)if(nm[v]<mn){mn=nm[v];best=v;}
    for(int v=0;v<ns;v++)s->pm[v]=nm[v]-mn;
    s->head=(s->head+1)%MLSE_TB; if(s->filled<MLSE_TB)s->filled++;
    if(s->filled<MLSE_TB)return -1;
    int st=best,idx=(s->head+MLSE_TB-1)%MLSE_TB;
    for(int t=0;t<MLSE_TB-1;t++){st=s->pre[idx*ns+st];idx=(idx+MLSE_TB-1)%MLSE_TB;}
    return s->sym[idx*ns+st]; }

int main(void){
    /* fixed datapath */
    eq_state_t fs; eq_reset(&fs);
    int64_t pm[256]; int8_t tbs[MLSE_TB*256], tbp[MLSE_TB*256];
    mlse_t ml; mlse_init(&ml,pm,tbs,tbp,mlse_target,mlse_levels,MLSE_NS,MLSE_M,MLSE_L,MLSE_TB);
    /* float datapath */
    double rb[EQ_FFE_LEN]={0}, tf[EQ_FFE_LEN];
    for(int k=0;k<EQ_FFE_LEN;k++) tf[k]=eq_ffe_taps[k]/(double)(1<<EQ_TAP_Q);
    mlse_f mf; mlf_init(&mf);

    long mism=0;                     /* fixed vs float detector decisions     */
    long me=0,mt=0, se=0,st=0;       /* MLSE / slicer errors vs ground truth  */
    long mfe=0,mft=0;                /* float MLSE errors vs ground truth     */
    double sig=0,err=0; long inr=0;
    int df=0, dx=0;                  /* emission counters (float / fixed)     */

    for(int n=0;n<EQ_NSTIM;n++){
        /* fixed FFE */
        int16_t yq = eq_ffe_step(&fs, eq_ffe_taps, eq_stimulus[n]);
        /* float FFE */
        for(int i=EQ_FFE_LEN-1;i>0;i--) rb[i]=rb[i-1];
        rb[0]=eq_stimulus[n]/(double)(1<<EQ_SAMP_Q);
        double yf=0; for(int k=0;k<EQ_FFE_LEN;k++) yf+=tf[k]*rb[k];

        /* quant SNR of FFE output (in-range) */
        double yqf=yq/(double)(1<<EQ_SAMP_Q);
        if(fabs(yf)<0.999){ double d=yf-yqf; sig+=yf*yf; err+=d*d; inr++; }

        /* slicer on fixed FFE out */
        int sidx = (yq<0)?0:1;
        if(n>=MLSE_SKIP){ int t=n-7; if(t>=0){ st++; se += (sidx!=eq_expect[n]); } }

        /* detectors */
        int rf = mlf_step(&mf, yf);
        int rx = mlse_step(&ml, yq);
        if(rf>=0 && rx>=0){ if(rf!=rx) mism++; }
        if(rx>=0){ if(dx>=MLSE_SKIP){ mt++; me += (rx!=eq_expect[dx]); } dx++; }
        if(rf>=0){ if(df>=MLSE_SKIP){ mft++; mfe += (rf!=eq_expect[df]); } df++; }
    }
    double snr = 10.0*log10(sig/(err+1e-30));
    printf("cosim (FFE -> PR-target -> MLSE):\n");
    printf("  FFE-output quant SNR (in-range) = %.1f dB\n", snr);
    printf("  MLSE  symbol errors: integer = %ld/%ld (%.2e)   float = %ld/%ld (%.2e)\n",
           me,mt, mt?(double)me/mt:0, mfe,mft, mft?(double)mfe/mft:0);
    printf("  slicer symbol errors = %ld / %ld  (%.2e)\n", se, st, st?(double)se/st:0);
    printf("  (integer vs float path differences = %ld; benign if SERs match)\n", mism);
    int ser_ok  = (me <= mfe);                          /* integer no worse than float */
    int beats   = (mt && st && (double)me/mt < (double)se/st);
    printf("  -> fixed-point port %s the float design (SER); MLSE %s the slicer.\n",
           ser_ok ? "MATCHES" : "DEGRADES vs",
           beats  ? "beats"   : "does NOT beat");
    return 0;
}
