/* bridge_cosim.c — fixed-point sign-off for the PAM-M receiver datapath.
 * Integer FFE+MLSE vs a floating-point reference on the SAME generated data:
 * decision agreement, MLSE vs slicer SER, FFE-output quant SNR.
 *   gcc -O2 -Itasks -o cosim bridge/bridge_cosim.c -lm && ./cosim */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include "serdes_workload.h"
#include "eq_fixed.h"
#include "mlse_fixed.h"

static double lvl_f(int j){ return mlse_levels[j]/(double)(1<<MLSE_LVL_Q); }
static double tgt_f(int k){ return mlse_target[k]/(double)(1<<EQ_SAMP_Q); }
static int slice_f(double y){ int b=0; double bd=1e300;
    for(int j=0;j<MLSE_M;j++){double d=y-lvl_f(j);d*=d;if(d<bd){bd=d;b=j;}} return b; }

/* float reference MLSE (mirrors mlse_fixed.h; double metrics, exp table) */
typedef struct { double pm[256], et[256*8]; int8_t sym[MLSE_TB*MLSE_NS], pre[MLSE_TB*MLSE_NS];
                 int pw[8]; int head,filled; } mlse_f;
static void mlf_init(mlse_f*s){
    s->pw[0]=1; for(int d=1;d<MLSE_L;d++)s->pw[d]=s->pw[d-1]*MLSE_M;
    for(int v=0;v<MLSE_NS;v++)for(int j=0;j<MLSE_M;j++){
        double e=tgt_f(0)*lvl_f(j);
        for(int d=1;d<MLSE_L;d++){int idx=(v/s->pw[d-1])%MLSE_M; e+=tgt_f(d)*lvl_f(idx);}
        s->et[v*MLSE_M+j]=e; }
    for(int v=0;v<MLSE_NS;v++)s->pm[v]=0; s->head=0;s->filled=0;
}
static int mlf_step(mlse_f*s,double y){ int ns=MLSE_NS,m=MLSE_M;
    int8_t*bs=s->sym+s->head*ns,*bp=s->pre+s->head*ns; double nm[256];
    for(int v=0;v<ns;v++)nm[v]=1e300;
    for(int v=0;v<ns;v++){double b=s->pm[v]; if(b>=1e299)continue;
        for(int j=0;j<m;j++){double d=y-s->et[v*m+j],c=b+d*d; int vp=(m*v+j)%ns;
            if(c<nm[vp]){nm[vp]=c;bs[vp]=j;bp[vp]=v;}}}
    double mn=nm[0];int best=0; for(int v=1;v<ns;v++)if(nm[v]<mn){mn=nm[v];best=v;}
    for(int v=0;v<ns;v++)s->pm[v]=nm[v]-mn;
    s->head=(s->head+1)%MLSE_TB; if(s->filled<MLSE_TB)s->filled++;
    if(s->filled<MLSE_TB)return -1;
    int st=best,idx=(s->head+MLSE_TB-1)%MLSE_TB;
    for(int t=0;t<MLSE_TB-1;t++){st=s->pre[idx*ns+st];idx=(idx+MLSE_TB-1)%MLSE_TB;}
    return s->sym[idx*ns+st]; }

int main(void){
    eq_state_t fs; eq_reset(&fs);
    int64_t pm[256]; int8_t tbs[MLSE_TB*256],tbp[MLSE_TB*256]; int32_t et[256*8];
    mlse_t ml; mlse_init(&ml,pm,tbs,tbp,et,mlse_target,mlse_levels,
                         MLSE_NS,MLSE_M,MLSE_L,MLSE_TB,MLSE_LVL_Q);
    double rb[EQ_FFE_LEN]={0}, tf[EQ_FFE_LEN];
    for(int k=0;k<EQ_FFE_LEN;k++) tf[k]=eq_ffe_taps[k]/(double)(1<<EQ_TAP_Q);
    mlse_f mf; mlf_init(&mf);

    long mism=0, me=0,mt=0, mfe=0,mft=0, se=0,st=0;
    double sig=0,err=0;
    int df=0,dx=0;
    for(int n=0;n<EQ_NSTIM;n++){
        int16_t yq=eq_ffe_step(&fs,eq_ffe_taps,eq_stimulus[n]);
        for(int i=EQ_FFE_LEN-1;i>0;i--) rb[i]=rb[i-1];
        rb[0]=eq_stimulus[n]/(double)(1<<EQ_SAMP_Q);
        double yf=0; for(int k=0;k<EQ_FFE_LEN;k++) yf+=tf[k]*rb[k];
        double yqf=yq/(double)(1<<EQ_SAMP_Q);
        if(fabs(yf)<0.999){double d=yf-yqf; sig+=yf*yf; err+=d*d;}
        int sidx=slice_f(yqf);
        if(n>=MLSE_SKIP){ st++; se+=(sidx!=eq_expect[n]); }
        int rf=mlf_step(&mf,yf), rx=mlse_step(&ml,yq);
        if(rf>=0&&rx>=0&&rf!=rx) mism++;
        if(rx>=0){ if(dx>=MLSE_SKIP){mt++; me+=(rx!=eq_expect[dx]);} dx++; }
        if(rf>=0){ if(df>=MLSE_SKIP){mft++; mfe+=(rf!=eq_expect[df]);} df++; }
    }
    double snr=10.0*log10(sig/(err+1e-30));
    printf("cosim PAM-%d (FFE -> PR-target -> MLSE, %d-state):\n",MLSE_M,MLSE_NS);
    printf("  FFE-output quant SNR (in-range) = %.1f dB\n",snr);
    printf("  MLSE  symbol errors: integer = %ld/%ld (%.2e)   float = %ld/%ld (%.2e)\n",
           me,mt,mt?(double)me/mt:0, mfe,mft,mft?(double)mfe/mft:0);
    printf("  slicer symbol errors = %ld / %ld  (%.2e)\n",se,st,st?(double)se/st:0);
    printf("  (integer vs float path differences = %ld; benign if SERs match)\n",mism);
    printf("  -> fixed-point port %s the float design (SER); MLSE %s the slicer.\n",
           (me<=mfe)?"MATCHES":"DEGRADES vs",
           (mt&&st&&(double)me/mt<(double)se/st)?"beats":"does NOT beat");
    return 0;
}