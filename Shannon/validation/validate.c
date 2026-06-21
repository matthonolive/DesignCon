/* validate.c — BER/SER vs SNR sweep. Trains+freezes the FFE at a design SNR
 * (taps are baked in the real system), then sweeps noise and compares, on the
 * SAME equalised signal:
 *     MLSE (float)   MLSE (fixed, the firmware kernel)
 *     DFE  (float, target-cancellation)   memoryless slicer (float)
 * over multiple seeds. Emits CSV to stdout-redirect via the Makefile.
 *
 *   gcc -O2 -Isimulation -Itasks -o validate bridge/validate.c \
 *       simulation/serdes_synth.c -lm
 *   ./validate            # human table
 *   ./validate csv > ber_sweep.csv
 */
#include "serdes_synth.h"
#include "mlse_fixed.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define FFE_LEN 14
#define TAP_Q   14
#define SAMP_Q  15
#define LVL_Q   15
#define TB      16
#define FFE_DELAY 6
#define ADC_BITS  10
#define NSYM    6000
#define SEEDS   12

static int M,L,ns;
static double h[32]; static int hlen;
static double htf[8];
static double wf[FFE_LEN];           /* trained float FFE (frozen)      */
static int16_t tapq[FFE_LEN], tgtq[8], lvlq[8];

static int q_round(double x,int qb){ long v=lrint(x*(double)(1<<qb));
    return v>32767?32767:(v<-32768?-32768:(int)v); }
/* slice against the PR main-cursor-scaled constellation (target[0]*amp) */
static int slice_sc(double y){ int b=0; double bd=1e300;
    for(int j=0;j<M;j++){ double d=y-htf[0]*syn_amp(j); d*=d; if(d<bd){bd=d;b=j;} } return b; }
static int16_t fsat16(int64_t v){ return v>32767?32767:(v<-32768?-32768:(int16_t)v); }
static int bits_per_sym(void){ int b=0,m=M; while(m>>=1)b++; return b; }
static int biterr(int a,int b){ int x=syn_gray(a)^syn_gray(b),c=0; while(x){c+=x&1;x>>=1;} return c; }

/* train float FFE to PR target at a clean design SNR, then freeze */
static void train(double design_sigma){
    hlen=syn_channel(h,32); syn_pr_target(htf,&L);
    M=SYN_M; ns=1; for(int i=0;i<L-1;i++) ns*=M;
    int N=NSYM; int *sym=malloc(N*sizeof(int)); double *rx=malloc(N*sizeof(double));
    srand(777); for(int i=0;i<N;i++) sym[i]=rand()%M;
    syn_generate(sym,N,h,hlen,design_sigma,ADC_BITS,rx);
    memset(wf,0,sizeof(wf)); wf[FFE_DELAY]=1.0;
    double rb[FFE_LEN]={0}, mu=0.01; int DELAY=FFE_DELAY+1;
    for(int pass=0;pass<12;pass++){ memset(rb,0,sizeof(rb));
        for(int n=0;n<N;n++){
            for(int i=FFE_LEN-1;i>0;i--) rb[i]=rb[i-1]; rb[0]=rx[n];
            double y=0; for(int k=0;k<FFE_LEN;k++) y+=wf[k]*rb[k];
            double dref=0; for(int k=0;k<L;k++){int m=n-DELAY-k; if(m>=0) dref+=htf[k]*syn_amp(sym[m]);}
            double e=dref-y; for(int k=0;k<FFE_LEN;k++) wf[k]+=mu*e*rb[k];
        } }
    for(int k=0;k<FFE_LEN;k++) tapq[k]=q_round(wf[k],TAP_Q);
    for(int k=0;k<L;k++) tgtq[k]=q_round(htf[k],SAMP_Q);
    for(int j=0;j<M;j++) lvlq[j]=q_round(syn_amp(j),LVL_Q);
    free(sym); free(rx);
}

/* float MLSE (table) */
typedef struct{ double pm[256],et[256*8]; int8_t sy[TB*256],pr[TB*256]; int pw[8],head,fill; } mlf_t;
static void mlf_init(mlf_t*s){ s->pw[0]=1; for(int d=1;d<L;d++)s->pw[d]=s->pw[d-1]*M;
    for(int v=0;v<ns;v++)for(int j=0;j<M;j++){ double e=htf[0]*syn_amp(j);
        for(int d=1;d<L;d++){int idx=(v/s->pw[d-1])%M; e+=htf[d]*syn_amp(idx);} s->et[v*M+j]=e; }
    for(int v=0;v<ns;v++)s->pm[v]=0; s->head=0;s->fill=0; }
static int mlf_step(mlf_t*s,double y){ int8_t*bs=s->sy+s->head*ns,*bp=s->pr+s->head*ns;
    double nm[256]; for(int v=0;v<ns;v++)nm[v]=1e300;
    for(int v=0;v<ns;v++){double b=s->pm[v]; if(b>=1e299)continue;
        for(int j=0;j<M;j++){double d=y-s->et[v*M+j],c=b+d*d;int vp=(M*v+j)%ns;
            if(c<nm[vp]){nm[vp]=c;bs[vp]=j;bp[vp]=v;}}}
    double mn=nm[0];int best=0;for(int v=1;v<ns;v++)if(nm[v]<mn){mn=nm[v];best=v;}
    for(int v=0;v<ns;v++)s->pm[v]=nm[v]-mn;
    s->head=(s->head+1)%TB; if(s->fill<TB)s->fill++; if(s->fill<TB)return -1;
    int st=best,idx=(s->head+TB-1)%TB; for(int t=0;t<TB-1;t++){st=s->pr[idx*ns+st];idx=(idx+TB-1)%TB;}
    return s->sy[idx*ns+st]; }

/* one run at sigma/seed -> fills error counters (SER + bit errors) */
typedef struct{ long se,st, be,bt; } acc_t;
static int g_mlse, g_samp;             /* calibrated alignment lags */

static void run(double sigma,int seed,
                acc_t *mf,acc_t *mx,acc_t *dfe,acc_t *slc,int calib){
    int N=NSYM; int *sym=malloc(N*sizeof(int)); double *rx=malloc(N*sizeof(double));
    srand(seed); for(int i=0;i<N;i++) sym[i]=rand()%M;
    syn_generate(sym,N,h,hlen,sigma,ADC_BITS,rx);

    /* fixed kernel state */
    int64_t pm[256]; int8_t tbs[TB*256],tbp[TB*256]; int32_t et[256*8];
    mlse_t ml; mlse_init(&ml,pm,tbs,tbp,et,tgtq,lvlq,ns,M,L,TB,LVL_Q);
    int16_t fbuf[FFE_LEN]={0};
    mlf_t mlf; mlf_init(&mlf);
    double rbf[FFE_LEN]={0};
    int dhist[8]={0};                  /* DFE decision history (indices) */

    int *decF=malloc(N*sizeof(int)), nF=0;     /* float mlse emissions */
    int *decX=malloc(N*sizeof(int)), nX=0;     /* fixed mlse emissions */
    int *decD=malloc(N*sizeof(int));           /* dfe per sample */
    int *decS=malloc(N*sizeof(int));           /* slicer per sample */

    for(int n=0;n<N;n++){
        int16_t xq=q_round(rx[n],SAMP_Q);
        /* fixed FFE */
        int16_t yq=fsat16(({int64_t a=0; for(int i=FFE_LEN-1;i>0;i--)fbuf[i]=fbuf[i-1]; fbuf[0]=xq;
                            for(int k=0;k<FFE_LEN;k++)a+=(int32_t)tapq[k]*(int32_t)fbuf[k]; a>>TAP_Q;}));
        /* float FFE */
        for(int i=FFE_LEN-1;i>0;i--) rbf[i]=rbf[i-1]; rbf[0]=rx[n];
        double yf=0; for(int k=0;k<FFE_LEN;k++) yf+=wf[k]*rbf[k];
        /* detectors on the same equalised sample */
        decS[n]=slice_sc(yf);
        double yd=yf; for(int k=1;k<L;k++) yd-=htf[k]*syn_amp(dhist[k-1]);
        int dd=slice_sc(yd); decD[n]=dd;
        for(int k=L-2;k>0;k--) dhist[k]=dhist[k-1]; if(L>=2) dhist[0]=dd;
        int rf=mlf_step(&mlf,yf); if(rf>=0) decF[nF++]=rf;
        int rx2=mlse_step(&ml,yq);  if(rx2>=0) decX[nX++]=rx2;
    }

    if(calib){ /* find lags minimising errors vs ground truth */
        long bb=1<<30; for(int g=-(FFE_LEN+L+8);g<=4;g++){ long e=0,c=0;
            for(int d=0;d<nF;d++){int t=d+g; if(t<0||t>=N)continue; e+=(decF[d]!=sym[t]);c++;}
            if(c>nF/2&&e<bb){bb=e;g_mlse=g;} }
        bb=1<<30; for(int g=-(FFE_LEN+L+8);g<=4;g++){ long e=0,c=0;
            for(int n=0;n<N;n++){int t=n+g; if(t<0||t>=N)continue; e+=(decS[n]!=sym[t]);c++;}
            if(c>N/2&&e<bb){bb=e;g_samp=g;} }
    }
    int SK=(g_mlse<0?-g_mlse:0)+TB+4;
    for(int d=SK;d<nF;d++){int t=d+g_mlse; if(t<0||t>=N)continue; mf->st++; if(decF[d]!=sym[t]){mf->se++; mf->be+=biterr(decF[d],sym[t]);} mf->bt+=bits_per_sym();}
    for(int d=SK;d<nX;d++){int t=d+g_mlse; if(t<0||t>=N)continue; mx->st++; if(decX[d]!=sym[t]){mx->se++; mx->be+=biterr(decX[d],sym[t]);} mx->bt+=bits_per_sym();}
    int SS=(g_samp<0?-g_samp:0)+4;
    for(int n=SS;n<N;n++){int t=n+g_samp; if(t<0||t>=N)continue;
        dfe->st++; if(decD[n]!=sym[t]){dfe->se++; dfe->be+=biterr(decD[n],sym[t]);} dfe->bt+=bits_per_sym();
        slc->st++; if(decS[n]!=sym[t]){slc->se++; slc->be+=biterr(decS[n],sym[t]);} slc->bt+=bits_per_sym();}

    free(sym);free(rx);free(decF);free(decX);free(decD);free(decS);
}

/* measure rx signal power for an SNR estimate */
static double sig_power(void){
    int N=2000; int *sym=malloc(N*sizeof(int)); double *rx=malloc(N*sizeof(double));
    srand(1); for(int i=0;i<N;i++) sym[i]=rand()%M;
    syn_generate(sym,N,h,hlen,0.0,ADC_BITS,rx);
    double p=0; for(int n=64;n<N;n++) p+=rx[n]*rx[n]; p/=(N-64);
    free(sym);free(rx); return p;
}

int main(int argc,char**argv){
    int csv = (argc>=2 && !strcmp(argv[1],"csv"));
    train(0.008);                       /* design-SNR taps, frozen */
    double P=sig_power();
    double sig0=0.16; { acc_t a={0},b={0},c={0},d={0}; run(sig0,101,&a,&b,&c,&d,1); } /* calibrate lags */

    double sigmas[]={0.02,0.04,0.06,0.08,0.10,0.12,0.14,0.16,0.18,0.20};
    int NS=sizeof(sigmas)/sizeof(sigmas[0]);
    if(csv) printf("sigma,snr_db,ser_mlse_float,ser_mlse_fixed,ser_dfe,ser_slicer,ber_mlse_fixed,ber_dfe\n");
    else printf("PAM-%d, %d-state MLSE, %d seeds x %d sym/pt; FFE frozen at design SNR\n"
                "  sigma  SNR_dB    MLSEf     MLSEx      DFE      slicer | BER(MLSEx) BER(DFE)\n",
                M,ns,SEEDS,NSYM);
    for(int i=0;i<NS;i++){
        acc_t mf={0},mx={0},df={0},sl={0};
        for(int s=0;s<SEEDS;s++) run(sigmas[i],1000+s*7+i,&mf,&mx,&df,&sl,0);
        double snr=10.0*log10(P/(sigmas[i]*sigmas[i]));
        double serMf=(double)mf.se/mf.st, serMx=(double)mx.se/mx.st;
        double serD=(double)df.se/df.st, serS=(double)sl.se/sl.st;
        double berMx=(double)mx.be/mx.bt, berD=(double)df.be/df.bt;
        if(csv) printf("%.3f,%.2f,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e\n",
                       sigmas[i],snr,serMf,serMx,serD,serS,berMx,berD);
        else printf("  %.3f  %6.2f  %.2e  %.2e  %.2e  %.2e | %.2e  %.2e\n",
                    sigmas[i],snr,serMf,serMx,serD,serS,berMx,berD);
    }
    return 0;
}