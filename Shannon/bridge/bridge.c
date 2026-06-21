/* bridge.c — host mediator. Build via the Makefile 'workload' target.
 *   ./bridge gen [out.h] [n_export]
 *
 * Drives the synthetic link, trains the RX FFE to a partial-response target,
 * quantises, runs the EXACT fixed-point FFE+MLSE the firmware will run to
 * align ground truth, and emits the generated data header. */
#include "serdes_synth.h"
#include "bridge.h"
#include "mlse_fixed.h"          /* shared MLSE kernel (runtime-parametrised) */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define RX_FFE_LEN 14
#define TAP_Q      14
#define SAMP_Q     15
#define TB_DEPTH   16
#define FFE_DELAY  6            /* reference tap inside the FFE for training  */
#define NOISE_SIG  0.010
#define ADC_BITS   10

static int q_round(double x, int qb, int lo, int hi){
    long v = lrint(x * (double)(1 << qb));
    if (v < lo) v = lo; if (v > hi) v = hi; return (int)v;
}
static int16_t fsat16(int64_t v){ return v>32767?32767:(v<-32768?-32768:(int16_t)v); }

/* fixed FFE identical to eq_fixed.h::eq_ffe_step (runtime length) */
static int16_t ffe_fix(int16_t *buf, const int16_t *taps, int16_t x){
    for (int i = RX_FFE_LEN-1; i>0; i--) buf[i]=buf[i-1];
    buf[0]=x;
    int64_t acc=0;
    for (int k=0;k<RX_FFE_LEN;k++) acc += (int32_t)taps[k]*(int32_t)buf[k];
    return fsat16(acc >> TAP_Q);
}

int bridge_gen(const char *out, int n_export){
    /* ── physical model ──────────────────────────────────────────────── */
    double h[32]; int hlen = syn_channel(h, 32);
    double htf[8]; int L;  syn_pr_target(htf, &L);
    int M = SYN_M, ns = 1; for (int i=0;i<L-1;i++) ns*=M;

    int NSYM = n_export + 64;
    int *sym   = malloc(NSYM*sizeof(int));
    double *rx = malloc(NSYM*sizeof(double));
    srand(2024);
    for (int i=0;i<NSYM;i++) sym[i] = rand() & 1;
    syn_generate(sym, NSYM, h, hlen, NOISE_SIG, ADC_BITS, rx);

    /* ── LMS: train FFE so ffe(rx) ~= (PR target) conv (symbols), delayed ─ */
    double w[RX_FFE_LEN]={0}; w[FFE_DELAY]=1.0;
    double rb[RX_FFE_LEN]={0};
    double mu=0.01, mse=0;
    int DELAY = FFE_DELAY + 1;                 /* +channel main index */
    for (int pass=0; pass<8; pass++){
        memset(rb,0,sizeof(rb)); mse=0; int cnt=0;
        for (int n=0;n<NSYM;n++){
            for (int i=RX_FFE_LEN-1;i>0;i--) rb[i]=rb[i-1];
            rb[0]=rx[n];
            double y=0; for(int k=0;k<RX_FFE_LEN;k++) y+=w[k]*rb[k];
            /* desired = target convolved with the symbol stream, delayed */
            double dref=0;
            for (int k=0;k<L;k++){ int m=n-DELAY-k; if(m>=0) dref += htf[k]*syn_amp(sym[m]); }
            double e = dref - y;
            for (int k=0;k<RX_FFE_LEN;k++) w[k]+= mu*e*rb[k];
            if(n>64){ mse+=e*e; cnt++; }
        }
        mse/=cnt;
    }

    /* ── quantise ────────────────────────────────────────────────────── */
    int16_t taps[RX_FFE_LEN];
    for (int k=0;k<RX_FFE_LEN;k++) taps[k]=(int16_t)q_round(w[k],TAP_Q,-32768,32767);
    int16_t tgt[8]; for(int k=0;k<L;k++) tgt[k]=(int16_t)q_round(htf[k],SAMP_Q,-32768,32767);
    int16_t levels[8]; for(int j=0;j<M;j++) levels[j]=(int16_t)lround(syn_amp(j));
    int16_t *stim = malloc(n_export*sizeof(int16_t));
    for (int n=0;n<n_export;n++) stim[n]=(int16_t)q_round(rx[n],SAMP_Q,-32768,32767);

    /* ── run the EXACT fixed kernel to get the decision sequence ──────── */
    int64_t pm[256]; int8_t tbs[TB_DEPTH*256], tbp[TB_DEPTH*256];
    mlse_t ml; mlse_init(&ml, pm, tbs, tbp, tgt, levels, ns, M, L, TB_DEPTH);
    int16_t fbuf[RX_FFE_LEN]={0};
    int *dec = malloc(n_export*sizeof(int));     /* decided index per emission */
    int *sdec= malloc(n_export*sizeof(int));     /* memoryless slicer index    */
    int ndec=0;
    for (int n=0;n<n_export;n++){
        int16_t y = ffe_fix(fbuf, taps, stim[n]);
        sdec[n] = (y<0)?0:1;                     /* slicer index, per sample   */
        int d = mlse_step(&ml, y);
        if (d>=0) dec[ndec++]=d;
    }

    /* ── align decisions to ground truth: find lag g minimising errors ── */
    int bestg=0; long bestmis=ndec+1;
    for (int g=-(RX_FFE_LEN+L+8); g<=8; g++){
        long mis=0,c=0;
        for (int d=0; d<ndec; d++){ int t=d+g; if(t<0||t>=NSYM) continue; mis += (dec[d]!=sym[t]); c++; }
        if (c> ndec/2 && mis<bestmis){ bestmis=mis; bestg=g; }
    }
    int SKIP = (bestg<0? -bestg:0) + 4;

    /* eq_expect[d] = ground-truth symbol index for MLSE emission d */
    int8_t *expect = malloc(n_export*sizeof(int8_t));
    for (int d=0; d<n_export; d++){ int t=d+bestg; expect[d]=(t>=0&&t<NSYM)?(int8_t)sym[t]:0; }

    /* self-check: MLSE vs slicer error rate (slicer aligned the same way) */
    long mlse_err=0, mlse_tot=0;
    for (int d=SKIP; d<ndec; d++){ int t=d+bestg; if(t<0||t>=NSYM) continue; mlse_tot++; mlse_err+= (dec[d]!=sym[t]); }
    /* slicer needs its own lag (it has no traceback delay) */
    long slc_err=0, slc_tot=0; int sg=0; long sb=n_export+1;
    for (int g=-(RX_FFE_LEN+L+8); g<=8; g++){ long mis=0,c=0;
        for(int n=0;n<n_export;n++){int t=n+g; if(t<0||t>=NSYM)continue; mis+=(sdec[n]!=sym[t]); c++;}
        if(c>n_export/2 && mis<sb){sb=mis;sg=g;} }
    for (int n=SKIP;n<n_export;n++){int t=n+sg; if(t<0||t>=NSYM)continue; slc_tot++; slc_err+=(sdec[n]!=sym[t]);}

    /* ── emit header ─────────────────────────────────────────────────── */
    FILE *fp=fopen(out,"w"); if(!fp){perror("open");return 1;}
    fprintf(fp,
      "/* serdes_workload.h — AUTO-GENERATED by bridge (synthetic channel).\n"
      " * Receiver data for the firmware: FFE taps shaped to a PR target,\n"
      " * the PR target + symbol levels for the MLSE, the captured RX sample\n"
      " * stream, and the aligned ground-truth symbols for self-check.\n"
      " * Regenerate; do not hand-edit. */\n"
      "#ifndef SERDES_WORKLOAD_H\n#define SERDES_WORKLOAD_H\n#include <stdint.h>\n\n");
    fprintf(fp,"#define EQ_FFE_LEN  %d\n#define EQ_TAP_Q    %d\n#define EQ_SAMP_Q   %d\n",
            RX_FFE_LEN,TAP_Q,SAMP_Q);
    fprintf(fp,"#define EQ_NSTIM    %d\n\n",n_export);
    fprintf(fp,"#define MLSE_M      %d\n#define MLSE_L      %d\n#define MLSE_NS     %d\n",
            M,L,ns);
    fprintf(fp,"#define MLSE_TB     %d\n#define MLSE_SKIP   %d\n\n",TB_DEPTH,SKIP);

    fprintf(fp,"static const int16_t eq_ffe_taps[EQ_FFE_LEN] = {\n  ");
    for(int k=0;k<RX_FFE_LEN;k++){fprintf(fp,"%6d,",taps[k]); if((k%8)==7)fprintf(fp,"\n  ");}
    fprintf(fp,"\n};\n\nstatic const int16_t mlse_target[MLSE_L] = { ");
    for(int k=0;k<L;k++)fprintf(fp,"%d, ",tgt[k]);
    fprintf(fp,"};\nstatic const int16_t mlse_levels[MLSE_M] = { ");
    for(int j=0;j<M;j++)fprintf(fp,"%d, ",levels[j]);
    fprintf(fp,"};\n\nstatic const int16_t eq_stimulus[EQ_NSTIM] = {\n  ");
    for(int n=0;n<n_export;n++){fprintf(fp,"%6d,",stim[n]); if((n%12)==11)fprintf(fp,"\n  ");}
    fprintf(fp,"\n};\n\nstatic const int8_t eq_expect[EQ_NSTIM] = {\n  ");
    for(int d=0;d<n_export;d++){fprintf(fp,"%d,",expect[d]); if((d%32)==31)fprintf(fp,"\n  ");}
    fprintf(fp,"\n};\n\n#endif /* SERDES_WORKLOAD_H */\n");
    fclose(fp);

    printf("bridge gen: wrote %s\n", out);
    printf("  FFE train MSE=%.3e  lag g=%d  SKIP=%d  decisions=%d\n", mse, bestg, SKIP, ndec);
    printf("  self-check  MLSE errors = %ld / %ld   (%.2e)\n",
           mlse_err, mlse_tot, mlse_tot? (double)mlse_err/mlse_tot : 0.0);
    printf("  self-check  slicer errors = %ld / %ld (%.2e)  <- MLSE should be <= this\n",
           slc_err, slc_tot, slc_tot? (double)slc_err/slc_tot : 0.0);

    free(sym);free(rx);free(stim);free(dec);free(sdec);free(expect);
    return 0;
}

int main(int argc,char**argv){
    if(argc>=2 && !strcmp(argv[1],"gen")){
        const char *out = (argc>=3)? argv[2] : "tasks/serdes_workload.h";
        int n = (argc>=4)? atoi(argv[3]) : 2048;
        return bridge_gen(out,n);
    }
    fprintf(stderr,"usage: %s gen [out.h] [n_export]\n",argv[0]);
    return 1;
}
