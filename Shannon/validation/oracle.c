/* oracle.c — independent correctness check for the MLSE.
 * Compares a block Viterbi (the ACS+traceback under test, float) against an
 * EXHAUSTIVE maximum-likelihood search over every possible symbol sequence,
 * on identical noisy data with a forced known start state. If the Viterbi is
 * truly optimal they must produce identical decisions on every block.
 *   gcc -O2 -Isimulation -o oracle bridge/oracle.c simulation/serdes_synth.c -lm */
#include "serdes_synth.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

static int M,L,ns; static double htf[8];
#define K     10            /* block length (symbols)            */
#define PFX   (L-1)         /* known prefix => fixed start state */

static double expct(const int*s,int n){ double e=0;
    for(int k=0;k<L;k++){int m=n-k; if(m>=0) e+=htf[k]*syn_amp(s[m]);} return e; }

/* exhaustive ML over suffix s[PFX..K-1]; prefix fixed; returns best seq in out */
static void brute(const double*y,const int*pfx,int*out){
    int suf=K-PFX; long total=1; for(int i=0;i<suf;i++) total*=M;
    double best=1e300; int s[16];
    for(int i=0;i<PFX;i++) s[i]=pfx[i];
    for(long code=0; code<total; code++){
        long c=code; for(int i=0;i<suf;i++){ s[PFX+i]=c%M; c/=M; }
        double cost=0; for(int n=PFX;n<K;n++){ double d=y[n]-expct(s,n); cost+=d*d; }
        if(cost<best){ best=cost; memcpy(out,s,K*sizeof(int)); }
    }
}

/* block Viterbi (the logic under test), forced start state, processes the
 * post-prefix positions n=PFX..K-1 only (prefix defines the start state). */
static void viterbi(const double*y,int start,int*out){
    double pm[256], et[256*8]; int pre[K][256], sym[K][256]; int pw[8];
    pw[0]=1; for(int d=1;d<L;d++)pw[d]=pw[d-1]*M;
    for(int v=0;v<ns;v++)for(int j=0;j<M;j++){ double e=htf[0]*syn_amp(j);
        for(int d=1;d<L;d++){int idx=(v/pw[d-1])%M; e+=htf[d]*syn_amp(idx);} et[v*M+j]=e; }
    for(int v=0;v<ns;v++) pm[v]=(v==start)?0:1e300;
    for(int n=PFX;n<K;n++){ double nm[256]; for(int v=0;v<ns;v++)nm[v]=1e300;
        for(int v=0;v<ns;v++){ if(pm[v]>=1e299)continue;
            for(int j=0;j<M;j++){ double d=y[n]-et[v*M+j],c=pm[v]+d*d; int vp=(M*v+j)%ns;
                if(c<nm[vp]){nm[vp]=c; pre[n][vp]=v; sym[n][vp]=j;} } }
        for(int v=0;v<ns;v++) pm[v]=nm[v]; }
    int best=0; for(int v=1;v<ns;v++) if(pm[v]<pm[best]) best=v;
    int st=best; for(int n=K-1;n>=PFX;n--){ out[n]=sym[n][st]; st=pre[n][st]; }
}

int main(void){
    syn_pr_target(htf,&L); M=SYN_M; ns=1; for(int i=0;i<L-1;i++) ns*=M;
    int BLOCKS=2000; double sigma=0.06; long mismatch=0, compared=0; int diffblocks=0;
    srand(31337);
    for(int b=0;b<BLOCKS;b++){
        int s[16], pfx[8]; for(int i=0;i<K;i++) s[i]=rand()%M;
        for(int i=0;i<PFX;i++) pfx[i]=s[i];
        int start=0; for(int d=1;d<L;d++) start += s[PFX-d]* (int)pow(M,d-1); /* state from prefix */
        double y[16]; for(int n=0;n<K;n++){
            double nz=0; { double u1=(rand()+1.0)/(RAND_MAX+2.0),u2=(rand()+1.0)/(RAND_MAX+2.0);
                           nz=sqrt(-2*log(u1))*cos(2*M_PI*u2)*sigma; }
            y[n]=expct(s,n)+nz; }
        int ob[16], ov[16]; brute(y,pfx,ob); viterbi(y,start,ov);
        int diff=0;
        for(int n=PFX;n<K;n++){ compared++; if(ob[n]!=ov[n]){ mismatch++; diff=1; } }
        if(diff) diffblocks++;
    }
    printf("oracle: exhaustive ML vs Viterbi  (PAM-%d, %d-state, %d blocks, sigma=%.2f)\n",
           M,ns,BLOCKS,sigma);
    printf("  symbols compared = %ld   Viterbi != exhaustive-ML = %ld   (blocks differing %d)\n",
           compared, mismatch, diffblocks);
    printf("  -> Viterbi ACS+traceback is %s exhaustive ML.\n",
           mismatch==0 ? "PROVABLY OPTIMAL (identical to)" : "NOT identical to");
    return 0;
}