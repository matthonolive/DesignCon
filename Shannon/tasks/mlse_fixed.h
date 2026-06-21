/*
 * mlse_fixed.h — freestanding fixed-point streaming MLSE (Viterbi) detector.
 *
 * PAM-M over a short partial-response target h_t[0..L-1] (Q(EQ_SAMP_Q)) that
 * the RX FFE shapes the channel toward. Integer only; int64 path metrics
 * (free on RV64). One symbol decided per input sample, delayed by traceback
 * depth TB => bounded, data-independent inner loop => analysable WCET.
 *
 * Generalised for PAM-M with fractional, normalised symbol levels stored in
 * Q(LVL_Q): expected = (sum_d target[d]*level[sym_{n-d}]) >> LVL_Q. The
 * expected-sample table depends only on the trellis (not on the data), so it
 * is built once at init and the per-sample loop is division-free.
 *
 * State = last (L-1) symbol indices, base-M, newest in the LS digit;
 * n_states = M^(L-1). One source of truth for firmware, cosim and bridge.
 */
#ifndef MLSE_FIXED_H
#define MLSE_FIXED_H
#include <stdint.h>

#define MLSE_INF  ((int64_t)1 << 60)
#ifndef MLSE_LMAX
#define MLSE_LMAX 8
#endif
#ifndef MLSE_NSMAX
#define MLSE_NSMAX 256
#endif

typedef struct {
    int64_t       *pm;       /* [ns]      path metrics                        */
    int8_t        *tb_sym;   /* [tb*ns]   survivor: incoming symbol index     */
    int8_t        *tb_pre;   /* [tb*ns]   survivor: predecessor state         */
    int32_t       *exp_tab;  /* [ns*m]    precomputed expected sample, Q(SAMP)*/
    const int16_t *target;   /* [l]       PR target, Q(EQ_SAMP_Q)             */
    const int16_t *levels;   /* [m]       symbol amplitudes, Q(lvl_q)         */
    int32_t        pw[MLSE_LMAX];
    int  ns, m, l, tb, lvl_q, ns_mask;
    int  head, filled;
} mlse_t;

static inline void mlse_init(mlse_t *s, int64_t *pm, int8_t *tb_sym, int8_t *tb_pre,
                             int32_t *exp_tab, const int16_t *target,
                             const int16_t *levels, int ns, int m, int l,
                             int tb, int lvl_q){
    s->pm=pm; s->tb_sym=tb_sym; s->tb_pre=tb_pre; s->exp_tab=exp_tab;
    s->target=target; s->levels=levels;
    s->ns=ns; s->m=m; s->l=l; s->tb=tb; s->lvl_q=lvl_q;
    s->ns_mask = ((ns & (ns-1)) == 0) ? (ns-1) : 0;   /* fast wrap if pow2 */
    s->head=0; s->filled=0;
    s->pw[0]=1;
    for (int d=1; d<l && d<MLSE_LMAX; d++) s->pw[d]=s->pw[d-1]*m;
    /* expected[v][j] = sum_d target[d]*level[ symbol d steps back ], >> lvl_q */
    for (int v=0; v<ns; v++)
        for (int j=0; j<m; j++){
            int64_t e = (int64_t)target[0]*(int64_t)levels[j];
            for (int d=1; d<l; d++){
                int idx = (v / s->pw[d-1]) % m;
                e += (int64_t)target[d]*(int64_t)levels[idx];
            }
            exp_tab[v*m+j] = (int32_t)(e >> lvl_q);
        }
    for (int v=0; v<ns; v++) pm[v]=0;     /* cold start: equiprobable */
}

/* Feed one equalised sample y (Q(EQ_SAMP_Q)). Returns a decided symbol index,
 * or -1 while the traceback ring is still filling. */
static inline int mlse_step(mlse_t *s, int32_t y){
    int ns=s->ns, m=s->m, mask=s->ns_mask;
    int8_t  *bs = s->tb_sym + (long)s->head*ns;
    int8_t  *bp = s->tb_pre + (long)s->head*ns;
    const int32_t *et = s->exp_tab;

    int64_t nm[MLSE_NSMAX];
    for (int v=0; v<ns; v++) nm[v]=MLSE_INF;

    for (int v=0; v<ns; v++){
        int64_t base = s->pm[v];
        if (base >= MLSE_INF) continue;
        const int32_t *ev = et + v*m;
        int vb = m*v;
        for (int j=0; j<m; j++){
            int32_t d = y - ev[j];
            int64_t cand = base + (int64_t)d*d;
            int vp = mask ? ((vb+j) & mask) : ((vb+j) % ns);
            if (cand < nm[vp]){ nm[vp]=cand; bs[vp]=(int8_t)j; bp[vp]=(int8_t)v; }
        }
    }
    int64_t mn=nm[0]; int best=0;
    for (int v=1; v<ns; v++) if (nm[v]<mn){ mn=nm[v]; best=v; }
    for (int v=0; v<ns; v++) s->pm[v]=nm[v]-mn;

    s->head = (s->head+1) % s->tb;
    if (s->filled < s->tb) s->filled++;
    if (s->filled < s->tb) return -1;

    int state=best, idx=(s->head + s->tb - 1) % s->tb;
    for (int t=0; t<s->tb-1; t++){
        state = (s->tb_pre + (long)idx*ns)[state];
        idx = (idx + s->tb - 1) % s->tb;
    }
    return (s->tb_sym + (long)idx*ns)[state];
}

#endif /* MLSE_FIXED_H */