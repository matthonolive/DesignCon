/*
 * mlse_fixed.h — freestanding fixed-point streaming MLSE (Viterbi) detector.
 *
 * Detects PAM-M symbols against a short partial-response target h_t[0..L-1]
 * (Q(EQ_SAMP_Q)) produced by shaping the channel with the RX FFE. Integer
 * only; path metrics are int64 (free on RV64). One symbol decided per input
 * sample, delayed by the traceback depth TB. Bounded, data-independent work
 * per sample => analysable WCET.
 *
 * State layout is the last (L-1) symbol *indices* packed base-M, newest in the
 * least-significant digit. n_states = M^(L-1).
 *
 * The same code is used by the firmware (static buffers, macro-sized), the
 * cosim, and the bridge (heap buffers) — one source of truth, no drift.
 */
#ifndef MLSE_FIXED_H
#define MLSE_FIXED_H
#include <stdint.h>

#define MLSE_INF ((int64_t)1 << 60)
#ifndef MLSE_LMAX
#define MLSE_LMAX 8          /* max target length supported at runtime */
#endif

typedef struct {
    int64_t       *pm;       /* [ns]      path metrics                         */
    int8_t        *tb_sym;   /* [tb*ns]   survivor: incoming symbol index      */
    int8_t        *tb_pre;   /* [tb*ns]   survivor: predecessor state          */
    const int16_t *target;   /* [l]       PR target, Q(EQ_SAMP_Q)              */
    const int16_t *levels;   /* [m]       symbol amplitudes (unit ints)        */
    int32_t        pw[MLSE_LMAX]; /* m^(d) lookup for state-digit extraction   */
    int  ns, m, l, tb;
    int  head;               /* ring write position                            */
    int  filled;             /* steps stored (caps at tb)                      */
} mlse_t;

static inline void mlse_init(mlse_t *s, int64_t *pm, int8_t *tb_sym, int8_t *tb_pre,
                             const int16_t *target, const int16_t *levels,
                             int ns, int m, int l, int tb){
    s->pm = pm; s->tb_sym = tb_sym; s->tb_pre = tb_pre;
    s->target = target; s->levels = levels;
    s->ns = ns; s->m = m; s->l = l; s->tb = tb;
    s->head = 0; s->filled = 0;
    s->pw[0] = 1;
    for (int d = 1; d < l && d < MLSE_LMAX; d++) s->pw[d] = s->pw[d-1] * m;
    for (int v = 0; v < ns; v++) s->pm[v] = 0;   /* cold start: equiprobable */
}

/* Expected noiseless sample for current symbol index j in state v. Q(EQ_SAMP_Q). */
static inline int32_t mlse_expect(const mlse_t *s, int v, int j){
    int32_t e = (int32_t)s->target[0] * (int32_t)s->levels[j];
    for (int d = 1; d < s->l; d++){
        int idx = (v / s->pw[d-1]) % s->m;       /* symbol d steps in the past */
        e += (int32_t)s->target[d] * (int32_t)s->levels[idx];
    }
    return e;   /* target Q15 * unit-level int => stays Q15 */
}

/* Feed one equalized sample y (Q(EQ_SAMP_Q)). Returns a decided symbol index,
 * or -1 while the traceback ring is still filling. */
static inline int mlse_step(mlse_t *s, int32_t y){
    int ns = s->ns, m = s->m;
    int8_t *bs = s->tb_sym + (long)s->head * ns;
    int8_t *bp = s->tb_pre + (long)s->head * ns;

    int64_t nm[ /*ns*/ 256 ];                 /* ns <= 256 supported */
    for (int v = 0; v < ns; v++) nm[v] = MLSE_INF;

    for (int v = 0; v < ns; v++){
        int64_t base = s->pm[v];
        if (base >= MLSE_INF) continue;
        for (int j = 0; j < m; j++){
            int32_t e   = mlse_expect(s, v, j);
            int32_t d   = y - e;
            int64_t bm  = (int64_t)d * (int64_t)d;
            int64_t cand= base + bm;
            int vp = (j + m * v) % ns;
            if (cand < nm[vp]){ nm[vp] = cand; bs[vp] = (int8_t)j; bp[vp] = (int8_t)v; }
        }
    }
    /* normalise to keep metrics bounded */
    int64_t mn = nm[0]; int best = 0;
    for (int v = 1; v < ns; v++) if (nm[v] < mn){ mn = nm[v]; best = v; }
    for (int v = 0; v < ns; v++) s->pm[v] = nm[v] - mn;

    s->head = (s->head + 1) % s->tb;
    if (s->filled < s->tb) s->filled++;

    if (s->filled < s->tb) return -1;

    /* trace back tb steps from the current best state */
    int state = best;
    int idx   = (s->head + s->tb - 1) % s->tb;     /* most recent stored step */
    for (int t = 0; t < s->tb - 1; t++){
        int8_t *p = s->tb_pre + (long)idx * ns;
        state = p[state];
        idx = (idx + s->tb - 1) % s->tb;
    }
    int8_t *sym = s->tb_sym + (long)idx * ns;
    return sym[state];                             /* symbol decided tb ago */
}

#endif /* MLSE_FIXED_H */
