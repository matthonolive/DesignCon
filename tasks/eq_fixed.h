/*
 * eq_fixed.h  — freestanding fixed-point RX FFE kernel (header-only).
 *
 * No libc, no FP — integer MACs only (uses the RV64 'M' extension).
 * Include serdes_workload.h BEFORE this header (it defines EQ_FFE_LEN /
 * EQ_TAP_Q / EQ_SAMP_Q and the const tap + stimulus arrays).
 *
 * Format: samples are Q(EQ_SAMP_Q), taps are Q(EQ_TAP_Q). A product is
 * Q(EQ_SAMP_Q+EQ_TAP_Q); shifting right by EQ_TAP_Q returns a Q(EQ_SAMP_Q)
 * result. The 64-bit accumulator is free on RV64.
 */
#ifndef EQ_FIXED_H
#define EQ_FIXED_H
#include <stdint.h>

#ifndef EQ_FFE_LEN
#error "include serdes_workload.h before eq_fixed.h"
#endif

typedef struct { int16_t buf[EQ_FFE_LEN]; } eq_state_t;

static inline void eq_reset(eq_state_t *s) {
    for (int i = 0; i < EQ_FFE_LEN; i++) s->buf[i] = 0;
}

static inline int16_t eq_sat16(int64_t v) {
    if (v >  32767) return  32767;
    if (v < -32768) return -32768;
    return (int16_t)v;
}

/* Shift in one new sample (newest at buf[0]) and return the FFE output. */
static inline int16_t eq_ffe_step(eq_state_t *s, const int16_t *taps, int16_t x) {
    for (int i = EQ_FFE_LEN - 1; i > 0; i--) s->buf[i] = s->buf[i - 1];
    s->buf[0] = x;

    int64_t acc = 0;
    for (int k = 0; k < EQ_FFE_LEN; k++)
        acc += (int32_t)taps[k] * (int32_t)s->buf[k];

    return eq_sat16(acc >> EQ_TAP_Q);
}

/* One-bit slicer: returns +1 / -1 for the decision. */
static inline int eq_slice(int16_t y) { return (y < 0) ? -1 : +1; }

#endif /* EQ_FIXED_H */
