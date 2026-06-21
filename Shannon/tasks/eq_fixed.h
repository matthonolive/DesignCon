/*
 * eq_fixed.h — freestanding fixed-point RX FFE kernel (header-only).
 * Integer MACs only (RV64 'M'). Include serdes_workload.h first.
 * samples Q(EQ_SAMP_Q), taps Q(EQ_TAP_Q); product >> EQ_TAP_Q -> Q(EQ_SAMP_Q).
 */
#ifndef EQ_FIXED_H
#define EQ_FIXED_H
#include <stdint.h>

#ifndef EQ_FFE_LEN
#error "include serdes_workload.h before eq_fixed.h"
#endif

typedef struct { int16_t buf[EQ_FFE_LEN]; } eq_state_t;

static inline void eq_reset(eq_state_t *s){
    for (int i = 0; i < EQ_FFE_LEN; i++) s->buf[i] = 0;
}
static inline int16_t eq_sat16(int64_t v){
    if (v >  32767) return  32767;
    if (v < -32768) return -32768;
    return (int16_t)v;
}
/* Shift in one new sample (newest at buf[0]) and return the FFE output. */
static inline int16_t eq_ffe_step(eq_state_t *s, const int16_t *taps, int16_t x){
    for (int i = EQ_FFE_LEN - 1; i > 0; i--) s->buf[i] = s->buf[i - 1];
    s->buf[0] = x;
    int64_t acc = 0;
    for (int k = 0; k < EQ_FFE_LEN; k++)
        acc += (int32_t)taps[k] * (int32_t)s->buf[k];
    return eq_sat16(acc >> EQ_TAP_Q);
}
/* Memoryless slicer (kept for the MLSE-vs-slicer comparison). NRZ sign only. */
static inline int eq_slice(int16_t y){ return (y < 0) ? -1 : +1; }

#endif /* EQ_FIXED_H */
