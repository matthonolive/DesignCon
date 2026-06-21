/* tasks.c — the receiver datapath, run entirely as scheduler tasks.
 *
 *   task_rx_ffe : fixed-point RX FFE over a block of captured RX samples,
 *                 shaping the channel toward the partial-response target.
 *   task_mlse   : fixed-point streaming Viterbi MLSE over the equalised
 *                 block; emits symbol decisions and checks them against the
 *                 baked-in ground truth (errors counted, not just a sink).
 *   task_report : liveness + cumulative MLSE symbol-error report, 1 Hz.
 *
 * FFE and MLSE state run continuously within a pass over the stimulus and
 * reset together at each wrap, so the firmware reproduces the exact decision
 * sequence the bridge/cosim verified (0 symbol errors).
 *
 * Budgets are QEMU icount values for regression only; re-measure on C906
 * silicon (icount cycles != hardware cycles).
 */
#include "bare.h"
#include "scheduler.h"
#include "tasks.h"
#include "serdes_workload.h"
#include "eq_fixed.h"
#include "mlse_fixed.h"

/* ── tunables ─────────────────────────────────────────────────────────── */
#define FRAME_HZ     1000u
#define RUN_SECONDS  5u
#ifndef EQ_BLOCK
#define EQ_BLOCK     256u
#endif
#ifndef FFE_BUDGET
#define FFE_BUDGET   53000u
#endif
#ifndef MLSE_BUDGET
#define MLSE_BUDGET  230000u
#endif
#define HB_BUDGET    20000u

/* a block must not straddle the stimulus wrap */
typedef char eqblock_divides_nstim[(EQ_NSTIM % EQ_BLOCK == 0) ? 1 : -1];

/* ── shared datapath state ────────────────────────────────────────────── */
static eq_state_t eq_state;
static int16_t    ffe_block[EQ_BLOCK];      /* FFE -> MLSE handoff           */
static uint32_t   eq_pos;                    /* index into eq_stimulus        */

static mlse_t     mlse;
static int64_t    mlse_pm[MLSE_NS];
static int8_t     mlse_tbs[MLSE_TB * MLSE_NS];
static int8_t     mlse_tbp[MLSE_TB * MLSE_NS];
static uint32_t   dec_idx;                   /* decisions emitted this pass   */

static uint32_t   pass_errors, pass_decisions;
static uint64_t   total_errors, total_decisions;
static uint32_t   passes;
volatile uint32_t mlse_sink;                 /* defeat dead-code elimination  */

/* ── FFE task ─────────────────────────────────────────────────────────── */
static void task_rx_ffe(void){
    if(eq_pos == 0) eq_reset(&eq_state);     /* start of a fresh pass         */
    uint32_t p = eq_pos;
    for(uint32_t i = 0; i < EQ_BLOCK; i++)
        ffe_block[i] = eq_ffe_step(&eq_state, eq_ffe_taps, eq_stimulus[p++]);
    eq_pos = p;
}

/* ── MLSE task ────────────────────────────────────────────────────────── */
static void task_mlse(void){
    if(eq_pos - EQ_BLOCK == 0){              /* first block of the pass       */
        mlse_init(&mlse, mlse_pm, mlse_tbs, mlse_tbp,
                  mlse_target, mlse_levels, MLSE_NS, MLSE_M, MLSE_L, MLSE_TB);
        dec_idx = 0; pass_errors = 0; pass_decisions = 0;
    }
    int32_t acc = 0;
    for(uint32_t i = 0; i < EQ_BLOCK; i++){
        int d = mlse_step(&mlse, ffe_block[i]);
        if(d >= 0){
            if(dec_idx >= MLSE_SKIP){
                pass_decisions++;
                if(d != eq_expect[dec_idx]) pass_errors++;
            }
            dec_idx++;
            acc += d;
        }
    }
    mlse_sink += (uint32_t)acc;

    if(eq_pos == EQ_NSTIM){                   /* pass complete                 */
        total_errors    += pass_errors;
        total_decisions += pass_decisions;
        passes++;
        eq_pos = 0;
    }
}

/* ── report task ──────────────────────────────────────────────────────── */
static void task_report(void){
    uint32_t t = scheduler_ticks();
    if((t % FRAME_HZ) == 0){
        uputs("alive "); uputd(t / FRAME_HZ);
        uputs("s  passes="); uputd(passes);
        uputs("  mlse_err="); uputd(total_errors);
        uputc('/');           uputd(total_decisions);
        uputc('\n');
    }
}

/* ── schedule: one frame = FFE -> MLSE -> report ──────────────────────── */
static const work_item_t frame_main[] = {
    { task_rx_ffe, FFE_BUDGET,  "rx_ffe" },
    { task_mlse,   MLSE_BUDGET, "mlse"   },
    { task_report, HB_BUDGET,   "report" },
};
static const frame_t frames[] = {
    { frame_main, 3 },
};
static const schedule_t schedule = {
    frames, 1, FRAME_HZ, RUN_SECONDS
};

const schedule_t *tasks_schedule(void){ return &schedule; }
