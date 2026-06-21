/* scheduler.h — generic cyclic executive.
 *
 * The scheduler knows nothing about what tasks do. A caller supplies a
 * schedule (a major cycle of frames, each frame a list of work items with
 * WCET budgets) and the scheduler runs it phase-locked to a timer tick,
 * reporting frame-level misses and per-task budget overruns.
 */
#ifndef SCHEDULER_H
#define SCHEDULER_H
#include <stdint.h>

typedef void (*task_fn)(void);

typedef struct {
    task_fn      fn;
    uint32_t     budget;     /* WCET budget, in rdcycle units      */
    const char  *name;       /* shown in overrun reports           */
} work_item_t;

typedef struct {
    const work_item_t *items;
    uint8_t            n;
} frame_t;

typedef struct {
    const frame_t *frames;       /* major cycle = frames[0 .. n_frames-1] */
    uint8_t        n_frames;
    uint32_t       frame_hz;      /* minor-cycle (frame) rate              */
    uint32_t       run_seconds;   /* stop after N seconds (0 = run forever) */
} schedule_t;

/* Install the timer interrupt and run the executive. Returns only if
 * run_seconds elapses, via the platform exit. */
void     scheduler_run(const schedule_t *s);

/* Minor-cycle count since start (ISR-maintained). Tasks may read this. */
uint32_t scheduler_ticks(void);

#endif /* SCHEDULER_H */
