#include "bare.h"

static volatile uint64_t * const mtime    = (uint64_t *)MTIME;
static volatile uint64_t * const mtimecmp = (uint64_t *)MTIMECMP;
#define FRAME_HZ  1000u                       /* minor-cycle (frame) rate */
#define INTERVAL  (MTIME_HZ / FRAME_HZ)
volatile uint32_t ticks;                      /* ISR writes, executive reads */

void trap_handler(void) __attribute__((interrupt("machine"), aligned(4)));
void trap_handler(void){ *mtimecmp += INTERVAL; ticks++; }


static void equalize_block(void){
}

static void heartbeat(void){
    uint32_t t = ticks;
    if((t % FRAME_HZ) == 0){ uputs("alive "); uputd(t / FRAME_HZ); uputs("s\n"); }
}

typedef void (*work_fn)(void);
typedef struct { work_fn fn; uint32_t budget; } work_item_t;   /* budget in rdcycle units */
typedef struct { const work_item_t *items; uint8_t n; } frame_t;

/* placeholder budgets — set from measured WCET once there's real work */
static const work_item_t frame_a[] = { {equalize_block, 100000}, {heartbeat, 100000} };
static const work_item_t frame_b[] = { {equalize_block, 100000} };
static const frame_t schedule[] = { {frame_a, 2}, {frame_b, 1} };
#define MAJOR (sizeof(schedule)/sizeof(schedule[0]))

static void run_frame(uint32_t idx){
    const frame_t *fr = &schedule[idx];
    for(uint8_t i = 0; i < fr->n; i++){
        uint64_t t0 = rdcycle();
        fr->items[i].fn();
        uint64_t used = rdcycle() - t0;
        if(used > fr->items[i].budget){
            uputs("OVERRUN frame="); uputd(idx);
            uputs(" item="); uputd(i);
            uputs(" used="); uputd(used); uputc('\n');
        }
    }
}

static void exec_run(void){
    uint32_t last = ticks;
    for(;;){
        __asm__ volatile("wfi");              /* sleep until the next tick */
        uint32_t now = ticks;
        if(now == last) continue;             /* spurious wake */
        if(now - last > 1)                    /* fell behind = frame-level overrun */
            { uputs("MISSED "); uputd(now - last - 1); uputs(" frame(s)\n"); }
        run_frame(now % MAJOR);               /* phase locked to absolute frame number */
        last = now;
        if(now >= 5u * FRAME_HZ){ uputs("done\n"); qemu_exit(0); }   /* bounded test run */
    }
}

void main(void){
    uputs("scheduler.c running\n");
    write_csr(mtvec, (uint64_t)trap_handler);
    *mtimecmp = *mtime + INTERVAL;
    set_csr(mie, 1u << 7);                     /* MTIE */
    set_csr(mstatus, 1u << 3);                 /* MIE  */
    exec_run();
}