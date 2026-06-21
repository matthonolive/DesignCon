/* scheduler.c — cyclic executive. Generic: driven entirely by the
 * schedule_t passed to scheduler_run(). No task-specific code here. */
#include "bare.h"
#include "scheduler.h"

static volatile uint64_t * const mtime    = (uint64_t *)MTIME;
static volatile uint64_t * const mtimecmp = (uint64_t *)MTIMECMP;

static volatile uint32_t s_ticks;     /* ISR writes, executive reads */
static uint64_t          s_interval;  /* mtime ticks per frame       */

uint32_t scheduler_ticks(void){ return s_ticks; }

void trap_handler(void) __attribute__((interrupt("machine"), aligned(4)));
void trap_handler(void){ *mtimecmp += s_interval; s_ticks++; }

static void run_frame(const frame_t *fr, uint32_t idx){
    for(uint8_t i = 0; i < fr->n; i++){
        uint64_t t0   = rdcycle();
        fr->items[i].fn();
        uint64_t used = rdcycle() - t0;
        if(used > fr->items[i].budget){
            uputs("OVERRUN frame="); uputd(idx);
            uputs(" task=");   uputs(fr->items[i].name ? fr->items[i].name : "?");
            uputs(" used=");   uputd(used);
            uputs(" budget="); uputd(fr->items[i].budget);
            uputc('\n');
        }
    }
}

void scheduler_run(const schedule_t *s){
    s_interval = (uint64_t)MTIME_HZ / s->frame_hz;
    s_ticks    = 0;

    write_csr(mtvec, (uint64_t)trap_handler);
    *mtimecmp = *mtime + s_interval;
    set_csr(mie,     1u << 7);   /* MTIE */
    set_csr(mstatus, 1u << 3);   /* MIE  */

    uint32_t last = s_ticks;
    for(;;){
        __asm__ volatile("wfi");
        uint32_t now = s_ticks;
        if(now == last) continue;
        if(now - last > 1)
            { uputs("MISSED "); uputd(now - last - 1); uputs(" frame(s)\n"); }
        run_frame(&s->frames[now % s->n_frames], now % s->n_frames);
        last = now;
        if(s->run_seconds && now >= s->run_seconds * s->frame_hz){
            uputs("done\n");
            plat_exit(0);
        }
    }
}
