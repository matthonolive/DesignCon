/* bare.h — include in exactly one source file (scheduler.c) */
#ifndef BARE_H
#define BARE_H
#include <stdint.h>

/* QEMU virt addresses — lift into a board.h when you port to the C906 */
#define UART0_BASE 0x10000000UL
#define CLINT_BASE 0x02000000UL
#define MTIME      (CLINT_BASE + 0xBFF8)
#define MTIMECMP   (CLINT_BASE + 0x4000)
#define MTIME_HZ   10000000UL          /* virt timebase = 10 MHz */

static inline void uputc(char c){
    volatile uint8_t *u = (uint8_t *)UART0_BASE;
    while(!(u[5] & 0x20)) ;            /* wait for LSR.THRE */
    u[0] = c;
}
static inline void uputs(const char *s){ while(*s) uputc(*s++); }

static inline void uputd(uint64_t v){          /* unsigned decimal */
    char b[20]; int i = 0;
    if(!v){ uputc('0'); return; }
    while(v){ b[i++] = '0' + (v % 10); v /= 10; }
    while(i) uputc(b[--i]);
}
static inline void uputx(uint64_t v){          /* 64-bit hex */
    uputs("0x");
    for(int s = 60; s >= 0; s -= 4)
        uputc("0123456789abcdef"[(v >> s) & 0xf]);
}

#define read_csr(r)      ({ uint64_t v; __asm__ volatile("csrr %0, " #r : "=r"(v)); v; })
#define write_csr(r, x)  __asm__ volatile("csrw " #r ", %0" :: "rK"(x))
#define set_csr(r, x)    __asm__ volatile("csrs " #r ", %0" :: "rK"(x))

static inline uint64_t rdcycle(void){ uint64_t c; __asm__ volatile("rdcycle %0":"=r"(c)); return c; }

#define SIFIVE_TEST     0x100000UL
#define FINISHER_FAIL   0x3333
#define FINISHER_PASS   0x5555

static inline void qemu_exit(int code) __attribute__((noreturn));
static inline void qemu_exit(int code)
{
    volatile uint32_t *finisher = (volatile uint32_t *)SIFIVE_TEST;

    if (code == 0) {
        *finisher = FINISHER_PASS;
    } else {
        *finisher = ((uint32_t)code << 16) | FINISHER_FAIL;
    }

    for (;;) {
        __asm__ volatile("wfi");
    }
}

/* Entry: set stack, clear bss, call main(). RV64 (sd). For rv32 use sw / 4. */
__asm__(
".section .text.init   \n"
".globl _start         \n"
"_start:               \n"
"  la   sp, _stack_top \n"
"  la   t0, _bss_start \n"
"  la   t1, _bss_end   \n"
"1:bgeu t0, t1, 2f     \n"
"  sd   zero, 0(t0)    \n"
"  addi t0, t0, 8      \n"
"  j    1b             \n"
"2:call main           \n"
"3:wfi                 \n"
"  j    3b             \n");
#endif