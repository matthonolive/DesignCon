/* start.c — reset entry. Sets stack, clears .bss, calls main(). */
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
