OUTPUT_ARCH(riscv)
ENTRY(_start)
SECTIONS {
  . = 0x80000000;
  .text : { KEEP(*(.text.init)) *(.text*) }
  .rodata : { *(.rodata*) }
  .data : { *(.data*) }
  . = ALIGN(8); _bss_start = .;
  .bss : { *(.bss*) *(COMMON) }
  _bss_end = .;
  . = ALIGN(16); . += 0x4000; _stack_top = .;
}
