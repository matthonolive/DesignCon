CROSS  ?= riscv64-unknown-elf-
CC      = $(CROSS)gcc
ARCH    = -march=rv64imac_zicsr -mabi=lp64
CFLAGS = $(ARCH) -mcmodel=medany -O2 -g -Wall -ffreestanding -nostdlib -fno-common
ELF     = scheduler.elf
QEMU    = qemu-system-riscv64 -machine virt -cpu rv64 -bios none -nographic \
          -icount shift=0,sleep=off -kernel $(ELF)

define LINK_LD
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
endef
export LINK_LD

all: $(ELF)
$(ELF): scheduler.c bare.h link.ld
	$(CC) $(CFLAGS) -T link.ld -o $@ scheduler.c -lgcc
link.ld:                              # auto-generated; you only edit scheduler.c + bare.h
	@printf '%s' "$$LINK_LD" > $@
run: $(ELF) ; $(QEMU)                 # Ctrl-A X to quit
debug: $(ELF) ; $(QEMU) -s -S         # halts at reset, gdb on :1234
clean: ; rm -f $(ELF) link.ld
.PHONY: all run debug clean