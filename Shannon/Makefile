# ── serdes_rt ────────────────────────────────────────────────────────────
# scheduler/ + tasks/  -> cross-compiled firmware (RISC-V, Milk-V Duo)
# bridge/   + simulation/ -> host tools (sim -> firmware data, fixed/float check)
#
#   make workload   (re)generate tasks/serdes_workload.h from the synthetic link
#   make cosim      validate fixed-point FFE+MLSE vs float reference
#   make firmware   cross-build scheduler+tasks -> scheduler.elf
#   make run        run firmware in QEMU
#   make measure    run with tiny budgets to read each task's 'used=' WCET
#   make clean
# ───────────────────────────────────────────────────────────────────────────
CROSS  ?= riscv64-unknown-elf-
CC      = $(CROSS)gcc
HOSTCC ?= gcc

ARCH    = -march=rv64imac_zicsr -mabi=lp64
CFLAGS  = $(ARCH) -mcmodel=medany -O2 -g -Wall -ffreestanding -nostdlib -fno-common
INC     = -Ischeduler -Itasks

EQ_BLOCK    ?=
FFE_BUDGET  ?=
MLSE_BUDGET ?=
DEFS = $(if $(EQ_BLOCK),-DEQ_BLOCK=$(EQ_BLOCK)) \
       $(if $(FFE_BUDGET),-DFFE_BUDGET=$(FFE_BUDGET)) \
       $(if $(MLSE_BUDGET),-DMLSE_BUDGET=$(MLSE_BUDGET))

FW_SRC   = scheduler/scheduler.c scheduler/start.c tasks/tasks.c tasks/main.c
ELF      = scheduler.elf
WORKLOAD = tasks/serdes_workload.h
NSTIM   ?= 2048

KERNELS  = tasks/eq_fixed.h tasks/mlse_fixed.h
QEMU = qemu-system-riscv64 -machine virt -cpu rv64 -bios none -nographic \
       -icount shift=0,sleep=off -kernel $(ELF)

.PHONY: all firmware run measure workload cosim clean FORCE
all: firmware

# ── host: generate the boundary artifact ──────────────────────────────────
bridge/bridge: bridge/bridge.c simulation/serdes_synth.c simulation/serdes_synth.h \
               bridge/bridge.h tasks/mlse_fixed.h
	$(HOSTCC) -O2 -Isimulation -Itasks -o $@ bridge/bridge.c simulation/serdes_synth.c -lm

workload: bridge/bridge
	./bridge/bridge gen $(WORKLOAD) $(NSTIM)

# ── host: fixed-point vs float sign-off ────────────────────────────────────
bridge/cosim: bridge/bridge_cosim.c $(WORKLOAD) $(KERNELS)
	$(HOSTCC) -O2 -Itasks -o $@ bridge/bridge_cosim.c -lm

cosim: bridge/cosim
	./bridge/cosim

# ── firmware ───────────────────────────────────────────────────────────────
firmware: $(ELF)
$(ELF): $(FW_SRC) $(WORKLOAD) scheduler/bare.h scheduler/scheduler.h \
        tasks/tasks.h $(KERNELS) link.ld FORCE
	$(CC) $(CFLAGS) $(DEFS) $(INC) -T link.ld -o $@ $(FW_SRC) -lgcc
FORCE:

link.ld:
	@printf '%s\n' \
	  'OUTPUT_ARCH(riscv)' 'ENTRY(_start)' 'SECTIONS {' \
	  '  . = 0x80000000;' \
	  '  .text : { KEEP(*(.text.init)) *(.text*) }' \
	  '  .rodata : { *(.rodata*) }' '  .data : { *(.data*) }' \
	  '  . = ALIGN(8); _bss_start = .;' '  .bss : { *(.bss*) *(COMMON) }' \
	  '  _bss_end = .;' '  . = ALIGN(16); . += 0x4000; _stack_top = .;' \
	  '}' > link.ld

run: firmware
	$(QEMU)

measure: 
	$(MAKE) firmware FFE_BUDGET=1 MLSE_BUDGET=1
	$(QEMU) | head -4

clean:
	rm -f $(ELF) link.ld bridge/bridge bridge/cosim
