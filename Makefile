# ── serdes_rt ────────────────────────────────────────────────────────────
# Folders:
#   scheduler/  generic cyclic executive (+ bare.h)         ─┐ cross-compiled
#   tasks/      task set + FFE kernel + generated data       ─┘ to RISC-V (Duo)
#   bridge/     host: sim -> firmware data, fixed/float check ─┐ host tools
#   simulation/ host: s4p channel + link sim                  ─┘
#
# Common targets:
#   make firmware   cross-build scheduler+tasks  -> scheduler.elf
#   make run        run firmware in QEMU
#   make measure    run with budget=1 to read the equalizer's 'used=' WCET
#   make workload   (re)generate tasks/serdes_workload.h from the s4p
#   make cosim      validate fixed-point port vs float
#   make sim        build the host serdes simulation
#   make clean
# ───────────────────────────────────────────────────────────────────────────

CROSS ?= riscv64-unknown-elf-
CC     = $(CROSS)gcc
HOSTCC ?= gcc

ARCH   = -march=rv64imac_zicsr -mabi=lp64
CFLAGS = $(ARCH) -mcmodel=medany -O2 -g -Wall -ffreestanding -nostdlib -fno-common
INC    = -Ischeduler -Itasks

# firmware tunables (override on the command line, e.g. make run EQ_BLOCK=512)
EQ_BLOCK  ?=
EQ_BUDGET ?=
DEFS = $(if $(EQ_BLOCK),-DEQ_BLOCK=$(EQ_BLOCK)) $(if $(EQ_BUDGET),-DEQ_BUDGET=$(EQ_BUDGET))

FW_SRC = scheduler/scheduler.c scheduler/start.c tasks/tasks.c tasks/main.c
ELF    = scheduler.elf
WORKLOAD = tasks/serdes_workload.h
S4P    ?= simulation/Asic_Mezz_Retimer_L10_Thru.s4p

QEMU = qemu-system-riscv64 -machine virt -cpu rv64 -bios none -nographic \
       -icount shift=0,sleep=off -kernel $(ELF)

.PHONY: all firmware run measure workload cosim sim clean FORCE
all: firmware

# ── firmware (scheduler + tasks) ──────────────────────────────────────────
firmware: $(ELF)
$(ELF): $(FW_SRC) $(WORKLOAD) scheduler/bare.h scheduler/scheduler.h tasks/tasks.h \
        tasks/eq_fixed.h link.ld FORCE
	$(CC) $(CFLAGS) $(DEFS) $(INC) -T link.ld -o $@ $(FW_SRC) -lgcc

FORCE:

link.ld:
	@printf '%s\n' \
	  'OUTPUT_ARCH(riscv)' 'ENTRY(_start)' 'SECTIONS {' \
	  '  . = 0x80000000;' \
	  '  .text : { KEEP(*(.text.init)) *(.text*) }' \
	  '  .rodata : { *(.rodata*) }' '  .data : { *(.data*) }' \
	  '  . = ALIGN(8); _bss_start = .;' '  .bss : { *(.bss*) *(COMMON) }' \
	  '  _bss_end = .;' '  . = ALIGN(16); . += 0x4000; _stack_top = .;' '}' > $@

run: firmware ; $(QEMU)

# budget=1 forces every equalize call to overrun so 'used=' is printed
measure: $(WORKLOAD)
	$(CC) $(CFLAGS) -DEQ_BUDGET=1 $(if $(EQ_BLOCK),-DEQ_BLOCK=$(EQ_BLOCK)) $(INC) \
	  -T link.ld -o $(ELF) $(FW_SRC) -lgcc
	@echo "running 5s; peak 'used=' is your WCET ..."
	@timeout 15 $(QEMU) 2>/dev/null | grep 'task=equalize' | sort -u | tail -4 || true

# ── bridge / mediation ────────────────────────────────────────────────────
workload: $(WORKLOAD)
$(WORKLOAD): bridge/bridge $(S4P)
	cd bridge && ./bridge gen ../$(S4P) ../$(WORKLOAD) 2048

bridge/bridge: bridge/bridge.c simulation/s4p_channel.c simulation/s4p_channel.h \
               simulation/serdes_sim.h bridge/bridge.h
	$(HOSTCC) -O2 -Wall -Isimulation -o $@ bridge/bridge.c simulation/s4p_channel.c -lm

cosim: bridge/cosim ; cd bridge && ./cosim
bridge/cosim: bridge/bridge_cosim.c $(WORKLOAD) tasks/eq_fixed.h
	$(HOSTCC) -O2 -Wall -Itasks -o $@ bridge/bridge_cosim.c -lm

# ── host simulation ───────────────────────────────────────────────────────
sim: simulation/serdes_sim simulation/s4p2taps
simulation/serdes_sim: simulation/serdes_sim.c simulation/s4p_channel.c
	$(HOSTCC) -O2 -o $@ $^ -lm
simulation/s4p2taps: simulation/s4p2taps.c simulation/s4p_channel.c
	$(HOSTCC) -O2 -o $@ $^ -lm

clean:
	rm -f $(ELF) link.ld bridge/bridge bridge/cosim \
	      simulation/serdes_sim simulation/s4p2taps
	# note: tasks/serdes_workload.h is generated; 'make workload' to rebuild
