/* bridge.h — host-side mediator between the simulation and the firmware.
 *
 * The bridge is the ONLY place the two worlds meet. It never runs on the
 * Duo. Two modes:
 *
 *   gen   <s4p> [out.h] [n]   Drive the simulation channel, train the RX FFE,
 *                             and emit tasks/serdes_workload.h (fixed-point
 *                             taps + stimulus) for the firmware to compile in.
 *
 *   cosim                     Read the generated tasks/serdes_workload.h and
 *                             run the firmware's fixed-point kernel against a
 *                             floating-point reference, reporting decision
 *                             agreement and quantisation SNR.
 */
#ifndef BRIDGE_H
#define BRIDGE_H

int bridge_gen  (const char *s4p, const char *out_header, int n_export);
int bridge_cosim(void);

#endif /* BRIDGE_H */
