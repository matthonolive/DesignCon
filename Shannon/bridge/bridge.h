/* bridge.h — host mediator between simulation and firmware. Never runs on
 * the Duo. gen: synthesise/train -> emit tasks/serdes_workload.h. */
#ifndef BRIDGE_H
#define BRIDGE_H
int bridge_gen(const char *out_header, int n_export);
#endif
