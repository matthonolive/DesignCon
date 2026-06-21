/* tasks.h — the application's task set and its schedule.
 * tasks/ depends only on scheduler/ and the generated serdes_workload.h.
 * No reference to the host bridge or simulation. */
#ifndef TASKS_H
#define TASKS_H
#include "scheduler.h"
const schedule_t *tasks_schedule(void);
#endif
