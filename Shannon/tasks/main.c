/* main.c — firmware entry. Hands the schedule to the generic executive. */
#include "scheduler.h"
#include "tasks.h"
int main(void){ scheduler_run(tasks_schedule()); return 0; }
