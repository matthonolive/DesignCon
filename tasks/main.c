/* main.c — firmware entry. Assembles nothing itself; just hands the
 * task-defined schedule to the generic scheduler. */
#include "bare.h"
#include "scheduler.h"
#include "tasks.h"

void main(void){
    uputs("serdes-rt firmware\n");
    scheduler_run(tasks_schedule());
}
