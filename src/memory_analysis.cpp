//
// Created by sjf on 3/14/2022.
//

#include "memory_analysis.h"

unsigned int get_proc_mem() {
    int pid = getpid();
    char file_name[64] = { 0 };
    FILE* fd;
    char line_buff[512] = { 0 };
    sprintf(file_name, "/proc/%d/status", pid);
    fd = fopen(file_name, "r");
    if (fd == NULL) return 0;
    char name[64];
    int vmrss = 0;
    int i;
    for (i = 0; i < VMRSS_LINE - 1; i++) {
        fgets(line_buff, sizeof(line_buff), fd);
    }
    fgets(line_buff, sizeof(line_buff), fd);
    sscanf(line_buff, "%s %d", name, &vmrss);
    fclose(fd);
    return vmrss;
}

