#include "rtklib.h"

#include <stdlib.h>

extern int truth_open(rtk_t* rtk, const char* file);
extern int truth_read(rtk_t* rtk);

int truth_open(rtk_t* rtk, const char* file)
{
    char line[1024];
    int i;

    if (!(rtk->truth.fp = fopen(file, "r"))) {
        trace(1, "truth file open error: %s\n", file);
        return 0;
    }

    /* skip header lines (1–59) */
    for (i = 0; i < 59; i++) {
        if (!fgets(line, sizeof(line), rtk->truth.fp)) {
            fclose(rtk->truth.fp);
            rtk->truth.fp = NULL;
            return 0;
        }
    }

    return 1;
}

int truth_read(rtk_t* rtk)
{
    char line[1024];
    double x, y, z, tow;
    int week;

    if (!rtk->truth.fp) return 0;

    while (fgets(line, sizeof(line), rtk->truth.fp)) {

        /* parse only required fields */
        if (sscanf(line,
            "%*[^,],%*[^,],%d,%lf,%lf,%lf,%lf",
            &week, &tow, &x, &y, &z) < 5) {
            continue;
        }

        rtk->truth.week = week;
        rtk->truth.tow = tow;
        rtk->truth.rr[0] = x;
        rtk->truth.rr[1] = y;
        rtk->truth.rr[2] = z;

        return 1;
    }

    return 0;
}
