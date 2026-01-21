#include "rtklib.h"

#include <stdlib.h>

int truth_open(rtk_t* rtk, const char* file);
int truth_read(rtk_t* rtk);

int truth_open(rtk_t* rtk, const char* file)
{
    char line[1024];
    int i;

    if (!(rtk->truth.fp = fopen(file, "r"))) {
        trace(1, "truth file open error: %s\n", file);
        fprintf(stderr, "Truth file open error\n");

        return 0;
    }

    /* skip header lines (1�59) */
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
    double x, y, z, tow, week;

    if (!rtk->truth.fp) {
        fprintf(stderr, "No file\n");

        return 0;
    }
    while (fgets(line, sizeof(line), rtk->truth.fp)) {
        fprintf(stderr, "Scanning line\n");

        /* parse only required fields */
        if (sscanf(line,
            "%*[^,],%*[^,],%lf,%lf,%lf,%lf,%lf",
            &week, &tow, &x, &y, &z) == 5)
        {
            rtk->truth.week = (int)week;
            rtk->truth.tow = tow;
            rtk->truth.rr[0] = x;
            rtk->truth.rr[1] = y;
            rtk->truth.rr[2] = z;
            return 1;
        }
    }

    fprintf(stderr, "Done scanning\n");

    return 0;
}
