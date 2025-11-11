#include "rtklib.h"

extern void set_relative_origin(
    struct DTMData* DEM,
    double latitude,
    double longitude) {
    return; //TODO-TC
}

extern void get_relative_height(
    struct DTMData* DEM,
    int* E,
    int* N,
    double* h,
    int* out_of_bounds) {
    //TODO-TC Currently Fake data just for testing

    if (*E < -3 && *N < -3) {
        //fprintf(stderr, "out of bounds triggered\n");
        *out_of_bounds = 1;
        return;
    }

    *h = 1000.0f + 10 * *E + 20 * *N;

    return;
}