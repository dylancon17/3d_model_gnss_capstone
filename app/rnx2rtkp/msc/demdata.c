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
    float* h,
    boolean* out_of_bounds) {
    //TODO-TC Currently Fake data just for testing

    if (*E < -100 && *N < -100) {
        *out_of_bounds = 1;
        return;
    }

    *h = 1000.0f + *E + *N;

    return;
}