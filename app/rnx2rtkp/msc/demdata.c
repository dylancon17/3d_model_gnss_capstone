#include "rtklib.h"


/* set_relative_origin ---------------------------------------------------
* set an origin to traverse along the DEM from
* args   : DTMData   *DEM     I   DTM Object, see rtklib.h for definition and rnx2rtkp.c for initial setup
*          double    latitude I   latitude in radians
*          double    longitudeI   longitude in radians
* return : void (even if out of bounds, that will be handled in get relative height call, but we can change this if you need)*/
extern void set_relative_origin(
    struct DTMData* DEM,
    double latitude,
    double longitude) {
    return; //TODO-TC
}

/* get_relative_height ---------------------------------------------------
* Report height for a location with respect to the relative origin
* args   : DTMData   *DEM           I   DTM Object, see rtklib.h for definition and rnx2rtkp.c for initial setup
*          int       *E I           I   The number of steps East to take (distance = E * step_size) (negative indicates West)
*          int       *N             I   The number of steps North to take (distance = N * step_size) (negative indicates South)
*          double    *h             I   The height of the DEM at that point
*          int       *out_of_bounds I   1 if there is not DEM data at that point (h will not looked at if this is 1)
* return : void (data returned by setting h and out_of_bounds)*/
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