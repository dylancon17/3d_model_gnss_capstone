#include "rtklib.h"

#include <math.h>



typedef struct LineState {
    int E, N; //Current point on line
    double d; // Distance travelled along line (not Euclidean)
    int dE, dN; //Absolute step size to each point
    int sE, sN; // Unit direction to step along line 
    int err, e2; // Error tracking
} LineState;

void step_along_line(LineState* l);


//void los_update() #TODO-DC set imports
//
//for sat in RTKLibSatObject{// TODO-DC properly parse Sat
//    obstructed = check_los(sat_az, sat_elev, origin_lat, origin_long, origin_height, threshold, DTMClass);
//    if not check_los(sat_az, sat_elev, origin_lat, origin_long, origin_height, threshold, DTMClass) {
//        RTKLibObject.removeobservation(sat);
//    }
//}

void los_update(rtk_t* rtk, const obsd_t* obs, int* sat, int* iu, int* ir, int *ns, const double* rs) {
    // TODO confirm how rs is being passed
    int i;
    double rr[3] = rtk->x; // TODO properly declare this. Add in check for if RTK position doesn't yet exist! Confirm coordinate system. TODO, add filtering if initial position SD if terrible
    double *e, *azel; //warning get's overwritten each satellite. Should be fine?
    double r;
    boolean reject = 0;
    for (i = 0;i < ns && i < MAXOBS;i++) {
        geodist(rs + i * 6, rr, e); // TODO add logic here for if this fails
        satazel(rr, e, azel); // TODO add logic here for if this fails

        reject = !check_los(azel[0], azel[1], rr[0], rr[1], rr[2], DEM); //TODO figure out how to pass the DEM class. Maybe through opt?

        if (reject) {
            ns = ns - 1;
            //TODO remove and shift the indices of sat, iu, ir
        }
        
        // TODO add logging to report on what's happening
    }


extern boolean check_los(double sat_az, double sat_elev, double origin_lat, double origin_long, double origin_height, DTMData* DTM) {
    
    //Set up DTM call
    DTM->set_relative_origin(DTM, origin_lat, origin_long);
    boolean out_of_bounds = 0;
    

    //TODO use DEM height if > than origin_height
    // Set up height checks
    double max_checked_DTM_height = origin_height;
    double current_DTM_height = 0, sat_height = 0;
    double sat_vertical_slope = tanf(sat_elev);

    // Draws a line from (0,0) to (E1, N1) that's the maximum distance required to check
    int E1 = (int)roundf(DTM->max_distance_check * sinf(sat_az));
    int N1 = (int)roundf(DTM->max_distance_check * cosf(sat_az));

    // Set up Bresenham Algorithm
    int dE = abs(E1), sE = 0 < E1 ? 1 : -1;
    int dN = -abs(N1), sN = 0 < N1 ? 1 : -1;; // Defines directions of steps
    LineState line = {0, 0, 0, dE, dN, sE, sN, dE + dN, 0};


    while (1) {
        // Traverse the DTM
        step_along_line(&line);
        DTM->get_relative_height(DTM, &(line.E), &(line.N), &current_DTM_height, &out_of_bounds);

        // If fully traversed DTM, LOS is clear
        if (out_of_bounds) {return 1;}

        // If a higher height has already been checked (or lower than starting height), LOS is clear at that point. Continue traversing
        if (current_DTM_height < max_checked_DTM_height) {continue;}
        
        //Else, onto the new max height checked
        max_checked_DTM_height = current_DTM_height;

        // Check the height of the satellite based on the distance travelled along the satellite line
        sat_height = sat_vertical_slope * line.d + origin_height;

        // Satellite is lower than DTM, LOS is obstructed
        if (sat_height < current_DTM_height) { return 0; }

        // Have traversed the max distance required to trraverse, LOS is clear
        if (line.E == E1 && line.N == N1) { return 1; }

        // Otherwise continue traversing
    }
}

void step_along_line(LineState* l) // credit to https://zingl.github.io/bresenham.html
{
    l->e2 = 2 * l->err;

    if (l->e2 >= l->dN) {
        l->err += l->dN;
        l->E += l->sE;
    }
    
    if (l->e2 <= l->dE) {
        l->err += l->dE;
        l->N += l->sN;
    }

    if (l->e2 <= l->dE && l->e2 >= l->dN) {
        l->d += 1.41421356237f; //Root 2
    }
    else {
        l->d += 1;
    }
}