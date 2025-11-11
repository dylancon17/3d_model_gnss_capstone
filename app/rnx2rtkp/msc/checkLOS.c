#include "rtklib.h"

#include <stdlib.h>

#include <math.h>



typedef struct LineState {
    int E, N; //Current point on line
    double d; // Distance travelled along line (not Euclidean)
    int dE, dN; //Absolute step size to each point
    int sE, sN; // Unit direction to step along line 
    int err, e2; // Error tracking
} LineState;

void step_along_line(LineState* l);


/* los update ---------------------------------------------------
* check and update observations based off of line of sight
* args   : rtl_t     *rtk   I   rtk object. Assumed that position is propagated and initialized. Includes link to DEM
*          obsd_t    *level I   rover and base observations. Rover is all first then base. Can also check receiver flag.
*          int       *sat   I   PRN of common satellite[k]
*          int       *iu    I   index of rover sat[k] in obsd_t
*          int       *ir    I   index of base sat[k] in obsd_t
*          int       *ns    I   number of common satellites
*          double    *rs    I   satellite positions in ECEF
* return : int (number of rejected sats) */
extern int los_update(rtk_t* rtk, const obsd_t* obs, int* sat, int* iu, int* ir, int* ns, const double* rs) {
    // TODO confirm how rs is being passed
    int i, nrej=0;
    int rej_idx[MAXOBS];
    const double *rr = rtk->sol.rr; // TODO add filtering for if the estimated position is poor
    // TODO is a check needed if no initial position? Or does the Update step account for that?


    double e[3], azel[2]; //warning gets overwritten each satellite. Should be fine?
    double r;
    boolean reject = 0;

    for (i = 0;i < *ns && i < MAXOBS;i++) {
        r = geodist(rs + i * 6, rr, e); //TODO how is rs indexed
            /* geodist failure check */
        if (r <= 0) {
            /* Bad geometry → reject satellite */
            rej_idx[nrej++] = i;
            continue;
        }

        satazel(rr, e, azel);

        reject = !check_los(azel[0], azel[1], rr[0], rr[1], rr[2], &rtk->opt.dtm);


        if (reject) {
            rej_idx[nrej++] = i;
            trace(3, "Rejected satellite %d: az=%.2f°, el=%.2f°\n", sat[i], azel[0], azel[1]);
        }

        // TODO add logging to report on what's happening
    }

    // TODO add check if removing too many (indicates incorrect estimated position more likely)
    // Remove the rejected indices
    for (i = nrej - 1; i >= 0; i--) {
        int idx = rej_idx[i];
        memmove(sat + idx, sat + idx + 1, (*ns - idx - 1) * sizeof(int));
        memmove(iu + idx, iu + idx + 1, (*ns - idx - 1) * sizeof(int));
        memmove(ir + idx, ir + idx + 1, (*ns - idx - 1) * sizeof(int));

        (*ns)--;
    }

    return  nrej;
    
}


extern boolean check_los(double sat_az, double sat_elev, double origin_lat, double origin_long, double origin_height, struct DTMData* DTM) {
    
    //Set up DTM call
    set_relative_origin(DTM, origin_lat, origin_long);
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
    int dN = -abs(N1), sN = 0 < N1 ? 1 : -1; // Defines directions of steps
    LineState line = {0, 0, 0, dE, dN, sE, sN, dE + dN, 0};


    while (1) {
        // Traverse the DTM
        step_along_line(&line);
        get_relative_height(DTM, &(line.E), &(line.N), &current_DTM_height, &out_of_bounds);

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
        if (DTM->max_distance_check < line.d) { return 1; }

        // Otherwise continue traversing
    }
}

void step_along_line(LineState* l) // credit to https://zingl.github.io/bresenham.html
{
    l->e2 = 2 * l->err;

    // Traverse North/South
    if (l->e2 >= l->dN) {
        l->err += l->dN;
        l->E += l->sE;
    }
    
    // Traverse East/West
    if (l->e2 <= l->dE) {
        l->err += l->dE;
        l->N += l->sN;
    }

    // If traversed East/West and North/South add root 2 distance
    if (l->e2 <= l->dE && l->e2 >= l->dN) {
        l->d += 1.41421356237f; //Root 2
    }
    // Otherwise only travelled one axis
    else {
        l->d += 1;
    }
}