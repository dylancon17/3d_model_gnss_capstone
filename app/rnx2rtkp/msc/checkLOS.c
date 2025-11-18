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
    //fprintf(stderr, "%s\n", "LOS Update\n");

    int i, nrej=0;
    int rej_idx[MAXOBS];
    const double *rr = rtk->sol.rr; // TODO add filtering for if the estimated position is poor

    // Converts ECEF to LLH
    double pos[3];
    ecef2pos(rr, pos);
    
    // Set relative origin here as it is constant for the rest of the update
    set_relative_origin(&(rtk->opt.dtm), pos[0], pos[1]);

    double e[3], azel[2]; //warning gets overwritten each satellite. Should be fine?
    double r;
    boolean reject = 0;

    //fprintf(stderr, "%i %i", rtk->opt.dtm.max_distance_check, rtk->opt.dtm.reject_observations);

    for (i = 0;i < *ns && i < MAXOBS;i++) {
        r = geodist(rs + i * 6, rr, e); //TODO how is rs indexed
            /* geodist failure check */
        if (r <= 0) {
            /* Bad geometry → reject satellite */
            rej_idx[nrej++] = i;            continue;
        }

        satazel(pos, e, azel);
        //fprintf(stderr, "%d", sat[i]);
        reject = !check_los(azel[0], azel[1], pos[0], pos[1], pos[2], &(rtk->opt.dtm));

        if (reject) {
            rej_idx[nrej++] = i;
        }

        // TODO add logging to report on what's happening
    }

    //fprintf(stdout, "%d possible sats\n", *ns);


    // TODO add check if removing too many (indicates incorrect estimated position more likely)
    // Remove the rejected indices
    for (i = nrej - 1; i >= 0; i--) {
        int idx = rej_idx[i];
        memmove(sat + idx, sat + idx + 1, (*ns - idx - 1) * sizeof(int));
        memmove(iu + idx, iu + idx + 1, (*ns - idx - 1) * sizeof(int));
        memmove(ir + idx, ir + idx + 1, (*ns - idx - 1) * sizeof(int));

        (*ns)--;
    }

    //fprintf(stdout, "%d sats rejected\n", nrej);


    return  nrej;
    
}

//Assumes relative origin already set
extern boolean check_los(double sat_az, double sat_elev, double origin_lat, double origin_long, double origin_height, struct DTMData* DTM) {
    //fprintf(stderr, "Checking line of sight for: az: %lf elev: %lf at lat: %lf long: %lf height: %lf\n",
    //    sat_az * 180 / 3.14,
    //    sat_elev * 180 / 3.14,
    //    origin_lat * 180 / 3.14,
    //    origin_long * 180 / 3.14,
    //    origin_height);

    int out_of_bounds = 0;
    double current_DTM_height = 0, sat_height = 0;
    double sat_vertical_slope = tanf(sat_elev);
    int origin_x = 0, origin_y = 0;
    //If the starting height is below the DEM, use the DEM height, or if the option to use_dem_height_only is set
    get_relative_height(DTM, &origin_x, &origin_y, &current_DTM_height, &out_of_bounds);
    if (out_of_bounds != 1 && (current_DTM_height > origin_height || DTM->use_dem_height_only == 1)) {
        origin_height = current_DTM_height + DTM->antenna_dem_offset;
    }

    // Set up height checks
    double max_checked_DTM_height = origin_height;

    // Only search until a distance where you're guranteed to hit a building, or you've searched a ways
    double max_distance_steps = (DTM->max_dem_height - origin_height) / sat_vertical_slope;
    if (max_distance_steps > DTM->max_distance) {
        max_distance_steps = DTM->max_distance;
    }
    // Scale to units of steps instead of meters
    max_distance_steps = max_distance_steps / DTM->step_size;

    // Draws a line from (0,0) to (E1, N1) that's the maximum distance required to check
    int E1 = (int)roundf(max_distance_steps * sinf(sat_az)); // Find distance to end point and then reduce by step size to reduce the number of required steps
    int N1 = (int)roundf(max_distance_steps * cosf(sat_az));
    // fprintf(stderr, "Determining End Point max distance: %i sinf: %lf cosf: %lf float solution: %lf rounded float solution: %lf",
    //      DTM->max_distance_check, sinf(sat_az), cosf(sat_az), DTM->max_distance_check* sinf(sat_az), roundf(DTM->max_distance_check * sinf(sat_az)));
    // Set up Bresenham Algorithm
    int dE = abs(E1), sE = 0 < E1 ? 1 : -1; // Each step represents the step size
    int dN = -abs(N1), sN = 0 < N1 ? 1 : -1; // Defines directions of steps
    LineState line = {0, 0, 0, dE, dN, sE, sN, dE + dN, 0};


    while (1) {
        // Traverse the DTM
        step_along_line(&line);
        
        get_relative_height(DTM, &(line.E), &(line.N), &current_DTM_height, &out_of_bounds);

        //fprintf(stderr, "Line State: E: %d N: %d d: %lf dE: %d dN: %d sE: %d sN: %d err: %d e2: %d Height Comparison: %lf Boundary: %d\n",
        //    line.E, line.N, line.d, line.dE, line.dN, line.sE, line.sN, line.err, line.e2, current_DTM_height, out_of_bounds);

        // If fully traversed DTM, LOS is clear
        if (out_of_bounds == 1) {return 1;}

        // Have traversed the max distance required to trraverse, LOS is clear. Depends on the number of steps travelled
        if (max_distance_steps < line.d) {
            //fprintf(stderr, "covered max distance\n");
            return 1;
        }

        // If a higher height has already been checked (or lower than starting height), LOS is clear at that point. Continue traversing
        if (current_DTM_height < max_checked_DTM_height) {continue;}
        
        //Else, onto the new max height checked
        max_checked_DTM_height = current_DTM_height;

        // Check the height of the satellite based on the distance travelled along the satellite line
        
        // Calculate precise distance travelled for a height check. Only done here for efficiency
        // Get the distance in meters instead of the units of steps, required for slope uints
        line.d = sqrt(pow(line.E * DTM->step_size, 2) + pow(line.N * DTM->step_size, 2));
        
        sat_height = sat_vertical_slope * line.d + origin_height;

        // Change distance units back to step size instead of meters
        line.d = line.d / DTM->step_size; 

        // Satellite is lower than DTM, LOS is obstructed
        if (sat_height < current_DTM_height) { 
            //fprintf(stderr, "Reject sat, too low");
            return 0; 
        }

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