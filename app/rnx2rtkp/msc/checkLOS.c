#define _USE_MATH_DEFINES

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
void apply_probability(double* y, double* probability);
void determine_DTM_height_var(double* var, struct DSMData* DSM);
void determine_sat_height_var(double* var, double origin_horizontal_variance, double origin_vertical_variance, double sat_vertical_slope, struct DSMData* DTM);
void soltocov_rtk(sol_t* sol, double* P);

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
    fprintf(stderr, "%s", "LOS Update\n");

    int i, nrej=0;
    int rej_idx[MAXOBS];
    const double *rr = rtk->sol.rr; // TODO add filtering for if the estimated position is poor

    // Converts ECEF to LLH
    double pos[3];
    ecef2pos(rr, pos);
    
    // Converts ECEF covars to LLH
    double P[9]; // 3x3 ENU
    double Q[9];
    
    soltocov_rtk(&(rtk->sol), P);
    covenu(pos, P, Q);

    lat_long ll = { pos[0] * 180 / M_PI, pos[1] * 180 / M_PI};

    int out_of_bounds = 0;
    // Set relative origin here as it is constant for the rest of the update
    // TODO, out of bounds filtering may be required...can be used as an optimization. Must reset the scaling if done so
    set_relative_origin(&(rtk->opt.DSM),&ll,&(rtk->opt.UTM), &(rtk->opt.ellip), &out_of_bounds);

    double e[3], azel[2]; //warning gets overwritten each satellite. Should be fine?
    double r;
    double probability_of_obstruction = 0;

    double current_DTM_height;

    if (out_of_bounds == 1) { //If origin is out of bounds, don't search farther than that
        return 0;
    }

    for (i = 0;i < *ns && i < MAXOBS;i++) {
        fprintf(stderr, "Checking satellite %d\n", i);

        r = geodist(rs + i * 6, rr, e); //TODO how is rs indexed

            /* geodist failure check */
        if (r <= 0) {
            /* Bad geometry → reject satellite */
            rej_idx[nrej++] = i;
            // Reset scaling just in case...
            rtk->ssat[sat[i] - 1].obstruction_scaling = 1;
            fprintf(stderr, "Geodist Failed %d\n", i);

            continue;
        }

        satazel(pos, e, azel);
        //fprintf(stderr, "%d", sat[i]);
        fprintf(stderr, "Requesting probability: ");

        probability_of_obstruction = check_los(azel[0], azel[1], pos[0], pos[1], pos[2], Q[0] + Q[4], Q[8], &(rtk->opt.DSM));
        fprintf(stderr, "%lf\n", probability_of_obstruction);

        // Reject the signal if it's likely that it's obstructed (works for DEM processing options 1 and 2)
        if (probability_of_obstruction > rtk->opt.DSM.rejection_threshold) {

            rej_idx[nrej++] = i;
        }

        double scaling = 1 / (1 - probability_of_obstruction);
        if (scaling > rtk->opt.DSM.max_noise_scaling)
        {
            scaling = rtk->opt.DSM.max_noise_scaling;
        }

        // Save data. Warning! This will last across epochs unless overwritten (actually ssat is reset every epoch). Shouldn't be an issue unless implementation changed
        fprintf(stderr, "%lf", scaling);
        rtk->ssat[sat[i] - 1].obstruction_scaling = scaling;
    }

    //fprintf(stdout, "%d possible sats\n", *ns);

    // If deterministic or probabilistic rejection
    if (rtk->opt.DSM.processing_type > 0 && rtk->opt.DSM.processing_type != 4) {
        // Remove the rejected indices
        for (i = nrej - 1; i >= 0; i--) {
            int idx = rej_idx[i];
            memmove(sat + idx, sat + idx + 1, (*ns - idx - 1) * sizeof(int));
            memmove(iu + idx, iu + idx + 1, (*ns - idx - 1) * sizeof(int));
            memmove(ir + idx, ir + idx + 1, (*ns - idx - 1) * sizeof(int));

            (*ns)--;
        }
    }

    //fprintf(stdout, "%d sats rejected\n", nrej);


    return  nrej;
    
}

//Assumes relative origin already set
extern double check_los(double sat_az, double sat_elev, double origin_lat, double origin_long, double origin_height, double origin_horizontal_variance, double origin_vertical_variance, struct DSMData* DSM) {
    //fprintf(stderr, "Checking line of sight for: az: %lf elev: %lf at lat: %lf long: %lf height: %lf, result: ",
    //    sat_az * 180 / M_PI,
    //   sat_elev * 180 / M_PI,
    //    origin_lat * 180 / M_PI,
    //    origin_long * 180 / M_PI,
    //    origin_height);

    int out_of_bounds = 0;
    double current_DTM_height = 0, sat_height = 0;
    double sat_vertical_slope = tanf(sat_elev);
    int origin_x = 0, origin_y = 0;
    //If the starting height is below the DEM, use the DEM height, or if the option to use_dem_height_only is set
    get_relative_height(DSM, &origin_x, &origin_y, &current_DTM_height, &out_of_bounds);
    if (out_of_bounds != 1 && (current_DTM_height > origin_height || DSM->use_dem_height_only == 1)) {
        origin_height = current_DTM_height + DSM->antenna_dem_offset;
        origin_vertical_variance = DSM->vertical_point_variance + DSM->antenna_dem_offset_var;
    }

    // Set up height checks
    double max_checked_DTM_height = origin_height;

    // Only search until a distance where you're guranteed to hit a building, or you've searched a ways
    double max_distance_steps = (DSM->max_dsm_height - origin_height) / sat_vertical_slope;
    if (max_distance_steps > DSM->max_distance) {
        max_distance_steps = DSM->max_distance;
    }
    // Scale to units of steps instead of meters
    max_distance_steps = max_distance_steps / DSM->step_size;

    // Draws a line from (0,0) to (E1, N1) that's the maximum distance required to check
    int E1 = (int)roundf(max_distance_steps * sinf(sat_az)); // Find distance to end point and then reduce by step size to reduce the number of required steps
    int N1 = (int)roundf(max_distance_steps * cosf(sat_az));
    // fprintf(stderr, "Determining End Point max distance: %i sinf: %lf cosf: %lf float solution: %lf rounded float solution: %lf",
    //      DTM->max_distance_check, sinf(sat_az), cosf(sat_az), DTM->max_distance_check* sinf(sat_az), roundf(DTM->max_distance_check * sinf(sat_az)));
    // Set up Bresenham Algorithm
    int dE = abs(E1), sE = 0 < E1 ? 1 : -1; // Each step represents the step size
    int dN = -abs(N1), sN = 0 < N1 ? 1 : -1; // Defines directions of steps
    LineState line = {0, 0, 0, dE, dN, sE, sN, dE + dN, 0};

    double probability_of_obstruction = 0;
    double DTM_height_var = 0;
    double sat_height_var = 0;
    double y = 0;

    while (1) {
        // Traverse the DTM
        step_along_line(&line);
        
        get_relative_height(DSM, &(line.E), &(line.N), &current_DTM_height, &out_of_bounds);

        //fprintf(stderr, "---DTM Height: %lf, dE: %d: dN %d, out_of_bounds: %d---", current_DTM_height, line.E, line.N, out_of_bounds);

        //fprintf(stderr, "Line State: E: %d N: %d d: %lf dE: %d dN: %d sE: %d sN: %d err: %d e2: %d Height Comparison: %lf Boundary: %d\n",
        //    line.E, line.N, line.d, line.dE, line.dN, line.sE, line.sN, line.err, line.e2, current_DTM_height, out_of_bounds);

        // If fully traversed DTM, LOS is clear after that point, return current probability
        if (out_of_bounds == 1) {
            //fprintf(stderr, "out of bounds %lf\n", probability_of_obstruction);

            return probability_of_obstruction;
        }

        // Have traversed the max distance required to trraverse, LOS is clear after that point. Depends on the number of steps travelled
        if (max_distance_steps < line.d) {
            //fprintf(stderr, "covered max distance %lf\n", probability_of_obstruction);
            return probability_of_obstruction;
        }

        // If a higher height has already been checked, it's not needed to check it again, it's assumed to be low probability of obstruction. NOTE - this is a simplification for efficiency purposes. It also means only one probability per building is estiimated (assuming the building height is ~ constant)
        if (current_DTM_height <= max_checked_DTM_height) {
            continue;
        }
        
        //Else, onto the new max height checked
        max_checked_DTM_height = current_DTM_height;

        // Check the height of the satellite based on the distance travelled along the satellite line
        
        // Calculate precise distance travelled for a height check. Only done here for efficiency
        // Get the distance in meters instead of the units of steps, required for slope uints
        line.d = sqrt(pow(line.E * DSM->step_size, 2) + pow(line.N * DSM->step_size, 2));
        
        sat_height = sat_vertical_slope * line.d + origin_height;

        // Change distance units back to step size instead of meters
        line.d = line.d / DSM->step_size; 


        if (DSM->processing_type > 1) {
            //Calculate probability of obstruction
            determine_sat_height_var(&sat_height_var, origin_horizontal_variance, origin_vertical_variance, sat_vertical_slope, DSM);

            determine_DTM_height_var(&DTM_height_var, DSM);

            y = (current_DTM_height - sat_height) / sqrt(DTM_height_var + sat_height_var);

            apply_probability(&y, &probability_of_obstruction);

            // If hit 0.99 round up to 1.
            if (probability_of_obstruction > 0.99) {
                probability_of_obstruction = 1.0f;
                //fprintf(stderr, "hit max probability %lf\n", probability_of_obstruction);

                return probability_of_obstruction;
            }

            //fprintf(stderr, "Probability Updated To: %lf, Sat Height SD: %lf, DTM Height SD: %lf, y: %lf\n", probability_of_obstruction, sqrt(sat_height_var), sqrt(DTM_height_var), y);


            continue;
        }

        // if processing_type == 1 

        // Boolean rejection (old)
        // Satellite is lower than DTM, LOS is obstructed
        if (sat_height < current_DTM_height) { 
            //fprintf(stderr, "sat height is less than DTM height, %lf < %lf\n", sat_height, current_DTM_height);
            return 1;
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

void apply_probability(double* y, double* probability) {
    double prob = 0.5 * erfc(-*y / sqrt(2.0)); // CDF approximation

    // If probability has not yet been set
    if (*probability == 0.0f) {
        *probability = prob;
        return;
    }

    //If probability has been set multiply the probabilities. Probability should go up given additional probabilities
    *probability = 1 - (1 - prob) * (1 - *probability);
}

void determine_DTM_height_var(double* var, struct DSMData* DSM) {
    //DTM possible error sources:
    // Height accuracy
    // Coordinate Horizontal accuracy - insignificant as all coordinates along the line are called. Can make the assumption that surfaces are generally flat so slight horizontal error doesn't matter
    // Height errros due to calling the nearest height - insignficant as all coordinates along the line are called. Can make the assumption that surfaces are generally flat so no interpolation error would be introduced anyways
    *var = DSM->vertical_point_variance;
}

void determine_sat_height_var(double* var, double origin_horizontal_variance, double origin_vertical_variance, double sat_vertical_slope,DSMData* DSM) {
    //Formula for sat height variance: height = tan(elev) * horizontal distance traversed + origin_height
    //Possible Error Sources:
    // Elevation accuracy - function of sat err, pos err and baseline distance. Assumed to be insignificant due to extremely long baseline
    // Distance Traversed - euclidean distance of starting coordinate and traversed coordinate
    //  Starting Coordinate - Use the horizontal variance of the starting coordinate
    //  Traversed Coordinate - This coordinate can be incorrect due to not perfectly following the line. Bresenham is theoretically accurate to up 1/2 a pixel, arguably more accurate than that. 
    //      Standard deviation will be reflected as 0.68 * Step Size / 2. This is probably overly pessimistic.
    // Origin Height - Will have errors, depends on origin height source
    //  DEM Source - use determine_DTM_height_var + DEM offset var
    //  KF Source - use the estimated filter variance
    *var = origin_vertical_variance;
    double distance_var = pow(0.34 * DSM->step_size,2) + origin_horizontal_variance;
    *var = *var + pow(sat_vertical_slope, 2) * distance_var;

    fprintf(stderr, "Sat Height Var Calculated Using: Origin Vertical Variance: %lf, Origin Horizontal Variance: %lf, Distance Variance: %lf, Calced Variance, %lf\n", origin_vertical_variance, origin_horizontal_variance, distance_var, *var);


    return;
}

/* solution to covariance - copy from solution.c-------------------*/
void soltocov_rtk(sol_t* sol, double* P)
{
    P[0] = sol->qr[0]; /* xx or ee */
    P[4] = sol->qr[1]; /* yy or nn */
    P[8] = sol->qr[2]; /* zz or uu */
    P[1] = P[3] = sol->qr[3]; /* xy or en */
    P[5] = P[7] = sol->qr[4]; /* yz or nu */        
    P[2] = P[6] = sol->qr[5]; /* zx or ue */
}